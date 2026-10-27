"""Step 4 (analysis 4): taxonomy check of her names (suggestions only).

Beetles: her ReportedBeetleSpecies and AcceptedBeetleSpeciesName against the KG's species list (the lab's list, built
on the Johnson et al. 2026 Scolytinae checklist; kg.duckdb views species + names, synonyms included). Exact matching
of the binomial (subspecies trimmed); for names not found, the closest names in the list (Levenshtein ratio) are shown.
Genus / tribe given by her vs the list; authorship vs the list's authority.
Pathogens: every distinct pathogen name -> GBIF backbone species match (api.gbif.org/v1/species/match, strict=false,
kingdom hint from her category; the method the KG's partner resolver uses). Flags: synonym (GBIF accepted name differs),
fuzzy (GBIF matched a different spelling), higher rank only, no match. Also the KG resolver's own match for the same
written name when the KG has one. NOTE: she checked names against NCBI Taxonomy; GBIF and NCBI can disagree.
Outputs: out/taxonomy_beetles.csv, out/taxonomy_pathogens.csv, data/her_pathogen_gbif.parquet
"""
import json
import re
import time
from pathlib import Path

import duckdb
import pandas as pd
import requests
from rapidfuzz import fuzz, process
from unidecode import unidecode

A = Path(__file__).resolve().parents[2] / "data"
SNAP = Path(__file__).resolve().parents[2] / "data"
CACHE = A / "data/gbif_cache.jsonl"
S = requests.Session()
S.headers["User-Agent"] = "UF-Hulcr-lab pathogen dataset cross-check (research)"
KHINT = {"Fungi": "Fungi", "Microsporidia": "Fungi", "Bacteria": "Bacteria", "Nematode": "Animalia",
         "Protist": "Protozoa", "Virus": "Viruses"}


def binom(s):
    s = re.sub(r"\(.*?\)", " ", str(s or ""))
    w = re.findall(r"[A-Za-z][a-z-]+", s)
    return " ".join(w[:2]) if len(w) >= 2 else (w[0] if w else "")


def clean_pathogen(s):
    s = str(s or "").strip()
    s = re.sub(r"\b(cf|aff)\.\s*", "", s)
    s = re.sub(r"\s+(sp|spp)\.?(\s+\w+)?$", "", s)
    s = re.sub(r"\s+s\.\s*l\.?$", "", s)
    return re.sub(r"\s+", " ", s).strip()


cache = {}
if CACHE.exists():
    for line in open(CACHE):
        r = json.loads(line)
        cache[r["k"]] = r["v"]


def gbif(name, kingdom):
    k = f"{name}|{kingdom}"
    if k in cache:
        return cache[k]
    v = None
    for i in range(4):
        try:
            r = S.get("https://api.gbif.org/v1/species/match",
                      params={"name": name, "kingdom": kingdom, "strict": "false", "verbose": "false"}, timeout=60)
            if r.status_code == 200:
                v = r.json()
                break
        except Exception:
            pass
        time.sleep(3 * (i + 1))
    if v is not None:
        cache[k] = v
        with open(CACHE, "a") as fh:
            fh.write(json.dumps({"k": k, "v": v}) + "\n")
    time.sleep(0.1)
    return v or {}


def main():
    R = pd.read_parquet(A / "out/her_records.parquet")
    con = duckdb.connect(str(SNAP / "kg/kg.duckdb"), read_only=True)
    SP = con.sql("select species_id, name, authority, genus, tribe, subfamily from species").df()
    NM = con.sql("select species_id, valid_name, name, authority, year, is_valid_combination from names").df()
    NM["b"] = NM.name.map(binom)
    valid = dict(zip(SP.name.map(binom), SP.species_id))
    by_name = NM.groupby("b").species_id.agg(lambda s: sorted(set(s))).to_dict()
    sp_info = SP.set_index("species_id")
    allnames = sorted(set(NM.b) | set(valid))

    rows = []
    B = R.groupby(["ReportedBeetleSpecies", "AcceptedBeetleSpeciesName", "authorship", "genus", "tribe", "subfamily"],
                  dropna=False).agg(n_records=("pathogens", "size"),
                                    sources=("source", lambda s: "; ".join(sorted(set(s))[:5]))).reset_index()
    for b in B.itertuples():
        acc = binom(b.AcceptedBeetleSpeciesName)
        rep = binom(b.ReportedBeetleSpecies)
        row = dict(reported=b.ReportedBeetleSpecies, accepted_in_dataset=b.AcceptedBeetleSpeciesName,
                   authorship_in_dataset=b.authorship, genus_in_dataset=b.genus, tribe_in_dataset=b.tribe,
                   n_records=b.n_records, example_sources=b.sources)
        if acc in valid:
            sid = valid[acc]
            row.update(kg_status="valid name in the species list", kg_species_id=sid, kg_valid_name=sp_info.name[sid])
        elif acc in by_name:
            ids = by_name[acc]
            row.update(kg_status="synonym / earlier combination in the species list" if len(ids) == 1
                       else "ambiguous (name of several species)", kg_species_id=";".join(ids),
                       kg_valid_name=";".join(sp_info.name[i] for i in ids))
        else:
            near = process.extract(acc, allnames, scorer=fuzz.ratio, limit=3)
            ep = acc.split()[1] if len(acc.split()) == 2 else ""
            same_ep = NM[NM.b.str.split().str[-1].eq(ep) & (ep not in ("sp", "spp"))] if ep else NM.iloc[:0]
            # same epithet, possibly another genus (old combination) or a one-letter spelling variant
            if ep and len(ep) > 3:
                same_ep = pd.concat([same_ep, NM[NM.b.str.split().str[-1].map(
                    lambda e: fuzz.ratio(e, ep) >= 90 and e != ep)]])
            row.update(kg_status="not in the species list", kg_species_id="",
                       kg_valid_name="", closest_names="; ".join(f"{n} ({s:.0f})" for n, s, _ in near),
                       same_or_near_epithet_in_list="; ".join(sorted({f"{r.name} -> {r.valid_name}"
                                                                     for r in same_ep.itertuples()}))[:400])
        sid = str(row.get("kg_species_id", "")).split(";")[0]
        if sid in sp_info.index:
            row.update(kg_authority=sp_info.authority[sid], kg_genus=sp_info.genus[sid], kg_tribe=sp_info.tribe[sid],
                       kg_subfamily=sp_info.subfamily[sid])
            row["authority_differs"] = (unidecode(str(b.authorship)).replace(" ", "").strip("()") !=
                                        unidecode(str(sp_info.authority[sid])).replace(" ", "").strip("()"))
            row["tribe_differs"] = str(b.tribe) != str(sp_info.tribe[sid])
        # reported name -> same species?
        rep_ids = [valid[rep]] if rep in valid else by_name.get(rep, [])
        row["reported_resolves_to"] = ";".join(sp_info.name[i] for i in rep_ids if i in sp_info.index)
        row["reported_vs_accepted_consistent"] = (not rep_ids) or (sid in rep_ids)
        rows.append(row)
    TB = pd.DataFrame(rows)
    TB.to_csv(A / "out/taxonomy_beetles.csv", index=False)
    print("beetle name combinations", len(TB), "| accepted names", TB.accepted_in_dataset.nunique())
    print(TB.drop_duplicates("accepted_in_dataset").kg_status.value_counts().to_string())
    print(TB[TB.kg_status.ne("valid name in the species list")].to_string())
    print("authority differs:", int(TB.authority_differs.fillna(False).sum()), "| tribe differs:",
          int(TB.tribe_differs.fillna(False).sum()))
    print(TB[TB.tribe_differs.fillna(False) | TB.authority_differs.fillna(False)][
        ["accepted_in_dataset", "authorship_in_dataset", "kg_authority", "tribe_in_dataset", "kg_tribe"]].to_string())
    print("reported names resolving to another species than the accepted name:")
    print(TB[~TB.reported_vs_accepted_consistent][["reported", "accepted_in_dataset", "reported_resolves_to",
                                                   "kg_valid_name"]].to_string())

    # ---------------- pathogens
    PN = R.groupby(["pathogens", "categories"]).agg(n_records=("source", "size")).reset_index()
    KGP = pd.read_parquet(SNAP / "tables/pathogens.parquet",
                          columns=["partner_query_written", "partner", "gbif_name", "gbif_status", "gbif_species",
                                   "match_type", "source_kind"])
    kgmap = (KGP[KGP.gbif_name.notna()].assign(k=KGP.partner.str.lower().str.strip())
             .groupby("k").agg(kg_gbif_name=("gbif_name", "first"), kg_gbif_species=("gbif_species", "first"),
                               kg_gbif_status=("gbif_status", "first")).to_dict("index"))
    out = []
    for p in PN.itertuples():
        q = clean_pathogen(p.pathogens)
        single = len(q.split()) == 1
        g = gbif(q, KHINT.get(p.categories, "")) if q else {}
        acc = g.get("species") if g.get("rank") in ("SPECIES", "SUBSPECIES", "VARIETY", "FORM") else None
        flag = ""
        mt, st = g.get("matchType", "NONE"), g.get("status", "")
        if mt == "NONE":
            flag = "no GBIF match"
        elif mt == "HIGHERRANK":
            flag = "GBIF knows only a higher rank" if not single else ""
        elif st == "SYNONYM" and acc and binom(acc) != binom(q):
            flag = "GBIF: synonym of " + acc
        elif st == "DOUBTFUL":
            flag = "GBIF: doubtful name"
        elif mt == "FUZZY":
            flag = f"GBIF matched a different spelling: {g.get('canonicalName')}"
        if g.get("kingdom") and KHINT.get(p.categories) and g.get("kingdom") != KHINT[p.categories] and mt != "NONE":
            flag = (flag + "; " if flag else "") + f"GBIF kingdom {g.get('kingdom')} vs category {p.categories}"
        k = kgmap.get(str(p.pathogens).lower().strip(), {})
        out.append(dict(pathogen_in_dataset=p.pathogens, category=p.categories, n_records=p.n_records, queried=q,
                        gbif_match_type=mt, gbif_status=st, gbif_rank=g.get("rank"), gbif_canonical=g.get("canonicalName"),
                        gbif_accepted_species=acc, gbif_kingdom=g.get("kingdom"), gbif_phylum=g.get("phylum"),
                        gbif_class=g.get("class"), gbif_order=g.get("order"), gbif_family=g.get("family"),
                        gbif_genus=g.get("genus"), gbif_usage_key=g.get("usageKey"),
                        gbif_accepted_key=g.get("acceptedUsageKey") or g.get("usageKey"),
                        gbif_confidence=g.get("confidence"), flag=flag, kg_resolver_gbif_name=k.get("kg_gbif_name"),
                        kg_resolver_status=k.get("kg_gbif_status")))
    TP = pd.DataFrame(out)
    TP.to_parquet(A / "data/her_pathogen_gbif.parquet", index=False)
    TP.to_csv(A / "out/taxonomy_pathogens.csv", index=False)
    print("pathogen name x category rows", len(TP), "| distinct names", TP.pathogen_in_dataset.nunique(),
          "| distinct names case-insensitive", TP.pathogen_in_dataset.str.lower().str.strip().nunique())
    print("names in more than one category:", TP.groupby("pathogen_in_dataset").category.nunique().gt(1).sum())
    print(TP.groupby("pathogen_in_dataset").category.agg(list)[TP.groupby("pathogen_in_dataset").category.nunique() > 1])
    print("GBIF match types:", TP.gbif_match_type.value_counts().to_dict())
    print(TP[TP.flag.ne("")][["pathogen_in_dataset", "category", "n_records", "gbif_match_type", "gbif_status",
                              "gbif_canonical", "gbif_accepted_species", "flag"]].to_string())
    # near-duplicate names within her list (spelling variants of one taxon)
    names = sorted(TP.pathogen_in_dataset.unique())
    dup = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if a.lower() != b.lower() and fuzz.ratio(a.lower(), b.lower()) >= 90:
                dup.append((a, b, round(fuzz.ratio(a.lower(), b.lower()))))
            elif a.lower() == b.lower():
                dup.append((a, b, 100))
    print("near-identical pathogen names (possible spelling variants):")
    for d in dup:
        print("  ", d)
    pd.DataFrame(dup, columns=["name_a", "name_b", "similarity"]).to_csv(A / "out/taxonomy_pathogen_name_variants.csv",
                                                                         index=False)


if __name__ == "__main__":
    main()
