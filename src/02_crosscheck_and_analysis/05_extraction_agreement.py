"""Step 5 (analysis 2): agreement of the KG's automated reading with her expert extraction, for her publications the
KG also read. Her records are the gold standard.

KG sources of a publication: 'doi:<doi>'; sources whose title matches (source_meta; year +-1, title similarity >= 92);
and the KG ids of the same reference of her list (abs:pathref:<n>, pref_ copies). The reading is tagged
'independent' when the text came from the KG's own stores, 'her list' when it was fetched because it was on her list.
Two KG layers are compared:
  reader  = every pathogen / fungus claim the reading model (Gemma 4 26B) made for the source, before the KG's rule checks
            (claims_checked; a claim is beetle + organism + relation + verbatim quote);
  graph   = the edges that entered the automated literature and database pipeline from those claims after the rule checks (method llm_facts+checks;
            predicates has_natural_enemy, has_associated_organism, has_fungal_associate, not_associated_with).
Pairs are unique (beetle species, organism) per publication. Beetle: the species-list id of her accepted name.
Organism match levels: exact = same binomial (written name or GBIF name, 'cf.' / 'sp.' removed) or same GBIF taxon;
genus = same genus; category = KG names an organism of the same category (fungi, bacteria, viruses, nematodes, protists,
microsporidia) for that beetle in that paper (credits vague names such as 'hemolymph nematodes').
Recall = her pairs found by the KG / her pairs in the publications the KG read.
Precision = KG pairs in her six categories that match one of her pairs (genus level) / all such KG pairs. KG pairs she
did not record are not necessarily errors (she excluded associates and organisms not called pathogens); precision is
therefore a lower bound for 'correct reading', and an estimate of how far the automated filter matches her rules.
Outputs: out/agreement_by_publication.csv, out/agreement_pairs.csv, out/agreement_summary.csv
"""
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from rapidfuzz import fuzz
from unidecode import unidecode

A = Path(__file__).resolve().parents[2] / "data"
SNAP = Path(__file__).resolve().parents[2] / "data"
ST = Path(__file__).resolve().parents[2] / "data"
CATS = ["Fungi", "Bacteria", "Virus", "Nematode", "Protist", "Microsporidia"]
ORG_PRED = ("has_natural_enemy", "has_associated_organism", "has_fungal_associate", "not_associated_with")


def norm(s):
    s = unidecode(str(s or "")).lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def clean_name(s):
    s = unidecode(str(s or "")).strip()
    s = re.sub(r"\b(cf|aff)\.\s*", "", s)
    s = re.sub(r"\s+(sp|spp|sp\d+)\.?(\s+\w+)?$", "", s)
    w = re.findall(r"[A-Za-z][a-z-]+", s)
    if len(w) >= 2 and w[0][0].isupper() and w[1][0].islower():
        return f"{w[0]} {w[1]}".lower()
    return s.lower()


def genus_of(s):
    w = re.findall(r"[A-Za-z][a-z-]+", unidecode(str(s or "")))
    return w[0].lower() if w and w[0][0].isupper() else ""


GENERIC = {"unidentified", "unknown", "nematode", "nematodes", "rod", "baculovirus", "entomopoxvirus", "virus",
           "bacteria", "fungus", "fungi", "microsporidia", "microsporidium", "gut", "hemolymph", "intestinal",
           "nematodan", "un", "similar", "tylenchid", "flacherie"}


def cat_of(kingdom, phylum, group, name):
    k, p, g, n = (str(x or "").lower() for x in (kingdom, phylum, group, name))
    if "virus" in n or k == "viruses" or "epv" in n or g.startswith("virus"):
        return "Virus"
    if p == "microsporidia" or g.startswith("microsporid") or re.search(r"microsporid|nosema|unikaryon|chytridiopsis|"
                                                                       r"canningia|larssoniella|pleistophora|stempellia", n):
        return "Microsporidia"
    if k == "bacteria" or g.startswith("bacteri"):
        return "Bacteria"
    if p == "nematoda" or g.startswith("nematod") or "nematod" in n:
        return "Nematode"
    if p == "oomycota" or g.startswith("oomycet"):
        return "Other"
    if k in ("protozoa", "chromista") or g.startswith("protist") or re.search(r"gregarin|mattesia|menzbieria|"
                                                                             r"malamoeba|amoeb|protozo|ophryocystis", n):
        return "Protist"
    if k == "fungi" or g.startswith(("fung", "insect-pathogenic", "yeast", "ophiostom", "microascal")) or \
            re.search(r"beauveria|metarhizium|isaria|fung", n):
        return "Fungi"
    if k in ("animalia", "plantae") or g.startswith(("mite", "insect", "parasitoid", "predator", "bird", "fl", "ant",
                                                     "spider", "true bug", "other beetle", "snakefl", "mammal")):
        return "Other"
    return "Unknown"


def binom_label(lbl):
    w = re.findall(r"[A-Za-z][a-z-]+", unidecode(str(lbl or "")))
    if len(w) >= 2 and w[0][0].isupper() and w[1][0].islower():
        return f"{w[0]} {w[1]}".lower()
    return w[0].lower() if w else ""


def main():
    P = pd.read_parquet(A / "out/her_publications.parquet")
    R = pd.read_parquet(A / "out/her_records.parquet")
    RC = pd.read_csv(A / "out/recall_by_publication.csv")[["pub_id", "found_kg_search"]]
    TB = pd.read_csv(A / "out/taxonomy_beetles.csv")
    TP = pd.read_parquet(A / "data/her_pathogen_gbif.parquet")
    sid_of = {r.accepted_in_dataset: str(r.kg_species_id).split(";")[0] for r in TB.itertuples()
              if isinstance(r.kg_species_id, str) and r.kg_species_id}
    # synonyms resolved by epithet (Ips oregonis -> Ips pini) are left unmatched here on purpose (her accepted name)
    R = R.merge(P[["ref", "pub_id"]], on="ref")
    R["sid"] = R.AcceptedBeetleSpeciesName.map(sid_of)
    gb = TP.set_index("pathogen_in_dataset")
    R["p_clean"] = R.pathogens.map(clean_name)
    R["p_genus"] = R.pathogens.map(genus_of)
    R["p_gbif_species"] = R.pathogens.map(lambda n: str(gb.gbif_accepted_species.get(n) or "").lower())
    R["p_gbif_genus"] = R.pathogens.map(lambda n: str(gb.gbif_genus.get(n) or "").lower())
    R["p_gbif_key"] = R.pathogens.map(lambda n: gb.gbif_accepted_key.get(n) if gb.gbif_rank.get(n) in
                                      ("SPECIES", "SUBSPECIES", "VARIETY") else None)
    R["p_usage_key"] = R.pathogens.map(lambda n: gb.gbif_usage_key.get(n) if gb.gbif_rank.get(n) in
                                       ("SPECIES", "SUBSPECIES", "VARIETY") else None)
    R["vague"] = R.p_genus.isin(GENERIC) | R.p_genus.eq("")
    HP = R.drop_duplicates(["pub_id", "AcceptedBeetleSpeciesName", "pathogens"]).copy()
    print("her unique (publication, beetle, organism) pairs:", len(HP), "| with beetle id", int(HP.sid.notna().sum()))

    con = duckdb.connect(str(SNAP / "kg/kg.duckdb"), read_only=True)
    SM = con.sql("select source_id, title, year, doi from source_meta").df()
    SM["tn"] = SM.title.map(norm)
    SM["year"] = pd.to_numeric(SM.year, errors="coerce")
    PR = pd.read_parquet(SNAP / "tables/pathogen_refs.parquet")
    PRF = pd.read_parquet(SNAP / "tables/pathogen_refs_fetch.parquet")
    SP_ = pd.read_csv(ST / "work/sources_pathrefs.csv")
    SA = pd.read_csv(ST / "work/sources_abstracts.csv", usecols=["source_id", "doi", "pathogen_ref"])
    ab_pref = set(SA.loc[SA.pathogen_ref.astype(str).eq("True"), "source_id"])

    ids_of = {}
    for p in P.itertuples():
        ids = set()
        if p.doi:
            ids.add("doi:" + p.doi)
        t = norm(p.title_parsed)
        if len(t) >= 25:
            cand = SM[(SM.year - p.year).abs().le(1) | SM.year.isna()]
            for r in cand.itertuples():
                if len(r.tn) >= 20 and (fuzz.ratio(r.tn, t) >= 92 or
                                        (len(r.tn) >= len(t) and fuzz.partial_ratio(t, r.tn) >= 95)):
                    ids.add(r.source_id)
        for r in PR[(PR.year - p.year).abs().le(1)].itertuples():
            if norm(r.first_author) and norm(r.first_author) in norm(p.authors) and \
                    fuzz.token_set_ratio(norm(r.rest), norm(p.rest)) >= 90:
                ids.add(f"abs:pathref:{int(r.n)}")
                if r.held_source_id:
                    ids.add(r.held_source_id)
                for q in PRF[PRF.n.eq(r.n) & PRF.source_id.fillna("").ne("")].source_id:
                    ids.add(q)
        ids_of[p.pub_id] = ids
    allids = set().union(*ids_of.values())
    (A / "data/her_pub_kg_source_ids.json").write_text(json.dumps({k: sorted(v) for k, v in ids_of.items()}, indent=0))
    print("publications with >= 1 KG source id candidate:", sum(1 for v in ids_of.values() if v))

    CC = pd.read_parquet(SNAP / "tables/claims_checked.parquet",
                         columns=["source_id", "doc_id", "type", "valid_species_id", "taxon_written", "partner",
                                  "partner_full", "partner_group", "relation", "quote", "passes_checks", "reason",
                                  "claim_id", "subject"])
    CC = CC[CC.source_id.isin(allids)].copy()
    E = con.sql(f"""select e.edge_id, e.subject, e.predicate, e.object, e.source_id, e.claim_id, e.props, e.evidence,
                    n.label, n.props nprops, n.node_type from edges e left join nodes n on n.node_id = e.object
                    where e.method = 'llm_facts+checks' and e.predicate in {ORG_PRED}""").df()
    E = E[E.source_id.isin(allids)].copy()
    print("claims rows for these sources", len(CC), "| graph edges", len(E))

    def nprop(s, k):
        try:
            return json.loads(s or "{}").get(k)
        except Exception:
            return None
    E["kingdom"] = E.nprops.map(lambda s: nprop(s, "gbif_kingdom"))
    E["phylum"] = E.nprops.map(lambda s: nprop(s, "gbif_phylum"))
    E["group"] = E.nprops.map(lambda s: nprop(s, "group"))
    E["gbif_genus"] = E.nprops.map(lambda s: str(nprop(s, "gbif_genus") or "").lower())
    E["written"] = E.props.map(lambda s: nprop(s, "partner"))
    E["relation"] = E.props.map(lambda s: nprop(s, "relation"))
    E["cat"] = [cat_of(k, ph, g, f"{w} {l}") for k, ph, g, w, l in zip(E.kingdom, E.phylum, E.group, E.written, E.label)]
    E["b_label"] = E.label.map(binom_label)
    E["b_written"] = E.written.map(clean_name)
    E["g_written"] = E.written.map(genus_of)
    E["gkey"] = E.object.str.extract(r"^gbif:(\d+)")[0].astype(float)
    CC["cat"] = [cat_of("", "", g, f"{p} {pf}") for g, p, pf in zip(CC.partner_group, CC.partner, CC.partner_full)]
    CC["b_written"] = CC.partner.map(clean_name)
    CC["g_written"] = CC.partner.map(genus_of)
    CC["sid"] = "species:" + CC.valid_species_id.astype("Int64").astype(str)

    rows, pair_rows = [], []
    for p in P.itertuples():
        ids = ids_of[p.pub_id]
        c_all = CC[CC.source_id.isin(ids)]
        c_org = c_all[c_all.type.isin(["pathogen", "fungus"])]
        e = E[E.source_id.isin(ids)]
        docs = set(c_all.doc_id.dropna())
        read = len(c_all) > 0
        text = ("full text" if any(not d.startswith("abs_") for d in docs) else "abstract only") if read else ""
        via_her = any(d.startswith("pref_") for d in docs) or any(i.startswith("abs:pathref") and
                                                                  (c_all.source_id == i).any() for i in ids)
        if text == "abstract only" and not via_her:
            via_her = any(s in ab_pref for s in c_all.source_id) and not bool(RC.set_index("pub_id").found_kg_search[p.pub_id])
        hp = HP[HP.pub_id.eq(p.pub_id)]
        rows.append(dict(pub_id=p.pub_id, source=p.source, year=p.year, categories=p.categories, n_her_pairs=len(hp),
                         kg_source_ids=";".join(sorted(ids)), kg_read=read, kg_text=text,
                         kg_read_via="her list" if via_her else ("independent" if read else ""),
                         n_claims_all=len(c_all), n_claims_organism=len(c_org), n_graph_edges=len(e)))
        if not read:
            continue
        for h in hp.itertuples():
            def match(df, layer):
                if h.sid is None or (isinstance(h.sid, float)):
                    return "beetle not on list", ""
                d = df[df.sid.eq(h.sid)] if layer == "reader" else df[df.subject.eq(h.sid)]
                if layer == "reader":
                    ex = d[d.b_written.eq(h.p_clean) | (d.b_written.eq(h.p_gbif_species) & (h.p_gbif_species != ""))]
                else:
                    ex = d[d.b_written.eq(h.p_clean) | d.b_label.eq(h.p_clean) |
                           (d.b_label.eq(h.p_gbif_species) & (h.p_gbif_species != "")) |
                           (d.gkey.notna() & (d.gkey.isin([h.p_gbif_key, h.p_usage_key])))]
                if len(ex):
                    return "exact", ex.iloc[0].get("reason", "") if layer == "reader" else ""
                gg = {h.p_genus, h.p_gbif_genus} - {""}
                if not h.vague and gg:
                    ge = d[d.g_written.isin(gg)] if layer == "reader" else d[d.g_written.isin(gg) | d.gbif_genus.isin(gg)]
                    if len(ge):
                        return "genus", ""
                if (d.cat == h.categories).any():
                    return "category", ""
                return "none", ""
            m_reader, _ = match(c_org, "reader")
            m_graph, _ = match(e, "graph")
            # reader claim for the pair rejected by a rule?
            rej = c_org[(c_org.sid == h.sid) & (c_org.g_written.isin({h.p_genus, h.p_gbif_genus} - {""}))
                        & ~c_org.passes_checks.astype(bool)]
            pair_rows.append(dict(pub_id=p.pub_id, source=p.source, kg_text=text,
                                  kg_read_via="her list" if via_her else "independent", beetle=h.AcceptedBeetleSpeciesName,
                                  pathogen=h.pathogens, category=h.categories, vague_name=h.vague, reader=m_reader,
                                  graph=m_graph, rejected_reason="; ".join(sorted(set(rej.reason.astype(str)))),
                                  graph_predicates="; ".join(sorted(set(e[e.subject.eq(h.sid) & (
                                      e.g_written.isin({h.p_genus, h.p_gbif_genus} - {""}) |
                                      e.gbif_genus.isin({h.p_genus, h.p_gbif_genus} - {""}))].predicate)))))
        # precision: KG pairs (graph, species-level beetle, her six categories) in this paper
        # the KG's pathogen / parasite statements: has_natural_enemy, relation pathogen or parasite, her six categories
        ein = e[e.cat.isin(CATS) & e.subject.str.startswith("species:") & e.predicate.eq("has_natural_enemy")
                & e.relation.isin(["pathogen", "parasite"])].drop_duplicates(["subject", "object"])
        for x in ein.itertuples():
            hh = hp[hp.sid.eq(x.subject)]
            gg = {x.g_written, x.gbif_genus, genus_of(x.label)} - {""}
            kg_generic = x.object.startswith("unmatched:") and bool(re.search(
                r"fung|isolat|microsporid|gregarin|nematod|virus|bacteri|protozo|entomopath|pathogen", str(x.written), re.I))
            org_any = hp[(hp.p_genus.isin(gg) | hp.p_gbif_genus.isin(gg))] if gg and not kg_generic else hp.iloc[:0]
            if len(hh) and (hh.p_genus.isin(gg) | hh.p_gbif_genus.isin(gg)).any() and not kg_generic:
                st = "match: same beetle, same organism genus"
            elif len(hh) and (hh.categories == x.cat).any() and (kg_generic or hh[hh.categories == x.cat].vague.any()):
                st = "match: same beetle, same category, name vague"
            elif len(org_any):
                st = "organism in her records of this paper, other beetle"
            elif len(hh):
                st = "beetle in her records of this paper, organism not"
            else:
                st = "neither beetle nor organism in her records of this paper"
            pair_rows.append(dict(pub_id=p.pub_id, source=p.source, kg_text=text,
                                  kg_read_via="her list" if via_her else "independent", side="kg",
                                  beetle=x.subject, pathogen=f"{x.written} -> {x.label}", category=x.cat,
                                  predicate=x.predicate, relation=x.relation, precision_status=st))
    B = pd.DataFrame(rows)
    PP = pd.DataFrame(pair_rows)
    PP["side"] = PP.get("side", pd.Series(dtype=str)).fillna("her")
    B.to_csv(A / "out/agreement_by_publication.csv", index=False)
    PP.to_csv(A / "out/agreement_pairs.csv", index=False)
    print("publications read by the KG:", int(B.kg_read.sum()), "of", len(B), "|", B[B.kg_read].kg_text.value_counts().to_dict(),
          "| via", B[B.kg_read].kg_read_via.value_counts().to_dict())
    H = PP[PP.side.eq("her")]
    K = PP[PP.side.eq("kg")]

    def wilson(k, n, z=1.96):
        if n == 0:
            return (np.nan, np.nan)
        p = k / n
        d = 1 + z * z / n
        c = (p + z * z / (2 * n)) / d
        h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
        return f"{100 * (c - h):.1f}-{100 * (c + h):.1f}"
    out = []
    for subset, hm, km in [("all", H.index == H.index, K.index == K.index),
                           ("full text", H.kg_text.eq("full text"), K.kg_text.eq("full text")),
                           ("abstract only", H.kg_text.eq("abstract only"), K.kg_text.eq("abstract only")),
                           ("read independently", H.kg_read_via.eq("independent"), K.kg_read_via.eq("independent")),
                           ("read via her list", H.kg_read_via.eq("her list"), K.kg_read_via.eq("her list"))]:
        for cat in ["all"] + CATS:
            h = H[hm & (H.category.eq(cat) if cat != "all" else True)]
            k = K[km & (K.category.eq(cat) if cat != "all" else True)]
            row = dict(subset=subset, category=cat, her_pairs=len(h))
            for layer in ["reader", "graph"]:
                for lvl, ok in [("exact", ["exact"]), ("genus", ["exact", "genus"]),
                                ("category", ["exact", "genus", "category"])]:
                    n_ok = int(h[layer].isin(ok).sum())
                    row[f"recall_{layer}_{lvl}"] = round(100 * n_ok / len(h), 1) if len(h) else np.nan
                    if lvl == "genus":
                        row[f"recall_{layer}_{lvl}_ci"] = wilson(n_ok, len(h))
            row["kg_pairs_in_scope"] = len(k)
            for nm, ok in [("strict", ["match: same beetle, same organism genus"]),
                           ("with_vague", ["match: same beetle, same organism genus",
                                           "match: same beetle, same category, name vague"]),
                           ("organism_any_beetle", ["match: same beetle, same organism genus",
                                                    "match: same beetle, same category, name vague",
                                                    "organism in her records of this paper, other beetle"])]:
                okp = k.precision_status.isin(ok)
                row[f"precision_{nm}"] = round(100 * okp.mean(), 1) if len(k) else np.nan
                row[f"precision_{nm}_ci"] = wilson(int(okp.sum()), len(k))
            for s in ["beetle in her records of this paper, organism not",
                      "neither beetle nor organism in her records of this paper"]:
                row["n " + s] = int(k.precision_status.eq(s).sum())
            out.append(row)
    O = pd.DataFrame(out)
    O.to_csv(A / "out/agreement_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(O.to_string())
    print("her pairs, graph match:", H.graph.value_counts().to_dict(), "| reader:", H.reader.value_counts().to_dict())
    print("rejected reasons for missed pairs:", H[H.graph.isin(["none", "category"])].rejected_reason.value_counts().head(15).to_dict())
    print("KG predicate for genus-level matches:", H[H.graph.isin(["exact", "genus"])].graph_predicates.value_counts().head(10).to_dict())
    print("KG in-scope pairs precision status:", K.precision_status.value_counts().to_dict())
    print(K[K.precision_status.str.startswith(("beetle in", "neither"))][["source", "beetle", "pathogen", "category",
                                                                       "predicate", "relation"]].head(60).to_string())


if __name__ == "__main__":
    main()
