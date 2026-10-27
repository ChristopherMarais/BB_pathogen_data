"""Structured biodiversity and sequence database retrieval for bark and ambrosia beetle pathogens and parasites:
Part A: Global Biotic Interactions (GloBI)
Part B: NCBI GenBank (nuccore) host and isolation_source qualifiers
"""

# =====================================================================
# PART A: Global Biotic Interactions (GloBI)
# =====================================================================
"""Programmatic retrieval and taxonomic resolution of bark and ambrosia beetle pathogen and parasite
records from Global Biotic Interactions (GloBI).
"""
import argparse
import gzip
import hashlib
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import regex as re

KG = Path(__file__).resolve().parents[2] / "data"
sys.path.insert(0, str(KG / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "data"))
from kg_taxa import current_higher_of, higher_names  # noqa: E402

INVERSE = {"pathogenOf": "hasPathogen", "hasPathogen": "pathogenOf", "parasiteOf": "hasParasite",
           "hasParasite": "parasiteOf", "parasitoidOf": "hasParasitoid", "hasParasitoid": "parasitoidOf",
           "endoparasitoidOf": "hasParasitoid", "ectoparasitoidOf": "hasParasitoid",
           "ectoparasiteOf": "hasParasite", "endoparasiteOf": "hasParasite", "kleptoparasiteOf": "hasKleptoparasite",
           "preysOn": "preyedUponBy", "preyedUponBy": "preysOn", "eats": "eatenBy", "eatenBy": "eats",
           "kills": "killedBy", "killedBy": "kills", "vectorOf": "hasVector", "hasVector": "vectorOf",
           "dispersalVectorOf": "hasDispersalVector", "hasDispersalVector": "dispersalVectorOf", "hostOf": "hasHost",
           "hasHost": "hostOf", "livesOn": "hasEpibiont", "livesIn": "hasInhabitant", "livesInsideOf": "hasInhabitant",
           "symbiontOf": "symbiontOf", "mutualistOf": "mutualistOf", "commensalistOf": "commensalistOf",
           "interactsWith": "interactsWith", "visits": "visitedBy", "visitedBy": "visits", "laysEggsOn": "hasEggsLaidOnBy",
           "inhabits": "inhabitedBy", "ecologicallyRelatedTo": "ecologicallyRelatedTo"}
PLANT_HOST = {"eats", "hasHost", "livesIn", "livesOn", "livesInsideOf", "laysEggsOn", "kills", "parasiteOf",
              "pathogenOf", "inhabits", "endoparasiteOf", "ectoparasiteOf"}
PLANT_OTHER = {"interactsWith", "visits", "ecologicallyRelatedTo", "symbiontOf", "mutualistOf", "vectorOf"}
ENEMY = {"hasParasite": "parasite", "hasParasitoid": "parasitoid", "preyedUponBy": "predator", "eatenBy": "predator",
         "killedBy": "predator", "hasPathogen": "pathogen", "hasKleptoparasite": "parasite"}
ASSOC = {"hasEpibiont", "hasInhabitant", "commensalistOf", "symbiontOf", "mutualistOf", "interactsWith",
         "hasVector", "hasDispersalVector", "dispersalVectorOf", "vectorOf", "ecologicallyRelatedTo", "hostOf"}
SUBFAM = {"Scolytinae", "Platypodinae", "Scolytidae", "Platypodidae", "Ipidae"}
# datasets left out of the facts (kept in tables/globi_excluded.parquet with the reason)
EXCLUDE = {"fkeck/gpt_interactions": "interactions extracted by a third party's GPT model (Keck et al. 2025), not "
                                     "verified; the PMC articles it cites are read by our own pipeline instead",
           "globalbioticinteractions/trophiCH": "metaweb (trophiCH, Reji Chacko et al. 2023): links partly inferred "
                                                "from higher-taxon diet rules, not observations"}
BEETLE_PATH = re.compile(r"\b(?:Coleoptera|Curculionidae|Curculionoidea|Scolytinae|Platypodinae|Scolytidae|"
                         r"Platypodidae)\b")


def patterns(out):
    H = higher_names()
    pats = sorted({k[:1].upper() + k[1:] for k in H if len(k) >= 3} | SUBFAM)
    Path(out).write_text("\n".join(pats) + "\n")
    print("patterns", len(pats))


def kingdom(path):
    p = str(path or "")
    if re.search(r"\bMicrosporidi", p):
        return "Microsporidia"
    for k, rx in (("Fungi", r"\bFungi\b"), ("Plantae", r"\b(?:Plantae|Viridiplantae|Tracheophyta|Streptophyta)\b"),
                  ("Animalia", r"\b(?:Animalia|Metazoa)\b"), ("Bacteria", r"\bBacteria\b"),
                  ("Viruses", r"\b(?:Viruses|Riboviria|Duplodnaviria|Varidnaviria|Monodnaviria)\b"),
                  ("Oomycota", r"\b(?:Oomycota|Oomycetes|Peronosporomycetes)\b")):
        if re.search(rx, p):
            return k
    return "Protista" if re.search(r"\b(?:Apicomplexa|Amoebozoa|Ciliophora|Euglenozoa|Protozoa)\b", p) else ""


def is_mite(path):
    return bool(re.search(r"\b(?:Acari|Arachnida|Mesostigmata|Trombidiformes|Sarcoptiformes)\b", str(path or "")))


def make_resolver():
    from names_lib import Resolver
    R = Resolver()
    cache = {}

    def beetle(name):
        n = re.sub(r"\s+", " ", str(name or "")).strip()
        if n in cache:
            return cache[n]
        out = ("", None)
        w = n.split(" ")
        if len(w) >= 2 and re.match(r"^[A-Z][a-z]+$", w[0]) and re.match(r"^[a-z\-]+$", w[1]):
            r = R.resolve(w[0], w[1], subspecies=w[2] if len(w) > 2 and re.match(r"^[a-z\-]+$", w[2]) else None)
            if r["valid_species_id"] is not None:
                out = (f"species:{int(r['valid_species_id'])}", int(r["valid_species_id"]))
        elif len(w) == 1:
            h = current_higher_of(n)  # current names only (historical genera include homonyms)
            if h:
                out = (h, None)
        cache[n] = out
        return out
    return beetle


def map_row(t_b, kg, mite):
    """beetle-perspective GloBI verb + partner kingdom -> (type, relation) or None."""
    if kg == "Plantae":
        if t_b in PLANT_HOST:
            return "host", "host"
        if t_b in PLANT_OTHER:
            return "host", "other"
        return None
    if kg in ("Fungi", "Oomycota"):
        if t_b in ("hasPathogen", "hasParasite", "killedBy"):
            return "pathogen", "pathogen" if t_b != "hasParasite" else "parasite"
        if t_b in ASSOC or t_b in ("eats", "hasHost", "livesOn", "livesIn", "livesInsideOf", "inhabits"):
            return "fungus", "associated"
        return None
    if kg == "Animalia":
        if t_b in ENEMY:
            return "pathogen", ENEMY[t_b]
        if t_b in ("hasEpibiont", "commensalistOf", "symbiontOf", "hasDispersalVector", "dispersalVectorOf",
                   "hasInhabitant", "mutualistOf"):
            return "pathogen", "phoretic" if mite and t_b != "hasInhabitant" else "other"
        return None
    if kg in ("Bacteria", "Viruses", "Microsporidia", "Protista"):
        if t_b in ENEMY:
            return "pathogen", ENEMY[t_b]
        if t_b in ("symbiontOf", "mutualistOf", "hasInhabitant"):
            return "pathogen", "endosymbiont"
        if t_b in ASSOC:
            return "pathogen", "other"
    return None


def parse(path):
    hdr = (Path(path).parent / "header.tsv").read_text().rstrip("\n").split("\t")
    beetle = make_resolver()
    rows, seen, n = [], set(), 0
    for line in gzip.open(path, "rt", encoding="utf-8", errors="replace"):
        n += 1
        v = line.rstrip("\n").split("\t")
        if len(v) != len(hdr):
            continue
        d = dict(zip(hdr, v))
        sb, ss = beetle(d["sourceTaxonName"])
        tb, ts = beetle(d["targetTaxonName"])
        if sb and d["sourceTaxonPathNames"] and not BEETLE_PATH.search(d["sourceTaxonPathNames"]):
            sb = ""  # homonym outside the beetles (e.g. the mammal genus Platypus)
        if tb and d["targetTaxonPathNames"] and not BEETLE_PATH.search(d["targetTaxonPathNames"]):
            tb = ""
        if bool(sb) == bool(tb):  # neither side, or beetle-beetle
            continue
        it = d["interactionTypeName"]
        if sb:
            subj, vid, side, t_b = sb, ss, "target", it
        else:
            subj, vid, side, t_b = tb, ts, "source", INVERSE.get(it, it)
        pname, ppath = d[f"{side}TaxonName"], d[f"{side}TaxonPathNames"]
        kg = kingdom(ppath)
        m = map_row(t_b, kg, is_mite(ppath))
        if not m or not pname or pname.lower() in ("no:match", "none"):
            continue
        typ, rel = m
        if rel == "parasite" and re.search(r"\bHymenoptera\b", ppath):
            rel = "parasitoid"
        rdoi = d.get("referenceDoi", "").strip().lower().replace("https://doi.org/", "")
        rcit = d.get("referenceCitation", "").strip()
        ns = d.get("sourceNamespace", "")
        bside = "source" if sb else "target"
        eppo = ns.endswith("/eppo")
        if eppo:  # GloBI's EPPO import attaches references to host rows unreliably: cite the EPPO page itself
            code = next((x.split(":", 1)[1] for x in d[f"{bside}TaxonIds"].split("|") + [d[f"{bside}TaxonId"]]
                         if x.strip().startswith("EPPO:")), "")
            rdoi, rcit = "", f"EPPO Global Database, host plants of {d[f'{bside}TaxonName']}" + \
                (f" (EPPO code {code})" if code else "")
            sid = f"eppo_gd:{code or d[f'{bside}TaxonName']}"
        else:
            sid = f"doi:{rdoi}" if rdoi.startswith("10.") else \
                f"globi:{ns}:{hashlib.sha1(rcit.encode()).hexdigest()[:10]}" if rcit else f"globi:{ns}"
        key = (subj, typ, rel, pname, sid)
        if key in seen:
            continue
        seen.add(key)
        url = f"https://gd.eppo.int/taxon/{sid.split(':', 1)[1]}/hosts" if eppo and " " not in sid else \
            f"https://doi.org/{rdoi}" if rdoi.startswith("10.") else d.get("referenceUrl", "") or \
            f"https://www.globalbioticinteractions.org/?interactionType=interactsWith&sourceTaxon={d['sourceTaxonName']}"
        rows.append(dict(subject=subj, valid_species_id=vid, type=typ, relation=rel, setting=f"GloBI: {t_b}",
                         partner_group={"Plantae": "plant", "Fungi": "fungus", "Oomycota": "oomycete",
                                        "Animalia": "mite" if is_mite(ppath) else "animal", "Bacteria": "bacterium",
                                        "Viruses": "virus", "Microsporidia": "microsporidium",
                                        "Protista": "protist"}.get(kg, ""),
                         partner=pname, partner_full=pname, partner_query=pname,
                         hint={"Plantae": "Plantae", "Fungi": "Fungi", "Animalia": "Animalia",
                               "Bacteria": "Bacteria", "Viruses": "Viruses"}.get(kg, ""),
                         source_kind="dataset", source_id=sid,
                         citation=(rcit or "(no reference given)") + f"; via GloBI dataset: "
                                                                     f"{d.get('sourceCitation', '').strip()} ({ns})",
                         url=url, reference_doi=rdoi, reference_citation=rcit, globi_namespace=ns,
                         globi_interaction=it, beetle_written=d[f"{'source' if sb else 'target'}TaxonName"],
                         evidence=f"{d['sourceTaxonName']} | {it} | {d['targetTaxonName']} | {rcit[:300]}",
                         extraction=f"kg_globi:{ns}", retrieved=date.today().isoformat()))
    G = pd.DataFrame(rows)
    G["valid_species_id"] = G.valid_species_id.astype("Int64")
    ex = G.globi_namespace.isin(EXCLUDE)
    G[ex].assign(excluded_because=G.globi_namespace[ex].map(EXCLUDE)).to_parquet(KG / "tables/globi_excluded.parquet",
                                                                                 index=False)
    print("excluded rows", int(ex.sum()), G[ex].globi_namespace.value_counts().to_dict(), flush=True)
    G = G[~ex]
    G.to_parquet(KG / "tables/globi_interactions.parquet", index=False)
    sp = G[G.subject.str.startswith("species:")]
    print("rows read", n, "| kept", len(G), "| by type", G.type.value_counts().to_dict(), "| species by type",
          sp.groupby("type").valid_species_id.nunique().to_dict(), "| higher-taxon subjects",
          G.subject[~G.subject.str.startswith("species:")].nunique(), "| with reference DOI",
          int(G.reference_doi.str.startswith("10.").sum()), "| namespaces", G.globi_namespace.nunique(), flush=True)
    print(G.groupby(["type", "relation"]).size().to_string(), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--patterns")
    ap.add_argument("--parse")
    a = ap.parse_args()
    if a.patterns:
        patterns(a.patterns)
    if a.parse:
        parse(a.parse)


# =====================================================================
# PART B: NCBI GenBank (nuccore) Host / Isolation Qualifiers
# =====================================================================
"""Programmatic retrieval of Fungi, Nematoda, and Viruses sequence records isolated from bark and ambrosia
beetles in NCBI GenBank (nuccore) via Entrez E-utilities.
"""
import gzip
import json
import sys
from pathlib import Path

import pandas as pd
import regex as re

KG = Path(__file__).resolve().parents[2] / "data"
sys.path.insert(0, str(KG / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "data"))
from kg_ncbi import CACHE, eget, summaries  # noqa: E402
from kg_taxa import current_higher_of, higher_names  # noqa: E402

GROUPS = {"Fungi": "Fungi[Organism]", "Nematoda": "Nematoda[Organism]", "Viruses": "Viruses[Organism]"}
ENTOMOPATH = {"Beauveria", "Metarhizium", "Cordyceps", "Isaria", "Lecanicillium", "Akanthomyces", "Hirsutella",
              "Purpureocillium", "Conidiobolus", "Entomophthora", "Nosema", "Unikaryon", "Chytridiopsis", "Canningia",
              "Vairimorpha", "Endoreticulatus", "Larssoniella", "Microsporidium"}
WGS_CONTIG = re.compile(r"^[A-Z]{4,6}\d{8,}(?:\.\d+)?$")


def genus_words():
    H = higher_names()
    return sorted({k[:1].upper() + k[1:] for k, v in H.items() if v.startswith("genus:") and len(k) >= 3})


def search_all():
    gw = genus_words()
    chunks = [gw[i:i + 40] for i in range(0, len(gw), 40)]
    extra = '(Scolytinae[All Fields] OR Platypodinae[All Fields] OR "bark beetle"[All Fields] OR ' \
            '"ambrosia beetle"[All Fields] OR Scolytidae[All Fields] OR Platypodidae[All Fields])'
    recs = {}
    for g, gt in GROUPS.items():
        terms = [f"{gt} AND ({' OR '.join(w + '[All Fields]' for w in ch)})" for ch in chunks] + [f"{gt} AND {extra}"]
        for i, t in enumerate(terms):
            for r in summaries("nuccore", t, f"assoc_{g}_{i:02d}"):
                r["_group"] = g
                recs[r["uid"]] = r
    print("records found", len(recs), flush=True)
    return list(recs.values())


def make_finder():
    from names_lib import Resolver
    R = Resolver()
    gw = genus_words()
    rx = re.compile(r"\b(" + "|".join(sorted(gw, key=len, reverse=True)) + r")(?:\s+\(([A-Z][a-z]+)\))?"
                    r"(?:\s+([a-z][a-z\-]{2,}))?")
    cache = {}

    def find(text):
        out = set()
        for m in rx.finditer(text or ""):
            g, ep = m.group(1), m.group(3)
            if ep and ep not in ("sp", "spp", "and", "from", "of", "in", "on", "beetle", "beetles", "gallery"):
                k = (g, ep)
                if k not in cache:
                    r = R.resolve(g, ep)
                    cache[k] = f"species:{int(r['valid_species_id'])}" if r["valid_species_id"] is not None else ""
                if cache[k]:
                    out.add(cache[k])
                    continue
            h = current_higher_of(g)  # genus alone, "sp.", or an epithet not on the list; current names only
            if h:
                out.add(h)
        return out
    return find


def references(accs):
    """accession -> first non-direct-submission reference (title, journal, authors, doi, pmid); cached."""
    d = CACHE / "assoc_refs"
    d.mkdir(parents=True, exist_ok=True)
    out = {}
    for f in d.glob("*.jsonl.gz"):
        for line in gzip.open(f, "rt", encoding="utf-8"):
            x = json.loads(line)
            out[x["accession"]] = x
    todo = [a for a in accs if a not in out]
    print("references to fetch", len(todo), flush=True)
    for i in range(0, len(todo), 200):
        ch = todo[i:i + 200]
        x = eget("efetch.fcgi", db="nuccore", id=",".join(ch), rettype="gb", retmode="xml").text
        recs = []
        for seq in re.findall(r"<GBSeq>(.*?)</GBSeq>", x, re.S):
            acc = (re.search(r"<GBSeq_accession-version>(.*?)</", seq) or [None, ""])[1]
            ref = {}
            for rf in re.findall(r"<GBReference>(.*?)</GBReference>", seq, re.S):
                title = (re.search(r"<GBReference_title>(.*?)</", rf, re.S) or [None, ""])[1]
                if not title or title.strip().lower() == "direct submission":
                    continue
                doi = re.search(r"<GBXref_dbname>doi</GBXref_dbname>\s*<GBXref_id>(.*?)</", rf)
                au = re.findall(r"<GBAuthor>(.*?)</GBAuthor>", rf)
                ref = dict(title=title.strip(), journal=(re.search(r"<GBReference_journal>(.*?)</", rf, re.S) or
                                                         [None, ""])[1].strip(),
                           authors=au, doi=doi[1].lower() if doi else "",
                           pmid=(re.search(r"<GBReference_pubmed>(\d+)</", rf) or [None, ""])[1])
                break
            recs.append(dict(accession=acc, ref=ref))
        with gzip.open(d / f"b{len(list(d.glob('*.jsonl.gz'))):05d}.jsonl.gz", "wt", encoding="utf-8") as fh:
            for r in recs:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
                out[r["accession"]] = r
        if i % 4000 == 0:
            print(" ", i, "of", len(todo), flush=True)
    return out


def main():
    recs = search_all()
    find = make_finder()
    rows = []
    for r in recs:
        acc = r.get("accessionversion") or r.get("caption", "")
        if WGS_CONTIG.match(acc) and not re.search(r"0{6,}(?:\.\d+)?$", acc):
            continue
        st, sn = (r.get("subtype") or "").split("|"), (r.get("subname") or "").split("|")
        sub = dict(zip(st, sn)) if len(st) == len(sn) else {}
        host, iso = sub.get("host", ""), sub.get("isolation_source", "")
        subj = find(host) | find(iso) | find(r.get("title", ""))
        if any(x.startswith("species:") for x in subj):
            subj = {x for x in subj if x.startswith("species:")}  # a genus is the subject only when no species is named
        if not subj:
            continue
        org = r.get("organism", "")
        g = r["_group"]
        txt = f"{host} {iso}".lower()
        rel = "mycangial_symbiont" if "mycang" in txt else "gut_associate" if re.search(r"\bgut|intestin", txt) \
            else "isolated_from"
        typ = "pathogen" if g != "Fungi" or org.split(" ")[0] in ENTOMOPATH else "fungus"
        for s in subj:
            rows.append(dict(subject=s, accession=acc, organism=org, group=g, type=typ, relation=rel, host=host,
                             isolation_source=iso, country=sub.get("country", "") or sub.get("geo_loc_name", ""),
                             title=r.get("title", ""), createdate=r.get("createdate", ""),
                             retrieved=r.get("_retrieved", "")))
    A = pd.DataFrame(rows)
    print("records naming our beetles", A.accession.nunique(), "| rows", len(A), flush=True)
    refs = references(sorted(A.accession.unique()))
    ref = A.accession.map(lambda a: (refs.get(a) or {}).get("ref") or {})
    A["ref_doi"] = ref.map(lambda x: x.get("doi", ""))
    A["ref_pmid"] = ref.map(lambda x: x.get("pmid", ""))
    A["ref_cit"] = ref.map(lambda x: (f"{x['authors'][0].split(',')[0]} et al. " if len(x.get("authors", [])) > 2 else
                                      " & ".join(a.split(",")[0] for a in x.get("authors", [])) + " " if x else "") +
                           f"{x.get('title', '')}. {x.get('journal', '')}".strip() if x else "")
    A["source_id"] = ["doi:" + d if d else f"pmid:{p}" if p else f"ncbi:nuccore:{a}"
                      for d, p, a in zip(A.ref_doi, A.ref_pmid, A.accession)]
    G = A.groupby(["subject", "organism", "type", "relation", "source_id"]).agg(
        group=("group", "first"), n_records=("accession", "nunique"),
        accessions=("accession", lambda s: ",".join(sorted(set(s))[:10])),
        hosts_written=("host", lambda s: " | ".join(sorted({x for x in s if x})[:5])),
        isolation_sources=("isolation_source", lambda s: " | ".join(sorted({x for x in s if x})[:5])),
        countries=("country", lambda s: " | ".join(sorted({x.split(":")[0] for x in s if x})[:10])),
        ref_cit=("ref_cit", "first"), ref_doi=("ref_doi", "first"), retrieved=("retrieved", "first")).reset_index()
    G["valid_species_id"] = pd.to_numeric(G.subject.str.extract(r"^species:(\d+)$")[0]).astype("Int64")
    G["partner"] = G.organism
    G["partner_full"] = G.organism
    G["partner_query"] = G.organism.str.replace(r"\s+(?:sp\.|cf\.|aff\.).*$", "", regex=True)
    G["partner_group"] = G.group.map({"Fungi": "fungus", "Nematoda": "nematode", "Viruses": "virus"})
    G["hint"] = G.group.map({"Fungi": "Fungi", "Nematoda": "Animalia", "Viruses": "Viruses"})
    G["setting"] = "GenBank isolate"
    G["source_kind"] = "database"
    G["citation"] = [(c + "; " if c else "") + f"GenBank {a.split(',')[0]} ({o}; host / source: {h or i})"
                     for c, a, o, h, i in zip(G.ref_cit, G.accessions, G.organism, G.hosts_written,
                                              G.isolation_sources)]
    G["url"] = ["https://doi.org/" + d if d else f"https://www.ncbi.nlm.nih.gov/nuccore/{a.split(',')[0]}"
                for d, a in zip(G.ref_doi, G.accessions)]
    G["evidence"] = [f"{o} | host: {h} | isolation source: {i} | {a}" for o, h, i, a in
                     zip(G.organism, G.hosts_written, G.isolation_sources, G.accessions)]
    G["extraction"] = "kg_ncbi_assoc"
    G.to_parquet(KG / "tables/ncbi_assoc.parquet", index=False)
    sp = G[G.valid_species_id.notna()]
    print("rows", len(G), "| by type", G.type.value_counts().to_dict(), "| species by type",
          sp.groupby("type").valid_species_id.nunique().to_dict(), "| higher-taxon subjects",
          G.subject[G.valid_species_id.isna()].nunique(), "| with a paper (doi / pmid)",
          int((~G.source_id.str.startswith("ncbi:")).sum()),
          "| partner organisms", G.organism.nunique(), flush=True)


if __name__ == "__main__":
    main()
