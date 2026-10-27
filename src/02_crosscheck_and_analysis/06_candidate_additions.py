"""Step 6 (analyses 3 and 5): KG pathogen / parasite records resting on sources NOT in her dataset, filtered with HER
inclusion rules -> candidates for the authors to verify (not additions).

Starting set: automated-pipeline edges read from the literature by the automated pipeline (method llm_facts+checks), plus,
as a separate stratum, rows of other curated databases (GloBI datasets with their original reference; GenBank isolate
records). Excluded as derived from her list: the Notion copy of her database (published_dataset rows 'via the user's
pathogen database', notion_pathogen_db:*), texts fetched from her reference list (pref_ documents, abs:pathref:*).
Her inclusion rules, applied in this order:
  R1 relation: predicate has_natural_enemy with relation pathogen or parasite (predators, parasitoids excluded;
     has_associated_organism = no harm stated -> excluded: phoretic, gut associates, endosymbionts, 'isolated from');
  R2 category: the organism is a fungus, bacterium, virus, nematode, protist or microsporidian (GBIF classification of
     the resolved name, else name words); mites, insects, birds and other animals excluded;
  R3 fungi and bacteria only when called pathogens of the BEETLE: an insect-pathogenic fungus group (Cordycipitaceae,
     Clavicipitaceae, Ophiocordycipitaceae, Entomophthorales, Laboulbeniales ...) or a quote with harm wording
     (entomopathogenic, pathogenic to, mycosis, infected, mortality, killed, virulence ...); quotes about tree disease
     or symbiosis (wilt, canker, blue stain, vector, mycangium, ambrosia fungus, symbiont) without harm wording are
     excluded; ophiostomatoid fungi, Microascales, wood-decay / plant-pathogen orders are excluded;
  R4 primary sources only: reviews, catalogues, bibliographies, books of synthesis, websites and Wikipedia excluded
     (title words), and second-hand statements (the quote cites another work; KG field second_hand_of) are moved to a
     'pointer' list (the cited work may be a primary source she lacks);
  R5 the source is not one of her 280 publications (DOI, title match, or KG id of her reference; step 5 mapping).
Outputs (out/): candidates_summary.csv, candidates_by_category_decade.csv,
  INTERNAL_candidate_review_table.csv (with supporting quotes: lab only, NOT for the paper),
  INTERNAL_candidate_sample50.csv (random sample for manual precision check), candidates_pointers.csv,
  platypodinae_records.csv
"""
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from unidecode import unidecode

import importlib.util

A = Path(__file__).resolve().parents[2] / "data"
SNAP = Path(__file__).resolve().parents[2] / "data"
ST = Path(__file__).resolve().parents[2] / "data"
spec = importlib.util.spec_from_file_location("agr", A / "scripts/05_extraction_agreement.py")
agr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agr)
cat_of, clean_name, genus_of, binom_label = agr.cat_of, agr.clean_name, agr.genus_of, agr.binom_label
CATS = agr.CATS

HARM = re.compile(r"entomopatho|entomogenous|insect[- ]pathogen|pathogenic (to|for|against|on)|pathogen of|"
                  r"pathogens? (of|in|from|to) (the )?(bark|ambrosia|beetle|adult|larva|scolyt|ips|dendroctonus)|"
                  r"mycos[ie]s|infect|mortality|kill|died|dead (beetle|adult|larva)|cadaver|virulen|lethal|"
                  r"larvicid|insecticid|biocontrol|biological control|parasit|disease of (the )?(beetle|insect)|"
                  r"diseased (beetle|adult|larva|insect)|mummif|sporulat|LT50|LC50|LD50|bioassay", re.I)
TREE = re.compile(r"wilt|canker|blue[- ]?stain|sap[- ]?stain|root rot|dieback|die-back|decline of|vector|vectored|"
                  r"mycangi|symbio|ambrosia fung|tree (death|mortality)|phytopatho|plant pathogen|"
                  r"pathogenic to (the )?(tree|host tree|seedling|pine|spruce|oak|elm|avocado)|lesion|inoculat\w+ (into|in) (the )?(tree|seedling|log|stem)", re.I)
BAD_FUNGI = re.compile(r"ophiostom|microascal|ceratocyst|grosmann|leptograph|raffaelea|harringtonia|ambrosiella|"
                       r"dryadomyces|geosmithia|fusarium|neocosmospora|armillaria|heterobasidion|diplodia|botryosph|"
                       r"phytophthora|ophiostoma|endoconidiophora|sporothrix|graphium|entomocorticium|flavodon|"
                       r"polyporal|russulal|agarical|hymenochaet|diaporthal|helotial|pucciniales|mycosphaerel|"
                       r"corticial|eurotial", re.I)
INSECT_FUNGI = re.compile(r"insect-pathogenic|cordycip|clavicip|ophiocordycip|entomophthor|laboulben|"
                          r"beauveria|metarhizium|isaria|cordyceps|lecanicillium|akanthomyces|paecilomyces|"
                          r"verticillium lecanii|hirsutella|tolypocladium|entomophaga|zoophthora|metschnikow|"
                          r"conidiobolus|pandora|aschersonia|nomuraea|simplicillium|purpureocillium", re.I)
REVIEW = re.compile(r"\breview\b|overview|synthesis|catalog|catalogue|bibliograph|checklist|check-list|"
                    r"natural enemies of|pathogens in bark beetles|bark and wood boring insects in living trees|"
                    r"handbook|encyclopa?edia|textbook|world list|annotated list|state of knowledge|"
                    r"literature|a list of|insect pathogens? database", re.I)


def main():
    con = duckdb.connect(str(SNAP / "kg/kg.duckdb"), read_only=True)
    E = con.sql("""select e.edge_id, e.subject, e.predicate, e.object, e.source_id, e.citation, e.url, e.evidence,
                   e.locator, e.claim_id, e.method, e.props, n.label, n.props nprops,
                   ef.fact_id from edges e left join nodes n on n.node_id = e.object
                   left join edge_facts ef on ef.edge_id = e.edge_id
                   where e.predicate in ('has_natural_enemy','has_associated_organism')
                   and e.method in ('llm_facts+checks','published_dataset','genbank_isolate')""").df()
    print("organism edges (literature + databases)", len(E), E.method.value_counts().to_dict())

    def jp(s, k):
        try:
            return json.loads(s or "{}").get(k)
        except Exception:
            return None
    for k, c in [("gbif_kingdom", "kingdom"), ("gbif_phylum", "phylum"), ("gbif_class", "class"),
                 ("gbif_order", "order"), ("gbif_family", "family"), ("gbif_genus", "gbif_genus"), ("group", "group"),
                 ("gbif_rank", "rank")]:
        E[c] = E.nprops.map(lambda s: jp(s, k))
    E["written"] = E.props.map(lambda s: jp(s, "partner"))
    E["relation"] = E.props.map(lambda s: jp(s, "relation"))
    E["relation_read"] = E.props.map(lambda s: jp(s, "relation_read"))
    E["cat"] = [cat_of(k, p, g, f"{w} {l}") for k, p, g, w, l in zip(E.kingdom, E.phylum, E.group, E.written, E.label)]
    # her-derived rows
    CC = pd.read_parquet(SNAP / "tables/claims_checked.parquet", columns=["claim_id", "doc_id", "quote", "evidence",
                                                                          "relation", "passage_id"])
    cdoc = dict(zip(CC.claim_id, CC.doc_id))
    E["doc_id"] = E.claim_id.map(cdoc)
    her = (E.citation.fillna("").str.contains("via the user's pathogen database") | E.source_id.str.startswith(
        ("notion_pathogen_db", "pref:", "abs:pathref")) | E.doc_id.fillna("").str.startswith("pref_"))
    print("edges derived from her list (excluded):", int(her.sum()))
    E = E[~her].copy()
    E["stratum"] = np.where(E.method.eq("llm_facts+checks"), "literature read by the KG",
                            np.where(E.method.eq("genbank_isolate"), "GenBank isolate record", "other database (GloBI)"))
    n0 = E.groupby("stratum").size().to_dict()
    steps = [("start", E)]
    E1 = E[E.predicate.eq("has_natural_enemy") & E.relation.isin(["pathogen", "parasite"])]
    steps.append(("R1 natural enemy, relation pathogen or parasite", E1))
    E2 = E1[E1.cat.isin(CATS)]
    steps.append(("R2 one of her six categories", E2))
    q = E2.evidence.fillna("")
    taxo = E2.group.fillna("") + " " + E2.order.fillna("") + " " + E2.family.fillna("") + " " + E2.label.fillna("") + \
        " " + E2.written.fillna("")
    harm, tree = q.str.contains(HARM), q.str.contains(TREE)
    fungal_ok = (taxo.str.contains(INSECT_FUNGI) & ~taxo.str.contains(BAD_FUNGI)) | \
                (harm & ~tree & ~taxo.str.contains(BAD_FUNGI))
    bact_ok = harm | taxo.str.contains(r"bacillus|serratia|pseudomonas|paenibacillus|lysinibacillus|brevibacillus|"
                                       r"photorhabdus|xenorhabdus|rickettsiella", case=False)
    keep3 = np.where(E2.cat.eq("Fungi"), fungal_ok, np.where(E2.cat.eq("Bacteria"), bact_ok & ~(tree & ~harm), True))
    E3 = E2[keep3]
    steps.append(("R3 fungi / bacteria only as pathogens of the beetle", E3))
    # R4 sources: reviews and secondary
    SM = con.sql("select source_id, title, year, citation, doi from source_meta").df().drop_duplicates("source_id")
    E3 = E3.merge(SM.rename(columns={"title": "src_title", "year": "src_year", "citation": "src_citation",
                                     "doi": "src_doi"}), on="source_id", how="left")
    srctext = E3.src_title.fillna("") + " " + E3.citation.fillna("")
    is_review = srctext.str.contains(REVIEW) | E3.source_id.str.startswith(("wiki:", "web:")) | \
        E3.citation.fillna("").str.contains("Wood & Bright|Wood and Bright|A catalog of Scolytidae", case=False)
    F = con.sql("select fact_id, supporting_json from facts where fact_id in (select fact_id from edge_facts)").df()
    F = F[F.fact_id.isin(set(E3.fact_id))]
    sh = {}
    for r in F.itertuples():
        for s in json.loads(r.supporting_json or "[]"):
            if s.get("edge_id") is not None:
                sh[int(s["edge_id"])] = s.get("second_hand_of") or []
    E3["second_hand_of"] = E3.edge_id.map(lambda i: sh.get(int(i), []))
    is_second = E3.second_hand_of.map(len).gt(0)
    POINT = E3[~is_review & is_second].copy()
    E4 = E3[~is_review & ~is_second].copy()
    steps.append(("R4 primary source (no review / catalogue / website; not second-hand)", E4))
    # R5 not one of her publications
    ids = json.loads((A / "data/her_pub_kg_source_ids.json").read_text())
    src2pub = {s: p for p, ss in ids.items() for s in ss}
    P = pd.read_parquet(A / "out/her_publications.parquet")
    her_dois = set(P.doi[P.doi.ne("")])
    E4["her_pub"] = E4.source_id.map(src2pub)
    E4.loc[E4.her_pub.isna() & E4.src_doi.fillna("").str.lower().isin(her_dois), "her_pub"] = "doi match"
    # GloBI rows carry the original reference in the citation; match its DOI too
    E4["cit_doi"] = E4.citation.fillna("").str.extract(r"(10\.\d{4,9}/[^\s,;]+)")[0].str.lower().str.rstrip(".")
    E4.loc[E4.her_pub.isna() & E4.cit_doi.isin(her_dois), "her_pub"] = "doi match (citation)"
    # her 4 publications whose only records lack a taxon (not in the 276 analytical publications) also count as hers
    from rapidfuzz import fuzz
    FS = pd.read_csv(A / "figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv", encoding="utf-8-sig")
    FS["ref"] = FS["title"].fillna("").str.replace(r"\s+", " ", regex=True).str.strip()
    extra = FS[~FS.ref.isin(set(P.ref))].drop_duplicates("ref")
    for r in extra.itertuples():
        t = agr.norm(r.ref)
        hit = E4.her_pub.isna() & E4.src_title.fillna("").map(lambda s: len(s) > 25 and fuzz.partial_ratio(agr.norm(s), t) >= 92)
        E4.loc[hit, "her_pub"] = "her file (record without taxon): " + str(r.source)
        print("her no-taxon publication", r.source, "-> KG edges", int(hit.sum()))
    E4[E4.her_pub.fillna("").str.startswith("her file")][
        ["her_pub", "subject", "written", "label", "cat", "relation", "citation", "url", "evidence"]].rename(
        columns={"evidence": "quote_INTERNAL"}).to_csv(A / "out/INTERNAL_her_no_taxon_publications_kg_records.csv",
                                                        index=False)
    E5 = E4[E4.her_pub.isna()].copy()
    # every primary source with an in-scope pathogen / parasite record (her publications included): for the trend step
    E4[E4.stratum.eq("literature read by the KG")].groupby("source_id").agg(
        year=("src_year", "first"), her_pub=("her_pub", "first"), n_edges=("edge_id", "size"),
        categories=("cat", lambda s: "; ".join(sorted(set(s))))).reset_index().to_csv(
        A / "out/kg_pathogen_sources_primary.csv", index=False)
    # also the KG's species-level in-scope records (any source, her-derived excluded) for the surveillance step
    E4.assign(in_her=E4.her_pub.notna())[["subject", "object", "label", "cat", "source_id", "stratum", "in_her"]].to_csv(
        A / "data/kg_inscope_records.csv", index=False)
    steps.append(("R5 source not in her dataset", E5))
    S = []
    for name, d in steps:
        S.append(dict(step=name, **{f"edges: {k}": int((d.stratum == k).sum()) for k in n0},
                      sources=d.source_id.nunique(), beetle_taxa=d.subject.nunique()))
    print(pd.DataFrame(S).to_string())
    labels = {"start": "all organism records (literature + databases)",
              "R1 natural enemy, relation pathogen or parasite": "R1 called pathogen or parasite",
              "R2 one of her six categories": "R2 one of the six categories",
              "R3 fungi / bacteria only as pathogens of the beetle": "R3 fungi / bacteria harm the beetle",
              "R4 primary source (no review / catalogue / website; not second-hand)": "R4 primary source",
              "R5 source not in her dataset": "R5 source not in her dataset"}
    FU = pd.DataFrame(S)
    FU["label"] = FU.step.map(labels)
    FU.to_csv(A / "out/candidate_funnel.csv", index=False)

    # ---------------- compare with her dataset
    R = pd.read_parquet(A / "out/her_records.parquet")
    TB = pd.read_csv(A / "out/taxonomy_beetles.csv")
    sid_of = {r.accepted_in_dataset: str(r.kg_species_id).split(";")[0] for r in TB.itertuples()
              if isinstance(r.kg_species_id, str) and r.kg_species_id}
    # Ips oregoni / oregonis = Ips pini synonyms (taxonomy step) are mapped by epithet when unique
    for r in TB.itertuples():
        if r.accepted_in_dataset not in sid_of and isinstance(getattr(r, "same_or_near_epithet_in_list", ""), str):
            vals = {x.split("->")[-1].strip() for x in r.same_or_near_epithet_in_list.split(";") if "->" in x}
            if len(vals) == 1:
                sid_of[r.accepted_in_dataset] = "name:" + vals.pop()
    SP = con.sql("select species_id, name, genus, tribe, subfamily from species").df()
    name2sid = dict(zip(SP.name, SP.species_id))
    sid_of = {k: (name2sid.get(v[5:], v) if v.startswith("name:") else v) for k, v in sid_of.items()}
    R["sid"] = R.AcceptedBeetleSpeciesName.map(sid_of)
    TP = pd.read_parquet(A / "data/her_pathogen_gbif.parquet").set_index("pathogen_in_dataset")
    R["p_clean"] = R.pathogens.map(clean_name)
    R["p_gen"] = R.pathogens.map(genus_of)
    R["p_gbif_sp"] = R.pathogens.map(lambda n: str(TP.gbif_accepted_species.get(n) or "").lower())
    R["p_gbif_gen"] = R.pathogens.map(lambda n: str(TP.gbif_genus.get(n) or "").lower())
    her_sids = set(R.sid.dropna())
    her_names = set(R.p_clean) | set(R.p_gbif_sp) - {""}
    her_genera = (set(R.p_gen) | set(R.p_gbif_gen)) - {""}
    her_pairs_g = {(s, g) for s, a, b in zip(R.sid, R.p_gen, R.p_gbif_gen) for g in (a, b) if g}
    her_pairs_sp = {(s, n) for s, a, b in zip(R.sid, R.p_clean, R.p_gbif_sp) for n in (a, b) if n}
    sp_info = SP.set_index("species_id")
    for D in (E5, POINT):
        D["beetle"] = D.subject.map(lambda s: sp_info.name.get(s, s))
        D["beetle_subfamily"] = D.subject.map(lambda s: sp_info.subfamily.get(s, ""))
        D["beetle_genus"] = D.subject.map(lambda s: sp_info.genus.get(s, s.split(":")[-1] if s.startswith("genus:") else ""))
        D["p_name"] = [binom_label(l) if r in ("SPECIES", "SUBSPECIES", "VARIETY") else clean_name(w)
                       for l, r, w in zip(D.label, D["rank"], D.written)]
        D["p_genus"] = [str(g).lower() if g else genus_of(w) for g, w in zip(D.gbif_genus, D.written)]
        D["beetle_in_her_data"] = D.subject.isin(her_sids)
        D["pathogen_genus_in_her_data"] = D.p_genus.isin(her_genera)
        D["pathogen_name_in_her_data"] = D.p_name.isin(her_names)
        D["pair_in_her_data_species"] = [(s, n) in her_pairs_sp for s, n in zip(D.subject, D.p_name)]
        D["pair_in_her_data_genus"] = [(s, g) in her_pairs_g for s, g in zip(D.subject, D.p_genus)]
        D["year"] = pd.to_numeric(D.src_year, errors="coerce")
        D["decade"] = (D.year // 10 * 10).astype("Int64")
    E5["quote"] = E5.evidence
    # language of the source (text corpus field, else guess)
    SRC = pd.read_parquet(ST / "dataset/sources.parquet", columns=["source_id", "language"])
    lang = dict(zip(SRC.source_id, SRC.language))
    E5["language"] = E5.source_id.map(lang).fillna("")
    E5["language"] = np.where(E5.language.eq("") & E5.quote.fillna("").str.contains(r"[一-鿿]"), "zh",
                              E5.language)

    # ---------------- focused subset: the source is about bark / ambrosia beetles (its title names a beetle term or the
    # beetle's genus) and the quote is the source's own statement (no citation marker: '(Author 1990)', '[12]',
    # 'et al.'). Checked on the 50-record sample: removes passing mentions, reference-list titles and second-hand
    # statements while keeping most genuine records.
    BEETLE_T = re.compile(r"beetle|scolyt|platypod|ipidae|borer|broca|borkenk|kornik|k[uů]rovec|lýkožrout|escaraba|"
                          r"besouro|scolyte|xylebor|engraver|pinhole|shot[- ]?hole|ambrosia", re.I)
    CIT = re.compile(r"\(\s*[A-Z][^()]{0,80}?\d{4}[a-z]?\s*\)|\[\s*\d+(?:\s*[,–-]\s*\d+)*\s*\]|\bet al\b|"
                     r"\[\s*[A-Z][^\]]{0,80}\d{4}|\(\d+\)\s*(found|reported|showed)|cited|according to", re.I)
    ttl = (E5.src_title.fillna("") + " " + E5.citation.fillna("")).str.slice(0, 400)
    E5["title_on_topic"] = [bool(BEETLE_T.search(t)) or (isinstance(g, str) and len(g) > 3 and g.lower() in t.lower())
                            for t, g in zip(ttl, E5.beetle_genus)]
    E5["quote_cites_other_work"] = E5.quote.fillna("").str.contains(CIT)
    E5["focused"] = E5.title_on_topic & ~E5.quote_cites_other_work

    # ---------------- Platypodinae (analysis 5): all stages
    Pl = E3[E3.subject.map(lambda s: sp_info.subfamily.get(s, "") == "Platypodinae" if s.startswith("species:")
                           else s in ("subfamily:Platypodinae",) or (s.startswith("genus:") and s[6:] in
                                                                     set(SP[SP.subfamily.eq("Platypodinae")].genus)))]
    PlA = E[E.subject.map(lambda s: sp_info.subfamily.get(s, "") == "Platypodinae" if s.startswith("species:")
                          else (s.startswith("genus:") and s[6:] in set(SP[SP.subfamily.eq("Platypodinae")].genus)))]
    PlA = PlA.merge(SM[["source_id", "year"]].rename(columns={"year": "src_year"}), on="source_id", how="left")
    PlA["beetle"] = PlA.subject.map(lambda s: sp_info.name.get(s, s))
    PlA["passes_her_rules_R1_R3"] = PlA.edge_id.isin(set(E3.edge_id))
    PlA["passes_R4"] = PlA.edge_id.isin(set(E4.edge_id))
    PlA[["edge_id", "fact_id", "stratum", "beetle", "predicate", "relation", "relation_read", "written", "label",
         "cat", "group", "citation", "url", "src_year", "passes_her_rules_R1_R3", "passes_R4", "evidence"]].rename(
        columns={"evidence": "quote_INTERNAL"}).to_csv(A / "out/platypodinae_records.csv", index=False)
    print("Platypodinae organism edges (all, not her-derived):", len(PlA), "| by predicate/relation",
          PlA.groupby(["predicate", "relation"]).size().to_dict(), "| passing R1-R3", int(PlA.passes_her_rules_R1_R3.sum()),
          "| passing R1-R4", int(PlA.passes_R4.sum()))
    print(PlA[PlA.cat.isin(CATS)][["beetle", "predicate", "relation", "written", "cat", "citation"]].to_string()[:6000])

    # ---------------- review table
    E5["new_pair"] = ~E5.pair_in_her_data_genus
    E5["new_beetle_species"] = E5.subject.str.startswith("species:") & ~E5.beetle_in_her_data
    E5["new_pathogen_genus"] = ~E5.pathogen_genus_in_her_data & E5.p_genus.ne("")
    E5["beetle_rank"] = E5.subject.str.split(":").str[0]
    cols = ["edge_id", "fact_id", "stratum", "beetle", "beetle_rank", "beetle_genus", "beetle_subfamily", "written",
            "label", "p_name", "p_genus", "cat", "group", "relation", "relation_read", "citation", "src_doi", "url",
            "year", "language", "beetle_in_her_data", "pathogen_name_in_her_data", "pathogen_genus_in_her_data",
            "pair_in_her_data_species", "pair_in_her_data_genus", "new_pair", "new_beetle_species",
            "new_pathogen_genus", "title_on_topic", "quote_cites_other_work", "focused", "source_id", "locator", "quote"]
    RT = E5[cols].rename(columns={"written": "pathogen_as_written", "label": "pathogen_resolved_gbif",
                                  "cat": "category", "src_doi": "doi", "quote": "supporting_quote_INTERNAL_ONLY"})
    RT = RT.sort_values(["new_pair", "category", "beetle"], ascending=[False, True, True])
    RT.insert(0, "INTERNAL", "lab only - quotes are copyrighted; not for publication")
    RT.to_csv(A / "out/INTERNAL_candidate_review_table.csv", index=False)
    POINT[["edge_id", "fact_id", "beetle", "written", "label", "cat", "relation", "citation", "url",
           "second_hand_of", "evidence"]].rename(columns={"evidence": "quote_INTERNAL"}).to_csv(
        A / "out/INTERNAL_candidates_pointers.csv", index=False)

    # ---------------- summaries
    lit = E5[E5.stratum.eq("literature read by the KG")]
    out = []
    for stratum, d in [("literature read by the KG", lit), ("literature, focused subset", lit[lit.focused]),
                       ("other databases (GloBI, GenBank)", E5[~E5.stratum.eq("literature read by the KG")]),
                       ("all", E5)]:
        sp = d[d.beetle_rank.eq("species")]
        newp = sp[sp.new_pair]
        out.append(dict(stratum=stratum, edges=len(d), sources=d.source_id.nunique(),
                        sources_with_new_pair=newp.source_id.nunique(),
                        beetle_species=sp.subject.nunique(), new_beetle_species=sp[sp.new_beetle_species].subject.nunique(),
                        pathogen_names=d.p_name.nunique(), new_pathogen_genera=d[d.new_pathogen_genus].p_genus.nunique(),
                        pairs_species_level=sp.drop_duplicates(["subject", "p_name"]).shape[0],
                        new_pairs_genus_level=newp.drop_duplicates(["subject", "p_genus"]).shape[0],
                        new_pairs_name_level=sp[~sp.pair_in_her_data_species].drop_duplicates(["subject", "p_name"]).shape[0],
                        genus_level_beetle_edges=int(d.beetle_rank.eq("genus").sum())))
    O = pd.DataFrame(out)
    O.to_csv(A / "out/candidates_summary.csv", index=False)
    print(O.to_string())
    CD = []
    for (cat, dec), g in lit[lit.beetle_rank.eq("species")].groupby(["cat", "decade"]):
        CD.append(dict(category=cat, decade=dec, sources=g.source_id.nunique(),
                       new_pairs=g[g.new_pair].drop_duplicates(["subject", "p_genus"]).shape[0],
                       new_beetles=g[g.new_beetle_species].subject.nunique()))
    pd.DataFrame(CD).to_csv(A / "out/candidates_by_category_decade.csv", index=False)
    print(lit[lit.beetle_rank.eq("species")].groupby("cat").agg(
        sources=("source_id", "nunique"), species=("subject", "nunique"),
        new_pairs=("new_pair", "sum")).to_string())
    print("by decade:", lit.groupby("decade").source_id.nunique().to_dict())
    print("language:", lit.drop_duplicates("source_id").language.value_counts().to_dict())

    # ---------------- random sample of 50 for manual reading (literature stratum, species-level, unique source-pair)
    U = lit[lit.beetle_rank.eq("species")].drop_duplicates(["source_id", "subject", "p_name"])
    ctx = dict(zip(CC.claim_id, CC.evidence))
    smp = U.sample(n=min(50, len(U)), random_state=20261023)
    print("sample50: focused", int(smp.focused.sum()), "| focused population", int(U.focused.sum()), "of", len(U))
    # publication-level sample: 40 random sources of the focused subset; all their focused records are printed
    srcs = U[U.focused].drop_duplicates("source_id").sample(n=min(40, U[U.focused].source_id.nunique()),
                                                            random_state=20261024).source_id
    S2 = U[U.focused & U.source_id.isin(set(srcs))].assign(context_INTERNAL=lambda d: d.claim_id.map(ctx))
    S2[["source_id", "edge_id", "beetle", "written", "label", "cat", "relation", "citation", "src_title", "url", "year",
        "new_pair", "evidence", "context_INTERNAL"]].rename(columns={"evidence": "quote_INTERNAL"}).to_csv(
        A / "out/INTERNAL_candidate_sample40_publications.csv", index=False)
    print("publication-level sample: sources", S2.source_id.nunique(), "records", len(S2),
          "| focused population sources", U[U.focused].source_id.nunique())
    smp = smp.assign(context_INTERNAL=smp.claim_id.map(ctx))
    smp[["edge_id", "fact_id", "beetle", "written", "label", "cat", "relation", "relation_read", "citation", "url",
         "year", "new_pair", "evidence", "context_INTERNAL"]].rename(columns={"evidence": "quote_INTERNAL"}).to_csv(
        A / "out/INTERNAL_candidate_sample50.csv", index=False)
    print("sample population (unique source x beetle x organism):", len(U))


if __name__ == "__main__":
    main()
