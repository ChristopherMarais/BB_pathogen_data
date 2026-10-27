"""Step 3 (analysis 1): search recall. Of her 280 primary publications, how many did the KG's OWN literature search
find, without her list?

Independent search stores (none uses her reference list):
  y2025      2025 species-name harvest (BB_papers/download/doi_list.csv + processed_dois.txt; Oct 2025, before her
             list reached the KG) and its full texts (download/full_text_downloads, named by DOI)
  s2026      2026 search hits (BB_papers/search_2026/hits_openalex / hits_epmc / hits_scopus; queries = species names
             and synonyms, genus / tribe names + a beetle term, topic phrases, and from 2026-10-16 fungus / natural-enemy
             taxa + a beetle term ('fungi', 'enemy', 'gfung', 'genemy' passes), taxonomists' surnames ('author'))
  corpus     documents the KG holds and read, from the stores above or from name-driven harvests (BHL, Internet
             Archive, Zobodat, Plazi, Naturalis, author PDFs, EPPO Reporting Service). EXCLUDED as derived from her list:
             pref_ documents (free copies fetched from her list), abs:pathref:* abstracts (looked up from her list),
             notion_pathogen_db:* rows (the Notion copy of her database), and any abstract fetched only because it was
             on her list (sources_abstracts.pathogen_ref True and DOI in no independent store).
  paradigms  a separate OpenAlex title/abstract search of the bark / ambrosia beetle literature made for the
             paradigm analysis (species_kg/analysis/paradigms, 7 queries); reported as a second, simpler search.
Matching: DOI first (exact, lower case). Then titles: candidates published within +-1 year of her year;
retrieval by token-set similarity, accepted when the Levenshtein ratio of the normalised titles is >= 90, or when the
candidate title (>= 60 % of her title's length) occurs in her full reference with partial ratio >= 95 ('title match');
best scores 80-90 are listed as 'possible' for manual review and NOT counted.
Outputs: out/recall_by_publication.csv, out/recall_summary.csv, out/recall_possible_matches.csv
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from rapidfuzz import fuzz, process
from unidecode import unidecode

A = Path(__file__).resolve().parents[2] / "data"
SNAP = Path(__file__).resolve().parents[2] / "data"
BB = Path(__file__).resolve().parents[2] / "data" / "literature_corpus"
ST = Path(__file__).resolve().parents[2] / "data"
PAR = Path(__file__).resolve().parents[2] / "data"


def norm(s):
    s = unidecode(str(s or "")).lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def cdoi(d):
    d = str(d or "").strip().lower().rstrip(".")
    d = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:)", "", d)
    return d if d.startswith("10.") else ""


HER_DERIVED = re.compile(r"^(pref:|abs:pathref|notion_pathogen_db)")


def main():
    P = pd.read_parquet(A / "out/her_publications.parquet")
    print("her publications", len(P), "| with DOI", int(P.doi.ne("").sum()))
    # ------------------------------------------------------------- DOI stores
    doi_src = {}

    def add(dois, tag):
        for d in dois:
            d = cdoi(d)
            if d:
                doi_src.setdefault(d, set()).add(tag)
    dl = pd.read_csv(BB / "download/doi_list.csv").iloc[:, 0]
    add(dl, "y2025_doi_list")
    add(open(BB / "download/processed_dois.txt", errors="replace").read().split(), "y2025_doi_list")
    ft = [p.stem for p in (BB / "download/full_text_downloads").iterdir()]
    # file name = DOI with the first '/' (2025) or every '/' replaced by '_' -> match on the slug instead
    slug = {}
    for f in ft:
        slug[f.lower()] = f
    hits = []
    for nm, src in [("hits_openalex", "openalex"), ("hits_epmc", "europepmc"), ("hits_scopus", "scopus")]:
        H = pd.read_parquet(BB / f"search_2026/{nm}.parquet")
        H["store"] = "s2026_" + src
        hits.append(H[["doi", "title", "year", "pass_", "store", "name"]])
    H = pd.concat(hits, ignore_index=True)
    H["doi"] = H.doi.map(cdoi)
    print("2026 hits", len(H), "| passes", H.pass_.value_counts().to_dict())
    for (st, ps), g in H.groupby(["store", "pass_"]):
        add(g.doi, f"{st}|{ps}")
    C = pd.read_parquet(BB / "search_2026/candidates.parquet")
    add(C.doi, "s2026_candidates")
    # corpus documents
    SRC = pd.read_parquet(ST / "dataset/sources.parquet")
    SA = pd.read_csv(ST / "work/sources_abstracts.csv", usecols=["source_id", "doi", "pathogen_ref"])
    ab_pref = set(SA.loc[SA.pathogen_ref.astype(str).eq("True"), "doi"].map(cdoi))
    SRC["doi_c"] = SRC.doi.map(cdoi)
    SRC["her_derived"] = (SRC.kind.eq("pathrefs_pdf") | SRC.source_id.str.match(HER_DERIVED)
                          | (SRC.kind.eq("abstract") & SRC.doi_c.isin(ab_pref)
                             & ~SRC.doi_c.map(lambda d: bool(d) and bool(doi_src.get(d)))))
    print("corpus sources", len(SRC), "| her-derived", int(SRC.her_derived.sum()),
          SRC[SRC.her_derived].kind.value_counts().to_dict())
    for k, g in SRC[~SRC.her_derived].groupby("kind"):
        add(g.doi_c, "corpus_" + k)
    # her-derived corpus documents (for the reading analysis only, never counted as found)
    her_doc_doi = {}
    for r in SRC[SRC.her_derived].itertuples():
        if r.doi_c:
            her_doc_doi.setdefault(r.doi_c, set()).add("herlist_" + r.kind)
    # paradigms OpenAlex corpus
    PC = pd.read_parquet(PAR / "corpus_raw.parquet", columns=["id", "doi", "title", "year", "type", "language"])
    PC["doi"] = PC.doi.map(cdoi)
    add(PC.doi, "paradigms_openalex")

    # ------------------------------------------------------------- title stores
    T = [H.assign(src=H.store + "|" + H.pass_)[["title", "year", "src", "doi"]],
         C.assign(src="s2026_candidates")[["title", "year", "src", "doi"]],
         SRC[~SRC.her_derived].assign(src="corpus_" + SRC.kind, doi=SRC.doi_c).rename(columns={"pub_year": "year"})[
             ["title", "year", "src", "doi"]],
         PC.assign(src="paradigms_openalex")[["title", "year", "src", "doi"]]]
    # Zobodat / IA candidate lists: name-driven harvests
    try:
        Z = pd.read_csv(ST / "work/zobodat_articles.csv")
        tcol = next(c for c in Z.columns if "title" in c.lower())
        ycol = next((c for c in Z.columns if "year" in c.lower()), None)
        T.append(pd.DataFrame({"title": Z[tcol], "year": Z[ycol] if ycol else np.nan, "src": "zobodat_list", "doi": ""}))
    except Exception as e:
        print("zobodat list", e)
    SM = pd.read_parquet(SNAP / "tables/source_meta.parquet")
    SM_ind = SM[~SM.source_id.str.match(HER_DERIVED)].copy()
    SM_ind["doi"] = SM_ind.doi.map(cdoi)
    T.append(SM_ind.assign(src="kg_source_meta")[["title", "year", "src", "doi"]])
    TT = pd.concat(T, ignore_index=True)
    TT["tn"] = TT.title.map(norm)
    TT = TT[TT.tn.str.len() >= 20]
    TT["year"] = pd.to_numeric(TT.year, errors="coerce")
    TT = TT.drop_duplicates(["tn", "year", "src"])
    print("title store rows", len(TT), TT.src.str.split("|").str[0].value_counts().to_dict())
    by_year = {y: g for y, g in TT.groupby(TT.year.fillna(-1).astype(int))}

    rows, possible = [], []
    for p in P.itertuples():
        t_n, rest_n = norm(p.title_parsed), norm(p.rest)
        res = dict(pub_id=p.pub_id)
        tags = set()
        if p.doi:
            tags |= doi_src.get(p.doi, set())
            s = p.doi.replace("/", "_")
            if s in slug or p.doi.replace("/", "_", 1) in slug:
                tags.add("y2025_full_text")
        res["doi_tags"] = sorted(tags)
        # titles
        cand = pd.concat([by_year.get(y, TT.iloc[:0]) for y in (p.year - 1, p.year, p.year + 1, -1)])
        title_tags, best_sc, best_t = set(), 0, ""
        if len(t_n) >= 15 and len(cand):
            choices = cand.tn.tolist()
            top = process.extract(t_n, choices, scorer=fuzz.token_set_ratio, limit=25, score_cutoff=70)
            for _, _, idx in top:
                c = cand.iloc[idx]
                r1 = fuzz.ratio(c.tn, t_n)
                r2 = fuzz.partial_ratio(c.tn, rest_n) if len(c.tn) >= 0.6 * len(t_n) else 0
                sc = max(r1, r2 - 3)  # partial matches count slightly less
                if sc >= 90 or (r2 >= 95):
                    title_tags.add(c.src)
                    if c.doi and c.doi != p.doi:
                        res.setdefault("title_match_dois", set()).add(c.doi)
                if sc > best_sc:
                    best_sc, best_t = sc, c.title
            if 80 <= best_sc < 90 and not title_tags:
                possible.append(dict(pub_id=p.pub_id, her_reference=p.ref, best_candidate=best_t, score=best_sc))
        res["title_tags"] = sorted(title_tags)
        res["best_title_score"] = best_sc
        res["best_title"] = best_t
        res["her_list_docs"] = sorted(her_doc_doi.get(p.doi, set())) if p.doi else []
        rows.append(res)
    R = pd.DataFrame(rows)
    R["title_match_dois"] = R.get("title_match_dois", pd.Series([None] * len(R))).map(
        lambda x: ";".join(sorted(x)) if isinstance(x, set) else "")
    R = P.drop(columns=["rest"]).merge(R, on="pub_id")
    allt = pd.Series([set(a) | set(b) for a, b in zip(R.doi_tags, R.title_tags)], index=R.index)

    def has(prefixes, exclude=()):
        return allt.map(lambda s: any(t.startswith(prefixes) and not t.startswith(tuple(exclude)) for t in s))
    R["found_y2025"] = has(("y2025",))
    R["found_s2026_names_topics"] = allt.map(lambda s: any(t.startswith("s2026") and t.split("|")[-1] in
                                                           ("syn", "new", "topic", "group", "valid", "recent")
                                                           for t in s))
    R["found_s2026_any"] = has(("s2026",))
    R["found_s2026_enemy_fungus"] = allt.map(lambda s: any(t.startswith("s2026") and t.split("|")[-1] in
                                                           ("fungi", "enemy", "gfung", "genemy") for t in s))
    R["in_corpus_independent"] = has(("corpus_", "kg_source_meta", "zobodat_list"))
    R["found_kg_search"] = R.found_y2025 | R.found_s2026_any | R.in_corpus_independent
    R["found_kg_search_without_enemy_pass"] = R.found_y2025 | R.found_s2026_names_topics | allt.map(
        lambda s: any(t.startswith(("corpus_", "kg_source_meta", "zobodat_list")) for t in s))
    R["found_paradigms_openalex"] = has(("paradigms",))
    R["found_any"] = R.found_kg_search | R.found_paradigms_openalex
    R["match_quality"] = np.select(
        [R.doi_tags.str.len().gt(0), R.title_tags.str.len().gt(0) & R.best_title_score.ge(97),
         R.title_tags.str.len().gt(0)],
        ["DOI", "title >= 97", "title 90-97"], "none")
    print("pass names seen:", sorted({t.split('|')[-1] for s in allt for t in s if t.startswith('s2026')}))
    # why missed
    def why(r):
        if r.found_kg_search:
            return ""
        if not r.oa_found:
            if r.pub_type in ("thesis", "report / proceedings", "book / chapter"):
                return f"grey literature ({r.pub_type}), no DOI-indexed record"
            if r.language_guess != "English":
                return f"non-English ({r.language_guess}), no DOI-indexed record"
            return "journal article without a DOI-indexed record (old or regional journal)"
        if not r.oa_has_abstract and not r.oa_beetle_term:
            return "in OpenAlex but no abstract and title without a bark-beetle term or name"
        if not r.oa_beetle_term:
            return "in OpenAlex; title / abstract name no bark-beetle term searched"
        return "in OpenAlex with a beetle term, but not retrieved (query caps or index of the year)"
    R["why_missed"] = R.apply(why, axis=1)
    R["doi_tags"] = R.doi_tags.map(";".join)
    R["title_tags"] = R.title_tags.map(";".join)
    R["her_list_docs"] = R.her_list_docs.map(";".join)
    R.to_csv(A / "out/recall_by_publication.csv", index=False)
    pd.DataFrame(possible).to_csv(A / "out/recall_possible_matches.csv", index=False)

    # ------------------------------------------------------------- summaries
    def wilson(k, n, z=1.96):
        if n == 0:
            return (np.nan, np.nan)
        p = k / n
        d = 1 + z * z / n
        c = (p + z * z / (2 * n)) / d
        h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
        return (c - h, c + h)
    out = []
    cats = sorted({c for s in R.categories for c in s.split("; ")})
    groups = [("all", R.index == R.index)]
    groups += [(f"decade {d}", R.decade.eq(d)) for d in sorted(R.decade.unique())]
    groups += [(f"period {a}-{b}", R.year.between(a, b)) for a, b in [(1914, 1969), (1970, 1999), (2000, 2026)]]
    groups += [(f"language {l}", R.language_guess.eq(l)) for l in R.language_guess.value_counts().index]
    groups += [(f"type {t}", R.pub_type.eq(t)) for t in R.pub_type.value_counts().index]
    groups += [(f"category {c}", R.categories.str.contains(c)) for c in cats]
    groups += [(f"full text flag {f}", R.full_text.eq(f)) for f in R.full_text.value_counts().index]
    real_doi = R.doi.fillna("").ne("") & ~R.doi.fillna("").str.startswith("10.5555/")
    groups += [("has article DOI", real_doi), ("no article DOI (none or CAB Abstracts record only)", ~real_doi),
               ("OpenAlex record by DOI", R.oa_found.fillna(False)),
               ("no OpenAlex record by DOI", ~R.oa_found.fillna(False))]
    for name, m in groups:
        g = R[m]
        n = len(g)
        row = dict(group=name, n_publications=n)
        for col in ["found_kg_search", "found_kg_search_without_enemy_pass", "found_y2025", "found_s2026_any",
                    "found_s2026_enemy_fungus", "in_corpus_independent", "found_paradigms_openalex", "found_any"]:
            k = int(g[col].sum())
            lo, hi = wilson(k, n)
            row[col] = k
            row[col + "_pct"] = round(100 * k / n, 1) if n else np.nan
            if col == "found_kg_search":
                row["ci95_lo"], row["ci95_hi"] = round(100 * lo, 1), round(100 * hi, 1)
        out.append(row)
    O = pd.DataFrame(out)
    O.to_csv(A / "out/recall_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(O[["group", "n_publications", "found_kg_search", "found_kg_search_pct", "ci95_lo", "ci95_hi",
             "found_kg_search_without_enemy_pass", "found_y2025", "found_s2026_any", "in_corpus_independent",
             "found_paradigms_openalex", "found_any"]].to_string())
    print("match quality:", R.match_quality.value_counts().to_dict())
    print("why missed:", R.why_missed.value_counts().to_dict())
    print("possible (unconfirmed) title matches:", len(possible))
    print(R[~R.found_kg_search][["pub_id", "source", "year", "pub_type", "language_guess", "doi", "oa_found",
                                  "why_missed", "best_title_score"]].to_string())


if __name__ == "__main__":
    main()
