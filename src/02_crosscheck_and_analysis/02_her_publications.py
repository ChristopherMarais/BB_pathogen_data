"""Step 2: Dong & Hulcr's dataset (figshare v1.1) -> one row per primary publication, with identifiers.

Publication = distinct full reference ('title' column). For each: first author, year, title (parsed from the reference),
publication type and language (inferred from the reference text; rules below), categories, records, full-text flag.
DOI, in order of trust:
  1. a DOI written in the reference;
  2. the DOI the KG found for the same reference of her list (tables/pathogen_refs.parquet: Crossref, accepted only when
     the Crossref title occurs in the reference) - an identifier of HER paper, not KG search evidence;
  3. the DOI of the same reference in the lab's Notion copy of her database (notion_pathogens: record DOI);
  4. Crossref bibliographic query (keyless; accepted when title similarity >= 90, year +-1 and first-author surname).
OpenAlex lookup (keyless; DOI, else title search with year +-1 and surname): indexed or not, OpenAlex language, type,
whether OpenAlex has an abstract, and whether title / abstract name a bark or ambrosia beetle term (what the KG's
searches query). Cache: data/lookup_cache.jsonl.
Outputs: out/her_records.parquet (analytical records), out/her_publications.parquet / .csv.
"""
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests
from rapidfuzz import fuzz
from unidecode import unidecode

A = Path(__file__).resolve().parents[2] / "data"
SNAP = Path(__file__).resolve().parents[2] / "data"
CACHE = A / "data/lookup_cache.jsonl"
S = requests.Session()
S.headers["User-Agent"] = "UF-Hulcr-lab pathogen dataset cross-check (research)"
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s,;\"<>]+", re.I)


def norm(s):
    s = unidecode(str(s or "")).lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def clean_doi(d):
    d = str(d or "").strip().lower().rstrip(".").rstrip(")")
    d = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:)", "", d)
    return d if d.startswith("10.") else ""


def parse_ref(ref):
    """-> authors part, year, rest (title + venue), title (first sentence of rest)."""
    r = re.sub(r"\s+", " ", str(ref)).strip()
    m = re.match(r"^(.*?)\(\s*(\d{4})[a-z]?\s*\)\.?\s*(.*)$", r)  # APA: Authors (2013). Title.
    if not m:
        m = re.match(r"^(.*?)[,.]?\s(\d{4})[a-z]?[.,:]\s*(.*)$", r)  # Elsevier: Authors, 1970. Title.
    if m:
        auth, yr, rest = m.group(1), int(m.group(2)), m.group(3)
    else:
        auth, yr, rest = r[:40], None, r
    # title = up to the first '. ' that is followed by a capital word and not inside an abbreviation of a name
    t = rest
    for mm in re.finditer(r"[.?!]\s+(?=[A-Z0-9(\[])", rest):
        cand = rest[:mm.start()]
        if len(cand) >= 15 and not re.search(r"\b[A-Z][a-z]?$|\b(Bals|Vuill|sp|spp|var|nov|n|f|L|Ratz|Gyll)$", cand):
            t = cand
            break
    surname = re.split(r"[,\s]", unidecode(auth).strip())[0] if auth else ""
    return auth.strip(), yr, rest.strip(), t.strip()


THESIS = re.compile(r"thesis|dissertation|ph\.? ?d\b|m\.? ?sc\b|master|diplom|doctoral|habilitation|disertac", re.I)
GREY = re.compile(r"\breport|res\. ?pap|research paper|research note|gen(eral)?\.? tech|technical|bulletin|proceedings|"
                  r"proc\.|conference|symposium|workshop|meeting|congress|abstracts? of|iobc|wprs|usda|forest service|"
                  r"extension|circular|leaflet|newsletter|annual|verhandlung|tagung|kongress|sbornik|zbornik", re.I)
BOOK = re.compile(r"\bIn:|\bIn [A-Z][^.]{0,80}\(eds?\.?\)|\(eds?\.?\)|\beds?\.,|\bpress\b|publisher|verlag|\bbook\b|"
                  r"\bpp\. ?\d+\s*$|academic|elsevier\b(?!.*journal)|springer(?!.*journal)", re.I)
ZH_J = re.compile(r"sinica|sinicae|chinese|china|zhongguo|xuebao|kunchong|linye|forest pest and disease|"
                  r"northeast forestry|beijing forestry|nanjing forestry|southwest forestry|northwest|central south|"
                  r"fujian|zhejiang|sichuan|yunnan|guangxi|shandong|hebei|henan|hunan|xinjiang|inner mongolia|"
                  r"journal of forestry research|forest research|plant protection|entomological knowledge|"
                  r"chinese bulletin|scientia|anhui|jiangsu|jilin|heilongjiang|shaanxi|gansu|qinghai|ningxia|tibet", re.I)
DE_W = re.compile(r"\b(über|uber|und|der|die|das|des|bei|zur|zum|von|einer|eines|untersuchungen|borkenk|"
                  r"anzeiger|zeitschrift|forstwiss|beiträge|beitrage|mitteilungen|anz\.|schädlingsk)\b", re.I)
RU_W = re.compile(r"[а-яА-Я]|\b(lesnoe|zhurnal|izvestiya|entomologicheskoe|obozrenie|trudy|zoologicheskii|"
                  r"parazitologiya|lesovedenie|vestnik)\b", re.I)
OTHER_W = {"cs": re.compile(r"\b(lesnická|lesnicka|práce|prace|kůrovec|kurovec|lýkožrout|lykozrout|zprávy|ochrana)\b", re.I),
           "pl": re.compile(r"\b(sylwan|leśne|lesne|polskie|kornik|badania)\b", re.I),
           "tr": re.compile(r"\b(dergisi|üzerinde|uzerinde|orman|türkiye|turkiye)\b", re.I),
           "fr": re.compile(r"\b(sur|les|des|une|scolyte|annales des sciences forestières|revue)\b", re.I),
           "es": re.compile(r"\b(sobre|del|los|las|una|revista|escarabajo)\b", re.I),
           "it": re.compile(r"\b(sulla|degli|della|nelle|redia|bollettino)\b", re.I),
           "ro": re.compile(r"\b(revista pădurilor|padurilor|gândacul|silvicultur)\b", re.I),
           "bg": re.compile(r"\b(nauka za gorata|gorata)\b", re.I),
           "ka": re.compile(r"\b(georgian|sakartvelos)\b", re.I),
           "ja": re.compile(r"[぀-ヿ]|nihon|nippon|japanese", re.I),
           "ko": re.compile(r"[가-힯]|korean", re.I)}


def pub_type(ref):
    if THESIS.search(ref):
        return "thesis"
    if GREY.search(ref):
        return "report / proceedings"
    if BOOK.search(ref):
        return "book / chapter"
    return "journal article"


DE_T = re.compile(r"\b(über|uber|und|der|die|zur|zum|bei|von|einer|eines|einiger|untersuchungen|untersuchung|beitrag|"
                  r"borkenkäfer\w*|borkenkafer\w*|krankheit\w*|wirkung|auftreten|gattung|kenntnis|möglichkeiten|"
                  r"vorkommen|erhebungen|häufigkeit|verbreitung|gegenspieler\w*|naturgeschichte|parasiten|"
                  r"pathogene|fichtenborkenkäfer\w*|experimentelle|insektenpathogenen)\b", re.I)
FR_T = re.compile(r"\b(les|des|du|sur|une|étude|etude|ennemis|naturels|chez|variations|dans|forêt|"
                  r"épizooties|epizooties|bactériennes|élevages|insectes|encyclopedie|entomologique)\b", re.I)


def language_guess(ref, title, oa_lang=None):
    """Language of the publication, inferred: script / journal words of the reference, title words (>= 2 function
    words of one language), OpenAlex's language when it is not English, Chinese journal names."""
    if re.search(r"[一-鿿]", ref) or re.search(r"\(in chinese|chinese with english", ref, re.I):
        return "Chinese"
    if RU_W.search(ref):
        return "Russian"
    if re.search(r"[łśżźąęćńŁŚŻŹ]", title) or len(OTHER_W["pl"].findall(title)) >= 2:
        return "Polish"
    if re.search(r"[řůěŘ]", title):
        return "Czech"
    fr, de = len(FR_T.findall(title)), len(DE_T.findall(title))
    if re.search(r"rev\.?\s?path\.?\s?v[ée]g|encyclop[ée]die entomologique", ref, re.I) or (fr >= 2 and fr > de):
        return "French"
    if de >= 2:
        return "German"
    # OpenAlex's language field was checked and found unreliable for these titles (English papers labelled es / it);
    # it is not used.
    for code, rx in OTHER_W.items():
        if code == "fr":
            continue
        if len(rx.findall(title)) >= 2:
            return {"cs": "Czech", "pl": "Polish", "tr": "Turkish", "es": "Spanish", "it": "Italian",
                    "ro": "Romanian", "bg": "Bulgarian", "ka": "Georgian", "ja": "Japanese", "ko": "Korean"}[code]
    if ZH_J.search(ref) and not re.search(r"journal of (invertebrate|applied|economic)|insect science|"
                                          r"forests\b|insects\b|pest management|plant protection research|"
                                          r"romanian", ref, re.I):
        return "Chinese (journal; English title)"
    return "English"


def load_cache():
    C = {}
    if CACHE.exists():
        for line in open(CACHE):
            r = json.loads(line)
            if isinstance(r["v"], dict) and r["v"].get("_http") not in (None, 404):
                continue  # transient failures are retried
            C[r["k"]] = r["v"]
    return C


CC = load_cache()


def cached(k, fn):
    if k in CC:
        return CC[k]
    v = fn()
    CC[k] = v
    with open(CACHE, "a") as fh:
        fh.write(json.dumps({"k": k, "v": v}) + "\n")
    return v


def get_json(url, params=None, tries=4, timeout=60):
    for k in range(tries):
        try:
            r = S.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return {"_http": 404}
            if r.status_code in (429, 500, 502, 503):
                time.sleep(4 * (k + 1))
                continue
            return {"_http": r.status_code}
        except Exception:
            time.sleep(4 * (k + 1))
    return {"_http": "failed"}


def crossref_find(ref, title, yr, surname):
    def f():
        time.sleep(0.25)
        j = get_json("https://api.crossref.org/works", {"query.bibliographic": ref[:300], "rows": 5,
                                                       "select": "DOI,title,issued,author,container-title,type"})
        return (j.get("message") or {}).get("items", []) if isinstance(j, dict) else []
    items = cached("cr|" + ref, f)
    best = None
    for it in items:
        t = (it.get("title") or [""])[0]
        sc = max(fuzz.ratio(norm(t), norm(title)), fuzz.partial_ratio(norm(t), norm(ref)) if len(t) > 25 else 0)
        y = ((it.get("issued") or {}).get("date-parts") or [[None]])[0][0]
        fam = [norm(a.get("family", "")) for a in it.get("author", [])]
        ok_y = yr is None or (y is not None and abs(int(y) - yr) <= 1)
        ok_a = not surname or any(norm(surname) and norm(surname) in f_ for f_ in fam) or not fam
        if sc >= 90 and ok_y and ok_a and (best is None or sc > best[1]):
            best = (it["DOI"].lower(), sc, t)
    return best


BEETLE_RE = re.compile(r"bark[- ]beetle|ambrosia[- ]beetle|scolyt|platypod|ipidae|borkenk|pinhole|shot[- ]?hole|"
                       r"\bips\b|dendroctonus|pine beetle|engraver|coffee berry borer", re.I)


OA_FAIL = [1]  # consecutive failed title searches; after 6 the title search is skipped (circuit breaker)


def openalex_lookup(doi, title, yr, surname):
    def by_doi():
        return get_json(f"https://api.openalex.org/works/doi:{doi}",
                        {"select": "id,doi,title,publication_year,type,language,abstract_inverted_index,authorships"})
    def by_title():
        t = re.sub(r"\s+", " ", re.sub(r"[^\w\s]", " ", title)).strip()[:200]
        f = f"title.search:{t}" + (f",publication_year:{yr - 1}-{yr + 1}" if yr else "")
        return get_json("https://api.openalex.org/works", {"filter": f, "per-page": 5,
                        "select": "id,doi,title,publication_year,type,language,abstract_inverted_index,authorships"},
                        tries=2, timeout=30)
    if doi:
        j = cached("oa|doi|" + doi, by_doi)
        if isinstance(j, dict) and j.get("id"):
            return j, "doi"
    if len(title) >= 15 and OA_FAIL[0] < 1:  # keyless title search timed out on 2026-10-23: skipped
        j = cached("oa|t|" + title, by_title)
        OA_FAIL[0] = OA_FAIL[0] + 1 if (isinstance(j, dict) and j.get("_http") == "failed") else 0
        for w in (j.get("results") or []) if isinstance(j, dict) else []:
            sc = fuzz.ratio(norm(w.get("title")), norm(title))
            y = w.get("publication_year")
            names = " ".join(norm((a.get("author") or {}).get("display_name", "")) for a in w.get("authorships", []))
            if sc >= 90 and (yr is None or (y and abs(y - yr) <= 1)) and (not surname or norm(surname) in names):
                return w, f"title({sc:.0f})"
        if isinstance(j, dict) and j.get("_http"):
            return {"_http": j["_http"]}, "error"
    return None, "not_found"


def main():
    D = pd.read_csv(A / "figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv", encoding="utf-8-sig")
    D.columns = [c.strip() for c in D.columns]
    print("figshare rows", len(D))
    D["pathogens"] = D["pathogens"].fillna("").str.strip()
    no_taxon = D.pathogens.eq("") | D.categories.isna()
    print("rows without pathogen taxon or category:", int(no_taxon.sum()))
    print(D[no_taxon][["Records ID", "AcceptedBeetleSpeciesName", "categories", "pathogens", "source"]].to_string())
    R = D[~no_taxon].copy()
    R["ref"] = R["title"].fillna("").str.replace(r"\s+", " ", regex=True).str.strip()
    print("analytical records", len(R), "| publications (distinct references)", R.ref.nunique(),
          "| distinct 'source' short refs", R.source.nunique(), "| beetle species", R.AcceptedBeetleSpeciesName.nunique(),
          "| genera", R.genus.nunique(), "| pathogen taxa", R.pathogens.nunique(),
          "| pathogen taxa (case-insensitive, trimmed)", R.pathogens.str.lower().str.strip().nunique())
    print("subfamily", R.subfamily.value_counts().to_dict())
    print("categories x records / taxa:\n", R.groupby("categories").agg(records=("pathogens", "size"),
                                                                       taxa=("pathogens", "nunique")).to_string())
    print("year range", R.year.min(), R.year.max())
    print("full text flag", R["full text"].value_counts(dropna=False).to_dict())
    # same short source with different references (and vice versa): possible splits / duplicates
    g = R.groupby("source").ref.nunique()
    print("short refs mapping to >1 full reference:", g[g > 1].to_dict())
    g2 = R.groupby("ref").source.nunique()
    print("full references with >1 short ref:", int((g2 > 1).sum()))
    R.to_parquet(A / "out/her_records.parquet", index=False)

    P = (R.groupby("ref").agg(source=("source", "first"), year=("year", "first"), years=("year", lambda s: sorted(set(s))),
                              categories=("categories", lambda s: "; ".join(sorted(set(s)))),
                              n_records=("pathogens", "size"), n_beetles=("AcceptedBeetleSpeciesName", "nunique"),
                              n_pathogens=("pathogens", "nunique"),
                              full_text=("full text", lambda s: "; ".join(sorted(set(map(str, s))))),
                              record_ids=("Records ID", lambda s: list(map(int, s)))).reset_index())
    P["pub_id"] = [f"P{i + 1:03d}" for i in range(len(P))]
    parsed = P.ref.map(parse_ref)
    P["authors"] = [p[0] for p in parsed]
    P["ref_year"] = [p[1] for p in parsed]
    P["rest"] = [p[2] for p in parsed]
    P["title_parsed"] = [p[3] for p in parsed]
    P["first_author"] = P.authors.map(lambda a: re.split(r"[,\s]", unidecode(a).strip())[0] if a else "")
    P["pub_type"] = P.ref.map(pub_type)
    P["language_guess"] = [language_guess(r, t) for r, t in zip(P.ref, P.title_parsed)]
    P["decade"] = (P.year // 10 * 10).astype(int)
    # ---- DOIs
    P["doi_in_ref"] = P.ref.map(lambda r: clean_doi((DOI_RE.findall(r) or [""])[0]))
    PR = pd.read_parquet(SNAP / "tables/pathogen_refs.parquet")
    PR = PR[PR.doi.fillna("").ne("") & PR.author_ok & PR.year_ok]
    def pr_match(row):
        best = ("", 0)
        for r in PR[(PR.year - row.year).abs() <= 1].itertuples():
            if norm(r.first_author) and norm(r.first_author) not in norm(row.authors):
                continue
            sc = fuzz.token_set_ratio(norm(r.rest), norm(row.rest))
            if sc > best[1]:
                best = (clean_doi(r.doi), sc)
        return best if best[1] >= 90 else ("", best[1])
    m = P.apply(pr_match, axis=1)
    P["doi_kg_reflist"] = [x[0] for x in m]
    NP = pd.read_parquet(SNAP / "tables/notion_pathogens.parquet")
    NP = NP[NP.link_how.eq("record_doi") & NP.source_id.str.startswith("doi:")]
    NP["ref"] = NP.citation.str.replace(r"\s*\(via the user's pathogen database, record [\d.]+\)\s*$", "", regex=True)
    NP["ref_n"] = NP.ref.map(norm)
    nmap = NP.groupby("ref_n").source_id.agg(lambda s: s.value_counts().index[0]).to_dict()
    P["doi_notion"] = P.ref.map(lambda r: clean_doi(nmap.get(norm(r), "")))
    # Crossref for the rest
    need = P.doi_in_ref.eq("") & P.doi_kg_reflist.eq("") & P.doi_notion.eq("")
    print("publications needing a Crossref lookup:", int(need.sum()))
    cr = {}
    for row in P[need].itertuples():
        b = crossref_find(row.ref, row.title_parsed, row.year, row.first_author)
        if b:
            cr[row.pub_id] = b
    P["doi_crossref"] = P.pub_id.map(lambda k: cr.get(k, ("",))[0])
    P["crossref_title_score"] = P.pub_id.map(lambda k: cr.get(k, (None, None))[1])
    P["crossref_title"] = P.pub_id.map(lambda k: cr.get(k, (None, None, ""))[2])
    def pick(r):
        for c, how in [("doi_in_ref", "in reference"), ("doi_kg_reflist", "KG reference list (Crossref, verified)"),
                       ("doi_notion", "Notion record"), ("doi_crossref", "Crossref query")]:
            if r[c]:
                return pd.Series([r[c], how])
        return pd.Series(["", "none"])
    P[["doi", "doi_how"]] = P.apply(pick, axis=1)
    # disagreements between DOI routes
    dis = P[[c for c in ["doi_in_ref", "doi_kg_reflist", "doi_notion", "doi_crossref"]]].apply(
        lambda r: len({x for x in r if x}) > 1, axis=1)
    print("publications whose DOI routes disagree:", int(dis.sum()))
    print(P[dis][["pub_id", "source", "doi_in_ref", "doi_kg_reflist", "doi_notion", "doi_crossref"]].to_string())
    # ---- OpenAlex
    oa = []
    for row in P.itertuples():
        w, how = openalex_lookup(row.doi, row.title_parsed, row.year, row.first_author)
        if w and w.get("id"):
            inv = w.get("abstract_inverted_index") or {}
            ab = " ".join(k for k, _ in sorted(((k, p) for k, ps in inv.items() for p in ps), key=lambda x: x[1]))
            txt = (w.get("title") or "") + " " + ab
            oa.append(dict(pub_id=row.pub_id, oa_found=True, oa_how=how, oa_id=w["id"], oa_doi=clean_doi(w.get("doi")),
                           oa_title=w.get("title"), oa_year=w.get("publication_year"), oa_type=w.get("type"),
                           oa_language=w.get("language"), oa_has_abstract=len(ab) > 50,
                           oa_beetle_term=bool(BEETLE_RE.search(txt))))
        else:
            oa.append(dict(pub_id=row.pub_id, oa_found=False, oa_how=how if w is None else f"error {w.get('_http')}"))
    P = P.merge(pd.DataFrame(oa), on="pub_id", how="left")
    # OpenAlex DOI fills a missing DOI only when found by strict title match
    fill = P.doi.eq("") & P.oa_found.fillna(False) & P.oa_doi.fillna("").ne("")
    P.loc[fill, "doi"] = P.loc[fill, "oa_doi"]
    P.loc[fill, "doi_how"] = "OpenAlex title match"
    P["language_guess"] = [language_guess(r, t, l) for r, t, l in zip(P.ref, P.title_parsed, P.oa_language)]
    # Chinese-language journals: references give English titles; a first author with a pinyin surname and no
    # article DOI is taken as a Chinese-language paper (inference, labelled as such)
    PINYIN = {"huang", "liu", "li", "zhi", "fan", "lin", "xue", "xiao", "nie", "zhang", "hu", "tang", "du", "wu", "yang",
              "wang", "chen", "zhao", "sun", "zhou", "xu", "gao", "ma", "luo", "song", "zheng", "han", "feng", "deng",
              "cao", "peng", "zeng", "xie", "dong", "yuan", "pan", "cai", "jiang", "yu", "ye", "wei", "lu", "shen",
              "qin", "ren", "yan", "lei", "tian", "jin", "duan", "hao", "kong", "bai", "mao", "qiu", "qian", "yin",
              "zou", "shi", "wen", "fu", "su", "lv", "ding", "gu", "jia", "xia", "fang", "kang", "gong", "cheng",
              "yao", "tan", "liao", "zhu", "he", "guo", "ji", "hou", "meng", "bao", "lai", "zhong", "qi", "chang"}
    no_article_doi = P.doi.eq("") | P.doi.str.startswith("10.5555/")
    zh = P.first_author.str.lower().isin(PINYIN) & no_article_doi & P.language_guess.eq("English")
    P.loc[zh, "language_guess"] = "Chinese (inferred: pinyin author, no DOI)"
    P.loc[P.language_guess.eq("Chinese (journal; English title)"), "language_guess"] = "Chinese (inferred: journal name)"
    # CAB Abstracts record DOIs (10.5555/...) identify an abstract record, not the paper
    P["doi_is_cab_abstract"] = P.doi.str.startswith("10.5555/")
    print("DOI routes:", P.doi_how.value_counts().to_dict(), "| CAB Abstracts record DOIs:",
          int(P.doi_is_cab_abstract.sum()))
    print("OpenAlex:", P.oa_how.astype(str).str.split("(").str[0].value_counts().to_dict(),
          "| language", P.oa_language.value_counts(dropna=False).to_dict())
    print("pub types:", P.pub_type.value_counts().to_dict())
    print("language guess:", P.language_guess.value_counts().to_dict())
    print("decades:", P.decade.value_counts().sort_index().to_dict())
    P["years"] = P.years.map(lambda x: ",".join(map(str, x)))
    P["record_ids"] = P.record_ids.map(lambda x: ",".join(map(str, x)))
    P.to_parquet(A / "out/her_publications.parquet", index=False)
    P.drop(columns=["rest"]).to_csv(A / "out/her_publications.csv", index=False)
    print(P[["pub_id", "source", "year", "pub_type", "language_guess", "doi", "doi_how", "oa_how"]].to_string())


if __name__ == "__main__":
    main()
