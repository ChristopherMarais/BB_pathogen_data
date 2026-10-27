"""Two-round programmatic literature search across OpenAlex, Europe PMC, Semantic Scholar, and Scopus
for bark and ambrosia beetle (Curculionidae: Scolytinae and Platypodinae) taxa, synonyms, and pathogen/parasite terms.
"""
import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).parent))
from names_lib import load_names  # noqa: E402

BB = Path(__file__).resolve().parents[2] / "data" / "literature_corpus"
OUT = BB / "search_2026"
TODAY = date.today().isoformat()
SINCE = "2025-10-01"
S = requests.Session()
S.headers["User-Agent"] = "UF-Hulcr-lab bark beetle literature search (research)"


def get(url, params):
    for k in range(6):
        try:
            r = S.get(url, params=params, timeout=90)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 429 and "api.openalex.org" in url and "budget" in r.text.lower():
                return {"_http": 429}
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(5 * (k + 1))
                continue
            return {"_http": r.status_code}
        except (requests.RequestException, ValueError):
            time.sleep(5 * (k + 1))
    return {"_http": "failed"}


CTX = "(beetle OR beetles OR Scolytinae OR Platypodinae OR Scolytidae OR Platypodidae OR Curculionidae)"
TOPICS = ["bark beetle", "bark beetles", "ambrosia beetle", "ambrosia beetles", "Scolytinae", "Platypodinae",
          "Scolytidae", "Platypodidae", "Ipidae", "Borkenkäfer", "Borkenkaefer", "Ambrosiakäfer", "scolyte",
          "scolytes", "escarabajo de corteza", "escarabajos de corteza", "escarabajo ambrosial",
          "besouro de ambrosia", "coleobrocas", "pinhole borer", "shot hole borer", "shot-hole borer",
          "timber beetle", "twig beetle", "coffee berry borer", "mycangia", "mycangium", "xylomycetophagy"]


# Fungus and natural-enemy passes (user 2026-10-16: "there has to be more mention of fungi associated with these
# beetles especially by beetle species"): fungal / enemy taxa and phrases with a bark or ambrosia beetle term, and every
# beetle genus with fungus words (gfung) and with enemy words (genemy).
CTX_BB = ('("bark beetle" OR "bark beetles" OR "ambrosia beetle" OR "ambrosia beetles" OR Scolytinae OR Platypodinae '
          'OR Scolytidae OR Platypodidae OR Ipidae)')
FUNGI = ["Raffaelea", "Ambrosiella", "Dryadomyces", "Harringtonia", "Afroraffaelea", "Meredithiella",
         "Phialophoropsis", "Toshionella", "Wolfgangiella", "Ophiostoma", "Ceratocystiopsis", "Leptographium",
         "Grosmannia", "Graphilbum", "Hyalorhinocladiella", "Sporothrix", "Esteya", "Ceratocystis", "Endoconidiophora",
         "Huntiella", "Ambrosiozyma", "Geosmithia", "Fusarium", "Neocosmospora", "Graphium", "Entomocorticium",
         "Flavodon", "Irpex", "Ogataea", "Kuraishia", "Cyberlindnera", "Nakazawaea", "Wickerhamomyces", "Candida",
         "Pichia", "yeasts", "ophiostomatoid", "ophiostomatoid fungi", "blue-stain fungi", "bluestain fungi",
         "blue stain", "ambrosia fungi", "ambrosia fungus", "fungal symbiont", "fungal symbionts", "fungal associates",
         "fungal associate", "mycangial fungi", "mycangium", "mycangia", "fungus farming", "fungal community",
         "mycobiome", "Ophiostomatales", "Microascales", "symbiotic fungi"]
ENEMIES = ["entomopathogenic fungi", "entomopathogenic fungus", "Beauveria", "Metarhizium", "Isaria", "Cordyceps",
           "Lecanicillium", "microsporidia", "Nosema", "Unikaryon", "Chytridiopsis", "Canningia", "gregarines",
           "Mattesia", "Menzbieria", "nematodes", "Contortylenchus", "Parasitorhabditis", "Parasitylenchus",
           "Sphaerularioides", "Ektaphelenchus", "Bursaphelenchus", "Cryptaphelenchus", "mites", "phoretic mites",
           "Pyemotes", "Tarsonemus", "Iponemus", "Dendrolaelaps", "parasitoids", "parasitoid", "Roptrocerus",
           "Rhopalicus", "Dinotiscus", "Coeloides", "Spathius", "Dendrosoter", "Ecphylus", "Cheiropachus", "Tomicobia",
           "Heydenia", "predators", "Thanasimus", "Medetera", "Rhizophagus", "Temnochila", "Enoclerus", "Corticeus",
           "natural enemies", "pathogens", "virus", "Wolbachia", "endosymbiont", "bacteria", "gut microbiota"]
G_FUNG = "(fungus OR fungi OR fungal OR symbiont OR symbionts OR mycangia OR mycangium OR yeast OR yeasts OR ophiostomatoid)"
G_ENEMY = ("(pathogen OR pathogens OR parasite OR parasites OR parasitoid OR parasitoids OR predator OR predators OR "
           "nematode OR nematodes OR mite OR mites OR entomopathogenic OR microsporidia OR \"natural enemies\")")


# Author pass (2026-10-18, user: searches targeted at Platypodinae, which lag on every data type): the bibliographies
# of the taxonomists who described most Platypodinae (and many Scolytinae), surname + a bark / ambrosia beetle term.
AUTHORS = ["Schedl", "Beaver", "Browne", "Nunberg", "Jordal", "Wood", "Bright", "Atkinson", "Sittichaya", "Roberts",
           "Strohmeyer", "Hagedorn", "Eggers", "Chapuis", "Thompson", "Saha", "Maiti", "Thube", "Gomez", "Smith",
           "Cognato", "Hulcr", "Knizek", "Petrov", "Mandelshtam", "Kirkendall", "Faccoli", "Alonso-Zarazaga",
           "Lanfranco", "Johnson", "Bateman", "Storer", "Sanchez", "Equihua", "Menocal", "Mecke", "Wu", "Lin",
           "Wang", "Kambestad", "Pistone", "Gohli", "Hongsuwan", "Murayama", "Nobuchi", "Choo", "Gebhardt",
           "Kasson", "Biedermann", "Peris", "Pfeffer", "Reitter", "Hopkins", "Swaine", "Blackman", "Blandford"]
CTX_TAX = ('(Platypodinae OR Platypodidae OR Platypus OR Scolytinae OR Scolytidae OR "ambrosia beetle" OR '
           '"ambrosia beetles" OR "bark beetle" OR "bark beetles" OR "pinhole borer")')


def cap_of(kind):
    return {"topic": 50000, "group": 5000, "fungi": 5000, "enemy": 5000, "author": 5000}.get(kind, 2000)


def query_of(name, kind):
    if kind == "group":
        return f'"{name}" AND {CTX}'
    if kind in ("fungi", "enemy"):
        return f'"{name}" AND {CTX_BB}'
    if kind == "author":
        return CTX_TAX
    if kind == "gfung":
        return f'"{name}" AND {G_FUNG}'
    if kind == "genemy":
        return f'"{name}" AND {G_ENEMY}'
    return f'"{name}"'


OA_BUDGET = []


def openalex(name, since, kind="syn"):
    if OA_BUDGET:  # keyless daily budget of the IP used up (401/429): skip until a key is in the key file
        return [], "budget"
    flt = f"from_publication_date:{since}" if since else None
    cur, out = "*", []
    while cur and len(out) < cap_of(kind):
        p = dict(search=query_of(name, kind), per_page=200, cursor=cur, **({"api_key": keyed("OPENALEX_API_KEY")}
                                                                          if keyed("OPENALEX_API_KEY") else {}),
                 select="id,doi,title,publication_year,type,open_access,best_oa_location,primary_location")
        if kind == "author":
            flt = ",".join(x for x in (flt, f"raw_author_name.search:{name}") if x)
        if flt:
            p["filter"] = flt
        js = get("https://api.openalex.org/works", p)
        if "_http" in js:
            if js["_http"] in (401, 429) and not keyed("OPENALEX_API_KEY"):
                OA_BUDGET.append(1)
            return out, str(js["_http"])
        out += js.get("results", [])
        cur = js.get("meta", {}).get("next_cursor")
        if not js.get("results"):
            break
    return out, "ok"


def epmc(name, since, kind="syn"):
    q = (f'AUTH:"{name}" AND {CTX_TAX}' if kind == "author" else query_of(name, kind)) + \
        (f" AND FIRST_PDATE:[{since} TO 3000-01-01]" if since else "")
    cur, out = "*", []
    while cur and len(out) < cap_of(kind):
        js = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                 dict(query=q, format="json", resultType="lite", pageSize=1000, cursorMark=cur))
        if "_http" in js:
            return out, str(js["_http"])
        res = js.get("resultList", {}).get("result", [])
        out += res
        nxt = js.get("nextCursorMark")
        cur = nxt if nxt and nxt != cur and res else None
    return out, "ok"


def keyed(name):
    try:
        from api_keys import key
        return key(name)
    except Exception:
        return None


S2_REFUSED = []


def s2(name, since, kind="syn"):
    if kind == "author":  # no author filter in the bulk search
        return [], "ok"
    if S2_REFUSED:  # the key was refused (403 on 2026-10-15: expired?); keyless requests share a tiny pool
        return [], "key_refused"
    k = keyed("SEMANTIC_SCHOLAR_API_KEY")
    q = query_of(name, kind).replace(" AND ", " + ").replace(" OR ", " | ")
    p = dict(query=q, fields="title,year,externalIds,openAccessPdf,venue,publicationTypes")
    if since:
        p["publicationDateOrYear"] = f"{since}:"
    out, token = [], None
    while len(out) < cap_of(kind):
        if token:
            p["token"] = token
        time.sleep(1.1)
        try:
            r = S.get("https://api.semanticscholar.org/graph/v1/paper/search/bulk", params=p,
                      headers={"x-api-key": k} if k else {}, timeout=90)
        except requests.RequestException:
            return out, "failed"
        if r.status_code == 403:
            S2_REFUSED.append(1)
        if r.status_code != 200:
            return out, str(r.status_code)
        js = r.json()
        out += js.get("data") or []
        token = js.get("token")
        if not token or not js.get("data"):
            break
    return out, "ok"


def scopus(name, since, kind="syn"):
    k = keyed("SCOPUS_API_KEY")
    if not k:
        return [], "no_key"
    q = f'TITLE-ABS-KEY("{name}")' + {
        "group": " AND TITLE-ABS-KEY(beetle OR Scolytinae OR Platypodinae OR Scolytidae OR Platypodidae OR Curculionidae)",
        "fungi": f" AND TITLE-ABS-KEY{CTX_BB}", "enemy": f" AND TITLE-ABS-KEY{CTX_BB}",
        "gfung": f" AND TITLE-ABS-KEY{G_FUNG}", "genemy": f" AND TITLE-ABS-KEY{G_ENEMY}"}.get(kind, "")
    if kind == "author":
        q = f"AUTHLASTNAME({name}) AND TITLE-ABS-KEY{CTX_TAX}"
    if since:
        q += f" AND PUBYEAR > {int(since[:4]) - 1}"
    out, cur = [], "*"
    while cur and len(out) < cap_of(kind):
        try:
            r = S.get("https://api.elsevier.com/content/search/scopus", headers={"X-ELS-APIKey": k,
                                                                                "Accept": "application/json"},
                      params=dict(query=q, cursor=cur, count=200,
                                  field="dc:identifier,prism:doi,dc:title,prism:coverDate,prism:publicationName,"
                                        "openaccessFlag,subtypeDescription"), timeout=90)
        except requests.RequestException:
            return out, "failed"
        if r.status_code != 200:
            return out, str(r.status_code)
        sr = r.json().get("search-results", {})
        ent = [e for e in sr.get("entry", []) if "error" not in e]
        out += ent
        nxt = (sr.get("cursor") or {}).get("@next")
        cur = nxt if ent and nxt and nxt != cur else None
        time.sleep(0.2)
    return out, "ok"


class Cache:
    def __init__(self, name):
        self.f = OUT / f"cache/{name}.jsonl"
        self.d = {}
        if self.f.exists():
            for line in open(self.f, encoding="utf-8"):
                r = json.loads(line)
                if r.get("status") == "ok":
                    self.d[r["k"]] = r
        self.fh = open(self.f, "a", encoding="utf-8")

    def fill(self, jobs, fn, workers):
        todo = [j for j in jobs if j[0] not in self.d]
        print(self.f.name, "to search", len(todo), "| cached", len(self.d), flush=True)
        with ThreadPoolExecutor(workers) as ex:
            for i, ((k, name, since), (res, st)) in enumerate(zip(todo, ex.map(lambda j: fn(j[1], j[2], j[0].split("|")[0]),
                                                                                 todo))):
                r = dict(k=k, name=name, since=since, status=st, n=len(res), results=res, retrieved=TODAY)
                self.fh.write(json.dumps(r, ensure_ascii=False) + "\n")
                if st == "ok":
                    self.d[k] = r
                if i % 500 == 0:
                    self.fh.flush()
                    print(" ", self.f.name, i, "of", len(todo), flush=True)
        self.fh.flush()


def norm_doi(d):
    d = str(d or "").strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "http://dx.doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d


def have_dois():
    h = set()
    dl = BB / "download/doi_list.csv"
    if dl.exists():
        h |= {norm_doi(x) for x in pd.read_csv(dl).iloc[:, 0].astype(str)}
    pf = BB / "download/processed_dois.txt"
    if pf.exists():
        h |= {norm_doi(x) for x in open(pf, encoding="utf-8", errors="replace")}
    ft = BB / "download/full_text_downloads"
    got = {norm_doi(f.stem.replace("_", "/", 1)) for f in ft.iterdir()} if ft.exists() else set()
    return h, got


def main(limit=0, only=None):
    (OUT / "cache").mkdir(parents=True, exist_ok=True)
    v, d = load_names()
    v = v[~v.is_placeholder]
    valid = {r.scientificName: int(r.valid_species_id) for r in v.itertuples()}
    syn = {}
    for r in d.itertuples():
        n = str(r.describedScientificName).strip()
        if int(r.name_valid_species_id) in set(valid.values()) and n not in valid and len(n.split()) >= 2:
            syn.setdefault(n, int(r.name_valid_species_id))
    groups = sorted({x for c in ("genus", "tribe", "subfamily") for x in v[c].dropna().astype(str) if len(x) > 2})
    jobs = [(f"topic|{t}", t, "") for t in TOPICS] + [(f"group|{g}", g, "") for g in groups] + \
        [(f"syn|{n}", n, "") for n in syn] + [(f"new|{n}", n, SINCE) for n in valid] + \
        [(f"fungi|{t}", t, "") for t in FUNGI] + [(f"enemy|{t}", t, "") for t in ENEMIES] + \
        [(f"gfung|{g}", g, "") for g in sorted({x for x in v.genus.dropna().astype(str) if len(x) > 2})] + \
        [(f"genemy|{g}", g, "") for g in sorted({x for x in v.genus.dropna().astype(str) if len(x) > 2})] + \
        [(f"author|{a}", a, "") for a in AUTHORS]
    if only:
        jobs = [j for j in jobs if j[0].split("|")[0] in only]
    print("group names", len(groups), "| topic phrases", len(TOPICS), flush=True)
    if limit:
        jobs = jobs[:limit]
    print("names: synonyms / original combinations", len(syn), "(all years) | valid", len(valid), f"(since {SINCE})",
          flush=True)
    oa, ep, ss, sc = Cache("openalex"), Cache("epmc"), Cache("s2"), Cache("scopus")
    with ThreadPoolExecutor(4) as ex:  # the services in parallel; Semantic Scholar 1 request / s with the key
        fs = [ex.submit(oa.fill, jobs, openalex, 3), ex.submit(ep.fill, jobs, epmc, 3), ex.submit(ss.fill, jobs, s2, 1),
              ex.submit(sc.fill, jobs, scopus, 2)]
        [f.result() for f in fs]
    vid_of = {**{n: v_ for n, v_ in valid.items()}, **syn}
    rows = []
    for r in oa.d.values():
        for w in r["results"]:
            boa = w.get("best_oa_location") or {}
            rows.append(dict(name=r["name"], valid_species_id=vid_of.get(r["name"]) if r["k"][:4] in ("syn|", "new|") else None,
                             pass_=r["k"].split("|")[0], openalex_id=w.get("id"), doi=norm_doi(w.get("doi")), title=w.get("title"),
                             year=w.get("publication_year"), type=w.get("type"),
                             is_oa=(w.get("open_access") or {}).get("is_oa"),
                             oa_status=(w.get("open_access") or {}).get("oa_status"), oa_pdf_url=boa.get("pdf_url"),
                             oa_licence=boa.get("license"),
                             source=((w.get("primary_location") or {}).get("source") or {}).get("display_name"),
                             retrieved=r["retrieved"]))
    H1 = pd.DataFrame(rows)
    H1["year"] = pd.to_numeric(H1.year, errors="coerce").astype("Int64")
    H1.to_parquet(OUT / "hits_openalex.parquet", index=False)
    rows = []
    for r in ep.d.values():
        for w in r["results"]:
            rows.append(dict(name=r["name"], valid_species_id=vid_of.get(r["name"]) if r["k"][:4] in ("syn|", "new|") else None,
                             pass_=r["k"].split("|")[0], epmc_id=f"{w.get('source')}:{w.get('id')}", pmcid=w.get("pmcid"),
                             doi=norm_doi(w.get("doi")), title=w.get("title"), year=w.get("pubYear"),
                             is_oa=w.get("isOpenAccess") == "Y", journal=w.get("journalTitle"),
                             retrieved=r["retrieved"]))
    H2 = pd.DataFrame(rows)
    H2["year"] = pd.to_numeric(H2.year, errors="coerce").astype("Int64")
    H2.to_parquet(OUT / "hits_epmc.parquet", index=False)
    rows = []
    for r in ss.d.values():
        for w in r["results"]:
            ext = w.get("externalIds") or {}
            rows.append(dict(name=r["name"], valid_species_id=vid_of.get(r["name"]) if r["k"][:4] in ("syn|", "new|") else None,
                             pass_=r["k"].split("|")[0], doi=norm_doi(ext.get("DOI")),
                             pmcid=f"PMC{ext['PubMedCentral']}" if ext.get("PubMedCentral") else None,
                             title=w.get("title"), year=w.get("year"), is_oa=bool(w.get("openAccessPdf")),
                             oa_pdf_url=(w.get("openAccessPdf") or {}).get("url"), source=w.get("venue"),
                             retrieved=r["retrieved"]))
    H3 = pd.DataFrame(rows)
    if len(H3):
        H3["year"] = pd.to_numeric(H3.year, errors="coerce").astype("Int64")
        H3.to_parquet(OUT / "hits_s2.parquet", index=False)
    rows = []
    for r in sc.d.values():
        for w in r["results"]:
            rows.append(dict(name=r["name"], valid_species_id=vid_of.get(r["name"]) if r["k"][:4] in ("syn|", "new|") else None,
                             pass_=r["k"].split("|")[0], doi=norm_doi(w.get("prism:doi")), title=w.get("dc:title"),
                             year=str(w.get("prism:coverDate") or "")[:4], is_oa=w.get("openaccessFlag") in (True, "true"),
                             source=w.get("prism:publicationName"), retrieved=r["retrieved"]))
    H4 = pd.DataFrame(rows)
    if len(H4):
        H4["year"] = pd.to_numeric(H4.year, errors="coerce").astype("Int64")
        H4.to_parquet(OUT / "hits_scopus.parquet", index=False)
    known, got = have_dois()
    A = pd.concat([H1.assign(via="openalex"), H2.assign(via="europepmc"), H3.assign(via="semantic_scholar"),
                   H4.assign(via="scopus")], ignore_index=True)
    A["key"] = A.doi.where(A.doi.fillna("").ne(""), A.get("pmcid"))
    A = A[A.key.fillna("").ne("")]
    A["in_bb_full_text"] = A.doi.isin(got)
    A["in_2025_doi_list"] = A.doi.isin(known)
    C = A.groupby("key").agg(doi=("doi", "first"), pmcid=("pmcid", "first"), title=("title", "first"),
                             year=("year", "first"), names=("name", lambda s: json.dumps(sorted(set(s)))),
                             species=("valid_species_id", lambda s: json.dumps(sorted({int(x) for x in s if x == x}))),
                             is_oa=("is_oa", "max"), oa_pdf_url=("oa_pdf_url", "first"),
                             oa_licence=("oa_licence", "first"), via=("via", lambda s: ",".join(sorted(set(s)))),
                             in_bb_full_text=("in_bb_full_text", "max"),
                             in_2025_doi_list=("in_2025_doi_list", "max")).reset_index()
    C.to_parquet(OUT / "candidates.parquet", index=False)
    new = C[~C.in_bb_full_text.astype(bool)]
    print("hits", len(A), "| distinct works", len(C), "| already full text in BB_papers", int(C.in_bb_full_text.sum()),
          "| not in BB_papers", len(new), "(of which open access", int(new.is_oa.fillna(False).astype(bool).sum()),
          ") | not in the 2025 DOI list", int((~C.in_2025_doi_list.astype(bool)).sum()), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="", help="search only these passes, e.g. fungi,enemy,gfung,genemy")
    a = ap.parse_args()
    main(a.limit, set(a.only.split(",")) - {""})
