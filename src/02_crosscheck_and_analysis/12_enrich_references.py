"""Step 12: Enrich all output CSV tables with complete, traceable bibliographic references
(full citation, DOI, clickable URL, page/locator, and dataset record IDs) so collaborators
can immediately locate the source publication and verify every data claim.
"""
import csv
from pathlib import Path
import re
import json
import urllib.parse

import duckdb
import pycountry

A = Path("/home/gmarais/gmarais/PhD/BarkBeetle_KG/papers/pathogens/analysis")
OUT = A / "out"
REPO_OUT = Path("/home/gmarais/gmarais/PhD/Pathogen_Paper_Download/1_code_github_repo/data/out")
FS_PATH = Path("/home/gmarais/gmarais/PhD/Pathogen_Paper_Download/1_code_github_repo/data/figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv")
DB_PATH = Path("/home/gmarais/gmarais/PhD/BarkBeetle_KG/species_kg/snapshots/kg_v2.22_2026-10-22/species_kg/kg/kg.duckdb")

con = duckdb.connect(str(DB_PATH), read_only=True)
src_lookup = {}
for sid, cit, doi, yr in con.execute("SELECT source_id, citation, doi, year FROM source_meta").fetchall():
    url = f"https://doi.org/{doi}" if doi else ""
    src_lookup[sid] = {"citation": cit or "", "doi": doi or "", "url": url, "year": str(yr or "")}

for sid, cit, url in con.execute("SELECT DISTINCT source_id, citation, url FROM edges WHERE source_id IS NOT NULL").fetchall():
    info = src_lookup.setdefault(sid, {"citation": "", "doi": "", "url": "", "year": ""})
    if not info["citation"] and cit:
        info["citation"] = cit
    if not info["url"] and url:
        info["url"] = url

author_year_lookup = {}
for sid, info in src_lookup.items():
    cit = info["citation"]
    yr = info["year"]
    doi = info["doi"]
    url = info["url"]
    if cit and yr:
        m = re.match(r"^\s*(?:[A-Z]\.\s*)*([A-Za-zÀ-ÿ\-]+)", cit)
        if m:
            key = f"{m.group(1).lower()}|{yr}"
            suffix = f" (https://doi.org/{doi})" if doi else (f" ({url})" if url else "")
            author_year_lookup.setdefault(key, []).append(f"{cit}{suffix}")

place_refs_by_sp = {}
for sp_id, sj in con.execute("""
    SELECT taxon_id, supporting_json
    FROM places
    WHERE taxon_rank = 'species' AND predicate = 'occurs_in'
      AND statuses LIKE '%introduced%' AND level IN ('medium', 'high')
""").fetchall():
    try:
        items = json.loads(sj or "[]")
        for it in items:
            sid = it.get("source_id") or ""
            cit = it.get("citation") or src_lookup.get(sid, {}).get("citation", "")
            url = it.get("url") or src_lookup.get(sid, {}).get("url", "")
            ref_str = (cit or sid) + (f" [{url}]" if url and url not in (cit or "") else "")
            if ref_str:
                place_refs_by_sp.setdefault(sp_id, set()).add(ref_str)
    except Exception:
        pass

with open(FS_PATH, encoding="utf-8-sig") as f:
    fs_rows = list(csv.DictReader(f))

with open(OUT / "recall_by_publication.csv") as f:
    rc_rows = list(csv.DictReader(f))

pub_by_id = {r["pub_id"]: r for r in rc_rows}
norm = lambda s: re.sub(r"\s+", " ", str(s or "")).strip()
pub_by_ref = {norm(r["ref"]): r for r in rc_rows}

for r in rc_rows:
    fa = (r.get("first_author") or "").lower().strip()
    yr = (r.get("year") or "").strip()
    if fa and yr:
        doi = r.get("doi", "")
        oa = r.get("oa_id", "") or ""
        link = f"https://doi.org/{doi}" if doi and not doi.startswith("10.5555/") else oa
        pid = r["pub_id"]
        src_s = r["source"]
        ref_s = r["ref"]
        suffix = f" ({link})" if link else ""
        author_year_lookup.setdefault(f"{fa}|{yr}", []).insert(0, f"[In Manual Dataset: {pid}] {src_s}: {ref_s}{suffix}")

rec_by_pub_bp = {}
rec_by_pub = {}
rec_by_beetle = {}
rec_by_reported_acc = {}
rec_by_pathogen = {}
rec_by_country = {}

for r in fs_rows:
    rid = r["Records ID"].strip()
    ref_full = norm(r["title"])
    short_src = r["source"].strip()
    yr = r["year"].strip()
    pub = pub_by_ref.get(ref_full, {})
    pid = pub.get("pub_id", "")
    doi = pub.get("doi", "")
    oa = pub.get("oa_id", "") or ""
    link = f"https://doi.org/{doi}" if doi and not doi.startswith("10.5555/") else oa
    suffix = f" ({link})" if link else ""
    cite_tag = f"[Rec #{rid}, {pid}] {short_src}: {ref_full}{suffix}"
    short_tag = f"[Rec #{rid}] {short_src}{suffix}"

    b = r["AcceptedBeetleSpeciesName"].strip()
    rep_b = r["ReportedBeetleSpecies"].strip()
    p = r["pathogens"].strip()
    c = r["country or region"].strip()

    if pid:
        rec_by_pub.setdefault(pid, []).append(rid)
        rec_by_pub_bp.setdefault((pid, b, p), []).append(rid)
    if b:
        rec_by_beetle.setdefault(b, []).append((rid, short_tag, cite_tag))
    rec_by_reported_acc.setdefault((rep_b, b), []).append((rid, short_tag, cite_tag))
    if p:
        rec_by_pathogen.setdefault(p, []).append((rid, short_tag, cite_tag))
    if c:
        rec_by_country.setdefault(c, []).append((rid, short_tag, cite_tag))

# 1. Enrich INTERNAL_candidate_review_table.csv and candidate_records_public.csv
with open(OUT / "INTERNAL_candidate_review_table.csv") as f:
    cand_rows = list(csv.DictReader(f))

cand_by_beetle = {}
for r in cand_rows:
    if r["source_id"] == "bhl_item:105679" and "Wood" not in r["citation"]:
        r["citation"] = (
            "Wood, S.L. & Bright, D.E. (1987). A catalog of Scolytidae and Platypodidae (Coleoptera), "
            "Part 1: Bibliography. Great Basin Naturalist Memoirs, no. 11 (Brigham Young University; BHL item 105679)"
        )
    if r["stratum"] == "other database (GloBI)" and "targetTaxon=" not in r["url"] and r["beetle"]:
        r["url"] = r["url"] + "&targetTaxon=" + urllib.parse.quote(r["beetle"])
    if not r["locator"]:
        if r["stratum"] == "GenBank isolate record":
            r["locator"] = r["source_id"].replace("ncbi:nuccore:", "GenBank accession ")
        elif r["stratum"] == "other database (GloBI)":
            r["locator"] = r["source_id"]
    b = r["beetle"]
    doi = r["doi"]
    link = f"https://doi.org/{doi}" if doi else r["url"]
    loc_str = f", {r['locator']}" if r["locator"] else ""
    c_entry = f"{r['citation']}{loc_str} [{link}]"
    cand_by_beetle.setdefault(b, {})[r["source_id"]] = c_entry

with open(OUT / "INTERNAL_candidate_review_table.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(cand_rows[0].keys()))
    w.writeheader()
    w.writerows(cand_rows)

pub_cand_rows = []
for r in cand_rows:
    rc = dict(r)
    rc["INTERNAL"] = "public metadata version (verbatim quotes omitted)"
    rc["supporting_quote_INTERNAL_ONLY"] = ""
    pub_cand_rows.append(rc)

for dest in [OUT / "candidate_records_public.csv", REPO_OUT / "candidate_records_public.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(pub_cand_rows[0].keys()))
        w.writeheader()
        w.writerows(pub_cand_rows)

# 2. Enrich INTERNAL_candidates_pointers.csv with resolved cited references
with open(OUT / "INTERNAL_candidates_pointers.csv") as f:
    ptr_rows = list(csv.DictReader(f))
for r in ptr_rows:
    raw_sh = r.get("second_hand_of", "")
    keys = re.findall(r"[a-z0-9_\-]+\|\d{4}", raw_sh.lower())
    resolved = []
    for k in keys:
        matches = author_year_lookup.get(k, [])
        if matches:
            resolved.append(f"{k} -> {matches[0]}")
        else:
            resolved.append(k)
    r["resolved_cited_references"] = " ; ".join(resolved) if resolved else raw_sh

ptr_cols = list(ptr_rows[0].keys())
for dest in [OUT / "INTERNAL_candidates_pointers.csv", REPO_OUT / "INTERNAL_candidates_pointers.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=ptr_cols)
        w.writeheader()
        for r in ptr_rows:
            row_out = dict(r)
            if dest == REPO_OUT / "INTERNAL_candidates_pointers.csv":
                row_out["quote_INTERNAL"] = ""
            w.writerow(row_out)

# 3. Enrich agreement_by_publication.csv and agreement_pairs.csv
with open(OUT / "agreement_by_publication.csv") as f:
    ag_pub = list(csv.DictReader(f))

ag_pub_by_id = {}
for r in ag_pub:
    pid = r["pub_id"]
    pub = pub_by_id.get(pid, {})
    doi = pub.get("doi", "")
    oa = pub.get("oa_id", "") or ""
    link = f"https://doi.org/{doi}" if doi and not doi.startswith("10.5555/") else oa
    r["dataset_record_ids"] = pub.get("record_ids", "")
    r["full_reference"] = pub.get("ref", "")
    r["doi"] = doi if not doi.startswith("10.5555/") else ""
    r["doi_or_link"] = link
    sids = [s.strip() for s in (r.get("kg_source_ids") or "").split(";") if s.strip()]
    r["matched_automated_citations"] = " ; ".join(src_lookup.get(s, {}).get("citation", s) for s in sids)
    ag_pub_by_id[pid] = r

for dest in [OUT / "agreement_by_publication.csv", REPO_OUT / "agreement_by_publication.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(ag_pub[0].keys()))
        w.writeheader()
        w.writerows(ag_pub)

with open(OUT / "agreement_pairs.csv") as f:
    ag_pairs = list(csv.DictReader(f))

for r in ag_pairs:
    pid = r["pub_id"]
    pub = pub_by_id.get(pid, {})
    agp = ag_pub_by_id.get(pid, {})
    doi = pub.get("doi", "")
    oa = pub.get("oa_id", "") or ""
    link = f"https://doi.org/{doi}" if doi and not doi.startswith("10.5555/") else oa
    rids = rec_by_pub_bp.get((pid, r["beetle"], r["pathogen"]), rec_by_pub.get(pid, []))
    r["dataset_record_ids"] = ", ".join(rids)
    r["full_reference"] = pub.get("ref", "")
    r["doi"] = doi if not doi.startswith("10.5555/") else ""
    r["doi_or_link"] = link
    r["kg_source_ids"] = agp.get("kg_source_ids", "")

for dest in [OUT / "agreement_pairs.csv", REPO_OUT / "agreement_pairs.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(ag_pairs[0].keys()))
        w.writeheader()
        w.writerows(ag_pairs)

# 4. Enrich kg_pathogen_sources_primary.csv (803 rows)
with open(OUT / "kg_pathogen_sources_primary.csv") as f:
    kg_prim = list(csv.DictReader(f))

for r in kg_prim:
    sid = r["source_id"]
    sinfo = src_lookup.get(sid, {})
    cit = sinfo.get("citation", "")
    if sid == "bhl_item:105679" and "Wood" not in cit:
        cit = (
            "Wood, S.L. & Bright, D.E. (1987). A catalog of Scolytidae and Platypodidae (Coleoptera), "
            "Part 1: Bibliography. Great Basin Naturalist Memoirs, no. 11 (BHL item 105679)"
        )
    doi = sinfo.get("doi", "")
    url = f"https://doi.org/{doi}" if doi else sinfo.get("url", "")
    r["citation"] = cit
    r["doi"] = doi
    r["url"] = url
    hpub = r.get("her_pub", "").strip()
    r["dataset_full_reference"] = pub_by_id.get(hpub, {}).get("ref", "") if hpub else ""

for dest in [OUT / "kg_pathogen_sources_primary.csv", REPO_OUT / "kg_pathogen_sources_primary.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(kg_prim[0].keys()))
        w.writeheader()
        w.writerows(kg_prim)

# 5. Enrich surveillance_species.csv, surveillance_top.csv, surveillance_ranked_short.csv, and geography_countries.csv
with open(OUT / "surveillance_species.csv") as f:
    sv_sp = list(csv.DictReader(f))

sp_ref_map = {}
for r in sv_sp:
    b = r["name"]
    sp_id = r["species_id"]
    d_recs = rec_by_beetle.get(b, [])
    seen_c = []
    for _, _, ct in d_recs:
        if ct not in seen_c:
            seen_c.append(ct)
    r["dataset_record_ids"] = ", ".join(rid for rid, _, _ in d_recs)
    r["dataset_references"] = " ; ".join(seen_c)
    c_dict = cand_by_beetle.get(b, {})
    r["candidate_references"] = " ; ".join(c_dict.values())
    p_refs = sorted(place_refs_by_sp.get(sp_id, set()))
    r["distribution_references"] = " ; ".join(p_refs[:15]) + (f" (+{len(p_refs)-15} more)" if len(p_refs) > 15 else "")
    sp_ref_map[b] = r

for dest in [OUT / "surveillance_species.csv", REPO_OUT / "surveillance_species.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(sv_sp[0].keys()))
        w.writeheader()
        w.writerows(sv_sp)

for fn in ["surveillance_top.csv", "surveillance_ranked_short.csv"]:
    with open(OUT / fn) as f:
        srows = list(csv.DictReader(f))
    for r in srows:
        b = r["name"]
        info = sp_ref_map.get(b, {})
        r["her_records"] = info.get("her_records", r.get("her_records", "0"))
        r["her_publications"] = info.get("her_publications", "0")
        r["kg_pathogen_sources_not_in_her_data"] = info.get("kg_pathogen_sources_not_in_her_data", r.get("kg_pathogen_sources_not_in_her_data", "0"))
        r["dataset_record_ids"] = info.get("dataset_record_ids", "")
        r["dataset_references"] = info.get("dataset_references", "")
        r["candidate_references"] = info.get("candidate_references", "")
        r["distribution_references"] = info.get("distribution_references", "")
    for dest in [OUT / fn, REPO_OUT / fn]:
        with open(dest, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(srows[0].keys()))
            w.writeheader()
            w.writerows(srows)

alias = {
    "usa": "US", "united states": "US", "u.s.a.": "US", "us": "US",
    "uk": "GB", "united kingdom": "GB", "england": "GB", "great britain": "GB", "scotland": "GB",
    "czech republic": "CZ", "czechia": "CZ", "czechoslovakia": "CZ",
    "russia": "RU", "ussr": "RU", "soviet union": "RU",
    "south korea": "KR", "korea": "KR", "republic of korea": "KR",
    "turkey": "TR", "türkiye": "TR",
    "ivory coast": "CI", "côte d'ivoire": "CI", "côte d’ivoire": "CI",
    "british columbia": "CA",
    "méxico": "MX", "mexico": "MX",
}
def to_iso2(name):
    s = str(name or "").strip()
    if not s or s.lower() in ("<na>", "nan", "none", "not stated", "unknown"):
        return ""
    low = s.lower()
    if low in alias:
        return alias[low]
    for fn_lookup in (lambda x: pycountry.countries.get(name=x), lambda x: pycountry.countries.search_fuzzy(x)[0]):
        try:
            return fn_lookup(s).alpha_2
        except Exception:
            pass
    return ""

iso_recs = {}
for c_raw, rlist in rec_by_country.items():
    iso = to_iso2(c_raw)
    if iso:
        iso_recs.setdefault(iso, []).extend(rlist)

with open(OUT / "geography_countries.csv") as f:
    gc_rows = list(csv.DictReader(f))
for r in gc_rows:
    iso = r["iso2"]
    rlist = iso_recs.get(iso, [])
    r["dataset_record_ids"] = ", ".join(rid for rid, _, _ in rlist)
    seen_c = []
    for _, _, ct in rlist:
        if ct not in seen_c:
            seen_c.append(ct)
    r["dataset_references"] = " ; ".join(seen_c)

for dest in [OUT / "geography_countries.csv", REPO_OUT / "geography_countries.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(gc_rows[0].keys()))
        w.writeheader()
        w.writerows(gc_rows)

# 6. Enrich taxonomy_beetles.csv and taxonomy_pathogens.csv
with open(OUT / "taxonomy_beetles.csv") as f:
    tb_rows = list(csv.DictReader(f))
for r in tb_rows:
    key = (r["reported"].strip(), r["accepted_in_dataset"].strip())
    rlist = rec_by_reported_acc.get(key, [])
    r["dataset_record_ids"] = ", ".join(rid for rid, _, _ in rlist)
    seen_c = []
    for _, _, ct in rlist:
        if ct not in seen_c:
            seen_c.append(ct)
    r["dataset_references"] = " ; ".join(seen_c)
for dest in [OUT / "taxonomy_beetles.csv", REPO_OUT / "taxonomy_beetles.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(tb_rows[0].keys()))
        w.writeheader()
        w.writerows(tb_rows)

with open(OUT / "taxonomy_pathogens.csv") as f:
    tp_rows = list(csv.DictReader(f))
for r in tp_rows:
    p = r["pathogen_in_dataset"].strip()
    rlist = rec_by_pathogen.get(p, [])
    r["dataset_record_ids"] = ", ".join(rid for rid, _, _ in rlist)
    seen_c = []
    for _, _, ct in rlist:
        if ct not in seen_c:
            seen_c.append(ct)
    r["dataset_references"] = " ; ".join(seen_c)
for dest in [OUT / "taxonomy_pathogens.csv", REPO_OUT / "taxonomy_pathogens.csv"]:
    with open(dest, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(tp_rows[0].keys()))
        w.writeheader()
        w.writerows(tp_rows)

print("Successfully enriched all CSV tables with complete source references, DOIs/URLs, locators, and dataset record IDs!")
