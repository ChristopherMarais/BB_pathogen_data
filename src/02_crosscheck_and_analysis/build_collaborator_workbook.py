"""Comprehensive review and verification workbook for collaborators (Y. Dong & J. Hulcr):
Combines:
  - Sheets 1–5: Candidate pathogen/parasite records found by the automated literature & database search (with interactive
    drop-down columns for collaborator verification);
  - Sheet 6: Secondary citation pointers (statements in primary papers that cite another work);
  - Sheet 7: The full current Figshare dataset (1,014 rows) enriched side-by-side with automated quality-control (QC) flags
    for beetle taxonomy, GBIF pathogen taxonomy, ecological relationship standardization, full-text flags, and automated
    search/extraction cross-check status;
  - Sheet 8: Automated search recall across all 276 record-contributing primary publications;
  - Sheet 9: Automated text-extraction agreement across all 895 publication x beetle x pathogen pairs;
  - Sheet 10: Pathogen surveillance gaps (introduced/regulated beetle species and countries) and rescaled decadal trends.

LAB ONLY: Sheets 1–6 quote verbatim supporting sentences from (often copyrighted) sources.
Usage: python build_review_sheet.py <INTERNAL_candidate_review_table.csv> <out.xlsx>"""
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.datavalidation import DataValidation

SEED, N_SAMPLE = 20261026, 100
src, out = sys.argv[1], sys.argv[2]
A = Path(__file__).resolve().parents[2] / "data"

def clean_df(df):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object or str(df[c].dtype).startswith("str"):
            df[c] = df[c].fillna("").astype(str).str.replace(ILLEGAL_CHARACTERS_RE, " ", regex=True)
    return df

T = clean_df(pd.read_csv(src, dtype=str).fillna(""))
B = lambda s: s.astype(str).str.lower().isin(["true", "1", "yes"])
T["focused_b"], T["new_pair_b"], T["new_beetle_b"] = B(T.focused), B(T.new_pair), B(T.new_beetle_species)
lit = T.stratum.str.contains("literature", case=False)
T["sheet"] = np.where(~lit | ~T.focused_b, "5", "4")
T.loc[lit & T.beetle_subfamily.eq("Platypodinae"), "sheet"] = "2"
T.loc[lit & T.focused_b & T.beetle.eq("Hypothenemus hampei"), "sheet"] = "3"
pool = T[lit & T.focused_b & T.sheet.isin(["3", "4"])]
sample_ids = pool.sample(n=min(N_SAMPLE, len(pool)), random_state=SEED).index
T.loc[sample_ids, "sheet"] = "1"
print("Candidate records per sheet:", T.sheet.value_counts().sort_index().to_dict(), "| total", len(T))

COLS = [
    ("Record ID", "fact_id", 16), ("Source Stratum", "stratum", 22), ("Beetle", "beetle", 24),
    ("Organism as written", "pathogen_as_written", 26), ("Organism (GBIF name)", "pathogen_resolved_gbif", 26),
    ("Category", "category", 12), ("Relation as read", "relation_read", 14), ("Year", "year", 7),
    ("Publication", "citation", 60), ("DOI / link", "link", 30), ("Page / Locator", "locator", 16),
    ("Source ID", "source_id", 24), ("Language", "language", 9),
    ("Beetle new to dataset", "new_beetle", 10), ("Pair new to dataset", "new_pair", 10),
    ("Supporting sentence (lab only)", "supporting_quote_INTERNAL_ONLY", 70),
    ("MEETS INCLUSION RULES?", None, 14), ("Beetle correct?", None, 10), ("Organism correct?", None, 10),
    ("If not included: reason", None, 26), ("Country / region (if stated)", None, 18), ("Notes", None, 30)
]
T["link"] = np.where(T.doi.ne(""), "https://doi.org/" + T.doi.str.replace("https://doi.org/", "", regex=False), T.url)
T["stratum"] = T["stratum"].str.replace("literature read by the KG", "automated literature reading", regex=False)
T["new_beetle"] = np.where(T.new_beetle_b, "yes", "")
T["new_pair"] = np.where(T.new_pair_b, "yes", "")
T["relation_read"] = T.relation_read.where(T.relation_read.ne(""), T.relation)

wb = Workbook()
ws = wb.active
ws.title = "Start here"
intro = [
    ("Bark & Ambrosia Beetle Pathogens — Collaborator Verification & Cross-Check Workbook", True),
    ("LAB ONLY: Sheets 1–6 quote verbatim sentences from sources, many of them copyrighted. Do not distribute publicly.", True),
    ("", False),
    ("1. Overview of this Workbook", True),
    ("This workbook accompanies the manuscript 'Candidate entomopathogens reported from bark and ambrosia beetles' and "
     "provides all underlying records for collaborator verification across two parts: (A) candidate pathogen and parasite "
     "records identified by the independent automated literature and database search (Sheets 1–6), and (B) row-by-row "
     "quality-control checks, search recall, extraction agreement, and surveillance gap tables for the existing 1,010-record "
     "dataset (Sheets 7–10).", False),
    ("", False),
    ("2. How the Automated Sources Were Obtained and Filtered", True),
    ("• Taxonomic query backbone: Built from 7,895 valid species and 10,861 described names (all synonyms and original "
     "combinations) from Johnson et al. (2026) and Wood & Bright (1992).", False),
    ("• Two-round bibliographic API search: Round 1 (2025) queried all 7,858 valid species names in PubMed, OpenAlex, "
     "Semantic Scholar, Crossref, Scopus, Web of Science, Springer Nature, and Dimensions (~136,000 DOIs). Round 2 (October "
     "2026) ran 13,922 structured queries per service in OpenAlex, Europe PMC, Scopus, and Semantic Scholar covering all "
     "5,654 synonyms/original combinations, 345 genus/tribe/subfamily names with beetle context terms, 28 multilingual topic "
     "phrases, 62 pathogen/parasite/natural-enemy terms and 58 fungal terms paired with beetle terms, every valid beetle "
     "genus paired with pathogen/enemy terms, and 56 specialist author surnames.", False),
    ("• Authorized full-text, archival, abstract, and database retrieval: Full texts (~32,800 documents; 466,051 pages) "
     "were obtained via Europe PMC, PMC Open Access, Unpaywall/OpenAlex OA links, Pensoft/ZooKeys XML, institutional "
     "Elsevier and Springer Nature TDM APIs, Biodiversity Heritage Library (113,503 name-indexed pages + 565 Internet "
     "Archive monographs), Plazi TreatmentBank (4,245 XML treatments), and Zobodat (313 central European journal PDFs). "
     "Abstracts (~113,000) were retrieved for remaining works via OpenAlex, Europe PMC, Crossref, Elsevier, Springer, and "
     "Scopus APIs. Structured records were also queried from Global Biotic Interactions (GloBI; 78,325 rows screened) and "
     "NCBI GenBank nuccore host/isolation qualifiers (54,289 Fungi/Nematoda/Viruses records screened; 14,205 matched).", False),
    ("• Automated text extraction & deterministic verification: Pages and abstracts naming a beetle and containing a "
     "pathogen/parasite/fungus cue (72,112 candidate passages/pages) were processed with an instruction-tuned open-weight "
     "language model (Gemma 4 26B-A4B-it, temperature 0, strict JSON schema). Claims had to pass four deterministic gates: "
     "(i) verbatim supporting sentence alignment >= 90% in the source text; (ii) exact beetle synonym/valid-name resolution; "
     "(iii) GBIF Backbone Taxonomy resolution of the partner organism; and (iv) rule-based rejection of reference-list lines, "
     "negated claims, generic group words, and plant-pathogenic/symbiotic fungi lacking explicit beetle-harm wording.", False),
    ("• Sequential inclusion filters for candidate records (R1–R6): R1 = reported as pathogen or parasite; R2 = belongs to "
     "one of the six target groups (Fungi, Bacteria, Viruses, Nematodes, Protists, Microsporidia); R3 = fungi and bacteria "
     "restricted to entomopathogenic groups or explicit beetle-harm statements; R4 = primary research sources only (reviews, "
     "catalogs, books, websites, and second-hand citations excluded); R5 = publication absent from the manual dataset; "
     "R6 (focused subset) = publication title names a bark/ambrosia beetle term or genus and supporting sentence contains "
     "no citation marker to another work.", False),
    ("", False),
    ("3. Guide to the Sheets & Order of Review", True),
    ("• Sheet '1 Random sample' (100 records, seed 20261026): CHECK FIRST. Checking all 100 rows in the yellow columns "
     "provides a stratified manual precision estimate (with 95% confidence interval) for the focused candidate records.", False),
    ("• Sheet '2 Platypodinae' (13 records): Candidate pathogen records for Platypodinae (backing the suggested Results sentence).", False),
    ("• Sheet '3 Coffee berry borer' (244 records): Focused candidate records for Hypothenemus hampei.", False),
    ("• Sheet '4 Other candidates' (499 records): All other focused primary-source candidate records.", False),
    ("• Sheet '5 Lower priority' (809 records): Candidates outside the title-focus filter plus GloBI and NCBI GenBank records.", False),
    ("• Sheet '6 Secondary pointers' (540 records): Statements in primary papers that cite another work (useful as literature pointers).", False),
    ("• Sheet '7 Current dataset + QC' (1,014 rows): Your complete Figshare dataset with side-by-side automated QC flags "
     "(beetle taxonomy & synonymy, GBIF pathogen name standardization & 'Nematodan' placeholder flag, ecological-relationship "
     "and full-text column cleanup flags, and whether each publication/record was matched in the automated cross-check).", False),
    ("• Sheet '8 Search recall (276 pubs)': Publication-by-publication comparison showing DOI/OpenAlex resolution and "
     "which automated search passes retrieved each of the 276 record-contributing publications.", False),
    ("• Sheet '9 Extraction agreement (895)': Pair-by-pair comparison showing automated text-extraction recall across all "
     "895 beetle–pathogen pairs in your dataset.", False),
    ("• Sheet '10 Surveillance & trends': Summary tables of invasive beetle species lacking pathogen records, country-level "
     "pathogen vs. introduced-beetle counts, and rescaled decadal publication trends (backing Figures 5–7 and Usage Notes).", False),
]
for i, (txt, bold) in enumerate(intro, 1):
    c = ws.cell(row=i, column=1, value=txt)
    c.font = Font(bold=bold, size=13 if i == 1 else 11)
    c.alignment = Alignment(wrap_text=True, vertical="top")
ws.column_dimensions["A"].width = 135

YELLOW = PatternFill("solid", fgColor="FFF2CC")
HEAD = PatternFill("solid", fgColor="D9E1F2")
FLAG_FILL = PatternFill("solid", fgColor="FCE4D6")

names = {
    "1": "1 Random sample",
    "2": "2 Platypodinae",
    "3": "3 Coffee berry borer",
    "4": "4 Other candidates",
    "5": "5 Lower priority",
}
reasons = ('"not called pathogen/parasite,associate only (no harm),cites another work,review/secondary source,'
           'wrong beetle,wrong organism,not one of the six groups,already in dataset,other"')

for key, title in names.items():
    D = T[T.sheet.eq(key)].sort_values(["citation", "beetle", "pathogen_as_written"])
    s = wb.create_sheet(title)
    for j, (h, _, w) in enumerate(COLS, 1):
        c = s.cell(row=1, column=j, value=h)
        c.font = Font(bold=True)
        c.fill = YELLOW if COLS[j - 1][1] is None else HEAD
        c.alignment = Alignment(wrap_text=True, vertical="top")
        s.column_dimensions[get_column_letter(j)].width = w
    for i, (_, r) in enumerate(D.iterrows(), 2):
        for j, (_, col, _) in enumerate(COLS, 1):
            if col is None:
                s.cell(row=i, column=j).fill = YELLOW
                continue
            v = r[col]
            c = s.cell(row=i, column=j, value=int(v) if col == "year" and v.isdigit() else v)
            c.alignment = Alignment(wrap_text=col in ("citation", "supporting_quote_INTERNAL_ONLY"), vertical="top")
            if col == "link" and v.startswith("http"):
                c.hyperlink = v
                c.font = Font(color="0563C1", underline="single")
            if col in ("beetle", "pathogen_as_written", "pathogen_resolved_gbif"):
                c.font = Font(italic=True)
    n = len(D) + 1
    # Reviewer columns are Q, R, S, T in the 22-column layout (17th to 20th columns)
    for col_letter, formula in (("Q", '"Yes,No,Unsure"'), ("R", '"Yes,No"'), ("S", '"Yes,No"'), ("T", reasons)):
        dv = DataValidation(type="list", formula1=formula, allow_blank=True)
        s.add_data_validation(dv)
        dv.add(f"{col_letter}2:{col_letter}{max(n, 2)}")
    s.freeze_panes = "D2"
    s.auto_filter.ref = f"A1:{get_column_letter(len(COLS))}{max(n, 1)}"


def write_table_sheet(wb, title, df, col_specs, freeze="C2", flag_cols=()):
    """Generic helper to write a formatted sheet from a DataFrame with (header, col_name, width, is_review)."""
    s = wb.create_sheet(title)
    df = clean_df(df)
    for j, (h, col, w, is_rev) in enumerate(col_specs, 1):
        c = s.cell(row=1, column=j, value=h)
        c.font = Font(bold=True)
        c.fill = YELLOW if is_rev else HEAD
        c.alignment = Alignment(wrap_text=True, vertical="top")
        s.column_dimensions[get_column_letter(j)].width = w
    for i, (_, r) in enumerate(df.iterrows(), 2):
        for j, (h, col, w, is_rev) in enumerate(col_specs, 1):
            if is_rev and (col is None or col not in r):
                s.cell(row=i, column=j).fill = YELLOW
                continue
            v = r[col] if col in r else ""
            if isinstance(v, str) and v.isdigit() and col in ("Records ID", "year", "n_records", "n_beetles", "n_pathogens"):
                v_out = int(v)
            else:
                v_out = v
            c = s.cell(row=i, column=j, value=v_out)
            c.alignment = Alignment(wrap_text=w >= 35, vertical="top")
            if is_rev:
                c.fill = YELLOW
            elif col in flag_cols and str(v).strip() != "":
                c.fill = FLAG_FILL
            if isinstance(v, str) and v.startswith("http"):
                c.hyperlink = v
                c.font = Font(color="0563C1", underline="single")
    n = len(df) + 1
    if freeze:
        s.freeze_panes = freeze
    s.auto_filter.ref = f"A1:{get_column_letter(len(col_specs))}{max(n, 1)}"
    return s


# ---------------- Sheet 6: Secondary pointers
pointers_csv = A / "out/INTERNAL_candidates_pointers.csv"
if "INTERNAL" in str(src) and Path("/home/gmarais/gmarais/PhD/BarkBeetle_KG/papers/pathogens/analysis/out/INTERNAL_candidates_pointers.csv").exists():
    pointers_csv = Path("/home/gmarais/gmarais/PhD/BarkBeetle_KG/papers/pathogens/analysis/out/INTERNAL_candidates_pointers.csv")
P_df = clean_df(pd.read_csv(pointers_csv, dtype=str).fillna(""))
p_specs = [
    ("Fact ID", "fact_id", 16, False),
    ("Beetle", "beetle", 24, False),
    ("Organism as written", "written", 25, False),
    ("Organism (GBIF)", "label", 25, False),
    ("Category", "cat", 12, False),
    ("Relation", "relation", 12, False),
    ("Citing Publication", "citation", 55, False),
    ("URL / DOI", "url", 30, False),
    ("Cited Work(s) (second_hand_of)", "second_hand_of", 30, False),
    ("Resolved Cited Reference(s)", "resolved_cited_references", 60, False),
    ("Supporting sentence (lab only)", "quote_INTERNAL", 65, False),
    ("Collaborator Notes", None, 30, True),
]
write_table_sheet(wb, "6 Secondary pointers", P_df, p_specs, freeze="C2")

# ---------------- Sheet 7: Current dataset (1,014 rows) + automated QC flags
FIGSHARE_PATH = A / "figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv" if (A / "figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv").exists() else A / "data/figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv"
FS = clean_df(pd.read_csv(FIGSHARE_PATH, encoding="utf-8-sig", dtype=str).fillna(""))
TB = clean_df(pd.read_csv(A / "out/taxonomy_beetles.csv", dtype=str).fillna(""))
TP = clean_df(pd.read_csv(A / "out/taxonomy_pathogens.csv", dtype=str).fillna(""))
RC = clean_df(pd.read_csv(A / "out/recall_by_publication.csv", dtype=str).fillna(""))
AG = clean_df(pd.read_csv(A / "out/agreement_pairs.csv", dtype=str).fillna(""))

# Map beetle taxonomy check by (ReportedBeetleSpecies, AcceptedBeetleSpeciesName)
tb_map = {}
for _, r in TB.iterrows():
    key = (r["reported"].strip(), r["accepted_in_dataset"].strip())
    flags = []
    if r["kg_status"] != "valid name in the species list":
        flags.append(f"Status: {r['kg_status']}")
    if r["authority_differs"].lower() == "true":
        flags.append(f"Authority in checklist: {r['kg_authority']}")
    if r["tribe_differs"].lower() == "true":
        flags.append(f"Tribe in checklist: {r['kg_tribe']}")
    if r["reported_vs_accepted_consistent"].lower() == "false":
        flags.append(f"Reported name resolves to: {r['reported_resolves_to']}")
    if r["same_or_near_epithet_in_list"]:
        flags.append(f"Near match: {r['same_or_near_epithet_in_list']}")
    tb_map[key] = {
        "kg_valid_name": r["kg_valid_name"],
        "kg_authority": r["kg_authority"],
        "kg_tribe": r["kg_tribe"],
        "beetle_flag": "; ".join(flags),
    }

# Map pathogen taxonomy check by pathogen_in_dataset
tp_map = {}
for _, r in TP.iterrows():
    p = r["pathogen_in_dataset"].strip()
    pflags = []
    if r["flag"]:
        pflags.append(r["flag"])
    if p in ("Nematodan", "hemolymph nematodes", "Entomopoxvirus", "similar to the ItEPV"):
        pflags.append("Descriptive/higher-rank placeholder (consider standardizing or excluding from Fig. 2 top taxa)")
    tp_map[p] = {
        "gbif_name": r["gbif_accepted_species"] or r["gbif_canonical"] or r["kg_resolver_gbif_name"],
        "gbif_rank": r["gbif_rank"],
        "pathogen_flag": "; ".join(pflags),
    }

# Map publication recall by normalized title
norm_ref = lambda s: re.sub(r"\s+", " ", str(s or "")).strip()
rc_map = {}
for _, r in RC.iterrows():
    rf = norm_ref(r["ref"])
    doi = r["doi"]
    link = f"https://doi.org/{doi}" if doi and not doi.startswith("10.5555/") else (r["oa_id"] or "")
    rc_map[rf] = {
        "pub_id": r["pub_id"],
        "doi_or_link": link,
        "found_auto_search": "Yes" if r["found_kg_search"].lower() == "true" else f"No ({r['why_missed']})",
    }

FS["kg_valid_name"] = [tb_map.get((r.strip(), a.strip()), {}).get("kg_valid_name", "") for r, a in zip(FS["ReportedBeetleSpecies"], FS["AcceptedBeetleSpeciesName"])]
FS["kg_authority"] = [tb_map.get((r.strip(), a.strip()), {}).get("kg_authority", "") for r, a in zip(FS["ReportedBeetleSpecies"], FS["AcceptedBeetleSpeciesName"])]
FS["kg_tribe"] = [tb_map.get((r.strip(), a.strip()), {}).get("kg_tribe", "") for r, a in zip(FS["ReportedBeetleSpecies"], FS["AcceptedBeetleSpeciesName"])]
FS["beetle_flag"] = [tb_map.get((r.strip(), a.strip()), {}).get("beetle_flag", "") for r, a in zip(FS["ReportedBeetleSpecies"], FS["AcceptedBeetleSpeciesName"])]

FS["gbif_pathogen_name"] = [tp_map.get(p.strip(), {}).get("gbif_name", "") for p in FS["pathogens"]]
FS["gbif_pathogen_rank"] = [tp_map.get(p.strip(), {}).get("gbif_rank", "") for p in FS["pathogens"]]
FS["pathogen_flag"] = [tp_map.get(p.strip(), {}).get("pathogen_flag", "") for p in FS["pathogens"]]

def eco_flag(val, p):
    if not str(p).strip():
        return "Excluded row (no pathogen taxon)"
    v = str(val).strip()
    if not v:
        return "Blank (411 rows have blank ecological relationship)"
    if v not in ("pathogen", "parasite", "endoparasite", "associate"):
        return f"Non-standard value ('{v}' mixes relationship with study type/notes)"
    return ""

def ft_flag(val):
    v = str(val).strip()
    if not v:
        return "Blank full-text status"
    if v in ("non full-text", "non full text"):
        return f"Spelling variant ('{v}')"
    return ""

FS["eco_flag"] = [eco_flag(e, p) for e, p in zip(FS["ecological relationship"], FS["pathogens"])]
FS["ft_flag"] = [ft_flag(f) for f in FS["full text"]]
FS["pub_doi_link"] = [rc_map.get(norm_ref(t), {}).get("doi_or_link", "") for t in FS["title"]]
FS["found_auto_search"] = [rc_map.get(norm_ref(t), {}).get("found_auto_search", "") for t in FS["title"]]

fs_specs = [
    ("Records ID", "Records ID", 10, False),
    ("Reported Beetle", "ReportedBeetleSpecies", 22, False),
    ("Accepted Beetle (Dataset)", "AcceptedBeetleSpeciesName", 24, False),
    ("Authorship (Dataset)", "authorship", 20, False),
    ("Subfamily", "subfamily", 12, False),
    ("Tribe (Dataset)", "tribe", 14, False),
    ("Genus", "genus", 14, False),
    ("Checklist Valid Name", "kg_valid_name", 24, False),
    ("Checklist Authority", "kg_authority", 20, False),
    ("Checklist Tribe", "kg_tribe", 14, False),
    ("QC: Beetle Taxonomy Flag", "beetle_flag", 35, False),
    ("Category", "categories", 13, False),
    ("Pathogen (Dataset)", "pathogens", 26, False),
    ("GBIF Resolved Name", "gbif_pathogen_name", 26, False),
    ("GBIF Rank", "gbif_pathogen_rank", 12, False),
    ("QC: Pathogen Taxonomy Flag", "pathogen_flag", 38, False),
    ("Country / Region", "country or region", 18, False),
    ("Ecological Relationship", "ecological relationship", 22, False),
    ("QC: Ecological Rel. Flag", "eco_flag", 36, False),
    ("Validation Type", "validation type", 20, False),
    ("Experimental Conditions", "experimental conditions", 18, False),
    ("Identification Method", "identification method", 18, False),
    ("Infection Site", "infection site", 18, False),
    ("Organism Source", "organism source", 18, False),
    ("Source Citation (Short)", "source", 22, False),
    ("Year", "year", 8, False),
    ("Full Text Status", "full text", 16, False),
    ("QC: Full Text Flag", "ft_flag", 24, False),
    ("Full Reference", "title", 55, False),
    ("Resolved DOI / Link", "pub_doi_link", 30, False),
    ("Retrieved by Automated Search?", "found_auto_search", 28, False),
    ("Collaborator Check / Action", None, 30, True),
]
write_table_sheet(
    wb, "7 Current dataset + QC", FS, fs_specs, freeze="D2",
    flag_cols=("beetle_flag", "pathogen_flag", "eco_flag", "ft_flag")
)

# ---------------- Sheet 8: Search recall (276 publications)
RC["doi_link"] = [f"https://doi.org/{d}" if d and not d.startswith("10.5555/") else (o or "") for d, o in zip(RC["doi"], RC["oa_id"])]
rc_specs = [
    ("Pub ID", "pub_id", 9, False),
    ("Short Citation", "source", 22, False),
    ("Dataset Record IDs", "record_ids", 18, False),
    ("Year", "year", 8, False),
    ("Categories", "categories", 18, False),
    ("Records in Dataset", "n_records", 12, False),
    ("Beetle Species", "n_beetles", 12, False),
    ("Pathogen Taxa", "n_pathogens", 12, False),
    ("Publication Type", "pub_type", 18, False),
    ("Language", "language_guess", 16, False),
    ("Full-Text Status (Dataset)", "full_text", 18, False),
    ("Resolved DOI / OpenAlex Link", "doi_link", 30, False),
    ("Retrieved by Automated Search?", "found_kg_search", 16, False),
    ("Found in 2025 Valid-Name Search?", "found_y2025", 16, False),
    ("Found in 2026 Synonym/Topic Search?", "found_s2026_names_topics", 16, False),
    ("Found in 2026 Enemy/Fungus Search?", "found_s2026_enemy_fungus", 16, False),
    ("Match Method", "match_quality", 16, False),
    ("Reason Missed (if not retrieved)", "why_missed", 38, False),
    ("Matched Automated-Search Title", "best_title", 45, False),
    ("Full Reference", "ref", 60, False),
    ("Collaborator Notes", None, 28, True),
]
write_table_sheet(wb, "8 Search recall (276 pubs)", RC, rc_specs, freeze="C2", flag_cols=("why_missed",))

# ---------------- Sheet 9: Extraction agreement (895 pairs)
AG["side_label"] = AG["side"].map({"her": "Dataset pair (Recall)", "kg": "Automated extraction (Precision)"}).fillna(AG["side"])
AG["pred_or_rel"] = [gp if s == "her" else (f"{p} ({r})" if p and r else (p or r)) for s, gp, p, r in zip(AG["side"], AG["graph_predicates"], AG.get("predicate", [""]*len(AG)), AG.get("relation", [""]*len(AG)))]
AG["doi_link"] = AG["doi_or_link"]
ag_specs = [
    ("Pub ID", "pub_id", 9, False),
    ("Short Citation", "source", 22, False),
    ("Dataset Record ID(s)", "dataset_record_ids", 18, False),
    ("Perspective / Evaluation", "side_label", 24, False),
    ("Beetle Species", "beetle", 24, False),
    ("Pathogen / Associate", "pathogen", 28, False),
    ("Category", "category", 14, False),
    ("Text Processed", "kg_text", 16, False),
    ("Read via", "kg_read_via", 16, False),
    ("Reader Match (Dataset Pairs)", "reader", 18, False),
    ("Graph Match (Dataset Pairs)", "graph", 18, False),
    ("Extraction Filter Reason", "rejected_reason", 32, False),
    ("Precision Status (Automated Pairs)", "precision_status", 36, False),
    ("Extracted Relation / Predicate", "pred_or_rel", 24, False),
    ("Resolved DOI / Link", "doi_link", 30, False),
    ("Full Reference", "full_reference", 60, False),
    ("Source ID", "kg_source_ids", 24, False),
    ("Collaborator Notes", None, 28, True),
]
write_table_sheet(wb, "9 Extraction agreement (895)", AG, ag_specs, freeze="D2")

# ---------------- Sheet 10: Surveillance & trends
s10 = wb.create_sheet("10 Surveillance & trends")
SV = clean_df(pd.read_csv(A / "out/surveillance_top.csv", dtype=str).fillna(""))
GC = clean_df(pd.read_csv(A / "out/geography_countries.csv", dtype=str).fillna(""))
TR = clean_df(pd.read_csv(A / "out/trend_rescaled_decade.csv", dtype=str).fillna(""))

row_cursor = 1
def write_subtable(ws, start_row, title, subtitle, df, cols):
    c = ws.cell(row=start_row, column=1, value=title)
    c.font = Font(bold=True, size=12)
    c2 = ws.cell(row=start_row + 1, column=1, value=subtitle)
    c2.font = Font(italic=True, size=10)
    hdr_row = start_row + 2
    for j, (h, col, w) in enumerate(cols, 1):
        cell = ws.cell(row=hdr_row, column=j, value=h)
        cell.font = Font(bold=True)
        cell.fill = HEAD
        cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.column_dimensions[get_column_letter(j)].width = max(ws.column_dimensions[get_column_letter(j)].width or 0, w)
    for i, (_, r) in enumerate(df.iterrows(), hdr_row + 1):
        for j, (_, col, _) in enumerate(cols, 1):
            v = r[col] if col in r else ""
            if isinstance(v, str) and re.match(r"^-?\d+$", v):
                v = int(v)
            elif isinstance(v, str) and re.match(r"^-?\d+\.\d+$", v):
                v = round(float(v), 2)
            cell = ws.cell(row=i, column=j, value=v)
            cell.alignment = Alignment(vertical="top")
    return hdr_row + len(df) + 3

sv_cols = [
    ("Beetle Species", "name", 26),
    ("Subfamily", "subfamily", 14),
    ("Tribe", "tribe", 16),
    ("Introduced Countries", "introduced_countries", 16),
    ("EPPO Listed?", "eppo_listed", 12),
    ("Top-50 Forecast Invader?", "next_invader_top50", 14),
    ("Total Literature Works", "study_effort_works", 16),
    ("Records in Manual Dataset", "her_records", 16),
    ("Candidate Sources in Automated Search", "kg_pathogen_sources_not_in_her_data", 20),
    ("Dataset Record IDs", "dataset_record_ids", 18),
    ("Dataset References", "dataset_references", 55),
    ("Candidate Pathogen References (Automated Search)", "candidate_references", 55),
    ("Introduced-Range References", "distribution_references", 55),
]
row_cursor = write_subtable(
    s10, row_cursor,
    "Table 10A. Top 40 Invasive / Regulated Bark & Ambrosia Beetle Species and Their Pathogen Coverage",
    "Ranked by number of countries where introduced (from global occurrence records) and overall study effort.",
    SV, sv_cols
)

tr_cols = [
    ("Decade", "period", 14),
    ("Dataset Publications", "her_publications", 16),
    ("OpenAlex Bark/Ambrosia Beetle Works", "openalex_bark_beetle_works", 20),
    ("Automated Corpus Works", "kg_reading_corpus_works", 18),
    ("Dataset Pubs per 1,000 OpenAlex Works", "per1000_openalex", 20),
    ("95% CI Low (OpenAlex)", "per1000_openalex_lo", 16),
    ("95% CI High (OpenAlex)", "per1000_openalex_hi", 16),
    ("Dataset Pubs per 1,000 Corpus Works", "per1000_kg", 20),
    ("All Primary Pathogen Sources (Dataset + Automated)", "kg_primary_pathogen_sources", 22),
]
row_cursor = write_subtable(
    s10, row_cursor,
    "Table 10B. Decadal Publication Trend Rescaled by Total Bark & Ambrosia Beetle Literature Volume (Figure 6)",
    "Compares raw dataset publication counts against per-1,000 rates in OpenAlex and the automated reading corpus.",
    TR, tr_cols
)

GC["her_rec_int"] = pd.to_numeric(GC["her_records"], errors="coerce").fillna(0).astype(int)
GC["intro_int"] = pd.to_numeric(GC["kg_introduced_species"], errors="coerce").fillna(0).astype(int)
for c in ("her_records", "her_beetle_species", "her_publications", "kg_species", "kg_introduced_species"):
    GC[c] = pd.to_numeric(GC[c], errors="coerce").fillna(0).astype(int).astype(str)
gc_sub = GC[(GC["her_rec_int"] > 0) | (GC["intro_int"] >= 5)].copy()
gc_sub = gc_sub.sort_values(["her_rec_int", "intro_int"], ascending=[False, False])
gc_cols = [
    ("ISO2", "iso2", 10),
    ("Country", "country", 26),
    ("Dataset Pathogen Records", "her_records", 16),
    ("Dataset Beetle Species", "her_beetle_species", 16),
    ("Dataset Publications", "her_publications", 16),
    ("Total Known Beetle Species", "kg_species", 18),
    ("Introduced Beetle Species", "kg_introduced_species", 18),
    ("Dataset Record IDs", "dataset_record_ids", 22),
    ("Dataset References", "dataset_references", 60),
]
row_cursor = write_subtable(
    s10, row_cursor,
    "Table 10C. Country-Level Pathogen Records vs. Native and Introduced Beetle Diversity",
    "Shows all countries with >= 1 dataset pathogen record or >= 5 introduced bark/ambrosia beetle species.",
    gc_sub, gc_cols
)

wb.save(out)
print("Successfully wrote", out)
