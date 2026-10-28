# Automated Literature & Database Cross-Check for Bark and Ambrosia Beetle Pathogens

---

## Overview

This repository contains the pathogen-specific code and verification outputs used to:
1. **Retrieve and extract pathogen and parasite records** of bark and ambrosia beetles (*Curculionidae: Scolytinae* and *Platypodinae*) from bibliographic APIs, historical taxonomic archives, and structured biological databases (`src/01_automated_source_collection/`).
2. **Independently cross-check and analyze** the manually curated dataset of 1,010 pathogen/parasite records from 280 publications (`src/02_crosscheck_and_analysis/`), evaluating:
   - **Search recall** across the 276 record-contributing primary publications by decade, publication type, language, and DOI availability (**Figure 5A, B**);
   - **Automated text-extraction agreement** across 895 `publication × beetle × pathogen` pairs from full texts vs. abstracts (**Figure 5C**);
   - **Candidate pathogen and parasite records** from sources absent from the manual dataset, filtered sequentially by the manual dataset's inclusion criteria `R1`–`R6` (**Figure 7**);
   - **Literature-rescaled publication trends** (1914–2026) normalized per 1,000 bark and ambrosia beetle works in OpenAlex and the automated reading corpus (**Figure 6**);
   - **Taxonomic, geographic, and invasive-host surveillance gaps** across 7,895 valid beetle species (including 186 globally introduced species and EPPO-regulated pests).

---

## Methodology Summary

### Stage 1: Automated Literature & Database Retrieval (`src/01_automated_source_collection/`)
1. **Taxonomic Backbone & Two-Round API Search (`01_literature_api_search.py`)**:
   - Queries were constructed from all **7,895 valid species** and **10,861 described names** (historical synonyms and original genus combinations) compiled from Johnson et al. (2026) and Wood & Bright (1992).
   - **Round 1 (2025)**: Queried 7,858 valid species names across PubMed, OpenAlex, Semantic Scholar, Crossref, Scopus, Web of Science, Springer Nature, and Dimensions (~136,000 unique DOIs).
   - **Round 2 (October 2026)**: Executed 13,922 structured queries per service across OpenAlex, Europe PMC, Scopus (`TITLE-ABS-KEY`), and Semantic Scholar covering:
     - All 5,654 synonyms and original combinations (all years);
     - All 7,895 valid species names (2025–2026);
     - 345 genus-, subtribe-, tribe-, and subfamily-level names paired with beetle context terms;
     - 28 multilingual topic phrases (e.g., `"bark beetle"`, `"ambrosia beetle"`, `"pinhole borer"`, `"shot-hole borer"`, `"coffee berry borer"`, `"Borkenkäfer"`, `"scolyte"`, `"escarabajo de corteza"`);
     - 62 pathogen, parasite, and natural-enemy terms (including *Beauveria*, *Metarhizium*, *Isaria*, *Cordyceps*, *Lecanicillium*, *Nosema*, *Unikaryon*, *Chytridiopsis*, *Canningia*, *Mattesia*, *Menzbieria*, *Contortylenchus*, *Parasitorhabditis*, *Parasitylenchus*, *Sphaerularioides*, *Ektaphelenchus*, *Bursaphelenchus*, *Cryptaphelenchus*, *Wolbachia*, `"entomopathogenic fungi"`, `"microsporidia"`, `"gregarines"`, `"nematodes"`, `"virus"`, `"bacteria"`) and 58 fungal terms combined with bark/ambrosia beetle terms, plus every valid beetle genus paired with pathogen/enemy terms;
     - 56 specialist taxonomist surnames paired with beetle terms.
   - **Full Texts, Archives, and Abstracts**: Full texts (~32,800 documents; 466,051 pages) were obtained exclusively via authorized programmatic channels (Europe PMC, PubMed Central Open Access, Unpaywall/OpenAlex OA endpoints, Pensoft/ZooKeys XML, and Elsevier/Springer Nature TDM APIs), supplemented by the Biodiversity Heritage Library (113,503 name-indexed pages + 565 Internet Archive monographs), Plazi TreatmentBank (4,245 XML treatments), and Zobodat (313 central European journal PDFs). Abstracts (~113,000) were retrieved for remaining works via OpenAlex, Europe PMC, Crossref, Elsevier, Springer, and Scopus metadata APIs.

2. **Structured Biological Repositories (`02_structured_database_queries.py`)**:
   - **Global Biotic Interactions (GloBI)**: Filtered 78,325 interaction rows mentioning bark/ambrosia beetle taxa for exact beetle name resolution, pathogen/parasite interaction verbs (`pathogenOf`, `hasPathogen`, `parasiteOf`, `hasParasite`), and primary publication references (excluding unverified third-party LLM extractions and inferred metawebs).
   - **NCBI GenBank (`nuccore`)**: Screened 54,289 `Fungi`, `Nematoda`, and `Viruses` sequence records via Entrez E-utilities across all valid and historical beetle genus/subfamily names (14,205 matched in `/host`, `/isolation_source`, or title qualifiers), extracting linked primary publications (`GBReference` DOI/PMID).

3. **Aho–Corasick Indexing, Language Model Extraction, and Rule Verification (`03_text_indexing_and_llm_extraction.py`)**:
   - Every full-text page and abstract naming a bark or ambrosia beetle and containing a multilingual pathogen/parasite/fungus cue (72,112 candidate passages/pages) was processed with `Gemma 4 26B-A4B-it` (`vLLM`, temperature `0`, strict JSON schema).
   - Extracted records had to pass four deterministic verification gates:
     1. **Verbatim quote alignment** in the source text (`rapidfuzz` partial ratio >= 90%);
     2. **Beetle taxonomic resolution** to a valid species or genus via exact synonym/original-combination lookup or page-scoped genus abbreviation expansion;
     3. **Partner organism resolution** against the GBIF Backbone Taxonomy;
     4. **Context & negation filtering** excluding reference-list lines, negated statements, generic organism words, and plant-pathogenic/symbiotic fungi unless explicit entomopathogenic harm, infection, mycosis, or mortality to the beetle was stated.

---

### Stage 2: Cross-Check, Gap Analysis, and Figures (`src/02_crosscheck_and_analysis/`)
- **`R1`–`R6` Candidate Filtering (`06_candidate_additions.py`)**:
  - `R1`: Organism reported with relation `pathogen` or `parasite`;
  - `R2`: Organism belongs to one of the six target groups (Fungi, Bacteria, Viruses, Nematodes, Protists, Microsporidia);
  - `R3`: Fungi and bacteria restricted to entomopathogenic lineages or explicit statements of harm to the beetle;
  - `R4`: Primary research publications only (reviews, catalogs, books, websites, and second-hand citations excluded);
  - `R5`: Publication absent from the manual dataset;
  - `R6`: Focused subset (publication title contains a bark/ambrosia beetle term or genus, and supporting sentence contains no citation marker to another work).

---

## Repository Structure

```text
.
├── README.md
├── environment.yml
├── requirements.txt
├── run_all.sh
├── src/
│   ├── 01_automated_source_collection/
│   │   ├── 01_literature_api_search.py
│   │   ├── 02_structured_database_queries.py
│   │   └── 03_text_indexing_and_llm_extraction.py
│   └── 02_crosscheck_and_analysis/
│       ├── 01_profile.py
│       ├── 02_her_publications.py
│       ├── 03_search_recall.py
│       ├── 04_taxonomy.py
│       ├── 05_extraction_agreement.py
│       ├── 06_candidate_additions.py
│       ├── 06c_precision.py
│       ├── 07_trends.py
│       ├── 08_surveillance.py
│       ├── 09_check_manuscript_numbers.py
│       ├── 10_figures.py
│       ├── 11_summary_tables.py
│       ├── 12_enrich_references.py
│       └── build_collaborator_workbook.py
└── figures/
    ├── Figure5_crosscheck.pdf / .png
    ├── Figure6_rescaled_trend.pdf / .png
    └── Figure7_candidates.pdf / .png
```

---

## Reproducing the Analysis and Figures

```bash
git clone https://github.com/ChristopherMarais/BB_pathogen_data.git
cd BB_pathogen_data
conda env create -f environment.yml
conda activate barkbeetle-pathogens
bash run_all.sh
```

> **Data Availability**: The underlying manual dataset is hosted on Figshare (`bark_ambrosia_beetle_pathogen_database_v1.1.csv`). To run the full analysis pipeline locally, place data tables in the `data/` directory (git-ignored for proprietary and licensing integrity).
