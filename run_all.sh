#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

echo "=== Regenerating Figures 5-7 from data/out/ ==="
python src/02_crosscheck_and_analysis/10_figures.py

echo "=== Building Public Verification Excel Workbook ==="
python src/02_crosscheck_and_analysis/build_collaborator_workbook.py \
    data/out/candidate_records_public.csv \
    data/out/Pathogen_Dataset_and_Automated_Crosscheck_PUBLIC.xlsx

echo "Done."
