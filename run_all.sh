#!/usr/bin/env bash
# Reproduce All_figures.pdf from the exported cell tables.
#
# Usage:  bash run_all.sh <folder with *_cells.csv> [results folder]
#
# Step 0 (export of the curated Matilda trees to *_cells.csv) needs the raw data on the cluster and is run separately:
#   python export/export_lineage_cells.py --output <tables folder>
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
A="$HERE/analysis"
TABLES="$(cd "$1" && pwd)"
OUT="$(mkdir -p "${2:-results}" && cd "${2:-results}" && pwd)"
export PYTHONPATH="$A:${PYTHONPATH:-}"

run_analysis () {   # $1 = results folder
  cd "$1"
  python "$A/01_build_dataset.py" "$TABLES"
  # independent analyses, run in parallel (permutations and bootstraps take a few minutes)
  python "$A/02_birth_vs_division.py" received &
  python "$A/02_birth_vs_division.py" lin3 &
  python "$A/02_birth_vs_division.py" lin2 &
  python "$A/03_lineage_analyses.py" basic &
  python "$A/03_lineage_analyses.py" patterns_mouse &
  python "$A/03_lineage_analyses.py" patterns_human &
  python "$A/03_lineage_analyses.py" sistcous &
  python "$A/04_published_conventions.py" &
  wait
}

echo "== main analysis (mouse generations 0-3)"
run_analysis "$OUT"
echo "== robustness run (mouse generations 0-4)"
mkdir -p "$OUT/check_gen4"
MAX_GEN_MOUSE=4 run_analysis "$OUT/check_gen4"

echo "== figures"
cd "$OUT"
FIG="$OUT/figures"; mkdir -p "$FIG"
python "$A/fig_cartoon.py" "$FIG"
python "$A/fig_supplementary.py" "$FIG" fig_depth fig_pctdiff fig_generations fig_patterns fig_sisters
python "$A/fig_main.py" "$FIG"
python "$A/fig_grids.py" "$FIG"
python "$A/fig_published.py" "$FIG"
python "$A/fig_generation4_check.py" "$FIG" "$OUT/check_gen4"
python "$A/assemble_pdf.py" "$FIG" "$OUT/All_figures.pdf"
