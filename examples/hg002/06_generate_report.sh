#!/usr/bin/env bash
# 06_generate_report.sh — generate the final end-to-end HTML report.
#
# Reads all outputs from example_output/ (alignment plots, VCF, ClinVar
# annotations, liftover databases, SV candidates) and the step-timing log
# produced by run_all.sh to produce a single self-contained HTML report.
#
# This script is idempotent: re-running it regenerates the report from
# whatever output files are currently present.
#
# Requires:
#   helpers/generate_liftover_report.py
#   example_output/run_all_timings.tsv   (optional — from run_all.sh)
#
# Output:
#   example_output/e2e_report.html
#   example_output/e2e_report_lite.html
#
# Usage:
#   bash examples/hg002/06_generate_report.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT="example_output"
mkdir -p "$OUT"

TIMELOG="$OUT/run_all_timings.tsv"
REPORT="$OUT/e2e_report.html"

echo "=== Generating $REPORT ==="
python3 helpers/generate_liftover_report.py \
    --base-dir "$OUT" \
    --timelog  "$TIMELOG" \
    --out      "$REPORT"

echo "Report:      $REPORT"
echo "Lite report: ${REPORT%.html}_lite.html"
