#!/usr/bin/env bash
# 02_variant_vcf.sh — call diploid variants and generate chromosome alignment plots.
#
# Inputs:  example_output/hg002_hap{0,1}.alndb  (from 01_align_alnmap.sh)
# Outputs in example_output/:
#   hg002.vcf
#   hg002_hap{0,1}_aln_plot.html  (whole-genome alignment dot-plots)
#
# Usage:
#   bash examples/hg002/02_variant_vcf.sh

set -euo pipefail
cd "$(dirname "$0")"

# NOTE: This script requires pgr to be installed in your PATH.
# Install from the repository root with:
#   cargo install --path ../../pgr-bin
PGR="${PGR:-pgr}"
OUT="example_output"

if ! command -v "$PGR" &>/dev/null && [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr not found in PATH" >&2
    echo "       Install with: cargo install --path <repo-root>/pgr-bin" >&2
    exit 1
fi

for f in "$OUT/hg002_hap0.alndb" "$OUT/hg002_hap1.alndb"; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found — run 01_align_alnmap.sh first" >&2; exit 1; }
done

mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# Diploid VCF
# ---------------------------------------------------------------------------
echo "=== pgr variant diploid-vcf ==="
"$PGR" variant diploid-vcf \
    --hap0-path "$OUT/hg002_hap0.alndb" \
    --hap1-path "$OUT/hg002_hap1.alndb" \
    --output-prefix "$OUT/hg002" \
    --sample-name hg002
echo "VCF: $OUT/hg002.vcf"

echo

# ---------------------------------------------------------------------------
# Whole-genome alignment dot-plots (HTML)
# ---------------------------------------------------------------------------
echo "=== pgr plot chr-aln: hap0 ==="
"$PGR" plot chr-aln \
    --alndb-path "$OUT/hg002_hap0.alndb" \
    --output-prefix "$OUT/hg002_hap0_aln_plot"
echo "Plot: $OUT/hg002_hap0_aln_plot.html"

echo "=== pgr plot chr-aln: hap1 ==="
"$PGR" plot chr-aln \
    --alndb-path "$OUT/hg002_hap1.alndb" \
    --output-prefix "$OUT/hg002_hap1_aln_plot"
echo "Plot: $OUT/hg002_hap1_aln_plot.html"

echo
echo "Next: bash 03_annotate_vcf.sh"
