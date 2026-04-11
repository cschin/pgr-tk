#!/usr/bin/env bash
# 02_variant_vcf.sh — call diploid variants and generate chromosome alignment plots.
#
# Inputs:  hg002_hap0.alndb, hg002_hap1.alndb  (from 01_align_alnmap.sh)
# Outputs: hg002.vcf
#          hg002_hap{0,1}_aln_plot.html  (whole-genome alignment dot-plots)
#
# Usage:
#   bash examples/hg002/02_variant_vcf.sh

set -euo pipefail
cd "$(dirname "$0")"

PGR="${PGR:-../../target/release/pgr}"

if [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr binary not found at $PGR" >&2
    exit 1
fi

for f in hg002_hap0.alndb hg002_hap1.alndb; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found — run 01_align_alnmap.sh first" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# Diploid VCF
# ---------------------------------------------------------------------------
echo "=== pgr variant diploid-vcf ==="
"$PGR" variant diploid-vcf \
    --hap0-path hg002_hap0.alndb \
    --hap1-path hg002_hap1.alndb \
    --output-prefix hg002 \
    --sample-name hg002
echo "VCF: hg002.vcf"

echo

# ---------------------------------------------------------------------------
# Whole-genome alignment dot-plots (HTML)
# ---------------------------------------------------------------------------
echo "=== pgr plot chr-aln: hap0 ==="
"$PGR" plot chr-aln \
    --ctgmap-json-path hg002_hap0.ctgmap.json \
    --output-prefix hg002_hap0_aln_plot
echo "Plot: hg002_hap0_aln_plot.html"

echo "=== pgr plot chr-aln: hap1 ==="
"$PGR" plot chr-aln \
    --ctgmap-json-path hg002_hap1.ctgmap.json \
    --output-prefix hg002_hap1_aln_plot
echo "Plot: hg002_hap1_aln_plot.html"

echo
echo "Next: bash 03_annotate_vcf.sh"
