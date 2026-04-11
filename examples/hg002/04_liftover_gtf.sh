#!/usr/bin/env bash
# 04_liftover_gtf.sh — lift RefSeq GTF annotations from GRCh38 onto HG002
# haplotype contigs and generate an HTML summary report.
#
# For each haplotype this runs pgr align liftover-gtf twice:
#   pass 1 — build the liftover SQLite database + anomaly TSVs
#   pass 2 — add --ref-fa / --query-fa to also write transcript FASTA files
#
# Requires:
#   example_output/hg002_hap{0,1}.alndb  (from 01_align_alnmap.sh)
#   .manifest.sh                          (from 00_download.sh)
#
# Output is written to example_output/
#
# Usage:
#   bash examples/hg002/04_liftover_gtf.sh

set -euo pipefail
cd "$(dirname "$0")"

PGR="${PGR:-../../target/release/pgr}"
OUT="example_output"

if [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr binary not found at $PGR" >&2; exit 1
fi
for f in "$OUT/hg002_hap0.alndb" "$OUT/hg002_hap1.alndb" .manifest.sh; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found"; exit 1; }
done
source .manifest.sh

mkdir -p "$OUT"

liftover_hap() {
    local hap="$1" alndb="$2" query_fa="$3" db="$4"
    echo "=== Haplotype ${hap} liftover ==="
    "$PGR" align liftover-gtf \
        --alndb-path   "$alndb" \
        --gtf-path     "$GTF" \
        --output-db    "$db" \
        --target-chr-prefix "GRCh38#0#" \
        --min-coverage 0.5 \
        --full-coverage 0.9 \
        --ref-fa   "$REF_FA" \
        --query-fa "$query_fa"
    echo "    Database: $db"
    echo "    Transcripts: $(sqlite3 "$db" 'SELECT COUNT(*) FROM transcript_summary;')"
}

liftover_hap 0 "$OUT/hg002_hap0.alndb" "$HAP0_FA" "$OUT/hg002_hap0_liftover.db"
echo
liftover_hap 1 "$OUT/hg002_hap1.alndb" "$HAP1_FA" "$OUT/hg002_hap1_liftover.db"
echo

# ---------------------------------------------------------------------------
# HTML liftover report (scatter plots, NM gene funnel, per-chrom coverage)
# ---------------------------------------------------------------------------
echo "=== Generating $OUT/liftover_report.html ==="
python3 helpers/generate_liftover_report.py \
    --base-dir "$OUT" \
    --out      "$OUT/liftover_report.html"

echo "Liftover report: $OUT/liftover_report.html"
echo "See lift_over_examples.md for SQL queries against the liftover databases."
echo "Next: bash 05_query_mhc.sh  or  bash 06_generate_report.sh"
