#!/usr/bin/env bash
# query_mhc.sh — fetch GRCh38 MHC region and query a pgr-tk AGC database
#
# Usage:
#   bash scripts/query_mhc.sh <db_prefix> <output_prefix>
#
# <db_prefix>    — prefix of the pgr-tk database (.agcrs, .mdbi, .mdbv, .midx)
# <output_prefix> — prefix for output files (hit table, FASTA)
#
# Example:
#   bash scripts/query_mhc.sh /data/hprc/HPRC-y1 /tmp/mhc_query

set -euo pipefail

DB_PREFIX="${1:?usage: $0 <db_prefix> <output_prefix>}"
OUT_PREFIX="${2:?usage: $0 <db_prefix> <output_prefix>}"

PGR_QUERY="${PGR_QUERY:-$(dirname "$0")/../target/release/pgr-query}"
if [[ ! -x "$PGR_QUERY" ]]; then
    PGR_QUERY="$(dirname "$0")/../target/debug/pgr-query"
fi
if [[ ! -x "$PGR_QUERY" ]]; then
    echo "ERROR: pgr-query binary not found. Build with: cargo build -p pgr-bin --bin pgr-query" >&2
    exit 1
fi

MHC_FA="${OUT_PREFIX}.mhc_query.fa"

# ---------------------------------------------------------------------------
# 1. Fetch GRCh38 MHC region from UCSC REST API
#    chr6:28,510,120–33,480,577  (classical MHC + extended MHC boundaries)
#    The UCSC API uses 0-based half-open coordinates, matching BED convention.
# ---------------------------------------------------------------------------
echo "[1] Fetching GRCh38 MHC region (chr6:28510120-33480577) from UCSC ..."
UCSC_URL="https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr6;start=28510120;end=33480577"
curl -fsSL "$UCSC_URL" \
    | python3 -c "
import sys, json
data = json.load(sys.stdin)
seq  = data['dna'].upper()
# wrap at 80 chars
chrom = data.get('chrom', 'chr6')
start = data.get('start', 28510120)
end   = data.get('end',   33480577)
print(f'>{chrom}:{start}-{end}')
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" > "$MHC_FA"

echo "    Written: $MHC_FA  ($(wc -c < "$MHC_FA") bytes)"

# ---------------------------------------------------------------------------
# 2. Run pgr-query
# ---------------------------------------------------------------------------
echo "[2] Running pgr-query ..."
"$PGR_QUERY" \
    "$DB_PREFIX" \
    "$MHC_FA" \
    "$OUT_PREFIX" \
    --merge-range-tol 100000 \
    --max-count 128 \
    --max-query-count 128 \
    --max-target-count 128 \
    --min-anchor-count 10

echo ""
echo "Output files:"
ls -lh "${OUT_PREFIX}".* 2>/dev/null || true
echo ""
echo "Hit summary (first 20 lines):"
head -20 "${OUT_PREFIX}.000.hit" 2>/dev/null || true
