#!/usr/bin/env bash
# query_c9orf72.sh — query the C9orf72 locus against the HG002 pangenome,
# then run principal bundle decomposition and generate an interactive HTML
# visualisation.
#
# C9orf72 (chromosome 9 open reading frame 72) — GRCh38 gene span:
#   chr9:27,546,542–27,573,863  (NCBI RefSeq NM_145005)
# Query region with 100 kb flanks on both sides:
#   chr9:27,446,542–27,673,863
#
# This script reuses the pangenome archive and shimmer index built by
# 05_query_mhc.sh (hg002_pan.agcrs / hg002_chr6_pan.agcrs). Run that
# script first if the index does not exist yet.
#
# Steps:
#   1. Fetch the C9orf72 region (+10 kb flanks) from the UCSC REST API
#   2. pgr query seqs   — find HG002 contigs covering the locus
#   3. pgr bundle decomp — principal bundle decomposition of hit sequences
#   4. pgr bundle svg --html — interactive HTML visualisation
#
# Output is written to example_output/
#
# Usage:
#   bash examples/hg002/query_c9orf72.sh

set -euo pipefail
cd "$(dirname "$0")"

PGR="${PGR:-../../target/release/pgr}"
AGC_RS="${AGC_RS:-../../target/release/agc-rs}"
OUT="example_output"

for bin in "$PGR" "$AGC_RS"; do
    [[ -x "$bin" ]] || { echo "ERROR: $bin not found"; exit 1; }
done

mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# Detect available RAM and pick the same archive/index as 05_query_mhc.sh
# ---------------------------------------------------------------------------
_ram_gb=999
if [[ "$(uname)" == "Darwin" ]]; then
    _ram_bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
    _ram_gb=$(( _ram_bytes / 1073741824 ))
fi

if [[ "$(uname)" == "Darwin" && $_ram_gb -lt 32 ]]; then
    ARCHIVE="$OUT/hg002_chr6_pan.agcrs"
    DB_PREFIX="$OUT/hg002_chr6_pan"
    echo "NOTE: macOS with ${_ram_gb} GB RAM — using chr6-only pangenome."
    echo "      C9orf72 is on chr9; chr6-only mode will return no hits."
    echo "      Re-run on a machine with >= 32 GB RAM for full results."
else
    ARCHIVE="$OUT/hg002_pan.agcrs"
    DB_PREFIX="$OUT/hg002_pan"
fi

# Verify the index exists
for f in "${DB_PREFIX}.mdbi" "${DB_PREFIX}.mdbv" "${DB_PREFIX}.midx"; do
    [[ -f "$f" ]] || {
        echo "ERROR: $f not found — run 05_query_mhc.sh first to build the pangenome index."
        exit 1
    }
done

# C9orf72 query parameters
# Gene span (GRCh38 1-based): chr9:27,546,542–27,573,863
# Flank: 10,000 bp each side
# UCSC REST API uses 0-based half-open coordinates
CHROM="chr9"
GENE_START=27546542
GENE_END=27573863
FLANK=100000
QUERY_START=$(( GENE_START - FLANK ))   # 27536542 (1-based) → 27536541 (0-based)
QUERY_END=$(( GENE_END + FLANK ))       # 27583863

C9_FA="$OUT/c9orf72_query.fa"
OUT_PREFIX="$OUT/c9orf72_hits"
BUNDLE_PREFIX="$OUT/c9orf72_bundle"

# ---------------------------------------------------------------------------
# 1. Fetch C9orf72 region from UCSC REST API
# ---------------------------------------------------------------------------
if [[ ! -f "$C9_FA" ]]; then
    echo "=== [1] Fetching C9orf72 region (${CHROM}:${QUERY_START}-${QUERY_END}) from UCSC ==="
    UCSC_URL="https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=${CHROM};start=$(( QUERY_START - 1 ));end=${QUERY_END}"
    curl -fsSL "$UCSC_URL" \
        | python3 -c "
import sys, json
data = json.load(sys.stdin)
seq  = data['dna'].upper()
print('>chr9:${QUERY_START}-${QUERY_END}_C9orf72_10kb_flank')
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" > "$C9_FA"
    echo "    Written: $C9_FA  ($(wc -c < "$C9_FA" | tr -d ' ') bytes)"
else
    echo "[SKIP] $C9_FA already exists"
fi

# ---------------------------------------------------------------------------
# 2. Query C9orf72 region against the pangenome index
# ---------------------------------------------------------------------------
echo
echo "=== [2] pgr query seqs — C9orf72 vs HG002 pangenome ==="
"$PGR" query seqs \
    --pgr-db-prefix "$DB_PREFIX" \
    --query-fastx-path "$C9_FA" \
    --output-prefix "$OUT_PREFIX" \
    --memory-mode moderate \
    --merge-range-tol 50000 \
    --max-count 128 \
    --max-query-count 128 \
    --max-target-count 128 \
    --min-anchor-count 10

echo
echo "=== Hit summary (${OUT_PREFIX}.000.hit) ==="
echo "# columns: idx  query  q_bgn  q_end  q_len  anchors  src  contig  bgn  end  orient  name"
cat "${OUT_PREFIX}.000.hit" || echo "(no hits file — query may have returned no results)"

# ---------------------------------------------------------------------------
# 3. Principal bundle decomposition of hit sequences
# ---------------------------------------------------------------------------
HIT_FA="${OUT_PREFIX}.000.fa"

if [[ ! -s "$HIT_FA" ]]; then
    echo
    echo "NOTE: no hit sequences found — skipping bundle decomposition."
    exit 0
fi

echo
echo "=== [3] pgr bundle decomp — principal bundle decomposition of C9orf72 hits ==="
"$PGR" bundle decomp \
    --fastx-path "$HIT_FA" \
    --output-prefix "$BUNDLE_PREFIX"

# ---------------------------------------------------------------------------
# 4. Generate interactive HTML bundle visualisation
# ---------------------------------------------------------------------------
echo
echo "=== [4] pgr bundle svg --html — interactive bundle visualisation ==="
"$PGR" bundle svg \
    --bed-file-path "${BUNDLE_PREFIX}.bed" \
    --output-prefix "$BUNDLE_PREFIX" \
    --html

echo
echo "=== Done ==="
echo "  Query hits:  ${OUT_PREFIX}.000.hit"
echo "  Hit FASTA:   ${HIT_FA}"
echo "  Bundle BED:  ${BUNDLE_PREFIX}.bed"
echo "  Bundle HTML: ${BUNDLE_PREFIX}.html"
