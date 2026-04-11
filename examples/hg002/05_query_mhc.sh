#!/usr/bin/env bash
# 05_query_mhc.sh — build a 2-sample pangenome index (GRCh38 + HG002 mat+pat)
# and query the MHC region (chr6:28.5–33.5 Mb) against it.
#
# This demonstrates:
#   agc-rs create / append   — pack two assemblies into one archive
#   pgr index mdb            — build the shimmer index
#   pgr query seqs           — find HG002 contigs covering the MHC region
#
# The query sequence is fetched live from the UCSC REST API (GRCh38).
#
# Requires:
#   .manifest.sh  (from 00_download.sh)
#   agc-rs, pgr binaries
#
# Usage:
#   bash examples/hg002/05_query_mhc.sh

set -euo pipefail
cd "$(dirname "$0")"

PGR="${PGR:-../../target/release/pgr}"
AGC_RS="${AGC_RS:-../../target/release/agc-rs}"

for bin in "$PGR" "$AGC_RS"; do
    [[ -x "$bin" ]] || { echo "ERROR: $bin not found"; exit 1; }
done
[[ -f ".manifest.sh" ]] || { echo "ERROR: run 00_download.sh first"; exit 1; }
source .manifest.sh

ARCHIVE="hg002_pangenome.agcrs"
DB_PREFIX="hg002_pangenome"
MHC_FA="mhc_query.fa"
OUT_PREFIX="mhc_hits"

# ---------------------------------------------------------------------------
# 1. Build 2-sample archive (skip if already built)
# ---------------------------------------------------------------------------
if [[ ! -f "$ARCHIVE" ]]; then
    echo "=== [1] Building pangenome archive: $ARCHIVE ==="

    # Extract GRCh38 PanSN sample name from the FASTA header
    REF_SAMPLE=$(python3 -c "
import gzip, sys
with gzip.open(sys.argv[1], 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            # PanSN format: sample#hap#contig — take the sample field
            print(line[1:].split('#')[0].strip())
            break
" "$REF_FA")
    HAP0_SAMPLE=$(python3 -c "
import gzip, sys
with gzip.open(sys.argv[1], 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            print(line[1:].split('#')[0].strip())
            break
" "$HAP0_FA")
    HAP1_SAMPLE=$(python3 -c "
import gzip, sys
with gzip.open(sys.argv[1], 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            print(line[1:].split('#')[0].strip())
            break
" "$HAP1_FA")

    echo "    Reference: $REF_SAMPLE ($REF_FA)"
    echo "    Hap0:      $HAP0_SAMPLE ($HAP0_FA)"
    echo "    Hap1:      $HAP1_SAMPLE ($HAP1_FA)"

    "$AGC_RS" create --output "$ARCHIVE" --sample "$REF_SAMPLE" "$REF_FA"
    "$AGC_RS" append "$ARCHIVE" --sample "$HAP0_SAMPLE" "$HAP0_FA"
    "$AGC_RS" append "$ARCHIVE" --sample "$HAP1_SAMPLE" "$HAP1_FA"
    "$AGC_RS" info "$ARCHIVE"
else
    echo "[SKIP] $ARCHIVE already exists"
fi

# ---------------------------------------------------------------------------
# 2. Build shimmer index (skip if already built)
# ---------------------------------------------------------------------------
if [[ ! -f "${DB_PREFIX}.mdbi" ]]; then
    echo
    echo "=== [2] Building shimmer index: $DB_PREFIX ==="
    "$PGR" index mdb \
        --agcrs-input "$ARCHIVE" \
        --prefix "$DB_PREFIX" \
        --batch-size 1
else
    echo "[SKIP] ${DB_PREFIX}.mdbi already exists"
fi

# ---------------------------------------------------------------------------
# 3. Fetch GRCh38 MHC region from UCSC REST API
# ---------------------------------------------------------------------------
if [[ ! -f "$MHC_FA" ]]; then
    echo
    echo "=== [3] Fetching GRCh38 MHC region (chr6:28510120-33480577) from UCSC ==="
    UCSC_URL="https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr6;start=28510120;end=33480577"
    curl -fsSL "$UCSC_URL" \
        | python3 -c "
import sys, json
data = json.load(sys.stdin)
seq  = data['dna'].upper()
print('>chr6:28510120-33480577_MHC')
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" > "$MHC_FA"
    echo "    Written: $MHC_FA  ($(wc -c < "$MHC_FA" | tr -d ' ') bytes)"
else
    echo "[SKIP] $MHC_FA already exists"
fi

# ---------------------------------------------------------------------------
# 4. Query MHC against the pangenome index
# ---------------------------------------------------------------------------
echo
echo "=== [4] pgr query seqs — MHC vs HG002 pangenome ==="
"$PGR" query seqs \
    --pgr-db-prefix "$DB_PREFIX" \
    --query-fastx-path "$MHC_FA" \
    --output-prefix "$OUT_PREFIX" \
    --memory-mode moderate \
    --merge-range-tol 100000 \
    --max-count 128 \
    --max-query-count 128 \
    --max-target-count 128 \
    --min-anchor-count 10

echo
echo "=== Hit summary (${OUT_PREFIX}.000.hit) ==="
echo "# columns: idx  query  q_bgn  q_end  q_len  anchors  src  contig  bgn  end  orient  name"
cat "${OUT_PREFIX}.000.hit"
