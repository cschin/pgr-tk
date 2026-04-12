#!/usr/bin/env bash
# query_c9orf72.sh — query the C9orf72 locus against a chr9-specific HG002
# pangenome, then run principal bundle decomposition and generate an
# interactive HTML visualisation.
#
# C9orf72 (chromosome 9 open reading frame 72) — GRCh38 gene span:
#   chr9:27,546,542–27,573,863  (NCBI RefSeq NM_145005)
# Query region with 100 kb flanks on both sides:
#   chr9:27,446,542–27,673,863
#
# This script builds a dedicated chr9 pangenome archive and shimmer index
# (hg002_chr9_pan.agcrs / hg002_chr9_pan.*) from the same input FASTAs
# used by the main pipeline. Run 00_download.sh first.
#
# Steps:
#   1. Extract chr9 sequences from GRCh38 reference and both HG002 haplotypes
#   2. Build chr9 pangenome archive (agc-rs create / append)
#   3. Build shimmer index (pgr index mdb)
#   4. Fetch the C9orf72 region (+100 kb flanks) from the UCSC REST API
#   5. pgr query seqs   — find HG002 contigs covering the locus
#   6. pgr bundle decomp — principal bundle decomposition of hit sequences
#   7. pgr bundle svg --html — interactive HTML visualisation
#
# Output is written to example_output/
#
# Usage:
#   bash examples/hg002/query_c9orf72.sh

set -euo pipefail
cd "$(dirname "$0")"

# NOTE: This script requires pgr and agc-rs to be installed in your PATH.
# Install from the repository root with:
#   cargo install --path ../../pgr-bin
#   cargo install --path ../../agc-rs
PGR="${PGR:-pgr}"
AGC_RS="${AGC_RS:-agc-rs}"
OUT="example_output"

for bin in "$PGR" "$AGC_RS"; do
    if ! command -v "$bin" &>/dev/null && [[ ! -x "$bin" ]]; then
        echo "ERROR: $bin not found in PATH" >&2
        echo "       Install with: cargo install --path <repo-root>/${bin%-*}" >&2
        exit 1
    fi
done
[[ -f ".manifest.sh" ]] || { echo "ERROR: run 00_download.sh first"; exit 1; }
source .manifest.sh

mkdir -p "$OUT"

ARCHIVE="$OUT/hg002_chr9_pan.agcrs"
DB_PREFIX="$OUT/hg002_chr9_pan"

CHROM="chr9"
GENE_START=27546542
GENE_END=27573863
FLANK=10000
QUERY_START=$(( GENE_START - FLANK ))
QUERY_END=$(( GENE_END + FLANK ))

C9_FA="$OUT/c9orf72_query.fa"
OUT_PREFIX="$OUT/c9orf72_hits"
BUNDLE_PREFIX="$OUT/c9orf72_bundle"

# ---------------------------------------------------------------------------
# Helper: extract sequences matching a regex from a (possibly gzipped) FASTA
# ---------------------------------------------------------------------------
extract_seqs() {
    local src="$1" pattern="$2" dest="$3"
    python3 - "$src" "$pattern" "$dest" <<'PYEOF'
import sys, gzip, re
src, pattern, dest = sys.argv[1], sys.argv[2], sys.argv[3]
rx = re.compile(pattern)
keep = False
written = 0
opener = gzip.open if src.endswith('.gz') else open
with opener(src, 'rt') as fh, open(dest, 'w') as out:
    for line in fh:
        if line.startswith('>'):
            keep = bool(rx.search(line))
            if keep:
                out.write(line)
                written += 1
        elif keep:
            out.write(line)
print(f"  extracted {written} sequences → {dest}", flush=True)
PYEOF
}

# ---------------------------------------------------------------------------
# 1. Extract chr9 sequences
# ---------------------------------------------------------------------------
REF_CHR9="$OUT/ref_chr9.fa"
HAP0_CHR9="$OUT/hap0_chr9.fa"
HAP1_CHR9="$OUT/hap1_chr9.fa"

if [[ ! -s "$REF_CHR9" ]]; then
    echo "=== [1a] Extracting chr9 from GRCh38 ==="
    extract_seqs "$REF_FA" "^>\\S*#chr9(\\s|$)" "$REF_CHR9"
else
    echo "[SKIP] $REF_CHR9 already exists"
fi

if [[ ! -s "$HAP0_CHR9" ]]; then
    echo "=== [1b] Extracting chr9 from HG002 hap0 ==="
    extract_seqs "$HAP0_FA" "^>#*\\S*#chr9(\\s|$)" "$HAP0_CHR9"
else
    echo "[SKIP] $HAP0_CHR9 already exists"
fi

if [[ ! -s "$HAP1_CHR9" ]]; then
    echo "=== [1c] Extracting chr9 from HG002 hap1 ==="
    extract_seqs "$HAP1_FA" "^>#*\\S*#chr9(\\s|$)" "$HAP1_CHR9"
else
    echo "[SKIP] $HAP1_CHR9 already exists"
fi

# ---------------------------------------------------------------------------
# 2. Build chr9 pangenome archive
# ---------------------------------------------------------------------------
_archive_ok=false
if [[ -f "$ARCHIVE" ]]; then
    _sample_count=$("$AGC_RS" list "$ARCHIVE" 2>/dev/null | wc -l || echo 0)
    [[ "$_sample_count" -ge 2 ]] && _archive_ok=true
fi

if ! $_archive_ok; then
    echo
    echo "=== [2] Building chr9 pangenome archive: $ARCHIVE ==="

    _get_sample() {
        python3 -c "
import sys, gzip
opener = __import__('gzip').open if sys.argv[1].endswith('.gz') else open
with opener(sys.argv[1], 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            print(line[1:].split('#')[0].strip()); break
" "$1"
    }

    REF_SAMPLE=$(_get_sample "$REF_CHR9")
    HAP0_SAMPLE="$(_get_sample "$HAP0_CHR9")_mat"
    HAP1_SAMPLE="$(_get_sample "$HAP1_CHR9")_pat"
    echo "    Reference: $REF_SAMPLE  Hap0: $HAP0_SAMPLE  Hap1: $HAP1_SAMPLE"

    rm -f "$ARCHIVE" "$ARCHIVE-wal" "$ARCHIVE-shm"
    "$AGC_RS" create --output "$ARCHIVE" --sample "$REF_SAMPLE"  "$REF_CHR9"
    "$AGC_RS" append "$ARCHIVE"          --sample "$HAP0_SAMPLE" "$HAP0_CHR9"
    "$AGC_RS" append "$ARCHIVE"          --sample "$HAP1_SAMPLE" "$HAP1_CHR9"
    "$AGC_RS" info "$ARCHIVE"
else
    echo "[SKIP] $ARCHIVE already exists ($("$AGC_RS" list "$ARCHIVE" | wc -l | tr -d ' ') samples)"
fi

# ---------------------------------------------------------------------------
# 3. Build shimmer index
# ---------------------------------------------------------------------------
_index_ok=false
if [[ -f "${DB_PREFIX}.mdbi" && -s "${DB_PREFIX}.mdbi" && \
      -f "${DB_PREFIX}.mdbv" && -s "${DB_PREFIX}.mdbv" ]]; then
    _index_ok=true
fi

if ! $_index_ok; then
    echo
    echo "=== [3] Building shimmer index: $DB_PREFIX ==="
    "$PGR" index mdb \
        --agcrs-input "$ARCHIVE" \
        --prefix "$DB_PREFIX" \
        --batch-size 0
else
    echo "[SKIP] ${DB_PREFIX}.mdbi already exists"
fi

# ---------------------------------------------------------------------------
# 4. Fetch C9orf72 region from UCSC REST API
# ---------------------------------------------------------------------------
if [[ ! -f "$C9_FA" ]]; then
    echo
    echo "=== [4] Fetching C9orf72 region (${CHROM}:${QUERY_START}-${QUERY_END}) from UCSC ==="
    UCSC_URL="https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=${CHROM};start=$(( QUERY_START - 1 ));end=${QUERY_END}"
    curl -fsSL "$UCSC_URL" \
        | python3 -c "
import sys, json
data = json.load(sys.stdin)
seq  = data['dna'].upper()
print('>chr9:${QUERY_START}-${QUERY_END}_C9orf72_100kb_flank')
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" > "$C9_FA"
    echo "    Written: $C9_FA  ($(wc -c < "$C9_FA" | tr -d ' ') bytes)"
else
    echo "[SKIP] $C9_FA already exists"
fi

# ---------------------------------------------------------------------------
# 5. Query C9orf72 region against the chr9 pangenome index
# ---------------------------------------------------------------------------
echo
echo "=== [5] pgr query seqs — C9orf72 vs HG002 chr9 pangenome ==="
"$PGR" query seqs \
    --pgr-db-prefix "$DB_PREFIX" \
    --query-fastx-path "$C9_FA" \
    --output-prefix "$OUT_PREFIX" \
    --memory-mode moderate \
    --merge-range-tol 500 \
    --max-count 128 \
    --max-query-count 128 \
    --max-target-count 128 \
    --min-anchor-count 1

echo
echo "=== Hit summary (${OUT_PREFIX}.000.hit) ==="
echo "# columns: idx  query  q_bgn  q_end  q_len  anchors  src  contig  bgn  end  orient  name"
cat "${OUT_PREFIX}.000.hit" || echo "(no hits file — query may have returned no results)"

# ---------------------------------------------------------------------------
# 6. Principal bundle decomposition of hit sequences
# ---------------------------------------------------------------------------
HIT_FA="${OUT_PREFIX}.000.fa"

if [[ ! -s "$HIT_FA" ]]; then
    echo
    echo "NOTE: no hit sequences found — skipping bundle decomposition."
    exit 0
fi

echo
echo "=== [6] pgr bundle decomp — principal bundle decomposition of C9orf72 hits ==="
"$PGR" bundle decomp \
    --fastx-path "$HIT_FA" \
    --output-prefix "$BUNDLE_PREFIX"

# ---------------------------------------------------------------------------
# 7. Generate interactive HTML bundle visualisation
# ---------------------------------------------------------------------------
echo
echo "=== [7] pgr bundle svg --html — interactive bundle visualisation ==="
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
