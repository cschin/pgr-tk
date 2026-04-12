#!/usr/bin/env bash
# 05_query_mhc.sh — build a pangenome index (GRCh38 + HG002 mat+pat)
# and query the MHC region (chr6:28.5–33.5 Mb) against it.
#
# On macOS with < 32 GB RAM, only chr6 sequences are used to keep build
# time and memory reasonable (~250 Mbp vs ~12 Gbp for full GRCh38).
# On Linux or macOS with >= 32 GB RAM, the full genomes are used.
#
# This demonstrates:
#   agc-rs create / append   — pack sequences into one archive
#   pgr index mdb            — build the shimmer index
#   pgr query seqs           — find HG002 contigs covering the MHC region
#   pgr bundle decomp        — principal bundle decomposition of hit sequences
#   pgr bundle svg --html    — interactive HTML visualisation of the bundles
#
# The query sequence is fetched live from the UCSC REST API (GRCh38).
#
# Requires:
#   .manifest.sh                          (from 00_download.sh)
#   example_output/hg002_hap{0,1}.ctgmap.bed  (from 01_align_alnmap.sh,
#                                              needed only in chr6 mode)
#   agc-rs, pgr binaries
#
# Output is written to example_output/
#
# Usage:
#   bash examples/hg002/05_query_mhc.sh

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

# ---------------------------------------------------------------------------
# Detect available RAM and choose full-genome vs chr6-only mode
# ---------------------------------------------------------------------------
_ram_gb=999
if [[ "$(uname)" == "Darwin" ]]; then
    _ram_bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 0)
    _ram_gb=$(( _ram_bytes / 1073741824 ))
fi

if [[ "$(uname)" == "Darwin" && $_ram_gb -lt 32 ]]; then
    CHR6_ONLY=true
    echo "NOTE: macOS with ${_ram_gb} GB RAM detected."
    echo "      Using chr6-only mode to limit memory and build time."
    echo "      (Full-genome mode requires >= 32 GB RAM)"
    ARCHIVE="$OUT/hg002_chr6_pan.agcrs"
    DB_PREFIX="$OUT/hg002_chr6_pan"
else
    CHR6_ONLY=false
    ARCHIVE="$OUT/hg002_pan.agcrs"
    DB_PREFIX="$OUT/hg002_pan"
fi

MHC_FA="$OUT/mhc_query.fa"
OUT_PREFIX="$OUT/mhc_hits"

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
# 1. Prepare input FASTAs (chr6 subset or full genome)
# ---------------------------------------------------------------------------
if $CHR6_ONLY; then
    REF_INPUT="$OUT/ref_chr6.fa"
    HAP0_INPUT="$OUT/hap0_chr6.fa"
    HAP1_INPUT="$OUT/hap1_chr6.fa"

    # GRCh38 header format: ">GRCh38#0#chr6  AC:... description"
    # match contig name before any whitespace/description
    if [[ ! -s "$REF_INPUT" ]]; then
        echo "=== [1a] Extracting chr6 from GRCh38 ==="
        extract_seqs "$REF_FA" "^>\\S*#chr6(\\s|$)" "$REF_INPUT"
    else
        echo "[SKIP] $REF_INPUT already exists"
    fi

    # HG002 header format: ">HG002#2#chr6" (clean, no description)
    if [[ ! -s "$HAP0_INPUT" ]]; then
        echo "=== [1b] Extracting chr6 contig from HG002 hap0 ==="
        extract_seqs "$HAP0_FA" "^>#*\\S*#chr6(\\s|$)" "$HAP0_INPUT"
    else
        echo "[SKIP] $HAP0_INPUT already exists"
    fi

    if [[ ! -s "$HAP1_INPUT" ]]; then
        echo "=== [1c] Extracting chr6 contig from HG002 hap1 ==="
        extract_seqs "$HAP1_FA" "^>#*\\S*#chr6(\\s|$)" "$HAP1_INPUT"
    else
        echo "[SKIP] $HAP1_INPUT already exists"
    fi
else
    REF_INPUT="$REF_FA"
    HAP0_INPUT="$HAP0_FA"
    HAP1_INPUT="$HAP1_FA"
    echo "=== [1] Using full genomes (${_ram_gb} GB RAM available) ==="
fi

# ---------------------------------------------------------------------------
# 2. Build pangenome archive (skip if valid)
# ---------------------------------------------------------------------------
_archive_ok=false
if [[ -f "$ARCHIVE" ]]; then
    _sample_count=$("$AGC_RS" list "$ARCHIVE" 2>/dev/null | wc -l || echo 0)
    [[ "$_sample_count" -ge 2 ]] && _archive_ok=true
fi

if ! $_archive_ok; then
    echo
    echo "=== [2] Building pangenome archive: $ARCHIVE ==="

    _get_sample() {
        local fa="$1"
        python3 -c "
import sys, gzip
opener = __import__('gzip').open if sys.argv[1].endswith('.gz') else open
with opener(sys.argv[1], 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            print(line[1:].split('#')[0].strip()); break
" "$fa"
    }
    REF_SAMPLE=$(_get_sample "$REF_INPUT")
    HAP0_SAMPLE=$(_get_sample "$HAP0_INPUT")
    HAP1_SAMPLE=$(_get_sample "$HAP1_INPUT")

    echo "    Reference: $REF_SAMPLE"
    echo "    Hap0:      $HAP0_SAMPLE"
    echo "    Hap1:      $HAP1_SAMPLE"

    # Force distinct sample names: hap0 = maternal, hap1 = paternal
    HAP0_SAMPLE="${HAP0_SAMPLE}_mat"
    HAP1_SAMPLE="${HAP1_SAMPLE}_pat"
    echo "    (using sample names: $REF_SAMPLE / $HAP0_SAMPLE / $HAP1_SAMPLE)"

    rm -f "$ARCHIVE" "$ARCHIVE-wal" "$ARCHIVE-shm"
    "$AGC_RS" create --output "$ARCHIVE" --sample "$REF_SAMPLE"  "$REF_INPUT"
    "$AGC_RS" append "$ARCHIVE"          --sample "$HAP0_SAMPLE" "$HAP0_INPUT"
    "$AGC_RS" append "$ARCHIVE"          --sample "$HAP1_SAMPLE" "$HAP1_INPUT"
    "$AGC_RS" info "$ARCHIVE"
else
    echo "[SKIP] $ARCHIVE already exists ($("$AGC_RS" list "$ARCHIVE" | wc -l | tr -d ' ') samples)"
fi

# ---------------------------------------------------------------------------
# 3. Build shimmer index (skip if already built and non-empty)
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
# 4. Fetch GRCh38 MHC region from UCSC REST API
# ---------------------------------------------------------------------------
if [[ ! -f "$MHC_FA" ]]; then
    echo
    echo "=== [4] Fetching GRCh38 MHC region (chr6:28510120-33480577) from UCSC ==="
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
# 5. Query MHC against the pangenome index
# ---------------------------------------------------------------------------
echo
echo "=== [5] pgr query seqs — MHC vs HG002 pangenome ==="
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

# ---------------------------------------------------------------------------
# 6. Principal bundle decomposition of hit sequences
# ---------------------------------------------------------------------------
BUNDLE_PREFIX="$OUT/mhc_bundle"
HIT_FA="${OUT_PREFIX}.000.fa"

echo
echo "=== [6] pgr bundle decomp — principal bundle decomposition of MHC hits ==="
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
echo "Bundle HTML: ${BUNDLE_PREFIX}.html"
