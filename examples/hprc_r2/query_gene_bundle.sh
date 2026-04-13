#!/usr/bin/env bash
# query_gene_bundle.sh — query any gene locus from a pangenome index,
# run principal bundle decomposition on all haplotype hits, and generate
# an interactive HTML visualisation.
#
# Gene coordinates are derived from a GTF annotation file (GRCh38 RefSeq
# or Ensembl).  The reference sequence for the query is fetched directly
# from the agcrs archive via pgr query fetch — no UCSC API call needed.
#
# Steps:
#   1. Parse GTF → gene chromosome, start, end
#   2. Add flanks on both sides
#   3. Find the matching GRCh38 contig in the pangenome DB (pgr query fetch --list)
#   4. Fetch the reference sequence for the locus (pgr query fetch)
#   5. Query all haplotype contigs covering the locus (pgr query seqs)
#   6. Principal bundle decomposition of hit sequences (pgr bundle decomp)
#   7. Interactive HTML bundle visualisation (pgr bundle svg --html)
#
# Requirements: pgr in PATH
#   cargo install --path <repo-root>/pgr-bin
#
# Usage:
#   bash query_gene_bundle.sh <gene_name> <gtf_file> [db_prefix] [flank_bp] [out_dir]
#
#   gene_name  — gene name or gene_id to search in the GTF (e.g. C9orf72, BRCA1)
#   gtf_file   — path to a GRCh38 GTF file (plain or .gz)
#   db_prefix  — prefix shared by .agcrs, .mdbi, .mdbv, .midx  (default: hprc_r2)
#   flank_bp   — flanking bases added on each side               (default: 100000)
#   out_dir    — output directory                                 (default: <gene_name>_bundle)
#
# Environment overrides:
#   PGR=<path>           — path to pgr binary
#   REF_SAMPLE_HINT=<s>  — substring used to identify the reference sample
#                          in the DB (default: GRCh38)

set -euo pipefail

GENE_NAME="${1:?usage: $0 <gene_name> <gtf_file> [db_prefix] [flank_bp] [out_dir]}"
GTF_FILE="${2:?usage: $0 <gene_name> <gtf_file> [db_prefix] [flank_bp] [out_dir]}"
DB_PREFIX="${3:-hprc_r2}"
FLANK="${4:-100000}"
OUT_DIR="${5:-${GENE_NAME}_bundle}"
PGR="${PGR:-pgr}"
REF_SAMPLE_HINT="${REF_SAMPLE_HINT:-GRCh38}"

if ! command -v "$PGR" &>/dev/null && [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr not found in PATH" >&2
    echo "       Install with: cargo install --path <repo-root>/pgr-bin" >&2
    exit 1
fi
[[ -f "$GTF_FILE" ]] || { echo "ERROR: GTF file not found: $GTF_FILE" >&2; exit 1; }

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# 1. Parse GTF → gene chromosome, start (0-based), end
# ---------------------------------------------------------------------------
echo "=== [1] Parsing GTF for gene: ${GENE_NAME} ==="
GENE_COORDS="$OUT_DIR/gene_coords.tsv"

python3 - "$GTF_FILE" "$GENE_NAME" "$GENE_COORDS" <<'PYEOF'
import sys, gzip, re, os

gtf_path, gene_name, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

opener = gzip.open if gtf_path.endswith(".gz") else open
pattern = re.compile(
    r'gene_name\s+"([^"]+)"|gene_id\s+"([^"]+)"'
)

found = []
with opener(gtf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) < 9:
            continue
        feature = fields[2]
        if feature != "gene":
            continue
        attrs = fields[8]
        for m in pattern.finditer(attrs):
            name = m.group(1) or m.group(2)
            if name == gene_name:
                chrom  = fields[0]
                # GTF is 1-based inclusive; convert to 0-based half-open
                start  = int(fields[3]) - 1
                end    = int(fields[4])
                strand = fields[6]
                found.append((chrom, start, end, strand))
                break

if not found:
    print(f"ERROR: gene '{gene_name}' not found in GTF", file=sys.stderr)
    sys.exit(1)

# Use the broadest span if multiple records
chrom = found[0][0]
start = min(r[1] for r in found)
end   = max(r[2] for r in found)
strand = found[0][3]

with open(out_path, "w") as out:
    out.write(f"{chrom}\t{start}\t{end}\t{strand}\n")

print(f"  {chrom}:{start}-{end}  strand={strand}  ({end-start:,} bp)")
PYEOF

read CHROM GENE_START GENE_END STRAND < "$GENE_COORDS"
QUERY_START=$(( GENE_START - FLANK ))
QUERY_END=$(( GENE_END   + FLANK ))
(( QUERY_START < 0 )) && QUERY_START=0
echo "    Query region with ${FLANK} bp flanks: ${CHROM}:${QUERY_START}-${QUERY_END}"

# ---------------------------------------------------------------------------
# 2. List all sequences in the DB and find the matching reference contig
# ---------------------------------------------------------------------------
SEQ_LIST="$OUT_DIR/seq_list.tsv"
if [[ ! -s "$SEQ_LIST" ]]; then
    echo
    echo "=== [2] Listing sequences in ${DB_PREFIX}.midx ==="
    "$PGR" query fetch \
        --pgr-db-prefix "$DB_PREFIX" \
        --list \
        > "$SEQ_LIST"
    echo "    $(wc -l < "$SEQ_LIST") sequences indexed"
else
    echo
    echo "[SKIP] seq_list.tsv already exists"
fi

# Find reference sample and contig for the target chromosome
# seq_list columns: sid  src  ctg  length
# Match: src contains REF_SAMPLE_HINT AND ctg ends with the chromosome name
REF_REGION=$(python3 - "$SEQ_LIST" "$REF_SAMPLE_HINT" "$CHROM" "$QUERY_START" "$QUERY_END" <<'PYEOF'
import sys

seq_list, hint, chrom, q_start, q_end = \
    sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])

with open(seq_list) as fh:
    for line in fh:
        sid, src, ctg, length = line.rstrip().split("\t")
        length = int(length)
        # src must contain the hint; contig must end with the chromosome
        if hint.lower() in src.lower() and ctg.split("#")[-1] == chrom:
            # Clamp coordinates to contig bounds
            bgn = max(0, q_start)
            end = min(length, q_end)
            # label  src  ctg  bgn  end  strand(0=fwd)
            print(f"{chrom}:{bgn}-{end}_{src}\t{src}\t{ctg}\t{bgn}\t{end}\t0")
            break
    else:
        print(f"ERROR: no reference contig matching '{hint}' / '{chrom}' found",
              file=sys.stderr)
        sys.exit(1)
PYEOF
)

echo "    Reference region: $REF_REGION"

# ---------------------------------------------------------------------------
# 3. Fetch the reference sequence for the locus from the agcrs
# ---------------------------------------------------------------------------
QUERY_FA="$OUT_DIR/${GENE_NAME}_query.fa"
if [[ ! -s "$QUERY_FA" ]]; then
    echo
    echo "=== [3] Fetching reference sequence from agcrs ==="
    REF_REGION_FILE="$OUT_DIR/ref_region.tsv"
    echo "$REF_REGION" > "$REF_REGION_FILE"
    "$PGR" query fetch \
        --pgr-db-prefix "$DB_PREFIX" \
        --region-file   "$REF_REGION_FILE" \
        --output-file   "$QUERY_FA"
    echo "    Written: $QUERY_FA  ($(wc -c < "$QUERY_FA" | tr -d ' ') bytes)"
else
    echo
    echo "[SKIP] $QUERY_FA already exists"
fi

# ---------------------------------------------------------------------------
# 4. Query all haplotype contigs covering the locus
# ---------------------------------------------------------------------------
HIT_PREFIX="$OUT_DIR/${GENE_NAME}_hits"
echo
echo "=== [4] pgr query seqs — ${GENE_NAME} vs pangenome ==="
"$PGR" query seqs \
    --pgr-db-prefix    "$DB_PREFIX" \
    --query-fastx-path "$QUERY_FA" \
    --output-prefix    "$HIT_PREFIX" \
    --memory-mode      moderate \
    --merge-range-tol  100000 \
    --max-count        128 \
    --max-query-count  128 \
    --max-target-count 128 \
    --min-anchor-count 10

echo
echo "=== Hit summary (${HIT_PREFIX}.000.hit) ==="
echo "# columns: idx  query  q_bgn  q_end  q_len  anchors  src  contig  bgn  end  orient  name"
cat "${HIT_PREFIX}.000.hit"

HIT_FA="${HIT_PREFIX}.000.fa"
if [[ ! -s "$HIT_FA" ]]; then
    echo
    echo "NOTE: no hit sequences found — skipping bundle decomposition."
    exit 0
fi

# ---------------------------------------------------------------------------
# 5. Principal bundle decomposition
# ---------------------------------------------------------------------------
BUNDLE_PREFIX="$OUT_DIR/${GENE_NAME}_bundle"
echo
echo "=== [5] pgr bundle decomp ==="
"$PGR" bundle decomp \
    --fastx-path    "$HIT_FA" \
    --output-prefix "$BUNDLE_PREFIX"

# ---------------------------------------------------------------------------
# 6. Interactive HTML bundle visualisation
# ---------------------------------------------------------------------------
echo
echo "=== [6] pgr bundle svg --html ==="
"$PGR" bundle svg \
    --bed-file-path "${BUNDLE_PREFIX}.bed" \
    --output-prefix "$BUNDLE_PREFIX" \
    --html

echo
echo "=== Done ==="
echo "  Query FASTA : $QUERY_FA"
echo "  Hit FASTA   : $HIT_FA"
echo "  Bundle BED  : ${BUNDLE_PREFIX}.bed"
echo "  Bundle HTML : ${BUNDLE_PREFIX}.html"
