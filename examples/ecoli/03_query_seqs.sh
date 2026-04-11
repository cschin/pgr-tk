#!/usr/bin/env bash
# 03_query_seqs.sh — query a 50 kb window of MG1655 against the 3-strain
# E. coli index and show hits in all matching strains.
#
# Also demonstrates the three --memory-mode values so you can see the
# tradeoff between RAM use and query speed on a small example.
#
# Requires: example_output/ecoli_demo.agcrs + .mdbi/.mdbv/.midx
#           (created by 01_agcrs_basics.sh and 02_build_index.sh)
#
# Output is written to example_output/
#
# Usage:
#   bash examples/ecoli/03_query_seqs.sh

set -euo pipefail
cd "$(dirname "$0")"

PGR="${PGR:-../../target/release/pgr}"
TESTDATA="../../test_data/ecoli"
OUT="example_output"
DB_PREFIX="$OUT/ecoli_demo"
QUERY_FA="$OUT/query_mg1655_50k.fa"

if [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr binary not found at $PGR" >&2
    echo "       Build with: cargo build --release -p pgr-bin" >&2
    exit 1
fi

for f in "${DB_PREFIX}.mdbi" "${DB_PREFIX}.mdbv" "${DB_PREFIX}.midx"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: $f not found — run 02_build_index.sh first" >&2
        exit 1
    fi
done

mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# Extract a 50 kb query window from MG1655 (positions 1,500,000–1,550,000)
# ---------------------------------------------------------------------------
echo "=== Extracting 50 kb query region from MG1655 ==="
python3 - "$TESTDATA/ecoli_k12_mg1655.fna.gz" "$QUERY_FA" <<'PYEOF'
import sys, gzip

fasta_gz, out_fa = sys.argv[1], sys.argv[2]
BGN, END = 1_500_000, 1_550_000          # 0-based, half-open

seq_buf = []
in_target = False
with gzip.open(fasta_gz, 'rt') as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            in_target = 'NC_000913' in line
        elif in_target:
            seq_buf.append(line)

seq = ''.join(seq_buf)[BGN:END]
with open(out_fa, 'w') as out:
    out.write(f'>MG1655_NC_000913.3:{BGN}-{END}\n')
    for i in range(0, len(seq), 80):
        out.write(seq[i:i+80] + '\n')

print(f"  Query: {out_fa}  ({END-BGN} bp)")
PYEOF

echo

# ---------------------------------------------------------------------------
# Run pgr query seqs with each memory mode
# ---------------------------------------------------------------------------
run_query() {
    local mode="$1"
    local out="$OUT/query_result_${mode}"
    echo "--- memory-mode=${mode} ---"
    "$PGR" query seqs \
        --pgr-db-prefix "$DB_PREFIX" \
        --query-fastx-path "$QUERY_FA" \
        --output-prefix "$out" \
        --memory-mode "$mode" \
        --merge-range-tol 5000 \
        --min-anchor-count 3
    echo "  Hit file: ${out}.000.hit"
    echo "  Hits:"
    grep -v '^#' "${out}.000.hit" | awk '{printf "    %-30s %s:%s-%s orient=%s\n",$2,$8,$9,$10,$11}' || true
    echo
}

echo "=== Query: 50 kb MG1655 window vs. 3-strain archive ==="
echo "(Expected: hits in MG1655, W3110, and Sakai)"
echo

run_query moderate
run_query high
run_query low

echo "All three memory modes produced results above."
echo "Done."
