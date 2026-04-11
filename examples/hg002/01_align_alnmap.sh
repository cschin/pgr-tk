#!/usr/bin/env bash
# 01_align_alnmap.sh — align HG002 maternal and paternal assemblies to GRCh38.
#
# Produces per-haplotype alignment databases (.alndb) used by all downstream
# variant-calling and liftover steps.
#
# Runtime: ~1 h per haplotype on a 32-core machine with 64 GB RAM.
#
# Requires: files downloaded by 00_download.sh
# Output is written to example_output/
#
# Usage:
#   bash examples/hg002/01_align_alnmap.sh

set -euo pipefail
cd "$(dirname "$0")"

PGR="${PGR:-../../target/release/pgr}"
OUT="example_output"

if [[ ! -x "$PGR" ]]; then
    echo "ERROR: pgr binary not found at $PGR" >&2
    echo "       Build with: cargo build --release -p pgr-bin" >&2
    exit 1
fi

if [[ ! -f ".manifest.sh" ]]; then
    echo "ERROR: .manifest.sh not found — run 00_download.sh first" >&2
    exit 1
fi
source .manifest.sh

for f in "$REF_FA" "$HAP0_FA" "$HAP1_FA"; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found"; exit 1; }
done

mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# Haplotype 0 (maternal / mat)
# ---------------------------------------------------------------------------
echo "=== [hap0] pgr align alnmap: GRCh38 vs HG002 maternal ==="
"$PGR" align alnmap \
    --reference-fasta-path "$REF_FA" \
    --assembly-contig-path "$HAP0_FA" \
    --output-prefix "$OUT/hg002_hap0" \
    --preset default
echo "[hap0] Done. Output: $OUT/hg002_hap0.alndb"

echo

# ---------------------------------------------------------------------------
# Haplotype 1 (paternal / pat)
# ---------------------------------------------------------------------------
echo "=== [hap1] pgr align alnmap: GRCh38 vs HG002 paternal ==="
"$PGR" align alnmap \
    --reference-fasta-path "$REF_FA" \
    --assembly-contig-path "$HAP1_FA" \
    --output-prefix "$OUT/hg002_hap1" \
    --preset default
echo "[hap1] Done. Output: $OUT/hg002_hap1.alndb"

echo
echo "Both haplotypes aligned. Next: bash 02_variant_vcf.sh"
