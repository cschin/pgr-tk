#!/usr/bin/env bash
# build_liftover_agc.sh — build per-haplotype AGC archives for pgr align liftover-gtf --agc
#
# Each archive holds two samples:
#   1. reference  (GRCh38) — first sample, used as the AGC reference stream
#   2. haplotype  (HG002)  — stored as delta vs. the reference
#
# Sample names are taken directly from the PanSN FASTA headers
# (the part before the first '#').  The alndb infers the same names at
# runtime, so no extra flags are needed when running pgr align liftover-gtf.
#
# Outputs:
#   hg002_hap0_liftover.agcrs   GRCh38 + HG002 (mat / hap0)
#   hg002_hap1_liftover.agcrs   GRCh38 + HG002 (pat / hap1)
#
# Usage:
#   cd test_genomes
#   bash build_liftover_agc.sh

set -euo pipefail

AGC=../target/release/agc-rs
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz
HAP0=hg002v1.1.mat_MT.PanSN.fa.gz
HAP1=hg002v1.1.pat.PanSN.fa.gz

# ── helper: extract sample name from first header of a (possibly gzipped) FASTA
sample_name() {
    python3 -c "
import gzip, sys
opener = gzip.open if sys.argv[1].endswith('.gz') else open
with opener(sys.argv[1], 'rt') as f:
    for line in f:
        if line.startswith('>'):
            print(line[1:].split('#')[0].split()[0])
            break
" "$1"
}

REF_SAMPLE=$(sample_name "$REF")
HAP0_SAMPLE=$(sample_name "$HAP0")
HAP1_SAMPLE=$(sample_name "$HAP1")

echo "Reference sample  : $REF_SAMPLE  ($REF)"
echo "Haplotype-0 sample: $HAP0_SAMPLE ($HAP0)"
echo "Haplotype-1 sample: $HAP1_SAMPLE ($HAP1)"
echo

# ── haplotype 0 archive ──────────────────────────────────────────────────────
OUT0=hg002_hap0_liftover.agcrs
echo "=== Building $OUT0 ==="
rm -f "$OUT0" "${OUT0}-shm" "${OUT0}-wal"

echo "  adding reference ($REF_SAMPLE) ..."
$AGC create --output "$OUT0" --sample "$REF_SAMPLE" "$REF"

echo "  appending haplotype-0 ($HAP0_SAMPLE) ..."
$AGC append "$OUT0" --sample "$HAP0_SAMPLE" "$HAP0"

echo "  done."
$AGC info "$OUT0"
echo

# ── haplotype 1 archive ──────────────────────────────────────────────────────
OUT1=hg002_hap1_liftover.agcrs
echo "=== Building $OUT1 ==="
rm -f "$OUT1" "${OUT1}-shm" "${OUT1}-wal"

echo "  adding reference ($REF_SAMPLE) ..."
$AGC create --output "$OUT1" --sample "$REF_SAMPLE" "$REF"

echo "  appending haplotype-1 ($HAP1_SAMPLE) ..."
$AGC append "$OUT1" --sample "$HAP1_SAMPLE" "$HAP1"

echo "  done."
$AGC info "$OUT1"
echo

# ── usage reminder ───────────────────────────────────────────────────────────
echo "Archives ready.  Run liftover with sequence fetching:"
echo
echo "  ../target/release/pgr align liftover-gtf \\"
echo "      --alndb-path hg002_hap0.alndb --gtf-path hg38.ncbiRefSeq.gtf.gz \\"
echo "      --output-db hg002_hap0_liftover.db \\"
echo "      --target-chr-prefix \"GRCh38#0#\" --min-coverage 0.5 \\"
echo "      --agc $OUT0"
echo
echo "  ../target/release/pgr align liftover-gtf \\"
echo "      --alndb-path hg002_hap1.alndb --gtf-path hg38.ncbiRefSeq.gtf.gz \\"
echo "      --output-db hg002_hap1_liftover.db \\"
echo "      --target-chr-prefix \"GRCh38#0#\" --min-coverage 0.5 \\"
echo "      --agc $OUT1"
