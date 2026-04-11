#!/usr/bin/env bash
# 00_download.sh — download GRCh38, HG002 assemblies, RefSeq GTF, and ClinVar.
#
# Assembly URLs are resolved from the same HPRC Release-2 index CSV used by
# examples/hprc_r2/01_build_archive.sh, so both examples always stay in sync
# with the canonical HPRC data release.
#
# Files downloaded (~2.5 GB total):
#   GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz  (~900 MB)
#   hg002v1.1.mat_MT.PanSN.fa.gz                             (~800 MB)
#   hg002v1.1.pat.PanSN.fa.gz                                (~770 MB)
#   hg38.ncbiRefSeq.gtf.gz                                   (~42 MB)
#   clinvar.vcf.gz + .tbi                                    (~190 MB)
#
# Usage:
#   bash examples/hg002/00_download.sh

set -euo pipefail
cd "$(dirname "$0")"

S3_HTTP_BASE="https://human-pangenomics.s3.us-west-2.amazonaws.com"
CSV_URL="https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv"
CSV_CACHE="hprc_r2_index.csv"

GTF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_TBI_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"

dl() {
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        echo "[SKIP] $dest already exists"
    else
        echo "[DL]   $dest"
        curl -fL --retry 3 --retry-delay 5 -o "$dest" "$url"
    fi
}

# ---------------------------------------------------------------------------
# 1. Download the HPRC R2 assembly index (shared with hprc_r2/01_build_archive.sh)
# ---------------------------------------------------------------------------
dl "$CSV_URL" "$CSV_CACHE"

# ---------------------------------------------------------------------------
# 2. Resolve assembly URLs from the CSV for GRCh38, HG002-mat, HG002-pat
# ---------------------------------------------------------------------------
echo
echo "[INDEX] Resolving HG002 + GRCh38 URLs from $CSV_CACHE ..."
python3 - "$CSV_CACHE" "$S3_HTTP_BASE" <<'PYEOF'
import csv, sys, os

csv_path, http_base = sys.argv[1], sys.argv[2]

targets = {
    "GCA_000001405.15_GRCh38_no_alt_analysis_set": "GRCh38",
    "hg002v1.1.mat_MT":                            "HG002_mat",
    "hg002v1.1.pat":                               "HG002_pat",
}

found = {}
with open(csv_path, newline="") as f:
    for row in csv.DictReader(f):
        name  = row["assembly_name"].strip()
        s3url = row["assembly"].strip()
        for key, label in targets.items():
            if name.startswith(key) and label not in found:
                s3_prefix = "s3://human-pangenomics/"
                key_path  = s3url[len(s3_prefix):]
                http_url  = f"{http_base}/{key_path}"
                found[label] = (http_url, os.path.basename(key_path))

missing = [l for l in targets.values() if l not in found]
if missing:
    print(f"ERROR: could not resolve URLs for: {missing}", flush=True)
    sys.exit(1)

for label, (url, fname) in sorted(found.items()):
    # Write as shell assignments so the parent script can source them
    print(f'{label}_URL="{url}"')
    print(f'{label}_FILE="{fname}"')
PYEOF

# Source the resolved URLs into this shell
eval "$(python3 - "$CSV_CACHE" "$S3_HTTP_BASE" <<'PYEOF'
import csv, sys, os

csv_path, http_base = sys.argv[1], sys.argv[2]
targets = {
    "GCA_000001405.15_GRCh38_no_alt_analysis_set": "GRCh38",
    "hg002v1.1.mat_MT":                            "HG002_mat",
    "hg002v1.1.pat":                               "HG002_pat",
}
found = {}
with open(csv_path, newline="") as f:
    for row in csv.DictReader(f):
        name  = row["assembly_name"].strip()
        s3url = row["assembly"].strip()
        for key, label in targets.items():
            if name.startswith(key) and label not in found:
                s3_prefix = "s3://human-pangenomics/"
                key_path  = s3url[len(s3_prefix):]
                http_url  = f"{http_base}/{key_path}"
                found[label] = (http_url, os.path.basename(key_path))
for label, (url, fname) in sorted(found.items()):
    print(f'{label}_URL="{url}"')
    print(f'{label}_FILE="{fname}"')
PYEOF
)"

# ---------------------------------------------------------------------------
# 3. Download assemblies
# ---------------------------------------------------------------------------
echo
echo "[DL] Downloading assemblies (~2.5 GB) ..."
dl "$GRCh38_URL"   "$GRCh38_FILE"
dl "$HG002_mat_URL" "$HG002_mat_FILE"
dl "$HG002_pat_URL" "$HG002_pat_FILE"

# ---------------------------------------------------------------------------
# 4. Download annotation files
# ---------------------------------------------------------------------------
echo
echo "[DL] Downloading annotation files ..."
dl "$GTF_URL"          "hg38.ncbiRefSeq.gtf.gz"
dl "$CLINVAR_URL"      "clinvar.vcf.gz"
dl "$CLINVAR_TBI_URL"  "clinvar.vcf.gz.tbi"

# ---------------------------------------------------------------------------
# 5. Write a manifest for downstream scripts to source
# ---------------------------------------------------------------------------
cat > .manifest.sh <<EOF
# Auto-generated by 00_download.sh — do not edit by hand
REF_FA="${GRCh38_FILE}"
HAP0_FA="${HG002_mat_FILE}"
HAP1_FA="${HG002_pat_FILE}"
GTF="hg38.ncbiRefSeq.gtf.gz"
CLINVAR="clinvar.vcf.gz"
EOF

echo
echo "[DONE] Manifest written to .manifest.sh"
echo "       Source it in downstream scripts with: source .manifest.sh"
echo
ls -lh "$GRCh38_FILE" "$HG002_mat_FILE" "$HG002_pat_FILE" \
        hg38.ncbiRefSeq.gtf.gz clinvar.vcf.gz
