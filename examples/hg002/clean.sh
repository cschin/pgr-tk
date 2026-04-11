#!/usr/bin/env bash
# clean.sh — remove all generated outputs so run_all.sh starts fresh.
# Downloaded files in downloads/ are preserved.
#
# Usage:
#   bash examples/hg002/clean.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "Removing example_output/ ..."
rm -rf example_output/

echo "Removing .manifest.sh ..."
rm -f .manifest.sh

echo "Done. Run 00_download.sh then run_all.sh to regenerate everything."
