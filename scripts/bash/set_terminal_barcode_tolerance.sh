#!/bin/bash
# Set the mismatch tolerance of every terminal (Y) barcode in a barcode config file (fix 5).
# The 24 NYBot terminal barcodes are >= 5 mismatches apart over their shared length, so a
# tolerance of 1 (or 2) cannot assign a read to the wrong terminal barcode.
# Usage: scripts/bash/set_terminal_barcode_tolerance.sh <barcode_config.txt> [mismatches=1]
# Edits the file in place; keeps its line endings (the SPIDR configs use CRLF).
set -euo pipefail
cfg=$1; mm=${2:-1}
[ -f "$cfg" ] || { echo "no such file: $cfg" >&2; exit 1; }
before=$(grep -cP '^Y\t' "$cfg" || true)
sed -i -E "s/^(Y\t[^\t]+\t[ACGTN]+\t)[0-9]+(\r?)$/\1${mm}\2/" "$cfg"
after=$(grep -cP "^Y\t[^\t]+\t[ACGTN]+\t${mm}\r?$" "$cfg" || true)
echo "$cfg: $after of $before Y barcodes now allow $mm mismatch(es)"
[ "$after" -eq "$before" ] || { echo "warning: not every Y row was updated" >&2; exit 2; }
