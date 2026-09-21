#!/usr/bin/env bash
#
# Merge per-library QC summary tables into a single TSV.
# Input tables must have metrics as rows and one library column each.
#
# Usage:
#   bash scripts/bash/merge_summary_tables.sh \
#       --input workup/qc/summary/*.tsv \
#       --output workup/qc/summary/merged_experiment.tsv
#
set -euo pipefail

INPUT_FILES=()
OUTPUT=""

usage() {
    echo "Usage: $0 --input <file> [<file> ...] --output <merged.tsv>" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                INPUT_FILES+=("$1")
                shift
            done
            ;;
        --output)
            [[ $# -ge 2 ]] || usage
            OUTPUT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            ;;
    esac
done

if [[ ${#INPUT_FILES[@]} -eq 0 || -z "${OUTPUT}" ]]; then
    usage
fi

for file in "${INPUT_FILES[@]}"; do
    if [[ ! -f "${file}" ]]; then
        echo "Error: input file not found: ${file}" >&2
        exit 1
    fi
done

mkdir -p "$(dirname "${OUTPUT}")"

python - <<'PY' "${INPUT_FILES[@]}" "${OUTPUT}"
import sys

import pandas as pd

input_files = sys.argv[1:-1]
output_path = sys.argv[-1]

merged = None
for path in input_files:
    table = pd.read_csv(path, sep="\t")
    if "metric" not in table.columns:
        raise SystemExit(f"Error: missing 'metric' column in {path}")
    if table.shape[1] != 2:
        raise SystemExit(f"Error: expected one library column in {path}, found {table.shape[1] - 1}")

    if merged is None:
        merged = table.copy()
        continue

    library_column = [col for col in table.columns if col != "metric"][0]
    if library_column in merged.columns:
        raise SystemExit(f"Error: duplicate library column '{library_column}'")

    merged = merged.merge(table, on="metric", how="outer")

merged.to_csv(output_path, sep="\t", index=False)
print(f"Merged {len(input_files)} tables into {output_path}")
PY
