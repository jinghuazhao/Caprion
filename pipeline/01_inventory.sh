#!/usr/bin/env bash
# 01_inventory.sh
# Inventory RAW and mzML files without modifying them.

set -eo pipefail

BASE="${HOME}/Caprion/analysis/ZWK"
RAW_DIR="${BASE}/raw"
MZML_DIR="${BASE}/mzML"
REPORT_DIR="${HOME}/Caprion/analysis/pipeline/reports"

mkdir -p "$REPORT_DIR"

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: RAW directory not found: $RAW_DIR" >&2
    exit 1
fi

if [[ ! -d "$MZML_DIR" ]]; then
    echo "ERROR: mzML directory not found: $MZML_DIR" >&2
    exit 1
fi

RAW_LIST="${REPORT_DIR}/raw_files.txt"
MZML_LIST="${REPORT_DIR}/mzml_files.txt"
MISSING_LIST="${REPORT_DIR}/missing_mzml.txt"
EXTRA_LIST="${REPORT_DIR}/extra_mzml.txt"
SUMMARY="${REPORT_DIR}/inventory_summary.txt"

find "$RAW_DIR" -maxdepth 1 -type f \
    -iname '*.raw' -printf '%f\n' |
    sort -f > "$RAW_LIST"

find "$MZML_DIR" -maxdepth 1 -type f \
    -iname '*.mzML' -printf '%f\n' |
    sort -f > "$MZML_LIST"

# Compare filenames case-insensitively, replacing only the extension.
awk '{
    name = tolower($0)
    sub(/\.raw$/, "", name)
    print name
}' "$RAW_LIST" | sort -u > "${REPORT_DIR}/raw_stems.txt"

awk '{
    name = tolower($0)
    sub(/\.mzml$/, "", name)
    print name
}' "$MZML_LIST" | sort -u > "${REPORT_DIR}/mzml_stems.txt"

comm -23 \
    "${REPORT_DIR}/raw_stems.txt" \
    "${REPORT_DIR}/mzml_stems.txt" > "$MISSING_LIST"

comm -13 \
    "${REPORT_DIR}/raw_stems.txt" \
    "${REPORT_DIR}/mzml_stems.txt" > "$EXTRA_LIST"

{
    echo "Caprion RAW/mzML inventory"
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "RAW directory: $RAW_DIR"
    echo "mzML directory: $MZML_DIR"
    echo
    echo "RAW files: $(wc -l < "$RAW_LIST")"
    echo "mzML files: $(wc -l < "$MZML_LIST")"
    echo "RAW files without matching mzML: $(wc -l < "$MISSING_LIST")"
    echo "mzML files without matching RAW: $(wc -l < "$EXTRA_LIST")"
    echo
    echo "Missing mzML stems:"
    cat "$MISSING_LIST"
    echo
    echo "Extra mzML stems:"
    cat "$EXTRA_LIST"
} > "$SUMMARY"

cat "$SUMMARY"

echo
echo "Detailed reports: $REPORT_DIR"
