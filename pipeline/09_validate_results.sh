#!/usr/bin/env bash
# =============================================================================
# 09_validate_results.sh
#
# Purpose:
#   Validate the completed FragPipe production run and inventory its outputs.
#
# Inputs:
#   Latest production_runs/production_98_*/
#   production_v2.fp-manifest
#   FragPipe production log under reports/
#
# Outputs:
#   reports/results_validation_<timestamp>.txt
#
# Safety:
#   Read-only. Does not rerun FragPipe, modify results, or delete files.
#
# Execution:
#   bash 09_validate_results.sh
#   bash 09_validate_results.sh production_runs/production_98_<timestamp>
# =============================================================================

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
RUN_ROOT="${PIPELINE}/production_runs"
MANIFEST="${PIPELINE}/production_v2.fp-manifest"
REPORT_DIR="${PIPELINE}/reports"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
REPORT="${REPORT_DIR}/results_validation_${TIMESTAMP}.txt"

mkdir -p "${REPORT_DIR}"

if [[ $# -gt 0 ]]; then
    RUN_DIR="$1"
    [[ "${RUN_DIR}" = /* ]] || RUN_DIR="${PIPELINE}/${RUN_DIR}"
else
    shopt -s nullglob
    RUNS=("${RUN_ROOT}"/production_98_*)
    shopt -u nullglob

    if [[ ${#RUNS[@]} -eq 0 ]]; then
        echo "ERROR: no production run directories found." >&2
        exit 1
    fi

    RUN_DIR="${RUNS[0]}"
    for candidate in "${RUNS[@]}"; do
        [[ "${candidate}" -nt "${RUN_DIR}" ]] && RUN_DIR="${candidate}"
    done
fi

if [[ ! -d "${RUN_DIR}" ]]; then
    echo "ERROR: run directory not found: ${RUN_DIR}" >&2
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${MANIFEST}" >&2
    exit 1
fi

{
    echo "FRAGPIPE RESULTS VALIDATION"
    echo "==========================="
    echo "Run directory: ${RUN_DIR}"
    echo "Manifest: ${MANIFEST}"
    echo "Report time: $(date --iso-8601=seconds)"
    echo

    echo "1. MANIFEST"
    python - "${MANIFEST}" <<'PY'
import os
import sys

path = sys.argv[1]
rows = []
errors = []

with open(path, encoding="utf-8") as handle:
    for lineno, line in enumerate(handle, 1):
        line = line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 4:
            errors.append(f"Line {lineno}: expected 4 columns, got {len(fields)}")
            continue
        mzml, experiment, third, acquisition = fields
        rows.append(fields)
        if not os.path.isfile(mzml):
            errors.append(f"Line {lineno}: missing mzML: {mzml}")
        if third != "":
            errors.append(f"Line {lineno}: third column is not empty")
        if acquisition != "DDA":
            errors.append(f"Line {lineno}: acquisition is not DDA")

groups = {}
for _, experiment, _, _ in rows:
    groups[experiment] = groups.get(experiment, 0) + 1

print(f"Manifest rows: {len(rows)}")
print(f"Experiment groups: {len(groups)}")
for group, count in sorted(groups.items()):
    print(f"  {group}: {count}")

if len(rows) != 98:
    errors.append(f"Expected 98 manifest rows; found {len(rows)}")

if errors:
    print("Manifest validation: FAIL")
    for error in errors:
        print(f"  ERROR: {error}")
else:
    print("Manifest validation: PASS")
PY

    echo
    echo "2. RUN METADATA"
    if [[ -f "${RUN_DIR}/run_metadata.txt" ]]; then
        cat "${RUN_DIR}/run_metadata.txt"
    else
        echo "WARNING: run_metadata.txt not found."
    fi

    echo
    echo "3. OUTPUT FILE INVENTORY"
    find "${RUN_DIR}" -type f -printf '%s bytes\t%p\n' | sort -k2
    echo
    echo "Total files: $(find "${RUN_DIR}" -type f | wc -l)"
    echo "Empty files: $(find "${RUN_DIR}" -type f -empty | wc -l)"

    echo
    echo "4. OUTPUT COUNTS"
    for pattern in '*.pepXML' '*.pin' '*.tsv' '*.txt' '*.log' '*.mzML'; do
        count=$(find "${RUN_DIR}" -type f -name "${pattern}" | wc -l)
        printf '%-12s %s\n' "${pattern}" "${count}"
    done

    echo
    echo "5. SEARCH FOR FDR-RELATED REPORTS"
    find "${RUN_DIR}" -type f \
        \( -iname '*psm*' -o -iname '*peptide*' -o \
           -iname '*protein*' -o -iname '*report*' -o \
           -iname '*combined*' \) -print | sort

    echo
    echo "6. FDR CONFIGURATION WARNING"
    grep -nE 'ionquant\.(peptidefdr|proteinfdr|ionfdr)|With (peptide|protein) FDR' \
        "${RUN_DIR}/workflow_used.workflow" \
        "${REPORT_DIR}"/fragpipe_production_*.log 2>/dev/null || true

    echo
    echo "IMPORTANT INTERPRETATION"
    echo "This script inventories outputs; it does not establish FDR control."
    echo "The run log reported peptide and protein FDR values of 1.0."
    echo "Inspect the actual Philosopher/Percolator reports before claiming"
    echo "1% peptide or protein FDR. Do not infer FDR from file counts alone."
    echo
    echo "VALIDATION REPORT GENERATED"
    echo "Detailed scientific validation remains pending."
} | tee "${REPORT}"

echo
echo "Report saved to: ${REPORT}"
