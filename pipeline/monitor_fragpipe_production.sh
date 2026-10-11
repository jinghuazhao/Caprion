#!/usr/bin/env bash
# =============================================================================
# 10_monitor_fragpipe_production.sh
#
# Purpose:
#   Monitor FragPipe production runs and inspect their logs and outputs.
#
# Inputs:
#   production_runs/production_98_*/
#   reports/fragpipe_production_*.log
#
# Outputs:
#   Terminal report only.
#
# Safety:
#   Read-only. Does not stop processes, alter files, or delete results.
#
# Execution:
#   bash 10_monitor_fragpipe_production.sh
#   bash 10_monitor_fragpipe_production.sh production_runs/production_98_<timestamp>
# =============================================================================

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
RUN_ROOT="${PIPELINE}/production_runs"

if [[ ! -d "${RUN_ROOT}" ]]; then
    echo "No production_runs directory found: ${RUN_ROOT}"
    exit 1
fi

if [[ $# -gt 0 ]]; then
    RUN_DIR="$1"
    [[ "${RUN_DIR}" = /* ]] || RUN_DIR="${PIPELINE}/${RUN_DIR}"
    if [[ ! -d "${RUN_DIR}" ]]; then
        echo "ERROR: run directory not found: ${RUN_DIR}" >&2
        exit 1
    fi
    RUNS=("${RUN_DIR}")
else
    shopt -s nullglob
    RUNS=("${RUN_ROOT}"/production_98_*)
    shopt -u nullglob
fi

if [[ ${#RUNS[@]} -eq 0 ]]; then
    echo "No production run directories found."
    exit 0
fi

for run in "${RUNS[@]}"; do
    echo
    echo "============================================================"
    echo "Run: ${run}"
    echo "============================================================"

    if [[ -f "${run}/run_metadata.txt" ]]; then
        echo "--- Run metadata ---"
        cat "${run}/run_metadata.txt"
    else
        echo "Metadata file not found."
    fi

    echo
    echo "--- FragPipe log candidates ---"
    find "${PIPELINE}/reports" -maxdepth 1 -type f \
        -name 'fragpipe_production_*.log' -print 2>/dev/null || true

    echo
    echo "--- Run directory summary ---"
    find "${run}" -maxdepth 2 -type f -printf '%TY-%Tm-%Td %TH:%TM\t%s bytes\t%p\n' \
        2>/dev/null | sort | tail -n 40 || true

    echo
    echo "--- Output file counts ---"
    for pattern in '*.pepXML' '*.pin' '*.tsv' '*.mzML' '*.log'; do
        count=$(find "${run}" -type f -name "${pattern}" 2>/dev/null | wc -l)
        printf '%-12s %s\n' "${pattern}" "${count}"
    done

    echo
    echo "--- Recent FragPipe messages ---"
    if [[ -f "${run}/run_metadata.txt" ]]; then
        timestamp=$(sed -n 's/^Started: //p' "${run}/run_metadata.txt" | head -n 1)
        if [[ -n "${timestamp}" ]]; then
            echo "Started: ${timestamp}"
        fi
    fi

    echo "Note: the production log is stored under ${PIPELINE}/reports/."
done

echo
echo "Monitoring complete. No files were changed."
