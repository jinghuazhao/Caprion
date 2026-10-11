#!/usr/bin/env bash
# =============================================================================
# 09_run_fragpipe_production.sh
#
# Purpose:
#   Run FragPipe 22.0 on the validated 98-file Caprion DDA dataset.
#
# Inputs:
#   production_v2.fp-manifest
#   LFQ-MBR-Sage-comparison.workflow
#   Validated canonical FASTA database
#
# Outputs:
#   production_runs/production_98_<timestamp>/
#   reports/production_preflight_<timestamp>.txt
#   reports/fragpipe_production_<timestamp>.log
#
# Safety:
#   - Refuses to reuse an existing output directory.
#   - Does not modify the source workflow or manifest.
#   - Does not delete or modify the six-file pilot.
#   - Stops on failed preflight checks.
#
# Execution:
#   Run only in a suitable allocated compute job.
#   bash 09_run_fragpipe_production.sh
# =============================================================================

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
MANIFEST="${PIPELINE}/production_v2.fp-manifest"
WORKFLOW="${PIPELINE}/LFQ-MBR-Sage-comparison.workflow"
DATABASE="${HOME}/Caprion/analysis/ZWK/discovery_fresh/database_canonical/2026-10-09-decoys-UP000005640-canonical-proteins-20420.fasta.fas"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RUN_ROOT="${PIPELINE}/production_runs"
RUN_DIR="${RUN_ROOT}/production_98_${TIMESTAMP}"
REPORT="${PIPELINE}/reports/production_preflight_${TIMESTAMP}.txt"
LOG="${PIPELINE}/reports/fragpipe_production_${TIMESTAMP}.log"

# Adjust these to the resources actually allocated by the scheduler.
THREADS="${SLURM_CPUS_PER_TASK:-8}"
RAM_GB="${FRAGPIPE_RAM_GB:-32}"

mkdir -p "${RUN_ROOT}" "${PIPELINE}/reports"

if [[ -e "${RUN_DIR}" ]]; then
    echo "ERROR: output directory already exists: ${RUN_DIR}" >&2
    exit 1
fi

for file in "${MANIFEST}" "${WORKFLOW}" "${DATABASE}"; do
    if [[ ! -f "${file}" || ! -r "${file}" ]]; then
        echo "ERROR: required file missing or unreadable: ${file}" >&2
        exit 1
    fi
done

module load ceuadmin/FragPipe/22.0

FRAGPIPE="$(command -v fragpipe)"
if [[ -z "${FRAGPIPE}" ]]; then
    echo "ERROR: fragpipe command unavailable after module load." >&2
    exit 1
fi

# Validate the manifest before creating a production run.
python - "${MANIFEST}" "${DATABASE}" "${WORKFLOW}" <<'PY' | tee "${REPORT}"
import os
import sys

manifest, database, workflow = sys.argv[1:]
rows = []
errors = []

with open(manifest, encoding="utf-8") as handle:
    for lineno, line in enumerate(handle, 1):
        line = line.rstrip("\r\n")
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) != 4:
            errors.append(f"Line {lineno}: expected 4 columns, got {len(fields)}")
            continue

        mzml, experiment, third, acquisition = fields
        rows.append((mzml, experiment, third, acquisition))

        if not os.path.isabs(mzml):
            errors.append(f"Line {lineno}: non-absolute mzML path: {mzml}")
        if not os.path.isfile(mzml):
            errors.append(f"Line {lineno}: mzML file missing: {mzml}")
        if third != "":
            errors.append(f"Line {lineno}: third column must be empty")
        if acquisition != "DDA":
            errors.append(f"Line {lineno}: expected DDA, got {acquisition!r}")

print("FRAGPIPE PRODUCTION PREFLIGHT")
print("=============================")
print(f"Manifest: {manifest}")
print(f"Workflow: {workflow}")
print(f"Database: {database}")
print(f"Manifest rows: {len(rows)}")

groups = {}
for _, exp, _, _ in rows:
    groups[exp] = groups.get(exp, 0) + 1

print("\nExperiment counts:")
for exp, count in sorted(groups.items()):
    print(f"  {exp}: {count}")

if len(rows) != 98:
    errors.append(f"Expected 98 rows; found {len(rows)}")
if len(groups) != 7:
    errors.append(f"Expected 7 experiment groups; found {len(groups)}")
if not os.path.isfile(database) or os.path.getsize(database) == 0:
    errors.append("Database missing or empty")
if not os.path.isfile(workflow) or os.path.getsize(workflow) == 0:
    errors.append("Workflow missing or empty")

print("\nErrors:")
for error in errors:
    print(f"  FAIL: {error}")

if errors:
    print("\nPREFLIGHT: FAIL")
    sys.exit(1)

print("  None")
print("\nPREFLIGHT: PASS")
PY

grep -Fq "database.db-path=" "${WORKFLOW}" || {
    echo "ERROR: database.db-path setting not found." >&2
    exit 1
}
grep -Fq "database.decoy-tag=rev_" "${WORKFLOW}" || {
    echo "ERROR: expected decoy tag rev_ not found." >&2
    exit 1
}

if ! [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer." >&2
    exit 1
fi
if ! [[ "${RAM_GB}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: FRAGPIPE_RAM_GB must be a non-negative integer." >&2
    exit 1
fi

echo
echo "FragPipe: ${FRAGPIPE}"
echo "Threads: ${THREADS}"
echo "RAM limit: ${RAM_GB} GB"
echo "Output: ${RUN_DIR}"
echo "Log: ${LOG}"

# Create the run directory only after preflight has passed.
mkdir "${RUN_DIR}"

printf 'Started: %s\n' "$(date --iso-8601=seconds)" > "${RUN_DIR}/run_metadata.txt"
printf 'Manifest: %s\nWorkflow: %s\nDatabase: %s\n' \
    "${MANIFEST}" "${WORKFLOW}" "${DATABASE}" >> "${RUN_DIR}/run_metadata.txt"
printf 'Threads: %s\nRAM_GB: %s\nFragPipe: %s\n' \
    "${THREADS}" "${RAM_GB}" "${FRAGPIPE}" >> "${RUN_DIR}/run_metadata.txt"

# Preserve a record of the effective workflow and manifest.
cp -p "${WORKFLOW}" "${RUN_DIR}/workflow_used.workflow"
cp -p "${MANIFEST}" "${RUN_DIR}/manifest_used.fp-manifest"

CMD=(
    "${FRAGPIPE}"
    --headless
    --workflow "${WORKFLOW}"
    --manifest "${MANIFEST}"
    --workdir "${RUN_DIR}"
    --threads "${THREADS}"
    --ram "${RAM_GB}"
)

printf '\nCommand:' | tee -a "${LOG}"
printf ' %q' "${CMD[@]}" | tee -a "${LOG}"
printf '\n\n' | tee -a "${LOG}"

set +e
"${CMD[@]}" 2>&1 | tee -a "${LOG}"
STATUS=${PIPESTATUS[0]}
set -e

printf '\nFinished: %s\nExit status: %s\n' \
    "$(date --iso-8601=seconds)" "${STATUS}" >> "${RUN_DIR}/run_metadata.txt"

if [[ "${STATUS}" -ne 0 ]]; then
    echo "FragPipe failed with exit status ${STATUS}." | tee -a "${LOG}"
    echo "Preserving output and logs for diagnosis: ${RUN_DIR}" | tee -a "${LOG}"
    exit "${STATUS}"
fi

echo "FragPipe completed with exit status 0." | tee -a "${LOG}"
echo "Review the output reports before interpreting identifications."
