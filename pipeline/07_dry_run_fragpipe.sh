#!/usr/bin/env bash
# 07_dry_run_fragpipe.sh
#
# Purpose:
#   Validate the 98-file production manifest and prepared FragPipe workflow,
#   then execute FragPipe 22.0 in headless dry-run mode.
#
# Inputs:
#   Workflow:  LFQ-MBR-Sage-comparison.workflow
#   Manifest:  production_v2.fp-manifest
#   Metadata:  production_sample_metadata_v2.tsv
#   Database:  validated canonical UniProt FASTA with reverse decoys
#
# Outputs:
#   Workdir:   dry_run_98_v2/
#   Log:       reports/fragpipe_dry_run_v2.log
#   Report:    reports/fragpipe_dry_run_validation_v2.txt
#
# Prerequisites:
#   - FragPipe 22.0 module and dependencies are available.
#   - Java 11 or newer is active.
#   - Manifest contains 98 rows and four tab-separated fields.
#   - The third field is empty and the fourth field is DDA.
#   - All referenced mzML files exist.
#   - Workflow references the validated FASTA and decoy tag rev_.
#
# Safety:
#   - Does not modify the source workflow, prepared workflow, manifest,
#     metadata, FASTA database, or successful six-file pilot.
#   - Refuses to overwrite the versioned workdir, log, or report.
#   - Runs in dry-run mode; this is not the full identification analysis.
#
# Execution:
#   cd ~/Caprion/analysis/pipeline
#   bash 07_dry_run_fragpipe.sh
#
# Exit status:
#   0 = dry run and validation passed
#   Non-zero = preflight, FragPipe execution, or validation failed
# =============================================================================

set -eo pipefail

module load ceuadmin/FragPipe/22.0

PIPELINE="${HOME}/Caprion/analysis/pipeline"
FP="/usr/local/Cluster-Apps/ceuadmin/FragPipe/22.0"

WORKFLOW="${PIPELINE}/LFQ-MBR-Sage-comparison.workflow"
MANIFEST="${PIPELINE}/production_v2.fp-manifest"
METADATA="${PIPELINE}/production_sample_metadata_v2.tsv"

DATABASE="${HOME}/Caprion/analysis/ZWK/discovery_fresh/database_canonical/2026-10-09-decoys-UP000005640-canonical-proteins-20420.fasta.fas"

DRYRUN_DIR="${PIPELINE}/dry_run_98_v2"
JAVA_HOME_DIR="${PIPELINE}/fragpipe_test_home"
REPORT_DIR="${PIPELINE}/reports"

LOG="${REPORT_DIR}/fragpipe_dry_run_v2.log"
REPORT="${REPORT_DIR}/fragpipe_dry_run_validation_v2.txt"

mkdir -p "$REPORT_DIR"

# Refuse to overwrite evidence or output from a previous attempt.
for path in "$DRYRUN_DIR" "$LOG" "$REPORT"; do
    if [[ -e "$path" ]]; then
        printf 'ERROR: refusing to overwrite existing path:\n%s\n' "$path" >&2
        exit 1
    fi
done

# Check required files.
for path in \
    "$FP/bin/fragpipe" \
    "$WORKFLOW" \
    "$MANIFEST" \
    "$METADATA" \
    "$DATABASE"
do
    if [[ ! -f "$path" ]]; then
        printf 'ERROR: required file does not exist:\n%s\n' "$path" >&2
        exit 1
    fi
done

# FragPipe 22.0 requires Java 11 or newer.
JAVA_VERSION="$(java -version 2>&1 | head -n 1)"
printf 'Java: %s\n' "$JAVA_VERSION"

JAVA_MAJOR="$(
    java -version 2>&1 |
    awk -F '"' '/version/ {
        split($2, v, ".");
        if (v[1] == "1") print v[2];
        else print v[1];
        exit
    }'
)"

if [[ -z "$JAVA_MAJOR" ]] || (( JAVA_MAJOR < 11 )); then
    printf 'ERROR: FragPipe 22.0 requires Java 11 or newer.\n' >&2
    exit 1
fi

# Validate manifest structure, acquisition type, paths, and uniqueness.
python - "$MANIFEST" <<'PY'
import os
import sys

manifest = sys.argv[1]
rows = []
errors = []

with open(manifest, encoding="utf-8", newline="") as handle:
    for line_number, line in enumerate(handle, start=1):
        line = line.rstrip("\r\n")
        fields = line.split("\t")

        if len(fields) != 4:
            errors.append(
                f"Line {line_number}: expected 4 fields, found {len(fields)}"
            )
            continue

        mzml, experiment, third, acquisition = fields

        if not mzml:
            errors.append(f"Line {line_number}: empty mzML path")
        elif not os.path.isfile(mzml):
            errors.append(f"Line {line_number}: missing mzML file: {mzml}")

        if not experiment:
            errors.append(f"Line {line_number}: empty experiment name")

        if third != "":
            errors.append(
                f"Line {line_number}: third field must be empty"
            )

        if acquisition != "DDA":
            errors.append(
                f"Line {line_number}: acquisition type must be DDA"
            )

        rows.append(mzml)

if len(rows) != 98:
    errors.append(f"Expected 98 manifest rows, found {len(rows)}")

if len(set(rows)) != len(rows):
    errors.append("Manifest contains duplicate mzML paths")

if errors:
    print("MANIFEST VALIDATION: FAIL")
    for error in errors:
        print(f"  - {error}")
    sys.exit(1)

print("MANIFEST VALIDATION: PASS")
print("Rows: 98")
print("Columns: 4")
print("Third column: empty")
print("Acquisition type: DDA")
print("All mzML paths exist and are unique")
PY

# Validate workflow database path using canonical paths.
# $HOME resolves to the RDS-mounted path on this cluster.
python - "$WORKFLOW" "$DATABASE" <<'PY'
import os
import sys

workflow, database = sys.argv[1:]
database = os.path.realpath(database)

with open(workflow, encoding="utf-8", errors="replace") as handle:
    content = handle.read()

entries = [
    line.split("=", 1)[1].strip()
    for line in content.splitlines()
    if line.startswith("database.db-path=")
]

if len(entries) != 1:
    print(
        "ERROR: expected exactly one database.db-path entry; "
        f"found {len(entries)}."
    )
    sys.exit(1)

workflow_database = os.path.realpath(entries[0])

if workflow_database != database:
    print("ERROR: workflow database path does not resolve to the validated FASTA.")
    print(f"Workflow path: {entries[0]}")
    print(f"Resolved path: {workflow_database}")
    print(f"Expected path: {database}")
    sys.exit(1)

if "database.decoy-tag=rev_" not in content.splitlines():
    print("ERROR: workflow decoy tag is not rev_.")
    sys.exit(1)

print("WORKFLOW DATABASE VALIDATION: PASS")
print(f"Workflow path: {entries[0]}")
print(f"Resolved FASTA: {database}")
print("Decoy tag: rev_")
PY

# Run FragPipe in a new, isolated dry-run directory.
mkdir "$DRYRUN_DIR"

printf 'Starting FragPipe dry run at %s\n' "$(date -Is)" | tee "$LOG"
printf 'Workflow: %s\n' "$WORKFLOW" | tee -a "$LOG"
printf 'Manifest: %s\n' "$MANIFEST" | tee -a "$LOG"
printf 'Workdir: %s\n' "$DRYRUN_DIR" | tee -a "$LOG"

set +e
JAVA_OPTS="-Duser.home=$JAVA_HOME_DIR" \
    "$FP/bin/fragpipe" \
    --headless \
    --workflow "$WORKFLOW" \
    --manifest "$MANIFEST" \
    --workdir "$DRYRUN_DIR" \
    --threads 8 \
    --ram 24 \
    --dry-run 2>&1 | tee -a "$LOG"

FP_STATUS=${PIPESTATUS[0]}
set -e

if (( FP_STATUS != 0 )); then
    {
        printf 'FragPipe dry run: FAIL\n'
        printf 'Exit status: %s\n' "$FP_STATUS"
        printf 'Log: %s\n' "$LOG"
        printf 'Workdir: %s\n' "$DRYRUN_DIR"
    } | tee "$REPORT"
    exit "$FP_STATUS"
fi

# Validate the workflow FragPipe generated.
GENERATED_WORKFLOW="${DRYRUN_DIR}/fragpipe.workflow"

if [[ ! -f "$GENERATED_WORKFLOW" ]]; then
    {
        printf 'FragPipe dry run: FAIL\n'
        printf 'Generated workflow not found: %s\n' "$GENERATED_WORKFLOW"
        printf 'Log: %s\n' "$LOG"
    } | tee "$REPORT"
    exit 1
fi

python - "$GENERATED_WORKFLOW" "$DATABASE" <<'PY'
import os
import sys

workflow, database = sys.argv[1:]
database = os.path.realpath(database)

with open(workflow, encoding="utf-8", errors="replace") as handle:
    content = handle.read()

entries = [
    line.split("=", 1)[1].strip()
    for line in content.splitlines()
    if line.startswith("database.db-path=")
]

if len(entries) != 1:
    print("GENERATED WORKFLOW VALIDATION: FAIL")
    print(f"Expected one database.db-path entry; found {len(entries)}")
    sys.exit(1)

if os.path.realpath(entries[0]) != database:
    print("GENERATED WORKFLOW VALIDATION: FAIL")
    print(f"Generated path: {entries[0]}")
    print(f"Expected path: {database}")
    sys.exit(1)

if "database.decoy-tag=rev_" not in content.splitlines():
    print("GENERATED WORKFLOW VALIDATION: FAIL")
    print("Missing database.decoy-tag=rev_")
    sys.exit(1)

print("GENERATED WORKFLOW VALIDATION: PASS")
print(f"Resolved FASTA: {database}")
print("Decoy tag: rev_")
PY

{
    printf 'FragPipe dry run: PASS\n'
    printf 'Java: %s\n' "$JAVA_VERSION"
    printf 'Manifest: %s\n' "$MANIFEST"
    printf 'Manifest rows: 98\n'
    printf 'Manifest columns: 4 (third field empty)\n'
    printf 'Workflow: %s\n' "$WORKFLOW"
    printf 'Database: %s\n' "$(readlink -f "$DATABASE")"
    printf 'Decoy tag: rev_\n'
    printf 'Dry-run directory: %s\n' "$DRYRUN_DIR"
    printf 'Log: %s\n' "$LOG"
} | tee "$REPORT"

printf '\nDry-run validation completed successfully.\n'
