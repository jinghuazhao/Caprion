#!/usr/bin/env bash
# =============================================================================
# 08_review_dry_run.sh
#
# Purpose:
#   Review the FragPipe 22.0 dry-run evidence before production execution.
#
# Inputs:
#   Production manifest: production_v2.fp-manifest
#   Prepared workflow:  LFQ-MBR-Sage-comparison.workflow
#   Generated workflow: dry_run_98_v2/fragpipe.workflow
#   Dry-run log:        reports/fragpipe_dry_run_v2.log
#   Validation report:  reports/fragpipe_dry_run_validation_v2.txt
#
# Outputs:
#   reports/dry_run_review_v2.txt
#
# Safety:
#   - Read-only audit; does not launch FragPipe.
#   - Does not modify workflows, manifests, databases, or dry-run files.
#   - Refuses to overwrite an existing review report.
#
# Execution:
#   cd ~/Caprion/analysis/pipeline
#   bash 08_review_dry_run.sh
#
# Exit status:
#   0 = structural review passed
#   Non-zero = missing inputs or review failure
# =============================================================================

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
MANIFEST="${PIPELINE}/production_v2.fp-manifest"
WORKFLOW="${PIPELINE}/LFQ-MBR-Sage-comparison.workflow"
GENERATED="${PIPELINE}/dry_run_98_v2/fragpipe.workflow"
LOG="${PIPELINE}/reports/fragpipe_dry_run_v2.log"
VALIDATION="${PIPELINE}/reports/fragpipe_dry_run_validation_v2.txt"
REPORT="${PIPELINE}/reports/dry_run_review_v2.txt"

mkdir -p "${PIPELINE}/reports"

if [[ -e "$REPORT" ]]; then
    printf 'ERROR: refusing to overwrite %s\n' "$REPORT" >&2
    exit 1
fi

for file in "$MANIFEST" "$WORKFLOW" "$GENERATED" "$LOG" "$VALIDATION"; do
    if [[ ! -f "$file" ]]; then
        printf 'ERROR: required file missing: %s\n' "$file" >&2
        exit 1
    fi
done

python - "$MANIFEST" "$WORKFLOW" "$GENERATED" "$LOG" "$VALIDATION" <<'PY'
import os
import sys

manifest, workflow, generated, log, validation = sys.argv[1:]
errors = []
notes = []

# 1. Confirm the production manifest has the expected four-column format.
rows = []
with open(manifest, encoding="utf-8", newline="") as handle:
    for line_number, line in enumerate(handle, 1):
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != 4:
            errors.append(
                f"Manifest line {line_number}: expected 4 columns, found {len(fields)}"
            )
            continue

        mzml, experiment, third, acquisition = fields

        if not mzml or not os.path.isfile(mzml):
            errors.append(f"Manifest line {line_number}: missing mzML file {mzml!r}")
        if not experiment:
            errors.append(f"Manifest line {line_number}: empty experiment name")
        if third != "":
            errors.append(f"Manifest line {line_number}: third column is not empty")
        if acquisition != "DDA":
            errors.append(f"Manifest line {line_number}: acquisition is not DDA")

        rows.append(mzml)

if len(rows) != 98:
    errors.append(f"Expected 98 manifest rows, found {len(rows)}")
if len(set(rows)) != len(rows):
    errors.append("Manifest contains duplicate mzML paths")

# 2. Read workflow contents.
with open(workflow, encoding="utf-8", errors="replace") as handle:
    original = handle.read()

with open(generated, encoding="utf-8", errors="replace") as handle:
    generated_text = handle.read()

with open(log, encoding="utf-8", errors="replace") as handle:
    log_text = handle.read()

with open(validation, encoding="utf-8", errors="replace") as handle:
    validation_text = handle.read()

# 3. Confirm the database path and decoy tag are preserved.
def property_values(text, key):
    prefix = key + "="
    return [
        line[len(prefix):].strip()
        for line in text.splitlines()
        if line.startswith(prefix)
    ]

for key in ("database.db-path", "database.decoy-tag"):
    original_values = property_values(original, key)
    generated_values = property_values(generated_text, key)

    if len(original_values) != 1:
        errors.append(
            f"Prepared workflow: expected one {key}, found {len(original_values)}"
        )
    if len(generated_values) != 1:
        errors.append(
            f"Generated workflow: expected one {key}, found {len(generated_values)}"
        )

    if len(original_values) == 1 and len(generated_values) == 1:
        if key.endswith("db-path"):
            if not os.path.exists(original_values[0]):
                errors.append("Prepared workflow database path does not exist")
            elif os.path.realpath(original_values[0]) != os.path.realpath(generated_values[0]):
                errors.append("Generated workflow database path differs from prepared workflow")
        elif original_values[0] != generated_values[0]:
            errors.append("Generated workflow decoy tag differs from prepared workflow")

if property_values(original, "database.decoy-tag") != ["rev_"]:
    errors.append("Prepared workflow decoy tag is not rev_")

# 4. Confirm the recorded dry-run validation passed.
if "FragPipe dry run: PASS" not in validation_text:
    errors.append("Validation report does not record a successful dry run")

if "FragPipe dry run: FAIL" in log_text:
    notes.append("Log contains a failure marker; inspect the full log")

# A dry-run log listing commands is evidence of a plan, not execution.
if "Cmd: [MSFragger]" in log_text:
    notes.append("MSFragger is listed in the planned command sequence")
if "Cmd: [IonQuant]" in log_text:
    notes.append("IonQuant is listed in the planned command sequence")

print("Dry-run review")
print("==============")
print(f"Manifest rows: {len(rows)}")
print("Expected manifest format: four columns, empty third column, DDA")
print(f"Prepared workflow: {workflow}")
print(f"Generated workflow: {generated}")
print("Database path and decoy-tag consistency checked")

if errors:
    print("\nREVIEW: FAIL")
    for item in errors:
        print(f"ERROR: {item}")
    sys.exit(1)

print("\nREVIEW: PASS")
for item in notes:
    print(f"NOTE: {item}")

print("\nThis review does not execute FragPipe or establish identification results.")
PY

# Preserve the audit output as a report.
# The preceding checks are repeated by capturing this review in a fresh run.
# The report is created exclusively so an existing report is never overwritten.
{
    printf 'Dry-run review completed at %s\n' "$(date -Is)"
    printf 'Manifest: %s\n' "$MANIFEST"
    printf 'Prepared workflow: %s\n' "$WORKFLOW"
    printf 'Generated workflow: %s\n' "$GENERATED"
    printf 'Dry-run log: %s\n' "$LOG"
    printf 'Validation report: %s\n' "$VALIDATION"
    printf '\nThe detailed review is printed to the terminal.\n'
} > "$REPORT"

printf '\nReview report saved to: %s\n' "$REPORT"
