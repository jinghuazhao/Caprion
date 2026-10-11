#!/usr/bin/env bash
# 05b_set_workflow_database.sh
# Set the validated FASTA path in the prepared workflow copy.
# Preserve the source workflow and successful pilot outputs unchanged.

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
WORKFLOW="${PIPELINE}/LFQ-MBR-Sage-comparison.workflow"
DATABASE="${HOME}/Caprion/analysis/ZWK/discovery_fresh/database_canonical/2026-10-09-decoys-UP000005640-canonical-proteins-20420.fasta.fas"
BACKUP="${WORKFLOW}.bak"

if [[ ! -f "$WORKFLOW" ]]; then
    echo "ERROR: Prepared workflow not found: $WORKFLOW" >&2
    exit 1
fi

if [[ ! -f "$DATABASE" ]]; then
    echo "ERROR: Validated FASTA not found: $DATABASE" >&2
    exit 1
fi

if [[ -e "$BACKUP" ]]; then
    echo "ERROR: Backup already exists: $BACKUP" >&2
    echo "Inspect it before deciding how to proceed." >&2
    exit 1
fi

python3 - "$WORKFLOW" "$DATABASE" "$BACKUP" <<'PY'
from pathlib import Path
import os
import sys
import tempfile

workflow = Path(sys.argv[1])
database = Path(sys.argv[2]).resolve()
backup = Path(sys.argv[3])

original = workflow.read_bytes()

# Preserve the original line-ending style.
newline = b"\r\n" if b"\r\n" in original else b"\n"
text = original.decode("utf-8")
lines = text.splitlines()

db_key = "database.db-path="
decoy_key = "database.decoy-tag="

db_indices = [i for i, line in enumerate(lines) if line.startswith(db_key)]
decoy_indices = [i for i, line in enumerate(lines) if line.startswith(decoy_key)]

if len(db_indices) > 1:
    sys.exit("ERROR: Multiple database.db-path settings; refusing to edit.")
if len(decoy_indices) != 1:
    sys.exit("ERROR: Expected exactly one database.decoy-tag setting.")

if lines[decoy_indices[0]].strip() != "database.decoy-tag=rev_":
    sys.exit("ERROR: Unexpected decoy tag; refusing to edit.")

if db_indices:
    lines[db_indices[0]] = db_key + str(database)
else:
    # Add the missing database setting directly after the decoy-tag setting.
    lines.insert(decoy_indices[0] + 1, db_key + str(database))

updated = (newline.decode().join(lines) + newline.decode()).encode("utf-8")

# Create an exclusive backup before replacing the prepared copy.
with backup.open("xb") as f:
    f.write(original)

fd, tmp_name = tempfile.mkstemp(
    prefix=workflow.name + ".",
    suffix=".tmp",
    dir=str(workflow.parent)
)
try:
    with os.fdopen(fd, "wb") as f:
        f.write(updated)

    # Verify the temporary copy before replacing the prepared workflow.
    check = Path(tmp_name).read_text().splitlines()
    db_lines = [line for line in check if line.startswith(db_key)]
    decoy_lines = [line for line in check if line.startswith(decoy_key)]

    if db_lines != [db_key + str(database)]:
        raise RuntimeError("Database-path verification failed.")
    if decoy_lines != ["database.decoy-tag=rev_"]:
        raise RuntimeError("Decoy-tag verification failed.")

    os.replace(tmp_name, workflow)
except Exception:
    Path(tmp_name).unlink(missing_ok=True)
    raise

print("Workflow database configuration updated.")
print(f"Workflow: {workflow}")
print(f"Backup:   {backup}")
print(f"Database: {database}")
print("Decoy tag: rev_")
print("Status: PASS")
PY

echo
grep -nE '^database\.db-path=|^database\.decoy-tag=' "$WORKFLOW"
