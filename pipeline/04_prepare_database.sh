#!/usr/bin/env bash
# 04_prepare_database.sh
# Verify the existing target/decoy FASTA without modifying it.

set -eo pipefail

DB="${HOME}/Caprion/analysis/ZWK/discovery_fresh/database_canonical/2026-10-09-decoys-UP000005640-canonical-proteins-20420.fasta.fas"
REPORT_DIR="${HOME}/Caprion/analysis/pipeline/reports"
REPORT="${REPORT_DIR}/database_validation.txt"

mkdir -p "$REPORT_DIR"

if [[ ! -s "$DB" ]]; then
    echo "ERROR: Database is missing or empty: $DB" >&2
    exit 1
fi

python3 - "$DB" "$REPORT" <<'PY'
import sys
from pathlib import Path

db = Path(sys.argv[1])
report = Path(sys.argv[2])

targets = 0
decoys = 0
invalid_headers = 0
empty_sequences = 0
duplicate_headers = 0
seen = set()
current_header = None
current_length = 0
first_headers = []

def finish_sequence():
    global empty_sequences, current_length
    if current_header is not None and current_length == 0:
        empty_sequences += 1

with db.open("r", encoding="utf-8") as fh:
    for line_no, line in enumerate(fh, 1):
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            finish_sequence()
            current_header = line[1:].strip()
            current_length = 0

            if not current_header:
                invalid_headers += 1
                continue

            if current_header in seen:
                duplicate_headers += 1
            seen.add(current_header)

            first_headers.append(current_header) if len(first_headers) < 5 else None

            if current_header.startswith("rev_"):
                decoys += 1
            else:
                targets += 1
        else:
            if current_header is None:
                invalid_headers += 1
            else:
                current_length += len(line)

finish_sequence()

ok = (
    invalid_headers == 0
    and empty_sequences == 0
    and duplicate_headers == 0
    and targets == 20420
    and decoys == 20420
)

lines = [
    "Database validation report",
    f"File: {db}",
    f"Size: {db.stat().st_size:,} bytes",
    f"Target sequences: {targets}",
    f"Decoy sequences (rev_ prefix): {decoys}",
    f"Invalid/orphan sequence lines: {invalid_headers}",
    f"Empty sequences: {empty_sequences}",
    f"Duplicate headers: {duplicate_headers}",
    f"First headers: {first_headers}",
    f"Status: {'PASS' if ok else 'FAIL'}",
]

report.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))

if not ok:
    sys.exit(1)
PY

echo "Report: $REPORT"
echo "Database was not modified."
