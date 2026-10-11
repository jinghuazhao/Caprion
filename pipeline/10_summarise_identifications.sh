#!/usr/bin/env bash
# =============================================================================
# 10_summarise_identifications.sh
#
# Purpose:
#   Summarise FragPipe PSM, peptide, ion and protein identification outputs.
#
# Inputs:
#   production_runs/production_98_20261011_152549/
#
# Outputs:
#   reports/identification_summary_<timestamp>.txt
#   reports/identification_table_inventory_<timestamp>.tsv
#
# Safety:
#   Read-only. Does not modify or delete production outputs.
#
# Execution:
#   bash 10_summarise_identifications.sh
#   bash 10_summarise_identifications.sh production_runs/production_98_20261011_152549
# =============================================================================

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
RUN_ROOT="${PIPELINE}/production_runs"
REPORT_DIR="${PIPELINE}/reports"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

mkdir -p "${REPORT_DIR}"

if [[ $# -gt 0 ]]; then
    RUN_DIR="$1"
    [[ "${RUN_DIR}" = /* ]] || RUN_DIR="${PIPELINE}/${RUN_DIR}"
else
    RUN_DIR="${RUN_ROOT}/production_98_20261011_152549"
fi

if [[ ! -d "${RUN_DIR}" ]]; then
    echo "ERROR: production directory not found: ${RUN_DIR}" >&2
    exit 1
fi

SUMMARY="${REPORT_DIR}/identification_summary_${TIMESTAMP}.txt"
INVENTORY="${REPORT_DIR}/identification_table_inventory_${TIMESTAMP}.tsv"

python - "${RUN_DIR}" "${SUMMARY}" "${INVENTORY}" <<'PY'
import csv
import os
import sys
from collections import Counter

run_dir, summary_path, inventory_path = sys.argv[1:]

combined_names = [
    "combined_ion.tsv",
    "combined_peptide.tsv",
    "combined_modified_peptide.tsv",
    "combined_protein.tsv",
]

per_experiment_names = ["psm.tsv", "peptide.tsv", "protein.tsv"]

summary = []
inventory = []

def note(message=""):
    summary.append(message)

def read_table(path):
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        rows = list(reader)
    return headers, rows

def column_kind(header):
    lower = header.lower()
    if any(term in lower for term in ("q-value", "qvalue", "q_value", "fdr")):
        return "FDR/q-value"
    if "probability" in lower or lower.endswith("_prob"):
        return "probability"
    if "intensity" in lower or "abundance" in lower:
        return "intensity"
    return ""

def summarize_table(path, label):
    if not os.path.isfile(path):
        note(f"WARNING: missing table: {path}")
        inventory.append([label, path, "MISSING", 0, ""])
        return

    try:
        headers, rows = read_table(path)
    except Exception as exc:
        note(f"WARNING: cannot read {path}: {exc}")
        inventory.append([label, path, "READ_ERROR", 0, ""])
        return

    inventory.append([
        label, path, "OK", len(rows), "; ".join(headers)
    ])

    note(f"\n{label}")
    note("-" * len(label))
    note(f"File: {path}")
    note(f"Data rows: {len(rows)}")
    note(f"Columns: {len(headers)}")
    note("Headers: " + "; ".join(headers))

    special = [(h, column_kind(h)) for h in headers if column_kind(h)]
    if special:
        note("Potential QC columns:")
        for header, kind in special:
            note(f"  {header}: {kind}")

        for header, kind in special:
            values = []
            for row in rows:
                value = row.get(header, "")
                try:
                    number = float(value)
                except (TypeError, ValueError):
                    continue
                if number == number and abs(number) != float("inf"):
                    values.append(number)

            if values:
                note(
                    f"  Numeric range for {header}: "
                    f"min={min(values):.6g}, max={max(values):.6g}"
                )
    else:
        note("No obvious FDR/q-value, probability or intensity column detected.")

note("FRAGPIPE IDENTIFICATION SUMMARY")
note("===============================")
note(f"Run directory: {run_dir}")
note("Read-only summary; no production data were modified.")

note("\n1. COMBINED TABLES")
note("==================")
for name in combined_names:
    summarize_table(os.path.join(run_dir, name), name)

note("\n2. PER-EXPERIMENT TABLES")
note("========================")
experiment_dirs = sorted(
    entry for entry in os.listdir(run_dir)
    if os.path.isdir(os.path.join(run_dir, entry))
    and not entry.startswith(".")
)

for experiment in experiment_dirs:
    folder = os.path.join(run_dir, experiment)
    available = [name for name in per_experiment_names
                 if os.path.isfile(os.path.join(folder, name))]
    if not available:
        continue

    note(f"\nEXPERIMENT: {experiment}")
    for name in per_experiment_names:
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            summarize_table(path, f"{experiment}/{name}")

note("\n3. FDR INTERPRETATION")
note("=====================")
note(
    "This script reports table dimensions, headers and numeric ranges. "
    "It does not certify FDR control."
)
note(
    "Inspect the actual q-value/FDR columns and the corresponding Philosopher "
    "and Percolator reports before claiming 1% PSM, peptide or protein FDR."
)
note(
    "The production log reported peptide FDR=1.0 and protein FDR=1.0. "
    "This must be resolved before treating protein identifications as "
    "filtered at 1% FDR."
)

with open(summary_path, "w", encoding="utf-8") as handle:
    handle.write("\n".join(summary) + "\n")

with open(inventory_path, "w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["table", "path", "status", "data_rows", "headers"])
    writer.writerows(inventory)

print("\n".join(summary))
print(f"\nSummary report: {summary_path}")
print(f"Table inventory: {inventory_path}")
PY
