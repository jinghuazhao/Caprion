#!/usr/bin/env bash
# 11_report_sample_qc.sh
#
# Purpose:
#   Summarise per-file PSM, peptide and protein identification counts from
#   the completed FragPipe production run, alongside experiment-level
#   Philosopher filtering summaries.
#
# Inputs:
#   production_runs/production_98_20261011_152549/
#
# Outputs:
#   reports/sample_qc_<timestamp>.tsv
#   reports/sample_qc_<timestamp>.txt
#
# Safety:
#   Read-only with respect to production outputs. Creates new reports only.
#   Does not rerun FragPipe or change any existing result files.
#
# Run:
#   bash 11_report_sample_qc.sh
#   bash 11_report_sample_qc.sh /path/to/production_run

set -eo pipefail

ROOT="$HOME/Caprion/analysis/pipeline"
RUN="${1:-$ROOT/production_runs/production_98_20261011_152549}"
REPORTS="$ROOT/reports"
STAMP="$(date +%Y%m%d_%H%M%S)"

mkdir -p "$REPORTS"

if [ ! -d "$RUN" ]; then
    echo "ERROR: run directory does not exist: $RUN" >&2
    exit 1
fi

TSV="$REPORTS/sample_qc_${STAMP}.tsv"
TXT="$REPORTS/sample_qc_${STAMP}.txt"

python - "$RUN" "$TSV" "$TXT" <<'PY'
import csv
import os
import re
import sys
from collections import defaultdict

run, tsv_path, txt_path = sys.argv[1:4]

experiments = [
    "IgY14_Pool4",
    "IgY14_Pool5",
    "IgY14_Pool6",
    "MARS14_Pool4",
    "MARS14_Pool5",
    "MARS14_Pool6",
    "PQC_3_QC",
]

def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames or [], rows

def unique_nonempty(rows, key):
    return len({r.get(key, "").strip() for r in rows
                if r.get(key, "").strip()})

def get_spectrum_file(value):
    value = (value or "").strip()
    if not value:
        return ""
    return os.path.basename(value)

def numeric(value):
    try:
        x = float(value)
        if x == x and abs(x) != float("inf"):
            return x
    except (TypeError, ValueError):
        pass
    return None

# Collect PSM, peptide and protein tables for each experiment.
records = []
warnings = []
totals = defaultdict(int)

for exp in experiments:
    exp_dir = os.path.join(run, exp)
    psm_path = os.path.join(exp_dir, "psm.tsv")
    peptide_path = os.path.join(exp_dir, "peptide.tsv")
    protein_path = os.path.join(exp_dir, "protein.tsv")

    if not os.path.isdir(exp_dir):
        warnings.append(f"Missing experiment directory: {exp}")
        continue

    tables = {}
    for label, path in (
        ("psm", psm_path),
        ("peptide", peptide_path),
        ("protein", protein_path),
    ):
        if not os.path.isfile(path):
            warnings.append(f"Missing {label} table: {path}")
            tables[label] = ([], [])
            continue
        try:
            tables[label] = read_tsv(path)
        except Exception as exc:
            warnings.append(f"Cannot read {path}: {exc}")
            tables[label] = ([], [])

    psm_header, psm_rows = tables["psm"]
    pep_header, pep_rows = tables["peptide"]
    prot_header, prot_rows = tables["protein"]

    # Per-file PSM counts use the actual Spectrum File field.
    psm_counts = defaultdict(int)
    file_probabilities = defaultdict(list)
    file_intensities = defaultdict(list)

    for row in psm_rows:
        file_name = get_spectrum_file(row.get("Spectrum File"))
        if not file_name:
            continue
        psm_counts[file_name] += 1

        p = numeric(row.get("Probability"))
        if p is not None:
            file_probabilities[file_name].append(p)

        intensity = numeric(row.get("Intensity"))
        if intensity is not None:
            file_intensities[file_name].append(intensity)

    # Per-file metrics only exist for files represented in psm.tsv.
    for file_name in sorted(psm_counts):
        probs = file_probabilities[file_name]
        intensities = file_intensities[file_name]
        records.append({
            "Experiment": exp,
            "Spectrum_File": file_name,
            "PSM_Count": psm_counts[file_name],
            "Peptide_Table_Rows": len(pep_rows),
            "Protein_Table_Rows": len(prot_rows),
            "Distinct_PSM_Peptides": unique_nonempty(psm_rows, "Peptide"),
            "Distinct_PSM_Proteins": unique_nonempty(psm_rows, "Protein ID"),
            "Mean_PSM_Probability": (
                sum(probs) / len(probs) if probs else ""
            ),
            "Median_PSM_Probability": (
                sorted(probs)[len(probs)//2] if probs else ""
            ),
            "Nonzero_PSM_Intensity_Rows": sum(x > 0 for x in intensities),
            "PSM_Intensity_Rows": len(intensities),
            "QC_Note": "PSM probability is not an FDR estimate",
        })

    totals[exp + "|PSM"] = len(psm_rows)
    totals[exp + "|Peptide"] = len(pep_rows)
    totals[exp + "|Protein"] = len(prot_rows)

fields = [
    "Experiment", "Spectrum_File", "PSM_Count",
    "Peptide_Table_Rows", "Protein_Table_Rows",
    "Distinct_PSM_Peptides", "Distinct_PSM_Proteins",
    "Mean_PSM_Probability", "Median_PSM_Probability",
    "Nonzero_PSM_Intensity_Rows", "PSM_Intensity_Rows", "QC_Note",
]

with open(tsv_path, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for row in records:
        writer.writerow(row)

lines = [
    "FRAGPIPE SAMPLE QC SUMMARY",
    "=" * 80,
    f"Run directory: {run}",
    f"Per-file rows: {len(records)}",
    "",
    "Experiment-level filtered table counts",
    "-" * 80,
    "Experiment\tPSMs\tPeptides\tProteins",
]

for exp in experiments:
    lines.append(
        f"{exp}\t{totals.get(exp+'|PSM', 0)}\t"
        f"{totals.get(exp+'|Peptide', 0)}\t"
        f"{totals.get(exp+'|Protein', 0)}"
    )

lines += [
    "",
    "Interpretation notes",
    "-" * 80,
    "PSM counts are assigned using Spectrum File in each experiment's psm.tsv.",
    "Peptide/protein table row counts are experiment-level, not per-file counts.",
    "Distinct PSM peptides/proteins are calculated from the observed PSM table.",
    "PSM probability is not itself a PSM FDR or q-value.",
    "Per-experiment filtering FDR is recorded in each experiment's filter.log.",
    "Quantification and match-between-runs FDR require separate assessment.",
    "",
    "Warnings",
    "-" * 80,
]
lines.extend(warnings or ["None detected by this script."])
lines += ["", f"Per-file table: {tsv_path}"]

with open(txt_path, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")

print("\n".join(lines))
print(f"\nQC TSV written: {tsv_path}")
PY
