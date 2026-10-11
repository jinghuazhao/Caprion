#!/usr/bin/env bash
# 12_summarise_batch.sh
#
# Purpose:
#   Produce a consolidated, auditable summary of the completed Caprion
#   98-file FragPipe production batch.
#
# Inputs:
#   production_runs/production_98_20261011_152549/
#   reports/fragpipe_production_20261011_152549.log
#   reports/results_validation_20261011_170819.txt
#   reports/identification_summary_20261011_171039.txt
#   reports/sample_qc_20261011_171510.tsv
#
# Outputs:
#   reports/batch_summary_<timestamp>.txt
#   reports/batch_experiment_counts_<timestamp>.tsv
#   reports/batch_sample_qc_<timestamp>.tsv
#
# Safety:
#   Reads existing outputs only. Does not rerun FragPipe, alter workflow
#   settings, overwrite previous reports, or modify production result files.
#
# Run:
#   bash 12_summarise_batch.sh
#   bash 12_summarise_batch.sh /path/to/production_run
#
# Requirements:
#   Bash, Python 3.11 (or compatible Python 3), standard library only.
#
# Important interpretation:
#   Philosopher's final filtered counts and TSV row counts are reported
#   separately because they may represent different processing stages.
#   PSM probability is not an FDR estimate.
#   IonQuant quantification and match-between-runs settings are reviewed
#   separately from Philosopher identification filtering.

set -eo pipefail

ROOT="$HOME/Caprion/analysis/pipeline"
RUN="${1:-$ROOT/production_runs/production_98_20261011_152549}"
REPORTS="$ROOT/reports"
STAMP="$(date +%Y%m%d_%H%M%S)"

mkdir -p "$REPORTS"

if [ ! -d "$RUN" ]; then
    echo "ERROR: production run directory not found: $RUN" >&2
    exit 1
fi

PROD_LOG="$REPORTS/fragpipe_production_20261011_152549.log"
VALIDATION="$REPORTS/results_validation_20261011_170819.txt"
IDENT_SUMMARY="$REPORTS/identification_summary_20261011_171039.txt"
SAMPLE_QC="$REPORTS/sample_qc_20261011_171510.tsv"

OUT_TXT="$REPORTS/batch_summary_${STAMP}.txt"
OUT_EXP="$REPORTS/batch_experiment_counts_${STAMP}.tsv"
OUT_SAMPLE="$REPORTS/batch_sample_qc_${STAMP}.tsv"

python - "$ROOT" "$RUN" "$PROD_LOG" "$VALIDATION" \
    "$IDENT_SUMMARY" "$SAMPLE_QC" "$OUT_TXT" "$OUT_EXP" "$OUT_SAMPLE" <<'PY'
import csv
import os
import re
import sys
from collections import Counter

(
    root, run, prod_log, validation, ident_summary, sample_qc,
    out_txt, out_exp, out_sample
) = sys.argv[1:10]

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
        return reader.fieldnames or [], list(reader)

def write_tsv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

def read_text(path):
    if not os.path.isfile(path):
        return ""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()

def count_data_rows(path):
    if not os.path.isfile(path):
        return None
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return max(0, sum(1 for _ in fh) - 1)

def number(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

warnings = []
lines = []

# ---------------------------------------------------------------------
# Production outputs and per-experiment TSV row counts.
# ---------------------------------------------------------------------
experiment_rows = []
experiment_fdr = {}

for exp in experiments:
    exp_dir = os.path.join(run, exp)
    psm_path = os.path.join(exp_dir, "psm.tsv")
    peptide_path = os.path.join(exp_dir, "peptide.tsv")
    protein_path = os.path.join(exp_dir, "protein.tsv")
    filter_path = os.path.join(exp_dir, "filter.log")

    psm_rows = count_data_rows(psm_path)
    peptide_rows = count_data_rows(peptide_path)
    protein_rows = count_data_rows(protein_path)

    if not os.path.isdir(exp_dir):
        warnings.append(f"Missing experiment directory: {exp_dir}")
    for label, path, value in (
        ("PSM", psm_path, psm_rows),
        ("peptide", peptide_path, peptide_rows),
        ("protein", protein_path, protein_rows),
    ):
        if value is None:
            warnings.append(f"Missing {label} table: {path}")

    fdr_lines = []
    if os.path.isfile(filter_path):
        for line in read_text(filter_path).splitlines():
            if (
                "Final report numbers after FDR filtering" in line
                or "Converged to" in line
            ):
                fdr_lines.append(line.strip())
    else:
        warnings.append(f"Missing filter log: {filter_path}")

    final_counts = {}
    for line in fdr_lines:
        if "Final report numbers after FDR filtering" in line:
            for key in ("ions", "peptides", "proteins", "psms"):
                match = re.search(rf"\b{key}=(\d+)\b", line)
                if match:
                    final_counts[key] = int(match.group(1))

    row = {
        "Experiment": exp,
        "PSM_TSV_Rows": "" if psm_rows is None else psm_rows,
        "Peptide_TSV_Rows": "" if peptide_rows is None else peptide_rows,
        "Protein_TSV_Rows": "" if protein_rows is None else protein_rows,
        "Philosopher_Final_PSMs": final_counts.get("psms", ""),
        "Philosopher_Final_Peptides": final_counts.get("peptides", ""),
        "Philosopher_Final_Proteins": final_counts.get("proteins", ""),
        "Philosopher_Final_Ions": final_counts.get("ions", ""),
        "Filter_Log": os.path.relpath(filter_path, run),
    }
    experiment_rows.append(row)
    experiment_fdr[exp] = fdr_lines

exp_fields = [
    "Experiment",
    "PSM_TSV_Rows",
    "Peptide_TSV_Rows",
    "Protein_TSV_Rows",
    "Philosopher_Final_PSMs",
    "Philosopher_Final_Peptides",
    "Philosopher_Final_Proteins",
    "Philosopher_Final_Ions",
    "Filter_Log",
]
write_tsv(out_exp, exp_fields, experiment_rows)

# ---------------------------------------------------------------------
# Sample QC table: preserve original columns and add a filename caution.
# ---------------------------------------------------------------------
sample_rows = []
if os.path.isfile(sample_qc):
    sample_fields, sample_rows = read_tsv(sample_qc)
    for row in sample_rows:
        spectrum = row.get("Spectrum_File", "")
        if spectrum.lower().endswith(".pep.xml"):
            row["Filename_Interpretation"] = (
                "PSM Spectrum_File is a pepXML filename; map to original "
                "mzML via manifest/workflow metadata before treating it "
                "as the raw-file name."
            )
        else:
            row["Filename_Interpretation"] = (
                "Confirm mapping to the original mzML manifest if needed."
            )
    sample_fields = sample_fields + ["Filename_Interpretation"]
    write_tsv(out_sample, sample_fields, sample_rows)
else:
    warnings.append(f"Missing sample QC TSV: {sample_qc}")
    sample_fields = []
    sample_rows = []

psm_values = [
    number(row.get("PSM_Count"))
    for row in sample_rows
    if number(row.get("PSM_Count")) is not None
]
file_counts = Counter(row.get("Experiment", "") for row in sample_rows)
sample_by_exp = Counter()
for row in sample_rows:
    sample_by_exp[row.get("Experiment", "")] += 1

# ---------------------------------------------------------------------
# Production log and auxiliary reports.
# ---------------------------------------------------------------------
prod_text = read_text(prod_log)
validation_text = read_text(validation)
ident_text = read_text(ident_summary)

if not prod_text:
    warnings.append(f"Production log not found/readable: {prod_log}")
if not validation_text:
    warnings.append(f"Validation report not found/readable: {validation}")
if not ident_text:
    warnings.append(f"Identification summary not found/readable: {ident_summary}")

# Look for high-level production completion evidence without assuming
# that a generic word such as "complete" proves every analysis stage.
completion_lines = [
    line.strip() for line in prod_text.splitlines()
    if re.search(
        r"all jobs done|ionquant.*exit code|finished successfully|"
        r"completed successfully",
        line,
        re.IGNORECASE,
    )
]
completion_lines = completion_lines[-20:]

# Find IonQuant settings separately from Philosopher filtering.
ionquant_lines = [
    line.strip() for line in prod_text.splitlines()
    if re.search(
        r"ion FDR|peptide FDR|protein FDR|--ionfdr|--peptidefdr|"
        r"--proteinfdr|mbrmincorr|mbrrttol|mbrtoprun|use-labeling",
        line,
        re.IGNORECASE,
    )
]
ionquant_lines = ionquant_lines[-30:]

# ---------------------------------------------------------------------
# Compose human-readable report.
# ---------------------------------------------------------------------
lines += [
    "CAPRION FRAGPIPE PRODUCTION BATCH SUMMARY",
    "=" * 88,
    f"Production run: {run}",
    f"Report generated: {__import__('datetime').datetime.now().astimezone().isoformat(timespec='seconds')}",
    "",
    "1. BATCH STRUCTURE",
    "-" * 88,
    f"Experiment groups expected: {len(experiments)}",
    f"Per-file QC rows found: {len(sample_rows)}",
    f"Expected input files: 98",
    f"Per-file QC rows match 98 expected files: {'YES' if len(sample_rows) == 98 else 'NO'}",
    "",
    "Per-file QC rows by experiment:",
]

for exp in experiments:
    lines.append(f"  {exp}: {sample_by_exp.get(exp, 0)}")

lines += [
    "",
    "2. PER-EXPERIMENT IDENTIFICATION COUNTS",
    "-" * 88,
    "PSM_TSV_Rows/Peptide_TSV_Rows/Protein_TSV_Rows are row counts in the",
    "post-run TSV files. Philosopher_Final_* are extracted from each",
    "experiment's filter.log final post-processing summary. They are kept",
    "separate because their generation stages or definitions may differ.",
    "",
    "Experiment\tPSM TSV\tPeptide TSV\tProtein TSV\tFinal PSMs\t"
    "Final peptides\tFinal proteins\tFinal ions",
]

for row in experiment_rows:
    lines.append(
        "\t".join(str(row.get(k, "")) for k in (
            "Experiment",
            "PSM_TSV_Rows",
            "Peptide_TSV_Rows",
            "Protein_TSV_Rows",
            "Philosopher_Final_PSMs",
            "Philosopher_Final_Peptides",
            "Philosopher_Final_Proteins",
            "Philosopher_Final_Ions",
        ))
    )

lines += [
    "",
    "3. PHILOSOPHER FILTERING EVIDENCE",
    "-" * 88,
    "The per-experiment filter.log files report target-decoy FDR estimates.",
    "The previously inspected logs show approximately 1% FDR at the initial",
    "PSM, peptide, ion and protein filtering stages, followed by sequential",
    "FDR estimation and post-processing. Final counts are included above.",
    "See each experiment's filter.log for thresholds and decoy counts.",
    "",
    "4. PER-FILE QC",
    "-" * 88,
]

if psm_values:
    lines.append(f"Per-file PSM count range: {min(psm_values)}–{max(psm_values)}")
    lines.append(f"Median per-file PSM count: {sorted(psm_values)[len(psm_values)//2]}")
else:
    lines.append("Per-file PSM count statistics unavailable.")

lines += [
    "Per-file PSM probabilities are search-engine confidence scores, not FDRs.",
    "A low PSM count is a review flag, not an automatic exclusion criterion.",
    "The Spectrum_File field in the QC table contains pepXML filenames in the",
    "observed output. Do not assume these are the original mzML filenames;",
    "use the manifest/workflow metadata to confirm the mapping.",
    "",
    "5. QUANTIFICATION AND MATCH-BETWEEN-RUNS CAVEAT",
    "-" * 88,
    "Philosopher identification filtering and IonQuant quantification are",
    "distinct stages. The inspected production log shows IonQuant configured",
    "with ion FDR 0.01, while peptide and protein FDR arguments were 1.0.",
    "Do not interpret those IonQuant settings as evidence that Philosopher's",
    "earlier identification FDR filtering was absent. Conversely, do not",
    "claim that all quantitative or match-between-runs transfers have been",
    "validated at 1% peptide/protein FDR solely from the Philosopher logs.",
    "Review IonQuant's quantification/MBR output and intended analysis policy",
    "before making that claim.",
    "",
    "Relevant IonQuant settings/log lines (selected):",
]

lines.extend(ionquant_lines or ["No matching IonQuant lines found in production log."])

lines += [
    "",
    "6. INPUT VALIDATION AND IDENTIFICATION REPORTS",
    "-" * 88,
    f"Validation report: {validation}",
    f"Identification summary: {ident_summary}",
    f"Sample QC report: {sample_qc}",
    "These existing reports were read, not modified.",
    "",
    "7. PRODUCTION COMPLETION EVIDENCE",
    "-" * 88,
]
lines.extend(completion_lines or [
    "No recognised completion marker found by the script; inspect the full log."
])

lines += [
    "",
    "8. WARNINGS / ITEMS TO REVIEW",
    "-" * 88,
]

if warnings:
    lines.extend(f"WARNING: {item}" for item in warnings)
else:
    lines.append("No missing paths detected by this script.")

lines += [
    "",
    "Review item: TSV peptide row counts differ from Philosopher final peptide",
    "counts in at least some experiments. This report preserves both figures;",
    "the difference should be explained from table provenance before choosing",
    "one count as the canonical batch-wide peptide count.",
    "",
    "Review item: verify the original mzML-to-pepXML mapping before using the",
    "QC Spectrum_File column as a raw-sample identifier.",
    "",
    "OUTPUT FILES",
    "-" * 88,
    f"Batch report: {out_txt}",
    f"Experiment count table: {out_exp}",
    f"Sample QC copy with filename caution: {out_sample}",
    "",
]

with open(out_txt, "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines))

print("\n".join(lines))
PY

echo
echo "Batch report written: $OUT_TXT"
echo "Experiment counts written: $OUT_EXP"
echo "Sample QC copy written: $OUT_SAMPLE"
