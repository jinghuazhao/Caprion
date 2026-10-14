#!/usr/bin/env python3
"""Audit adjacent isotope-like trace spacing against theoretical charge-2 spacing."""

from pathlib import Path
import csv
import statistics
from collections import defaultdict

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
INPUT = REPORTS / "full_trace_chromatography_v1.tsv"

FEATURE_ID = "10275997325121874585"
THEORETICAL_SPACING_DA = 1.00335483507 / 2
TOLERANCE_MDA = 3.0

DETAIL_OUT = REPORTS / "adjacent_theoretical_spacing_v1.tsv"
SUMMARY_OUT = REPORTS / "adjacent_theoretical_spacing_summary_v1.tsv"


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def median_abs(values):
    return statistics.median(abs(x) for x in values) if values else None


if not INPUT.exists():
    raise FileNotFoundError(
        f"Missing {INPUT}. Run the script that creates "
        "reports/full_trace_chromatography_v1.tsv first."
    )

rows = read_tsv(INPUT)
by_scan = defaultdict(dict)

for row in rows:
    if row["feature_id"] != FEATURE_ID:
        continue
    scan = int(row["scan_index_zero_based"])
    trace = int(row["trace_index"])
    by_scan[scan][trace] = {
        "rt_sec": float(row["rt_sec"]),
        "observed_mz": float(row["observed_mz"]),
        "peak_intensity": float(row["peak_intensity"]),
    }

detail = []
pairs = [(i, i + 1) for i in range(1, 6)]

for scan in sorted(by_scan):
    traces = by_scan[scan]
    for lower, upper in pairs:
        if lower not in traces or upper not in traces:
            continue

        lo = traces[lower]
        hi = traces[upper]
        spacing = hi["observed_mz"] - lo["observed_mz"]
        residual_mda = (spacing - THEORETICAL_SPACING_DA) * 1000

        detail.append({
            "feature_id": FEATURE_ID,
            "scan_index_zero_based": scan,
            "rt_sec": f'{lo["rt_sec"]:.6f}',
            "lower_trace_index": lower,
            "upper_trace_index": upper,
            "lower_observed_mz": f'{lo["observed_mz"]:.9f}',
            "upper_observed_mz": f'{hi["observed_mz"]:.9f}',
            "observed_spacing_Da": f"{spacing:.9f}",
            "theoretical_spacing_Da": f"{THEORETICAL_SPACING_DA:.12f}",
            "spacing_residual_mDa": f"{residual_mda:.6f}",
            "abs_spacing_residual_mDa": f"{abs(residual_mda):.6f}",
            "within_3mDa": abs(residual_mda) <= TOLERANCE_MDA,
            "lower_intensity": f'{lo["peak_intensity"]:.6f}',
            "upper_intensity": f'{hi["peak_intensity"]:.6f}',
        })

with DETAIL_OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(detail[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(detail)

summary = []
for lower, upper in pairs:
    subset = [
        r for r in detail
        if int(r["lower_trace_index"]) == lower
        and int(r["upper_trace_index"]) == upper
    ]
    residuals = [float(r["spacing_residual_mDa"]) for r in subset]
    spacings = [float(r["observed_spacing_Da"]) for r in subset]
    within = sum(bool(r["within_3mDa"]) for r in subset)

    summary.append({
        "lower_trace_index": lower,
        "upper_trace_index": upper,
        "n_matched_scans": len(subset),
        "median_observed_spacing_Da": (
            f"{statistics.median(spacings):.9f}" if spacings else ""
        ),
        "median_spacing_residual_mDa": (
            f"{statistics.median(residuals):.6f}" if residuals else ""
        ),
        "median_abs_spacing_residual_mDa": (
            f"{median_abs(residuals):.6f}" if residuals else ""
        ),
        "max_abs_spacing_residual_mDa": (
            f"{max(abs(x) for x in residuals):.6f}" if residuals else ""
        ),
        "within_3mDa_count": within,
        "within_3mDa_fraction": (
            f"{within / len(subset):.4f}" if subset else ""
        ),
    })

with SUMMARY_OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

print(f"Feature: {FEATURE_ID}")
print(f"Theoretical charge-2 spacing: {THEORETICAL_SPACING_DA:.12f} Da")
print(f"Per-scan rows: {len(detail)}")
print(f"Detailed report: {DETAIL_OUT}")
print(f"Summary report: {SUMMARY_OUT}")
print("\nPair summary:")
for row in summary:
    print(
        f'{row["lower_trace_index"]}-{row["upper_trace_index"]}: '
        f'n={row["n_matched_scans"]}, '
        f'median spacing={row["median_observed_spacing_Da"]} Da, '
        f'median |residual|={row["median_abs_spacing_residual_mDa"]} mDa, '
        f'within ±3 mDa={row["within_3mDa_count"]}/'
        f'{row["n_matched_scans"]}'
    )

print(
    "\nInterpretation note: this is a spacing diagnostic, not proof that "
    "the traces belong to the same isotope envelope."
)
