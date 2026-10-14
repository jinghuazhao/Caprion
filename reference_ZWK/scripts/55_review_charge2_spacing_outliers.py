#!/usr/bin/env python3
"""Review scan-level isotope-spacing outliers for one provisional feature.

Diagnostic only: no IG assignments or monoisotopic anchors are changed.
"""
from pathlib import Path
import pandas as pd

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
INPUT = REPORTS / "charge2_spacing_residual_outliers_v2.tsv"
SUMMARY = REPORTS / "charge2_spacing_residual_summary_v2.tsv"

if not INPUT.is_file():
    raise FileNotFoundError(f"Missing Script 54 outlier report: {INPUT}")

outliers = pd.read_csv(INPUT, sep="\t")
summary = (
    pd.read_csv(SUMMARY, sep="\t")
    if SUMMARY.is_file() else pd.DataFrame()
)

if outliers.empty:
    print("No outlier scans were recorded.")
    raise SystemExit(0)

# Retain only columns available in this report version.
preferred = [
    "trace_a", "trace_b", "scan_index",
    "rt_sec_a", "rt_sec_b",
    "observed_mz_a", "observed_mz_b",
    "signed_delta_mDa", "spacing_residual_mDa",
    "abs_residual_mDa", "peak_intensity_a", "peak_intensity_b",
]
columns = [c for c in preferred if c in outliers.columns]
outliers = outliers[columns].copy()

for col in ["scan_index", "rt_sec_a", "spacing_residual_mDa",
            "abs_residual_mDa"]:
    if col in outliers.columns:
        outliers[col] = pd.to_numeric(outliers[col], errors="coerce")

sort_cols = [
    c for c in ["trace_a", "trace_b", "rt_sec_a", "scan_index"]
    if c in outliers.columns
]
outliers = outliers.sort_values(sort_cols)

detail_path = REPORTS / "charge2_spacing_outlier_review_v2.tsv"
outliers.to_csv(detail_path, sep="\t", index=False)

print(f"Outlier records: {len(outliers):,}")
print("\nOutlier counts by adjacent trace pair:")
print(
    outliers.groupby(["trace_a", "trace_b"])
    .agg(
        outlier_scans=("scan_index", "nunique"),
        median_abs_residual_mDa=("abs_residual_mDa", "median"),
        max_abs_residual_mDa=("abs_residual_mDa", "max"),
    ).round(4).to_string()
)

print("\nLargest absolute residuals:")
print(
    outliers.sort_values("abs_residual_mDa", ascending=False)
    .head(15).round(5).to_string(index=False)
)

print("\nOutliers ordered by retention time:")
print(outliers.round(5).to_string(index=False))

if not summary.empty:
    print("\nScript 54 pair summary:")
    print(summary.round(4).to_string(index=False))

print("\nWrote:")
print(detail_path)
print("\nDiagnostic only: no IG assignments or anchors were changed.")
