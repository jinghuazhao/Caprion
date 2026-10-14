#!/usr/bin/env python3
"""Inspect raw trace-peak records at selected charge-2 spacing outliers.

Diagnostic only: does not change IG assignments or assign monoisotopic anchors.
"""
from pathlib import Path
import pandas as pd

ROOT = Path.cwd()
INPUT = ROOT / "openms_pilot" / "isotope_trace_peaks.tsv"
OUT = ROOT / "reports"
OUT.mkdir(parents=True, exist_ok=True)

FEATURE_ID = "10275997325121874585"
CHARGE = 2
TRACE_IDS = [1, 2, 3, 4, 5, 6]
SCAN_IDS = [5954, 5981, 5988]

if not INPUT.is_file():
    raise FileNotFoundError(f"Missing input: {INPUT}")

df = pd.read_csv(INPUT, sep="\t", dtype={"feature_id": str})
required = {
    "feature_id", "feature_charge", "trace_index", "scan_index",
    "rt_sec", "observed_mz", "peak_intensity"
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

for col in [
    "feature_charge", "trace_index", "scan_index",
    "rt_sec", "observed_mz", "peak_intensity"
]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

x = df.loc[
    (df["feature_id"] == FEATURE_ID)
    & (df["feature_charge"] == CHARGE)
    & df["trace_index"].isin(TRACE_IDS)
    & df["scan_index"].isin(SCAN_IDS)
].copy()

if x.empty:
    raise ValueError("No raw records found for the selected feature and scans.")

sort_cols = ["scan_index", "trace_index"]
if "trace_index" in x.columns:
    sort_cols.append("observed_mz")
x = x.sort_values(sort_cols)

detail_path = OUT / "charge2_raw_outlier_peak_records_v2.tsv"
x.to_csv(detail_path, sep="\t", index=False)

print(f"Feature: {FEATURE_ID}; charge: {CHARGE}")
print(f"Selected scans: {SCAN_IDS}")
print(f"Raw records found: {len(x)}")

print("\nRecord counts by scan and trace:")
print(
    x.groupby(["scan_index", "trace_index"])
    .agg(
        records=("observed_mz", "size"),
        distinct_mz=("observed_mz", "nunique"),
        median_mz=("observed_mz", "median"),
        min_mz=("observed_mz", "min"),
        max_mz=("observed_mz", "max"),
        summed_intensity=("peak_intensity", "sum"),
    ).round(6).to_string()
)

print("\nRaw records:")
print(x.to_string(index=False))
print("\nWrote:")
print(detail_path)
print("\nDiagnostic only: no IG assignments or anchors were changed.")
