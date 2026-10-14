#!/usr/bin/env python3
"""Audit scan-level spacing residuals for a provisional charge-2 ladder.

Diagnostic only: does not alter IG assignments or assign monoisotopic anchors.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.cwd()
INPUT = ROOT / "openms_pilot" / "isotope_trace_peaks.tsv"
OUT = ROOT / "reports"
OUT.mkdir(parents=True, exist_ok=True)

FEATURE_ID = "10275997325121874585"
CHARGE = 2
ADJACENT_PAIRS = [(1, 2), (2, 3), (3, 4), (4, 5)]
EXPECTED_STEP = 1.00335483507 / CHARGE
TOLERANCE_MDA = 2.0

df = pd.read_csv(INPUT, sep="\t", dtype={"feature_id": str})
required = {
    "feature_id", "feature_charge", "trace_index", "scan_index",
    "rt_sec", "observed_mz", "peak_intensity"
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

for col in ["feature_charge", "trace_index", "scan_index",
            "rt_sec", "observed_mz", "peak_intensity"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

x = df.loc[
    (df.feature_id == FEATURE_ID)
    & (df.feature_charge == CHARGE)
    & df.trace_index.isin(sorted({t for p in ADJACENT_PAIRS for t in p}))
].copy()

if x.empty:
    raise ValueError(f"No observations found for feature {FEATURE_ID}")

# Collapse duplicate records within trace/scan without treating duplicates
# as independent observations.
scan = (
    x.groupby(["trace_index", "scan_index"], as_index=False)
    .agg(
        rt_sec=("rt_sec", "median"),
        observed_mz=("observed_mz", "median"),
        peak_intensity=("peak_intensity", "sum"),
        source_rows=("observed_mz", "size"),
    )
)

details = []
summaries = []
expected_mda = EXPECTED_STEP * 1000

for a, b in ADJACENT_PAIRS:
    aa = scan.loc[
        scan.trace_index == a,
        ["scan_index", "rt_sec", "observed_mz", "peak_intensity"]
    ]
    bb = scan.loc[
        scan.trace_index == b,
        ["scan_index", "rt_sec", "observed_mz", "peak_intensity"]
    ]
    m = aa.merge(bb, on="scan_index", suffixes=("_a", "_b"))

    if m.empty:
        summaries.append({
            "trace_a": a, "trace_b": b, "common_scans": 0
        })
        continue

    m["trace_a"] = a
    m["trace_b"] = b
    m["signed_delta_mDa"] = (
        m["observed_mz_b"] - m["observed_mz_a"]
    ) * 1000
    m["spacing_residual_mDa"] = (
        m["signed_delta_mDa"] - expected_mda
    )
    m["abs_residual_mDa"] = m["spacing_residual_mDa"].abs()
    m["within_2mDa"] = m["abs_residual_mDa"] <= TOLERANCE_MDA
    m["positive_order"] = m["signed_delta_mDa"] > 0

    details.append(m)

    residual = m["spacing_residual_mDa"]
    abs_residual = residual.abs()

    # Retention-time trend: residual in mDa per second.
    drift = np.nan
    if m["rt_sec_a"].nunique() >= 3:
        drift = float(
            np.polyfit(m["rt_sec_a"], residual, 1)[0]
        )

    summaries.append({
        "trace_a": a,
        "trace_b": b,
        "common_scans": len(m),
        "positive_order_fraction": float(m["positive_order"].mean()),
        "median_residual_mDa": float(residual.median()),
        "residual_MAD_mDa": float(
            (residual - residual.median()).abs().median()
        ),
        "p90_abs_residual_mDa": float(abs_residual.quantile(0.90)),
        "max_abs_residual_mDa": float(abs_residual.max()),
        "fraction_within_2mDa": float(m["within_2mDa"].mean()),
        "residual_drift_mDa_per_sec": drift,
    })

detail_df = (
    pd.concat(details, ignore_index=True)
    if details else pd.DataFrame()
)
summary_df = pd.DataFrame(summaries)

detail_path = OUT / "charge2_spacing_residual_scan_detail_v2.tsv"
summary_path = OUT / "charge2_spacing_residual_summary_v2.tsv"
outlier_path = OUT / "charge2_spacing_residual_outliers_v2.tsv"

detail_df.to_csv(detail_path, sep="\t", index=False)
summary_df.to_csv(summary_path, sep="\t", index=False)

if not detail_df.empty:
    outliers = detail_df.loc[~detail_df["within_2mDa"]].copy()
else:
    outliers = pd.DataFrame()
outliers.to_csv(outlier_path, sep="\t", index=False)

print(f"Feature: {FEATURE_ID}; charge: {CHARGE}")
print(f"Expected adjacent spacing: {EXPECTED_STEP:.9f} m/z")
print("\nResidual summary:")
print(summary_df.round(4).to_string(index=False))

print("\nOutlier counts by pair:")
if outliers.empty:
    print("No scan-level residuals exceeded the ±2 mDa screen.")
else:
    print(
        outliers.groupby(["trace_a", "trace_b"])
        .size().rename("outlier_scans").to_string()
    )

print("\nWrote:")
print(summary_path)
print(detail_path)
print(outlier_path)
print("\nDiagnostic only: no IG assignments or anchors were changed.")
