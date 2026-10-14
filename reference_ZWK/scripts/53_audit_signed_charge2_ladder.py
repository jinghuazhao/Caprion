#!/usr/bin/env python3
"""Audit signed scan-level spacing for adjacent traces in a provisional ladder.

Diagnostic only: no candidate IG assignments or monoisotopic anchors are changed.
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
TRACES = [1, 2, 3, 4, 5]
EXPECTED = 1.00335483507 / CHARGE
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
    & (df.trace_index.isin(TRACES))
].copy()

if x.empty:
    raise ValueError(f"No trace observations found for feature {FEATURE_ID}")

# One observation per trace and scan; preserve duplicate-row counts.
scan = (
    x.groupby(["trace_index", "scan_index"], as_index=False)
    .agg(
        rt_sec=("rt_sec", "median"),
        observed_mz=("observed_mz", "median"),
        peak_intensity=("peak_intensity", "sum"),
        source_rows=("observed_mz", "size"),
    )
)

rows = []
for a, b in zip(TRACES[:-1], TRACES[1:]):
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
        rows.append({
            "trace_a": a, "trace_b": b, "common_scans": 0,
            "positive_order_fraction": np.nan,
            "median_signed_delta_mDa": np.nan,
            "mad_signed_delta_mDa": np.nan,
            "median_spacing_residual_mDa": np.nan,
            "fraction_within_2mDa": np.nan,
            "median_rt_delta_sec": np.nan,
            "intensity_correlation": np.nan,
        })
        continue

    delta = (m.observed_mz_b - m.observed_mz_a) * 1000
    expected_mda = EXPECTED * 1000
    residual = delta - expected_mda
    corr = (
        m.peak_intensity_a.corr(m.peak_intensity_b)
        if len(m) >= 3 else np.nan
    )

    rows.append({
        "trace_a": a,
        "trace_b": b,
        "common_scans": len(m),
        "positive_order_fraction": float((delta > 0).mean()),
        "median_signed_delta_mDa": float(delta.median()),
        "mad_signed_delta_mDa": float((delta - delta.median()).abs().median()),
        "median_spacing_residual_mDa": float(residual.median()),
        "fraction_within_2mDa": float(residual.abs().le(TOLERANCE_MDA).mean()),
        "median_rt_delta_sec": float(
            (m.rt_sec_b - m.rt_sec_a).median()
        ),
        "intensity_correlation": corr,
    })

    m = m.assign(
        trace_a=a,
        trace_b=b,
        signed_delta_mDa=delta,
        expected_delta_mDa=expected_mda,
        spacing_residual_mDa=residual,
        positive_order=(delta > 0),
    )
    # Retain all scan-level observations for reviewing outliers.
    path = OUT / f"charge2_adjacent_trace_{a}_{b}_scan_audit_v2.tsv"
    m.to_csv(path, sep="\t", index=False)

summary = pd.DataFrame(rows)
summary_path = OUT / "charge2_signed_adjacent_ladder_summary_v2.tsv"
scan_path = OUT / "charge2_signed_ladder_scan_observations_v2.tsv"

scan.to_csv(scan_path, sep="\t", index=False)
summary.to_csv(summary_path, sep="\t", index=False)

print(f"Feature: {FEATURE_ID}; charge: {CHARGE}")
print(f"Expected adjacent isotope spacing: {EXPECTED:.9f} m/z")
print("\nAdjacent-trace audit:")
print(summary.round(4).to_string(index=False))
print("\nWrote:")
print(summary_path)
print(scan_path)
print("\nIndividual scan-level pair audits were written for each adjacent pair.")
print("Diagnostic only: no IG assignments or anchors were changed.")
