#!/usr/bin/env python3
"""Audit scan-level isotope-like spacing for one provisional charge-2 feature.

Diagnostic only:
- does not modify candidate isotope groups;
- does not assign monoisotopic anchors;
- does not use the original XLSX assignments.
"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd

ROOT = Path.cwd()
INPUT = ROOT / "openms_pilot" / "isotope_trace_peaks.tsv"
OUT = ROOT / "reports"
OUT.mkdir(parents=True, exist_ok=True)

FEATURE_ID = "10275997325121874585"
CHARGE = 2
TRACES = [1, 2, 3, 4, 5, 6]
C13 = 1.00335483507
EXPECTED_STEP = C13 / CHARGE
TOLERANCE_MDA = 2.0

if not INPUT.is_file():
    raise FileNotFoundError(f"Missing input: {INPUT}")

df = pd.read_csv(INPUT, sep="\t", dtype={"feature_id": str})
required = {
    "feature_id", "feature_charge", "trace_index", "scan_index",
    "rt_sec", "observed_mz", "peak_intensity"
}
missing = required.difference(df.columns)
if missing:
    raise ValueError(f"Input is missing required columns: {sorted(missing)}")

df["feature_charge"] = pd.to_numeric(df["feature_charge"], errors="coerce")
df["trace_index"] = pd.to_numeric(df["trace_index"], errors="coerce")
df["scan_index"] = pd.to_numeric(df["scan_index"], errors="coerce")
for col in ["rt_sec", "observed_mz", "peak_intensity"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

x = df.loc[
    (df["feature_id"] == FEATURE_ID)
    & (df["feature_charge"] == CHARGE)
    & (df["trace_index"].isin(TRACES))
].copy()

if x.empty:
    raise ValueError(
        f"No rows found for feature {FEATURE_ID}, charge {CHARGE}. "
        "Check the input file and trace-index convention."
    )

# Collapse duplicate rows only if they share the same trace and scan.
# Median m/z avoids silently treating duplicate records as independent evidence.
scan = (
    x.groupby(["trace_index", "scan_index"], as_index=False)
    .agg(
        rt_sec=("rt_sec", "median"),
        observed_mz=("observed_mz", "median"),
        peak_intensity=("peak_intensity", "sum"),
        source_rows=("observed_mz", "size"),
    )
)

scan.to_csv(OUT / "charge2_envelope_scan_observations_v2.tsv",
            sep="\t", index=False)

pair_rows = []
for a, b in itertools.combinations(TRACES, 2):
    aa = scan.loc[scan.trace_index == a,
                  ["scan_index", "rt_sec", "observed_mz", "peak_intensity"]]
    bb = scan.loc[scan.trace_index == b,
                  ["scan_index", "rt_sec", "observed_mz", "peak_intensity"]]
    m = aa.merge(bb, on="scan_index", suffixes=("_a", "_b"))
    if m.empty:
        pair_rows.append({
            "trace_a": a, "trace_b": b, "common_scans": 0,
            "median_delta_mDa": np.nan, "mad_delta_mDa": np.nan,
            "nearest_step": np.nan, "spacing_residual_mDa": np.nan,
            "fraction_within_2mDa_of_step": np.nan,
            "median_rt_delta_sec": np.nan,
            "intensity_correlation": np.nan,
        })
        continue

    # Use the positive mass difference irrespective of trace order.
    delta = (m.observed_mz_b - m.observed_mz_a) * 1000
    delta = delta.where(delta >= 0, -delta)
    med = float(delta.median())
    mad = float((delta - med).abs().median())
    nearest_step = max(1, int(round(med / (EXPECTED_STEP * 1000))))
    residual = abs(med - nearest_step * EXPECTED_STEP * 1000)
    within = float(
        (delta - nearest_step * EXPECTED_STEP * 1000).abs()
        .le(TOLERANCE_MDA).mean()
    )

    corr = np.nan
    if len(m) >= 3:
        corr = m["peak_intensity_a"].corr(m["peak_intensity_b"])

    pair_rows.append({
        "trace_a": a,
        "trace_b": b,
        "common_scans": len(m),
        "median_delta_mDa": med,
        "mad_delta_mDa": mad,
        "nearest_step": nearest_step,
        "spacing_residual_mDa": residual,
        "fraction_within_2mDa_of_step": within,
        "median_rt_delta_sec": float(
            (m.rt_sec_b - m.rt_sec_a).median()
        ),
        "intensity_correlation": corr,
    })

pairs = pd.DataFrame(pair_rows)
pairs.to_csv(OUT / "charge2_envelope_pair_audit_v2.tsv",
             sep="\t", index=False)

print(f"Feature: {FEATURE_ID}; charge: {CHARGE}")
print(f"Rows selected: {len(x):,}")
print(f"Unique trace-scan observations: {len(scan):,}")
print(f"Expected one-isotope spacing: {EXPECTED_STEP:.9f} m/z")
print("\nTrace coverage:")
print(
    scan.groupby("trace_index")
    .agg(
        scans=("scan_index", "nunique"),
        mz_median=("observed_mz", "median"),
        intensity_sum=("peak_intensity", "sum"),
    ).to_string()
)
print("\nPairwise spacing audit:")
print(pairs.round(4).to_string(index=False))
print("\nWrote:")
print(OUT / "charge2_envelope_scan_observations_v2.tsv")
print(OUT / "charge2_envelope_pair_audit_v2.tsv")
print("\nDiagnostic only: no IG assignments or anchors were changed.")
