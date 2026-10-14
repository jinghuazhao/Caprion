#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
data = root / "openms_pilot"

source = data / "isotope_trace_peaks.tsv"
groups_file = data / "candidate_isotope_groups_v2.tsv"

out_long = data / "candidate_charge2_feature_scan_trajectories_v2.tsv"
out_pairs = data / "candidate_charge2_feature_scan_pair_deltas_v2.tsv"
out_summary = data / "candidate_charge2_feature_trace_summary_v2.tsv"

feature_id = "10275997325121874585"
charge = 2
C13 = 1.00335483507
EXPECTED_SPACING = C13 / charge

raw = pd.read_csv(source, sep="\t", low_memory=False)

required = {
    "feature_id", "feature_charge", "trace_index",
    "scan_index", "rt_sec", "observed_mz", "peak_intensity"
}
missing = required - set(raw.columns)
if missing:
    raise ValueError(f"Input lacks columns: {sorted(missing)}")

raw["feature_id"] = raw["feature_id"].astype(str)
for col in [
    "feature_charge", "trace_index", "scan_index",
    "rt_sec", "observed_mz", "peak_intensity"
]:
    raw[col] = pd.to_numeric(raw[col], errors="coerce")

raw = raw.loc[
    (raw["feature_id"] == feature_id)
    & (raw["feature_charge"] == charge)
].copy()

if raw.empty:
    raise RuntimeError(f"No rows found for feature {feature_id}, charge {charge}")

# Collapse multiple raw rows within a trace/scan, if present.
long = (
    raw.dropna(subset=["trace_index", "scan_index"])
    .groupby(["trace_index", "scan_index"], as_index=False)
    .agg(
        rt_sec=("rt_sec", "median"),
        observed_mz=("observed_mz", "median"),
        peak_intensity=("peak_intensity", "sum"),
        n_raw_rows=("scan_index", "size"),
    )
)
long["feature_id"] = feature_id
long["charge"] = charge
long["expected_isotope_spacing_Da"] = EXPECTED_SPACING

# Attach existing candidate IG membership for context only.
groups = pd.read_csv(groups_file, sep="\t", low_memory=False)
groups["feature_id"] = groups["feature_id"].astype(str)
groups["trace_index"] = pd.to_numeric(groups["trace_index"], errors="coerce")
groups["charge"] = pd.to_numeric(groups["charge"], errors="coerce")

membership = groups.loc[
    (groups["feature_id"] == feature_id)
    & (groups["charge"] == charge),
    ["trace_index", "candidate_ig_id", "candidate_status"]
].drop_duplicates()

long = long.merge(membership, on="trace_index", how="left")
long = long.sort_values(["trace_index", "scan_index"])
long.to_csv(out_long, sep="\t", index=False)

summary_rows = []
for trace_index, g in long.groupby("trace_index"):
    mz = g["observed_mz"].dropna()
    rt = g["rt_sec"].dropna()
    intensity = g["peak_intensity"].dropna()
    summary_rows.append({
        "feature_id": feature_id,
        "charge": charge,
        "trace_index": int(trace_index),
        "candidate_ig_id": ";".join(
            sorted(g["candidate_ig_id"].dropna().astype(str).unique())
        ),
        "n_scans": g["scan_index"].nunique(),
        "rt_min_sec": rt.min() if len(rt) else np.nan,
        "rt_median_sec": rt.median() if len(rt) else np.nan,
        "rt_max_sec": rt.max() if len(rt) else np.nan,
        "mz_min": mz.min() if len(mz) else np.nan,
        "mz_median": mz.median() if len(mz) else np.nan,
        "mz_max": mz.max() if len(mz) else np.nan,
        "mz_range_mDa": (
            (mz.max() - mz.min()) * 1000 if len(mz) else np.nan
        ),
        "intensity_sum": intensity.sum() if len(intensity) else np.nan,
        "intensity_max": intensity.max() if len(intensity) else np.nan,
    })

summary = pd.DataFrame(summary_rows).sort_values("trace_index")
summary.to_csv(out_summary, sep="\t", index=False)

# Compare all trace pairs scan by scan.
pair_rows = []
traces = {
    int(t): g.set_index("scan_index").sort_index()
    for t, g in long.groupby("trace_index")
}

for trace_a, trace_b in combinations(sorted(traces), 2):
    a = traces[trace_a]
    b = traces[trace_b]
    common = a.index.intersection(b.index)

    if len(common) == 0:
        continue

    aa = a.loc[common]
    bb = b.loc[common]

    delta = (
        bb["observed_mz"].to_numpy(dtype=float)
        - aa["observed_mz"].to_numpy(dtype=float)
    )
    valid = np.isfinite(delta)
    delta = delta[valid]

    if len(delta) == 0:
        continue

    median_delta = float(np.median(delta))
    steps = np.rint(delta / EXPECTED_SPACING)
    residual = delta - steps * EXPECTED_SPACING

    ia = aa["peak_intensity"].to_numpy(dtype=float)
    ib = bb["peak_intensity"].to_numpy(dtype=float)
    good = np.isfinite(ia) & np.isfinite(ib)
    corr = np.nan
    if good.sum() >= 3 and np.std(ia[good]) > 0 and np.std(ib[good]) > 0:
        corr = float(np.corrcoef(ia[good], ib[good])[0, 1])

    rt_delta = (
        bb["rt_sec"].to_numpy(dtype=float)
        - aa["rt_sec"].to_numpy(dtype=float)
    )
    rt_delta = rt_delta[np.isfinite(rt_delta)]

    pair_rows.append({
        "feature_id": feature_id,
        "charge": charge,
        "trace_A": trace_a,
        "trace_B": trace_b,
        "n_common_scans": len(common),
        "median_mz_delta_mDa": median_delta * 1000,
        "mz_delta_MAD_mDa": (
            np.median(np.abs(delta - median_delta)) * 1000
        ),
        "mz_delta_range_mDa": (np.max(delta) - np.min(delta)) * 1000,
        "nearest_isotope_step": int(np.rint(median_delta / EXPECTED_SPACING)),
        "median_abs_spacing_residual_mDa": (
            np.median(np.abs(residual)) * 1000
        ),
        "fraction_spacing_residual_le_2mDa": (
            np.mean(np.abs(residual) <= 0.002)
        ),
        "intensity_correlation": corr,
        "median_rt_delta_sec": (
            float(np.median(rt_delta)) if len(rt_delta) else np.nan
        ),
    })

pairs = pd.DataFrame(pair_rows)
if not pairs.empty:
    pairs = pairs.sort_values(
        ["n_common_scans", "trace_A", "trace_B"],
        ascending=[False, True, True]
    )
pairs.to_csv(out_pairs, sep="\t", index=False)

print(f"Scan trajectories: {out_long}")
print(f"Trace summary: {out_summary}")
print(f"Scan-level pair deltas: {out_pairs}")

print("\n=== Trace summary ===")
print(summary.to_string(index=False))

print("\n=== Trace-pair scan audit ===")
if pairs.empty:
    print("No trace pairs had common scans.")
else:
    cols = [
        "trace_A", "trace_B", "n_common_scans",
        "median_mz_delta_mDa", "mz_delta_MAD_mDa",
        "mz_delta_range_mDa", "nearest_isotope_step",
        "median_abs_spacing_residual_mDa",
        "fraction_spacing_residual_le_2mDa",
        "intensity_correlation",
    ]
    print(pairs[cols].to_string(index=False))
