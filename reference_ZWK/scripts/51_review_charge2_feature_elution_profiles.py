#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
data = root / "openms_pilot"

source = data / "candidate_charge2_feature_scan_trajectories_v2.tsv"
out_aligned = data / "candidate_charge2_feature_aligned_intensity_v2.tsv"
out_corr = data / "candidate_charge2_feature_elution_correlation_v2.tsv"
out_summary = data / "candidate_charge2_feature_elution_profile_summary_v2.tsv"

feature_id = "10275997325121874585"
trace_ids = [1, 2, 3, 4, 5, 6]

df = pd.read_csv(source, sep="\t", low_memory=False)

required = {
    "feature_id", "charge", "trace_index", "scan_index",
    "rt_sec", "peak_intensity", "observed_mz"
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Input lacks columns: {sorted(missing)}")

df["feature_id"] = df["feature_id"].astype(str)
for col in ["charge", "trace_index", "scan_index", "rt_sec",
            "peak_intensity", "observed_mz"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.loc[
    (df["feature_id"] == feature_id)
    & (df["charge"] == 2)
    & (df["trace_index"].isin(trace_ids))
].copy()

if df.empty:
    raise RuntimeError("No rows found for the selected feature and traces.")

# Use one row per trace and scan.
df = (
    df.groupby(["scan_index", "trace_index"], as_index=False)
    .agg(
        rt_sec=("rt_sec", "median"),
        intensity=("peak_intensity", "sum"),
        observed_mz=("observed_mz", "median"),
    )
)

# Wide scan-aligned tables. Missing values mean that a trace was not observed
# at that scan; they are not silently replaced with zero for correlations.
intensity = df.pivot(index="scan_index", columns="trace_index", values="intensity")
rt = df.pivot(index="scan_index", columns="trace_index", values="rt_sec")
mz = df.pivot(index="scan_index", columns="trace_index", values="observed_mz")

for t in trace_ids:
    if t not in intensity.columns:
        intensity[t] = np.nan
    if t not in rt.columns:
        rt[t] = np.nan
    if t not in mz.columns:
        mz[t] = np.nan

intensity = intensity.reindex(columns=trace_ids).sort_index()
rt = rt.reindex(index=intensity.index, columns=trace_ids)
mz = mz.reindex(index=intensity.index, columns=trace_ids)

# Normalise each trace by its own maximum. This compares profile shape, not
# absolute abundance.
normalised = intensity.copy()
for t in trace_ids:
    vmax = normalised[t].max(skipna=True)
    if pd.notna(vmax) and vmax > 0:
        normalised[t] = normalised[t] / vmax

aligned = pd.DataFrame({
    "scan_index": intensity.index,
    "rt_sec_median": rt.median(axis=1, skipna=True).to_numpy(),
})

for t in trace_ids:
    aligned[f"trace_{t}_intensity"] = intensity[t].to_numpy()
    aligned[f"trace_{t}_normalised_intensity"] = normalised[t].to_numpy()
    aligned[f"trace_{t}_mz"] = mz[t].to_numpy()
    aligned[f"trace_{t}_observed"] = intensity[t].notna().to_numpy()

aligned.to_csv(out_aligned, sep="\t", index=False)

# Pairwise correlations over scans where both traces are present.
corr_rows = []
for i, a in enumerate(trace_ids):
    for b in trace_ids[i + 1:]:
        both = normalised[[a, b]].dropna()
        corr = np.nan
        if (
            len(both) >= 3
            and both[a].nunique() > 1
            and both[b].nunique() > 1
        ):
            corr = both[a].corr(both[b])

        corr_rows.append({
            "trace_A": a,
            "trace_B": b,
            "n_common_scans": len(both),
            "normalised_intensity_correlation": corr,
        })

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(out_corr, sep="\t", index=False)

# Summary includes where each trace is observed and its intensity profile.
summary_rows = []
for t in trace_ids:
    g = df.loc[df["trace_index"] == t]
    observed = g.dropna(subset=["scan_index", "rt_sec", "intensity"])
    vals = observed["intensity"]
    peak = vals.max() if len(vals) else np.nan
    half = peak / 2 if pd.notna(peak) else np.nan
    above_half = observed.loc[observed["intensity"] >= half] if pd.notna(half) else observed.iloc[0:0]

    summary_rows.append({
        "feature_id": feature_id,
        "charge": 2,
        "trace_index": t,
        "n_observed_scans": observed["scan_index"].nunique(),
        "scan_min": observed["scan_index"].min() if len(observed) else np.nan,
        "scan_max": observed["scan_index"].max() if len(observed) else np.nan,
        "rt_min_sec": observed["rt_sec"].min() if len(observed) else np.nan,
        "rt_median_sec": observed["rt_sec"].median() if len(observed) else np.nan,
        "rt_max_sec": observed["rt_sec"].max() if len(observed) else np.nan,
        "intensity_sum": vals.sum() if len(vals) else np.nan,
        "intensity_max": peak,
        "scans_at_or_above_half_max": len(above_half),
        "rt_half_max_min_sec": above_half["rt_sec"].min() if len(above_half) else np.nan,
        "rt_half_max_max_sec": above_half["rt_sec"].max() if len(above_half) else np.nan,
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(out_summary, sep="\t", index=False)

print(f"Aligned intensity table: {out_aligned}")
print(f"Pairwise correlation table: {out_corr}")
print(f"Elution profile summary: {out_summary}")

print("\n=== Elution profile summary ===")
print(summary.to_string(index=False))

print("\n=== Pairwise normalised-intensity correlations ===")
print(corr_df.to_string(index=False))
