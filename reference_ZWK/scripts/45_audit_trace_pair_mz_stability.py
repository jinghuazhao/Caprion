#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
data = root / "openms_pilot"
src = data / "candidate_isotope_group_selected_scan_trajectories_v2.tsv"
detail_out = data / "candidate_isotope_group_trace_pair_stability_detail_v2.tsv"
summary_out = data / "candidate_isotope_group_trace_pair_stability_summary_v2.tsv"

print(f"Reading {src}", flush=True)
if not src.exists():
    raise FileNotFoundError(src)

df = pd.read_csv(src, sep="\t", low_memory=False)
required = [
    "pair_id", "feature_id", "scan_index", "trace_A", "trace_B",
    "mz_A", "mz_B", "rt_A_sec", "rt_B_sec"
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}; available: {df.columns.tolist()}")

for c in ["scan_index", "trace_A", "trace_B", "mz_A", "mz_B",
          "rt_A_sec", "rt_B_sec"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=[
    "pair_id", "scan_index", "trace_A", "trace_B", "mz_A", "mz_B"
]).copy()
df["scan_index"] = df["scan_index"].astype(int)
df["trace_A"] = df["trace_A"].astype(int)
df["trace_B"] = df["trace_B"].astype(int)

keys = ["pair_id", "trace_A", "trace_B", "scan_index"]
dups = df[df.duplicated(keys, keep=False)]
if not dups.empty:
    dup_path = data / "candidate_isotope_group_trace_pair_duplicate_rows_v2.tsv"
    dups.to_csv(dup_path, sep="\t", index=False)
    raise ValueError(f"Duplicate trace-pair/scan rows: {len(dups)}; inspect {dup_path}")

df["mz_delta"] = df["mz_B"] - df["mz_A"]
df["mz_delta_mDa"] = df["mz_delta"] * 1000.0
df["rt_delta_sec"] = df["rt_B_sec"] - df["rt_A_sec"]

group_keys = ["pair_id", "trace_A", "trace_B"]
df = df.sort_values(group_keys + ["scan_index"]).copy()
df["adjacent_delta_change_mDa"] = (
    df.groupby(group_keys)["mz_delta_mDa"].diff()
)
df.to_csv(detail_out, sep="\t", index=False)

rows = []
spacing_mDa = 1000.0 * 1.00335483507
for (pair_id, trace_a, trace_b), g in df.groupby(group_keys, sort=False):
    delta = g["mz_delta_mDa"].to_numpy(dtype=float)
    adjacent = g["adjacent_delta_change_mDa"].dropna().to_numpy(dtype=float)
    med = float(np.median(delta))
    nearest_n = int(round(med / spacing_mDa))
    rt_delta = g["rt_delta_sec"].dropna().abs()

    rows.append({
        "pair_id": pair_id,
        "feature_id": g["feature_id"].iloc[0],
        "trace_A": trace_a,
        "trace_B": trace_b,
        "n_scans": len(g),
        "scan_index_min": int(g["scan_index"].min()),
        "scan_index_max": int(g["scan_index"].max()),
        "median_delta_mDa": med,
        "delta_MAD_mDa": float(np.median(np.abs(delta - med))),
        "delta_range_mDa": float(np.ptp(delta)),
        "fraction_within_2mDa_of_median": float(np.mean(np.abs(delta - med) <= 2.0)),
        "median_abs_adjacent_change_mDa": float(np.median(np.abs(adjacent))) if len(adjacent) else np.nan,
        "max_abs_adjacent_change_mDa": float(np.max(np.abs(adjacent))) if len(adjacent) else np.nan,
        "median_abs_rt_delta_sec": float(np.median(rt_delta)) if len(rt_delta) else np.nan,
        "nearest_integer_isotope_spacing": nearest_n,
        "residual_to_nearest_integer_spacing_mDa": med - nearest_n * spacing_mDa,
    })

summary = pd.DataFrame(rows)
summary.to_csv(summary_out, sep="\t", index=False)

print(f"Input rows: {len(df)}")
print(f"Trace pairs: {len(summary)}")
print(f"Trace pairs with >=5 scans: {(summary['n_scans'] >= 5).sum()}")
print(f"Detail: {detail_out}")
print(f"Summary: {summary_out}")
print("\n=== Trace-pair-specific stability ===")
print(summary.sort_values(["pair_id", "trace_A", "trace_B"]).to_string(index=False))
