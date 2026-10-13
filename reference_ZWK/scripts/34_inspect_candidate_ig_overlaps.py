#!/usr/bin/env python3
"""Expand potentially overlapping candidate isotope-group pairs for review."""

from pathlib import Path
import pandas as pd

BASE = Path("openms_pilot")
PAIR_FILE = BASE / "candidate_isotope_group_overlap_v2.tsv"
MEMBER_FILE = BASE / "candidate_isotope_groups_v2.tsv"
OUT = BASE / "candidate_isotope_group_overlap_review_v2.tsv"
SUMMARY = BASE / "candidate_isotope_group_overlap_review_summary_v2.tsv"

pairs = pd.read_csv(PAIR_FILE, sep="\t")
members = pd.read_csv(MEMBER_FILE, sep="\t")

if pairs.empty:
    print("No candidate-group overlap pairs to inspect.")
    raise SystemExit(0)

key = ["acquisition", "feature_id", "charge"]
lookup = {}

for k, g in members.groupby(key, sort=False):
    lookup[k] = {
        ig_id: x.sort_values("trace_mz_median")
        for ig_id, x in g.groupby("candidate_ig_id", sort=False)
    }

rows = []

for _, p in pairs.iterrows():
    k = (p["acquisition"], p["feature_id"], int(p["charge"]))
    groups = lookup.get(k, {})
    ga = groups.get(p["candidate_ig_id_a"])
    gb = groups.get(p["candidate_ig_id_b"])

    if ga is None or gb is None:
        continue

    ma = ga["trace_mz_median"].to_numpy(float)
    mb = gb["trace_mz_median"].to_numpy(float)
    spacing = 1.00335483507 / int(p["charge"])

    nearest = min(
        (abs(a - b), a, b) for a in ma for b in mb
    )
    distance, mz_a, mz_b = nearest

    rows.append({
        **p.to_dict(),
        "n_traces_a": len(ma),
        "n_traces_b": len(mb),
        "mz_positions_a": ",".join(f"{x:.6f}" for x in ma),
        "mz_positions_b": ",".join(f"{x:.6f}" for x in mb),
        "nearest_trace_mz_a": mz_a,
        "nearest_trace_mz_b": mz_b,
        "nearest_trace_distance_mz": distance,
        "nearest_distance_as_isotope_steps": distance / spacing,
        "nearest_distance_mDa": distance * 1000,
        "status_a": ga["candidate_status"].iloc[0],
        "status_b": gb["candidate_status"].iloc[0],
    })

out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False)

if not out.empty:
    summary = (
        out.groupby(["acquisition", "status_a", "status_b"])
        .agg(
            pairs=("candidate_ig_id_a", "size"),
            median_nearest_distance_mDa=(
                "nearest_distance_mDa", "median"
            ),
            min_nearest_distance_mDa=(
                "nearest_distance_mDa", "min"
            ),
            max_nearest_distance_mDa=(
                "nearest_distance_mDa", "max"
            ),
            median_traces_a=("n_traces_a", "median"),
            median_traces_b=("n_traces_b", "median"),
        )
        .reset_index()
    )
    summary.to_csv(SUMMARY, sep="\t", index=False)

print(f"Expanded overlap review: {OUT}")
print(f"Overlap review summary: {SUMMARY}")
print(f"Pairs expanded: {len(out)}")

if not out.empty:
    print("\n=== Nearest trace distance, mDa ===")
    print(out["nearest_distance_mDa"].describe().to_string())
    print("\n=== Pair types ===")
    print(
        out.groupby(["acquisition", "status_a", "status_b"])
        .size().to_string()
    )
    print("\n=== Closest 20 pairs ===")
    cols = [
        "acquisition", "feature_id", "charge",
        "candidate_ig_id_a", "candidate_ig_id_b",
        "nearest_trace_distance_mz", "nearest_distance_mDa",
        "nearest_distance_as_isotope_steps",
    ]
    print(
        out.sort_values("nearest_distance_mDa")[cols]
        .head(20).to_string(index=False)
    )
