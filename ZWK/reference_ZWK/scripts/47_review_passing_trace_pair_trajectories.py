#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

root = Path(__file__).resolve().parents[1]
data = root / "openms_pilot"

detail_file = data / "candidate_isotope_group_charge_aware_spacing_detail_v2.tsv"
summary_file = data / "candidate_isotope_group_charge_aware_spacing_summary_v2.tsv"
groups_file = data / "candidate_isotope_groups_v2.tsv"

out_detail = data / "candidate_isotope_group_passing_pair_trajectory_review_v2.tsv"
out_summary = data / "candidate_isotope_group_passing_pair_trajectory_summary_v2.tsv"
out_members = data / "candidate_isotope_group_passing_pair_membership_v2.tsv"

print("Reading charge-aware spacing audit...", flush=True)
summary = pd.read_csv(summary_file, sep="\t", low_memory=False)
detail = pd.read_csv(detail_file, sep="\t", low_memory=False)
groups = pd.read_csv(groups_file, sep="\t", low_memory=False)

required_summary = {
    "pair_id", "feature_id", "trace_A", "trace_B", "charge",
    "passes_both_screening_rules"
}
missing = required_summary - set(summary.columns)
if missing:
    raise ValueError(f"Summary is missing columns: {sorted(missing)}")

passing = summary.loc[
    summary["passes_both_screening_rules"].astype(str).str.lower().isin(
        ["true", "1", "yes"]
    )
].copy()

if passing.empty:
    raise ValueError("No trace pairs passed both screening rules.")

print(f"Passing trace pairs: {len(passing)}", flush=True)

# Select the exact scan-level rows belonging to the passing trace pairs.
keys = ["pair_id", "trace_A", "trace_B"]
passing_keys = passing[keys].drop_duplicates()
review = detail.merge(passing_keys, on=keys, how="inner", validate="many_to_one")

if review.empty:
    raise ValueError("No scan-level rows matched the passing trace pairs.")

# Add the charge-aware expected spacing and isotope-index difference.
review["expected_spacing_mDa"] = (
    1000.0 * 1.00335483507 / review["charge"]
)
review["delta_mDa"] = review["mz_delta"] * 1000.0
review["spacing_step"] = np.rint(
    review["delta_mDa"] / review["expected_spacing_mDa"]
).astype(int)
review["spacing_residual_mDa"] = (
    review["delta_mDa"]
    - review["spacing_step"] * review["expected_spacing_mDa"]
)

review = review.sort_values(keys + ["scan_index"])
review["adjacent_delta_change_mDa"] = (
    review.groupby(keys)["delta_mDa"].diff()
)

review.to_csv(out_detail, sep="\t", index=False)

# Summarise each passing trace pair without pooling different trace combinations.
rows = []
for (pair_id, trace_a, trace_b), g in review.groupby(keys, sort=False):
    g = g.sort_values("scan_index")
    delta = g["delta_mDa"].to_numpy(dtype=float)
    residual = g["spacing_residual_mDa"].to_numpy(dtype=float)
    adj = g["adjacent_delta_change_mDa"].dropna().to_numpy(dtype=float)

    rows.append({
        "pair_id": pair_id,
        "feature_id": g["feature_id"].iloc[0],
        "charge": int(g["charge"].iloc[0]),
        "trace_A": int(trace_a),
        "trace_B": int(trace_b),
        "n_scans": len(g),
        "median_delta_mDa": float(np.median(delta)),
        "delta_MAD_mDa": float(np.median(np.abs(delta - np.median(delta)))),
        "delta_range_mDa": float(np.ptp(delta)),
        "fraction_within_2mDa_of_median": float(
            np.mean(np.abs(delta - np.median(delta)) <= 2)
        ),
        "median_abs_spacing_residual_mDa": float(np.median(np.abs(residual))),
        "fraction_spacing_residual_le_2mDa": float(np.mean(np.abs(residual) <= 2)),
        "median_abs_adjacent_change_mDa": (
            float(np.median(np.abs(adj))) if len(adj) else np.nan
        ),
        "median_abs_rt_delta_sec": float(
            np.median(np.abs(g["rt_delta_sec"].dropna()))
        ),
    })

pair_summary = pd.DataFrame(rows)
pair_summary.to_csv(out_summary, sep="\t", index=False)

# Map the reviewed trace indices back to candidate IG IDs, where possible.
membership_columns = [
    "acquisition", "candidate_ig_id", "feature_id",
    "charge", "trace_index", "candidate_status",
    "monoisotopic_anchor_assigned"
]
available = [c for c in membership_columns if c in groups.columns]
if {"feature_id", "charge", "trace_index"}.issubset(available):
    member_rows = []
    for _, r in passing.iterrows():
        feature = str(r["feature_id"])
        charge = int(r["charge"])
        for side in ["A", "B"]:
            trace_col = f"trace_{side}"
            trace_index = int(r[trace_col])
            match = groups.loc[
                (groups["feature_id"].astype(str) == feature)
                & (pd.to_numeric(groups["charge"], errors="coerce") == charge)
                & (pd.to_numeric(groups["trace_index"], errors="coerce") == trace_index)
            ].copy()
            if not match.empty:
                match["review_pair_id"] = r["pair_id"]
                match["review_trace_side"] = side
                member_rows.append(match[available + ["review_pair_id", "review_trace_side"]])

    if member_rows:
        pd.concat(member_rows, ignore_index=True).drop_duplicates().to_csv(
            out_members, sep="\t", index=False
        )
        print(f"Membership report: {out_members}")

print(f"Scan-level review: {out_detail}")
print(f"Pair summary: {out_summary}")
print(f"Scan rows reviewed: {len(review)}")
print("\n=== Passing-pair trajectory review ===")
print(pair_summary.to_string(index=False))
