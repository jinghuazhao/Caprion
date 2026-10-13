#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

root = Path(__file__).resolve().parents[1]
data = root / "openms_pilot"

trace_files = {
    "first_pilot": data / "isotope_trace_peaks.tsv",
    "second_acquisition": data / "isotope_trace_peaks_second.tsv",
}
groups_file = data / "candidate_isotope_groups_v2.tsv"
passing_file = data / "candidate_isotope_group_charge_aware_spacing_summary_v2.tsv"

out_trace = data / "candidate_isotope_group_passing_feature_trace_inventory_v2.tsv"
out_membership = data / "candidate_isotope_group_passing_feature_membership_audit_v2.tsv"
out_pairs = data / "candidate_isotope_group_passing_feature_pair_membership_v2.tsv"

target_features = {
    ("first_pilot", "10275997325121874585"),
    ("second_acquisition", "12167939761234023158"),
}

def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", low_memory=False)

print("Reading candidate IG memberships...", flush=True)
groups = read_tsv(groups_file)
passing = read_tsv(passing_file)

required_groups = {
    "acquisition", "candidate_ig_id", "feature_id",
    "charge", "trace_index"
}
missing = required_groups - set(groups.columns)
if missing:
    raise ValueError(f"Candidate group table lacks columns: {sorted(missing)}")

# Keep all candidate IGs in the two features, not only the passing trace pairs.
selected_groups = groups.loc[
    [
        (str(a), str(f)) in target_features
        for a, f in zip(groups["acquisition"], groups["feature_id"])
    ]
].copy()

if selected_groups.empty:
    raise ValueError("No candidate groups matched the two target features.")

selected_groups["feature_id"] = selected_groups["feature_id"].astype(str)
selected_groups["trace_index"] = pd.to_numeric(
    selected_groups["trace_index"], errors="coerce"
).astype("Int64")

# Export membership inventory, retaining every original group column.
selected_groups.to_csv(out_membership, sep="\t", index=False)

print("Summarising raw scan-level trace data...", flush=True)
trace_summaries = []

for acquisition, path in trace_files.items():
    if not path.exists():
        print(f"Skipping missing trace file: {path}", flush=True)
        continue

    raw = read_tsv(path)
    required_raw = {
        "feature_id", "feature_charge", "trace_index",
        "scan_index", "rt_sec", "observed_mz", "peak_intensity"
    }
    missing = required_raw - set(raw.columns)
    if missing:
        raise ValueError(f"{path.name} lacks columns: {sorted(missing)}")

    feature_ids = {
        feature for acq, feature in target_features if acq == acquisition
    }
    raw["feature_id"] = raw["feature_id"].astype(str)
    raw = raw.loc[raw["feature_id"].isin(feature_ids)].copy()

    if raw.empty:
        print(f"No raw rows for target feature(s) in {path.name}", flush=True)
        continue

    for col in [
        "feature_charge", "trace_index", "scan_index",
        "rt_sec", "observed_mz", "peak_intensity"
    ]:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")

    for (feature_id, charge, trace_index), g in raw.groupby(
        ["feature_id", "feature_charge", "trace_index"], dropna=False
    ):
        mz = g["observed_mz"].dropna()
        rt = g["rt_sec"].dropna()
        intensity = g["peak_intensity"].dropna()

        trace_summaries.append({
            "acquisition": acquisition,
            "feature_id": feature_id,
            "charge": charge,
            "trace_index": trace_index,
            "n_rows": len(g),
            "n_unique_scans": g["scan_index"].nunique(),
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
            "intensity_median": intensity.median() if len(intensity) else np.nan,
            "intensity_max": intensity.max() if len(intensity) else np.nan,
        })

trace_summary = pd.DataFrame(trace_summaries)
if trace_summary.empty:
    raise ValueError("No trace summaries were generated.")

trace_summary = trace_summary.sort_values(
    ["acquisition", "feature_id", "charge", "trace_index"]
)

# Add the candidate IG IDs and status associated with each trace.
join_keys = ["acquisition", "feature_id", "charge", "trace_index"]
group_lookup_cols = join_keys + [
    c for c in [
        "candidate_ig_id", "candidate_status",
        "monoisotopic_anchor_assigned"
    ] if c in selected_groups.columns
]
group_lookup = selected_groups[group_lookup_cols].drop_duplicates()

trace_summary = trace_summary.merge(
    group_lookup,
    on=join_keys,
    how="left",
    validate="many_to_many"
)
trace_summary.to_csv(out_trace, sep="\t", index=False)

# Annotate whether a trace is in either side of one of the passing pairs.
passing_required = {"feature_id", "trace_A", "trace_B", "charge", "pair_id"}
missing = passing_required - set(passing.columns)
if missing:
    raise ValueError(f"Passing-pair table lacks columns: {sorted(missing)}")

pair_rows = []
for _, p in passing.iterrows():
    feature_id = str(p["feature_id"])
    charge = int(p["charge"])

    for side in ["A", "B"]:
        trace_index = int(p[f"trace_{side}"])
        match = trace_summary.loc[
            (trace_summary["feature_id"].astype(str) == feature_id)
            & (pd.to_numeric(trace_summary["charge"], errors="coerce") == charge)
            & (pd.to_numeric(trace_summary["trace_index"], errors="coerce") == trace_index)
        ].copy()

        if match.empty:
            pair_rows.append({
                "pair_id": p["pair_id"],
                "feature_id": feature_id,
                "charge": charge,
                "trace_side": side,
                "trace_index": trace_index,
                "membership_found": False,
            })
        else:
            match["pair_id"] = p["pair_id"]
            match["trace_side"] = side
            match["membership_found"] = True
            pair_rows.append(match)

pd.concat(
    [x if isinstance(x, pd.DataFrame) else pd.DataFrame([x]) for x in pair_rows],
    ignore_index=True
).to_csv(out_pairs, sep="\t", index=False)

print(f"Membership audit: {out_membership}")
print(f"Trace inventory: {out_trace}")
print(f"Passing-pair membership: {out_pairs}")
print(f"Candidate IG membership rows: {len(selected_groups)}")
print(f"Raw trace summaries: {len(trace_summary)}")

print("\n=== Trace inventory ===")
display_cols = [
    c for c in [
        "acquisition", "feature_id", "charge", "trace_index",
        "candidate_ig_id", "candidate_status",
        "n_unique_scans", "rt_min_sec", "rt_max_sec",
        "mz_median", "mz_range_mDa", "intensity_sum"
    ] if c in trace_summary.columns
]
print(trace_summary[display_cols].to_string(index=False))
