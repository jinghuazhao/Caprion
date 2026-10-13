#!/usr/bin/env python3
"""Export per-scan trace profiles for the closest candidate IG pairs."""

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path("openms_pilot")
PAIR_FILE = BASE / "candidate_isotope_group_pair_relationships_rt_v2.tsv"
MEMBER_FILE = BASE / "candidate_isotope_groups_v2.tsv"

TRACE_FILES = {
    "first_pilot": BASE / "isotope_trace_peaks.tsv",
    "second_acquisition": BASE / "isotope_trace_peaks_second.tsv",
}

OUT = BASE / "candidate_isotope_group_pair_profiles_v2.tsv"
SELECTED = BASE / "candidate_isotope_group_pair_profiles_selected_v2.tsv"

N_PAIRS = 20


def main():
    pairs = pd.read_csv(PAIR_FILE, sep="\t", dtype={"feature_id": str})
    members = pd.read_csv(
        MEMBER_FILE, sep="\t",
        dtype={"feature_id": str, "trace_index": str},
    )

    required_pairs = {
        "acquisition", "feature_id", "charge",
        "candidate_ig_id_a", "candidate_ig_id_b",
        "nearest_distance_mDa",
    }
    # Prefer the actual RT-audit table's available m/z metric if needed.
    if "nearest_distance_mDa" not in pairs.columns:
        if "best_median_nearest_error_mDa" in pairs.columns:
            pairs["nearest_distance_mDa"] = pairs[
                "best_median_nearest_error_mDa"
            ]
        else:
            pairs["nearest_distance_mDa"] = np.nan

    missing = required_pairs - set(pairs.columns)
    if missing:
        raise ValueError(f"Pair table lacks columns: {sorted(missing)}")

    # Select the closest pairs by their original nearest trace distance,
    # if available; otherwise use the diagnostic m/z error.
    pairs = pairs.sort_values(
        ["nearest_distance_mDa", "acquisition"],
        na_position="last",
    ).head(N_PAIRS).copy()

    pairs.to_csv(SELECTED, sep="\t", index=False)

    # Map each selected candidate group to its trace indices.
    wanted = []
    for _, p in pairs.iterrows():
        for label, ig_col in [
            ("A", "candidate_ig_id_a"),
            ("B", "candidate_ig_id_b"),
        ]:
            g = members[
                (members["acquisition"].astype(str) == str(p["acquisition"]))
                & (members["feature_id"].astype(str) == str(p["feature_id"]))
                & (members["charge"].astype(int) == int(p["charge"]))
                & (members["candidate_ig_id"] == p[ig_col])
            ]
            for _, m in g.iterrows():
                wanted.append({
                    "acquisition": str(p["acquisition"]),
                    "feature_id": str(p["feature_id"]),
                    "charge": int(p["charge"]),
                    "pair_id": (
                        f"{p['acquisition']}:{p['feature_id']}:"
                        f"z{int(p['charge'])}:"
                        f"{p['candidate_ig_id_a']}__{p['candidate_ig_id_b']}"
                    ),
                    "pair_side": label,
                    "candidate_ig_id": p[ig_col],
                    "trace_index": str(m["trace_index"]),
                    "trace_mz_median": m["trace_mz_median"],
                    "trace_n_scans": m["trace_n_scans"],
                })

    wanted_df = pd.DataFrame(wanted)
    if wanted_df.empty:
        print("No trace memberships found for selected pairs.")
        return

    profile_parts = []

    for acquisition, path in TRACE_FILES.items():
        w = wanted_df[wanted_df["acquisition"] == acquisition]
        if w.empty:
            continue

        peaks = pd.read_csv(
            path, sep="\t",
            usecols=[
                "feature_id", "trace_index", "scan_index",
                "rt_sec", "observed_mz", "peak_intensity",
            ],
            dtype={"feature_id": str, "trace_index": str},
        )
        peaks["rt_sec"] = pd.to_numeric(peaks["rt_sec"], errors="coerce")
        peaks["observed_mz"] = pd.to_numeric(
            peaks["observed_mz"], errors="coerce"
        )
        peaks["peak_intensity"] = pd.to_numeric(
            peaks["peak_intensity"], errors="coerce"
        )

        w = w.copy()
        w["feature_id"] = w["feature_id"].astype(str)
        peaks = peaks.merge(
            w,
            on=["feature_id", "trace_index"],
            how="inner",
            validate="many_to_one",
        )
        peaks["acquisition"] = acquisition
        profile_parts.append(peaks)

    if not profile_parts:
        print("No per-scan profiles matched the selected traces.")
        return

    profiles = pd.concat(profile_parts, ignore_index=True)
    profiles = profiles.sort_values(
        ["pair_id", "pair_side", "trace_index", "rt_sec", "scan_index"]
    )
    profiles.to_csv(OUT, sep="\t", index=False)

    print(f"Selected pairs: {SELECTED}")
    print(f"Per-scan profiles: {OUT}")
    print(f"Pairs selected: {pairs.shape[0]}")
    print(f"Trace memberships selected: {len(wanted_df)}")
    print(f"Per-scan rows exported: {len(profiles)}")

    print("\n=== Rows by acquisition and pair side ===")
    print(
        profiles.groupby(["acquisition", "pair_side"])
        .size().to_string()
    )

    print("\n=== Trace profile summary ===")
    summary = profiles.groupby(
        ["acquisition", "pair_id", "pair_side", "trace_index"],
        as_index=False,
    ).agg(
        scans=("scan_index", "nunique"),
        rt_min_sec=("rt_sec", "min"),
        rt_max_sec=("rt_sec", "max"),
        mz_median=("observed_mz", "median"),
        mz_min=("observed_mz", "min"),
        mz_max=("observed_mz", "max"),
        intensity_sum=("peak_intensity", "sum"),
    )
    print(summary.head(30).to_string(index=False))


if __name__ == "__main__":
    main()
