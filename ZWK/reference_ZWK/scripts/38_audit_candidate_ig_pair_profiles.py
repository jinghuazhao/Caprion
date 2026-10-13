#!/usr/bin/env python3
"""Audit m/z ladder structure and trace-level co-elution in exported IG pairs."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
PROFILE_FILE = BASE / "candidate_isotope_group_pair_profiles_v2.tsv"

OUT_PAIRS = BASE / "candidate_isotope_group_pair_profile_audit_v2.tsv"
OUT_TRACES = BASE / "candidate_isotope_group_pair_trace_audit_v2.tsv"

ISOTOPE_MASS = 1.00335483507
MZ_TOLERANCE_MDA = 2.0
RT_MATCH_TOLERANCE_SEC = 2.0
MIN_COMMON_SCANS = 3


def numeric(df, columns):
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def main():
    if not PROFILE_FILE.exists():
        raise FileNotFoundError(
            f"Missing {PROFILE_FILE}; run script 37 first."
        )

    df = pd.read_csv(
        PROFILE_FILE,
        sep="\t",
        dtype={
            "feature_id": str,
            "trace_index": str,
            "pair_id": str,
            "candidate_ig_id": str,
        },
    )
    required = {
        "acquisition", "feature_id", "charge", "pair_id",
        "pair_side", "candidate_ig_id", "trace_index",
        "scan_index", "rt_sec", "observed_mz", "peak_intensity",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Profile file lacks columns: {sorted(missing)}")

    df = numeric(
        df,
        ["charge", "scan_index", "rt_sec", "observed_mz", "peak_intensity"],
    )
    df = df.dropna(
        subset=["charge", "scan_index", "rt_sec", "observed_mz"]
    ).copy()
    df["charge"] = df["charge"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)

    # One summary row per trace in each selected pair.
    trace_rows = []
    for keys, g in df.groupby(
        [
            "acquisition", "pair_id", "pair_side",
            "candidate_ig_id", "trace_index",
        ],
        dropna=False,
    ):
        acquisition, pair_id, side, ig_id, trace_index = keys
        g = g.sort_values(["rt_sec", "scan_index"])
        trace_rows.append({
            "acquisition": acquisition,
            "feature_id": g["feature_id"].iloc[0],
            "pair_id": pair_id,
            "pair_side": side,
            "candidate_ig_id": ig_id,
            "trace_index": trace_index,
            "charge": int(g["charge"].iloc[0]),
            "scans": g["scan_index"].nunique(),
            "rt_min_sec": g["rt_sec"].min(),
            "rt_max_sec": g["rt_sec"].max(),
            "mz_median": g["observed_mz"].median(),
            "mz_min": g["observed_mz"].min(),
            "mz_max": g["observed_mz"].max(),
            "intensity_sum": g["peak_intensity"].sum(),
        })

    traces = pd.DataFrame(trace_rows)
    if traces.empty:
        print("No trace profiles to audit.")
        return

    traces["rt_mid_sec"] = (
        traces["rt_min_sec"] + traces["rt_max_sec"]
    ) / 2
    traces.to_csv(OUT_TRACES, sep="\t", index=False)

    # Pair-level summaries.
    pair_rows = []
    for pair_id, all_g in traces.groupby("pair_id"):
        a = all_g[all_g["pair_side"] == "A"].copy()
        b = all_g[all_g["pair_side"] == "B"].copy()

        if a.empty or b.empty:
            pair_rows.append({
                "pair_id": pair_id,
                "audit_status": "missing_side",
            })
            continue

        charge_values = all_g["charge"].dropna().unique()
        if len(charge_values) != 1:
            charge = np.nan
            expected = np.nan
        else:
            charge = int(charge_values[0])
            expected = ISOTOPE_MASS / charge

        # Adjacent trace spacing within each side.
        def side_spacing(side_df):
            mzs = np.sort(side_df["mz_median"].dropna().to_numpy())
            if len(mzs) < 2 or not np.isfinite(expected):
                return {
                    "n_adjacent_gaps": 0,
                    "adjacent_gap_median_mda": np.nan,
                    "adjacent_gap_p95_mda": np.nan,
                    "adjacent_gaps_within_2mDa": np.nan,
                }
            gaps = np.diff(mzs)
            steps = np.maximum(1, np.rint(gaps / expected).astype(int))
            errors = np.abs(gaps - steps * expected) * 1000
            return {
                "n_adjacent_gaps": len(gaps),
                "adjacent_gap_median_mda": float(np.median(errors)),
                "adjacent_gap_p95_mda": float(np.quantile(errors, 0.95)),
                "adjacent_gaps_within_2mDa": float(
                    np.mean(errors <= MZ_TOLERANCE_MDA)
                ),
            }

        spacing_a = side_spacing(a)
        spacing_b = side_spacing(b)

        # Trace-level RT overlap: compare every trace in A with every trace in B.
        rt_overlaps = []
        rt_midpoint_diffs = []
        trace_pairs_common_scans = []

        for _, ta in a.iterrows():
            ga = df[
                (df["pair_id"] == pair_id)
                & (df["pair_side"] == "A")
                & (df["trace_index"].astype(str) == str(ta["trace_index"]))
            ]
            scans_a = set(ga["scan_index"].unique())

            for _, tb in b.iterrows():
                gb = df[
                    (df["pair_id"] == pair_id)
                    & (df["pair_side"] == "B")
                    & (df["trace_index"].astype(str) == str(tb["trace_index"]))
                ]
                scans_b = set(gb["scan_index"].unique())

                lo = max(ta["rt_min_sec"], tb["rt_min_sec"])
                hi = min(ta["rt_max_sec"], tb["rt_max_sec"])
                overlap = max(0.0, hi - lo)
                union_span = max(
                    ta["rt_max_sec"], tb["rt_max_sec"]
                ) - min(ta["rt_min_sec"], tb["rt_min_sec"])
                overlap_fraction = (
                    overlap / union_span if union_span > 0 else np.nan
                )

                common_scans = len(scans_a & scans_b)
                rt_mid_diff = abs(
                    ta["rt_mid_sec"] - tb["rt_mid_sec"]
                )

                rt_overlaps.append(overlap_fraction)
                rt_midpoint_diffs.append(rt_mid_diff)
                trace_pairs_common_scans.append(common_scans)

        # Assess whether traces across both groups can form a single ladder.
        combined = all_g.sort_values("mz_median")
        mzs = combined["mz_median"].dropna().to_numpy()
        combined_errors = []
        if len(mzs) >= 2 and np.isfinite(expected):
            gaps = np.diff(mzs)
            steps = np.maximum(1, np.rint(gaps / expected).astype(int))
            combined_errors = (
                np.abs(gaps - steps * expected) * 1000
            ).tolist()

        combined_fraction = (
            float(np.mean(np.asarray(combined_errors) <= MZ_TOLERANCE_MDA))
            if combined_errors else np.nan
        )
        combined_median_error = (
            float(np.median(combined_errors))
            if combined_errors else np.nan
        )

        # Conservative, explicitly heuristic labels.
        if not np.isfinite(combined_fraction):
            label = "inconclusive"
        elif (
            combined_fraction >= 0.8
            and np.isfinite(combined_median_error)
            and combined_median_error <= MZ_TOLERANCE_MDA
            and np.nanmedian(rt_overlaps) >= 0.5
            and max(trace_pairs_common_scans, default=0) >= MIN_COMMON_SCANS
        ):
            label = "potential_split_ladder"
        elif (
            spacing_a["n_adjacent_gaps"] >= 1
            and spacing_b["n_adjacent_gaps"] >= 1
            and spacing_a["adjacent_gaps_within_2mDa"] >= 0.8
            and spacing_b["adjacent_gaps_within_2mDa"] >= 0.8
            and combined_fraction < 0.8
        ):
            label = "potential_separate_ladders"
        else:
            label = "inconclusive"

        pair_rows.append({
            "pair_id": pair_id,
            "acquisition": all_g["acquisition"].iloc[0],
            "feature_id": all_g["feature_id"].iloc[0],
            "charge": charge,
            "expected_isotope_spacing_mz": expected,
            "n_traces_A": len(a),
            "n_traces_B": len(b),
            "mz_range_A": a["mz_median"].max() - a["mz_median"].min(),
            "mz_range_B": b["mz_median"].max() - b["mz_median"].min(),
            "rt_overlap_fraction_median": (
                float(np.nanmedian(rt_overlaps)) if rt_overlaps else np.nan
            ),
            "rt_midpoint_difference_median_sec": (
                float(np.nanmedian(rt_midpoint_diffs))
                if rt_midpoint_diffs else np.nan
            ),
            "common_scan_count_max": (
                max(trace_pairs_common_scans)
                if trace_pairs_common_scans else 0
            ),
            "combined_ladder_gap_count": len(combined_errors),
            "combined_ladder_gap_error_median_mDa": combined_median_error,
            "combined_ladder_gaps_within_2mDa_fraction": combined_fraction,
            **{f"A_{k}": v for k, v in spacing_a.items()},
            **{f"B_{k}": v for k, v in spacing_b.items()},
            "heuristic_relationship": label,
            "interpretation_limit": (
                "Diagnostic only; not a validated IG merge or mono-anchor call."
            ),
        })

    result = pd.DataFrame(pair_rows)
    result.to_csv(OUT_PAIRS, sep="\t", index=False)

    print(f"Trace audit: {OUT_TRACES}")
    print(f"Pair audit:  {OUT_PAIRS}")
    print(f"Traces audited: {len(traces)}")
    print(f"Pairs audited: {len(result)}")

    print("\n=== Heuristic relationship counts ===")
    print(
        result["heuristic_relationship"]
        .value_counts(dropna=False)
        .to_string()
    )

    cols = [
        "acquisition", "feature_id", "charge",
        "n_traces_A", "n_traces_B",
        "rt_overlap_fraction_median",
        "common_scan_count_max",
        "combined_ladder_gap_error_median_mDa",
        "combined_ladder_gaps_within_2mDa_fraction",
        "heuristic_relationship",
    ]
    print("\n=== Pair audit preview ===")
    print(result[[c for c in cols if c in result.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
