#!/usr/bin/env python3
"""Export aligned scan-level trajectories for four selected candidate IG pairs."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "candidate_isotope_group_pair_profiles_v2.tsv"
OUT_LONG = BASE / "candidate_isotope_group_selected_scan_trajectories_v2.tsv"
OUT_SUMMARY = BASE / "candidate_isotope_group_selected_scan_trajectory_summary_v2.tsv"

ISOTOPE_MASS = 1.00335483507

SELECTED_FEATURES = {
    "10275997325121874585": "high_intensity_correlation",
    "7118767039353031980": "low_intensity_correlation",
    "8339255300597246646": "many_common_scans",
    "12167939761234023158": "low_spacing_error",
}


def main():
    df = pd.read_csv(
        INPUT,
        sep="\t",
        dtype={
            "pair_id": str,
            "feature_id": str,
            "trace_index": str,
            "candidate_ig_id": str,
        },
    )

    required = {
        "pair_id", "acquisition", "feature_id", "pair_side",
        "trace_index", "charge", "scan_index", "rt_sec",
        "observed_mz", "peak_intensity",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    for col in [
        "charge", "scan_index", "rt_sec", "observed_mz", "peak_intensity"
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(
        subset=["charge", "scan_index", "rt_sec", "observed_mz"]
    ).copy()
    df["charge"] = df["charge"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)
    df["peak_intensity"] = df["peak_intensity"].fillna(0.0)

    df["feature_id"] = df["feature_id"].astype(str)
    selected = df[df["feature_id"].isin(SELECTED_FEATURES)].copy()

    found = set(selected["feature_id"].unique())
    missing_features = set(SELECTED_FEATURES) - found
    if missing_features:
        raise ValueError(
            "Selected feature IDs were not found in the profile file: "
            + ", ".join(sorted(missing_features))
        )

    # Aggregate multiple peak records from the same trace and scan.
    scan = (
        selected.groupby(
            [
                "pair_id", "acquisition", "feature_id", "pair_side",
                "trace_index", "candidate_ig_id", "charge", "scan_index",
            ],
            dropna=False,
            as_index=False,
        )
        .agg(
            rt_sec=("rt_sec", "median"),
            observed_mz=("observed_mz", "median"),
            peak_intensity=("peak_intensity", "sum"),
        )
    )

    scan["selection_reason"] = scan["feature_id"].map(SELECTED_FEATURES)

    # Build scan-level comparisons for every A/B trace combination.
    trajectory_rows = []
    summary_rows = []

    for pair_id, pg in scan.groupby("pair_id", sort=False):
        feature_id = str(pg["feature_id"].iloc[0])
        acquisition = pg["acquisition"].iloc[0]
        charge = int(pg["charge"].iloc[0])
        spacing = ISOTOPE_MASS / charge

        side_profiles = {}
        for side in ["A", "B"]:
            side_profiles[side] = []
            sg = pg[pg["pair_side"].astype(str).str.upper() == side]

            for trace_index, tg in sg.groupby("trace_index", sort=False):
                side_profiles[side].append(
                    tg.sort_values("scan_index").set_index("scan_index")
                )

        pair_comparisons = []

        for ia, a in enumerate(side_profiles["A"], start=1):
            for ib, b in enumerate(side_profiles["B"], start=1):
                common = a.index.intersection(b.index)
                if len(common) == 0:
                    continue

                a = a.sort_index()
                b = b.sort_index()
                common = common.sort_values()

                # Align by scan index. Retention times are retained separately
                # so possible timing differences remain visible.
                aa = a.loc[common]
                bb = b.loc[common]

                delta = (
                    bb["observed_mz"].to_numpy(dtype=float)
                    - aa["observed_mz"].to_numpy(dtype=float)
                )
                delta_mDa = delta * 1000.0
                nearest_step = int(round(np.nanmedian(delta) / spacing))
                spacing_residual = (delta - nearest_step * spacing) * 1000.0

                intensity_a = aa["peak_intensity"].to_numpy(dtype=float)
                intensity_b = bb["peak_intensity"].to_numpy(dtype=float)

                corr = np.nan
                if (
                    len(common) >= 3
                    and np.std(intensity_a) > 0
                    and np.std(intensity_b) > 0
                ):
                    corr = float(np.corrcoef(intensity_a, intensity_b)[0, 1])

                trace_a = str(aa["trace_index"].iloc[0])
                trace_b = str(bb["trace_index"].iloc[0])

                for k, scan_index in enumerate(common):
                    trajectory_rows.append({
                        "pair_id": pair_id,
                        "feature_id": feature_id,
                        "acquisition": acquisition,
                        "selection_reason": SELECTED_FEATURES[feature_id],
                        "charge": charge,
                        "scan_index": int(scan_index),
                        "trace_A": trace_a,
                        "trace_B": trace_b,
                        "rt_A_sec": float(aa.loc[scan_index, "rt_sec"]),
                        "rt_B_sec": float(bb.loc[scan_index, "rt_sec"]),
                        "rt_delta_B_minus_A_sec": float(
                            bb.loc[scan_index, "rt_sec"]
                            - aa.loc[scan_index, "rt_sec"]
                        ),
                        "mz_A": float(aa.loc[scan_index, "observed_mz"]),
                        "mz_B": float(bb.loc[scan_index, "observed_mz"]),
                        "mz_delta_B_minus_A": float(delta[k]),
                        "mz_delta_mDa": float(delta_mDa[k]),
                        "nearest_integer_spacing_steps": nearest_step,
                        "spacing_residual_mDa": float(spacing_residual[k]),
                        "intensity_A": float(intensity_a[k]),
                        "intensity_B": float(intensity_b[k]),
                        "intensity_correlation_common_scans": corr,
                    })

                pair_comparisons.append({
                    "trace_A": trace_a,
                    "trace_B": trace_b,
                    "common_scan_count": len(common),
                    "median_mz_delta_mDa": float(np.nanmedian(delta_mDa)),
                    "mad_mz_delta_mDa": float(
                        np.nanmedian(np.abs(delta_mDa - np.nanmedian(delta_mDa)))
                    ),
                    "range_mz_delta_mDa": float(
                        np.nanmax(delta_mDa) - np.nanmin(delta_mDa)
                    ),
                    "median_abs_spacing_residual_mDa": float(
                        np.nanmedian(np.abs(spacing_residual))
                    ),
                    "median_abs_rt_delta_sec": float(
                        np.nanmedian(np.abs(
                            bb["rt_sec"].to_numpy()
                            - aa["rt_sec"].to_numpy()
                        ))
                    ),
                    "intensity_correlation": corr,
                })

        comp = pd.DataFrame(pair_comparisons)
        if comp.empty:
            continue

        summary_rows.append({
            "pair_id": pair_id,
            "feature_id": feature_id,
            "acquisition": acquisition,
            "selection_reason": SELECTED_FEATURES[feature_id],
            "charge": charge,
            "n_trace_pairs_compared": len(comp),
            "max_common_scan_count": int(comp["common_scan_count"].max()),
            "median_common_scan_count": float(comp["common_scan_count"].median()),
            "median_abs_mz_delta_mad_mDa": float(comp["mad_mz_delta_mDa"].median()),
            "median_abs_spacing_residual_mDa": float(
                comp["median_abs_spacing_residual_mDa"].median()
            ),
            "median_abs_rt_delta_sec": float(
                comp["median_abs_rt_delta_sec"].median()
            ),
            "median_intensity_correlation": float(
                comp["intensity_correlation"].median()
            ),
            "interpretation_limit": (
                "Scan-level descriptive comparison only; not an IG merge rule."
            ),
        })

    trajectories = pd.DataFrame(trajectory_rows)
    summaries = pd.DataFrame(summary_rows)

    trajectories.to_csv(OUT_LONG, sep="\t", index=False)
    summaries.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"Scan trajectories: {OUT_LONG}")
    print(f"Pair summary:      {OUT_SUMMARY}")
    print(f"Selected features found: {len(found)}")
    print(f"Aligned scan rows: {len(trajectories)}")
    print(f"Pairs summarized:  {len(summaries)}")

    if not summaries.empty:
        print("\n=== Selected pair trajectory summary ===")
        print(summaries.to_string(index=False))


if __name__ == "__main__":
    main()
