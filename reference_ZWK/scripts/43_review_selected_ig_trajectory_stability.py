#!/usr/bin/env python3
"""Review scan-level m/z and intensity stability for selected candidate IG pairs."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "candidate_isotope_group_selected_scan_trajectories_v2.tsv"
OUT_SCAN = BASE / "candidate_isotope_group_trajectory_stability_scans_v2.tsv"
OUT_PAIR = BASE / "candidate_isotope_group_trajectory_stability_summary_v2.tsv"

REQUIRED = {
    "pair_id", "feature_id", "acquisition", "charge", "scan_index",
    "trace_A", "trace_B", "rt_A_sec", "rt_B_sec",
    "mz_A", "mz_B", "mz_delta_B_minus_A", "intensity_A", "intensity_B",
}


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def robust_mad(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return np.nan
    med = np.median(x)
    return float(np.median(np.abs(x - med)))


def main():
    df = pd.read_csv(INPUT, sep="\t", dtype={
        "pair_id": str, "feature_id": str, "trace_A": str, "trace_B": str
    })

    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    numeric = [
        "charge", "scan_index", "rt_A_sec", "rt_B_sec",
        "mz_A", "mz_B", "mz_delta_B_minus_A",
        "intensity_A", "intensity_B",
    ]
    for col in numeric:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=[
        "charge", "scan_index", "mz_A", "mz_B", "intensity_A", "intensity_B"
    ]).copy()
    df["charge"] = df["charge"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)

    scan_results = []
    pair_results = []

    for pair_id, pg in df.groupby("pair_id", sort=False):
        pg = pg.sort_values("scan_index").copy()
        pg = pg.drop_duplicates("scan_index", keep="first")

        if len(pg) < 3:
            print(f"WARNING: {pair_id}: only {len(pg)} unique scans")

        delta = pg["mz_delta_B_minus_A"].to_numpy(dtype=float)
        intensity_a = pg["intensity_A"].to_numpy(dtype=float)
        intensity_b = pg["intensity_B"].to_numpy(dtype=float)

        delta_med = float(np.nanmedian(delta))
        delta_mad = robust_mad(delta)
        delta_min = float(np.nanmin(delta))
        delta_max = float(np.nanmax(delta))
        delta_range = delta_max - delta_min

        # Fractional deviations from the median m/z separation.
        # These describe stability and do not imply isotope identity.
        deviation = delta - delta_med

        pg["delta_mz_median_for_pair"] = delta_med
        pg["delta_mz_deviation_mDa"] = deviation * 1000.0
        pg["delta_mz_abs_deviation_mDa"] = np.abs(deviation) * 1000.0
        pg["intensity_A_fraction_of_max"] = (
            intensity_a / np.nanmax(intensity_a)
            if np.nanmax(intensity_a) > 0 else np.nan
        )
        pg["intensity_B_fraction_of_max"] = (
            intensity_b / np.nanmax(intensity_b)
            if np.nanmax(intensity_b) > 0 else np.nan
        )

        # Compare relative intensity profiles to reduce the effect of scale.
        corr_raw = safe_corr(intensity_a, intensity_b)
        corr_relative = safe_corr(
            pg["intensity_A_fraction_of_max"],
            pg["intensity_B_fraction_of_max"],
        )

        # Scan-wise log intensity ratio, restricted to positive intensities.
        positive = (intensity_a > 0) & (intensity_b > 0)
        if positive.any():
            log_ratio = np.log2(intensity_b[positive] / intensity_a[positive])
            log_ratio_median = float(np.median(log_ratio))
            log_ratio_mad = robust_mad(log_ratio)
        else:
            log_ratio_median = np.nan
            log_ratio_mad = np.nan

        # A simple descriptive peak-shape comparison. No biological cutoff
        # is applied; values are exported for inspection.
        pair_results.append({
            "pair_id": pair_id,
            "feature_id": str(pg["feature_id"].iloc[0]),
            "acquisition": pg["acquisition"].iloc[0],
            "charge": int(pg["charge"].iloc[0]),
            "trace_A": pg["trace_A"].iloc[0],
            "trace_B": pg["trace_B"].iloc[0],
            "n_unique_common_scans": len(pg),
            "rt_min_sec": float(
                min(pg["rt_A_sec"].min(), pg["rt_B_sec"].min())
            ),
            "rt_max_sec": float(
                max(pg["rt_A_sec"].max(), pg["rt_B_sec"].max())
            ),
            "median_mz_delta_B_minus_A": delta_med,
            "median_abs_mz_delta_deviation_mDa": delta_mad * 1000.0,
            "range_mz_delta_mDa": delta_range * 1000.0,
            "p90_abs_mz_delta_deviation_mDa": float(
                np.nanquantile(np.abs(deviation), 0.90) * 1000.0
            ),
            "fraction_scans_within_1mDa_of_median_delta": float(
                np.nanmean(np.abs(deviation) * 1000.0 <= 1.0)
            ),
            "fraction_scans_within_2mDa_of_median_delta": float(
                np.nanmean(np.abs(deviation) * 1000.0 <= 2.0)
            ),
            "raw_intensity_correlation": corr_raw,
            "normalised_intensity_correlation": corr_relative,
            "median_log2_intensity_ratio_B_over_A": log_ratio_median,
            "mad_log2_intensity_ratio_B_over_A": log_ratio_mad,
            "interpretation_limit": (
                "Descriptive trajectory stability; no identity or merge decision."
            ),
        })

        scan_results.append(pg)

    scan_df = pd.concat(scan_results, ignore_index=True) if scan_results else pd.DataFrame()
    pair_df = pd.DataFrame(pair_results)

    scan_df.to_csv(OUT_SCAN, sep="\t", index=False)
    pair_df.to_csv(OUT_PAIR, sep="\t", index=False)

    print(f"Scan-level review: {OUT_SCAN}")
    print(f"Pair summary:      {OUT_PAIR}")
    print(f"Pairs reviewed:    {len(pair_df)}")
    print(f"Unique scan rows:  {len(scan_df)}")

    if not pair_df.empty:
        cols = [
            "feature_id", "charge", "n_unique_common_scans",
            "median_mz_delta_B_minus_A",
            "median_abs_mz_delta_deviation_mDa",
            "range_mz_delta_mDa",
            "fraction_scans_within_2mDa_of_median_delta",
            "raw_intensity_correlation",
            "normalised_intensity_correlation",
            "mad_log2_intensity_ratio_B_over_A",
        ]
        print("\n=== Trajectory stability summary ===")
        print(pair_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
