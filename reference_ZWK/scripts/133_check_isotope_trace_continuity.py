#!/usr/bin/env python3
"""Check intensity co-movement and m/z continuity of four fitted MS1 peaks."""

from pathlib import Path
import numpy as np
import pandas as pd

FIT = Path(
    "reports/feature10275997325121874585_joint_fit_scans_v1.tsv"
)
OUT_SCAN = Path(
    "reports/feature10275997325121874585_trace_continuity_v1.tsv"
)
OUT_CORR = Path(
    "reports/feature10275997325121874585_trace_intensity_correlations_v1.tsv"
)
OUT_SUMMARY = Path(
    "reports/feature10275997325121874585_trace_continuity_summary_v1.tsv"
)

PEAKS = range(1, 5)
WINDOW = 5  # centred rolling median window, in observed scans


def rolling_median(values, window=WINDOW):
    """Robust local trend; interpolate endpoints for complete comparison."""
    return (
        pd.Series(values)
        .rolling(window=window, center=True, min_periods=3)
        .median()
        .interpolate(limit_direction="both")
        .to_numpy(dtype=float)
    )


def main():
    if not FIT.exists():
        raise SystemExit(f"Missing input: {FIT}")

    df = pd.read_csv(FIT, sep="\t").sort_values("rt_seconds").reset_index(drop=True)

    required = {"native_scan_number", "rt_seconds"}
    for i in PEAKS:
        required |= {f"peak{i}_mz", f"peak{i}_intensity"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing columns: {sorted(missing)}")

    # Calculate a local robust trend and residual for each peak position.
    for i in PEAKS:
        mz_col = f"peak{i}_mz"
        intensity_col = f"peak{i}_intensity"
        trend_col = f"peak{i}_smooth_mz"

        df[trend_col] = rolling_median(df[mz_col].to_numpy())
        df[f"peak{i}_mz_residual_mDa"] = (
            df[mz_col] - df[trend_col]
        ) * 1000

        # Relative intensity helps compare isotope profiles with different
        # absolute intensities.
        positive = df[intensity_col].where(df[intensity_col] > 0)
        df[f"peak{i}_log10_intensity"] = np.log10(positive)

    # Correlations between the four isotope-like intensity traces.
    corr_rows = []
    for i in PEAKS:
        for j in range(i + 1, 5):
            a = df[f"peak{i}_intensity"]
            b = df[f"peak{j}_intensity"]
            valid = a.notna() & b.notna() & (a > 0) & (b > 0)

            pearson = a[valid].corr(b[valid], method="pearson")
            spearman = a[valid].corr(b[valid], method="spearman")
            corr_rows.append({
                "peak_pair": f"{i}-{j}",
                "n_scans": int(valid.sum()),
                "pearson_intensity_correlation": pearson,
                "spearman_intensity_correlation": spearman,
            })

    corr = pd.DataFrame(corr_rows)

    # Per-track continuity summary.
    summary_rows = []
    for i in PEAKS:
        residual = df[f"peak{i}_mz_residual_mDa"].dropna()
        mz = df[f"peak{i}_mz"].dropna()
        intensity = df[f"peak{i}_intensity"].dropna()

        summary_rows.append({
            "track": i,
            "n_scans": int(len(mz)),
            "median_mz": float(mz.median()),
            "min_mz": float(mz.min()),
            "max_mz": float(mz.max()),
            "median_abs_smooth_residual_mDa": float(residual.abs().median()),
            "max_abs_smooth_residual_mDa": float(residual.abs().max()),
            "scans_abs_residual_gt_3mDa": int((residual.abs() > 3).sum()),
            "scans_abs_residual_gt_5mDa": int((residual.abs() > 5).sum()),
            "median_intensity": float(intensity.median()),
            "max_intensity": float(intensity.max()),
        })

    # A simple per-scan flag if any track deviates substantially from its
    # rolling-median trend. This is a screening flag, not a rejection rule.
    residual_cols = [f"peak{i}_mz_residual_mDa" for i in PEAKS]
    df["max_abs_track_residual_mDa"] = df[residual_cols].abs().max(axis=1)
    df["any_track_deviation_gt_5mDa"] = (
        df["max_abs_track_residual_mDa"] > 5
    )

    OUT_SCAN.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_SCAN, sep="\t", index=False, na_rep="NA")
    corr.to_csv(OUT_CORR, sep="\t", index=False, na_rep="NA")
    pd.DataFrame(summary_rows).to_csv(
        OUT_SUMMARY, sep="\t", index=False, na_rep="NA"
    )

    print(f"Fitted scans: {len(df)}")
    print("\nIntensity correlations:")
    print(corr.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nTrack continuity summary:")
    print(pd.DataFrame(summary_rows).to_string(
        index=False, float_format=lambda x: f"{x:.4f}"
    ))

    print("\nScans with any track >5 mDa from its local median trend:")
    show = df.loc[df["any_track_deviation_gt_5mDa"], [
        "native_scan_number", "rt_seconds",
        *residual_cols, "max_abs_track_residual_mDa"
    ]]
    if show.empty:
        print("None")
    else:
        print(show.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    print("\nOutputs:")
    print(f"  {OUT_SCAN}")
    print(f"  {OUT_CORR}")
    print(f"  {OUT_SUMMARY}")
    print(
        "\nCaution: only the 14 scans selected by the joint-fit script are "
        "assessed. Correlation and smoothness are supporting diagnostics, "
        "not proof that the four peaks share a precursor."
    )


if __name__ == "__main__":
    main()
