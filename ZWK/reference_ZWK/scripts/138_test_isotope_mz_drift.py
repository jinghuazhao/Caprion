#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

INPUT = Path(
    "reports/feature10275997325121874585_raw_isotope_spacing_v1.tsv"
)
OUT_PEAKS = Path(
    "reports/feature10275997325121874585_mz_drift_peak_fits_v1.tsv"
)
OUT_SCAN = Path(
    "reports/feature10275997325121874585_mz_drift_residuals_v1.tsv"
)
OUT_SUMMARY = Path(
    "reports/feature10275997325121874585_mz_drift_summary_v1.tsv"
)

RT_MIN, RT_MAX = 755.0, 769.0
MZ_COLS = [f"peak{i}_mz" for i in range(1, 5)]


def fit_line(x, y):
    """OLS fit y = intercept + slope*x; slope reported per second."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return None
    slope, intercept = np.polyfit(x, y, 1)
    pred = intercept + slope * x
    resid = y - pred
    ss_res = np.sum(resid**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {
        "n": len(x),
        "intercept": intercept,
        "slope_Da_per_s": slope,
        "slope_mDa_per_s": slope * 1000,
        "slope_mDa_per_min": slope * 60000,
        "r_squared": r2,
        "residual_sd_mDa": np.std(resid, ddof=1) * 1000,
        "residual_mad_mDa": np.median(
            np.abs(resid - np.median(resid))
        ) * 1000,
        "max_abs_residual_mDa": np.max(np.abs(resid)) * 1000,
        "rt_min": x.min(),
        "rt_max": x.max(),
    }


def main():
    d = pd.read_csv(INPUT, sep="\t")
    d["rt_seconds"] = pd.to_numeric(d["rt_seconds"], errors="coerce")
    for col in MZ_COLS:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    d = d[d["rt_seconds"].between(RT_MIN, RT_MAX)].copy()
    d = d.sort_values("rt_seconds")

    if d.empty:
        raise SystemExit("No rows found in requested RT range.")

    x = d["rt_seconds"].to_numpy(float)
    fits = []
    fitted_values = {}

    for i, col in enumerate(MZ_COLS, start=1):
        y = d[col].to_numpy(float)
        fit = fit_line(x, y)
        if fit is None:
            continue
        fit["peak"] = i
        fit["mz_column"] = col
        fits.append(fit)

        predicted = fit["intercept"] + fit["slope_Da_per_s"] * x
        residual = (y - predicted) * 1000
        fitted_values[i] = (predicted, residual)

        d[f"peak{i}_fitted_mz"] = predicted
        d[f"peak{i}_linear_residual_mDa"] = residual

    fits_df = pd.DataFrame(fits)
    fits_df.to_csv(OUT_PEAKS, sep="\t", index=False, na_rep="NA")
    d.to_csv(OUT_SCAN, sep="\t", index=False, na_rep="NA")

    summary = []
    if len(fits_df):
        slopes = fits_df["slope_mDa_per_s"].to_numpy(float)
        summary.extend([
            {"metric": "n_peaks_fitted", "value": len(fits_df)},
            {"metric": "mean_slope_mDa_per_s", "value": np.mean(slopes)},
            {"metric": "slope_sd_across_peaks_mDa_per_s",
             "value": np.std(slopes, ddof=1) if len(slopes) > 1 else np.nan},
            {"metric": "slope_range_mDa_per_s",
             "value": np.ptp(slopes)},
            {"metric": "common_drift_slope_mDa_per_s",
             "value": np.median(slopes)},
        ])

    # Relative drift: subtract the mean of the four peak positions in each
    # scan. This removes shared movement and reveals differential movement.
    matrix = d[MZ_COLS].to_numpy(float)
    row_mean = np.nanmean(matrix, axis=1)
    rel_rows = []
    for i, col in enumerate(MZ_COLS):
        rel = (matrix[:, i] - row_mean) * 1000
        ok = np.isfinite(x) & np.isfinite(rel)
        if ok.sum() >= 3:
            slope, intercept = np.polyfit(x[ok], rel[ok], 1)
            rel_rows.append({
                "peak": i + 1,
                "relative_slope_mDa_per_s": slope,
                "relative_slope_mDa_per_min": slope * 60,
                "relative_position_start_mDa": (
                    np.nanmedian(rel[x <= x.min() + 1.5])
                ),
                "relative_position_end_mDa": (
                    np.nanmedian(rel[x >= x.max() - 1.5])
                ),
            })

    rel_df = pd.DataFrame(rel_rows)
    if not rel_df.empty:
        rel_df.to_csv(
            OUT_SUMMARY, sep="\t", index=False, na_rep="NA"
        )
    else:
        pd.DataFrame(columns=[
            "peak", "relative_slope_mDa_per_s",
            "relative_slope_mDa_per_min",
            "relative_position_start_mDa",
            "relative_position_end_mDa",
        ]).to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"RT window: {RT_MIN}-{RT_MAX} s")
    print(f"Scans analysed: {len(d)}")
    print("\nIndividual linear fits:")
    print(fits_df[[
        "peak", "n", "slope_mDa_per_s", "slope_mDa_per_min",
        "r_squared", "residual_sd_mDa", "max_abs_residual_mDa"
    ]].to_string(index=False, float_format=lambda z: f"{z:.5f}"))

    print("\nRelative drift after subtracting each scan's mean m/z:")
    print(rel_df.to_string(index=False, float_format=lambda z: f"{z:.5f}"))

    print("\nOutputs:")
    print(f"  {OUT_PEAKS}")
    print(f"  {OUT_SCAN}")
    print(f"  {OUT_SUMMARY}")
    print("\nInterpretation: similar slopes suggest shared m/z drift; "
          "they do not alone prove a common precursor.")


if __name__ == "__main__":
    main()
