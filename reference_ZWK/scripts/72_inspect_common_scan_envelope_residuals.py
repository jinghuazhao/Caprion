#!/usr/bin/env python3
"""Inspect nested isotope-envelope fits on scans shared by models 4, 5 and 6.

Diagnostic only: assumes trace_index 1..N corresponds to consecutive
isotope positions. Does not assign isotope groups or use original XLSX data.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
INPUT = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
OUT = ROOT / "reports"

FEATURE_ID = "10275997325121874585"
EXPECTED_SPACING = 1.00335483507 / 2.0
MODELS = (4, 5, 6)
THRESHOLD_MDA = 3.0

PER_SCAN_OUT = OUT / "common_scan_envelope_residuals_v2_per_scan.tsv"
PER_TRACE_OUT = OUT / "common_scan_envelope_residuals_v2_per_trace.tsv"
SUMMARY_OUT = OUT / "common_scan_envelope_residuals_v2_summary.tsv"


def fit_model(mz, isotope_index):
    """Return fixed-spacing and free-spacing fits."""
    fixed_offsets = mz - isotope_index * EXPECTED_SPACING
    fixed_intercept = float(np.mean(fixed_offsets))
    fixed_predicted = fixed_intercept + isotope_index * EXPECTED_SPACING
    fixed_residual = (mz - fixed_predicted) * 1000.0

    design = np.column_stack([np.ones(len(isotope_index)), isotope_index])
    coef, _, _, _ = np.linalg.lstsq(design, mz, rcond=None)
    free_intercept = float(coef[0])
    free_spacing = float(coef[1])
    free_predicted = free_intercept + free_spacing * isotope_index
    free_residual = (mz - free_predicted) * 1000.0

    return {
        "fixed_intercept": fixed_intercept,
        "fixed_predicted": fixed_predicted,
        "fixed_residual_mDa": fixed_residual,
        "free_intercept": free_intercept,
        "free_spacing": free_spacing,
        "free_predicted": free_predicted,
        "free_residual_mDa": free_residual,
    }


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t", low_memory=False)

    required = {"feature_id", "trace_index", "scan_index", "observed_mz"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    df = df.loc[df["feature_id"].astype(str).eq(FEATURE_ID)].copy()
    if df.empty:
        raise ValueError(f"Feature not found: {FEATURE_ID}")

    for col in ("trace_index", "scan_index", "observed_mz"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["trace_index", "scan_index", "observed_mz"])
    df["trace_index"] = df["trace_index"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)

    duplicate = df.duplicated(["scan_index", "trace_index"], keep=False)
    if duplicate.any():
        raise ValueError(
            "Duplicate scan/trace observations exist; resolve them before fitting."
        )

    by_scan = {
        int(scan): group.sort_values("trace_index")
        for scan, group in df.groupby("scan_index")
    }

    eligible = {}
    for n in MODELS:
        eligible[n] = {
            scan for scan, group in by_scan.items()
            if set(range(1, n + 1)).issubset(set(group["trace_index"]))
        }

    common_scans = sorted(set.intersection(*(eligible[n] for n in MODELS)))
    if not common_scans:
        raise ValueError("No scans are common to all three models.")

    scan_rows = []
    trace_rows = []

    for scan in common_scans:
        group = by_scan[scan].set_index("trace_index")

        for n in MODELS:
            model = group.loc[list(range(1, n + 1))]
            trace_indices = model.index.to_numpy(dtype=int)
            mz = model["observed_mz"].to_numpy(dtype=float)
            isotope_index = trace_indices - 1
            fit = fit_model(mz, isotope_index)

            fixed_abs = np.abs(fit["fixed_residual_mDa"])
            free_abs = np.abs(fit["free_residual_mDa"])

            scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index": scan,
                "n_traces_model": n,
                "n_traces": len(model),
                "expected_spacing_Da": EXPECTED_SPACING,
                "fixed_median_abs_residual_mDa": float(np.median(fixed_abs)),
                "fixed_max_abs_residual_mDa": float(np.max(fixed_abs)),
                "fixed_all_within_3mDa": bool(np.all(fixed_abs <= THRESHOLD_MDA)),
                "fixed_traces_over_3mDa": ",".join(
                    str(int(trace_indices[i]))
                    for i in np.where(fixed_abs > THRESHOLD_MDA)[0]
                ),
                "free_fitted_spacing_Da": fit["free_spacing"],
                "free_spacing_delta_mDa": (
                    fit["free_spacing"] - EXPECTED_SPACING
                ) * 1000.0,
                "free_median_abs_residual_mDa": float(np.median(free_abs)),
                "free_max_abs_residual_mDa": float(np.max(free_abs)),
                "free_all_within_3mDa": bool(np.all(free_abs <= THRESHOLD_MDA)),
                "free_traces_over_3mDa": ",".join(
                    str(int(trace_indices[i]))
                    for i in np.where(free_abs > THRESHOLD_MDA)[0]
                ),
            })

            for i, trace_index in enumerate(trace_indices):
                trace_rows.append({
                    "feature_id": FEATURE_ID,
                    "scan_index": scan,
                    "n_traces_model": n,
                    "trace_index": int(trace_index),
                    "isotope_index_assumed": int(trace_index - 1),
                    "observed_mz": mz[i],
                    "fixed_predicted_mz": fit["fixed_predicted"][i],
                    "fixed_residual_mDa": fit["fixed_residual_mDa"][i],
                    "fixed_within_3mDa": (
                        abs(fit["fixed_residual_mDa"][i]) <= THRESHOLD_MDA
                    ),
                    "free_predicted_mz": fit["free_predicted"][i],
                    "free_residual_mDa": fit["free_residual_mDa"][i],
                    "free_within_3mDa": (
                        abs(fit["free_residual_mDa"][i]) <= THRESHOLD_MDA
                    ),
                    "free_fitted_spacing_Da": fit["free_spacing"],
                    "free_spacing_delta_mDa": (
                        fit["free_spacing"] - EXPECTED_SPACING
                    ) * 1000.0,
                })

    scan_df = pd.DataFrame(scan_rows).sort_values(
        ["scan_index", "n_traces_model"]
    )
    trace_df = pd.DataFrame(trace_rows).sort_values(
        ["scan_index", "n_traces_model", "trace_index"]
    )

    scan_df.to_csv(PER_SCAN_OUT, sep="\t", index=False)
    trace_df.to_csv(PER_TRACE_OUT, sep="\t", index=False)

    summary_rows = []
    for n in MODELS:
        part = scan_df.loc[scan_df["n_traces_model"].eq(n)]
        summary_rows.append({
            "n_traces_model": n,
            "n_common_scans": len(part),
            "median_scan_median_abs_fixed_mDa": (
                part["fixed_median_abs_residual_mDa"].median()
            ),
            "median_scan_max_abs_fixed_mDa": (
                part["fixed_max_abs_residual_mDa"].median()
            ),
            "max_fixed_abs_residual_mDa": (
                part["fixed_max_abs_residual_mDa"].max()
            ),
            "scans_all_fixed_within_3mDa": int(
                part["fixed_all_within_3mDa"].sum()
            ),
            "median_fitted_spacing_Da": (
                part["free_fitted_spacing_Da"].median()
            ),
            "median_spacing_delta_mDa": (
                part["free_spacing_delta_mDa"].median()
            ),
            "scans_all_free_within_3mDa": int(
                part["free_all_within_3mDa"].sum()
            ),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SUMMARY_OUT, sep="\t", index=False)

    print(f"Feature: {FEATURE_ID}")
    print(f"Expected charge-2 spacing: {EXPECTED_SPACING:.12f} Da")
    print(f"Common scans: {len(common_scans)}")
    print(f"Common scan indices: {', '.join(map(str, common_scans))}")
    print(f"Per-scan results: {PER_SCAN_OUT}")
    print(f"Per-trace results: {PER_TRACE_OUT}")
    print(f"Summary: {SUMMARY_OUT}")
    print("\nPer-scan model comparison:")
    print(scan_df.to_string(index=False))
    print("\nSummary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
