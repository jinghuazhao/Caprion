#!/usr/bin/env python3
"""Fit nested isotope-spacing models for one candidate MS1 feature.

This is a diagnostic only. It does not assign isotope groups and does not
use the original XLSX isotope-group assignments.
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
TRACE_COLUMN = "trace_index"

PER_SCAN_FILE = OUT / "whole_envelope_spacing_v2_per_scan.tsv"
PER_TRACE_FILE = OUT / "whole_envelope_spacing_v2_per_trace.tsv"
SUMMARY_FILE = OUT / "whole_envelope_spacing_v2_model_summary.tsv"


def fit_fixed_spacing(mz, isotope_indices):
    """Fit intercept with theoretical spacing held fixed."""
    offsets = mz - isotope_indices * EXPECTED_SPACING
    intercept = float(np.mean(offsets))
    residuals = mz - (intercept + isotope_indices * EXPECTED_SPACING)
    return intercept, residuals


def fit_free_spacing(mz, isotope_indices):
    """Fit intercept and spacing by ordinary least squares."""
    design = np.column_stack(
        [np.ones(len(isotope_indices)), isotope_indices]
    )
    coefficients, _, _, _ = np.linalg.lstsq(design, mz, rcond=None)
    intercept = float(coefficients[0])
    spacing = float(coefficients[1])
    residuals = mz - (intercept + spacing * isotope_indices)
    return intercept, spacing, residuals


def residual_stats(residuals):
    abs_residuals = np.abs(residuals) * 1000.0
    return {
        "mean_abs_residual_mDa": float(np.mean(abs_residuals)),
        "median_abs_residual_mDa": float(np.median(abs_residuals)),
        "max_abs_residual_mDa": float(np.max(abs_residuals)),
        "fraction_trace_residuals_within_3mDa": float(
            np.mean(abs_residuals <= 3.0)
        ),
    }


def main():
    if not INPUT.exists():
        raise FileNotFoundError(f"Missing input: {INPUT}")

    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT, sep="\t", low_memory=False)

    required = {
        "feature_id", "trace_index", "scan_index", "observed_mz"
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(
            "Input is missing required columns: "
            + ", ".join(sorted(missing))
        )

    df = df.loc[
        df["feature_id"].astype(str).eq(FEATURE_ID)
    ].copy()

    if df.empty:
        raise ValueError(f"Feature not found: {FEATURE_ID}")

    df["trace_index"] = pd.to_numeric(df["trace_index"], errors="coerce")
    df["scan_index"] = pd.to_numeric(df["scan_index"], errors="coerce")
    df["observed_mz"] = pd.to_numeric(df["observed_mz"], errors="coerce")
    df = df.dropna(subset=["trace_index", "scan_index", "observed_mz"])

    # This analysis treats trace 1 as isotope index 0, trace 2 as index 1,
    # etc. That is a hypothesis to test, not an established assignment.
    df["trace_index"] = df["trace_index"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)
    df["isotope_index_assumed"] = df["trace_index"] - 1

    duplicate_keys = df.duplicated(
        subset=["scan_index", "trace_index"], keep=False
    )
    if duplicate_keys.any():
        duplicate_path = OUT / "whole_envelope_spacing_v2_duplicate_keys.tsv"
        df.loc[duplicate_keys].sort_values(
            ["scan_index", "trace_index"]
        ).to_csv(duplicate_path, sep="\t", index=False)
        raise ValueError(
            f"Duplicate scan/trace keys found; inspect {duplicate_path}"
        )

    by_scan = {
        int(scan): group.sort_values("trace_index")
        for scan, group in df.groupby("scan_index")
    }

    eligible_scans = {}
    for n in MODELS:
        eligible_scans[n] = sorted(
            scan
            for scan, group in by_scan.items()
            if set(range(1, n + 1)).issubset(
                set(group["trace_index"].astype(int))
            )
        )

    common_scans = set(eligible_scans[4])
    common_scans.intersection_update(eligible_scans[5])
    common_scans.intersection_update(eligible_scans[6])
    common_scans = sorted(common_scans)

    per_scan_rows = []
    per_trace_rows = []

    for n in MODELS:
        scans = eligible_scans[n]

        for scan in scans:
            group = by_scan[scan]
            group = group.loc[
                group["trace_index"].between(1, n)
            ].sort_values("trace_index")

            # Require exactly traces 1..N; eligibility was checked above.
            trace_numbers = group["trace_index"].to_numpy(dtype=int)
            mz = group["observed_mz"].to_numpy(dtype=float)
            isotope_indices = trace_numbers - 1

            fixed_intercept, fixed_residuals = fit_fixed_spacing(
                mz, isotope_indices
            )
            free_intercept, free_spacing, free_residuals = fit_free_spacing(
                mz, isotope_indices
            )

            fixed_stats = residual_stats(fixed_residuals)
            free_stats = residual_stats(free_residuals)
            all_fixed_within = bool(
                np.all(np.abs(fixed_residuals * 1000.0) <= 3.0)
            )
            all_free_within = bool(
                np.all(np.abs(free_residuals * 1000.0) <= 3.0)
            )

            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index": scan,
                "n_traces_model": n,
                "common_scan_for_all_models": scan in common_scans,
                "n_trace_observations": len(group),
                "expected_spacing_Da": EXPECTED_SPACING,
                "fixed_intercept_mz": fixed_intercept,
                **{
                    "fixed_" + key: value
                    for key, value in fixed_stats.items()
                },
                "fixed_all_traces_within_3mDa": all_fixed_within,
                "free_intercept_mz": free_intercept,
                "free_fitted_spacing_Da": free_spacing,
                "free_spacing_delta_mDa": (
                    free_spacing - EXPECTED_SPACING
                ) * 1000.0,
                **{
                    "free_" + key: value
                    for key, value in free_stats.items()
                },
                "free_all_traces_within_3mDa": all_free_within,
            })

            for position, (_, row) in enumerate(group.iterrows()):
                per_trace_rows.append({
                    "feature_id": FEATURE_ID,
                    "scan_index": scan,
                    "n_traces_model": n,
                    "common_scan_for_all_models": scan in common_scans,
                    "trace_index": int(row["trace_index"]),
                    "isotope_index_assumed": int(
                        row["trace_index"] - 1
                    ),
                    "observed_mz": float(row["observed_mz"]),
                    "expected_spacing_Da": EXPECTED_SPACING,
                    "fixed_predicted_mz": (
                        fixed_intercept
                        + isotope_indices[position] * EXPECTED_SPACING
                    ),
                    "fixed_residual_mDa": (
                        fixed_residuals[position] * 1000.0
                    ),
                    "free_predicted_mz": (
                        free_intercept
                        + isotope_indices[position] * free_spacing
                    ),
                    "free_residual_mDa": (
                        free_residuals[position] * 1000.0
                    ),
                })

    scan_df = pd.DataFrame(per_scan_rows)
    trace_df = pd.DataFrame(per_trace_rows)

    if scan_df.empty:
        raise ValueError(
            "No scans contain all required traces for any model."
        )

    scan_df.to_csv(PER_SCAN_FILE, sep="\t", index=False)
    trace_df.to_csv(PER_TRACE_FILE, sep="\t", index=False)

    summary_rows = []
    for n in MODELS:
        model_all = scan_df.loc[scan_df["n_traces_model"].eq(n)]
        model_common = model_all.loc[
            model_all["common_scan_for_all_models"]
        ]

        for subset_name, subset in (
            ("all_eligible_scans", model_all),
            ("common_scans", model_common),
        ):
            if subset.empty:
                summary_rows.append({
                    "n_traces_model": n,
                    "scan_subset": subset_name,
                    "n_scans": 0,
                })
                continue

            model_trace = trace_df.loc[
                trace_df["n_traces_model"].eq(n)
                & trace_df["common_scan_for_all_models"].eq(
                    subset_name == "common_scans"
                )
            ].copy()

            # For the all-scans summary, include every eligible scan.
            # For the common subset, include only scans eligible for all models.
            if subset_name == "all_eligible_scans":
                model_trace = trace_df.loc[
                    trace_df["n_traces_model"].eq(n)
                    & ~trace_df["common_scan_for_all_models"]
                ].copy()
                # Add the common scans too: all-eligible means all scans.
                model_trace = pd.concat([
                    model_trace,
                    trace_df.loc[
                        trace_df["n_traces_model"].eq(n)
                        & trace_df["common_scan_for_all_models"]
                    ],
                ], ignore_index=True)

            fixed_abs = np.abs(
                model_trace["fixed_residual_mDa"].to_numpy(dtype=float)
            )
            free_abs = np.abs(
                model_trace["free_residual_mDa"].to_numpy(dtype=float)
            )

            summary_rows.append({
                "n_traces_model": n,
                "scan_subset": subset_name,
                "n_scans": len(subset),
                "n_trace_residuals": len(model_trace),
                "median_scan_median_abs_fixed_mDa": float(
                    subset["fixed_median_abs_residual_mDa"].median()
                ),
                "median_scan_max_abs_fixed_mDa": float(
                    subset["fixed_max_abs_residual_mDa"].median()
                ),
                "max_fixed_abs_residual_mDa": float(fixed_abs.max()),
                "fraction_fixed_trace_residuals_within_3mDa": float(
                    np.mean(fixed_abs <= 3.0)
                ),
                "fraction_scans_all_fixed_residuals_within_3mDa": float(
                    subset["fixed_all_traces_within_3mDa"].mean()
                ),
                "median_fitted_spacing_Da": float(
                    subset["free_fitted_spacing_Da"].median()
                ),
                "median_fitted_spacing_delta_mDa": float(
                    subset["free_spacing_delta_mDa"].median()
                ),
                "median_scan_median_abs_free_mDa": float(
                    subset["free_median_abs_residual_mDa"].median()
                ),
                "max_free_abs_residual_mDa": float(free_abs.max()),
                "fraction_free_trace_residuals_within_3mDa": float(
                    np.mean(free_abs <= 3.0)
                ),
                "fraction_scans_all_free_residuals_within_3mDa": float(
                    subset["free_all_traces_within_3mDa"].mean()
                ),
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SUMMARY_FILE, sep="\t", index=False)

    print(f"Feature: {FEATURE_ID}")
    print(f"Input trace observations: {len(df)}")
    print(f"Distinct scans: {df['scan_index'].nunique()}")
    print(f"Expected charge-2 spacing: {EXPECTED_SPACING:.12f} Da")
    print(f"Scans eligible for traces 1-4: {len(eligible_scans[4])}")
    print(f"Scans eligible for traces 1-5: {len(eligible_scans[5])}")
    print(f"Scans eligible for traces 1-6: {len(eligible_scans[6])}")
    print(f"Common scans across all models: {len(common_scans)}")
    print(f"Per-scan results: {PER_SCAN_FILE}")
    print(f"Per-trace results: {PER_TRACE_FILE}")
    print(f"Model summary: {SUMMARY_FILE}")
    print()
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
