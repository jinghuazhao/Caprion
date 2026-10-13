#!/usr/bin/env python3
"""Audit isotope-spacing consistency for winning v2 anchor assignments.

This is a diagnostic only. It does not modify or validate monoisotopic
assignments and does not use the original XLSX isotope-group assignments.
"""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
BEST_FILE = BASE / "anchor_one_to_one_v2_best.tsv"
SCAN_FILE = BASE / "anchor_one_to_one_v2_per_scan.tsv"

OUT_DETAIL = BASE / "winning_anchor_spacing_audit_v2.tsv"
OUT_SUMMARY = BASE / "winning_anchor_spacing_summary_v2.tsv"

ISOTOPE_SPACING = 1.00335483507
PPM_TOL = 10.0


def main():
    best = pd.read_csv(BEST_FILE, sep="\t")
    scans = pd.read_csv(SCAN_FILE, sep="\t")

    candidate_cols = [
        "acquisition", "Spectrum", "feature_id", "charge",
        "nominated_isotope_offset", "best_anchor_offset",
    ]

    # Restrict the per-scan table to winning offsets only.
    winners = best[candidate_cols].drop_duplicates()

    detail = scans.merge(
        winners,
        on=[
            "acquisition", "Spectrum", "feature_id", "charge",
            "nominated_isotope_offset",
        ],
        how="inner",
        suffixes=("", "_winner"),
    )

    detail = detail[
        detail["anchor_offset"] == detail["best_anchor_offset"]
    ].copy()

    spacing_records = []

    for row in detail.itertuples(index=False):
        z = int(row.charge)
        if z <= 0:
            continue

        expected = ISOTOPE_SPACING / z
        observed_mzs = []

        for k in range(5):
            mz = getattr(row, f"M{k}_observed_mz")
            observed_mzs.append(
                float(mz) if pd.notna(mz) else np.nan
            )

        errors = []
        for k in range(4):
            a, b = observed_mzs[k], observed_mzs[k + 1]
            if np.isfinite(a) and np.isfinite(b):
                gap = b - a
                errors.append(
                    (gap - expected) / expected * 1e6
                )
            else:
                errors.append(np.nan)

        finite_errors = [
            abs(x) for x in errors if np.isfinite(x)
        ]

        spacing_records.append({
            "acquisition": row.acquisition,
            "Spectrum": row.Spectrum,
            "feature_id": row.feature_id,
            "charge": z,
            "best_anchor_offset": int(row.best_anchor_offset),
            "scan_index": int(row.scan_index),
            "isotope_positions_found": int(row.isotope_positions_found),
            "consecutive_isotopes_from_anchor":
                int(row.consecutive_isotopes_from_anchor),
            "median_abs_adjacent_spacing_error_ppm":
                float(np.median(finite_errors))
                if finite_errors else np.nan,
            "max_abs_adjacent_spacing_error_ppm":
                float(np.max(finite_errors))
                if finite_errors else np.nan,
            "n_adjacent_gaps_checked": len(finite_errors),
            **{
                f"M{k}_observed_mz": observed_mzs[k]
                for k in range(5)
            },
            **{
                f"gap_{k}_{k+1}_error_ppm": errors[k]
                for k in range(4)
            },
        })

    result = pd.DataFrame(spacing_records)
    result.to_csv(OUT_DETAIL, sep="\t", index=False)

    # Summarise scan-level spacing diagnostics by winning offset.
    summary = result.groupby(
        "best_anchor_offset", dropna=False
    ).agg(
        scan_rows=("Spectrum", "size"),
        candidates=("Spectrum", "nunique"),
        median_positions=("isotope_positions_found", "median"),
        median_consecutive=("consecutive_isotopes_from_anchor", "median"),
        median_abs_spacing_error_ppm=(
            "median_abs_adjacent_spacing_error_ppm", "median"
        ),
        median_max_gap_error_ppm=(
            "max_abs_adjacent_spacing_error_ppm", "median"
        ),
    ).reset_index()

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print("Detail:", OUT_DETAIL)
    print("Summary:", OUT_SUMMARY)
    print("Detail rows:", len(result))
    print("\nSummary by winning offset:")
    print(summary.to_string(index=False))

    # Diagnostic proportions: adjacent gaps near the expected isotope spacing.
    if len(result):
        valid = result["n_adjacent_gaps_checked"] > 0
        for tol in (10, 50, 100, 500):
            fraction = (
                result.loc[valid, "max_abs_adjacent_spacing_error_ppm"]
                .le(tol).mean()
            )
            print(
                f"Fraction of scans with all checked gaps within "
                f"{tol} ppm: {fraction:.4f}"
            )


if __name__ == "__main__":
    main()
