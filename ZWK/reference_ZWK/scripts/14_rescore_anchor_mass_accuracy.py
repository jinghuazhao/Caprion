#!/usr/bin/env python3
"""Rescore isotope anchors using individual peak mass accuracy."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")

INPUT = BASE / "monoisotopic_anchor_comparison_per_scan_fast.tsv"
OUT_PER_SCAN = BASE / "anchor_mass_accuracy_per_scan.tsv"
OUT_SUMMARY = BASE / "anchor_mass_accuracy_summary.tsv"
OUT_BEST = BASE / "anchor_mass_accuracy_best.tsv"

# Gaussian scale for individual peak errors, in ppm of theoretical m/z.
# This is a tunable scoring parameter, not a validated probability.
MASS_ERROR_SCALE_PPM = 5.0


def main():
    df = pd.read_csv(INPUT, sep="\t", low_memory=False)

    required = {
        "anchor_offset",
        "mean_anchor_score",
        "mean_isotope_positions_found",
        "mean_consecutive_isotopes",
        "median_abs_spacing_error_ppm",
    }
    missing = required.difference(df.columns)
    if missing:
        raise SystemExit(
            "Missing expected columns: " + ", ".join(sorted(missing))
        )

    # Individual peak errors are not present in the original summary.
    # Reuse the original coverage/continuity information, but do not
    # pretend that spacing error is individual peak mass accuracy.
    positions = pd.to_numeric(
        df["mean_isotope_positions_found"], errors="coerce"
    )
    consecutive = pd.to_numeric(
        df["mean_consecutive_isotopes"], errors="coerce"
    )

    df["coverage_fraction"] = positions / 5.0
    df["continuity_fraction"] = consecutive / 5.0

    # Explicitly label this as a provisional score. Individual peak
    # mass accuracy cannot be reconstructed from the summary columns.
    df["provisional_coverage_continuity_score"] = (
        0.55 * df["coverage_fraction"]
        + 0.30 * df["continuity_fraction"]
    )

    # Record why this is not yet an individual-mass-accuracy score.
    df["individual_mass_accuracy_available"] = False

    df.to_csv(OUT_PER_SCAN, sep="\t", index=False)

    group_cols = [
        c for c in [
            "acquisition",
            "Spectrum",
            "feature_id",
            "charge",
            "nominated_isotope_offset",
            "anchor_offset",
        ]
        if c in df.columns
    ]

    summary = (
        df.groupby(group_cols, dropna=False)
        .agg(
            mean_coverage_fraction=("coverage_fraction", "mean"),
            mean_continuity_fraction=("continuity_fraction", "mean"),
            provisional_score=(
                "provisional_coverage_continuity_score", "mean"
            ),
            scans=("provisional_coverage_continuity_score", "size"),
        )
        .reset_index()
    )

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    candidate_cols = [
        c for c in [
            "acquisition",
            "Spectrum",
            "feature_id",
            "charge",
            "nominated_isotope_offset",
        ]
        if c in summary.columns
    ]

    best_rows = []
    for _, group in summary.groupby(candidate_cols, dropna=False):
        best_score = group["provisional_score"].max()
        winners = group[
            np.isclose(group["provisional_score"], best_score)
        ].copy()
        winners["number_of_tied_winners"] = len(winners)
        best_rows.append(winners)

    best = (
        pd.concat(best_rows, ignore_index=True)
        if best_rows else pd.DataFrame()
    )
    best.to_csv(OUT_BEST, sep="\t", index=False)

    print(f"Input rows: {len(df):,}")
    print(f"Anchor summary rows: {len(summary):,}")
    print(f"Best rows, retaining ties: {len(best):,}")

    if not best.empty:
        print("\nTied-winner counts:")
        print(best["number_of_tied_winners"].value_counts().sort_index())

    print("\nImportant: this is not yet an individual-mass-accuracy score.")
    print("The existing input table does not contain individual peak ppm errors.")
    print("Outputs:")
    print(f"  {OUT_PER_SCAN}")
    print(f"  {OUT_SUMMARY}")
    print(f"  {OUT_BEST}")


if __name__ == "__main__":
    main()
