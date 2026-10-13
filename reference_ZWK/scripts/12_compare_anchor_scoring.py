#!/usr/bin/env python3
"""Compare original and Gaussian-penalty isotope-anchor scores."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "monoisotopic_anchor_comparison_per_scan_fast.tsv"
OUT_SCAN = BASE / "anchor_scoring_comparison_per_scan.tsv"
OUT_SUMMARY = BASE / "anchor_scoring_comparison_summary.tsv"
OUT_BEST = BASE / "anchor_scoring_comparison_best.tsv"
SPACING_SCALE_PPM = 10.0

KEYS = [
    "acquisition", "Spectrum", "feature_id", "charge",
    "nominated_isotope_offset",
]
REQUIRED = KEYS + [
    "scan_index", "anchor_offset", "isotope_positions_found",
    "consecutive_isotopes_from_anchor",
    "median_abs_adjacent_spacing_error_ppm", "anchor_score",
]


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    df = pd.read_csv(INPUT, sep="\t", low_memory=False)
    missing = set(REQUIRED) - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    error = pd.to_numeric(
        df["median_abs_adjacent_spacing_error_ppm"], errors="coerce"
    ).abs()
    positions = pd.to_numeric(df["isotope_positions_found"], errors="coerce")
    consecutive = pd.to_numeric(
        df["consecutive_isotopes_from_anchor"], errors="coerce"
    )

    spacing_quality = np.exp(
        -0.5 * (error / SPACING_SCALE_PPM) ** 2
    )
    spacing_quality = spacing_quality.where(np.isfinite(error), 0.0)

    df["spacing_quality_gaussian"] = spacing_quality
    df["anchor_score_gaussian"] = (
        0.55 * positions / 5.0
        + 0.30 * consecutive / 5.0
        + 0.15 * spacing_quality
    )
    df.to_csv(OUT_SCAN, sep="\t", index=False)

    group_keys = KEYS + ["anchor_offset"]
    summary = df.groupby(group_keys, dropna=False).agg(
        scans_examined=("anchor_score_gaussian", "size"),
        original_mean_score=("anchor_score", "mean"),
        gaussian_mean_score=("anchor_score_gaussian", "mean"),
        median_abs_spacing_error_ppm=(
            "median_abs_adjacent_spacing_error_ppm",
            lambda x: pd.to_numeric(x, errors="coerce").abs().median(),
        ),
        mean_isotope_positions_found=("isotope_positions_found", "mean"),
        mean_consecutive_isotopes=(
            "consecutive_isotopes_from_anchor", "mean"
        ),
    ).reset_index()

    summary["original_rank"] = summary.groupby(KEYS)[
        "original_mean_score"
    ].rank(method="min", ascending=False)
    summary["gaussian_rank"] = summary.groupby(KEYS)[
        "gaussian_mean_score"
    ].rank(method="min", ascending=False)

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    best = summary[summary["gaussian_rank"] == 1].copy()
    best.to_csv(OUT_BEST, sep="\t", index=False)

    original_ties = (
        summary[summary["original_rank"] == 1]
        .groupby(KEYS, dropna=False).size()
    )
    gaussian_ties = best.groupby(KEYS, dropna=False).size()

    print(f"Scan-level rows: {len(df):,}")
    print(f"Candidate/anchor groups: {len(summary):,}")
    print(f"Candidates with original best-score ties: "
          f"{int((original_ties > 1).sum()):,}")
    print(f"Candidates with Gaussian-score ties: "
          f"{int((gaussian_ties > 1).sum()):,}")
    print(f"Candidates with unique Gaussian winner: "
          f"{int((gaussian_ties == 1).sum()):,}")

    mask = (
        (summary["acquisition"] == "second_acquisition")
        & (summary["Spectrum"]
           == "szwk021704i19101xms1.20054.20054.4")
        & (summary["feature_id"].astype(str)
           == "10796541564565835142")
    )
    cols = [
        "anchor_offset", "original_mean_score", "gaussian_mean_score",
        "median_abs_spacing_error_ppm", "original_rank", "gaussian_rank",
    ]
    print("\nExample candidate:")
    print(summary.loc[mask, cols].sort_values("anchor_offset").to_string(index=False))

    print(f"\nWrote {OUT_SCAN}")
    print(f"Wrote {OUT_SUMMARY}")
    print(f"Wrote {OUT_BEST}")


if __name__ == "__main__":
    main()
