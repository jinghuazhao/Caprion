#!/usr/bin/env python3
"""Audit mass-accuracy anchor winners and compare with nominated offsets."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
BEST = BASE / "anchor_mass_accuracy_best.tsv"
SUMMARY = BASE / "anchor_mass_accuracy_summary.tsv"

OUT = BASE / "mass_accuracy_winner_audit.tsv"
OUT_GROUPS = BASE / "mass_accuracy_winner_audit_summary.tsv"

def main():
    best = pd.read_csv(BEST, sep="\t", low_memory=False)
    summary = pd.read_csv(SUMMARY, sep="\t", low_memory=False)

    keys = [
        "acquisition",
        "Spectrum",
        "feature_id",
        "charge",
        "nominated_isotope_offset",
    ]

    required = keys + [
        "anchor_offset",
        "mean_score",
        "median_abs_individual_peak_error_ppm",
        "mean_coverage",
        "mean_continuity",
        "n_tied_best_offsets",
    ]

    missing = set(required) - set(best.columns)
    if missing:
        raise SystemExit(
            "Missing columns in best table: " + ", ".join(sorted(missing))
        )

    best["agrees_with_nominated"] = (
        best["anchor_offset"] == best["nominated_isotope_offset"]
    )

    # Recalculate a unique-candidate count from the actual grouping keys.
    unique_candidates = best[keys].drop_duplicates().shape[0]

    # A candidate can have multiple winner rows if its best offsets tie.
    candidate_winners = (
        best.groupby(keys, dropna=False)
        .agg(
            n_winner_rows=("anchor_offset", "size"),
            n_tied_best_offsets=("n_tied_best_offsets", "max"),
            agrees_with_nominated=("agrees_with_nominated", "any"),
            best_score=("mean_score", "max"),
            best_median_abs_peak_error_ppm=(
                "median_abs_individual_peak_error_ppm", "min"
            ),
            best_mean_coverage=("mean_coverage", "max"),
            best_mean_continuity=("mean_continuity", "max"),
        )
        .reset_index()
    )

    candidate_winners["unique_winner"] = (
        candidate_winners["n_winner_rows"] == 1
    )

    candidate_winners.to_csv(OUT, sep="\t", index=False)

    print(f"Rows in best table: {len(best):,}")
    print(f"Unique candidates: {unique_candidates:,}")
    print(f"Audited candidate rows: {len(candidate_winners):,}")
    print(
        "Candidates with unique winner:",
        int(candidate_winners["unique_winner"].sum()),
    )
    print(
        "Candidates with ties:",
        int((~candidate_winners["unique_winner"]).sum()),
    )
    print(
        "Candidates whose best offset agrees with nominated offset:",
        int(candidate_winners["agrees_with_nominated"].sum()),
    )

    print("\nTie counts per candidate:")
    print(
        candidate_winners["n_tied_best_offsets"]
        .value_counts()
        .sort_index()
    )

    print("\nAgreement by tie status:")
    print(
        candidate_winners.groupby("unique_winner")[
            "agrees_with_nominated"
        ].agg(["count", "sum", "mean"])
    )

    print("\nWrote:", OUT)

if __name__ == "__main__":
    main()
