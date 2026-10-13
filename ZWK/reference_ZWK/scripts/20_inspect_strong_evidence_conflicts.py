#!/usr/bin/env python3
"""Inspect competing offsets for candidates with stronger isotope evidence."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
SUMMARY = BASE / "anchor_mass_accuracy_summary.tsv"
STABILITY = BASE / "stability_by_evidence.tsv"

OUT = BASE / "strong_evidence_offset_conflicts.tsv"
OUT_SUMMARY = BASE / "strong_evidence_offset_conflicts_summary.tsv"

KEYS = [
    "acquisition",
    "Spectrum",
    "feature_id",
    "charge",
    "nominated_isotope_offset",
]

def main():
    summary = pd.read_csv(SUMMARY, sep="\t", low_memory=False)
    stability = pd.read_csv(STABILITY, sep="\t", low_memory=False)

    strong = stability[
        stability["evidence_stratum"]
        == "stronger_mass_and_series_evidence"
    ].copy()

    if strong.empty:
        raise SystemExit("No candidates in the stronger-evidence stratum.")

    merged = summary.merge(
        strong[KEYS + [
            "modal_offset_fraction",
            "nominated_offset_fraction",
            "scan_tie_fraction",
            "evidence_stratum",
        ]],
        on=KEYS,
        how="inner",
        validate="many_to_one",
    )

    # Rank offsets within each candidate, preserving near-equal scores.
    merged["rank"] = merged.groupby(KEYS)["mean_score"].rank(
        method="min", ascending=False
    )

    top = merged[merged["rank"] == 1].copy()
    top["agrees_with_nominated"] = (
        top["anchor_offset"] == top["nominated_isotope_offset"]
    )

    top.to_csv(OUT, sep="\t", index=False)

    result = (
        top.groupby("agrees_with_nominated")
        .agg(
            winning_offset_rows=("anchor_offset", "size"),
            median_score=("mean_score", "median"),
            median_peak_error_ppm=(
                "median_abs_individual_peak_error_ppm", "median"
            ),
            median_coverage=("mean_coverage", "median"),
            median_continuity=("mean_continuity", "median"),
            median_modal_offset_fraction=("modal_offset_fraction", "median"),
            median_scan_tie_fraction=("scan_tie_fraction", "median"),
        )
        .reset_index()
    )

    result.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print("Strong-evidence candidates:", len(strong))
    print("Winning offset rows:", len(top))
    print("\nWinning offsets by agreement with nominated offset:")
    print(result.to_string(index=False))
    print("\nWrote:", OUT)
    print("Wrote:", OUT_SUMMARY)

if __name__ == "__main__":
    main()
