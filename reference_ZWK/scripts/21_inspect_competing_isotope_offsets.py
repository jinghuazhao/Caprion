#!/usr/bin/env python3
"""Inspect per-isotope mass errors for competing monoisotopic anchor offsets."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
BEST_FILE = BASE / "strong_evidence_offset_conflicts.tsv"
SCAN_FILE = BASE / "anchor_mass_accuracy_per_scan.tsv"

OUT_DETAIL = BASE / "competing_offset_isotope_detail.tsv"
OUT_SUMMARY = BASE / "competing_offset_isotope_summary.tsv"

MAX_CANDIDATES = 50
ISOTOPES = range(5)

KEYS = [
    "acquisition",
    "Spectrum",
    "feature_id",
    "charge",
    "nominated_isotope_offset",
]


def main():
    best = pd.read_csv(BEST_FILE, sep="\t", low_memory=False)
    scans = pd.read_csv(SCAN_FILE, sep="\t", low_memory=False)

    for name, df in [("best", best), ("scans", scans)]:
        missing = [c for c in KEYS if c not in df.columns]
        if missing:
            raise ValueError(f"{name} table is missing keys: {missing}")

    if "anchor_offset" not in best.columns:
        raise ValueError("Winning-offset table lacks anchor_offset.")

    # Select candidates with multiple winning offsets.
    offset_counts = (
        best.groupby(KEYS, dropna=False)["anchor_offset"]
        .nunique()
        .rename("n_winning_offsets")
        .reset_index()
    )
    competing = offset_counts.query("n_winning_offsets > 1")

    if competing.empty:
        print("No candidates with competing winning offsets.")
        return

    # Prefer candidates with high scan-level ambiguity, when available.
    if "scan_tie_fraction" in best.columns:
        tie_col = "scan_tie_fraction"
    elif "mean_scan_tie_fraction" in best.columns:
        tie_col = "mean_scan_tie_fraction"
    else:
        tie_col = None

    if tie_col:
        tie_values = best.groupby(KEYS, dropna=False)[tie_col].max().reset_index()
        competing = competing.merge(tie_values, on=KEYS, how="left")
        competing = competing.sort_values(
            [tie_col, "n_winning_offsets"], ascending=False
        )
    else:
        competing = competing.sort_values("n_winning_offsets", ascending=False)

    competing = competing.head(MAX_CANDIDATES)

    selected = scans.merge(competing[KEYS], on=KEYS, how="inner")

    required = ["anchor_offset", "psm_mz", "charge"]
    missing = [c for c in required if c not in selected.columns]
    if missing:
        raise ValueError(f"Scan table is missing required columns: {missing}")

    # Reshape M0...M4 observed m/z and ppm-error columns to long format.
    records = []
    for i in ISOTOPES:
        mz_col = f"M{i}_observed_mz"
        err_col = f"M{i}_ppm_error"

        if mz_col not in selected.columns or err_col not in selected.columns:
            continue

        part = selected[KEYS + [
            "scan_index", "rt_sec", "psm_rt_sec", "psm_mz",
            "anchor_offset", "anchor_mz", "anchor_score",
            "n_individual_peak_matches", "median_abs_individual_peak_error_ppm",
            mz_col, err_col,
        ]].copy()

        part = part.rename(columns={
            mz_col: "observed_mz",
            err_col: "individual_peak_error_ppm",
        })
        part["isotope_index"] = i
        part["theoretical_mz"] = (
            part["psm_mz"]
            - part["anchor_offset"] * (
                1.00335483507 / part["charge"].astype(float)
            )
            + i * (
                1.00335483507 / part["charge"].astype(float)
            )
        )
        records.append(part)

    if not records:
        raise ValueError("No M0...M4 observed-m/z and ppm-error columns found.")

    detail = pd.concat(records, ignore_index=True)
    detail["observed_mz"] = pd.to_numeric(detail["observed_mz"], errors="coerce")
    detail["individual_peak_error_ppm"] = pd.to_numeric(
        detail["individual_peak_error_ppm"], errors="coerce"
    )
    detail["absolute_error_ppm"] = detail["individual_peak_error_ppm"].abs()
    detail["anchor_agrees_with_nominated"] = (
        detail["anchor_offset"] == detail["nominated_isotope_offset"]
    )

    detail.to_csv(OUT_DETAIL, sep="\t", index=False)

    group_keys = KEYS + ["anchor_offset"]
    summary = (
        detail.groupby(group_keys, dropna=False)
        .agg(
            scans=("scan_index", "nunique"),
            isotope_observations=("observed_mz", "count"),
            median_abs_error_ppm=("absolute_error_ppm", "median"),
            max_abs_error_ppm=("absolute_error_ppm", "max"),
            median_anchor_score=("anchor_score", "median"),
            median_individual_matches=("n_individual_peak_matches", "median"),
            median_m0_observed_mz=("observed_mz", "median"),
        )
        .reset_index()
    )

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print("Competing candidates selected:", len(competing))
    print("Long-format isotope observations:", len(detail))
    print("Candidate-offset summaries:", len(summary))
    print("\nFirst 20 summaries:")
    print(summary.head(20).to_string(index=False))
    print("\nWrote:", OUT_DETAIL)
    print("Wrote:", OUT_SUMMARY)


if __name__ == "__main__":
    main()
