#!/usr/bin/env python3
"""Stratify scan-level anchor stability by isotope evidence quality."""

from pathlib import Path
import pandas as pd

BASE = Path("openms_pilot")
STABILITY = BASE / "scan_level_anchor_stability.tsv"
BEST = BASE / "anchor_mass_accuracy_best.tsv"
PER_SCAN = BASE / "anchor_mass_accuracy_per_scan.tsv"

OUT = BASE / "stability_by_evidence.tsv"
OUT_SUMMARY = BASE / "stability_by_evidence_summary.tsv"


def main():
    stability = pd.read_csv(STABILITY, sep="\t", low_memory=False)
    per_scan = pd.read_csv(PER_SCAN, sep="\t", low_memory=False)

    keys = [
        "acquisition",
        "Spectrum",
        "feature_id",
        "charge",
        "nominated_isotope_offset",
    ]

    required = keys + [
        "scan_index",
        "anchor_offset",
        "isotope_positions_found",
        "consecutive_isotopes_from_anchor",
        "n_individual_peak_matches",
        "median_abs_individual_peak_error_ppm",
        "mass_accuracy_score",
    ]
    missing = set(required) - set(per_scan.columns)
    if missing:
        raise SystemExit(
            "Missing columns in per-scan table: "
            + ", ".join(sorted(missing))
        )

    # Summarise evidence quality across candidate/scan/offset rows.
    per_scan["isotope_positions_found"] = pd.to_numeric(
        per_scan["isotope_positions_found"], errors="coerce"
    )
    per_scan["consecutive_isotopes_from_anchor"] = pd.to_numeric(
        per_scan["consecutive_isotopes_from_anchor"], errors="coerce"
    )
    per_scan["n_individual_peak_matches"] = pd.to_numeric(
        per_scan["n_individual_peak_matches"], errors="coerce"
    )
    per_scan["median_abs_individual_peak_error_ppm"] = pd.to_numeric(
        per_scan["median_abs_individual_peak_error_ppm"], errors="coerce"
    )

    evidence = (
        per_scan.groupby(keys, dropna=False)
        .agg(
            mean_positions_found=("isotope_positions_found", "mean"),
            mean_consecutive_isotopes=(
                "consecutive_isotopes_from_anchor", "mean"
            ),
            mean_individual_peak_matches=(
                "n_individual_peak_matches", "mean"
            ),
            median_abs_peak_error_ppm=(
                "median_abs_individual_peak_error_ppm", "median"
            ),
            mean_score=("mass_accuracy_score", "mean"),
        )
        .reset_index()
    )

    out = stability.merge(evidence, on=keys, how="left", validate="one_to_one")

    # Descriptive evidence strata; these are not validated confidence classes.
    out["evidence_stratum"] = "weak_or_incomplete"
    strong = (
        (out["mean_positions_found"] >= 4)
        & (out["mean_consecutive_isotopes"] >= 3)
        & (out["mean_individual_peak_matches"] >= 4)
        & (out["median_abs_peak_error_ppm"] <= 5)
    )
    out.loc[strong, "evidence_stratum"] = "stronger_mass_and_series_evidence"

    intermediate = (
        (out["mean_positions_found"] >= 3)
        & (out["mean_consecutive_isotopes"] >= 2)
        & (out["mean_individual_peak_matches"] >= 3)
        & ~strong
    )
    out.loc[intermediate, "evidence_stratum"] = "intermediate_evidence"

    out.to_csv(OUT, sep="\t", index=False)

    summary = (
        out.groupby("evidence_stratum", dropna=False)
        .agg(
            candidates=("Spectrum", "size"),
            median_scans=("scans_with_winner", "median"),
            median_modal_offset_fraction=("modal_offset_fraction", "median"),
            median_nominated_offset_fraction=(
                "nominated_offset_fraction", "median"
            ),
            median_scan_tie_fraction=("scan_tie_fraction", "median"),
            median_abs_peak_error_ppm=("median_abs_peak_error_ppm", "median"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print("Candidates:", len(out))
    print("\nSummary by evidence stratum:")
    print(summary.to_string(index=False))
    print("\nWrote:", OUT)
    print("Wrote:", OUT_SUMMARY)


if __name__ == "__main__":
    main()
