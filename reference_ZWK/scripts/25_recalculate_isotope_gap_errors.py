#!/usr/bin/env python3
"""Express isotope-gap deviations in mDa and precursor-m/z ppm."""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "winning_isotope_sequence_examples_v2.tsv"
OUT = BASE / "winning_isotope_gap_errors_v2.tsv"
SUMMARY = BASE / "winning_isotope_gap_error_summary_v2.tsv"

ISOTOPE_SPACING = 1.00335483507


def main():
    d = pd.read_csv(INPUT, sep="\t")

    valid = (
        d["observed_gap_mz"].notna()
        & d["expected_gap_mz"].notna()
        & d["observed_mz"].notna()
        & d["next_observed_mz"].notna()
    )
    g = d.loc[valid].copy()

    g["gap_deviation_mz"] = (
        g["observed_gap_mz"] - g["expected_gap_mz"]
    )
    g["gap_deviation_mDa"] = 1000 * g["gap_deviation_mz"]
    g["gap_deviation_ppm_at_peak_mz"] = (
        g["gap_deviation_mz"]
        / ((g["observed_mz"] + g["next_observed_mz"]) / 2)
        * 1e6
    )
    g["abs_gap_deviation_mDa"] = g["gap_deviation_mDa"].abs()
    g["abs_gap_deviation_ppm_at_peak_mz"] = (
        g["gap_deviation_ppm_at_peak_mz"].abs()
    )

    g.to_csv(OUT, sep="\t", index=False)

    summary = g.groupby(
        ["acquisition", "Spectrum", "feature_id", "charge",
         "best_anchor_offset"],
        dropna=False,
    ).agg(
        n_gaps=("gap_deviation_mDa", "size"),
        median_abs_gap_mDa=("abs_gap_deviation_mDa", "median"),
        max_abs_gap_mDa=("abs_gap_deviation_mDa", "max"),
        median_abs_gap_ppm_at_peak_mz=(
            "abs_gap_deviation_ppm_at_peak_mz", "median"
        ),
        median_abs_individual_peak_error_ppm=(
            "individual_peak_error_ppm",
            lambda x: np.nanmedian(np.abs(x)),
        ),
    ).reset_index()

    summary.to_csv(SUMMARY, sep="\t", index=False)

    print("Gap details:", OUT)
    print("Candidate summary:", SUMMARY)
    print("\nGap deviations by example:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
