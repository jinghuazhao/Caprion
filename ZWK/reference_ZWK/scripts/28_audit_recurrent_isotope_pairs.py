#!/usr/bin/env python3
"""Audit recurrence of candidate isotope relationships across scans."""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUTS = [
    ("first_pilot", BASE / "isotope_sequence_pair_audit_first_pilot.tsv"),
    (
        "second_acquisition",
        BASE / "isotope_sequence_pair_audit_second_acquisition.tsv",
    ),
]

OUT = BASE / "recurrent_isotope_pairs_v2.tsv"
SUMMARY = BASE / "recurrent_isotope_pair_summary_v2.tsv"


def main():
    results = []
    summaries = []

    for acquisition, path in INPUTS:
        print(f"Reading {path}")
        d = pd.read_csv(path, sep="\t")

        required = {
            "feature_id", "charge", "scan_index",
            "lower_mz", "upper_mz", "isotope_steps",
            "abs_gap_deviation_mDa",
            "both_peaks_positive_intensity",
        }
        missing = required - set(d.columns)
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")

        # Group observations by feature and candidate peak pair.
        # Round m/z to 5 decimal places to absorb minor scan-to-scan
        # centroid fluctuations without using a broad matching tolerance.
        d["lower_mz_key"] = d["lower_mz"].round(5)
        d["upper_mz_key"] = d["upper_mz"].round(5)

        keys = [
            "feature_id", "charge", "lower_mz_key",
            "upper_mz_key", "isotope_steps",
        ]

        agg = d.groupby(keys, dropna=False).agg(
            n_pair_observations=("scan_index", "size"),
            n_scans=("scan_index", "nunique"),
            median_abs_gap_mDa=("abs_gap_deviation_mDa", "median"),
            p95_abs_gap_mDa=(
                "abs_gap_deviation_mDa",
                lambda x: x.quantile(0.95),
            ),
            positive_intensity_fraction=(
                "both_peaks_positive_intensity", "mean"
            ),
            median_lower_mz=("lower_mz", "median"),
            median_upper_mz=("upper_mz", "median"),
        ).reset_index()

        agg.insert(0, "acquisition", acquisition)

        # Recurrence categories are descriptive, not confidence scores.
        agg["recurrence_class"] = pd.cut(
            agg["n_scans"],
            bins=[0, 1, 2, 4, 9, np.inf],
            labels=["1 scan", "2 scans", "3–4 scans", "5–9 scans", "10+ scans"],
        )

        results.append(agg)

        one_step = agg[agg["isotope_steps"] == 1]
        summary = {
            "acquisition": acquisition,
            "candidate_relationships": len(agg),
            "one_step_relationships": len(one_step),
            "one_step_seen_in_1_scan": int((one_step["n_scans"] == 1).sum()),
            "one_step_seen_in_2plus_scans": int((one_step["n_scans"] >= 2).sum()),
            "one_step_seen_in_3plus_scans": int((one_step["n_scans"] >= 3).sum()),
            "one_step_seen_in_5plus_scans": int((one_step["n_scans"] >= 5).sum()),
            "median_scans_per_one_step_relationship": one_step["n_scans"].median(),
            "median_gap_mDa_recurrent_2plus": one_step.loc[
                one_step["n_scans"] >= 2, "median_abs_gap_mDa"
            ].median(),
            "median_gap_mDa_recurrent_3plus": one_step.loc[
                one_step["n_scans"] >= 3, "median_abs_gap_mDa"
            ].median(),
        }
        summaries.append(summary)
        print(pd.Series(summary).to_string())

    combined = pd.concat(results, ignore_index=True)
    combined.to_csv(OUT, sep="\t", index=False)

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(SUMMARY, sep="\t", index=False)

    print(f"\nPair recurrence table: {OUT}")
    print(f"Summary: {SUMMARY}")
    print(f"Relationship rows: {len(combined)}")

    print("\n=== One-step recurrence classes ===")
    print(
        combined[combined["isotope_steps"] == 1]
        .groupby(["acquisition", "recurrence_class"], observed=False)
        .size()
        .rename("relationships")
        .reset_index()
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
