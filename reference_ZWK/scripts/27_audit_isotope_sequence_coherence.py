#!/usr/bin/env python3
"""Audit coherent isotope sequences in existing OpenMS trace-peak tables.

This is a diagnostic audit, not a final monoisotopic-anchor assignment.
It does not read or use the original XLSX isotope-group assignments.
"""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUTS = [
    ("first_pilot", BASE / "isotope_trace_peaks.tsv"),
    ("second_acquisition", BASE / "isotope_trace_peaks_second.tsv"),
]

ISOTOPE_SPACING = 1.00335483507
GAP_TOLERANCE_MDA = 2.0


def main():
    summaries = []
    detail_frames = []

    for acquisition, path in INPUTS:
        print(f"Reading {path}")
        d = pd.read_csv(path, sep="\t")

        required = {
            "feature_id", "feature_charge", "scan_index",
            "observed_mz", "peak_intensity",
        }
        missing = required - set(d.columns)
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")

        d = d.dropna(
            subset=[
                "feature_id", "feature_charge", "scan_index",
                "observed_mz", "peak_intensity",
            ]
        ).copy()

        d["feature_charge"] = pd.to_numeric(
            d["feature_charge"], errors="coerce"
        )
        d["observed_mz"] = pd.to_numeric(
            d["observed_mz"], errors="coerce"
        )
        d["peak_intensity"] = pd.to_numeric(
            d["peak_intensity"], errors="coerce"
        )
        d = d.dropna(
            subset=["feature_charge", "observed_mz", "peak_intensity"]
        )
        d = d[d["feature_charge"] > 0]

        # Deduplicate exact peak observations before sequence construction.
        d = d.drop_duplicates(
            ["feature_id", "scan_index", "observed_mz"]
        )

        # Build within-feature, within-scan candidate peak pairs.
        group_cols = ["feature_id", "feature_charge", "scan_index"]
        records = []

        for (feature_id, charge, scan_index), g in d.groupby(
            group_cols, sort=False
        ):
            z = int(charge)
            if z <= 0 or len(g) < 2:
                continue

            g = g.sort_values("observed_mz")
            mz = g["observed_mz"].to_numpy(dtype=float)
            intensity = g["peak_intensity"].to_numpy(dtype=float)
            expected = ISOTOPE_SPACING / z

            for i in range(len(mz) - 1):
                for j in range(i + 1, len(mz)):
                    gap = mz[j] - mz[i]
                    # Only adjacent-isotope-sized gaps are candidates here.
                    # Larger gaps may represent missing isotope positions.
                    isotope_steps = int(round(gap / expected))
                    if isotope_steps < 1 or isotope_steps > 4:
                        continue

                    deviation_mda = (
                        gap - isotope_steps * expected
                    ) * 1000

                    records.append({
                        "acquisition": acquisition,
                        "feature_id": feature_id,
                        "charge": z,
                        "scan_index": scan_index,
                        "lower_mz": mz[i],
                        "upper_mz": mz[j],
                        "observed_gap_mz": gap,
                        "isotope_steps": isotope_steps,
                        "expected_gap_mz": isotope_steps * expected,
                        "gap_deviation_mDa": deviation_mda,
                        "abs_gap_deviation_mDa": abs(deviation_mda),
                        "lower_intensity": intensity[i],
                        "upper_intensity": intensity[j],
                        "both_peaks_positive_intensity": (
                            intensity[i] > 0 and intensity[j] > 0
                        ),
                    })

        pairs = pd.DataFrame(records)
        detail_path = BASE / f"isotope_sequence_pair_audit_{acquisition}.tsv"
        pairs.to_csv(detail_path, sep="\t", index=False)
        detail_frames.append(pairs)

        if pairs.empty:
            print(f"No candidate pairs found for {acquisition}")
            continue

        direct = pairs[pairs["isotope_steps"] == 1]
        summary = {
            "acquisition": acquisition,
            "input_peak_rows": len(d),
            "feature_scan_groups": d.groupby(group_cols).ngroups,
            "candidate_peak_pairs": len(pairs),
            "one_step_pairs": len(direct),
            "one_step_pairs_within_1mDa": int(
                (direct["abs_gap_deviation_mDa"] <= 1).sum()
            ),
            "one_step_pairs_within_2mDa": int(
                (direct["abs_gap_deviation_mDa"] <= 2).sum()
            ),
            "one_step_median_abs_gap_mDa": direct[
                "abs_gap_deviation_mDa"
            ].median(),
            "one_step_p95_abs_gap_mDa": direct[
                "abs_gap_deviation_mDa"
            ].quantile(0.95),
        }
        summaries.append(summary)
        print(f"Pair detail: {detail_path}")
        print(pd.Series(summary).to_string())

    if summaries:
        summary_path = BASE / "isotope_sequence_coherence_summary.tsv"
        pd.DataFrame(summaries).to_csv(
            summary_path, sep="\t", index=False
        )
        print(f"\nSummary: {summary_path}")

    if detail_frames:
        all_pairs = pd.concat(detail_frames, ignore_index=True)
        all_path = BASE / "isotope_sequence_pair_audit_all.tsv"
        all_pairs.to_csv(all_path, sep="\t", index=False)
        print(f"Combined pair audit: {all_path}")
        print(f"Combined pair rows: {len(all_pairs)}")


if __name__ == "__main__":
    main()
