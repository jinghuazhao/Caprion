#!/usr/bin/env python3
"""Sensitivity analysis of recurrence to m/z rounding resolution."""

from pathlib import Path

import pandas as pd

BASE = Path("openms_pilot")
INPUTS = [
    ("first_pilot", BASE / "isotope_sequence_pair_audit_first_pilot.tsv"),
    (
        "second_acquisition",
        BASE / "isotope_sequence_pair_audit_second_acquisition.tsv",
    ),
]
OUT = BASE / "isotope_pair_mz_resolution_sensitivity.tsv"

# Decimal places tested. Coarser rounding merges more nearby m/z pairs.
DECIMALS = [6, 5, 4, 3]


def main():
    rows = []

    for acquisition, path in INPUTS:
        print(f"Reading {path}")
        d = pd.read_csv(path, sep="\t")

        required = {
            "feature_id", "charge", "scan_index",
            "lower_mz", "upper_mz", "isotope_steps",
            "abs_gap_deviation_mDa",
        }
        missing = required - set(d.columns)
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")

        d = d[d["isotope_steps"] == 1].copy()

        for decimals in DECIMALS:
            x = d.copy()
            x["lower_key"] = x["lower_mz"].round(decimals)
            x["upper_key"] = x["upper_mz"].round(decimals)

            keys = [
                "feature_id", "charge", "lower_key",
                "upper_key",
            ]

            grouped = x.groupby(keys, dropna=False).agg(
                n_scans=("scan_index", "nunique"),
                median_abs_gap_mDa=("abs_gap_deviation_mDa", "median"),
            ).reset_index()

            n_relationships = len(grouped)
            row = {
                "acquisition": acquisition,
                "m_z_decimal_places": decimals,
                "candidate_observations": len(x),
                "unique_relationships": n_relationships,
                "relationships_2plus_scans": int(
                    (grouped["n_scans"] >= 2).sum()
                ),
                "relationships_3plus_scans": int(
                    (grouped["n_scans"] >= 3).sum()
                ),
                "relationships_5plus_scans": int(
                    (grouped["n_scans"] >= 5).sum()
                ),
                "fraction_2plus": (
                    (grouped["n_scans"] >= 2).mean()
                    if n_relationships else float("nan")
                ),
                "fraction_3plus": (
                    (grouped["n_scans"] >= 3).mean()
                    if n_relationships else float("nan")
                ),
                "median_scans_per_relationship": grouped["n_scans"].median(),
                "median_gap_mDa_for_2plus": grouped.loc[
                    grouped["n_scans"] >= 2, "median_abs_gap_mDa"
                ].median(),
                "median_gap_mDa_for_3plus": grouped.loc[
                    grouped["n_scans"] >= 3, "median_abs_gap_mDa"
                ].median(),
            }
            rows.append(row)

    result = pd.DataFrame(rows)
    result.to_csv(OUT, sep="\t", index=False)

    print(f"\nResolution sensitivity results: {OUT}")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
