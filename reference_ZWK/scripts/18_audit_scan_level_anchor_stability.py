#!/usr/bin/env python3
"""Measure how consistently each isotope anchor wins across scans."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "anchor_mass_accuracy_per_scan.tsv"

OUT_CANDIDATES = BASE / "scan_level_anchor_stability.tsv"
OUT_SUMMARY = BASE / "scan_level_anchor_stability_summary.tsv"

KEYS = [
    "acquisition",
    "Spectrum",
    "feature_id",
    "charge",
    "nominated_isotope_offset",
]

def main():
    df = pd.read_csv(INPUT, sep="\t", low_memory=False)

    required = KEYS + [
        "scan_index",
        "anchor_offset",
        "mass_accuracy_score",
        "median_abs_individual_peak_error_ppm",
        "isotope_positions_found",
        "consecutive_isotopes_from_anchor",
    ]
    missing = set(required) - set(df.columns)
    if missing:
        raise SystemExit(
            "Missing columns: " + ", ".join(sorted(missing))
        )

    df["mass_accuracy_score"] = pd.to_numeric(
        df["mass_accuracy_score"], errors="coerce"
    )
    df["anchor_offset"] = pd.to_numeric(
        df["anchor_offset"], errors="coerce"
    )
    df["scan_index"] = pd.to_numeric(
        df["scan_index"], errors="coerce"
    )

    scan_keys = KEYS + ["scan_index"]
    df = df.dropna(subset=["mass_accuracy_score", "anchor_offset", "scan_index"])

    # Select the best offset separately within each candidate and scan.
    # Preserve ties instead of arbitrarily selecting one.
    scan_winners = []
    for key, group in df.groupby(scan_keys, dropna=False, sort=False):
        best_score = group["mass_accuracy_score"].max()
        winners = group[
            np.isclose(group["mass_accuracy_score"], best_score)
        ]

        for _, winner in winners.iterrows():
            scan_winners.append({
                **dict(zip(scan_keys, key)),
                "winning_offset": int(winner["anchor_offset"]),
                "winning_score": float(winner["mass_accuracy_score"]),
                "n_tied_offsets_in_scan": len(winners),
            })

    winners = pd.DataFrame(scan_winners)

    results = []
    for key, group in winners.groupby(KEYS, dropna=False, sort=False):
        counts = group["winning_offset"].value_counts()
        total_scans = len(group)
        top_count = int(counts.max())
        top_offsets = sorted(counts[counts == top_count].index.astype(int))

        nominated = int(key[-1])
        nominated_count = int(
            (group["winning_offset"] == nominated).sum()
        )

        results.append({
            **dict(zip(KEYS, key)),
            "scans_with_winner": total_scans,
            "modal_winner_offsets": ",".join(map(str, top_offsets)),
            "modal_offset_fraction": top_count / total_scans,
            "nominated_offset_fraction": nominated_count / total_scans,
            "mean_scan_winner_score": group["winning_score"].mean(),
            "scan_tie_fraction": (
                group["n_tied_offsets_in_scan"].gt(1).mean()
            ),
        })

    out = pd.DataFrame(results)
    out.to_csv(OUT_CANDIDATES, sep="\t", index=False)

    print(f"Candidates audited: {len(out):,}")
    print("\nScans with winner:")
    print(out["scans_with_winner"].describe())

    print("\nModal offset fraction:")
    print(out["modal_offset_fraction"].describe())

    print("\nNominated offset fraction:")
    print(out["nominated_offset_fraction"].describe())

    print("\nScan-level tie fraction:")
    print(out["scan_tie_fraction"].describe())

    # Overall summary stratified by number of scans.
    out["scan_count_bin"] = pd.cut(
        out["scans_with_winner"],
        bins=[0, 1, 2, 5, 10, np.inf],
        labels=["1", "2", "3-5", "6-10", ">10"],
    )

    summary = (
        out.groupby("scan_count_bin", observed=False)
        .agg(
            candidates=("scans_with_winner", "size"),
            median_modal_offset_fraction=("modal_offset_fraction", "median"),
            median_nominated_offset_fraction=("nominated_offset_fraction", "median"),
            median_scan_tie_fraction=("scan_tie_fraction", "median"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print("\nSummary by scan count:")
    print(summary.to_string(index=False))
    print("\nWrote:", OUT_CANDIDATES)
    print("Wrote:", OUT_SUMMARY)

if __name__ == "__main__":
    main()
