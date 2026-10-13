#!/usr/bin/env python3
"""Audit isotope-anchor winner stability across scans."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "anchor_mass_accuracy_summary.tsv"

OUT_CANDIDATES = BASE / "anchor_stability_per_candidate.tsv"
OUT_SUMMARY = BASE / "anchor_stability_summary.tsv"

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
        "anchor_offset",
        "mean_score",
        "median_abs_individual_peak_error_ppm",
        "mean_coverage",
        "mean_continuity",
        "scans",
    ]
    missing = set(required) - set(df.columns)
    if missing:
        raise SystemExit(
            "Missing columns: " + ", ".join(sorted(missing))
        )

    results = []

    for key, group in df.groupby(KEYS, dropna=False, sort=False):
        group = group.copy()
        scores = pd.to_numeric(group["mean_score"], errors="coerce")
        finite = scores.notna()

        if not finite.any():
            continue

        group = group.loc[finite].copy()
        group["mean_score"] = pd.to_numeric(
            group["mean_score"], errors="coerce"
        )

        best_score = group["mean_score"].max()
        winners = group[
            np.isclose(group["mean_score"], best_score)
        ]

        # Relative score margin between the best and second-best offset.
        ordered = np.sort(group["mean_score"].unique())[::-1]
        margin = (
            float(ordered[0] - ordered[1])
            if len(ordered) > 1 else np.nan
        )

        winner_offsets = sorted(winners["anchor_offset"].astype(int).unique())
        nominated = int(key[-1])

        results.append({
            **dict(zip(KEYS, key)),
            "n_offsets_evaluated": int(group["anchor_offset"].nunique()),
            "n_best_offsets": len(winner_offsets),
            "best_offsets": ",".join(map(str, winner_offsets)),
            "best_score": float(best_score),
            "best_vs_second_score_margin": margin,
            "best_agrees_with_nominated": nominated in winner_offsets,
            "best_median_abs_peak_error_ppm": float(
                pd.to_numeric(
                    winners["median_abs_individual_peak_error_ppm"],
                    errors="coerce",
                ).min()
            ),
            "best_mean_coverage": float(
                pd.to_numeric(winners["mean_coverage"], errors="coerce").max()
            ),
            "best_mean_continuity": float(
                pd.to_numeric(
                    winners["mean_continuity"], errors="coerce"
                ).max()
            ),
            "scans_for_best_offset": int(
                pd.to_numeric(winners["scans"], errors="coerce").max()
            ),
        })

    out = pd.DataFrame(results)
    out["unique_winner"] = out["n_best_offsets"] == 1
    out.to_csv(OUT_CANDIDATES, sep="\t", index=False)

    print(f"Candidates audited: {len(out):,}")
    print("Unique winners:", int(out["unique_winner"].sum()))
    print("Tied winners:", int((~out["unique_winner"]).sum()))

    print("\nAgreement with nominated offset:")
    print(
        out.groupby("unique_winner")["best_agrees_with_nominated"]
        .agg(["count", "sum", "mean"])
    )

    print("\nScore-margin distribution:")
    print(out["best_vs_second_score_margin"].describe())

    print("\nWrote:", OUT_CANDIDATES)

if __name__ == "__main__":
    main()
