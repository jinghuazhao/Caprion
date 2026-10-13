#!/usr/bin/env python3
"""Audit internal consistency and overlap of provisional isotope groups."""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "candidate_isotope_groups_v2.tsv"
OUT = BASE / "candidate_isotope_group_consistency_v2.tsv"
OVERLAP_OUT = BASE / "candidate_isotope_group_overlap_v2.tsv"

ISOTOPE_SPACING = 1.00335483507
MAX_STEPS = 4
GAP_TOLERANCE_MDA = 2.0


def main():
    d = pd.read_csv(INPUT, sep="\t")
    required = {
        "acquisition", "candidate_ig_id", "feature_id", "charge",
        "trace_index", "trace_mz_median", "trace_n_scans",
        "candidate_status",
    }
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"Input lacks columns: {sorted(missing)}")

    detail = []
    overlap_rows = []

    group_keys = [
        "acquisition", "candidate_ig_id", "feature_id", "charge"
    ]

    for (acquisition, ig_id, feature_id, charge), g in d.groupby(
        group_keys, sort=False
    ):
        charge = int(charge)
        if charge <= 0:
            continue

        g = g.sort_values("trace_mz_median").reset_index(drop=True)
        mz = g["trace_mz_median"].to_numpy(dtype=float)
        expected_step = ISOTOPE_SPACING / charge

        gaps = np.diff(mz)
        step_numbers = np.rint(gaps / expected_step).astype(int)
        deviations = (
            gaps - step_numbers * expected_step
        ) * 1000.0

        adjacent_valid = (
            (step_numbers >= 1)
            & (step_numbers <= MAX_STEPS)
            & (np.abs(deviations) <= GAP_TOLERANCE_MDA)
        )

        detail.append({
            "acquisition": acquisition,
            "candidate_ig_id": ig_id,
            "feature_id": feature_id,
            "charge": charge,
            "n_traces": len(g),
            "n_adjacent_gaps": len(gaps),
            "n_adjacent_gaps_valid": int(adjacent_valid.sum()),
            "fraction_adjacent_gaps_valid": (
                float(adjacent_valid.mean()) if len(gaps) else np.nan
            ),
            "median_abs_adjacent_gap_deviation_mDa": (
                float(np.median(np.abs(deviations)))
                if len(deviations) else np.nan
            ),
            "max_abs_adjacent_gap_deviation_mDa": (
                float(np.max(np.abs(deviations)))
                if len(deviations) else np.nan
            ),
            "min_trace_mz": float(np.min(mz)),
            "max_trace_mz": float(np.max(mz)),
            "candidate_status": g["candidate_status"].iloc[0],
            "all_adjacent_gaps_valid": bool(
                len(gaps) > 0 and adjacent_valid.all()
            ),
            "monoisotopic_anchor_assigned": False,
        })

    # Detect m/z overlap between different candidate groups within a feature.
    for (acquisition, feature_id, charge), g in d.groupby(
        ["acquisition", "feature_id", "charge"], sort=False
    ):
        group_ranges = g.groupby("candidate_ig_id").agg(
            mz_min=("trace_mz_median", "min"),
            mz_max=("trace_mz_median", "max"),
            n_traces=("trace_index", "nunique"),
        ).reset_index()

        if len(group_ranges) < 2:
            continue

        ranges = group_ranges.to_dict("records")
        for i in range(len(ranges) - 1):
            a = ranges[i]
            for j in range(i + 1, len(ranges)):
                b = ranges[j]

                overlap = max(
                    0.0,
                    min(a["mz_max"], b["mz_max"])
                    - max(a["mz_min"], b["mz_min"]),
                )

                # Also flag groups with very close boundary traces.
                boundary_gap = max(
                    0.0,
                    max(a["mz_min"], b["mz_min"])
                    - min(a["mz_max"], b["mz_max"]),
                )

                if overlap > 0 or boundary_gap <= expected_step * 0.1:
                    overlap_rows.append({
                        "acquisition": acquisition,
                        "feature_id": feature_id,
                        "charge": int(charge),
                        "candidate_ig_id_a": a["candidate_ig_id"],
                        "candidate_ig_id_b": b["candidate_ig_id"],
                        "n_traces_a": a["n_traces"],
                        "n_traces_b": b["n_traces"],
                        "mz_min_a": a["mz_min"],
                        "mz_max_a": a["mz_max"],
                        "mz_min_b": b["mz_min"],
                        "mz_max_b": b["mz_max"],
                        "range_overlap_mz": overlap,
                        "boundary_gap_mz": boundary_gap,
                    })

    result = pd.DataFrame(detail)
    result.to_csv(OUT, sep="\t", index=False)

    overlaps = pd.DataFrame(overlap_rows)
    overlaps.to_csv(OVERLAP_OUT, sep="\t", index=False)

    print(f"Consistency audit: {OUT}")
    print(f"Candidate groups audited: {len(result)}")
    print(f"Group-overlap audit: {OVERLAP_OUT}")
    print(f"Potential overlapping/nearby group pairs: {len(overlaps)}")

    if not result.empty:
        print("\n=== Internal consistency ===")
        print(
            result.groupby(
                ["acquisition", "candidate_status"]
            ).agg(
                groups=("candidate_ig_id", "size"),
                all_gaps_valid=("all_adjacent_gaps_valid", "sum"),
                median_valid_gap_fraction=(
                    "fraction_adjacent_gaps_valid", "median"
                ),
                median_abs_gap_deviation_mDa=(
                    "median_abs_adjacent_gap_deviation_mDa", "median"
                ),
            ).to_string()
        )


if __name__ == "__main__":
    main()
