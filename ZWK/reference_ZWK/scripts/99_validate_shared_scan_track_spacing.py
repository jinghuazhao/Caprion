#!/usr/bin/env python3
"""Validate adjacent centroid-track spacing using shared scans only."""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

REPORTS = Path("reports")
INPUT = REPORTS / "feature10275997325121874585_centroid_tracks_v1.tsv"
OUT_ROWS = REPORTS / "feature10275997325121874585_shared_scan_track_spacing_v1.tsv"
OUT_SUMMARY = REPORTS / "feature10275997325121874585_shared_scan_track_spacing_summary_v1.tsv"

EXPECTED = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    rows = read_tsv(INPUT)
    if not rows:
        raise RuntimeError(f"No rows in {INPUT}")

    by_track_scan = defaultdict(dict)
    for row in rows:
        track = int(row["track_id"])
        scan = int(row["scan_index_zero_based"])
        by_track_scan[track][scan] = {
            "mz": float(row["centroid_mz"]),
            "intensity": float(row["centroid_intensity"]),
            "native_scan": row["native_scan_number"],
            "rt": float(row["retention_time_seconds"]),
        }

    tracks = sorted(by_track_scan)
    detail = []
    summaries = []

    # Compare every pair of tracks whose median m/z difference is close
    # to one isotope spacing. Require at least 3 shared scans for summary.
    medians = {
        track: median(obs["mz"] for obs in by_track_scan[track].values())
        for track in tracks
    }

    candidate_pairs = []
    for i, lower in enumerate(tracks):
        for upper in tracks[i + 1:]:
            median_spacing = medians[upper] - medians[lower]
            if abs(median_spacing - EXPECTED) <= 0.030:
                candidate_pairs.append((lower, upper))

    for lower, upper in candidate_pairs:
        shared = sorted(
            set(by_track_scan[lower]) & set(by_track_scan[upper])
        )
        residuals = []
        for scan in shared:
            low = by_track_scan[lower][scan]
            high = by_track_scan[upper][scan]
            spacing = high["mz"] - low["mz"]
            residual_mda = (spacing - EXPECTED) * 1000.0
            residuals.append(residual_mda)

            detail.append({
                "lower_track_id": lower,
                "upper_track_id": upper,
                "scan_index_zero_based": scan,
                "native_scan_number": low["native_scan"],
                "retention_time_seconds": low["rt"],
                "lower_centroid_mz": low["mz"],
                "upper_centroid_mz": high["mz"],
                "lower_intensity": low["intensity"],
                "upper_intensity": high["intensity"],
                "spacing_da": spacing,
                "spacing_residual_mDa": residual_mda,
                "within_3_mDa": abs(residual_mda) <= TOLERANCE_MDA,
            })

        within = sum(abs(value) <= TOLERANCE_MDA for value in residuals)
        summaries.append({
            "lower_track_id": lower,
            "upper_track_id": upper,
            "lower_track_median_mz": medians[lower],
            "upper_track_median_mz": medians[upper],
            "median_track_spacing_da": medians[upper] - medians[lower],
            "shared_scans": len(shared),
            "median_shared_scan_spacing_da": (
                median(row["spacing_da"] for row in detail
                       if row["lower_track_id"] == lower
                       and row["upper_track_id"] == upper)
                if shared else None
            ),
            "median_abs_residual_mDa": (
                median(abs(value) for value in residuals)
                if residuals else None
            ),
            "max_abs_residual_mDa": (
                max(abs(value) for value in residuals)
                if residuals else None
            ),
            "within_3_mDa": within,
            "within_3_mDa_fraction": within / len(shared) if shared else None,
        })

    detail_fields = [
        "lower_track_id", "upper_track_id", "scan_index_zero_based",
        "native_scan_number", "retention_time_seconds",
        "lower_centroid_mz", "upper_centroid_mz",
        "lower_intensity", "upper_intensity", "spacing_da",
        "spacing_residual_mDa", "within_3_mDa",
    ]
    with OUT_ROWS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_fields, delimiter="\t")
        writer.writeheader()
        for row in detail:
            writer.writerow({
                key: f"{value:.9f}" if isinstance(value, float) and "mz" in key
                else f"{value:.6f}" if isinstance(value, float)
                else value
                for key, value in row.items()
            })

    summary_fields = [
        "lower_track_id", "upper_track_id",
        "lower_track_median_mz", "upper_track_median_mz",
        "median_track_spacing_da", "shared_scans",
        "median_shared_scan_spacing_da", "median_abs_residual_mDa",
        "max_abs_residual_mDa", "within_3_mDa",
        "within_3_mDa_fraction",
    ]
    with OUT_SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)

    print(f"Input: {INPUT}")
    print(f"Shared-scan detail: {OUT_ROWS}")
    print(f"Pair summary: {OUT_SUMMARY}")
    print(f"Candidate track pairs: {len(candidate_pairs)}")
    print(f"Theoretical spacing: {EXPECTED:.12f} Da")
    print("\nShared-scan pair results:")
    for row in summaries:
        print(
            f'{row["lower_track_id"]}-{row["upper_track_id"]}: '
            f'shared={row["shared_scans"]}, '
            f'median_spacing={row["median_shared_scan_spacing_da"] if row["median_shared_scan_spacing_da"] is not None else float("nan"):.6f} Da, '
            f'median_abs_residual={row["median_abs_residual_mDa"] if row["median_abs_residual_mDa"] is not None else float("nan"):.3f} mDa, '
            f'within_3_mDa={row["within_3_mDa"]}/{row["shared_scans"]}'
        )
    print("\nCaution: track persistence and spacing are supportive evidence, not proof of common precursor identity.")


if __name__ == "__main__":
    main()
