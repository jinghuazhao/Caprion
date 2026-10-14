#!/usr/bin/env python3
"""Track native centroid peaks across scans and assess adjacent m/z spacing."""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

REPORTS = Path("reports")
INPUT = REPORTS / "feature10275997325121874585_all_trace_scan_centroids_v1.tsv"
OUT_DETAIL = REPORTS / "feature10275997325121874585_centroid_tracks_v1.tsv"
OUT_SUMMARY = REPORTS / "feature10275997325121874585_centroid_track_summary_v1.tsv"
OUT_PAIRS = REPORTS / "feature10275997325121874585_centroid_track_spacing_v1.tsv"

EXPECTED_SPACING = 1.00335483507 / 2.0
TRACK_TOLERANCE_DA = 0.020
SPACING_TOLERANCE_MDA = 3.0


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    rows = read_tsv(INPUT)
    by_scan = defaultdict(list)

    for row in rows:
        by_scan[int(row["scan_index_zero_based"])].append({
            "scan": int(row["scan_index_zero_based"]),
            "native_scan": row["native_scan_number"],
            "rt": float(row["rt_seconds"]),
            "peak_index": int(row["peak_index"]),
            "mz": float(row["centroid_mz"]),
            "intensity": float(row["centroid_intensity"]),
        })

    # Greedy nearest-m/z tracking across scans. Each scan can contribute
    # at most one centroid to a track.
    tracks = []
    for scan in sorted(by_scan):
        peaks = sorted(by_scan[scan], key=lambda p: p["mz"])
        used_tracks = set()

        for peak in peaks:
            candidates = []
            for track_id, track in enumerate(tracks):
                if track_id in used_tracks:
                    continue
                if track[-1]["scan"] >= scan:
                    continue
                delta = abs(peak["mz"] - track[-1]["mz"])
                if delta <= TRACK_TOLERANCE_DA:
                    candidates.append((delta, track_id))

            if candidates:
                _, track_id = min(candidates)
                tracks[track_id].append(peak)
                used_tracks.add(track_id)
            else:
                tracks.append([peak])
                used_tracks.add(len(tracks) - 1)

    # Sort tracks by median m/z, then assign stable IDs for this report.
    tracks.sort(key=lambda track: median(p["mz"] for p in track))
    track_summaries = []
    detail_rows = []

    for track_id, track in enumerate(tracks, 1):
        mzs = [p["mz"] for p in track]
        intensities = [p["intensity"] for p in track]
        scans = [p["scan"] for p in track]
        summary = {
            "track_id": track_id,
            "observations": len(track),
            "scans_present": len(set(scans)),
            "scans_total": len(by_scan),
            "persistence_fraction": len(set(scans)) / len(by_scan),
            "median_mz": median(mzs),
            "min_mz": min(mzs),
            "max_mz": max(mzs),
            "mz_range_mDa": (max(mzs) - min(mzs)) * 1000.0,
            "median_intensity": median(intensities),
            "max_intensity": max(intensities),
        }
        track_summaries.append(summary)

        for peak in track:
            detail_rows.append({
                "track_id": track_id,
                "scan_index_zero_based": peak["scan"],
                "native_scan_number": peak["native_scan"],
                "retention_time_seconds": peak["rt"],
                "peak_index": peak["peak_index"],
                "centroid_mz": peak["mz"],
                "centroid_intensity": peak["intensity"],
            })

    detail_fields = [
        "track_id", "scan_index_zero_based", "native_scan_number",
        "retention_time_seconds", "peak_index", "centroid_mz",
        "centroid_intensity",
    ]
    with OUT_DETAIL.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_fields, delimiter="\t")
        writer.writeheader()
        for row in detail_rows:
            writer.writerow({
                key: f"{value:.9f}" if isinstance(value, float) and "mz" in key
                else f"{value:.6f}" if isinstance(value, float)
                else value
                for key, value in row.items()
            })

    summary_fields = [
        "track_id", "observations", "scans_present", "scans_total",
        "persistence_fraction", "median_mz", "min_mz", "max_mz",
        "mz_range_mDa", "median_intensity", "max_intensity",
    ]
    with OUT_SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for row in track_summaries:
            writer.writerow(row)

    pair_rows = []
    for lower, upper in zip(track_summaries, track_summaries[1:]):
        spacing = upper["median_mz"] - lower["median_mz"]
        residual_mda = (spacing - EXPECTED_SPACING) * 1000.0
        pair_rows.append({
            "lower_track_id": lower["track_id"],
            "upper_track_id": upper["track_id"],
            "lower_median_mz": lower["median_mz"],
            "upper_median_mz": upper["median_mz"],
            "median_spacing_da": spacing,
            "spacing_residual_mDa": residual_mda,
            "within_3_mDa": abs(residual_mda) <= SPACING_TOLERANCE_MDA,
            "lower_scans_present": lower["scans_present"],
            "upper_scans_present": upper["scans_present"],
            "shared_scans": len(
                set(p["scan"] for p in tracks[lower["track_id"] - 1])
                & set(p["scan"] for p in tracks[upper["track_id"] - 1])
            ),
        })

    pair_fields = [
        "lower_track_id", "upper_track_id", "lower_median_mz",
        "upper_median_mz", "median_spacing_da", "spacing_residual_mDa",
        "within_3_mDa", "lower_scans_present", "upper_scans_present",
        "shared_scans",
    ]
    with OUT_PAIRS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_fields, delimiter="\t")
        writer.writeheader()
        for row in pair_rows:
            writer.writerow({
                key: f"{value:.9f}" if isinstance(value, float) and "mz" in key
                else f"{value:.6f}" if isinstance(value, float)
                else value
                for key, value in row.items()
            })

    persistent = [
        row for row in track_summaries
        if row["scans_present"] >= 8
    ]

    print(f"Centroid input: {INPUT}")
    print(f"Track details: {OUT_DETAIL}")
    print(f"Track summary: {OUT_SUMMARY}")
    print(f"Adjacent-track spacing: {OUT_PAIRS}")
    print(f"Scans represented: {len(by_scan)}")
    print(f"Centroid observations: {len(rows)}")
    print(f"Tracks generated: {len(tracks)}")
    print(f"Tracks present in at least 8 scans: {len(persistent)}")
    print("\nMost persistent tracks (up to 20):")
    for row in sorted(
        track_summaries,
        key=lambda r: (-r["scans_present"], r["median_mz"])
    )[:20]:
        print(
            f'track {row["track_id"]:>3}: '
            f'{row["scans_present"]}/{row["scans_total"]} scans, '
            f'm/z={row["median_mz"]:.6f}, '
            f'range={row["mz_range_mDa"]:.2f} mDa, '
            f'median intensity={row["median_intensity"]:.1f}'
        )

    print(
        "\nCaution: greedy m/z tracking can merge unrelated peaks or split "
        "a drifting peak; inspect tracks and shared-scan spacing before "
        "interpreting them as isotope groups."
    )


if __name__ == "__main__":
    main()
