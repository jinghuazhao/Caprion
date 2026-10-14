#!/usr/bin/env python3
"""Compare candidate isotope envelopes using same-scan spacing and intensity."""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

REPORTS = Path("reports")
INPUT = REPORTS / "feature10275997325121874585_centroid_tracks_v1.tsv"
OUT = REPORTS / "feature10275997325121874585_envelope_candidate_scores_v1.tsv"
DETAIL = REPORTS / "feature10275997325121874585_envelope_candidate_scan_detail_v1.tsv"

EXPECTED_SPACING = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0

CANDIDATES = {
    "tracks_2_5": [2, 3, 4, 5],
    "tracks_2_6": [2, 3, 4, 5, 6],
    "tracks_2_7": [2, 3, 4, 5, 6, 7],
}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    by_scan = defaultdict(dict)
    for row in read_tsv(INPUT):
        track = int(row["track_id"])
        scan = int(row["scan_index_zero_based"])
        by_scan[scan][track] = {
            "mz": float(row["centroid_mz"]),
            "intensity": float(row["centroid_intensity"]),
            "native_scan": row["native_scan_number"],
            "rt": float(row["retention_time_seconds"]),
        }

    detail_rows = []
    summary_rows = []

    for candidate_name, track_ids in CANDIDATES.items():
        eligible = []
        for scan in sorted(by_scan):
            peaks = by_scan[scan]
            if not all(track in peaks for track in track_ids):
                continue

            spacings = [
                peaks[b]["mz"] - peaks[a]["mz"]
                for a, b in zip(track_ids, track_ids[1:])
            ]
            residuals = [
                (spacing - EXPECTED_SPACING) * 1000.0
                for spacing in spacings
            ]
            abs_residuals = [abs(value) for value in residuals]
            all_within = all(value <= TOLERANCE_MDA for value in abs_residuals)

            intensities = [peaks[track]["intensity"] for track in track_ids]
            monotonic_decrease = all(
                a >= b for a, b in zip(intensities, intensities[1:])
            )

            row = {
                "candidate": candidate_name,
                "track_sequence": "-".join(map(str, track_ids)),
                "scan_index_zero_based": scan,
                "native_scan_number": peaks[track_ids[0]]["native_scan"],
                "retention_time_seconds": peaks[track_ids[0]]["rt"],
                "n_peaks": len(track_ids),
                "median_abs_spacing_residual_mDa": median(abs_residuals),
                "max_abs_spacing_residual_mDa": max(abs_residuals),
                "all_adjacent_spacings_within_3_mDa": all_within,
                "intensity_monotonically_decreasing": monotonic_decrease,
                "intensity_sequence": ";".join(f"{x:.1f}" for x in intensities),
            }
            eligible.append(row)
            detail_rows.append(row)

        summary_rows.append({
            "candidate": candidate_name,
            "track_sequence": "-".join(map(str, track_ids)),
            "n_peaks": len(track_ids),
            "scans_with_all_peaks": len(eligible),
            "scans_all_adjacent_spacings_within_3_mDa": sum(
                row["all_adjacent_spacings_within_3_mDa"] for row in eligible
            ),
            "fraction_scans_all_spacings_within_3_mDa": (
                sum(row["all_adjacent_spacings_within_3_mDa"] for row in eligible)
                / len(eligible) if eligible else None
            ),
            "median_scan_median_abs_residual_mDa": (
                median(row["median_abs_spacing_residual_mDa"] for row in eligible)
                if eligible else None
            ),
            "scans_intensity_monotonic_decrease": sum(
                row["intensity_monotonically_decreasing"] for row in eligible
            ),
        })

    def write_tsv(path, rows):
        if not rows:
            raise RuntimeError(f"No rows to write to {path}")
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(OUT, summary_rows)
    write_tsv(DETAIL, detail_rows)

    print(f"Candidate summary: {OUT}")
    print(f"Per-scan detail: {DETAIL}")
    print(f"Theoretical spacing: {EXPECTED_SPACING:.12f} Da")
    print()
    for row in summary_rows:
        frac = row["fraction_scans_all_spacings_within_3_mDa"]
        med = row["median_scan_median_abs_residual_mDa"]
        print(
            f'{row["candidate"]}: '
            f'complete_scans={row["scans_with_all_peaks"]}, '
            f'all_spacings_within_3_mDa='
            f'{row["scans_all_adjacent_spacings_within_3_mDa"]}/'
            f'{row["scans_with_all_peaks"]}, '
            f'fraction={frac if frac is not None else float("nan"):.3f}, '
            f'median_scan_median_abs_residual='
            f'{med if med is not None else float("nan"):.3f} mDa, '
            f'monotonic_intensity='
            f'{row["scans_intensity_monotonic_decrease"]}/'
            f'{row["scans_with_all_peaks"]}'
        )

    print(
        "\nCaution: these are candidate-envelope diagnostics, not a calibrated "
        "probability that the peaks belong to one precursor."
    )


if __name__ == "__main__":
    main()
