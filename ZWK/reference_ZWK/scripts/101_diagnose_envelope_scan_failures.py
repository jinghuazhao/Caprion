#!/usr/bin/env python3
"""Identify adjacent isotope-spacing failures in candidate envelopes."""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

REPORTS = Path("reports")
INPUT = REPORTS / "feature10275997325121874585_centroid_tracks_v1.tsv"
OUT_SCAN = REPORTS / "feature10275997325121874585_envelope_spacing_failures_v1.tsv"
OUT_PAIR = REPORTS / "feature10275997325121874585_envelope_pair_diagnostics_v1.tsv"

EXPECTED = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0

CANDIDATES = {
    "tracks_2_5": [2, 3, 4, 5],
    "tracks_2_6": [2, 3, 4, 5, 6],
    "tracks_2_7": [2, 3, 4, 5, 6, 7],
}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    by_scan = defaultdict(dict)

    for row in read_tsv(INPUT):
        scan = int(row["scan_index_zero_based"])
        track = int(row["track_id"])
        by_scan[scan][track] = {
            "mz": float(row["centroid_mz"]),
            "intensity": float(row["centroid_intensity"]),
            "native_scan": row["native_scan_number"],
            "rt": float(row["retention_time_seconds"]),
        }

    scan_rows = []
    pair_observations = defaultdict(list)

    for candidate, track_ids in CANDIDATES.items():
        for scan in sorted(by_scan):
            peaks = by_scan[scan]
            if not all(track in peaks for track in track_ids):
                continue

            residuals = []
            failed_pairs = []
            pair_text = []

            for lower, upper in zip(track_ids, track_ids[1:]):
                low = peaks[lower]
                high = peaks[upper]
                spacing = high["mz"] - low["mz"]
                residual = (spacing - EXPECTED) * 1000.0
                passed = abs(residual) <= TOLERANCE_MDA

                residuals.append(residual)
                pair_text.append(
                    f"{lower}-{upper}:{residual:+.3f}mDa"
                )
                if not passed:
                    failed_pairs.append(f"{lower}-{upper}")

                pair_observations[(candidate, lower, upper)].append({
                    "scan": scan,
                    "residual": residual,
                    "spacing": spacing,
                    "passed": passed,
                })

            scan_rows.append({
                "candidate": candidate,
                "track_sequence": "-".join(map(str, track_ids)),
                "scan_index_zero_based": scan,
                "native_scan_number": peaks[track_ids[0]]["native_scan"],
                "retention_time_seconds": peaks[track_ids[0]]["rt"],
                "adjacent_pair_residuals_mDa": ";".join(pair_text),
                "failed_adjacent_pairs": ",".join(failed_pairs),
                "all_pairs_within_3_mDa": not failed_pairs,
                "median_absolute_spacing_residual_mDa": median(
                    abs(value) for value in residuals
                ),
            })

    pair_rows = []
    for (candidate, lower, upper), observations in sorted(pair_observations.items()):
        residuals = [obs["residual"] for obs in observations]
        passed = sum(obs["passed"] for obs in observations)

        pair_rows.append({
            "candidate": candidate,
            "lower_track": lower,
            "upper_track": upper,
            "n_shared_scans": len(observations),
            "median_spacing_da": median(
                obs["spacing"] for obs in observations
            ),
            "median_residual_mDa": median(residuals),
            "median_absolute_residual_mDa": median(
                abs(value) for value in residuals
            ),
            "max_absolute_residual_mDa": max(
                abs(value) for value in residuals
            ),
            "within_3_mDa": passed,
            "outside_3_mDa": len(observations) - passed,
            "fraction_within_3_mDa": passed / len(observations),
        })

    write_tsv(OUT_SCAN, scan_rows)
    write_tsv(OUT_PAIR, pair_rows)

    print(f"Per-scan diagnosis: {OUT_SCAN}")
    print(f"Pair diagnostics: {OUT_PAIR}")
    print(f"Per-scan candidate rows: {len(scan_rows)}")
    print("\nFailure counts by candidate:")
    for candidate in CANDIDATES:
        subset = [row for row in scan_rows if row["candidate"] == candidate]
        failures = sum(not row["all_pairs_within_3_mDa"] for row in subset)
        print(
            f"{candidate}: complete_scans={len(subset)}, "
            f"scans_with_any_spacing_failure={failures}, "
            f"all_pass={len(subset) - failures}"
        )

    print("\nPair diagnostics:")
    for row in pair_rows:
        print(
            f'{row["candidate"]} pair {row["lower_track"]}-{row["upper_track"]}: '
            f'n={row["n_shared_scans"]}, '
            f'median_abs_residual={row["median_absolute_residual_mDa"]:.3f} mDa, '
            f'within_3_mDa={row["within_3_mDa"]}/{row["n_shared_scans"]}'
        )


if __name__ == "__main__":
    main()
