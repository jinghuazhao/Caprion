#!/usr/bin/env python3
"""Inspect near-theoretical native centroid pairs not matching adjacent trace labels."""

import csv
from pathlib import Path

REPORTS = Path("reports")
PAIR_FILE = REPORTS / "feature10275997325121874585_alternative_pairs_v1.tsv"
CENTROID_FILE = REPORTS / "feature10275997325121874585_all_trace_scan_centroids_v1.tsv"

def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

pairs = read_tsv(PAIR_FILE)
centroids = read_tsv(CENTROID_FILE)

exceptions = [
    row for row in pairs
    if row["within_3_mDa"].strip().lower() in {"true", "1", "yes"}
    and row["trace_adjacent_pair"].strip().lower() not in {"true", "1", "yes"}
]

print(f"Non-adjacent-label pairs within ±3 mDa: {len(exceptions)}")

for number, pair in enumerate(exceptions, 1):
    scan = pair["scan_index_zero_based"]
    low = float(pair["lower_centroid_mz"])
    high = float(pair["upper_centroid_mz"])

    print(f"\nException {number}")
    for key in (
        "scan_index_zero_based", "native_scan_number",
        "retention_time_seconds", "lower_centroid_mz",
        "upper_centroid_mz", "spacing_da",
        "spacing_residual_mDa", "lower_trace_labels",
        "upper_trace_labels",
    ):
        print(f"  {key}: {pair[key]}")

    nearby = [
        row for row in centroids
        if row["scan_index_zero_based"] == scan
        and low - 0.03 <= float(row["centroid_mz"]) <= high + 0.03
    ]
    print("  Nearby native centroids:")
    for row in sorted(nearby, key=lambda r: float(r["centroid_mz"])):
        print(
            f'    m/z={float(row["centroid_mz"]):.9f} '
            f'intensity={float(row["centroid_intensity"]):.1f} '
            f'peak_index={row["peak_index"]}'
        )
