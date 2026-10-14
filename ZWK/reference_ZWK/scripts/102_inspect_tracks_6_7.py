#!/usr/bin/env python3
"""Inspect scan-level spacing between centroid tracks 6 and 7."""

import csv
from pathlib import Path

INPUT = Path(
    "reports/feature10275997325121874585_centroid_tracks_v1.tsv"
)
OUTPUT = Path(
    "reports/feature10275997325121874585_tracks_6_7_scan_inspection_v1.tsv"
)

EXPECTED_SPACING = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0


def main():
    by_scan = {}

    with INPUT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            track = int(row["track_id"])
            if track not in (6, 7):
                continue

            scan = int(row["scan_index_zero_based"])
            by_scan.setdefault(scan, {})[track] = {
                "mz": float(row["centroid_mz"]),
                "intensity": float(row["centroid_intensity"]),
                "native_scan": row["native_scan_number"],
                "rt": float(row["retention_time_seconds"]),
            }

    results = []
    for scan in sorted(by_scan):
        peaks = by_scan[scan]
        if 6 not in peaks or 7 not in peaks:
            continue

        p6, p7 = peaks[6], peaks[7]
        spacing = p7["mz"] - p6["mz"]
        residual_mda = (spacing - EXPECTED_SPACING) * 1000.0

        results.append({
            "scan_index_zero_based": scan,
            "native_scan_number": p6["native_scan"],
            "retention_time_seconds": f'{p6["rt"]:.5f}',
            "track6_mz": f'{p6["mz"]:.6f}',
            "track6_intensity": f'{p6["intensity"]:.1f}',
            "track7_mz": f'{p7["mz"]:.6f}',
            "track7_intensity": f'{p7["intensity"]:.1f}',
            "observed_spacing_da": f"{spacing:.9f}",
            "expected_spacing_da": f"{EXPECTED_SPACING:.9f}",
            "residual_mDa": f"{residual_mda:+.3f}",
            "within_3_mDa": abs(residual_mda) <= TOLERANCE_MDA,
            "track7_to_track6_intensity_ratio": (
                f'{p7["intensity"] / p6["intensity"]:.5f}'
                if p6["intensity"] else ""
            ),
        })

    if not results:
        raise SystemExit("No scans containing both tracks 6 and 7 were found.")

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(results[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"Output: {OUTPUT}")
    print(f"Shared scans: {len(results)}")
    print(
        "Scan  NativeScan       RT(s)     Track6_mz   Track7_mz"
        "   Spacing(Da)   Residual(mDa)  Pass  I7/I6"
    )
    for row in results:
        print(
            f'{row["scan_index_zero_based"]:4d} '
            f'{row["native_scan_number"]:>10s} '
            f'{row["retention_time_seconds"]:>11s} '
            f'{row["track6_mz"]:>12s} '
            f'{row["track7_mz"]:>12s} '
            f'{row["observed_spacing_da"]:>13s} '
            f'{row["residual_mDa"]:>13s} '
            f'{str(row["within_3_mDa"]):>5s} '
            f'{row["track7_to_track6_intensity_ratio"]:>7s}'
        )


if __name__ == "__main__":
    main()
