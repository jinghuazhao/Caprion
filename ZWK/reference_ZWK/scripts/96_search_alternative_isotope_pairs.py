#!/usr/bin/env python3
"""Search native centroid data for alternative charge-2 isotope-spaced pairs."""

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path.cwd()
REPORTS = ROOT / "reports"

CENTROID_FILE = REPORTS / (
    "feature10275997325121874585_all_trace_scan_centroids_v1.tsv"
)
TRACE_MATCH_FILE = REPORTS / (
    "feature10275997325121874585_trace_centroid_matches_v1.tsv"
)
OUT = REPORTS / "feature10275997325121874585_alternative_pairs_v1.tsv"
SUMMARY = REPORTS / "feature10275997325121874585_alternative_pairs_summary_v1.tsv"

THEORETICAL_SPACING = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    centroids = read_tsv(CENTROID_FILE)
    trace_matches = read_tsv(TRACE_MATCH_FILE)

    by_scan = defaultdict(list)
    for row in centroids:
        by_scan[int(row["scan_index_zero_based"])].append({
            "native_scan": row["native_scan_number"],
            "rt": float(row["rt_seconds"]),
            "peak_index": int(row["peak_index"]),
            "mz": float(row["centroid_mz"]),
            "intensity": float(row["centroid_intensity"]),
        })

    trace_mz_by_scan = defaultdict(dict)
    for row in trace_matches:
        idx = row["trace_index"].strip()
        if idx:
            trace_mz_by_scan[int(row["scan_index_zero_based"])][int(idx)] = (
                float(row["trace_mz"])
            )

    out_rows = []
    for scan, peaks in sorted(by_scan.items()):
        peaks.sort(key=lambda p: p["mz"])
        trace_mzs = trace_mz_by_scan.get(scan, {})

        for i, lower in enumerate(peaks):
            for upper in peaks[i + 1:]:
                spacing = upper["mz"] - lower["mz"]

                # Only consider pairs near the theoretical adjacent isotope spacing.
                residual_mda = (spacing - THEORETICAL_SPACING) * 1000.0
                if abs(residual_mda) > 10.0:
                    continue

                lower_trace = [
                    idx for idx, mz in trace_mzs.items()
                    if abs(mz - lower["mz"]) <= 0.015
                ]
                upper_trace = [
                    idx for idx, mz in trace_mzs.items()
                    if abs(mz - upper["mz"]) <= 0.015
                ]

                out_rows.append({
                    "scan_index_zero_based": scan,
                    "native_scan_number": lower["native_scan"],
                    "retention_time_seconds": lower["rt"],
                    "lower_peak_index": lower["peak_index"],
                    "upper_peak_index": upper["peak_index"],
                    "lower_centroid_mz": lower["mz"],
                    "upper_centroid_mz": upper["mz"],
                    "lower_intensity": lower["intensity"],
                    "upper_intensity": upper["intensity"],
                    "spacing_da": spacing,
                    "spacing_residual_mDa": residual_mda,
                    "within_3_mDa": abs(residual_mda) <= TOLERANCE_MDA,
                    "lower_trace_labels": ",".join(map(str, sorted(lower_trace))),
                    "upper_trace_labels": ",".join(map(str, sorted(upper_trace))),
                    "trace_adjacent_pair": (
                        len(lower_trace) == 1
                        and len(upper_trace) == 1
                        and upper_trace[0] == lower_trace[0] + 1
                    ),
                })

    fields = [
        "scan_index_zero_based", "native_scan_number",
        "retention_time_seconds", "lower_peak_index", "upper_peak_index",
        "lower_centroid_mz", "upper_centroid_mz", "lower_intensity",
        "upper_intensity", "spacing_da", "spacing_residual_mDa",
        "within_3_mDa", "lower_trace_labels", "upper_trace_labels",
        "trace_adjacent_pair",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    good = [r for r in out_rows if r["within_3_mDa"]]
    trace_good = [r for r in good if r["trace_adjacent_pair"]]
    summary_rows = [
        {"metric": "unique_scans_with_centroids", "value": len(by_scan)},
        {"metric": "all_pairs_within_10_mDa", "value": len(out_rows)},
        {"metric": "pairs_within_3_mDa", "value": len(good)},
        {"metric": "within_3_mDa_pairs_matching_adjacent_trace_labels", "value": len(trace_good)},
        {"metric": "within_3_mDa_pairs_not_matching_adjacent_trace_labels", "value": len(good) - len(trace_good)},
    ]
    with SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Alternative-pair table: {OUT}")
    print(f"Summary: {SUMMARY}")
    for row in summary_rows:
        print(f'{row["metric"]}: {row["value"]}')
    print("Note: spacing alone does not establish common precursor identity.")


if __name__ == "__main__":
    main()
