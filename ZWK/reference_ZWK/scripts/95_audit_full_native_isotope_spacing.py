#!/usr/bin/env python3
"""Audit adjacent isotope-trace spacing using native mzML centroid matches."""

import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

ROOT = Path.cwd()
REPORTS = ROOT / "reports"

INPUT = REPORTS / "feature10275997325121874585_trace_centroid_matches_v1.tsv"
OUT_ROWS = REPORTS / "feature10275997325121874585_native_spacing_full_v1.tsv"
OUT_SUMMARY = REPORTS / "feature10275997325121874585_native_spacing_summary_full_v1.tsv"

THEORETICAL_SPACING = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    rows = read_tsv(INPUT)
    if not rows:
        raise RuntimeError(f"No rows in {INPUT}")

    by_scan = defaultdict(dict)
    for row in rows:
        scan = int(row["scan_index_zero_based"])
        trace_index = numeric(row["trace_index"])
        if trace_index is None:
            continue
        trace_index = int(trace_index)

        observed_mz = numeric(row["trace_mz"])
        centroid_mz = numeric(row["nearest_centroid_mz"])
        intensity = numeric(row["nearest_centroid_intensity"])
        matched = row["within_0_015_Da"].strip().lower() in {
            "true", "1", "yes"
        }

        by_scan[scan][trace_index] = {
            "native_scan": row["native_scan_number"],
            "rt": numeric(row["spectrum_rt_seconds"]),
            "trace_mz": observed_mz,
            "centroid_mz": centroid_mz,
            "intensity": intensity,
            "matched": matched,
        }

    pair_rows = []
    for scan in sorted(by_scan):
        traces = by_scan[scan]
        for lower_index in sorted(traces):
            upper_index = lower_index + 1
            if upper_index not in traces:
                continue

            lower = traces[lower_index]
            upper = traces[upper_index]

            spacing = (
                upper["centroid_mz"] - lower["centroid_mz"]
                if lower["centroid_mz"] is not None
                and upper["centroid_mz"] is not None
                else None
            )
            residual_mda = (
                (spacing - THEORETICAL_SPACING) * 1000.0
                if spacing is not None
                else None
            )
            both_matched = lower["matched"] and upper["matched"]

            pair_rows.append({
                "scan_index_zero_based": scan,
                "native_scan_number": lower["native_scan"],
                "retention_time_seconds": lower["rt"],
                "lower_trace_index": lower_index,
                "upper_trace_index": upper_index,
                "lower_trace_mz": lower["trace_mz"],
                "upper_trace_mz": upper["trace_mz"],
                "lower_native_centroid_mz": lower["centroid_mz"],
                "upper_native_centroid_mz": upper["centroid_mz"],
                "lower_centroid_intensity": lower["intensity"],
                "upper_centroid_intensity": upper["intensity"],
                "both_trace_peaks_matched": both_matched,
                "native_centroid_spacing_da": spacing,
                "theoretical_spacing_da": THEORETICAL_SPACING,
                "spacing_residual_mDa": residual_mda,
                "within_3_mDa": (
                    abs(residual_mda) <= TOLERANCE_MDA
                    if residual_mda is not None
                    else False
                ),
            })

    pair_fields = [
        "scan_index_zero_based", "native_scan_number",
        "retention_time_seconds", "lower_trace_index", "upper_trace_index",
        "lower_trace_mz", "upper_trace_mz",
        "lower_native_centroid_mz", "upper_native_centroid_mz",
        "lower_centroid_intensity", "upper_centroid_intensity",
        "both_trace_peaks_matched", "native_centroid_spacing_da",
        "theoretical_spacing_da", "spacing_residual_mDa", "within_3_mDa",
    ]

    with OUT_ROWS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_fields, delimiter="\t")
        writer.writeheader()
        for row in pair_rows:
            writer.writerow({
                key: (
                    f"{value:.9f}"
                    if isinstance(value, float) and math.isfinite(value)
                    and ("mz" in key or "spacing" in key)
                    else f"{value:.6f}"
                    if isinstance(value, float) and math.isfinite(value)
                    else value
                )
                for key, value in row.items()
            })

    summary_rows = []
    pair_labels = sorted({
        (row["lower_trace_index"], row["upper_trace_index"])
        for row in pair_rows
    })

    for lower_index, upper_index in pair_labels:
        subset = [
            row for row in pair_rows
            if row["lower_trace_index"] == lower_index
            and row["upper_trace_index"] == upper_index
        ]
        matched = [
            row for row in subset
            if row["both_trace_peaks_matched"]
            and row["spacing_residual_mDa"] is not None
        ]
        residuals = [row["spacing_residual_mDa"] for row in matched]
        within = sum(abs(value) <= TOLERANCE_MDA for value in residuals)

        summary_rows.append({
            "lower_trace_index": lower_index,
            "upper_trace_index": upper_index,
            "pair_rows": len(subset),
            "both_peaks_matched": len(matched),
            "median_native_spacing_da": (
                median(row["native_centroid_spacing_da"] for row in matched)
                if matched else None
            ),
            "median_spacing_residual_mDa": median(residuals) if residuals else None,
            "median_absolute_residual_mDa": (
                median(abs(value) for value in residuals) if residuals else None
            ),
            "max_absolute_residual_mDa": (
                max(abs(value) for value in residuals) if residuals else None
            ),
            "within_3_mDa": within,
            "within_3_mDa_fraction": within / len(matched) if matched else None,
        })

    summary_fields = [
        "lower_trace_index", "upper_trace_index", "pair_rows",
        "both_peaks_matched", "median_native_spacing_da",
        "median_spacing_residual_mDa", "median_absolute_residual_mDa",
        "max_absolute_residual_mDa", "within_3_mDa",
        "within_3_mDa_fraction",
    ]
    with OUT_SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=summary_fields, delimiter="\t"
        )
        writer.writeheader()
        for row in summary_rows:
            writer.writerow({
                key: (
                    f"{value:.9f}"
                    if isinstance(value, float) and math.isfinite(value)
                    and "spacing_da" in key
                    else f"{value:.6f}"
                    if isinstance(value, float) and math.isfinite(value)
                    else value
                )
                for key, value in row.items()
            })

    print(f"Input trace-centroid matches: {INPUT}")
    print(f"Full adjacent-pair audit: {OUT_ROWS}")
    print(f"Pair summary: {OUT_SUMMARY}")
    print(f"Adjacent-pair observations: {len(pair_rows)}")
    print(f"Theoretical spacing: {THEORETICAL_SPACING:.12f} Da")
    print()
    print("Pair summary:")
    for row in summary_rows:
        med = row["median_native_spacing_da"]
        absres = row["median_absolute_residual_mDa"]
        fraction = row["within_3_mDa_fraction"]
        print(
            f'{row["lower_trace_index"]}-{row["upper_trace_index"]}: '
            f'n={row["pair_rows"]}, matched={row["both_peaks_matched"]}, '
            f'median_spacing={med if med is not None else float("nan"):.6f} Da, '
            f'median_abs_residual={absres if absres is not None else float("nan"):.3f} mDa, '
            f'within_3_mDa={row["within_3_mDa"]}/'
            f'{row["both_peaks_matched"]}, '
            f'fraction={fraction if fraction is not None else float("nan"):.3f}'
        )


if __name__ == "__main__":
    main()
