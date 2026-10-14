#!/usr/bin/env python3
"""Rank alternative native MS1 centroids for isotope-spacing outliers.

Uses Script 67's exported local-centroid table and pair diagnostics.
This is a diagnostic ranking only; it does not assign isotope groups.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PAIR_FILE = (
    ROOT / "reports/native_centroid_outlier_inspection_v2_pair_diagnostics.tsv"
)
CENTROID_FILE = (
    ROOT / "reports/native_centroid_outlier_inspection_v2_local_centroids.tsv"
)

OUTPUT_CANDIDATES = (
    ROOT / "reports/alternative_native_isotope_centroid_candidates_v2.tsv"
)
OUTPUT_SUMMARY = (
    ROOT / "reports/alternative_native_isotope_centroid_summary_v2.tsv"
)

ISOTOPE_SPACING = 1.00335483507 / 2.0

# Candidate must be within this distance of the predicted isotope position.
MAX_EXPECTED_ERROR_MDA = 20.0

# Identify centroids very close to any trace target in the same scan.
TRACE_TARGET_TOLERANCE_MDA = 5.0

# Report this many best-ranked candidates per pair.
TOP_N = 10


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def num(value, default=float("nan")) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def int_or_none(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def main() -> int:
    for path in (PAIR_FILE, CENTROID_FILE):
        if not path.exists():
            print(
                f"ERROR: missing input {path}. Run Script 67 first.",
                file=sys.stderr,
            )
            return 1

    pair_rows = read_tsv(PAIR_FILE)
    centroid_rows = read_tsv(CENTROID_FILE)

    # Centroids are duplicated in the local-centroid export when they fall
    # within the windows around both trace targets. Deduplicate by scan/mz/
    # intensity so the same physical centroid is not ranked multiple times.
    centroids_by_scan: dict[int, dict[tuple[float, float], dict]] = {}
    for row in centroid_rows:
        scan = int_or_none(row.get("scan_index"))
        mz = num(row.get("native_centroid_mz"))
        intensity = num(row.get("native_centroid_intensity"))
        if scan is None or not math.isfinite(mz):
            continue

        key = (round(mz, 8), round(intensity, 3))
        centroids_by_scan.setdefault(scan, {})[key] = {
            "scan_index": scan,
            "native_centroid_mz": mz,
            "native_centroid_intensity": intensity,
        }

    candidates_out = []
    summary_out = []

    for pair in pair_rows:
        scan = int_or_none(pair.get("scan_index"))
        trace_a = int_or_none(pair.get("trace_a"))
        trace_b = int_or_none(pair.get("trace_b"))

        if scan is None or trace_a is None or trace_b is None:
            continue

        mz_a = num(pair.get("trace_a_native_mz"))
        mz_b = num(pair.get("trace_b_native_mz"))
        if not math.isfinite(mz_a) or not math.isfinite(mz_b):
            continue

        expected_mz_b = mz_a + ISOTOPE_SPACING
        original_spacing = mz_b - mz_a
        original_error_mda = (original_spacing - ISOTOPE_SPACING) * 1000.0

        scan_centroids = list(centroids_by_scan.get(scan, {}).values())

        # Candidate trace-target m/z values from this scan are useful for
        # identifying whether an alternative centroid overlaps another trace.
        trace_targets = []
        for candidate_pair in pair_rows:
            if int_or_none(candidate_pair.get("scan_index")) != scan:
                continue
            for field in ("trace_a_mz", "trace_b_mz"):
                target = num(candidate_pair.get(field))
                if math.isfinite(target):
                    trace_targets.append(target)

        ranked = []
        for centroid in scan_centroids:
            candidate_mz = centroid["native_centroid_mz"]
            candidate_intensity = centroid["native_centroid_intensity"]

            error_mda = (candidate_mz - expected_mz_b) * 1000.0
            abs_error_mda = abs(error_mda)

            if abs_error_mda > MAX_EXPECTED_ERROR_MDA:
                continue

            trace_distances = [
                abs(candidate_mz - target) * 1000.0
                for target in trace_targets
            ]
            nearest_trace_distance = (
                min(trace_distances) if trace_distances else float("nan")
            )
            near_trace_target = (
                math.isfinite(nearest_trace_distance)
                and nearest_trace_distance <= TRACE_TARGET_TOLERANCE_MDA
            )

            # Rank by mass agreement first, intensity second.
            ranked.append(
                {
                    "scan_index": scan,
                    "trace_a": trace_a,
                    "trace_b": trace_b,
                    "trace_a_native_mz": mz_a,
                    "trace_b_native_mz_original": mz_b,
                    "expected_trace_b_mz": expected_mz_b,
                    "original_spacing_da": original_spacing,
                    "original_spacing_residual_mda": original_error_mda,
                    "candidate_native_mz": candidate_mz,
                    "candidate_intensity": candidate_intensity,
                    "candidate_spacing_da": candidate_mz - mz_a,
                    "candidate_spacing_residual_mda": error_mda,
                    "candidate_abs_residual_mda": abs_error_mda,
                    "nearest_trace_target_distance_mda": nearest_trace_distance,
                    "near_existing_trace_target": near_trace_target,
                    "candidate_is_original_upper_centroid": (
                        abs(candidate_mz - mz_b) < 1e-7
                    ),
                }
            )

        ranked.sort(
            key=lambda row: (
                row["candidate_abs_residual_mda"],
                -row["candidate_intensity"],
            )
        )

        for rank, candidate in enumerate(ranked[:TOP_N], start=1):
            candidate["candidate_rank"] = rank
            candidates_out.append(candidate)

        best = ranked[0] if ranked else None
        summary_out.append(
            {
                "scan_index": scan,
                "trace_pair": f"{trace_a}-{trace_b}",
                "expected_upper_mz": expected_mz_b,
                "original_upper_mz": mz_b,
                "original_spacing_residual_mda": original_error_mda,
                "local_centroids_in_export": len(scan_centroids),
                "candidates_within_20mDa": len(ranked),
                "best_candidate_mz": (
                    best["candidate_native_mz"] if best else float("nan")
                ),
                "best_candidate_intensity": (
                    best["candidate_intensity"] if best else float("nan")
                ),
                "best_candidate_residual_mda": (
                    best["candidate_spacing_residual_mda"]
                    if best
                    else float("nan")
                ),
                "best_candidate_is_original": (
                    best["candidate_is_original_upper_centroid"]
                    if best
                    else False
                ),
                "best_candidate_near_trace_target": (
                    best["near_existing_trace_target"] if best else False
                ),
            }
        )

    candidate_fields = [
        "scan_index",
        "trace_a",
        "trace_b",
        "candidate_rank",
        "trace_a_native_mz",
        "trace_b_native_mz_original",
        "expected_trace_b_mz",
        "original_spacing_da",
        "original_spacing_residual_mda",
        "candidate_native_mz",
        "candidate_intensity",
        "candidate_spacing_da",
        "candidate_spacing_residual_mda",
        "candidate_abs_residual_mda",
        "nearest_trace_target_distance_mda",
        "near_existing_trace_target",
        "candidate_is_original_upper_centroid",
    ]

    summary_fields = [
        "scan_index",
        "trace_pair",
        "expected_upper_mz",
        "original_upper_mz",
        "original_spacing_residual_mda",
        "local_centroids_in_export",
        "candidates_within_20mDa",
        "best_candidate_mz",
        "best_candidate_intensity",
        "best_candidate_residual_mda",
        "best_candidate_is_original",
        "best_candidate_near_trace_target",
    ]

    write_tsv(OUTPUT_CANDIDATES, candidates_out, candidate_fields)
    write_tsv(OUTPUT_SUMMARY, summary_out, summary_fields)

    print(f"Outlier pairs inspected: {len(summary_out)}")
    print(f"Candidate rows exported: {len(candidates_out)}")
    print(f"Candidate table: {OUTPUT_CANDIDATES}")
    print(f"Summary table: {OUTPUT_SUMMARY}")
    print()
    print(
        "scan\tpair\toriginal_residual_mDa\t"
        "candidate_count\tbest_mz\tbest_residual_mDa\t"
        "best_is_original\tnear_other_trace"
    )

    for row in summary_out:
        best_mz = num(row["best_candidate_mz"])
        best_error = num(row["best_candidate_residual_mda"])
        print(
            f"{row['scan_index']}\t{row['trace_pair']}\t"
            f"{num(row['original_spacing_residual_mda']):.3f}\t"
            f"{row['candidates_within_20mDa']}\t"
            f"{best_mz:.6f}\t{best_error:.3f}\t"
            f"{row['best_candidate_is_original']}\t"
            f"{row['best_candidate_near_trace_target']}"
        )

    print()
    print(
        "Caution: the candidate list is limited to centroids exported by "
        "Script 67. A missing alternative may simply fall outside its "
        "local-centroid export window."
    )
    print(
        "A candidate with a smaller mass residual is not automatically the "
        "correct isotope; inspect intensity, trace overlap, and neighbouring "
        "scans before changing assignments."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
