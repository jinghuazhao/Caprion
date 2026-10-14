#!/usr/bin/env python3
"""Expand the native MS1 centroid search at isotope-spacing outliers.

Reads the original mzML directly. It does not use the original XLSX IG
assignments and does not automatically reassign isotope traces.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import pyopenms as oms

ROOT = Path(__file__).resolve().parents[1]

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
PAIR_FILE = (
    ROOT / "reports/native_centroid_outlier_inspection_v2_pair_diagnostics.tsv"
)
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"

OUTPUT_CANDIDATES = (
    ROOT / "reports/expanded_native_isotope_candidates_v2.tsv"
)
OUTPUT_SUMMARY = (
    ROOT / "reports/expanded_native_isotope_candidates_summary_v2.tsv"
)

FEATURE_ID = "10275997325121874585"
ISOTOPE_SPACING = 1.00335483507 / 2.0

# Search a wider region around the predicted upper isotope position.
SEARCH_HALF_WIDTH_DA = 0.15
TRACE_TARGET_TOLERANCE_MDA = 5.0
TOP_N = 15


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


def to_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def main() -> int:
    for path in (MZML, PAIR_FILE, TRACE_FILE):
        if not path.exists():
            print(f"ERROR: missing input: {path}", file=sys.stderr)
            return 1

    pair_rows = read_tsv(PAIR_FILE)
    trace_rows = [
        row for row in read_tsv(TRACE_FILE)
        if row.get("feature_id") == FEATURE_ID
    ]

    if not trace_rows:
        print(f"ERROR: no trace observations for {FEATURE_ID}", file=sys.stderr)
        return 1

    traces_by_scan: dict[int, list[dict]] = {}
    for row in trace_rows:
        scan = to_int(row.get("scan_index"))
        if scan is not None:
            traces_by_scan.setdefault(scan, []).append(row)

    # Build the requested scan/pair list from Script 67's diagnostics.
    requests = []
    for row in pair_rows:
        scan = to_int(row.get("scan_index"))
        trace_a = to_int(row.get("trace_a"))
        trace_b = to_int(row.get("trace_b"))
        mz_a = num(row.get("trace_a_native_mz"))
        mz_b = num(row.get("trace_b_native_mz"))

        if (
            scan is None or trace_a is None or trace_b is None
            or not math.isfinite(mz_a) or not math.isfinite(mz_b)
        ):
            continue

        requests.append({
            "scan_index": scan,
            "trace_a": trace_a,
            "trace_b": trace_b,
            "trace_a_mz": mz_a,
            "trace_b_mz": mz_b,
            "expected_upper_mz": mz_a + ISOTOPE_SPACING,
            "original_spacing_residual_mda": (
                (mz_b - mz_a) - ISOTOPE_SPACING
            ) * 1000.0,
        })

    if not requests:
        print("ERROR: no valid requests in Script 67 output.", file=sys.stderr)
        return 1

    requested_scans = {row["scan_index"] for row in requests}

    # Read the native mzML. The scan index here is the zero-based position
    # in the MSExperiment, matching the indexing used in Script 67.
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(MZML), exp)

    spectra: dict[int, list[tuple[float, float]]] = {}
    for scan_index, spectrum in enumerate(exp):
        if scan_index not in requested_scans or spectrum.getMSLevel() != 1:
            continue
        mzs, intensities = spectrum.get_peaks()
        spectra[scan_index] = [
            (float(mz), float(intensity))
            for mz, intensity in zip(mzs, intensities)
        ]

    candidate_output = []
    summary_output = []

    for request in requests:
        scan = request["scan_index"]
        expected_mz = request["expected_upper_mz"]
        lower_mz = expected_mz - SEARCH_HALF_WIDTH_DA
        upper_mz = expected_mz + SEARCH_HALF_WIDTH_DA
        peaks = spectra.get(scan, [])

        scan_trace_targets = []
        for trace in traces_by_scan.get(scan, []):
            target = num(trace.get("observed_mz"))
            if math.isfinite(target):
                scan_trace_targets.append({
                    "trace_index": to_int(trace.get("trace_index")),
                    "mz": target,
                })

        candidates = []
        for mz, intensity in peaks:
            if not lower_mz <= mz <= upper_mz:
                continue

            signed_error_mda = (mz - expected_mz) * 1000.0
            trace_distances = [
                {
                    "trace_index": item["trace_index"],
                    "distance_mda": abs(mz - item["mz"]) * 1000.0,
                }
                for item in scan_trace_targets
            ]
            nearest_trace = (
                min(trace_distances, key=lambda item: item["distance_mda"])
                if trace_distances else None
            )

            candidates.append({
                **request,
                "candidate_native_mz": mz,
                "candidate_intensity": intensity,
                "candidate_spacing_da": mz - request["trace_a_mz"],
                "candidate_spacing_residual_mda": signed_error_mda,
                "candidate_abs_residual_mda": abs(signed_error_mda),
                "nearest_trace_index": (
                    nearest_trace["trace_index"] if nearest_trace else None
                ),
                "nearest_trace_target_distance_mda": (
                    nearest_trace["distance_mda"] if nearest_trace
                    else float("nan")
                ),
                "near_trace_target": bool(
                    nearest_trace
                    and nearest_trace["distance_mda"] <= TRACE_TARGET_TOLERANCE_MDA
                ),
                "is_original_upper_centroid": (
                    abs(mz - request["trace_b_mz"]) < 1e-7
                ),
            })

        # Rank by isotope-spacing agreement first, then intensity.
        candidates.sort(
            key=lambda row: (
                row["candidate_abs_residual_mda"],
                -row["candidate_intensity"],
            )
        )

        for rank, candidate in enumerate(candidates[:TOP_N], start=1):
            candidate["candidate_rank"] = rank
            candidate_output.append(candidate)

        best = candidates[0] if candidates else None
        summary_output.append({
            "scan_index": scan,
            "trace_pair": f"{request['trace_a']}-{request['trace_b']}",
            "expected_upper_mz": expected_mz,
            "original_upper_mz": request["trace_b_mz"],
            "original_spacing_residual_mda": (
                request["original_spacing_residual_mda"]
            ),
            "search_lower_mz": lower_mz,
            "search_upper_mz": upper_mz,
            "native_spectrum_found": scan in spectra,
            "candidates_in_window": len(candidates),
            "best_candidate_mz": (
                best["candidate_native_mz"] if best else float("nan")
            ),
            "best_candidate_intensity": (
                best["candidate_intensity"] if best else float("nan")
            ),
            "best_candidate_residual_mda": (
                best["candidate_spacing_residual_mda"]
                if best else float("nan")
            ),
            "best_candidate_is_original": (
                best["is_original_upper_centroid"] if best else False
            ),
            "best_candidate_near_trace_target": (
                best["near_trace_target"] if best else False
            ),
        })

    candidate_fields = [
        "scan_index", "trace_a", "trace_b",
        "trace_a_mz", "trace_b_mz", "expected_upper_mz",
        "original_spacing_residual_mda",
        "candidate_rank", "candidate_native_mz", "candidate_intensity",
        "candidate_spacing_da", "candidate_spacing_residual_mda",
        "candidate_abs_residual_mda", "nearest_trace_index",
        "nearest_trace_target_distance_mda", "near_trace_target",
        "is_original_upper_centroid",
    ]

    summary_fields = [
        "scan_index", "trace_pair", "expected_upper_mz",
        "original_upper_mz", "original_spacing_residual_mda",
        "search_lower_mz", "search_upper_mz", "native_spectrum_found",
        "candidates_in_window", "best_candidate_mz",
        "best_candidate_intensity", "best_candidate_residual_mda",
        "best_candidate_is_original", "best_candidate_near_trace_target",
    ]

    write_tsv(OUTPUT_CANDIDATES, candidate_output, candidate_fields)
    write_tsv(OUTPUT_SUMMARY, summary_output, summary_fields)

    print(f"Feature: {FEATURE_ID}")
    print(f"Native MS1 spectra found: {len(spectra)}")
    print(f"Pairs inspected: {len(summary_output)}")
    print(f"Candidate rows exported: {len(candidate_output)}")
    print(f"Candidate table: {OUTPUT_CANDIDATES}")
    print(f"Summary table: {OUTPUT_SUMMARY}")
    print()
    print(
        "scan\tpair\toriginal_residual_mDa\tcandidates\t"
        "best_mz\tbest_residual_mDa\tbest_is_original"
    )

    for row in summary_output:
        print(
            f"{row['scan_index']}\t{row['trace_pair']}\t"
            f"{num(row['original_spacing_residual_mda']):.3f}\t"
            f"{row['candidates_in_window']}\t"
            f"{num(row['best_candidate_mz']):.6f}\t"
            f"{num(row['best_candidate_residual_mda']):.3f}\t"
            f"{row['best_candidate_is_original']}"
        )

    print()
    print(
        "Diagnostic only: mass agreement does not prove isotope identity. "
        "Review candidate intensity, trace overlap and neighbouring scans "
        "before changing any assignments."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
