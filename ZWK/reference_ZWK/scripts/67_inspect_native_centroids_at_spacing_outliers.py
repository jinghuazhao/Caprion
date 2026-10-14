#!/usr/bin/env python3
"""Inspect native MS1 centroids at isotope-spacing outliers.

Compare trace-guided m/z spacing with the nearest native MS1 centroid
masses in selected scans. This is a diagnostic audit, not a definitive
isotope-group assignment.
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
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
SPACING_FILE = (
    ROOT / "reports/trace_guided_isotope_spacing_v2_same_scan_pairs.tsv"
)

FEATURE_ID = "10275997325121874585"
CHARGE = 2
ISOTOPE_SPACING = 1.00335483507 / CHARGE
MATCH_TOLERANCE_DA = 0.015
LOCAL_HALF_WINDOW_DA = 0.08

OUTPUT_PAIRS = (
    ROOT / "reports/native_centroid_outlier_inspection_v2_pair_diagnostics.tsv"
)
OUTPUT_CENTROIDS = (
    ROOT / "reports/native_centroid_outlier_inspection_v2_local_centroids.tsv"
)

# Scan index, first trace index, second trace index.
TARGET_PAIRS = [
    (5981, 2, 3),
    (5988, 2, 3),
    (5891, 5, 6),
    (5904, 5, 6),
    (5865, 5, 6),
    (5917, 5, 6),
    (5813, 3, 4),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: list[dict],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def fnum(value, default=float("nan")) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def nearest_centroid(
    peaks: list[tuple[float, float]],
    target_mz: float,
    tolerance: float,
) -> tuple[float, float, float] | None:
    """Return (centroid m/z, intensity, signed error in mDa)."""
    candidates = [
        (mz, intensity)
        for mz, intensity in peaks
        if abs(mz - target_mz) <= tolerance
    ]
    if not candidates:
        return None

    mz, intensity = min(candidates, key=lambda item: abs(item[0] - target_mz))
    return mz, intensity, (mz - target_mz) * 1000.0


def main() -> int:
    for path in (MZML, TRACE_FILE, SPACING_FILE):
        if not path.exists():
            print(f"ERROR: required input not found: {path}", file=sys.stderr)
            return 1

    trace_rows = read_tsv(TRACE_FILE)
    spacing_rows = read_tsv(SPACING_FILE)

    feature_rows = [
        row
        for row in trace_rows
        if row.get("feature_id") == FEATURE_ID
    ]
    if not feature_rows:
        print(
            f"ERROR: feature {FEATURE_ID} not found in {TRACE_FILE}",
            file=sys.stderr,
        )
        return 1

    observations: dict[tuple[int, int], dict] = {}
    for row in feature_rows:
        scan_index = int(row["scan_index"])
        trace_index = int(row["trace_index"])
        observations[(scan_index, trace_index)] = row

    # Load native MS1 spectra, retaining the selected scan indices.
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(MZML), exp)

    requested_scans = {scan for scan, _, _ in TARGET_PAIRS}
    native_spectra: dict[int, list[tuple[float, float]]] = {}

    for scan_index, spectrum in enumerate(exp):
        if scan_index not in requested_scans:
            continue
        if spectrum.getMSLevel() != 1:
            continue

        mzs, intensities = spectrum.get_peaks()
        native_spectra[scan_index] = [
            (float(mz), float(intensity))
            for mz, intensity in zip(mzs, intensities)
        ]

    diagnostics = []
    centroid_rows = []

    for scan_index, trace_a, trace_b in TARGET_PAIRS:
        obs_a = observations.get((scan_index, trace_a))
        obs_b = observations.get((scan_index, trace_b))
        peaks = native_spectra.get(scan_index, [])

        row = {
            "scan_index": scan_index,
            "feature_id": FEATURE_ID,
            "charge": CHARGE,
            "trace_a": trace_a,
            "trace_b": trace_b,
            "expected_spacing_da": ISOTOPE_SPACING,
            "native_spectrum_found": bool(peaks),
        }

        for label, obs in (("a", obs_a), ("b", obs_b)):
            trace_mz = fnum(obs.get("observed_mz")) if obs else float("nan")
            trace_intensity = (
                fnum(obs.get("peak_intensity")) if obs else float("nan")
            )
            row[f"trace_{label}_mz"] = trace_mz
            row[f"trace_{label}_intensity"] = trace_intensity

            matched = (
                nearest_centroid(peaks, trace_mz, MATCH_TOLERANCE_DA)
                if peaks and math.isfinite(trace_mz)
                else None
            )

            if matched:
                native_mz, native_intensity, native_error_mda = matched
                row[f"trace_{label}_native_mz"] = native_mz
                row[f"trace_{label}_native_intensity"] = native_intensity
                row[f"trace_{label}_native_error_mda"] = native_error_mda
            else:
                row[f"trace_{label}_native_mz"] = float("nan")
                row[f"trace_{label}_native_intensity"] = float("nan")
                row[f"trace_{label}_native_error_mda"] = float("nan")

            if math.isfinite(trace_mz):
                for mz, intensity in peaks:
                    if abs(mz - trace_mz) <= LOCAL_HALF_WINDOW_DA:
                        centroid_rows.append(
                            {
                                "scan_index": scan_index,
                                "feature_id": FEATURE_ID,
                                "trace_index": (
                                    trace_a if label == "a" else trace_b
                                ),
                                "trace_label": label,
                                "trace_target_mz": trace_mz,
                                "native_centroid_mz": mz,
                                "native_centroid_intensity": intensity,
                                "delta_from_trace_mda": (
                                    mz - trace_mz
                                ) * 1000.0,
                            }
                        )

        trace_a_mz = row["trace_a_mz"]
        trace_b_mz = row["trace_b_mz"]
        native_a_mz = row["trace_a_native_mz"]
        native_b_mz = row["trace_b_native_mz"]

        if math.isfinite(trace_a_mz) and math.isfinite(trace_b_mz):
            trace_spacing = trace_b_mz - trace_a_mz
            row["trace_spacing_da"] = trace_spacing
            row["trace_spacing_residual_mda"] = (
                trace_spacing - ISOTOPE_SPACING
            ) * 1000.0
        else:
            row["trace_spacing_da"] = float("nan")
            row["trace_spacing_residual_mda"] = float("nan")

        if math.isfinite(native_a_mz) and math.isfinite(native_b_mz):
            native_spacing = native_b_mz - native_a_mz
            row["native_spacing_da"] = native_spacing
            row["native_spacing_residual_mda"] = (
                native_spacing - ISOTOPE_SPACING
            ) * 1000.0
        else:
            row["native_spacing_da"] = float("nan")
            row["native_spacing_residual_mda"] = float("nan")

        diagnostics.append(row)

    pair_fields = [
        "scan_index",
        "feature_id",
        "charge",
        "trace_a",
        "trace_b",
        "native_spectrum_found",
        "expected_spacing_da",
        "trace_a_mz",
        "trace_a_intensity",
        "trace_a_native_mz",
        "trace_a_native_intensity",
        "trace_a_native_error_mda",
        "trace_b_mz",
        "trace_b_intensity",
        "trace_b_native_mz",
        "trace_b_native_intensity",
        "trace_b_native_error_mda",
        "trace_spacing_da",
        "trace_spacing_residual_mda",
        "native_spacing_da",
        "native_spacing_residual_mda",
    ]

    centroid_fields = [
        "scan_index",
        "feature_id",
        "trace_index",
        "trace_label",
        "trace_target_mz",
        "native_centroid_mz",
        "native_centroid_intensity",
        "delta_from_trace_mda",
    ]

    write_tsv(OUTPUT_PAIRS, diagnostics, pair_fields)
    write_tsv(OUTPUT_CENTROIDS, centroid_rows, centroid_fields)

    print(f"Feature: {FEATURE_ID}")
    print(f"Charge: {CHARGE}")
    print(f"Expected isotope spacing: {ISOTOPE_SPACING:.9f} Da")
    print(f"Target pairs: {len(TARGET_PAIRS)}")
    print(f"Native MS1 spectra found: {len(native_spectra)}")
    print(f"Pair diagnostics: {OUTPUT_PAIRS}")
    print(f"Local centroids: {OUTPUT_CENTROIDS}")
    print()
    print(
        "scan\tpair\ttrace_spacing\ttrace_residual_mDa\t"
        "native_spacing\tnative_residual_mDa\t"
        "trace_a_native_error_mDa\ttrace_b_native_error_mDa"
    )

    # Fix: precompute errors before formatting, avoiding nested f-strings
    # with conflicting quote characters.
    for row in diagnostics:
        trace_a = row["trace_a"]
        trace_b = row["trace_b"]
        error_a = fnum(row.get("trace_a_native_error_mda"))
        error_b = fnum(row.get("trace_b_native_error_mda"))

        print(
            f"{row['scan_index']}\t{trace_a}-{trace_b}\t"
            f"{fnum(row.get('trace_spacing_da')):.6f}\t"
            f"{fnum(row.get('trace_spacing_residual_mda')):.3f}\t"
            f"{fnum(row.get('native_spacing_da')):.6f}\t"
            f"{fnum(row.get('native_spacing_residual_mda')):.3f}\t"
            f"{error_a:.3f}\t{error_b:.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
