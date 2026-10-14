#!/usr/bin/env python3
"""Audit OpenMS trace observations against native mzML MS1 centroids.

For each observation in the selected feature's existing trace table:
  * locate the native spectrum by scan_index;
  * verify the spectrum is MS1;
  * find the nearest native centroid to observed_mz;
  * record signed/absolute m/z error and intensity diagnostics;
  * summarise match coverage by trace.

This is a correspondence audit, not isotope-group reconstruction.
No original XLSX isotope-group assignments are read.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from pyteomics import mzml


PROJECT = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
DEFAULT_MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
DEFAULT_TRACE = PROJECT / "openms_pilot/isotope_trace_peaks.tsv"
DEFAULT_FEATURE = "10275997325121874585"


def as_float(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def percentile(values, p):
    return float(np.percentile(values, p)) if values else math.nan


def median(values):
    return float(statistics.median(values)) if values else math.nan


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t",
            extrasaction="ignore", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_trace_observations(trace_path, feature_id, trace_indices):
    observations = []
    required = {
        "feature_id", "trace_index", "scan_index", "rt_sec",
        "observed_mz", "peak_intensity",
    }

    with trace_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Trace table is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            if str(row["feature_id"]) != str(feature_id):
                continue

            try:
                trace_index = int(row["trace_index"])
                scan_index = int(row["scan_index"])
            except (TypeError, ValueError):
                continue

            if trace_indices and trace_index not in trace_indices:
                continue

            observed_mz = as_float(row["observed_mz"])
            rt_sec = as_float(row["rt_sec"])
            trace_intensity = as_float(row["peak_intensity"])

            if not math.isfinite(observed_mz):
                continue

            observations.append({
                "feature_id": str(feature_id),
                "trace_index": trace_index,
                "scan_index": scan_index,
                "trace_rt_sec": rt_sec,
                "trace_observed_mz": observed_mz,
                "trace_peak_intensity": trace_intensity,
                "trace_feature_mz": as_float(row.get("feature_mz")),
                "trace_feature_charge": row.get("feature_charge", ""),
            })

    observations.sort(key=lambda x: (x["scan_index"], x["trace_index"]))
    return observations


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--mzml", type=Path, default=DEFAULT_MZML)
    parser.add_argument("--trace-table", type=Path, default=DEFAULT_TRACE)
    parser.add_argument("--feature-id", default=DEFAULT_FEATURE)
    parser.add_argument(
        "--traces", default="1,2,3,4,5,6",
        help="Comma-separated trace indices; default: 1,2,3,4,5,6",
    )
    parser.add_argument(
        "--tolerance-da", type=float, default=0.015,
        help="Maximum absolute m/z difference for a match; default 0.015 Da",
    )
    parser.add_argument(
        "--out-prefix",
        default="trace_native_centroid_audit_v2",
        help="Basename for reports under reports/",
    )
    args = parser.parse_args()

    if args.tolerance_da <= 0:
        parser.error("--tolerance-da must be positive")
    if not args.mzml.is_file():
        parser.error(f"mzML file not found: {args.mzml}")
    if not args.trace_table.is_file():
        parser.error(f"Trace table not found: {args.trace_table}")

    trace_indices = {
        int(item.strip()) for item in args.traces.split(",") if item.strip()
    }

    observations = load_trace_observations(
        args.trace_table, args.feature_id, trace_indices
    )
    if not observations:
        raise SystemExit(
            "No observations found. Check --feature-id, --trace-table "
            "and --traces."
        )

    # Read only spectra needed by the trace table. The spectrum list index is
    # the scan_index used by the validated OpenMS trace-table mapping.
    needed_indices = {x["scan_index"] for x in observations}
    max_index = max(needed_indices)
    spectra_by_index = {}

    print(f"mzML: {args.mzml}")
    print(f"Trace table: {args.trace_table}")
    print(f"Feature ID: {args.feature_id}")
    print(f"Trace indices: {sorted(trace_indices)}")
    print(f"Trace observations: {len(observations):,}")
    print(f"Distinct scan indices: {len(needed_indices):,}")
    print(f"Centroid match tolerance: ±{args.tolerance_da:.6f} Da")
    print("Reading native spectra ...", flush=True)

    with mzml.MzML(str(args.mzml), use_index=False) as reader:
        for spectrum_index, spectrum in enumerate(reader):
            if spectrum_index > max_index:
                break
            if spectrum_index not in needed_indices:
                continue

            ms_level = int(spectrum.get("ms level", 0))
            rt = as_float(
                spectrum.get("scanList", {})
                .get("scan", [{}])[0]
                .get("scan start time")
            )
            if math.isfinite(rt):
                # mzML scan start time is normally in minutes unless converted
                # by pyteomics; normalise to seconds.
                unit = (
                    spectrum.get("scanList", {})
                    .get("scan", [{}])[0]
                    .get("scan start time unitName", "")
                )
                if "minute" in str(unit).lower():
                    rt *= 60.0

            mz_array = np.asarray(spectrum.get("m/z array", []), dtype=float)
            intensity_array = np.asarray(
                spectrum.get("intensity array", []), dtype=float
            )

            if len(mz_array) != len(intensity_array):
                raise ValueError(
                    f"m/z and intensity array lengths differ at spectrum "
                    f"index {spectrum_index}"
                )

            if len(mz_array):
                order = np.argsort(mz_array)
                mz_array = mz_array[order]
                intensity_array = intensity_array[order]

            spectra_by_index[spectrum_index] = {
                "ms_level": ms_level,
                "native_rt_sec": rt,
                "mz": mz_array,
                "intensity": intensity_array,
                "native_peak_count": len(mz_array),
                "spectrum_id": spectrum.get("id", ""),
            }

    missing_spectra = sorted(needed_indices - set(spectra_by_index))
    print(f"Native spectra found: {len(spectra_by_index):,}")
    print(f"Requested indices absent from mzML: {len(missing_spectra):,}")

    detailed_rows = []
    per_trace = defaultdict(list)

    for obs in observations:
        idx = obs["scan_index"]
        native = spectra_by_index.get(idx)

        result = dict(obs)
        result.update({
            "native_spectrum_found": bool(native),
            "native_ms_level": native["ms_level"] if native else "",
            "native_spectrum_id": native["spectrum_id"] if native else "",
            "native_rt_sec": native["native_rt_sec"] if native else math.nan,
            "native_peak_count": native["native_peak_count"] if native else 0,
            "nearest_native_mz": math.nan,
            "nearest_native_intensity": math.nan,
            "signed_error_mDa": math.nan,
            "abs_error_mDa": math.nan,
            "matched_within_tolerance": False,
            "trace_to_native_intensity_ratio": math.nan,
            "log10_intensity_ratio": math.nan,
            "status": "",
        })

        if native is None:
            result["status"] = "spectrum_index_not_found"
        elif native["ms_level"] != 1:
            result["status"] = "spectrum_is_not_MS1"
        elif native["native_peak_count"] == 0:
            result["status"] = "empty_native_spectrum"
        else:
            mzs = native["mz"]
            intensities = native["intensity"]
            target = obs["trace_observed_mz"]

            insertion = int(np.searchsorted(mzs, target))
            candidates = []
            if insertion > 0:
                candidates.append(insertion - 1)
            if insertion < len(mzs):
                candidates.append(insertion)

            nearest_i = min(candidates, key=lambda i: abs(mzs[i] - target))
            nearest_mz = float(mzs[nearest_i])
            nearest_intensity = float(intensities[nearest_i])
            signed_error_mDa = (nearest_mz - target) * 1000.0
            abs_error_mDa = abs(signed_error_mDa)
            matched = abs(nearest_mz - target) <= args.tolerance_da

            result.update({
                "nearest_native_mz": nearest_mz,
                "nearest_native_intensity": nearest_intensity,
                "signed_error_mDa": signed_error_mDa,
                "abs_error_mDa": abs_error_mDa,
                "matched_within_tolerance": matched,
                "status": "matched" if matched else "nearest_peak_outside_tolerance",
            })

            trace_intensity = obs["trace_peak_intensity"]
            if (
                matched and trace_intensity > 0 and nearest_intensity > 0
            ):
                ratio = trace_intensity / nearest_intensity
                result["trace_to_native_intensity_ratio"] = ratio
                result["log10_intensity_ratio"] = math.log10(ratio)

        detailed_rows.append(result)
        per_trace[obs["trace_index"]].append(result)

    report_dir = PROJECT / "reports"
    prefix = args.out_prefix
    detailed_path = report_dir / f"{prefix}_observations.tsv"
    summary_path = report_dir / f"{prefix}_summary.tsv"
    missing_path = report_dir / f"{prefix}_unmatched.tsv"
    absent_path = report_dir / f"{prefix}_missing_spectrum_indices.tsv"

    detail_fields = [
        "feature_id", "trace_index", "scan_index",
        "trace_rt_sec", "native_rt_sec", "trace_observed_mz",
        "nearest_native_mz", "signed_error_mDa", "abs_error_mDa",
        "trace_peak_intensity", "nearest_native_intensity",
        "trace_to_native_intensity_ratio", "log10_intensity_ratio",
        "native_spectrum_found", "native_ms_level", "native_spectrum_id",
        "native_peak_count", "matched_within_tolerance", "status",
        "trace_feature_mz", "trace_feature_charge",
    ]
    write_tsv(detailed_path, detailed_rows, detail_fields)

    summary_rows = []
    for trace_index in sorted(per_trace):
        rows = per_trace[trace_index]
        found = [r for r in rows if r["native_spectrum_found"]]
        ms1_rows = [r for r in found if r["native_ms_level"] == 1]
        matched = [r for r in ms1_rows if r["matched_within_tolerance"]]
        errors = [r["abs_error_mDa"] for r in matched]
        signed = [r["signed_error_mDa"] for r in matched]
        ratios = [
            r["trace_to_native_intensity_ratio"] for r in matched
            if math.isfinite(r["trace_to_native_intensity_ratio"])
        ]

        summary_rows.append({
            "feature_id": args.feature_id,
            "trace_index": trace_index,
            "observations": len(rows),
            "native_spectrum_found": len(found),
            "native_ms1_spectra": len(ms1_rows),
            "matched_within_tolerance": len(matched),
            "match_fraction_of_observations": (
                len(matched) / len(rows) if rows else math.nan
            ),
            "match_fraction_of_native_MS1_observations": (
                len(matched) / len(ms1_rows) if ms1_rows else math.nan
            ),
            "median_abs_error_mDa": median(errors),
            "p95_abs_error_mDa": percentile(errors, 95),
            "max_abs_error_mDa": max(errors) if errors else math.nan,
            "median_signed_error_mDa": median(signed),
            "median_trace_to_native_intensity_ratio": median(ratios),
            "matched_intensity_ratio_count": len(ratios),
        })

    summary_fields = [
        "feature_id", "trace_index", "observations",
        "native_spectrum_found", "native_ms1_spectra",
        "matched_within_tolerance", "match_fraction_of_observations",
        "match_fraction_of_native_MS1_observations",
        "median_abs_error_mDa", "p95_abs_error_mDa",
        "max_abs_error_mDa", "median_signed_error_mDa",
        "median_trace_to_native_intensity_ratio",
        "matched_intensity_ratio_count",
    ]
    write_tsv(summary_path, summary_rows, summary_fields)

    unmatched_rows = [
        row for row in detailed_rows
        if not row["matched_within_tolerance"]
    ]
    write_tsv(missing_path, unmatched_rows, detail_fields)

    write_tsv(
        absent_path,
        [{"scan_index": i} for i in missing_spectra],
        ["scan_index"],
    )

    print("\nPer-trace correspondence summary")
    print(
        "trace observations native_MS1 matched match_fraction "
        "median_abs_error_mDa p95_abs_error_mDa"
    )
    for row in summary_rows:
        median_error = row["median_abs_error_mDa"]
        p95_error = row["p95_abs_error_mDa"]
        print(
            f"{row['trace_index']:>5} "
            f"{row['observations']:>12} "
            f"{row['native_ms1_spectra']:>9} "
            f"{row['matched_within_tolerance']:>7} "
            f"{row['match_fraction_of_observations']:.3f} "
            f"{median_error:.3f} "
            f"{p95_error:.3f}"
        )

    print("\nOutputs:")
    print(f"  {detailed_path}")
    print(f"  {summary_path}")
    print(f"  {missing_path}")
    print(f"  {absent_path}")

    print("\nInterpretation cautions:")
    print("- A match means an observed trace m/z has a nearby native centroid.")
    print("- It does not prove that the centroid is an isotope of this feature.")
    print("- Intensity ratios are diagnostic only; processing/scaling may differ.")
    print("- Do not infer monoisotopic identity or merge isotope groups from this audit.")


if __name__ == "__main__":
    main()
