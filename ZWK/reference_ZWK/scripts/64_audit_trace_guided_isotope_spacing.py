#!/usr/bin/env python3
"""Audit same-scan isotope-like spacing between existing OpenMS traces.

This script:
  1. Reads the existing trace observations for one feature.
  2. Matches observations by scan_index.
  3. Computes adjacent-trace m/z spacing only within the same scan.
  4. Compares spacing with the expected isotope spacing for a charge-2 ion.
  5. Reports trace-level and scan-level diagnostics.

It does not assign isotope indices, identify a monoisotopic peak, merge traces,
or use the original XLSX isotope-group assignments.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np


PROJECT = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
DEFAULT_TRACE = PROJECT / "openms_pilot/isotope_trace_peaks.tsv"
DEFAULT_FEATURE = "10275997325121874585"

# Difference between adjacent isotope peaks for z=2, using 13C spacing.
C13_SPACING_DA = 1.00335483507
EXPECTED_Z2_SPACING_DA = C13_SPACING_DA / 2.0


def finite_float(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return math.nan
    return x if math.isfinite(x) else math.nan


def median(values):
    return float(statistics.median(values)) if values else math.nan


def percentile(values, p):
    return float(np.percentile(values, p)) if values else math.nan


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-table", type=Path, default=DEFAULT_TRACE)
    parser.add_argument("--feature-id", default=DEFAULT_FEATURE)
    parser.add_argument("--traces", default="1,2,3,4,5,6")
    parser.add_argument("--charge", type=int, default=2)
    parser.add_argument(
        "--spacing-tolerance-mda",
        type=float,
        default=3.0,
        help="Diagnostic tolerance around expected spacing in mDa; default 3",
    )
    parser.add_argument(
        "--out-prefix",
        default="trace_guided_isotope_spacing_v2",
    )
    args = parser.parse_args()

    if not args.trace_table.is_file():
        parser.error(f"Trace table not found: {args.trace_table}")
    if args.charge < 1:
        parser.error("--charge must be at least 1")
    if args.spacing_tolerance_mda < 0:
        parser.error("--spacing-tolerance-mda cannot be negative")

    selected_traces = {
        int(x.strip()) for x in args.traces.split(",") if x.strip()
    }
    expected_spacing = C13_SPACING_DA / args.charge
    tolerance_da = args.spacing_tolerance_mda / 1000.0

    by_scan = defaultdict(dict)
    feature_rows = 0
    duplicate_observations = []

    with args.trace_table.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "feature_id", "trace_index", "scan_index",
            "rt_sec", "observed_mz", "peak_intensity",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing trace-table columns: {sorted(missing)}")

        for row in reader:
            if str(row["feature_id"]) != str(args.feature_id):
                continue

            try:
                trace_index = int(row["trace_index"])
                scan_index = int(row["scan_index"])
            except (TypeError, ValueError):
                continue

            if trace_index not in selected_traces:
                continue

            mz = finite_float(row["observed_mz"])
            rt = finite_float(row["rt_sec"])
            intensity = finite_float(row["peak_intensity"])
            if not math.isfinite(mz):
                continue

            feature_rows += 1
            observation = {
                "feature_id": str(args.feature_id),
                "scan_index": scan_index,
                "trace_index": trace_index,
                "rt_sec": rt,
                "observed_mz": mz,
                "peak_intensity": intensity,
            }

            if trace_index in by_scan[scan_index]:
                duplicate_observations.append({
                    "scan_index": scan_index,
                    "trace_index": trace_index,
                    "previous_mz": by_scan[scan_index][trace_index]["observed_mz"],
                    "duplicate_mz": mz,
                    "previous_rt_sec": by_scan[scan_index][trace_index]["rt_sec"],
                    "duplicate_rt_sec": rt,
                })
                # Retain the first observation rather than silently replacing it.
                continue

            by_scan[scan_index][trace_index] = observation

    if feature_rows == 0:
        raise SystemExit(
            "No observations found. Check feature ID and trace indices."
        )

    # Same-scan adjacent trace spacing.
    spacing_rows = []
    pair_values = defaultdict(list)
    pair_missing_counts = defaultdict(int)
    trace_presence = defaultdict(int)

    ordered_traces = sorted(selected_traces)

    for scan_index in sorted(by_scan):
        scan = by_scan[scan_index]
        present = sorted(scan)
        for trace_index in present:
            trace_presence[trace_index] += 1

        # Evaluate each adjacent numeric trace pair, even when an intermediate
        # trace is absent, so missingness is explicit rather than skipped.
        for lower_trace, upper_trace in zip(ordered_traces, ordered_traces[1:]):
            lower = scan.get(lower_trace)
            upper = scan.get(upper_trace)

            row = {
                "feature_id": args.feature_id,
                "scan_index": scan_index,
                "lower_trace": lower_trace,
                "upper_trace": upper_trace,
                "lower_mz": lower["observed_mz"] if lower else math.nan,
                "upper_mz": upper["observed_mz"] if upper else math.nan,
                "lower_intensity": lower["peak_intensity"] if lower else math.nan,
                "upper_intensity": upper["peak_intensity"] if upper else math.nan,
                "rt_sec_lower": lower["rt_sec"] if lower else math.nan,
                "rt_sec_upper": upper["rt_sec"] if upper else math.nan,
                "observed_spacing_da": math.nan,
                "spacing_residual_mda": math.nan,
                "abs_spacing_residual_mda": math.nan,
                "within_spacing_tolerance": False,
                "status": "",
            }

            if lower is None or upper is None:
                pair_missing_counts[(lower_trace, upper_trace)] += 1
                row["status"] = "one_or_both_traces_absent"
            else:
                spacing = upper["observed_mz"] - lower["observed_mz"]
                residual_mda = (spacing - expected_spacing) * 1000.0
                abs_residual = abs(residual_mda)
                within = abs_residual <= args.spacing_tolerance_mda

                row.update({
                    "observed_spacing_da": spacing,
                    "spacing_residual_mda": residual_mda,
                    "abs_spacing_residual_mda": abs_residual,
                    "within_spacing_tolerance": within,
                    "status": (
                        "within_tolerance"
                        if within else "outside_tolerance"
                    ),
                })
                pair_values[(lower_trace, upper_trace)].append(row)

            spacing_rows.append(row)

    # Summary by adjacent pair.
    pair_summary = []
    for lower_trace, upper_trace in zip(ordered_traces, ordered_traces[1:]):
        rows = pair_values[(lower_trace, upper_trace)]
        residuals = [r["spacing_residual_mda"] for r in rows]
        abs_residuals = [abs(x) for x in residuals]
        within = sum(
            r["within_spacing_tolerance"] for r in rows
        )
        pair_summary.append({
            "feature_id": args.feature_id,
            "charge_hypothesis": args.charge,
            "lower_trace": lower_trace,
            "upper_trace": upper_trace,
            "expected_spacing_da": expected_spacing,
            "scans_with_both_traces": len(rows),
            "scans_missing_one_or_both": pair_missing_counts[
                (lower_trace, upper_trace)
            ],
            "within_tolerance_count": within,
            "within_tolerance_fraction": (
                within / len(rows) if rows else math.nan
            ),
            "median_observed_spacing_da": median([
                r["observed_spacing_da"] for r in rows
            ]),
            "median_signed_residual_mda": median(residuals),
            "median_abs_residual_mda": median(abs_residuals),
            "p95_abs_residual_mda": percentile(abs_residuals, 95),
            "max_abs_residual_mda": (
                max(abs_residuals) if abs_residuals else math.nan
            ),
        })

    # Trace coverage summary.
    trace_summary = []
    for trace_index in ordered_traces:
        observations = [
            by_scan[s][trace_index]
            for s in sorted(by_scan)
            if trace_index in by_scan[s]
        ]
        rts = [x["rt_sec"] for x in observations if math.isfinite(x["rt_sec"])]
        mzs = [x["observed_mz"] for x in observations]
        intensities = [
            x["peak_intensity"] for x in observations
            if math.isfinite(x["peak_intensity"])
        ]
        trace_summary.append({
            "feature_id": args.feature_id,
            "trace_index": trace_index,
            "observation_count": len(observations),
            "distinct_scan_count": len({x["scan_index"] for x in observations}),
            "median_observed_mz": median(mzs),
            "min_observed_mz": min(mzs) if mzs else math.nan,
            "max_observed_mz": max(mzs) if mzs else math.nan,
            "median_intensity": median(intensities),
            "intensity_sum": sum(intensities),
            "min_rt_sec": min(rts) if rts else math.nan,
            "max_rt_sec": max(rts) if rts else math.nan,
        })

    report_dir = PROJECT / "reports"
    spacing_path = report_dir / f"{args.out_prefix}_same_scan_pairs.tsv"
    pair_summary_path = report_dir / f"{args.out_prefix}_pair_summary.tsv"
    trace_summary_path = report_dir / f"{args.out_prefix}_trace_summary.tsv"
    duplicates_path = report_dir / f"{args.out_prefix}_duplicate_observations.tsv"

    spacing_fields = [
        "feature_id", "scan_index", "lower_trace", "upper_trace",
        "lower_mz", "upper_mz", "lower_intensity", "upper_intensity",
        "rt_sec_lower", "rt_sec_upper", "observed_spacing_da",
        "spacing_residual_mda", "abs_spacing_residual_mda",
        "within_spacing_tolerance", "status",
    ]
    pair_summary_fields = [
        "feature_id", "charge_hypothesis", "lower_trace", "upper_trace",
        "expected_spacing_da", "scans_with_both_traces",
        "scans_missing_one_or_both", "within_tolerance_count",
        "within_tolerance_fraction", "median_observed_spacing_da",
        "median_signed_residual_mda", "median_abs_residual_mda",
        "p95_abs_residual_mda", "max_abs_residual_mda",
    ]
    trace_summary_fields = [
        "feature_id", "trace_index", "observation_count",
        "distinct_scan_count", "median_observed_mz",
        "min_observed_mz", "max_observed_mz", "median_intensity",
        "intensity_sum", "min_rt_sec", "max_rt_sec",
    ]

    write_tsv(spacing_path, spacing_rows, spacing_fields)
    write_tsv(pair_summary_path, pair_summary, pair_summary_fields)
    write_tsv(trace_summary_path, trace_summary, trace_summary_fields)
    write_tsv(
        duplicates_path,
        duplicate_observations,
        [
            "scan_index", "trace_index", "previous_mz", "duplicate_mz",
            "previous_rt_sec", "duplicate_rt_sec",
        ],
    )

    print(f"Feature ID: {args.feature_id}")
    print(f"Trace-table rows read for selected traces: {feature_rows}")
    print(f"Distinct scan indices: {len(by_scan)}")
    print(f"Charge hypothesis: {args.charge}")
    print(f"Expected adjacent isotope spacing: {expected_spacing:.9f} Da")
    print(f"Spacing tolerance: ±{args.spacing_tolerance_mda:.3f} mDa")
    print(f"Duplicate scan/trace observations: {len(duplicate_observations)}")

    print("\nTrace coverage")
    print("trace observations distinct_scans median_mz")
    for row in trace_summary:
        print(
            f"{row['trace_index']:>5} "
            f"{row['observation_count']:>11} "
            f"{row['distinct_scan_count']:>13} "
            f"{row['median_observed_mz']:.6f}"
        )

    print("\nAdjacent-trace same-scan spacing")
    print(
        "pair both_present missing median_spacing "
        "median_abs_residual_mDa within_tolerance"
    )
    for row in pair_summary:
        print(
            f"{row['lower_trace']}-{row['upper_trace']} "
            f"{row['scans_with_both_traces']:>12} "
            f"{row['scans_missing_one_or_both']:>7} "
            f"{row['median_observed_spacing_da']:.6f} "
            f"{row['median_abs_residual_mda']:.3f} "
            f"{row['within_tolerance_count']}/"
            f"{row['scans_with_both_traces']}"
        )

    print("\nOutputs:")
    print(f"  {spacing_path}")
    print(f"  {pair_summary_path}")
    print(f"  {trace_summary_path}")
    print(f"  {duplicates_path}")

    print("\nCaution: spacing consistency is necessary but not sufficient to")
    print("establish isotope identity, monoisotopic status, or peptide linkage.")


if __name__ == "__main__":
    main()
