#!/usr/bin/env python3
"""Combine same-scan spacing and intensity evidence for candidate MS1 traces.

Inputs are the existing Script 64 and Script 65 reports. The script creates:
  * one per-scan, adjacent-pair evidence table;
  * a pair-level diagnostic summary;
  * a ranked table of the largest absolute spacing residuals.

Missing trace observations remain missing. No original XLSX assignments are
read, and no isotope-group or monoisotopic assignments are made.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

PROJECT = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
REPORTS = PROJECT / "reports"


def number(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return math.nan
    return x if math.isfinite(x) else math.nan


def integer(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def median(values):
    return float(statistics.median(values)) if values else math.nan


def read_tsv(path):
    if not path.is_file():
        raise FileNotFoundError(f"Required report not found: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
    parser.add_argument(
        "--spacing",
        type=Path,
        default=REPORTS / "trace_guided_isotope_spacing_v2_same_scan_pairs.tsv",
    )
    parser.add_argument(
        "--per-scan-intensity",
        type=Path,
        default=REPORTS / "trace_intensity_coherence_v2_per_scan.tsv",
    )
    parser.add_argument(
        "--pair-intensity",
        type=Path,
        default=REPORTS / "trace_intensity_coherence_v2_adjacent_pair_summary.tsv",
    )
    parser.add_argument("--feature-id", default="10275997325121874585")
    parser.add_argument("--spacing-tolerance-mda", type=float, default=3.0)
    parser.add_argument("--worst-n", type=int, default=20)
    parser.add_argument("--out-prefix", default="combined_trace_evidence_v2")
    args = parser.parse_args()

    if args.spacing_tolerance_mda < 0:
        parser.error("--spacing-tolerance-mda cannot be negative")
    if args.worst_n < 1:
        parser.error("--worst-n must be positive")

    spacing_rows = read_tsv(args.spacing)
    intensity_scan_rows = read_tsv(args.per_scan_intensity)
    intensity_pair_rows = read_tsv(args.pair_intensity)

    # Restrict spacing records to the requested feature.
    spacing_rows = [
        r for r in spacing_rows
        if str(r.get("feature_id", "")) == str(args.feature_id)
    ]
    if not spacing_rows:
        raise SystemExit(
            f"No spacing rows for feature {args.feature_id} in {args.spacing}"
        )

    scan_intensity = {}
    for row in intensity_scan_rows:
        if str(row.get("feature_id", "")) != str(args.feature_id):
            continue
        idx = integer(row.get("scan_index"))
        if idx is not None:
            scan_intensity[idx] = row

    pair_intensity = {}
    for row in intensity_pair_rows:
        if str(row.get("feature_id", "")) != str(args.feature_id):
            continue
        a = integer(row.get("trace_a"))
        b = integer(row.get("trace_b"))
        if a is not None and b is not None:
            pair_intensity[(a, b)] = row

    combined = []
    for s in spacing_rows:
        idx = integer(s.get("scan_index"))
        a = integer(s.get("lower_trace"))
        b = integer(s.get("upper_trace"))
        if idx is None or a is None or b is None:
            continue

        intensity_row = scan_intensity.get(idx, {})
        mz_a = number(s.get("lower_mz"))
        mz_b = number(s.get("upper_mz"))
        spacing = number(s.get("observed_spacing_da"))
        residual = number(s.get("spacing_residual_mda"))
        abs_residual = number(s.get("abs_spacing_residual_mda"))
        ia = number(intensity_row.get(f"trace_{a}_intensity"))
        ib = number(intensity_row.get(f"trace_{b}_intensity"))

        intensity_ratio = math.nan
        log10_ratio = math.nan
        if math.isfinite(ia) and math.isfinite(ib) and ia > 0 and ib > 0:
            intensity_ratio = ib / ia
            log10_ratio = math.log10(intensity_ratio)

        pair_summary = pair_intensity.get((a, b), {})
        within = (
            math.isfinite(abs_residual)
            and abs_residual <= args.spacing_tolerance_mda
        )

        if not (math.isfinite(mz_a) and math.isfinite(mz_b)):
            evidence_flag = "missing_one_or_both_traces"
        elif not math.isfinite(abs_residual):
            evidence_flag = "spacing_not_calculated"
        elif within:
            evidence_flag = "spacing_within_tolerance"
        else:
            evidence_flag = "spacing_outside_tolerance"

        combined.append({
            "feature_id": args.feature_id,
            "scan_index": idx,
            "trace_a": a,
            "trace_b": b,
            "rt_sec": (
                number(s.get("rt_sec_lower"))
                if math.isfinite(number(s.get("rt_sec_lower")))
                else number(s.get("rt_sec_upper"))
            ),
            "mz_a": mz_a,
            "mz_b": mz_b,
            "intensity_a": ia,
            "intensity_b": ib,
            "observed_spacing_da": spacing,
            "spacing_residual_mda": residual,
            "abs_spacing_residual_mda": abs_residual,
            "within_spacing_tolerance": within,
            "intensity_ratio_b_over_a": intensity_ratio,
            "log10_intensity_ratio_b_over_a": log10_ratio,
            "pair_common_scan_count": integer(
                pair_summary.get("common_positive_intensity_scans")
            ),
            "pearson_raw_intensity": number(
                pair_summary.get("pearson_raw_intensity")
            ),
            "spearman_intensity": number(
                pair_summary.get("spearman_intensity")
            ),
            "pearson_log10_intensity": number(
                pair_summary.get("pearson_log10_intensity")
            ),
            "evidence_flag": evidence_flag,
        })

    combined.sort(key=lambda r: (r["trace_a"], r["trace_b"], r["scan_index"]))

    # Pair-level summary from the combined scan-level evidence.
    grouped = defaultdict(list)
    for row in combined:
        grouped[(row["trace_a"], row["trace_b"])].append(row)

    pair_summary_out = []
    for (a, b), rows in sorted(grouped.items()):
        present = [
            r for r in rows
            if math.isfinite(r["abs_spacing_residual_mda"])
        ]
        within_count = sum(r["within_spacing_tolerance"] for r in present)
        abs_errors = [r["abs_spacing_residual_mda"] for r in present]
        signed_errors = [r["spacing_residual_mda"] for r in present]
        ratios = [
            r["intensity_ratio_b_over_a"] for r in rows
            if math.isfinite(r["intensity_ratio_b_over_a"])
        ]
        psummary = pair_intensity.get((a, b), {})

        pair_summary_out.append({
            "feature_id": args.feature_id,
            "trace_a": a,
            "trace_b": b,
            "scan_rows": len(rows),
            "scans_with_spacing": len(present),
            "scans_within_tolerance": within_count,
            "fraction_within_tolerance": (
                within_count / len(present) if present else math.nan
            ),
            "median_spacing_residual_mda": median(signed_errors),
            "median_abs_spacing_residual_mda": median(abs_errors),
            "max_abs_spacing_residual_mda": (
                max(abs_errors) if abs_errors else math.nan
            ),
            "median_intensity_ratio_b_over_a": median(ratios),
            "common_positive_intensity_scans": integer(
                psummary.get("common_positive_intensity_scans")
            ),
            "pearson_raw_intensity": number(
                psummary.get("pearson_raw_intensity")
            ),
            "spearman_intensity": number(
                psummary.get("spearman_intensity")
            ),
            "pearson_log10_intensity": number(
                psummary.get("pearson_log10_intensity")
            ),
        })

    # Worst residuals are ranked only when both traces are present.
    worst = [
        r for r in combined
        if math.isfinite(r["abs_spacing_residual_mda"])
    ]
    worst.sort(key=lambda r: r["abs_spacing_residual_mda"], reverse=True)
    worst = worst[:args.worst_n]

    out_fields = [
        "feature_id", "scan_index", "rt_sec", "trace_a", "trace_b",
        "mz_a", "mz_b", "intensity_a", "intensity_b",
        "observed_spacing_da", "spacing_residual_mda",
        "abs_spacing_residual_mda", "within_spacing_tolerance",
        "intensity_ratio_b_over_a", "log10_intensity_ratio_b_over_a",
        "pair_common_scan_count", "pearson_raw_intensity",
        "spearman_intensity", "pearson_log10_intensity", "evidence_flag",
    ]
    pair_fields = [
        "feature_id", "trace_a", "trace_b", "scan_rows",
        "scans_with_spacing", "scans_within_tolerance",
        "fraction_within_tolerance", "median_spacing_residual_mda",
        "median_abs_spacing_residual_mda", "max_abs_spacing_residual_mda",
        "median_intensity_ratio_b_over_a", "common_positive_intensity_scans",
        "pearson_raw_intensity", "spearman_intensity",
        "pearson_log10_intensity",
    ]

    prefix = args.out_prefix
    detailed_path = REPORTS / f"{prefix}_per_scan_pair.tsv"
    summary_path = REPORTS / f"{prefix}_pair_summary.tsv"
    worst_path = REPORTS / f"{prefix}_worst_spacing_residuals.tsv"

    write_tsv(detailed_path, combined, out_fields)
    write_tsv(summary_path, pair_summary_out, pair_fields)
    write_tsv(worst_path, worst, out_fields)

    print(f"Feature ID: {args.feature_id}")
    print(f"Combined scan/pair rows: {len(combined)}")
    print(f"Pair groups: {len(pair_summary_out)}")
    print(f"Spacing tolerance: ±{args.spacing_tolerance_mda:.3f} mDa")

    print("\nCombined pair summary")
    print(
        "pair scans_with_spacing within_tolerance median_abs_residual_mDa "
        "Pearson intensity median_intensity_ratio"
    )
    for row in pair_summary_out:
        print(
            f"{row['trace_a']}-{row['trace_b']} "
            f"{row['scans_with_spacing']:>5} "
            f"{row['scans_within_tolerance']}/"
            f"{row['scans_with_spacing']} "
            f"{row['median_abs_spacing_residual_mda']:.3f} "
            f"{row['pearson_raw_intensity']:.3f} "
            f"{row['median_intensity_ratio_b_over_a']:.3f}"
        )

    print("\nLargest spacing residuals")
    print("pair scan_index RT_sec abs_residual_mDa intensity_ratio")
    for row in worst:
        print(
            f"{row['trace_a']}-{row['trace_b']} "
            f"{row['scan_index']} {row['rt_sec']:.5f} "
            f"{row['abs_spacing_residual_mda']:.3f} "
            f"{row['intensity_ratio_b_over_a']:.3f}"
        )

    print("\nOutputs:")
    print(f"  {detailed_path}")
    print(f"  {summary_path}")
    print(f"  {worst_path}")

    print("\nCaution: these are joined diagnostics, not isotope assignments.")


if __name__ == "__main__":
    main()
