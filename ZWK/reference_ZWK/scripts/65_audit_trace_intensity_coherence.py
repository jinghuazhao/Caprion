#!/usr/bin/env python3
"""Audit intensity coherence between existing OpenMS traces.

Calculates pairwise intensity correlations only over scans where both traces
were observed. Missing trace observations remain missing and are never replaced
with zero.

Correlations are descriptive diagnostics; they do not prove isotope identity.
No original XLSX isotope-group assignments are read.
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


def finite_float(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return math.nan
    return x if math.isfinite(x) else math.nan


def median(values):
    return float(statistics.median(values)) if values else math.nan


def pearson(x, y):
    if len(x) < 3:
        return math.nan
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if np.std(x) == 0 or np.std(y) == 0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    if len(x) < 3:
        return math.nan
    rx = np.argsort(np.argsort(np.asarray(x, dtype=float)))
    ry = np.argsort(np.argsort(np.asarray(y, dtype=float)))
    return pearson(rx, ry)


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
    parser.add_argument(
        "--min-pairs",
        type=int,
        default=5,
        help="Minimum common scans for a pairwise summary flag; default 5",
    )
    parser.add_argument(
        "--out-prefix",
        default="trace_intensity_coherence_v2",
    )
    args = parser.parse_args()

    if not args.trace_table.is_file():
        parser.error(f"Trace table not found: {args.trace_table}")
    if args.min_pairs < 3:
        parser.error("--min-pairs must be at least 3")

    traces = sorted({
        int(x.strip()) for x in args.traces.split(",") if x.strip()
    })
    by_scan = defaultdict(dict)
    duplicate_rows = []

    with args.trace_table.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "feature_id", "trace_index", "scan_index",
            "rt_sec", "observed_mz", "peak_intensity",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")

        for row in reader:
            if str(row["feature_id"]) != str(args.feature_id):
                continue

            try:
                trace_index = int(row["trace_index"])
                scan_index = int(row["scan_index"])
            except (TypeError, ValueError):
                continue
            if trace_index not in traces:
                continue

            observation = {
                "scan_index": scan_index,
                "trace_index": trace_index,
                "rt_sec": finite_float(row["rt_sec"]),
                "observed_mz": finite_float(row["observed_mz"]),
                "peak_intensity": finite_float(row["peak_intensity"]),
            }

            if trace_index in by_scan[scan_index]:
                duplicate_rows.append({
                    "scan_index": scan_index,
                    "trace_index": trace_index,
                    "retained_intensity":
                        by_scan[scan_index][trace_index]["peak_intensity"],
                    "duplicate_intensity": observation["peak_intensity"],
                })
                continue

            by_scan[scan_index][trace_index] = observation

    if not by_scan:
        raise SystemExit("No observations found; check feature ID and trace table.")

    # Per-scan matrix: absent observations are blank/NaN, never zero-filled.
    scan_rows = []
    matrix = {trace: [] for trace in traces}
    for scan_index in sorted(by_scan):
        present = by_scan[scan_index]
        rts = [
            present[t]["rt_sec"] for t in traces
            if t in present and math.isfinite(present[t]["rt_sec"])
        ]
        row = {
            "feature_id": args.feature_id,
            "scan_index": scan_index,
            "rt_sec": median(rts),
            "traces_present": ",".join(map(str, sorted(present))),
        }
        for trace in traces:
            obs = present.get(trace)
            intensity = (
                obs["peak_intensity"]
                if obs and math.isfinite(obs["peak_intensity"])
                else math.nan
            )
            mz = (
                obs["observed_mz"]
                if obs and math.isfinite(obs["observed_mz"])
                else math.nan
            )
            row[f"trace_{trace}_mz"] = mz
            row[f"trace_{trace}_intensity"] = intensity
            matrix[trace].append(intensity)
        scan_rows.append(row)

    # Pairwise correlation over common observed scans only.
    pair_rows = []
    for lower, upper in zip(traces, traces[1:]):
        x, y, rts, scan_ids = [], [], [], []
        for scan_index in sorted(by_scan):
            left = by_scan[scan_index].get(lower)
            right = by_scan[scan_index].get(upper)
            if left is None or right is None:
                continue

            li = left["peak_intensity"]
            ri = right["peak_intensity"]
            if not (math.isfinite(li) and math.isfinite(ri)):
                continue
            if li <= 0 or ri <= 0:
                continue

            x.append(li)
            y.append(ri)
            scan_ids.append(scan_index)
            rt = left["rt_sec"]
            if math.isfinite(rt):
                rts.append(rt)

        n = len(x)
        p = pearson(x, y)
        s = spearman(x, y)
        log_p = pearson(np.log10(x), np.log10(y))
        status = (
            "sufficient_for_descriptive_summary"
            if n >= args.min_pairs
            else "limited_common_scans"
        )

        pair_rows.append({
            "feature_id": args.feature_id,
            "trace_a": lower,
            "trace_b": upper,
            "common_positive_intensity_scans": n,
            "first_common_scan_index": min(scan_ids) if scan_ids else "",
            "last_common_scan_index": max(scan_ids) if scan_ids else "",
            "pearson_raw_intensity": p,
            "spearman_intensity": s,
            "pearson_log10_intensity": log_p,
            "median_intensity_a": median(x),
            "median_intensity_b": median(y),
            "status": status,
        })

    # Overall per-trace coverage/intensity summary.
    trace_rows = []
    for trace in traces:
        obs = [
            by_scan[scan][trace]
            for scan in sorted(by_scan)
            if trace in by_scan[scan]
        ]
        valid = [
            x["peak_intensity"] for x in obs
            if math.isfinite(x["peak_intensity"])
        ]
        valid_positive = [x for x in valid if x > 0]
        rts = [
            x["rt_sec"] for x in obs
            if math.isfinite(x["rt_sec"])
        ]
        trace_rows.append({
            "feature_id": args.feature_id,
            "trace_index": trace,
            "observations": len(obs),
            "distinct_scan_indices": len({x["scan_index"] for x in obs}),
            "positive_intensity_observations": len(valid_positive),
            "median_intensity": median(valid),
            "median_positive_intensity": median(valid_positive),
            "intensity_sum": sum(valid),
            "min_rt_sec": min(rts) if rts else math.nan,
            "max_rt_sec": max(rts) if rts else math.nan,
        })

    report_dir = PROJECT / "reports"
    prefix = args.out_prefix
    per_scan_path = report_dir / f"{prefix}_per_scan.tsv"
    pair_path = report_dir / f"{prefix}_adjacent_pair_summary.tsv"
    trace_path = report_dir / f"{prefix}_trace_summary.tsv"
    duplicate_path = report_dir / f"{prefix}_duplicate_observations.tsv"

    per_scan_fields = [
        "feature_id", "scan_index", "rt_sec", "traces_present"
    ]
    for trace in traces:
        per_scan_fields.extend([
            f"trace_{trace}_mz",
            f"trace_{trace}_intensity",
        ])

    write_tsv(per_scan_path, scan_rows, per_scan_fields)
    write_tsv(
        pair_path,
        pair_rows,
        [
            "feature_id", "trace_a", "trace_b",
            "common_positive_intensity_scans",
            "first_common_scan_index", "last_common_scan_index",
            "pearson_raw_intensity", "spearman_intensity",
            "pearson_log10_intensity",
            "median_intensity_a", "median_intensity_b", "status",
        ],
    )
    write_tsv(
        trace_path,
        trace_rows,
        [
            "feature_id", "trace_index", "observations",
            "distinct_scan_indices", "positive_intensity_observations",
            "median_intensity", "median_positive_intensity",
            "intensity_sum", "min_rt_sec", "max_rt_sec",
        ],
    )
    write_tsv(
        duplicate_path,
        duplicate_rows,
        [
            "scan_index", "trace_index",
            "retained_intensity", "duplicate_intensity",
        ],
    )

    print(f"Feature ID: {args.feature_id}")
    print(f"Trace-table observations: {sum(map(len, by_scan.values()))}")
    print(f"Distinct scan indices: {len(by_scan)}")
    print(f"Duplicate scan/trace rows skipped: {len(duplicate_rows)}")

    print("\nAdjacent-trace intensity coherence")
    print("pair n_common Pearson Spearman log10_Pearson status")
    for row in pair_rows:
        print(
            f"{row['trace_a']}-{row['trace_b']} "
            f"{row['common_positive_intensity_scans']:>8} "
            f"{row['pearson_raw_intensity']:.3f} "
            f"{row['spearman_intensity']:.3f} "
            f"{row['pearson_log10_intensity']:.3f} "
            f"{row['status']}"
        )

    print("\nOutputs:")
    print(f"  {per_scan_path}")
    print(f"  {pair_path}")
    print(f"  {trace_path}")
    print(f"  {duplicate_path}")

    print("\nCautions:")
    print("- Pairwise correlations use only common scans with positive intensities.")
    print("- Correlation is descriptive and sensitive to the short RT window.")
    print("- Shared chromatographic trends do not prove isotope identity.")
    print("- No monoisotopic peak or independent IG is assigned.")


if __name__ == "__main__":
    main()
