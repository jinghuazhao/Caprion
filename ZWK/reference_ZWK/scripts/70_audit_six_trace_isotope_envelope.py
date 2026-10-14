#!/usr/bin/env python3
"""Audit the six-trace isotope-like envelope for one OpenMS feature.

Outputs:
  reports/six_trace_envelope_v2_per_scan.tsv
  reports/six_trace_envelope_v2_pair_summary.tsv
  reports/six_trace_envelope_v2_trace_summary.tsv
  reports/six_trace_envelope_v2_feature_summary.tsv

Mass-spacing and intensity evidence are reported separately. No automatic
isotope-group assignments are made.
"""

from __future__ import annotations

import csv
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
FEATURE_ID = "10275997325121874585"

ISOTOPE_SPACING = 1.00335483507 / 2.0
SPACING_TOLERANCE_MDA = 3.0

OUT_SCAN = ROOT / "reports/six_trace_envelope_v2_per_scan.tsv"
OUT_PAIR = ROOT / "reports/six_trace_envelope_v2_pair_summary.tsv"
OUT_TRACE = ROOT / "reports/six_trace_envelope_v2_trace_summary.tsv"
OUT_FEATURE = ROOT / "reports/six_trace_envelope_v2_feature_summary.tsv"


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
        value = float(value)
        return value if math.isfinite(value) else default
    except (TypeError, ValueError):
        return default


def median(values: list[float]) -> float:
    return statistics.median(values) if values else float("nan")


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) < 3 or len(x) != len(y):
        return float("nan")
    mx = statistics.mean(x)
    my = statistics.mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    denominator = math.sqrt(sum(v * v for v in dx) * sum(v * v for v in dy))
    if denominator == 0:
        return float("nan")
    return sum(a * b for a, b in zip(dx, dy)) / denominator


def fmt_bool(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def main() -> int:
    if not TRACE_FILE.exists():
        print(f"ERROR: missing trace file: {TRACE_FILE}", file=sys.stderr)
        return 1

    rows = [
        row for row in read_tsv(TRACE_FILE)
        if row.get("feature_id") == FEATURE_ID
    ]
    if not rows:
        print(f"ERROR: feature {FEATURE_ID} not found", file=sys.stderr)
        return 1

    # Key by scan and trace. Repeated observations are retained for auditing.
    by_scan: dict[int, dict[int, dict]] = {}
    duplicate_keys = []

    for row in rows:
        scan = int(row["scan_index"])
        trace = int(row["trace_index"])
        if trace < 1 or trace > 6:
            continue

        if trace in by_scan.setdefault(scan, {}):
            duplicate_keys.append({
                "scan_index": scan,
                "trace_index": trace,
                "issue": "duplicate scan/trace observation",
            })
            continue

        by_scan[scan][trace] = row

    scan_output = []
    for scan in sorted(by_scan):
        trace_map = by_scan[scan]
        out = {
            "feature_id": FEATURE_ID,
            "scan_index": scan,
            "n_traces_present": len(trace_map),
            "trace_indices_present": ",".join(
                str(t) for t in sorted(trace_map)
            ),
        }

        for trace in range(1, 7):
            row = trace_map.get(trace)
            out[f"trace_{trace}_mz"] = (
                num(row.get("observed_mz")) if row else float("nan")
            )
            out[f"trace_{trace}_intensity"] = (
                num(row.get("peak_intensity")) if row else float("nan")
            )

        for trace in range(1, 6):
            a = out[f"trace_{trace}_mz"]
            b = out[f"trace_{trace + 1}_mz"]
            ia = out[f"trace_{trace}_intensity"]
            ib = out[f"trace_{trace + 1}_intensity"]

            present = (
                math.isfinite(a) and math.isfinite(b)
                and math.isfinite(ia) and math.isfinite(ib)
                and ia > 0 and ib > 0
            )
            out[f"pair_{trace}_{trace + 1}_present"] = fmt_bool(present)

            if present:
                spacing = b - a
                residual = (spacing - ISOTOPE_SPACING) * 1000.0
                out[f"pair_{trace}_{trace + 1}_spacing_da"] = spacing
                out[f"pair_{trace}_{trace + 1}_residual_mda"] = residual
                out[f"pair_{trace}_{trace + 1}_within_3mda"] = fmt_bool(
                    abs(residual) <= SPACING_TOLERANCE_MDA
                )
                out[f"pair_{trace}_{trace + 1}_intensity_ratio"] = ib / ia
            else:
                for suffix in (
                    "spacing_da", "residual_mda", "intensity_ratio"
                ):
                    out[f"pair_{trace}_{trace + 1}_{suffix}"] = float("nan")
                out[f"pair_{trace}_{trace + 1}_within_3mda"] = ""

        scan_output.append(out)

    scan_fields = [
        "feature_id", "scan_index", "n_traces_present", "trace_indices_present"
    ]
    for trace in range(1, 7):
        scan_fields.extend([
            f"trace_{trace}_mz",
            f"trace_{trace}_intensity",
        ])
    for trace in range(1, 6):
        scan_fields.extend([
            f"pair_{trace}_{trace + 1}_present",
            f"pair_{trace}_{trace + 1}_spacing_da",
            f"pair_{trace}_{trace + 1}_residual_mda",
            f"pair_{trace}_{trace + 1}_within_3mda",
            f"pair_{trace}_{trace + 1}_intensity_ratio",
        ])

    write_tsv(OUT_SCAN, scan_output, scan_fields)

    pair_summary = []
    for trace in range(1, 6):
        pair = f"{trace}-{trace + 1}"
        residuals = []
        ratios = []
        x_intensities = []
        y_intensities = []

        for row in scan_output:
            if row[f"pair_{trace}_{trace + 1}_present"] != "TRUE":
                continue

            residual = num(row[f"pair_{trace}_{trace + 1}_residual_mda"])
            ratio = num(row[f"pair_{trace}_{trace + 1}_intensity_ratio"])
            ia = num(row[f"trace_{trace}_intensity"])
            ib = num(row[f"trace_{trace + 1}_intensity"])

            if math.isfinite(residual):
                residuals.append(residual)
            if math.isfinite(ratio):
                ratios.append(ratio)
            if ia > 0 and ib > 0:
                x_intensities.append(math.log10(ia))
                y_intensities.append(math.log10(ib))

        within = sum(
            abs(value) <= SPACING_TOLERANCE_MDA for value in residuals
        )
        pair_summary.append({
            "feature_id": FEATURE_ID,
            "trace_pair": pair,
            "n_common_positive_scans": len(residuals),
            "median_spacing_residual_mda": median(residuals),
            "median_absolute_spacing_residual_mda": median(
                [abs(value) for value in residuals]
            ),
            "n_within_3mda": within,
            "fraction_within_3mda": (
                within / len(residuals) if residuals else float("nan")
            ),
            "median_intensity_ratio_b_over_a": median(ratios),
            "log10_intensity_pearson": pearson(x_intensities, y_intensities),
        })

    pair_fields = [
        "feature_id", "trace_pair", "n_common_positive_scans",
        "median_spacing_residual_mda",
        "median_absolute_spacing_residual_mda",
        "n_within_3mda", "fraction_within_3mda",
        "median_intensity_ratio_b_over_a", "log10_intensity_pearson",
    ]
    write_tsv(OUT_PAIR, pair_summary, pair_fields)

    trace_summary = []
    for trace in range(1, 7):
        mz_values = []
        intensities = []
        scans = []

        for scan, trace_map in by_scan.items():
            row = trace_map.get(trace)
            if row is None:
                continue
            mz = num(row.get("observed_mz"))
            intensity = num(row.get("peak_intensity"))
            if not math.isfinite(mz) or not math.isfinite(intensity):
                continue
            mz_values.append(mz)
            intensities.append(intensity)
            scans.append(scan)

        trace_summary.append({
            "feature_id": FEATURE_ID,
            "trace_index": trace,
            "n_observations": len(mz_values),
            "n_distinct_scans": len(set(scans)),
            "median_mz": median(mz_values),
            "min_mz": min(mz_values) if mz_values else float("nan"),
            "max_mz": max(mz_values) if mz_values else float("nan"),
            "intensity_sum": sum(intensities),
            "median_intensity": median(intensities),
        })

    trace_fields = [
        "feature_id", "trace_index", "n_observations", "n_distinct_scans",
        "median_mz", "min_mz", "max_mz", "intensity_sum", "median_intensity",
    ]
    write_tsv(OUT_TRACE, trace_summary, trace_fields)

    # Feature-level descriptive summary. No acceptance rule is applied.
    feature_summary = [{
        "feature_id": FEATURE_ID,
        "n_scan_indices": len(by_scan),
        "n_observations": sum(len(v) for v in by_scan.values()),
        "n_duplicate_scan_trace_keys": len(duplicate_keys),
        "expected_charge2_spacing_da": ISOTOPE_SPACING,
        "spacing_tolerance_mda": SPACING_TOLERANCE_MDA,
        "interpretation": (
            "Descriptive six-trace audit; isotope identities remain provisional"
        ),
    }]
    write_tsv(
        OUT_FEATURE,
        feature_summary,
        [
            "feature_id", "n_scan_indices", "n_observations",
            "n_duplicate_scan_trace_keys", "expected_charge2_spacing_da",
            "spacing_tolerance_mda", "interpretation",
        ],
    )

    print(f"Feature: {FEATURE_ID}")
    print(f"Distinct scans: {len(by_scan)}")
    print(f"Trace observations: {sum(len(v) for v in by_scan.values())}")
    print(f"Duplicate scan/trace keys: {len(duplicate_keys)}")
    print(f"Per-scan table: {OUT_SCAN}")
    print(f"Adjacent-pair summary: {OUT_PAIR}")
    print(f"Per-trace summary: {OUT_TRACE}")
    print(f"Feature summary: {OUT_FEATURE}")
    print()
    print(
        "pair\tn\tmedian_abs_residual_mDa\twithin_3mDa\t"
        "fraction_within_3mDa\tlog10_intensity_r"
    )
    for row in pair_summary:
        print(
            f"{row['trace_pair']}\t"
            f"{row['n_common_positive_scans']}\t"
            f"{num(row['median_absolute_spacing_residual_mda']):.3f}\t"
            f"{row['n_within_3mda']}\t"
            f"{num(row['fraction_within_3mda']):.3f}\t"
            f"{num(row['log10_intensity_pearson']):.3f}"
        )

    print()
    print(
        "Caution: intensity correlation and co-elution are descriptive; "
        "they do not override inconsistent isotope spacing."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
