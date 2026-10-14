#!/usr/bin/env python3
"""Build a scan-aligned m/z and intensity audit for one candidate envelope."""

from pathlib import Path
import csv
import math
import statistics

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
INPUT = ROOT / "openms_pilot/isotope_trace_peaks.tsv"

FEATURE_ID = "10275997325121874585"
TRACE_INDICES = [1, 2, 3, 4, 5, 6]
CORE_TRACES = [1, 2, 3, 4]
SCANS = [5852, 5865, 5878, 5891, 5904, 5917]

# Carbon-13 isotope spacing divided by charge 2.
THEORETICAL_SPACING = 1.00335483507 / 2.0


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def fmt(value):
    return "" if value is None else f"{value:.9f}"


def median(values):
    return statistics.median(values) if values else None


if not INPUT.exists():
    raise FileNotFoundError(INPUT)

source = read_tsv(INPUT)

observations = {}
for row in source:
    if row["feature_id"] != FEATURE_ID:
        continue

    scan = int(row["scan_index"])
    trace = int(row["trace_index"])
    if scan not in SCANS or trace not in TRACE_INDICES:
        continue

    key = (scan, trace)
    if key in observations:
        raise ValueError(f"Duplicate scan/trace observation: {key}")

    observations[key] = {
        "rt": float(row["rt_sec"]),
        "mz": float(row["observed_mz"]),
        "intensity": float(row["peak_intensity"]),
    }

# Estimate m/z intercept independently within each scan from traces 1–4.
# The linear model is mz = intercept + isotope_index * spacing.
# The intercept is the median of the four core-adjusted m/z values.
scan_models = {}
for scan in SCANS:
    core = [
        observations[(scan, trace)]
        for trace in CORE_TRACES
        if (scan, trace) in observations
    ]

    if len(core) != len(CORE_TRACES):
        scan_models[scan] = {
            "intercept": None,
            "core_spacing_median": None,
            "core_spacing_mad": None,
        }
        continue

    adjusted = [
        observations[(scan, trace)]["mz"]
        - (trace - 1) * THEORETICAL_SPACING
        for trace in CORE_TRACES
    ]
    intercept = median(adjusted)

    spacings = [
        observations[(scan, trace + 1)]["mz"]
        - observations[(scan, trace)]["mz"]
        for trace in range(1, 4)
    ]

    scan_models[scan] = {
        "intercept": intercept,
        "core_spacing_median": median(spacings),
        "core_spacing_mad": median(
            [abs(v - median(spacings)) for v in spacings]
        ),
    }

# Scan-aligned long-form audit.
detail = []
for scan in SCANS:
    model = scan_models[scan]
    for trace in TRACE_INDICES:
        obs = observations.get((scan, trace))
        expected = (
            model["intercept"] + (trace - 1) * THEORETICAL_SPACING
            if model["intercept"] is not None else None
        )
        residual_mda = (
            (obs["mz"] - expected) * 1000.0
            if obs is not None and expected is not None else None
        )

        detail.append({
            "feature_id": FEATURE_ID,
            "scan_index_zero_based": scan,
            "trace_index": trace,
            "rt_sec": fmt(obs["rt"]) if obs else "",
            "observed_mz": fmt(obs["mz"]) if obs else "",
            "observed_intensity": fmt(obs["intensity"]) if obs else "",
            "expected_mz_from_core_model": fmt(expected),
            "residual_from_core_model_mDa": (
                f"{residual_mda:.6f}" if residual_mda is not None else ""
            ),
            "abs_residual_mDa": (
                f"{abs(residual_mda):.6f}"
                if residual_mda is not None else ""
            ),
            "within_3mDa": (
                abs(residual_mda) <= 3.0
                if residual_mda is not None else ""
            ),
            "core_intercept_mz": fmt(model["intercept"]),
            "core_median_spacing_Da": fmt(model["core_spacing_median"]),
            "core_spacing_MAD_Da": fmt(model["core_spacing_mad"]),
        })

detail_path = REPORTS / "scan_aligned_envelope_audit_v1.tsv"
write_tsv(detail_path, detail)

# Summarise residual behaviour by trace.
summary = []
for trace in TRACE_INDICES:
    rows = [
        row for row in detail
        if int(row["trace_index"]) == trace
        and row["abs_residual_mDa"] != ""
    ]
    residuals = [float(row["abs_residual_mDa"]) for row in rows]
    passing = sum(float(row["abs_residual_mDa"]) <= 3.0 for row in rows)

    summary.append({
        "trace_index": trace,
        "n_observations_with_core_model": len(rows),
        "n_within_3mDa": passing,
        "fraction_within_3mDa": (
            f"{passing / len(rows):.6f}" if rows else ""
        ),
        "median_abs_residual_mDa": (
            f"{median(residuals):.6f}" if residuals else ""
        ),
        "max_abs_residual_mDa": (
            f"{max(residuals):.6f}" if residuals else ""
        ),
    })

summary_path = REPORTS / "scan_aligned_envelope_summary_v1.tsv"
write_tsv(summary_path, summary)

# Adjacent trace spacings by scan.
spacing_rows = []
for scan in SCANS:
    for lower in range(1, 6):
        a = observations.get((scan, lower))
        b = observations.get((scan, lower + 1))
        spacing = b["mz"] - a["mz"] if a and b else None
        residual_mda = (
            (spacing - THEORETICAL_SPACING) * 1000.0
            if spacing is not None else None
        )
        spacing_rows.append({
            "scan_index_zero_based": scan,
            "lower_trace": lower,
            "upper_trace": lower + 1,
            "observed_spacing_Da": fmt(spacing),
            "theoretical_spacing_Da": fmt(THEORETICAL_SPACING),
            "spacing_residual_mDa": (
                f"{residual_mda:.6f}"
                if residual_mda is not None else ""
            ),
        })

spacing_path = REPORTS / "scan_aligned_adjacent_spacing_v1.tsv"
write_tsv(spacing_path, spacing_rows)

print(f"Feature: {FEATURE_ID}")
print(f"Theoretical charge-2 spacing: {THEORETICAL_SPACING:.12f} Da")
print(f"Core traces used for intercept: {CORE_TRACES}")
print("\nResidual summary against the core-derived model:")
for row in summary:
    print(
        f"trace {row['trace_index']}: "
        f"{row['n_within_3mDa']}/{row['n_observations_with_core_model']} "
        f"within 3 mDa; "
        f"median |residual|={row['median_abs_residual_mDa']} mDa; "
        f"max |residual|={row['max_abs_residual_mDa']} mDa"
    )

print("\nReports:")
print(detail_path)
print(summary_path)
print(spacing_path)

print(
    "\nCaution: traces 1–4 define the reference model, so their residuals "
    "are not independent validation. Traces 5–6 are diagnostic only; "
    "this script does not establish isotope identity or an independent IG."
)
