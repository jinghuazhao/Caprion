#!/usr/bin/env python3
"""Audit scan-wise charge-2 isotope spacing for traces 4–6."""

from pathlib import Path
import csv
import statistics

ROOT = Path.cwd()
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
REPORTS = ROOT / "reports"
REPORTS.mkdir(exist_ok=True)

FEATURE_ID = "10275997325121874585"
SCANS = [5852, 5865, 5878, 5891, 5904, 5917]
SPACING = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


with TRACE_FILE.open(newline="") as f:
    source = list(csv.DictReader(f, delimiter="\t"))

trace = {}
for row in source:
    if row["feature_id"] != FEATURE_ID:
        continue
    scan = int(row["scan_index"])
    idx = int(row["trace_index"])
    if scan not in SCANS or idx not in (1, 2, 3, 4, 5, 6):
        continue
    trace[(scan, idx)] = {
        "mz": float(row["observed_mz"]),
        "intensity": float(row["peak_intensity"]),
        "rt": float(row["rt_sec"]),
    }

per_scan = []
spacing_rows = []

for scan in SCANS:
    core = [
        trace[(scan, i)] for i in (1, 2, 3, 4)
        if (scan, i) in trace
    ]
    if len(core) != 4:
        continue

    # The core traces define the envelope intercept independently per scan.
    intercept = statistics.mean(
        trace[(scan, i)]["mz"] - (i - 1) * SPACING
        for i in (1, 2, 3, 4)
    )

    residuals = {}
    for idx in (1, 2, 3, 4, 5, 6):
        obs = trace.get((scan, idx))
        if obs is None:
            continue
        predicted = intercept + (idx - 1) * SPACING
        residuals[idx] = (obs["mz"] - predicted) * 1000.0

        per_scan.append({
            "feature_id": FEATURE_ID,
            "scan_index_zero_based": scan,
            "trace_index": idx,
            "rt_sec": f"{obs['rt']:.5f}",
            "observed_mz": f"{obs['mz']:.9f}",
            "core_intercept_mz": f"{intercept:.9f}",
            "predicted_mz_sequential_index": f"{predicted:.9f}",
            "residual_mDa": f"{residuals[idx]:.6f}",
            "within_3mDa": abs(residuals[idx]) <= TOLERANCE_MDA,
            "intensity": f"{obs['intensity']:.3f}",
        })

    for lower, upper in ((3, 4), (4, 5), (5, 6)):
        lo = trace.get((scan, lower))
        hi = trace.get((scan, upper))
        if lo is None or hi is None:
            continue
        observed_spacing = hi["mz"] - lo["mz"]
        residual = (observed_spacing - SPACING) * 1000.0
        spacing_rows.append({
            "feature_id": FEATURE_ID,
            "scan_index_zero_based": scan,
            "lower_trace_index": lower,
            "upper_trace_index": upper,
            "lower_mz": f"{lo['mz']:.9f}",
            "upper_mz": f"{hi['mz']:.9f}",
            "observed_spacing_Da": f"{observed_spacing:.9f}",
            "expected_charge2_spacing_Da": f"{SPACING:.12f}",
            "spacing_residual_mDa": f"{residual:.6f}",
            "within_3mDa": abs(residual) <= TOLERANCE_MDA,
        })

per_scan_path = REPORTS / "envelope_spacing_traces4_6_v1_per_scan.tsv"
spacing_path = REPORTS / "envelope_spacing_traces4_6_v1_adjacent.tsv"
write_tsv(per_scan_path, per_scan)
write_tsv(spacing_path, spacing_rows)

print(f"Feature: {FEATURE_ID}")
print(f"Expected charge-2 spacing: {SPACING:.12f} Da")
print(f"Per-trace residual report: {per_scan_path}")
print(f"Adjacent-spacing report: {spacing_path}")

print("\nAdjacent spacing summary:")
for pair in ((3, 4), (4, 5), (5, 6)):
    rows = [
        r for r in spacing_rows
        if (r["lower_trace_index"], r["upper_trace_index"]) == pair
    ]
    residuals = [float(r["spacing_residual_mDa"]) for r in rows]
    if not residuals:
        continue
    passing = sum(abs(x) <= TOLERANCE_MDA for x in residuals)
    print(
        f"Traces {pair[0]}–{pair[1]}: n={len(rows)}, "
        f"median residual={statistics.median(residuals):.3f} mDa, "
        f"median absolute residual={statistics.median(map(abs, residuals)):.3f} mDa, "
        f"within ±3 mDa={passing}/{len(rows)}"
    )

print("\nPer-scan residuals for traces 4–6:")
for scan in SCANS:
    selected = {
        int(r["trace_index"]): float(r["residual_mDa"])
        for r in per_scan
        if int(r["scan_index_zero_based"]) == scan
        and int(r["trace_index"]) in (4, 5, 6)
    }
    if selected:
        formatted = ", ".join(
            f"trace {i}: {selected[i]:+.3f} mDa"
            for i in sorted(selected)
        )
        print(f"scan {scan}: {formatted}")

print(
    "\nInterpretation: the core traces 1–4 define the intercept. "
    "Spacing and residuals are diagnostic evidence only; neither alone "
    "establishes isotope-group membership."
)
