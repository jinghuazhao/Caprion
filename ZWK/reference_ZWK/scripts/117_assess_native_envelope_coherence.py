#!/usr/bin/env python3
"""Assess mass-drift coherence of candidate native isotope tracks."""

from pathlib import Path
import csv
import numpy as np

INPUT = Path(
    "reports/feature10275997325121874585_native_track_validation_v1.tsv"
)
OUT_SCAN = Path(
    "reports/feature10275997325121874585_envelope_coherence_scans_v1.tsv"
)
OUT_SUMMARY = Path(
    "reports/feature10275997325121874585_envelope_coherence_summary_v1.tsv"
)

CORE = [2, 3, 4, 5]
ALL_TRACKS = [2, 3, 4, 5, 6, 7]
THEORETICAL_SPACING = 1.00335483507 / 2

with INPUT.open(newline="") as f:
    records = list(csv.DictReader(f, delimiter="\t"))

by_scan = {}
for r in records:
    scan = int(r["scan_index"])
    track = int(r["track"])
    by_scan.setdefault(scan, {})[track] = r

scan_rows = []
for scan in sorted(by_scan):
    d = by_scan[scan]
    row = {
        "scan_index": scan,
        "rt": float(next(iter(d.values()))["rt"]),
    }

    for t in ALL_TRACKS:
        r = d.get(t)
        found = r is not None and r["observed_mz"] != ""
        row[f"track{t}_found"] = int(found)
        row[f"track{t}_mz"] = float(r["observed_mz"]) if found else ""
        row[f"track{t}_intensity"] = float(r["intensity"]) if found else ""
        row[f"track{t}_error_mDa"] = float(r["error_mDa"]) if found else ""

    # Adjacent spacings for the conservative core.
    for a, b in zip(CORE[:-1], CORE[1:]):
        ma, mb = row[f"track{a}_mz"], row[f"track{b}_mz"]
        if ma != "" and mb != "":
            spacing = mb - ma
            row[f"spacing_{a}_{b}_Da"] = spacing
            row[f"spacing_{a}_{b}_residual_mDa"] = (
                spacing - THEORETICAL_SPACING
            ) * 1000
        else:
            row[f"spacing_{a}_{b}_Da"] = ""
            row[f"spacing_{a}_{b}_residual_mDa"] = ""

    # A scan is core-complete only when all four core peaks are observed.
    row["core_complete"] = int(
        all(row[f"track{t}_found"] for t in CORE)
    )

    # Intensity ordering is descriptive only, not a pass/fail criterion.
    if row["core_complete"]:
        intensities = [row[f"track{t}_intensity"] for t in CORE]
        row["core_intensity_order"] = ">".join(
            str(t) for t in sorted(
                CORE, key=lambda t: row[f"track{t}_intensity"], reverse=True
            )
        )
        row["core_intensity_strictly_decreasing"] = int(
            all(intensities[i] > intensities[i + 1]
                for i in range(len(intensities) - 1))
        )
    else:
        row["core_intensity_order"] = ""
        row["core_intensity_strictly_decreasing"] = ""

    scan_rows.append(row)

scan_fields = list(scan_rows[0])
OUT_SCAN.parent.mkdir(parents=True, exist_ok=True)
with OUT_SCAN.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=scan_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(scan_rows)

# Assess whether mass shifts are shared across core tracks.
summary = []
for t in CORE:
    vals = [
        r for r in scan_rows
        if r[f"track{t}_found"]
    ]
    errors = np.array(
        [r[f"track{t}_error_mDa"] for r in vals], dtype=float
    )
    mzs = np.array(
        [r[f"track{t}_mz"] for r in vals], dtype=float
    )
    summary.append({
        "metric": f"track{t}_mass_error",
        "n": len(vals),
        "median": float(np.median(errors)) if len(errors) else "",
        "mad": float(np.median(np.abs(errors - np.median(errors))))
                if len(errors) else "",
        "min": float(np.min(errors)) if len(errors) else "",
        "max": float(np.max(errors)) if len(errors) else "",
        "correlation_with_track2_mz": (
            float(np.corrcoef(
                mzs,
                np.array([
                    r["track2_mz"] for r in vals
                ], dtype=float)
            )[0, 1])
            if len(vals) >= 3 and t != 2 else ""
        ),
    })

core_complete = [r for r in scan_rows if r["core_complete"]]
for a, b in zip(CORE[:-1], CORE[1:]):
    residuals = np.array([
        r[f"spacing_{a}_{b}_residual_mDa"]
        for r in core_complete
        if r[f"spacing_{a}_{b}_residual_mDa"] != ""
    ], dtype=float)
    summary.append({
        "metric": f"spacing_{a}_{b}_residual",
        "n": len(residuals),
        "median": float(np.median(residuals)) if len(residuals) else "",
        "mad": float(np.median(
            np.abs(residuals - np.median(residuals))
        )) if len(residuals) else "",
        "min": float(np.min(residuals)) if len(residuals) else "",
        "max": float(np.max(residuals)) if len(residuals) else "",
        "correlation_with_track2_mz": "",
    })

summary_fields = [
    "metric", "n", "median", "mad", "min", "max",
    "correlation_with_track2_mz",
]
with OUT_SUMMARY.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

print("Input:", INPUT)
print("Scans examined:", len(scan_rows))
print("Core tracks:", CORE)
print("Core-complete scans:", len(core_complete))
print("\n=== Core envelope ===")
print("Track 2–5 all present:", len(core_complete), "/", len(scan_rows))

for a, b in zip(CORE[:-1], CORE[1:]):
    residuals = np.array([
        r[f"spacing_{a}_{b}_residual_mDa"]
        for r in core_complete
        if r[f"spacing_{a}_{b}_residual_mDa"] != ""
    ], dtype=float)
    if len(residuals):
        print(
            f"{a}-{b}: median residual={np.median(residuals):+.3f} mDa; "
            f"median absolute residual="
            f"{np.median(np.abs(residuals)):.3f} mDa; "
            f"within ±3 mDa={np.sum(np.abs(residuals) <= 3)}/{len(residuals)}"
        )

print("\nIntensity ordering (descriptive only):")
if core_complete:
    decreasing = sum(
        r["core_intensity_strictly_decreasing"] == 1
        for r in core_complete
    )
    print(
        "Strictly decreasing from track 2 through track 5:",
        decreasing, "/", len(core_complete)
    )

print("\nDetailed scan report:", OUT_SCAN)
print("Summary report:", OUT_SUMMARY)
print("No mzML files modified; no final IG assignments made.")
