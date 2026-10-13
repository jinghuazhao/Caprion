#!/usr/bin/env python3
"""Spatially indexed audit of MS1 peaks falling inside OpenMS feature hulls."""

from pathlib import Path
from collections import defaultdict
import bisect

import numpy as np
import pandas as pd
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
FEATUREXML = Path("openms_pilot/szwk900404i19201xms1.featureXML")
OUT = Path("openms_pilot")
OUT.mkdir(parents=True, exist_ok=True)

# Load feature map.
fm = oms.FeatureMap()
oms.FeatureXMLFile().load(str(FEATUREXML), fm)

meta = []
hulls = []

for f in fm:
    fid = int(f.getUniqueId())
    meta.append({
        "feature_id": fid,
        "feature_mz": float(f.getMZ()),
        "feature_rt_sec": float(f.getRT()),
        "feature_intensity": float(f.getIntensity()),
        "feature_charge": int(f.getCharge()),
        "feature_quality": float(f.getOverallQuality()),
    })

    for hi, hull in enumerate(f.getConvexHulls()):
        pts = hull.getHullPoints()
        if len(pts) < 3:
            continue

        polygon = [(float(p[0]), float(p[1])) for p in pts]
        xs = [p[0] for p in polygon]
        ys = [p[1] for p in polygon]

        hulls.append({
            "feature_id": fid,
            "hull_index": hi,
            "polygon": polygon,
            "rt_min": min(xs),
            "rt_max": max(xs),
            "mz_min": min(ys),
            "mz_max": max(ys),
        })

feature_meta = pd.DataFrame(meta)

# Index hulls by RT start. At a given RT, inspect only hulls that have started.
hulls.sort(key=lambda h: h["rt_min"])
rt_starts = [h["rt_min"] for h in hulls]

# Ray-casting point-in-polygon test.
def inside(x, y, polygon):
    inside_flag = False
    j = len(polygon) - 1

    for i, (xi, yi) in enumerate(polygon):
        xj, yj = polygon[j]

        if (yi > y) != (yj > y):
            crossing_x = xi + (y - yi) * (xj - xi) / (yj - yi)
            if x < crossing_x:
                inside_flag = not inside_flag

        j = i

    return inside_flag

# Load mzML.
exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

# Store compact records: spectrum index, within-spectrum peak index, RT,
# observed m/z, intensity, feature ID, hull index.
records = []
ms1_count = 0
peak_count = 0

for scan_index, spec in enumerate(exp):
    if spec.getMSLevel() != 1:
        continue

    ms1_count += 1
    rt = float(spec.getRT())
    mzs, intensities = spec.get_peaks()
    peak_count += len(mzs)

    # Only hulls whose RT range includes this spectrum can match.
    end = bisect.bisect_right(rt_starts, rt)
    active = [
        h for h in hulls[:end]
        if h["rt_max"] >= rt
    ]

    # Quickly reject peaks outside all active hull m/z ranges.
    if not active:
        continue

    mz_low = min(h["mz_min"] for h in active)
    mz_high = max(h["mz_max"] for h in active)
    mzs = np.asarray(mzs)
    intensities = np.asarray(intensities)

    candidate_indices = np.flatnonzero(
        (mzs >= mz_low) & (mzs <= mz_high)
    )

    # Further filter hulls by m/z bounds for each candidate peak.
    for pi in candidate_indices:
        mz = float(mzs[pi])
        intensity = float(intensities[pi])

        for h in active:
            if not (h["mz_min"] <= mz <= h["mz_max"]):
                continue

            if inside(rt, mz, h["polygon"]):
                records.append((
                    scan_index, int(pi), rt, mz, intensity,
                    h["feature_id"], h["hull_index"]
                ))

    if ms1_count % 250 == 0:
        print(
            f"MS1 spectra: {ms1_count:,}; centroid peaks: {peak_count:,}; "
            f"hull assignments: {len(records):,}",
            flush=True
        )

cols = [
    "scan_index", "peak_index", "rt_sec", "observed_mz",
    "peak_intensity", "feature_id", "hull_index"
]
assign = pd.DataFrame(records, columns=cols)

if assign.empty:
    raise SystemExit("No hull assignments found; check the input files.")

# Identify each centroid by its spectrum and within-spectrum peak index.
peak_counts = (
    assign.groupby(["scan_index", "peak_index"])["feature_id"]
    .nunique()
    .rename("n_features")
    .reset_index()
)
assign = assign.merge(
    peak_counts, on=["scan_index", "peak_index"], validate="many_to_one"
)
assign["shared_peak"] = assign["n_features"] > 1
assign["shared_intensity"] = assign["peak_intensity"].where(
    assign["shared_peak"], 0.0
)

assign.to_csv(
    OUT / "feature_overlap_assignments.tsv", sep="\t", index=False
)

# Feature-level audit.
summary = (
    assign.groupby("feature_id", as_index=False)
    .agg(
        assigned_peak_count=("peak_index", "size"),
        unique_peak_count=("peak_index", "nunique"),
        assigned_intensity=("peak_intensity", "sum"),
        shared_assignment_count=("shared_peak", "sum"),
        shared_intensity=("shared_intensity", "sum"),
    )
)

summary["shared_intensity_fraction"] = (
    summary["shared_intensity"] /
    summary["assigned_intensity"].replace(0, np.nan)
)

summary = feature_meta.merge(
    summary, on="feature_id", how="left", validate="one_to_one"
)

for c in [
    "assigned_peak_count", "unique_peak_count",
    "shared_assignment_count", "assigned_intensity", "shared_intensity"
]:
    summary[c] = summary[c].fillna(0)

summary["has_shared_peaks"] = summary["shared_assignment_count"] > 0
summary.to_csv(
    OUT / "feature_overlap_qc_geometry.tsv", sep="\t", index=False
)

print("\nFinished")
print("Features:", len(summary))
print("MS1 spectra:", ms1_count)
print("Centroid peaks:", peak_count)
print("Hull assignments:", len(assign))
print("Features with shared peaks:", int(summary["has_shared_peaks"].sum()))
print("Saved:", OUT / "feature_overlap_qc_geometry.tsv")
