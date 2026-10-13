#!/usr/bin/env python3
"""Audit MS1 centroid-peak assignments to overlapping OpenMS feature hulls.

This is an audit, not an intensity-correction procedure. FeatureXML hulls
define geometric regions; peaks inside multiple regions are flagged as
ambiguous. Original OpenMS feature intensities are not modified.
"""

from pathlib import Path
from collections import defaultdict

import pandas as pd
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
FEATUREXML = Path("openms_pilot/szwk900404i19201xms1.featureXML")
OUTDIR = Path("openms_pilot")
OUTDIR.mkdir(parents=True, exist_ok=True)

# -------------------- Load feature map --------------------

fm = oms.FeatureMap()
oms.FeatureXMLFile().load(str(FEATUREXML), fm)

features = []
hulls_by_feature = {}

for f in fm:
    fid = int(f.getUniqueId())
    hulls = []

    for hull in f.getConvexHulls():
        points = hull.getHullPoints()
        if len(points) < 3:
            continue

        # OpenMS hull coordinates: x=RT (seconds), y=m/z.
        polygon = [(float(p[0]), float(p[1])) for p in points]
        xs = [p[0] for p in polygon]
        ys = [p[1] for p in polygon]

        hulls.append({
            "polygon": polygon,
            "rt_min": min(xs),
            "rt_max": max(xs),
            "mz_min": min(ys),
            "mz_max": max(ys),
        })

    hulls_by_feature[fid] = hulls
    features.append({
        "feature_id": fid,
        "feature_mz": float(f.getMZ()),
        "feature_rt_sec": float(f.getRT()),
        "feature_intensity": float(f.getIntensity()),
        "feature_charge": int(f.getCharge()),
        "feature_quality": float(f.getOverallQuality()),
        "n_hulls": len(hulls),
    })

feature_meta = pd.DataFrame(features)
print("Features loaded:", len(feature_meta))
print("Features with hulls:", sum(bool(v) for v in hulls_by_feature.values()))

# -------------------- Polygon containment --------------------

def point_in_polygon(x, y, polygon):
    """Ray-casting containment test; boundary points may be ambiguous."""
    inside = False
    n = len(polygon)
    j = n - 1

    for i in range(n):
        xi, yi = polygon[i]
        xj, yj = polygon[j]

        # Check whether the horizontal ray crosses this polygon edge.
        if (yi > y) != (yj > y):
            x_cross = xi + (y - yi) * (xj - xi) / (yj - yi)
            if x < x_cross:
                inside = not inside

        j = i

    return inside


# -------------------- Read mzML and test peaks --------------------

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

# Index hulls by their RT bounds so each spectrum only checks nearby hulls.
hulls = []
for fid, fhulls in hulls_by_feature.items():
    for hidx, h in enumerate(fhulls):
        hulls.append({
            "feature_id": fid,
            "hull_index": hidx,
            **h,
        })

hulls.sort(key=lambda h: h["rt_min"])
hull_rt_mins = [h["rt_min"] for h in hulls]

import bisect

assignments = []
ms1_scans = 0
ms1_peaks = 0

for scan_index, spec in enumerate(exp):
    if spec.getMSLevel() != 1:
        continue

    ms1_scans += 1
    rt = float(spec.getRT())
    mzs, intensities = spec.get_peaks()
    ms1_peaks += len(mzs)

    # Hulls beginning after this RT cannot contain a peak from this scan.
    end = bisect.bisect_right(hull_rt_mins, rt)

    for peak_index, (mz, intensity) in enumerate(zip(mzs, intensities)):
        mz = float(mz)
        intensity = float(intensity)

        matching = []
        for h in hulls[:end]:
            if rt > h["rt_max"]:
                continue
            if not (h["mz_min"] <= mz <= h["mz_max"]):
                continue

            if point_in_polygon(rt, mz, h["polygon"]):
                matching.append(h)

        if not matching:
            continue

        feature_ids = sorted({h["feature_id"] for h in matching})

        for h in matching:
            assignments.append({
                "scan_index": scan_index,
                "peak_index": peak_index,
                "rt_sec": rt,
                "observed_mz": mz,
                "peak_intensity": intensity,
                "feature_id": h["feature_id"],
                "hull_index": h["hull_index"],
                "n_matching_features": len(feature_ids),
                "shared_across_features": len(feature_ids) > 1,
                "matching_feature_ids": ",".join(map(str, feature_ids)),
            })

    if ms1_scans % 500 == 0:
        print(f"Processed {ms1_scans} MS1 spectra; "
              f"{ms1_peaks:,} centroid peaks seen")

# -------------------- Save assignment-level audit --------------------

assign_df = pd.DataFrame(assignments)
assign_path = OUTDIR / "feature_overlap_assignments.tsv"
assign_df.to_csv(assign_path, sep="\t", index=False)

print("\nMS1 spectra:", ms1_scans)
print("MS1 centroid peaks:", f"{ms1_peaks:,}")
print("Hull assignments:", f"{len(assign_df):,}")

if assign_df.empty:
    print("No peaks were assigned to hulls; inspect input files and hull geometry.")
    raise SystemExit(0)

# A physical centroid peak is identified within this mzML by its spectrum
# index and within-spectrum peak index, rather than floating-point values.
peak_key = ["scan_index", "peak_index"]

peak_membership = (
    assign_df.groupby(peak_key, as_index=False)
    .agg(
        n_features=("feature_id", "nunique"),
        matching_feature_ids=("feature_id",
                              lambda s: ",".join(map(str, sorted(set(s)))),
        ),
        peak_intensity=("peak_intensity", "first"),
        observed_mz=("observed_mz", "first"),
        rt_sec=("rt_sec", "first"),
    )
)

shared_keys = peak_membership[peak_membership["n_features"] > 1]
shared_peak_set = set(
    zip(shared_keys["scan_index"], shared_keys["peak_index"])
)

assign_df["shared_peak"] = [
    (s, p) in shared_peak_set
    for s, p in zip(assign_df["scan_index"], assign_df["peak_index"])
]
assign_df["shared_peak_intensity"] = assign_df["peak_intensity"].where(
    assign_df["shared_peak"], 0.0
)

# -------------------- Feature-level summary --------------------

summary = (
    assign_df.groupby("feature_id", as_index=False)
    .agg(
        assigned_peak_count=("peak_index", "size"),
        unique_peak_count=("peak_index", "nunique"),
        assigned_intensity=("peak_intensity", "sum"),
        shared_assignment_count=("shared_peak", "sum"),
        shared_intensity=("shared_peak_intensity", "sum"),
        maximum_matching_features=("n_matching_features", "max"),
    )
)

summary["shared_intensity_fraction"] = (
    summary["shared_intensity"] / summary["assigned_intensity"].replace(0, pd.NA)
)

summary = feature_meta.merge(
    summary, on="feature_id", how="left", validate="one_to_one"
)

count_cols = [
    "assigned_peak_count", "unique_peak_count", "shared_assignment_count",
    "maximum_matching_features",
]
for col in count_cols:
    summary[col] = summary[col].fillna(0).astype(int)

for col in ["assigned_intensity", "shared_intensity"]:
    summary[col] = summary[col].fillna(0.0)

summary["has_shared_peaks"] = summary["shared_assignment_count"] > 0

summary_path = OUTDIR / "feature_overlap_qc_geometry.tsv"
summary.sort_values(
    "shared_intensity_fraction", ascending=False, na_position="last"
).to_csv(summary_path, sep="\t", index=False)

shared_path = OUTDIR / "shared_centroid_peaks_geometry.tsv"
shared_keys.to_csv(shared_path, sep="\t", index=False)

# -------------------- Report --------------------

print("\nFeature-level overlap summary")
print("-----------------------------")
print("Features:", len(summary))
print("Features with shared peaks:", int(summary["has_shared_peaks"].sum()))
print("Shared physical centroid peaks:", len(shared_keys))
print("Shared peak assignments:", int(assign_df["shared_peak"].sum()))

print("\nShared-intensity fraction:")
print(
    summary["shared_intensity_fraction"]
    .describe(percentiles=[0.5, 0.75, 0.9, 0.95, 0.99])
    .to_string()
)

print("\nTop 10 features by shared-intensity fraction:")
print(
    summary.sort_values("shared_intensity_fraction", ascending=False)
    .head(10)[[
        "feature_id", "feature_mz", "feature_rt_sec", "feature_charge",
        "assigned_peak_count", "shared_assignment_count",
        "shared_intensity_fraction",
    ]]
    .to_string(index=False)
)

print("\nSaved:")
print(assign_path)
print(summary_path)
print(shared_path)
