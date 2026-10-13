#!/usr/bin/env python3
import re
import numpy as np
import pandas as pd
import pyopenms as oms
from scipy.spatial import cKDTree

BASE = "search_98/szwk900404i19201"
FEATURE_XML = "openms_pilot/szwk900404i19201xms1.featureXML"
OUT = "openms_pilot/psm_feature_matches.tsv"

PPM = 10.0
RT_TOL = 10.0  # seconds

psm = pd.read_csv(f"{BASE}/psm.tsv", sep="\t")
psm["scan"] = pd.to_numeric(
    psm["Spectrum"].astype(str).str.extract(r"\.(\d+)\.\d+\.\d+$")[0],
    errors="coerce"
)

# Read independently detected MS1 features.
fm = oms.FeatureMap()
oms.FeatureXMLFile().load(FEATURE_XML, fm)

features = []
for f in fm:
    features.append({
        "feature_id": str(f.getUniqueId()),
        "feature_mz": float(f.getMZ()),
        "feature_rt_sec": float(f.getRT()),
        "feature_charge": int(f.getCharge()),
        "feature_intensity": float(f.getIntensity()),
        "feature_quality": float(f.getOverallQuality()),
    })

feat = pd.DataFrame(features)
if feat.empty:
    raise SystemExit("No features found in featureXML.")

# Index by charge and retention time, so only nearby features are examined.
indices = {}
for charge, group in feat.groupby("feature_charge"):
    arr = group[["feature_rt_sec", "feature_mz"]].to_numpy()
    # RT and m/z are queried separately below; RT index narrows candidates.
    indices[int(charge)] = (group.index.to_numpy(), cKDTree(arr[:, [0]]))

rows = []
for _, p in psm.iterrows():
    mz = pd.to_numeric(p.get("Calibrated Observed M/Z"), errors="coerce")
    if pd.isna(mz):
        mz = pd.to_numeric(p.get("Observed M/Z"), errors="coerce")
    rt = pd.to_numeric(p.get("Retention"), errors="coerce")
    charge = pd.to_numeric(p.get("Charge"), errors="coerce")

    if pd.isna(mz) or pd.isna(rt) or pd.isna(charge):
        continue

    charge = int(charge)
    if charge not in indices:
        continue

    global_idx, tree = indices[charge]
    local = tree.query_ball_point([float(rt)], r=RT_TOL)
    if not local:
        rows.append({
            "Spectrum": p["Spectrum"],
            "scan": p["scan"],
            "Peptide": p["Peptide"],
            "Probability": p["Probability"],
            "psm_mz": mz,
            "psm_rt_sec": rt,
            "charge": charge,
            "n_matches": 0,
            "feature_id": "",
            "feature_mz": np.nan,
            "feature_rt_sec": np.nan,
            "ppm_error": np.nan,
            "rt_error_sec": np.nan,
            "feature_intensity": np.nan,
            "feature_quality": np.nan,
        })
        continue

    candidates = feat.loc[global_idx[np.asarray(local)]].copy()
    candidates["ppm_error"] = (
        (candidates["feature_mz"] - float(mz)) / float(mz) * 1e6
    )
    candidates["rt_error_sec"] = (
        candidates["feature_rt_sec"] - float(rt)
    )
    candidates = candidates[candidates["ppm_error"].abs() <= PPM]

    if candidates.empty:
        rows.append({
            "Spectrum": p["Spectrum"],
            "scan": p["scan"],
            "Peptide": p["Peptide"],
            "Probability": p["Probability"],
            "psm_mz": mz,
            "psm_rt_sec": rt,
            "charge": charge,
            "n_matches": 0,
            "feature_id": "",
            "feature_mz": np.nan,
            "feature_rt_sec": np.nan,
            "ppm_error": np.nan,
            "rt_error_sec": np.nan,
            "feature_intensity": np.nan,
            "feature_quality": np.nan,
        })
        continue

    # Emit one row per candidate; do not hide ambiguity.
    n = len(candidates)
    for _, f in candidates.iterrows():
        rows.append({
            "Spectrum": p["Spectrum"],
            "scan": p["scan"],
            "Peptide": p["Peptide"],
            "Probability": p["Probability"],
            "psm_mz": mz,
            "psm_rt_sec": rt,
            "charge": charge,
            "n_matches": n,
            "feature_id": f["feature_id"],
            "feature_mz": f["feature_mz"],
            "feature_rt_sec": f["feature_rt_sec"],
            "ppm_error": f["ppm_error"],
            "rt_error_sec": f["rt_error_sec"],
            "feature_intensity": f["feature_intensity"],
            "feature_quality": f["feature_quality"],
        })

out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False)

per_psm = out.drop_duplicates("Spectrum")
print("PSMs:", len(psm))
print("PSMs with >=1 candidate feature:",
      (per_psm["n_matches"] > 0).sum())
print("PSMs with exactly 1 candidate:",
      (per_psm["n_matches"] == 1).sum())
print("PSMs with multiple candidates:",
      (per_psm["n_matches"] > 1).sum())
print("PSMs without a candidate:",
      (per_psm["n_matches"] == 0).sum())
print("Output:", OUT)
