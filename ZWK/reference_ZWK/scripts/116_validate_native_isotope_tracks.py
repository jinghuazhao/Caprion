#!/usr/bin/env python3
"""Validate candidate isotope tracks using native centroid MS1 peaks."""

from pathlib import Path
import csv
import numpy as np
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
OUT = Path("reports/feature10275997325121874585_native_track_validation_v1.tsv")

RT_MIN, RT_MAX = 750.0, 770.0
MZ_TOL = 0.015  # Da; report this tolerance explicitly
THEORETICAL_SPACING = 1.00335483507 / 2

TRACKS = {
    2: 1341.067749,
    3: 1341.570496,
    4: 1342.068359,
    5: 1342.571045,
    6: 1343.074463,
    7: 1343.573608,
}

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

spectra = []
for i, spectrum in enumerate(exp):
    if spectrum.getMSLevel() != 1:
        continue
    if not (RT_MIN <= spectrum.getRT() <= RT_MAX):
        continue
    if spectrum.getType() != oms.SpectrumSettings.SpectrumType.CENTROID:
        raise RuntimeError(
            f"Spectrum {i} at RT {spectrum.getRT():.5f} "
            f"is not marked CENTROID (type={spectrum.getType()})"
        )
    mz, intensity = spectrum.get_peaks()
    mz = np.asarray(mz, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    order = np.argsort(mz)
    spectra.append((i, float(spectrum.getRT()), mz[order], intensity[order]))

print("File:", MZML)
print("MS1 scans in RT window:", len(spectra))
print("Spectrum type: CENTROID")
print("Matching tolerance:", MZ_TOL, "Da")
print("Theoretical z=2 spacing:", f"{THEORETICAL_SPACING:.9f}", "Da")

if not spectra:
    raise RuntimeError("No MS1 spectra found in the selected RT window")

rows = []
for track_id, target_mz in TRACKS.items():
    for scan_index, rt, mz, intensity in spectra:
        lo = np.searchsorted(mz, target_mz - MZ_TOL, side="left")
        hi = np.searchsorted(mz, target_mz + MZ_TOL, side="right")
        if lo == hi:
            rows.append({
                "track": track_id, "scan_index": scan_index, "rt": rt,
                "target_mz": target_mz, "observed_mz": "",
                "intensity": "", "error_mDa": "", "n_candidates": 0,
            })
            continue

        # Choose the candidate closest in m/z; retain the number of candidates
        # so ambiguous windows remain visible in the report.
        candidates = np.arange(lo, hi)
        best = candidates[np.argmin(np.abs(mz[candidates] - target_mz))]
        rows.append({
            "track": track_id, "scan_index": scan_index, "rt": rt,
            "target_mz": target_mz,
            "observed_mz": float(mz[best]),
            "intensity": float(intensity[best]),
            "error_mDa": float((mz[best] - target_mz) * 1000),
            "n_candidates": int(hi - lo),
        })

OUT.parent.mkdir(parents=True, exist_ok=True)
fields = [
    "track", "scan_index", "rt", "target_mz", "observed_mz",
    "intensity", "error_mDa", "n_candidates",
]
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("\n=== Track summaries ===")
for track_id in TRACKS:
    found = [
        r for r in rows
        if r["track"] == track_id and r["observed_mz"] != ""
    ]
    errors = np.array([r["error_mDa"] for r in found], dtype=float)
    intensities = np.array([r["intensity"] for r in found], dtype=float)
    ambiguous = sum(r["n_candidates"] > 1 for r in found)
    print(
        f"Track {track_id}: detected {len(found)}/{len(spectra)} scans; "
        f"median error={np.median(errors):+.3f} mDa; "
        f"median intensity={np.median(intensities):.1f}; "
        f"ambiguous windows={ambiguous}"
        if found else
        f"Track {track_id}: detected 0/{len(spectra)} scans"
    )

print("\n=== Adjacent track spacing ===")
for lower, upper in zip(TRACKS, list(TRACKS)[1:]):
    low_by_scan = {
        r["scan_index"]: r for r in rows
        if r["track"] == lower and r["observed_mz"] != ""
    }
    high_by_scan = {
        r["scan_index"]: r for r in rows
        if r["track"] == upper and r["observed_mz"] != ""
    }
    shared = sorted(set(low_by_scan) & set(high_by_scan))
    spacings = np.array([
        high_by_scan[s]["observed_mz"] - low_by_scan[s]["observed_mz"]
        for s in shared
    ])
    if len(spacings):
        residuals = (spacings - THEORETICAL_SPACING) * 1000
        print(
            f"{lower}-{upper}: shared={len(spacings)}; "
            f"median spacing={np.median(spacings):.7f} Da; "
            f"median |residual|={np.median(np.abs(residuals)):.3f} mDa; "
            f"within 3 mDa={np.sum(np.abs(residuals) <= 3)}/{len(residuals)}"
        )
    else:
        print(f"{lower}-{upper}: no shared detections")

print("\nDetailed table:", OUT)
print("No peak picking performed; no mzML files modified.")
