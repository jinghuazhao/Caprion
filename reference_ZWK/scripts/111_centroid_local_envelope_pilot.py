#!/usr/bin/env python3
"""Centroid local MS1 spectra and inspect the candidate isotope envelope."""

from pathlib import Path
import csv
import numpy as np
import pyopenms as oms

ROOT = Path.cwd()
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
OUT = ROOT / "reports"
OUT.mkdir(exist_ok=True)

FEATURE_MZ = 1341.06775
CHARGE = 2
RT_MIN, RT_MAX = 750.0, 770.0
WINDOW_DA = 0.025
THEORETICAL_SPACING = 1.00335483507 / CHARGE

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

picker = oms.PeakPickerHiRes()
rows = []

for index, spectrum in enumerate(exp):
    if spectrum.getMSLevel() != 1:
        continue
    rt = float(spectrum.getRT())
    if not RT_MIN <= rt <= RT_MAX:
        continue

    centroid = oms.MSSpectrum()
    picker.pick(spectrum, centroid)
    mz, intensity = centroid.get_peaks()
    mz = np.asarray(mz, dtype=float)
    intensity = np.asarray(intensity, dtype=float)

    # Candidate monoisotopic track and four successive isotopes.
    for isotope_index in range(5):
        target_mz = FEATURE_MZ + isotope_index * THEORETICAL_SPACING
        mask = np.abs(mz - target_mz) <= WINDOW_DA
        if not np.any(mask):
            rows.append({
                "spectrum_index": index,
                "rt": rt,
                "isotope_index": isotope_index,
                "target_mz": target_mz,
                "observed_mz": "",
                "mass_error_mda": "",
                "intensity": 0.0,
                "n_peaks_in_window": 0,
            })
            continue

        candidates = np.where(mask)[0]
        best = candidates[np.argmax(intensity[candidates])]
        rows.append({
            "spectrum_index": index,
            "rt": rt,
            "isotope_index": isotope_index,
            "target_mz": target_mz,
            "observed_mz": mz[best],
            "mass_error_mda": (mz[best] - target_mz) * 1000,
            "intensity": intensity[best],
            "n_peaks_in_window": len(candidates),
        })

detail = OUT / "feature10275997325121874585_centroid_pilot_v1.tsv"
with detail.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("=== Local centroided envelope pilot ===")
print("MS1 spectra analysed:", len({r["spectrum_index"] for r in rows}))
print("Detail:", detail)
print("Theoretical isotope spacing:", f"{THEORETICAL_SPACING:.9f} Da")
print("m/z window:", f"+/- {WINDOW_DA:.3f} Da")

print("\n=== Per-isotope detection and mass-error summary ===")
for isotope_index in range(5):
    subset = [r for r in rows if r["isotope_index"] == isotope_index]
    detected = [r for r in subset if r["observed_mz"] != ""]
    errors = [float(r["mass_error_mda"]) for r in detected]
    print(
        f"isotope {isotope_index}: detected={len(detected)}/{len(subset)}",
        f"median_error_mDa={np.median(errors):.3f}" if errors else "",
        f"median_intensity={np.median([float(r['intensity']) for r in detected]):.1f}"
        if detected else "",
    )

print("\n=== Scan-level envelope observations ===")
for spectrum_index in sorted({r["spectrum_index"] for r in rows}):
    subset = sorted(
        [r for r in rows if r["spectrum_index"] == spectrum_index],
        key=lambda r: r["isotope_index"],
    )
    detected = [r for r in subset if r["observed_mz"] != ""]
    intensities = [float(r["intensity"]) for r in detected]
    print(
        f"index={spectrum_index} RT={subset[0]['rt']:.5f}",
        f"detected={len(detected)}/5",
        "intensities=" + ",".join(f"{x:.0f}" for x in intensities),
    )

print("\nInterpretation caution:")
print("- This is a fixed-m/z-window pilot, not a feature detector.")
print("- A peak near each theoretical isotope position does not prove common origin.")
print("- The nearest/highest peak in a window can switch between coeluting species.")
print("- Independent IG construction still requires chromatographic coherence and validation.")
