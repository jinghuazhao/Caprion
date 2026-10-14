#!/usr/bin/env python3
"""Validate spectrum-type interpretation and native centroid spacing."""

from pathlib import Path
import numpy as np
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
RT_TARGET = 761.13018
RT_TOL = 0.01
MZ_MIN, MZ_MAX = 1340.8, 1344.0
THEORETICAL_SPACING = 1.00335483507 / 2

print("PyOpenMS version:", oms.__version__)

print("\n=== Spectrum type enum ===")
for name in ("SpectrumSettings", "MSSpectrum"):
    cls = getattr(oms, name, None)
    if cls is not None:
        print(name, "available")
        for attr in dir(cls):
            if "TYPE" in attr.upper() or "CENTROID" in attr.upper() or "PROFILE" in attr.upper():
                try:
                    print(f"  {attr} = {getattr(cls, attr)}")
                except Exception:
                    pass

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

matches = [
    (i, s) for i, s in enumerate(exp)
    if s.getMSLevel() == 1 and abs(s.getRT() - RT_TARGET) <= RT_TOL
]
if len(matches) != 1:
    raise RuntimeError(f"Expected one spectrum near target RT; found {len(matches)}")

index, spectrum = matches[0]
mz, intensity = spectrum.get_peaks()
mz = np.asarray(mz, dtype=float)
intensity = np.asarray(intensity, dtype=float)

print("\n=== Target spectrum ===")
print("Index:", index)
print("RT:", spectrum.getRT())
print("MS level:", spectrum.getMSLevel())
print("Spectrum type enum:", spectrum.getType())
print("Number of m/z values:", len(mz))
print("Number of intensities:", len(intensity))
print("m/z array sorted:", bool(np.all(np.diff(mz) >= 0)))
print("Unique m/z values:", len(np.unique(mz)))

mask = (mz >= MZ_MIN) & (mz <= MZ_MAX)
local_mz = mz[mask]
local_i = intensity[mask]

print(f"\n=== Native peaks in {MZ_MIN}–{MZ_MAX} m/z ===")
for x, y in zip(local_mz, local_i):
    print(f"{x:.7f}\t{y:.1f}")

print("\n=== Adjacent local native-peak spacings ===")
if len(local_mz) > 1:
    for j in range(len(local_mz) - 1):
        spacing = local_mz[j + 1] - local_mz[j]
        residual_mda = (spacing - THEORETICAL_SPACING) * 1000
        print(
            f"{local_mz[j]:.7f} -> {local_mz[j + 1]:.7f}: "
            f"{spacing:.7f} Da; theoretical residual={residual_mda:+.3f} mDa"
        )
else:
    print("Fewer than two peaks in the selected interval.")

print("\n=== Workflow guard ===")
print(
    "Do not apply PeakPickerHiRes to this input until the source spectrum "
    "type has been verified from the mzML metadata and the intended "
    "processing strategy is documented."
)
print("No files were modified.")
