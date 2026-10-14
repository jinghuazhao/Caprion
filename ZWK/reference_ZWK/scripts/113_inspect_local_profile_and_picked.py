#!/usr/bin/env python3
"""Inspect native profile points and picked centroids around candidate isotope tracks."""

from pathlib import Path
import numpy as np
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
RT_TARGET = 761.13018
RT_TOL = 0.01
MZ_MIN, MZ_MAX = 1339.0, 1345.5

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

matches = [
    (i, s) for i, s in enumerate(exp)
    if s.getMSLevel() == 1 and abs(s.getRT() - RT_TARGET) <= RT_TOL
]
if len(matches) != 1:
    raise RuntimeError(f"Expected one spectrum near target RT; found {len(matches)}")

index, source = matches[0]
native_mz, native_i = source.get_peaks()
native_mz = np.asarray(native_mz, dtype=float)
native_i = np.asarray(native_i, dtype=float)

picker = oms.PeakPickerHiRes()
picked = oms.MSSpectrum()
picker.pick(source, picked)
picked_mz, picked_i = picked.get_peaks()
picked_mz = np.asarray(picked_mz, dtype=float)
picked_i = np.asarray(picked_i, dtype=float)

print(f"Spectrum index={index}, RT={source.getRT():.5f} s")
print(f"Native points: {len(native_mz)}; picked peaks: {len(picked_mz)}")
print(f"Inspecting m/z {MZ_MIN}–{MZ_MAX}")

print("\n=== Native profile points in the region ===")
mask = (native_mz >= MZ_MIN) & (native_mz <= MZ_MAX)
for mz, inten in zip(native_mz[mask], native_i[mask]):
    print(f"{mz:.7f}\t{inten:.1f}")

print("\n=== Picked centroid peaks in the region ===")
mask = (picked_mz >= MZ_MIN) & (picked_mz <= MZ_MAX)
for mz, inten in zip(picked_mz[mask], picked_i[mask]):
    print(f"{mz:.7f}\t{inten:.1f}")

print("\n=== Nearest picked peak to each native track position ===")
track_targets = {
    2: 1341.069458,
    3: 1341.571533,
    4: 1342.071533,
    5: 1342.573364,
    6: 1343.074707,
    7: 1343.576782,
}
for track, target in track_targets.items():
    if len(picked_mz) == 0:
        print(f"track {track}: no picked peaks")
        continue
    j = int(np.argmin(np.abs(picked_mz - target)))
    print(
        f"track {track}: native_target={target:.7f}, "
        f"nearest_picked={picked_mz[j]:.7f}, "
        f"delta_mDa={(picked_mz[j] - target) * 1000:.3f}, "
        f"intensity={picked_i[j]:.1f}"
    )

print("\nNo files were modified.")
