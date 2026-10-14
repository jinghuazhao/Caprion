#!/usr/bin/env python3
"""Survey spectrum types and local MS1 sampling in one ZWK mzML."""

from collections import Counter
from pathlib import Path
import pyopenms as oms

mzml = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)

if not mzml.is_file():
    raise FileNotFoundError(mzml)

print("Input:", mzml)
print("File size (GB):", round(mzml.stat().st_size / 1e9, 3))
print("Loading mzML metadata and spectra...")

exp = oms.MSExperiment()
oms.MzMLFile().load(str(mzml), exp)

print("\n=== Spectrum counts ===")
levels = Counter(s.getMSLevel() for s in exp)
print("Total spectra:", exp.size())
for level, count in sorted(levels.items()):
    print(f"MS{level}: {count}")

print("\n=== Spectrum data-point counts by MS level ===")
for level in sorted(levels):
    sizes = [
        len(s.get_peaks()[0])
        for s in exp if s.getMSLevel() == level
    ]
    if sizes:
        sizes_sorted = sorted(sizes)
        n = len(sizes_sorted)
        print(
            f"MS{level}: min={sizes_sorted[0]}, "
            f"median={sizes_sorted[n // 2]}, "
            f"max={sizes_sorted[-1]}"
        )

print("\n=== Local MS1 spectra, RT 750–770 s ===")
local = [
    (i, s) for i, s in enumerate(exp)
    if s.getMSLevel() == 1 and 750 <= s.getRT() <= 770
]
print("Number of local MS1 spectra:", len(local))

for i, s in local:
    mz, intensity = s.get_peaks()
    print(
        f"index={i:6d} RT={s.getRT():10.5f} "
        f"points={len(mz):7d} "
        f"type={s.getType()}"
    )

print("\n=== PeakPickerHiRes test on one local MS1 spectrum ===")
if local:
    index, source = local[len(local) // 2]
    print(
        f"Testing spectrum index={index}, RT={source.getRT():.5f} s, "
        f"input points={len(source.get_peaks()[0])}"
    )
    picker = oms.PeakPickerHiRes()
    output = oms.MSSpectrum()

    try:
        picker.pick(source, output)
        out_mz, out_intensity = output.get_peaks()
        print("Centroid peaks:", len(out_mz))
        print("Input points:", len(source.get_peaks()[0]))
        if len(out_mz):
            print("First five centroid m/z values:", out_mz[:5].tolist())
            print("First five centroid intensities:", out_intensity[:5].tolist())
    except Exception as exc:
        print("Peak-picking test failed:", repr(exc))
else:
    print("No MS1 spectra found in the requested RT interval.")

print("\nSurvey complete; source mzML was not modified.")
