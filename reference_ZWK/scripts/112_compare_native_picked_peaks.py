#!/usr/bin/env python3
"""Compare native centroid and PeakPickerHiRes output around candidate tracks."""

from pathlib import Path
import csv
import numpy as np
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
OUT = Path("reports")
OUT.mkdir(exist_ok=True)

# Median m/z values from the earlier native-centroid track audit.
TRACK_MZ = {
    2: 1341.067749,
    3: 1341.570496,
    4: 1342.068359,
    5: 1342.571045,
    6: 1343.074463,
    7: 1343.573608,
}
RT_MIN, RT_MAX = 750.0, 770.0
WINDOW_DA = 0.030

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

    native_mz, native_i = spectrum.get_peaks()
    native_mz = np.asarray(native_mz, dtype=float)
    native_i = np.asarray(native_i, dtype=float)

    picked = oms.MSSpectrum()
    picker.pick(spectrum, picked)
    picked_mz, picked_i = picked.get_peaks()
    picked_mz = np.asarray(picked_mz, dtype=float)
    picked_i = np.asarray(picked_i, dtype=float)

    for track, target in TRACK_MZ.items():
        for source_name, mz_array, i_array in (
            ("native", native_mz, native_i),
            ("picked", picked_mz, picked_i),
        ):
            mask = np.abs(mz_array - target) <= WINDOW_DA
            candidates = np.where(mask)[0]
            if len(candidates):
                best = candidates[np.argmax(i_array[candidates])]
                observed_mz = float(mz_array[best])
                intensity = float(i_array[best])
                error_mda = (observed_mz - target) * 1000
            else:
                observed_mz = ""
                intensity = 0.0
                error_mda = ""
            rows.append({
                "spectrum_index": index,
                "rt": rt,
                "track": track,
                "target_mz": target,
                "source": source_name,
                "observed_mz": observed_mz,
                "error_mda": error_mda,
                "intensity": intensity,
                "n_peaks_in_window": len(candidates),
            })

outfile = OUT / "feature10275997325121874585_native_vs_picked_v1.tsv"
with outfile.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("Output:", outfile)
print("Local MS1 spectra:", len({r["spectrum_index"] for r in rows}))
print("Target tracks:", len(TRACK_MZ))
print("Window: +/-", WINDOW_DA, "Da")

print("\n=== Detection by track and source ===")
for track in TRACK_MZ:
    for source in ("native", "picked"):
        subset = [
            r for r in rows
            if r["track"] == track and r["source"] == source
        ]
        found = [r for r in subset if r["observed_mz"] != ""]
        errors = [float(r["error_mda"]) for r in found]
        print(
            f"track={track} source={source}: "
            f"{len(found)}/{len(subset)} spectra",
            f"median_error_mDa={np.median(errors):.3f}" if errors else "",
            f"median_intensity={np.median([r['intensity'] for r in found]):.0f}"
            if found else "",
        )

print("\n=== Native versus picked m/z at each scan and track ===")
for index in sorted({r["spectrum_index"] for r in rows}):
    for track in TRACK_MZ:
        subset = [
            r for r in rows
            if r["spectrum_index"] == index and r["track"] == track
        ]
        native = next(r for r in subset if r["source"] == "native")
        picked = next(r for r in subset if r["source"] == "picked")
        if native["observed_mz"] != "" or picked["observed_mz"] != "":
            print(
                f"scan={index} RT={native['rt']:.5f} track={track} "
                f"native={native['observed_mz']} "
                f"picked={picked['observed_mz']}"
            )

print("\nThis is a diagnostic comparison, not an isotope-group assignment.")
