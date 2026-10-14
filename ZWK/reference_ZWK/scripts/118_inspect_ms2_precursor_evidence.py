#!/usr/bin/env python3
"""Inspect MS2 precursor metadata near a candidate MS1 feature."""

from pathlib import Path
import csv
import numpy as np
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
OUT = Path(
    "reports/feature10275997325121874585_ms2_precursor_evidence_v1.tsv"
)

RT_TARGET = 761.13018
RT_WINDOW = 10.0  # seconds on either side of the candidate RT
MZ_TOL = 0.015

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

rows = []

for scan_index, spectrum in enumerate(exp):
    if spectrum.getMSLevel() != 2:
        continue

    rt = float(spectrum.getRT())
    if abs(rt - RT_TARGET) > RT_WINDOW:
        continue

    precursors = spectrum.getPrecursors()
    if not precursors:
        rows.append({
            "scan_index": scan_index,
            "rt": rt,
            "rt_delta_s": rt - RT_TARGET,
            "precursor_number": "",
            "precursor_mz": "",
            "precursor_charge": "",
            "isolation_lower_mz": "",
            "isolation_upper_mz": "",
            "nearest_track": "",
            "nearest_track_error_mDa": "",
            "within_track_tolerance": "",
        })
        continue

    for pnum, precursor in enumerate(precursors, start=1):
        pmz = float(precursor.getMZ())
        charge = int(precursor.getCharge())

        # PyOpenMS reports isolation offsets relative to precursor m/z.
        lower = pmz - float(precursor.getIsolationWindowLowerOffset())
        upper = pmz + float(precursor.getIsolationWindowUpperOffset())

        nearest_track = min(
            TRACKS,
            key=lambda t: abs(TRACKS[t] - pmz)
        )
        error_mda = (pmz - TRACKS[nearest_track]) * 1000

        rows.append({
            "scan_index": scan_index,
            "rt": rt,
            "rt_delta_s": rt - RT_TARGET,
            "precursor_number": pnum,
            "precursor_mz": pmz,
            "precursor_charge": charge,
            "isolation_lower_mz": lower,
            "isolation_upper_mz": upper,
            "nearest_track": nearest_track,
            "nearest_track_error_mDa": error_mda,
            "within_track_tolerance": int(abs(error_mda) <= MZ_TOL * 1000),
        })

fields = [
    "scan_index", "rt", "rt_delta_s", "precursor_number",
    "precursor_mz", "precursor_charge",
    "isolation_lower_mz", "isolation_upper_mz",
    "nearest_track", "nearest_track_error_mDa",
    "within_track_tolerance",
]

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("Input:", MZML)
print("Target RT:", RT_TARGET)
print("RT window: ±", RT_WINDOW, "s")
print("MS2 spectra examined:", len({
    r["scan_index"] for r in rows
}))
print("Precursor records:", len(rows))
print("Output:", OUT)

with_precursor = [
    r for r in rows if r["precursor_mz"] != ""
]
near_track = [
    r for r in with_precursor
    if r["within_track_tolerance"] == 1
]

print("Precursor records within ±15 mDa of a candidate track:",
      len(near_track), "/", len(with_precursor))

if with_precursor:
    print("\nNearest-track summary:")
    for track in TRACKS:
        selected = [
            r for r in with_precursor
            if r["nearest_track"] == track
        ]
        if selected:
            print(f"Track {track}: {len(selected)} precursor records")

print("\nInterpretation caution:")
print("A precursor m/z match is supportive, not proof of a feature-to-PSM link.")
print("Check isolation bounds, charge, scan identity and peptide evidence separately.")
print("No mzML files modified.")
