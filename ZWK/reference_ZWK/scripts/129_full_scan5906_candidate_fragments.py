#!/usr/bin/env python3
"""Full-spectrum b/y screening for native MS2 scan 5906."""

import csv
from pathlib import Path
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
OUT = Path("reports/feature10275997325121874585_scan5906_full_fragment_matches_v1.tsv")

NATIVE_SCAN = 5906
TARGET_RT = 764.23206
TOLERANCE_DA = 0.020
PROTON = 1.007276466621
WATER = 18.010564684

PEPTIDE = "SSSYSKQFTSSTSYNRGDSTFESK"

AA = {
    "A": 71.037113805, "R": 156.101111050,
    "N": 114.042927470, "D": 115.026943065,
    "C": 103.009184505, "E": 129.042593135,
    "Q": 128.058577540, "G": 57.021463735,
    "H": 137.058911875, "I": 113.084063975,
    "L": 113.084063975, "K": 128.094963015,
    "M": 131.040484645, "F": 147.068413945,
    "P": 97.052763875, "S": 87.032028435,
    "T": 101.047678505, "W": 186.079312980,
    "Y": 163.063328575, "V": 99.068413945,
}

if not MZML.is_file():
    raise FileNotFoundError(MZML)

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

matches = [
    s for s in exp
    if s.getMSLevel() == 2
    and f"scan={NATIVE_SCAN}" in s.getNativeID()
]
if len(matches) != 1:
    raise RuntimeError(
        f"Expected one MS2 scan {NATIVE_SCAN}; found {len(matches)}"
    )

spectrum = matches[0]
rt = float(spectrum.getRT())
if abs(rt - TARGET_RT) > 0.05:
    raise RuntimeError(f"Unexpected RT for scan: {rt:.6f}")

mz, intensity = spectrum.get_peaks()
peaks = sorted(
    [(float(m), float(i)) for m, i in zip(mz, intensity)],
    key=lambda x: x[0],
)
print(f"Native scan={NATIVE_SCAN}; RT={rt:.5f}; peaks={len(peaks)}")

residue_masses = [AA[a] for a in PEPTIDE]
neutral_mass = sum(residue_masses) + WATER
precursor_mz = (neutral_mass + 2 * PROTON) / 2
print(f"Candidate unmodified [M+2H]2+ = {precursor_mz:.5f}")
print(f"Target precursor m/z = 1341.06800")
print(f"Difference = {(precursor_mz - 1341.068) * 1000:+.2f} mDa")

theoretical = []
for cut in range(1, len(PEPTIDE)):
    prefix = sum(residue_masses[:cut])
    suffix = sum(residue_masses[cut:])
    for series, neutral in (
        ("b", prefix),
        ("y", suffix + WATER),
    ):
        for charge in (1, 2, 3):
            value = (neutral + charge * PROTON) / charge
            theoretical.append((series, cut, charge, value))

rows = []
for series, cut, charge, theo in theoretical:
    candidates = [
        (obs, inten) for obs, inten in peaks
        if abs(obs - theo) <= TOLERANCE_DA
    ]
    for obs, inten in candidates:
        error = obs - theo
        rows.append({
            "ion": f"{series}{cut}",
            "charge": charge,
            "theoretical_mz": theo,
            "observed_mz": obs,
            "error_mDa": error * 1000,
            "error_ppm": error / theo * 1e6,
            "intensity": inten,
        })

rows.sort(key=lambda r: (r["ion"][0], int(r["ion"][1:]), r["charge"]))

OUT.parent.mkdir(parents=True, exist_ok=True)
fields = [
    "ion", "charge", "theoretical_mz", "observed_mz",
    "error_mDa", "error_ppm", "intensity",
]
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print("Matched theoretical ion observations:", len(rows))
print("Distinct matched ion labels:", len(set(r["ion"] for r in rows)))
print("Report:", OUT)
for row in rows:
    print(
        f"{row['ion']:5s} z={row['charge']} "
        f"theory={row['theoretical_mz']:.5f} "
        f"obs={row['observed_mz']:.5f} "
        f"error={row['error_mDa']:+.2f} mDa "
        f"intensity={row['intensity']:.1f}"
    )
print(
    "Interpretation caution: this tests only unmodified b/y ions. "
    "A glycosylated peptide may fragment differently; matches are "
    "not sufficient for identification without specificity and controls."
)
