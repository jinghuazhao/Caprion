#!/usr/bin/env python3
"""Screen native MS2 scan 5906 for peptide and HexNAc-bearing fragments."""

from pathlib import Path
import csv
import math
from pyopenms import MSExperiment, MzMLFile

ROOT = Path.cwd()
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
OUT = ROOT / "reports/feature10275997325121874585_scan5906_glycopeptide_fragments_v1.tsv"

TARGET_SCAN = 5906
TARGET_RT = 764.23206
RT_TOL = 0.05
TOL_DA = 0.020

PEPTIDE = "SSSYSKQFTSSTSYNRGDSTFESK"

# Monoisotopic residue masses, Da.
AA = {
    "A": 71.037113805, "R": 156.10111105, "N": 114.04292747,
    "D": 115.026943065, "C": 103.009184505, "E": 129.042593135,
    "Q": 128.05857754, "G": 57.021463735, "H": 137.058911875,
    "I": 113.084063975, "L": 113.084063975, "K": 128.094963015,
    "M": 131.040484645, "F": 147.068413945, "P": 97.052763875,
    "S": 87.032028435, "T": 101.047678505, "W": 186.07931298,
    "Y": 163.063328575, "V": 99.068413945,
}
PROTON = 1.007276466621
H2O = 18.010564684
HEXNAC = 203.07937252
HEX = 162.05282342
FUC = 146.05790880

if not MZML.exists():
    raise FileNotFoundError(MZML)

exp = MSExperiment()
MzMLFile().load(str(MZML), exp)

matches = []
for spec in exp:
    if spec.getMSLevel() != 2:
        continue
    if f"scan={TARGET_SCAN}" not in spec.getNativeID():
        continue
    rt = spec.getRT()
    if abs(rt - TARGET_RT) > RT_TOL:
        raise RuntimeError(f"Scan RT mismatch: {rt:.5f} s")
    mz_array, intensity_array = spec.get_peaks()
    peaks = [
        (float(mz), float(inten))
        for mz, inten in zip(mz_array, intensity_array)
    ]
    break
else:
    raise RuntimeError(f"MS2 scan {TARGET_SCAN} not found")

peaks.sort()
print(f"Scan {TARGET_SCAN}; RT={rt:.5f}; full peaks={len(peaks)}")

def mz_from_neutral(neutral_mass, z=1):
    return (neutral_mass + z * PROTON) / z

def match_peak(theory):
    candidates = [
        (abs(mz - theory), mz, intensity)
        for mz, intensity in peaks
        if abs(mz - theory) <= TOL_DA
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda x: x[0])

# Cumulative b-ion residue masses; y ions include water.
residues = [AA[a] for a in PEPTIDE]
total_residue_mass = sum(residues)
ion_defs = []

for n in range(1, len(PEPTIDE)):
    b_neutral = sum(residues[:n])
    y_neutral = sum(residues[n:]) + H2O

    ion_defs.append((f"b{n}", b_neutral, "unmodified"))
    ion_defs.append((f"y{len(PEPTIDE)-n}", y_neutral, "unmodified"))

    # Screen a single HexNAc attached to either fragment.
    ion_defs.append((f"b{n}+HexNAc", b_neutral + HEXNAC, "HexNAc-bearing"))
    ion_defs.append((f"y{len(PEPTIDE)-n}+HexNAc",
                     y_neutral + HEXNAC, "HexNAc-bearing"))

    # Common neutral losses from HexNAc-bearing fragments.
    ion_defs.append((f"b{n}+HexNAc-H2O",
                     b_neutral + HEXNAC - H2O, "HexNAc neutral loss"))
    ion_defs.append((f"y{len(PEPTIDE)-n}+HexNAc-H2O",
                     y_neutral + HEXNAC - H2O, "HexNAc neutral loss"))

# Diagnostic glycan ions/losses, included as reference hypotheses.
diagnostics = [
    ("HexNAc oxonium", HEXNAC + PROTON),
    ("HexNAc-H2O", HEXNAC - H2O + PROTON),
    ("HexNAc-2H2O", HEXNAC - 2 * H2O + PROTON),
    ("HexNAc-related 138", 138.05500),
    ("HexNAc-related 126", 126.05500),
    ("Hex oxonium", HEX + PROTON),
    ("Hex-H2O", HEX - H2O + PROTON),
    ("Fuc oxonium", FUC + PROTON),
]

rows = []
for label, theory in diagnostics:
    hit = match_peak(theory)
    if hit:
        err, obs, inten = hit
        rows.append({
            "ion": label, "class": "diagnostic",
            "charge": 1, "theoretical_mz": theory,
            "observed_mz": obs, "error_mDa": (obs-theory)*1000,
            "intensity": inten,
        })

for label, neutral, ion_class in ion_defs:
    for z in (1, 2):
        theory = mz_from_neutral(neutral, z)
        hit = match_peak(theory)
        if hit:
            err, obs, inten = hit
            rows.append({
                "ion": label, "class": ion_class,
                "charge": z, "theoretical_mz": theory,
                "observed_mz": obs, "error_mDa": (obs-theory)*1000,
                "intensity": inten,
            })

# Remove duplicate observations for the same theoretical ion/charge.
rows.sort(key=lambda r: (-r["intensity"], abs(r["error_mDa"])))
OUT.parent.mkdir(parents=True, exist_ok=True)
fields = [
    "ion", "class", "charge", "theoretical_mz",
    "observed_mz", "error_mDa", "intensity",
]
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"Matched theoretical ions: {len(rows)}")
print(f"Distinct matched labels: {len(set(r['ion'] for r in rows))}")
print(f"Report: {OUT}")
print("\nStrongest matches:")
for r in rows[:30]:
    print(
        f"{r['ion']:24s} {r['class']:20s} "
        f"z={r['charge']} theory={r['theoretical_mz']:.5f} "
        f"obs={r['observed_mz']:.5f} "
        f"error={r['error_mDa']:+.2f} mDa "
        f"intensity={r['intensity']:.1f}"
    )

print("\nCaution: this is a hypothesis screen, not a glycopeptide identification.")
print("HexNAc-bearing fragment matches need sequence coverage and competing-candidate controls.")
