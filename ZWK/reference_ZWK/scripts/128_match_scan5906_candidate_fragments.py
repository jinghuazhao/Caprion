#!/usr/bin/env python3
"""Test a candidate peptide against native MS2 scan 5906 product ions."""

import csv
from pathlib import Path

PEAKS_FILE = Path(
    "reports/feature10275997325121874585_ms2_scan5905_peaks_v1.tsv"
)
OUT_FILE = Path(
    "reports/feature10275997325121874585_scan5906_fragment_matches_v1.tsv"
)

PEPTIDE = "SSSYSKQFTSSTSYNRGDSTFESK"
PROTON = 1.007276466621
H2O = 18.010564684
CO = 27.994914620
TOLERANCE_DA = 0.02

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

GLYCAN_IONS = {
    "HexNAc oxonium": 204.086649,
    "HexNAc-H2O": 186.076084,
    "HexNAc-2H2O": 168.065519,
    "HexNAc-related 138": 138.0550,
    "HexNAc-related 126": 126.0550,
}

if not PEAKS_FILE.is_file():
    raise FileNotFoundError(PEAKS_FILE)

peaks = []
with PEAKS_FILE.open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        peaks.append({
            "mz": float(row["mz"]),
            "intensity": float(row["intensity"]),
        })

peaks.sort(key=lambda p: p["mz"])
aa_mass = [AA[a] for a in PEPTIDE]
total_residue_mass = sum(aa_mass)
neutral_peptide_mass = total_residue_mass + H2O
precursor_mz = (neutral_peptide_mass + 2 * PROTON) / 2

print("Candidate:", PEPTIDE)
print(f"Unmodified neutral peptide mass: {neutral_peptide_mass:.5f}")
print(f"Unmodified [M+2H]2+ m/z: {precursor_mz:.5f}")
print("Target precursor m/z: 1341.06800")
print(f"Unmodified precursor error: "
      f"{(precursor_mz - 1341.068) * 1000:+.2f} mDa")
print(f"Product-ion tolerance: ±{TOLERANCE_DA:.3f} Da")

def nearest_peak(theoretical):
    peak = min(peaks, key=lambda p: abs(p["mz"] - theoretical))
    error = peak["mz"] - theoretical
    if abs(error) > TOLERANCE_DA:
        return None
    return peak, error, error / theoretical * 1e6

results = []

# Standard singly charged b/y ions for the unmodified peptide.
for cut in range(1, len(PEPTIDE)):
    prefix_mass = sum(aa_mass[:cut])
    suffix_mass = sum(aa_mass[cut:])

    ions = [
        (f"b{cut}", prefix_mass + PROTON),
        (f"y{len(PEPTIDE)-cut}",
         suffix_mass + H2O + PROTON),
    ]

    for ion_name, theoretical in ions:
        match = nearest_peak(theoretical)
        if match:
            peak, error, ppm = match
            results.append({
                "ion": ion_name,
                "theoretical_mz": theoretical,
                "observed_mz": peak["mz"],
                "error_mDa": error * 1000,
                "error_ppm": ppm,
                "intensity": peak["intensity"],
                "kind": "unmodified peptide ion",
            })

# Common glycan-related diagnostic ions.
for ion_name, theoretical in GLYCAN_IONS.items():
    match = nearest_peak(theoretical)
    if match:
        peak, error, ppm = match
        results.append({
            "ion": ion_name,
            "theoretical_mz": theoretical,
            "observed_mz": peak["mz"],
            "error_mDa": error * 1000,
            "error_ppm": ppm,
            "intensity": peak["intensity"],
            "kind": "diagnostic ion",
        })

results.sort(key=lambda r: (r["kind"], r["theoretical_mz"]))

fields = [
    "kind", "ion", "theoretical_mz", "observed_mz",
    "error_mDa", "error_ppm", "intensity",
]
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
with OUT_FILE.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

print("\nMatched ions:")
for r in results:
    print(
        f"{r['kind']:25s} {r['ion']:10s} "
        f"theory={r['theoretical_mz']:.5f} "
        f"observed={r['observed_mz']:.5f} "
        f"error={r['error_mDa']:+.2f} mDa "
        f"intensity={r['intensity']:.1f}"
    )

peptide_matches = [
    r for r in results if r["kind"] == "unmodified peptide ion"
]
print("\nUnmodified b/y matches:", len(peptide_matches))
print("Report:", OUT_FILE)
print(
    "Caution: this is a preliminary screen using the top-50 peak list; "
    "it does not model glycan-bearing fragments, neutral losses, "
    "or multiply charged product ions."
)
