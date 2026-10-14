#!/usr/bin/env python3
"""Rank FragPipe peptide sequences against native MS2 scan 5906."""

from pathlib import Path
import csv
import math
import re
from collections import defaultdict
from pyopenms import MSExperiment, MzMLFile

ROOT = Path.cwd()
PSM_ROOT = ROOT / "search_98"
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
OUT = ROOT / "reports/feature10275997325121874585_scan5906_database_candidates_v1.tsv"

SCAN = 5906
TARGET_RT = 764.23206
TARGET_MZ = 1341.068
TARGET_Z = 2
FRAGMENT_TOL = 0.020
PRECURSOR_WINDOW = 0.100

PROTON = 1.007276466621
H2O = 18.010564684
AA = {
    "A": 71.037113805, "R": 156.10111105, "N": 114.04292747,
    "D": 115.026943065, "C": 103.009184505, "E": 129.042593135,
    "Q": 128.05857754, "G": 57.021463735, "H": 137.058911875,
    "I": 113.084063975, "L": 113.084063975, "K": 128.094963015,
    "M": 131.040484645, "F": 147.068413945, "P": 97.052763875,
    "S": 87.032028435, "T": 101.047678505, "W": 186.07931298,
    "Y": 163.063328575, "V": 99.068413945,
}

# Load the exact native MS2 spectrum.
exp = MSExperiment()
MzMLFile().load(str(MZML), exp)
peaks = None
scan_rt = None
for spec in exp:
    if spec.getMSLevel() == 2 and f"scan={SCAN}" in spec.getNativeID():
        scan_rt = spec.getRT()
        mzs, intensities = spec.get_peaks()
        peaks = sorted((float(m), float(i)) for m, i in zip(mzs, intensities))
        break

if peaks is None:
    raise RuntimeError(f"MS2 scan {SCAN} not found")
if abs(scan_rt - TARGET_RT) > 0.05:
    raise RuntimeError(f"Unexpected scan RT: {scan_rt}")

mz_values = [p[0] for p in peaks]
max_intensity = max((p[1] for p in peaks), default=1.0)
base_threshold = max_intensity * 0.001

def peptide_mass(seq):
    try:
        return sum(AA[a] for a in seq) + H2O
    except KeyError:
        return None

def fragment_mz(neutral, z=1):
    return (neutral + z * PROTON) / z

def fragment_matches(seq):
    """Return matched b/y labels, unique observed peaks and intensity score."""
    residues = [AA[a] for a in seq]
    hits = []
    used_peaks = set()
    n = len(seq)
    b_mass = 0.0
    suffix_mass = sum(residues) + H2O

    for cut in range(1, n):
        b_mass += residues[cut - 1]
        suffix_mass -= residues[cut - 1]

        ions = [
            (f"b{cut}", b_mass),
            (f"y{n-cut}", suffix_mass),
        ]
        for label, neutral in ions:
            for z in (1, 2):
                theory = fragment_mz(neutral, z)
                # Linear scan is small (62 peaks per spectrum).
                candidates = [
                    (abs(m - theory), idx, m, inten)
                    for idx, (m, inten) in enumerate(peaks)
                    if abs(m - theory) <= FRAGMENT_TOL
                ]
                if not candidates:
                    continue
                err, idx, observed, intensity = min(candidates)
                hits.append((label, z, theory, observed, err, intensity, idx))
                used_peaks.add(idx)

    strong_hits = sum(
        1 for h in hits if h[5] >= base_threshold
    )
    intensity_score = sum(
        math.sqrt(h[5] / max_intensity) for h in hits
    )
    return hits, len(used_peaks), strong_hits, intensity_score

def parse_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

# Read PSM tables and collect unique unmodified peptide sequences.
# We only use peptide sequences, not PSM scores, as a candidate list.
seq_sources = defaultdict(set)
psm_files = list(PSM_ROOT.glob("**/psm.tsv"))
if not psm_files:
    raise FileNotFoundError(f"No psm.tsv files under {PSM_ROOT}")

for file in psm_files:
    try:
        with file.open(newline="", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            seq_col = next(
                (c for c in ("Peptide", "peptide", "Peptide Sequence",
                             "Modified Peptide") if c in reader.fieldnames),
                None
            )
            if seq_col is None:
                continue
            for row in reader:
                seq = (row.get(seq_col) or "").strip()
                # Only standard amino-acid sequences; modified notation excluded.
                if not seq or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
                    continue
                mass = peptide_mass(seq)
                if mass is None:
                    continue
                mz = fragment_mz(mass, TARGET_Z)
                if abs(mz - TARGET_MZ) <= PRECURSOR_WINDOW:
                    seq_sources[seq].add(str(file.parent.name))
    except (OSError, csv.Error):
        continue

print(f"MS2 scan={SCAN}; RT={scan_rt:.5f}; peaks={len(peaks)}")
print(f"PSM tables scanned: {len(psm_files)}")
print(f"Unique unmodified sequences within ±{PRECURSOR_WINDOW:.3f} Da: "
      f"{len(seq_sources)}")

results = []
for seq, runs in seq_sources.items():
    neutral = peptide_mass(seq)
    calc_mz = fragment_mz(neutral, TARGET_Z)
    hits, n_peaks, n_strong, intensity_score = fragment_matches(seq)
    results.append({
        "peptide": seq,
        "length": len(seq),
        "calculated_precursor_mz": calc_mz,
        "precursor_error_mDa": (calc_mz - TARGET_MZ) * 1000,
        "precursor_error_ppm": (calc_mz - TARGET_MZ) / TARGET_MZ * 1e6,
        "matched_fragment_ions": len(hits),
        "unique_matched_peaks": n_peaks,
        "matched_ions_above_0.1pct_basepeak": n_strong,
        "fragment_intensity_score": intensity_score,
        "fragpipe_runs": ";".join(sorted(runs)),
        "fragment_matches": ";".join(
            f"{h[0]}(z{h[1]},err={h[4]*1000:+.2f}mDa,I={h[5]:.0f})"
            for h in sorted(hits, key=lambda x: -x[5])
        ),
    })

results.sort(
    key=lambda r: (
        -r["unique_matched_peaks"],
        -r["matched_ions_above_0.1pct_basepeak"],
        -r["fragment_intensity_score"],
        abs(r["precursor_error_ppm"]),
    )
)

OUT.parent.mkdir(parents=True, exist_ok=True)
fields = [
    "peptide", "length", "calculated_precursor_mz",
    "precursor_error_mDa", "precursor_error_ppm",
    "matched_fragment_ions", "unique_matched_peaks",
    "matched_ions_above_0.1pct_basepeak", "fragment_intensity_score",
    "fragpipe_runs", "fragment_matches",
]
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

print(f"Candidates ranked: {len(results)}")
print(f"Report: {OUT}")
print("\nTop 20 candidates:")
for r in results[:20]:
    print(
        f"{r['peptide']:32s} "
        f"precursor={r['precursor_error_ppm']:+7.1f} ppm "
        f"ions={r['matched_fragment_ions']:2d} "
        f"peaks={r['unique_matched_peaks']:2d} "
        f"strong={r['matched_ions_above_0.1pct_basepeak']:2d}"
    )

print("\nCaution: candidate generation uses unmodified sequences already present in PSMs.")
print("This cannot discover unseen peptides or correctly model glycan-bearing fragments.")
print("Precursor and fragment matches must be evaluated together; this is not an identification.")
