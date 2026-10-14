#!/usr/bin/env python3
"""Inspect one candidate MS2 scan and search existing FragPipe output."""

from pathlib import Path
import csv
import re
import numpy as np
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
FRAGPIPE = Path("search_98")
OUT_PEAKS = Path(
    "reports/feature10275997325121874585_ms2_scan5905_peaks_v1.tsv"
)
OUT_PSMS = Path(
    "reports/feature10275997325121874585_ms2_scan5905_psm_matches_v1.tsv"
)

TARGET_SCAN_INDEX = 5905
TARGET_RT = 764.23206
RT_TOL = 0.02
TOP_N = 50

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

# First locate by scan index, then validate by MS level and RT.
if TARGET_SCAN_INDEX >= exp.size():
    raise RuntimeError(f"Scan index {TARGET_SCAN_INDEX} is outside the experiment")

spectrum = exp[TARGET_SCAN_INDEX]
if spectrum.getMSLevel() != 2:
    raise RuntimeError(
        f"Index {TARGET_SCAN_INDEX} is MS level {spectrum.getMSLevel()}, not MS2"
    )
if abs(float(spectrum.getRT()) - TARGET_RT) > RT_TOL:
    raise RuntimeError(
        f"RT mismatch: observed {spectrum.getRT():.5f}, expected {TARGET_RT:.5f}"
    )

print("File:", MZML)
print("Spectrum index:", TARGET_SCAN_INDEX)
print("RT:", spectrum.getRT())
print("MS level:", spectrum.getMSLevel())
print("Native spectrum type:", spectrum.getType())

print("\n=== Spectrum metadata ===")
try:
    print("Native ID:", spectrum.getNativeID())
except Exception as exc:
    print("Native ID unavailable:", repr(exc))

precursors = spectrum.getPrecursors()
print("Number of precursor records:", len(precursors))

for i, p in enumerate(precursors, start=1):
    pmz = float(p.getMZ())
    charge = int(p.getCharge())
    lower_offset = float(p.getIsolationWindowLowerOffset())
    upper_offset = float(p.getIsolationWindowUpperOffset())
    print(f"Precursor {i}:")
    print(f"  m/z: {pmz:.7f}")
    print(f"  charge: {charge}")
    print(f"  isolation lower offset: {lower_offset:.7f}")
    print(f"  isolation upper offset: {upper_offset:.7f}")
    print(f"  isolation bounds: {pmz-lower_offset:.7f}–{pmz+upper_offset:.7f}")

mz, intensity = spectrum.get_peaks()
mz = np.asarray(mz, dtype=float)
intensity = np.asarray(intensity, dtype=float)

order = np.argsort(intensity)[::-1][:TOP_N]
peak_rows = [
    {"rank": rank, "mz": float(mz[j]), "intensity": float(intensity[j])}
    for rank, j in enumerate(order, start=1)
]

OUT_PEAKS.parent.mkdir(parents=True, exist_ok=True)
with OUT_PEAKS.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=["rank", "mz", "intensity"], delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(peak_rows)

print("\n=== Top product-ion peaks ===")
for r in peak_rows[:20]:
    print(f"{r['rank']:2d}\t{r['mz']:.6f}\t{r['intensity']:.1f}")
print("Peak table:", OUT_PEAKS)

# Search existing tabular FragPipe outputs. A scan identifier can be
# represented in different ways, so test common forms and retain matches.
patterns = [
    re.compile(r"\b5905\b"),
    re.compile(r"scan(?:=|:|_)?0*5905\b", re.IGNORECASE),
]
psm_matches = []

for path in FRAGPIPE.rglob("*"):
    if not path.is_file() or path.suffix.lower() not in {".tsv", ".txt", ".csv"}:
        continue
    try:
        with path.open("r", errors="replace") as f:
            header = f.readline().rstrip("\n")
            if not header:
                continue
            for line_number, line in enumerate(f, start=2):
                if any(p.search(line) for p in patterns):
                    psm_matches.append({
                        "file": str(path),
                        "line_number": line_number,
                        "header": header,
                        "row": line.rstrip("\n"),
                    })
    except (OSError, UnicodeError):
        continue

with OUT_PSMS.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["file", "line_number", "header", "row"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(psm_matches)

print("\n=== Existing FragPipe text-table matches ===")
print("Rows containing scan number 5905:", len(psm_matches))
for r in psm_matches[:20]:
    print(f"{r['file']}:{r['line_number']}")
    print("  ", r["row"][:400])
print("Search report:", OUT_PSMS)
print("\nNo source files modified; no peptide or IG assignment made.")
