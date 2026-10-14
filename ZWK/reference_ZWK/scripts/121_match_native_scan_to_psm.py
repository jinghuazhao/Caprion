#!/usr/bin/env python3
"""Match native mzML scan 5906 to PSMs for the exact acquisition."""

from pathlib import Path
import csv

ACQUISITION = "szwk900404i19201xms1"
NATIVE_SCAN = 5906

PSM_FILE = Path("search_98") / ACQUISITION / "psm.tsv"
OUT = Path(
    "reports/feature10275997325121874585_native_scan5906_psm_match_v1.tsv"
)

if not PSM_FILE.exists():
    raise FileNotFoundError(f"Expected PSM table not found: {PSM_FILE}")

with PSM_FILE.open(newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fields = reader.fieldnames
    if not fields:
        raise RuntimeError("PSM table has no header")

    # FragPipe normally stores names such as:
    # szwk900404i19201xms1.05876.05876.2
    # The middle integer is the spectrum/scan number used by FragPipe.
    name_candidates = [
        field for field in fields
        if field.lower() in {
            "spectrum", "spectrumname", "spectrum_name",
            "spec", "specid", "spec_id"
        }
    ]
    if not name_candidates:
        raise RuntimeError(
            "Could not identify the spectrum-name column. "
            f"Available columns: {fields}"
        )

    name_col = name_candidates[0]
    matched = []
    for row in reader:
        name = row.get(name_col, "")
        if not name.startswith(ACQUISITION + "."):
            continue

        parts = name.split(".")
        if len(parts) >= 4:
            try:
                scan_a = int(parts[1])
                scan_b = int(parts[2])
            except ValueError:
                continue

            if scan_a == NATIVE_SCAN or scan_b == NATIVE_SCAN:
                matched.append(row)

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(matched)

print("PSM table:", PSM_FILE)
print("Acquisition:", ACQUISITION)
print("Native mzML scan:", NATIVE_SCAN)
print("Spectrum-name column:", name_col)
print("Matching PSM rows:", len(matched))
print("Output:", OUT)

for i, row in enumerate(matched, start=1):
    print(f"\n--- Match {i} ---")
    for key in (
        name_col, "Peptide", "Modified Peptide", "Charge",
        "Precursor", "Hyperscore", "Expectation", "QValue",
        "Protein", "Gene"
    ):
        if key in row:
            print(f"{key}: {row[key]}")
