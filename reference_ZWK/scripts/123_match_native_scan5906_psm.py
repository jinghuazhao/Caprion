#!/usr/bin/env python3
"""Match native MS2 scan 5906 to a FragPipe PSM for one acquisition."""

import csv
import re
from pathlib import Path

ACQUISITION = "szwk900404i19201xms1"
NATIVE_SCAN = 5906

PSM_FILE = Path("search_98/szwk900404i19201/psm.tsv")
OUT_FILE = Path(
    "reports/feature10275997325121874585_native_scan5906_psm_match_v1.tsv"
)
NEAR_FILE = Path(
    "reports/feature10275997325121874585_native_scan5906_nearby_psms_v1.tsv"
)

if not PSM_FILE.is_file():
    raise FileNotFoundError(f"PSM table not found: {PSM_FILE}")

OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

exact_matches = []
nearby = []
same_acquisition = []
unparsed = 0

with PSM_FILE.open(newline="", errors="replace") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    if not reader.fieldnames or "Spectrum" not in reader.fieldnames:
        raise ValueError(
            f"Expected 'Spectrum' column; found: {reader.fieldnames}"
        )

    fields = reader.fieldnames

    for row in reader:
        spectrum = (row.get("Spectrum") or "").strip()
        if not spectrum.startswith(ACQUISITION + "."):
            continue

        match = re.match(
            rf"^{re.escape(ACQUISITION)}\.(\d+)\.(\d+)\.(\d+)$",
            spectrum,
        )

        if not match:
            unparsed += 1
            continue

        scan_a, scan_b, charge = map(int, match.groups())
        info = {
            "Spectrum": spectrum,
            "Parsed scan A": scan_a,
            "Parsed scan B": scan_b,
            "Parsed charge": charge,
            "Peptide": row.get("Peptide", ""),
            "Modified Peptide": row.get("Modified Peptide", ""),
            "Retention": row.get("Retention", ""),
            "Observed M/Z": row.get("Observed M/Z", ""),
            "Expectation": row.get("Expectation", ""),
            "Probability": row.get("Probability", ""),
        }
        same_acquisition.append((min(abs(scan_a - NATIVE_SCAN),
                                     abs(scan_b - NATIVE_SCAN)), info, row))

        if scan_a == NATIVE_SCAN or scan_b == NATIVE_SCAN:
            exact_matches.append(row)

        if (
            abs(scan_a - NATIVE_SCAN) <= 20
            or abs(scan_b - NATIVE_SCAN) <= 20
        ):
            nearby.append(info)

# Write exact PSM rows, including a header if there are no matches.
with OUT_FILE.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(exact_matches)

# Write nearby scan identifiers and their PSM details for diagnosis.
near_fields = [
    "Spectrum", "Parsed scan A", "Parsed scan B", "Parsed charge",
    "Peptide", "Modified Peptide", "Retention", "Observed M/Z",
    "Expectation", "Probability",
]
with NEAR_FILE.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=near_fields, delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(nearby)

print("PSM table:", PSM_FILE)
print("Acquisition:", ACQUISITION)
print("Spectrum rows parsed for acquisition:", len(same_acquisition))
print("Rows with scan number 5906:", len(exact_matches))
print("Rows within 20 scan numbers of 5906:", len(nearby))
print("Unparsed Spectrum names for acquisition:", unparsed)
print("Exact-match report:", OUT_FILE)
print("Nearby-scan report:", NEAR_FILE)

if exact_matches:
    print("\nExact PSM matches:")
    for row in exact_matches:
        print(
            "Spectrum =", row.get("Spectrum"),
            "| Peptide =", row.get("Peptide"),
            "| Modified Peptide =", row.get("Modified Peptide"),
            "| Charge =", row.get("Charge"),
            "| Retention =", row.get("Retention"),
        )
else:
    print("\nNo exact PSM match. Nearest spectrum IDs for this acquisition:")
    for distance, info, _ in sorted(same_acquisition, key=lambda x: x[0])[:10]:
        print(
            f"distance={distance}: {info['Spectrum']} | "
            f"peptide={info['Peptide']} | RT={info['Retention']}"
        )
