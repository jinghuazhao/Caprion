#!/usr/bin/env python3
"""Audit PSMs near a candidate feature's RT and precursor m/z."""

import csv
import math
import re
from pathlib import Path

ACQUISITION = "szwk900404i19201xms1"
PSM_FILE = Path("search_98/szwk900404i19201/psm.tsv")

TARGET_RT = 764.23206
RT_TOLERANCE = 5.0

# Provisional track references from the earlier native-centroid audit.
TRACKS = {
    2: 1341.067749,
    3: 1341.570496,
    4: 1342.068359,
    5: 1342.571045,
    6: 1343.074463,
    7: 1343.573608,
}

OUT_FILE = Path(
    "reports/feature10275997325121874585_psms_near_rt_764s_v1.tsv"
)

def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None

if not PSM_FILE.is_file():
    raise FileNotFoundError(f"PSM table not found: {PSM_FILE}")

OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

output_fields = [
    "Spectrum",
    "Native scan number in FragPipe name",
    "RT",
    "RT difference from target (s)",
    "Charge",
    "Peptide",
    "Modified Peptide",
    "Observed M/Z",
    "Calculated M/Z",
    "Closest provisional track",
    "Closest track reference m/z",
    "Observed m/z minus track (mDa)",
    "Within 15 mDa of provisional track",
    "Expectation",
    "Probability",
]

results = []
acquisition_rows = 0
unparsed_rt = 0
unparsed_mz = 0

with PSM_FILE.open(newline="", errors="replace") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    required = {"Spectrum", "Retention", "Observed M/Z"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(
            f"Missing expected PSM columns: {sorted(missing)}; "
            f"found {reader.fieldnames}"
        )

    for row in reader:
        spectrum = (row.get("Spectrum") or "").strip()
        if not spectrum.startswith(ACQUISITION + "."):
            continue

        acquisition_rows += 1
        rt = number(row.get("Retention"))
        observed_mz = number(row.get("Observed M/Z"))

        if rt is None:
            unparsed_rt += 1
            continue

        rt_delta = rt - TARGET_RT
        if abs(rt_delta) > RT_TOLERANCE:
            continue

        if observed_mz is None:
            unparsed_mz += 1
            closest_track = ""
            closest_reference = ""
            mz_error_mda = ""
            within_15_mda = ""
        else:
            track_id, reference_mz = min(
                TRACKS.items(),
                key=lambda item: abs(observed_mz - item[1]),
            )
            mz_error_mda = (observed_mz - reference_mz) * 1000
            closest_track = track_id
            closest_reference = reference_mz
            within_15_mda = abs(mz_error_mda) <= 15

        match = re.match(
            rf"^{re.escape(ACQUISITION)}\.(\d+)\.(\d+)\.(\d+)$",
            spectrum,
        )
        scan_id = match.group(1) if match else ""

        results.append({
            "Spectrum": spectrum,
            "Native scan number in FragPipe name": scan_id,
            "RT": rt,
            "RT difference from target (s)": round(rt_delta, 5),
            "Charge": row.get("Charge", ""),
            "Peptide": row.get("Peptide", ""),
            "Modified Peptide": row.get("Modified Peptide", ""),
            "Observed M/Z": observed_mz if observed_mz is not None else "",
            "Calculated M/Z": row.get("Calculated M/Z", ""),
            "Closest provisional track": closest_track,
            "Closest track reference m/z": closest_reference,
            "Observed m/z minus track (mDa)": (
                round(mz_error_mda, 3)
                if isinstance(mz_error_mda, (int, float))
                else ""
            ),
            "Within 15 mDa of provisional track": within_15_mda,
            "Expectation": row.get("Expectation", ""),
            "Probability": row.get("Probability", ""),
        })

results.sort(key=lambda row: abs(row["RT difference from target (s)"]))

with OUT_FILE.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=output_fields, delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("PSM file:", PSM_FILE)
print("Acquisition rows:", acquisition_rows)
print(f"Target RT: {TARGET_RT:.5f} s")
print(f"RT window: {TARGET_RT - RT_TOLERANCE:.5f}–"
      f"{TARGET_RT + RT_TOLERANCE:.5f} s")
print("PSMs in RT window:", len(results))
print("Rows with unparseable RT:", unparsed_rt)
print("Rows with unparseable m/z encountered:", unparsed_mz)
print("Output:", OUT_FILE)

for row in results:
    print(
        f"RT={row['RT']:.5f}, "
        f"delta_RT={row['RT difference from target (s)']:+.5f}s, "
        f"scan={row['Native scan number in FragPipe name']}, "
        f"observed_mz={row['Observed M/Z']}, "
        f"closest_track={row['Closest provisional track']}, "
        f"mz_error_mDa={row['Observed m/z minus track (mDa)']}, "
        f"peptide={row['Peptide']}"
    )
