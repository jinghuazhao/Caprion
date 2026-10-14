#!/usr/bin/env python3
"""Reconcile nearby FragPipe PSMs with native mzML MS2 spectra by RT."""

import csv
import re
from pathlib import Path

import pyopenms as oms

ACQ = "szwk900404i19201xms1"
TARGET_RT = 764.23206
RT_WINDOW = 5.0
MATCH_TOLERANCE = 0.05  # seconds

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
PSM_FILE = Path("search_98/szwk900404i19201/psm.tsv")
OUT = Path("reports/feature10275997325121874585_native_psm_rt_mapping_v1.tsv")


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def native_scan_number(native_id):
    match = re.search(r"(?:^|\s)scan=(\d+)", native_id)
    return int(match.group(1)) if match else None


if not MZML.is_file():
    raise FileNotFoundError(MZML)
if not PSM_FILE.is_file():
    raise FileNotFoundError(PSM_FILE)

# Load native MS2 spectra, retaining native IDs and precursor information.
exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

native_ms2 = []
target_scan = None

for index, spectrum in enumerate(exp):
    if spectrum.getMSLevel() != 2:
        continue

    rt = float(spectrum.getRT())
    native_id = spectrum.getNativeID()
    scan = native_scan_number(native_id)

    precursors = spectrum.getPrecursors()
    precursor = precursors[0] if precursors else None
    precursor_mz = float(precursor.getMZ()) if precursor else None
    precursor_charge = int(precursor.getCharge()) if precursor else None

    item = {
        "native_rt": rt,
        "native_id": native_id,
        "native_scan": scan,
        "spectrum_index": index,
        "precursor_mz": precursor_mz,
        "precursor_charge": precursor_charge,
    }
    native_ms2.append(item)

    if scan == 5906:
        target_scan = item

if target_scan is None:
    print("WARNING: no MS2 spectrum with native ID scan=5906 was found.")
else:
    print("Target native MS2 scan:")
    print(target_scan)

# Read PSMs belonging only to this acquisition and near the target RT.
psms = []

with PSM_FILE.open(newline="", errors="replace") as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        spectrum_name = (row.get("Spectrum") or "").strip()
        if not spectrum_name.startswith(ACQ + "."):
            continue

        match = re.match(
            rf"^{re.escape(ACQ)}\.(\d+)\.(\d+)\.(\d+)$",
            spectrum_name,
        )
        rt = number(row.get("Retention"))

        if rt is None or abs(rt - TARGET_RT) > RT_WINDOW:
            continue

        psm_scan = int(match.group(1)) if match else None

        nearest = min(
            native_ms2,
            key=lambda s: abs(s["native_rt"] - rt),
        ) if native_ms2 else None

        psm_mz = number(row.get("Observed M/Z"))
        rt_difference = (
            nearest["native_rt"] - rt if nearest else None
        )

        psms.append({
            "Spectrum": spectrum_name,
            "FragPipe scan": psm_scan,
            "PSM RT (s)": rt,
            "Nearest native scan": (
                nearest["native_scan"] if nearest else ""
            ),
            "Native spectrum ID": (
                nearest["native_id"] if nearest else ""
            ),
            "Native RT (s)": (
                nearest["native_rt"] if nearest else ""
            ),
            "Native minus PSM RT (s)": (
                round(rt_difference, 6)
                if rt_difference is not None else ""
            ),
            "RT match within tolerance": (
                abs(rt_difference) <= MATCH_TOLERANCE
                if rt_difference is not None else False
            ),
            "PSM observed precursor m/z": (
                psm_mz if psm_mz is not None else ""
            ),
            "Native precursor m/z": (
                nearest["precursor_mz"]
                if nearest and nearest["precursor_mz"] is not None
                else ""
            ),
            "Native precursor charge": (
                nearest["precursor_charge"]
                if nearest and nearest["precursor_charge"] is not None
                else ""
            ),
            "Peptide": row.get("Peptide", ""),
            "Modified Peptide": row.get("Modified Peptide", ""),
        })

psms.sort(key=lambda row: abs(row["PSM RT (s)"] - TARGET_RT))

OUT.parent.mkdir(parents=True, exist_ok=True)
fields = list(psms[0].keys()) if psms else [
    "Spectrum", "FragPipe scan", "PSM RT (s)",
    "Nearest native scan", "Native spectrum ID", "Native RT (s)",
    "Native minus PSM RT (s)", "RT match within tolerance",
    "PSM observed precursor m/z", "Native precursor m/z",
    "Native precursor charge", "Peptide", "Modified Peptide",
]

with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(psms)

print("\nNative MS2 spectra in file:", len(native_ms2))
print("Nearby PSMs:", len(psms))
print("PSMs matching a native MS2 by RT tolerance:",
      sum(r["RT match within tolerance"] is True for r in psms))
print("Output:", OUT)

for row in psms:
    print(
        f"PSM scan={row['FragPipe scan']}, "
        f"RT={row['PSM RT (s)']:.5f}, "
        f"native scan={row['Nearest native scan']}, "
        f"native RT={row['Native RT (s)']}, "
        f"RT difference={row['Native minus PSM RT (s)']}, "
        f"PSM m/z={row['PSM observed precursor m/z']}, "
        f"native precursor m/z={row['Native precursor m/z']}, "
        f"peptide={row['Peptide']}"
    )
