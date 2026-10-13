#!/usr/bin/env python3

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

if len(sys.argv) != 2 or sys.argv[1] not in {
    "xms2", "xms3", "xms4", "xms5"
}:
    sys.exit(f"Usage: {sys.argv[0]} xms2|xms3|xms4|xms5")

suffix = sys.argv[1]
run = "szwk021704i19101"

root = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
mzml = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    f"analysis/ZWK/mzML/{run}{suffix}.mzML"
)
psm_file = root / "search_98" / run / "psm.tsv"
report = root / "reports" / f"spectrum_scan_mapping_{suffix}.tsv"

if not mzml.is_file():
    sys.exit(f"mzML file not found: {mzml}")
if not psm_file.is_file():
    sys.exit(f"PSM file not found: {psm_file}")

# Read native MS2 scan numbers and retention times from mzML.
ns = {"mzml": "http://psi.hupo.org/ms/mzml"}
native = {}

for _, elem in ET.iterparse(mzml, events=("end",)):
    if elem.tag.endswith("}spectrum"):
        params = {
            p.attrib.get("accession"): p.attrib
            for p in elem.iter()
            if p.tag.endswith("}cvParam")
        }

        ms_level = params.get("MS:1000511", {}).get("value")
        if ms_level == "2":
            match = re.search(
                r"(?:^| )scan=(\d+)",
                elem.attrib.get("id", "")
            )
            if match:
                scan = int(match.group(1))
                rt = None

                for p in elem.iter():
                    if (
                        p.tag.endswith("}cvParam")
                        and p.attrib.get("accession") == "MS:1000016"
                    ):
                        rt = float(p.attrib["value"])
                        unit = p.attrib.get("unitAccession")
                        if unit == "UO:0000031":  # minute
                            rt *= 60
                        elif unit == "UO:0000010":  # second
                            pass
                        else:
                            raise ValueError(
                                f"Unrecognised RT unit for scan {scan}: "
                                f"{p.attrib}"
                            )
                        break

                if rt is not None:
                    native[scan] = rt

        elem.clear()

# Read PSMs and match spectrum scan numbers to native MS2 scans.
psms = pd.read_csv(psm_file, sep="\t")
psms = psms[
    psms["Spectrum File"].str.contains(
        f"{run}{suffix}.pep.xml", regex=False, na=False
    )
].copy()

print(f"PSMs selected for {suffix}: {len(psms)}")
pattern = re.compile(
    rf"{re.escape(suffix)}\.(\d+)\.(\d+)\.(\d+)$"
)

records = []
for _, row in psms.iterrows():
    spectrum = str(row["Spectrum"])
    match = pattern.search(spectrum)

    scan1 = int(match.group(1)) if match else None
    scan2 = int(match.group(2)) if match else None
    name_charge = int(match.group(3)) if match else None
    scan_found = scan1 in native if scan1 is not None else False
    native_rt = native.get(scan1) if scan_found else None

    fragpipe_rt_sec = (
        float(row["Retention"]) if pd.notna(row["Retention"]) else None
    )
    rt_difference_sec = (
        fragpipe_rt_sec - native_rt
        if fragpipe_rt_sec is not None and native_rt is not None
        else None
    )

    charge = int(row["Charge"]) if pd.notna(row["Charge"]) else None

    records.append({
        "Spectrum": spectrum,
        "Scan_1": scan1,
        "Scan_2": scan2,
        "Native_MS2_scan_found": scan_found,
        "PSM_charge": charge,
        "Spectrum_name_charge": name_charge,
        "Charge_agrees": charge == name_charge,
        "FragPipe_RT_sec": fragpipe_rt_sec,
        "Native_RT_sec": native_rt,
        "RT_difference_sec": rt_difference_sec,
    })

result = pd.DataFrame(records)
report.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(report, sep="\t", index=False)

n = len(result)
found = int(result["Native_MS2_scan_found"].sum())
charge_ok = int(result["Charge_agrees"].sum())
rt = result["RT_difference_sec"].dropna()

print(f"Acquisition: {run}{suffix}")
print(f"PSMs: {n}")
print(f"Native MS2 scan IDs found: {found}/{n}")
print(f"Charge agrees: {charge_ok}/{n}")

if not rt.empty:
    print("\nRetention-time difference (seconds):")
    print(rt.describe(percentiles=[0.5, 0.9, 0.95, 0.99]).to_string())
    print(f"\nAbsolute difference <= 0.1 sec: {(rt.abs() <= 0.1).sum()}")

print(f"\nSaved: {report}")
