#!/usr/bin/env python3
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
TRACE = REPORTS / "full_trace_chromatography_v1.tsv"
OUT = REPORTS / "native_spectrum_index_mapping_complete_v1.tsv"
TRACE_AUDIT = REPORTS / "trace_scan_mapping_complete_audit_v1.tsv"

def local_name(tag):
    return tag.rsplit("}", 1)[-1]

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

mapping = []
ordinal = -1

for event, elem in ET.iterparse(MZML, events=("end",)):
    if local_name(elem.tag) != "spectrum":
        continue

    ordinal += 1
    spectrum_id = elem.attrib.get("id", "")
    ms_level = ""
    rt_sec = ""

    for node in elem.iter():
        if local_name(node.tag) != "cvParam":
            continue
        accession = node.attrib.get("accession", "")
        if accession == "MS:1000511":
            ms_level = node.attrib.get("value", "")
        elif accession == "MS:1000016":
            value = float(node.attrib.get("value", "nan"))
            unit = node.attrib.get("unitAccession", "")
            rt_sec = value * 60.0 if unit == "UO:0000031" else value

    # Extract the native scan number from the mzML spectrum ID.
    native_scan = ""
    if "scan=" in spectrum_id:
        try:
            native_scan = int(spectrum_id.rsplit("scan=", 1)[1].split()[0])
        except ValueError:
            native_scan = ""

    mapping.append({
        "spectrum_ordinal_zero_based": ordinal,
        "spectrum_ordinal_one_based": ordinal + 1,
        "native_spectrum_id": spectrum_id,
        "native_scan_number": native_scan,
        "ms_level": ms_level,
        "rt_sec": f"{float(rt_sec):.6f}" if rt_sec != "" else "",
    })
    elem.clear()

with OUT.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(mapping[0]), delimiter="\t")
    w.writeheader()
    w.writerows(mapping)

# Audit every trace observation against the new mapping.
trace_rows = read_tsv(TRACE)
by_index = {int(r["spectrum_ordinal_zero_based"]): r for r in mapping}
audit = []

for r in trace_rows:
    idx = int(r["scan_index_zero_based"])
    m = by_index.get(idx)
    trace_rt = float(r["rt_sec"])
    mzml_rt = float(m["rt_sec"]) if m and m["rt_sec"] else float("nan")
    audit.append({
        "feature_id": r["feature_id"],
        "scan_index_zero_based": idx,
        "trace_index": r["trace_index"],
        "trace_rt_sec": f"{trace_rt:.6f}",
        "mapped_native_scan_number": m["native_scan_number"] if m else "",
        "mapped_ms_level": m["ms_level"] if m else "",
        "mapped_rt_sec": m["rt_sec"] if m else "",
        "rt_difference_sec": (
            f"{mzml_rt - trace_rt:.6f}"
            if mzml_rt == mzml_rt else ""
        ),
        "rt_match_within_0.02_sec": (
            abs(mzml_rt - trace_rt) <= 0.02
            if mzml_rt == mzml_rt else False
        ),
    })

with TRACE_AUDIT.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(audit[0]), delimiter="\t")
    w.writeheader()
    w.writerows(audit)

print(f"Complete mapping: {OUT}")
print(f"Trace mapping audit: {TRACE_AUDIT}")
print(f"Spectra mapped: {len(mapping)}")
print(f"Trace observations audited: {len(audit)}")
print(
    "Trace observations with retention-time match:",
    sum(str(r["rt_match_within_0.02_sec"]).lower() == "true" for r in audit),
    "/",
    len(audit),
)
print(
    "Trace observations mapped to MS1:",
    sum(str(r["mapped_ms_level"]) == "1" for r in audit),
    "/",
    len(audit),
)
