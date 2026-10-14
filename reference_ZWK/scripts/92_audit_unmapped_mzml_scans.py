#!/usr/bin/env python3
import csv
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
FEATURE = "10275997325121874585"

TRACE = REPORTS / "full_trace_chromatography_v1.tsv"
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
OUT = REPORTS / "unmapped_mzml_scan_audit_v1.tsv"

TARGETS = {5813, 5826, 5941, 5954, 5957, 5970, 5981, 5988}
RT_TOLERANCE_SEC = 0.02

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def local_name(tag):
    return tag.rsplit("}", 1)[-1]

trace_rows = read_tsv(TRACE)
trace_info = {}
for r in trace_rows:
    if r["feature_id"] != FEATURE:
        continue
    idx = int(r["scan_index_zero_based"])
    if idx in TARGETS:
        trace_info.setdefault(idx, {
            "rt_sec": float(r["rt_sec"]),
            "trace_indices": set(),
        })
        trace_info[idx]["trace_indices"].add(int(r["trace_index"]))

# Parse mzML spectra in document order; the ordinal is zero-based.
results = []
ordinal = -1

for event, elem in ET.iterparse(MZML, events=("end",)):
    if local_name(elem.tag) != "spectrum":
        continue

    ordinal += 1
    if ordinal not in TARGETS:
        elem.clear()
        continue

    spectrum_id = elem.attrib.get("id", "")
    ms_level = ""
    rt_sec = ""

    for node in elem.iter():
        name = local_name(node.tag)
        if name == "cvParam":
            accession = node.attrib.get("accession", "")
            if accession == "MS:1000511":
                ms_level = node.attrib.get("value", "")
            elif accession == "MS:1000016":
                value = float(node.attrib.get("value", "nan"))
                unit = node.attrib.get("unitAccession", "")
                # mzML retention time is usually seconds or minutes.
                rt_sec = value * 60.0 if unit == "UO:0000031" else value

    info = trace_info.get(ordinal, {})
    trace_rt = info.get("rt_sec", float("nan"))
    mzml_rt = float(rt_sec) if rt_sec != "" else float("nan")

    results.append({
        "scan_index_zero_based": ordinal,
        "trace_rt_sec": f"{trace_rt:.6f}" if trace_rt == trace_rt else "",
        "mzml_rt_sec": f"{mzml_rt:.6f}" if mzml_rt == mzml_rt else "",
        "rt_difference_sec": (
            f"{mzml_rt - trace_rt:.6f}"
            if mzml_rt == mzml_rt and trace_rt == trace_rt else ""
        ),
        "trace_indices": ",".join(
            str(x) for x in sorted(info.get("trace_indices", set()))
        ),
        "mzml_ms_level": ms_level,
        "mzml_spectrum_id": spectrum_id,
        "rt_match_within_tolerance": (
            abs(mzml_rt - trace_rt) <= RT_TOLERANCE_SEC
            if mzml_rt == mzml_rt and trace_rt == trace_rt else ""
        ),
    })
    elem.clear()

found = {int(r["scan_index_zero_based"]) for r in results}
for idx in sorted(TARGETS - found):
    info = trace_info.get(idx, {})
    rt = info.get("rt_sec", float("nan"))
    results.append({
        "scan_index_zero_based": idx,
        "trace_rt_sec": f"{rt:.6f}" if rt == rt else "",
        "mzml_rt_sec": "",
        "rt_difference_sec": "",
        "trace_indices": ",".join(
            str(x) for x in sorted(info.get("trace_indices", set()))
        ),
        "mzml_ms_level": "",
        "mzml_spectrum_id": "",
        "rt_match_within_tolerance": False,
    })

results.sort(key=lambda r: int(r["scan_index_zero_based"]))

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print(f"mzML: {MZML}")
print(f"Report: {OUT}")
print(f"Target indices: {len(TARGETS)}")
print(f"Spectra found: {len(found)}")
print("\nResults:")
for r in results:
    print(
        f'index={r["scan_index_zero_based"]} '
        f'ms_level={r["mzml_ms_level"]} '
        f'trace_RT={r["trace_rt_sec"]} '
        f'mzML_RT={r["mzml_rt_sec"]} '
        f'RT_match={r["rt_match_within_tolerance"]} '
        f'id={r["mzml_spectrum_id"]}'
    )
