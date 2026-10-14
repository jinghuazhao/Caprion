#!/usr/bin/env python3
"""Audit OpenMS spectrum indices against native mzML spectrum IDs."""

from pathlib import Path
import re
import xml.etree.ElementTree as ET
import pandas as pd

ROOT = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
OUT = ROOT / "reports"
FEATURE_ID = "10275997325121874585"
TARGET_INDICES = [5852, 5865, 5878, 5891, 5904, 5917]

NS = "{http://psi.hupo.org/ms/mzml}"


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    traces = pd.read_csv(TRACE_FILE, sep="\t", low_memory=False)
    traces = traces.loc[
        traces["feature_id"].astype(str).eq(FEATURE_ID)
        & pd.to_numeric(traces["scan_index"], errors="coerce").isin(
            TARGET_INDICES
        )
    ].copy()

    trace_rt = {}
    for _, row in traces.iterrows():
        scan = int(row["scan_index"])
        rt = pd.to_numeric(row.get("rt_sec"), errors="coerce")
        if pd.notna(rt):
            trace_rt.setdefault(scan, []).append(float(rt))

    rows = []
    with MZML.open("rb") as handle:
        ordinal = -1
        for _, elem in ET.iterparse(handle, events=("end",)):
            if elem.tag == NS + "spectrum":
                ordinal += 1
                sid = elem.attrib.get("id", "")
                match = re.search(r"(?:^|[ ;])scan=(\d+)", sid)
                native_scan = int(match.group(1)) if match else None

                ms_level = None
                for cv in elem.findall(NS + "cvParam"):
                    if cv.attrib.get("accession") == "MS:1000511":
                        ms_level = int(float(cv.attrib["value"]))
                        break

                rt_sec = None
                for cv in elem.findall(
                    ".//" + NS + "scan/" + NS + "cvParam"
                ):
                    if cv.attrib.get("accession") == "MS:1000016":
                        rt = float(cv.attrib["value"])
                        rt_sec = (
                            rt * 60
                            if cv.attrib.get("unitAccession") == "UO:0000031"
                            else rt
                        )
                        break

                if (
                    ordinal in TARGET_INDICES
                    or native_scan in TARGET_INDICES
                    or (
                        rt_sec is not None
                        and any(
                            abs(rt_sec - target_rt) <= 2.0
                            for values in trace_rt.values()
                            for target_rt in values
                        )
                    )
                ):
                    rows.append({
                        "spectrum_ordinal_zero_based": ordinal,
                        "spectrum_ordinal_one_based": ordinal + 1,
                        "native_spectrum_id": sid,
                        "native_scan_number": native_scan,
                        "ms_level": ms_level,
                        "rt_sec": rt_sec,
                    })

                elem.clear()

    result = pd.DataFrame(rows)
    output = OUT / "native_spectrum_index_mapping_v2.tsv"
    result.to_csv(output, sep="\t", index=False)

    print(f"mzML spectra indexed: {ordinal + 1}")
    print(f"Mapping report: {output}")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
