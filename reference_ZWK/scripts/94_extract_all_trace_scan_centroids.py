#!/usr/bin/env python3
"""Extract native MS1 centroid peaks for every scan in the candidate feature traces."""

import base64
import csv
import math
import re
import struct
import zlib
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
TRACE_FILE = REPORTS / "full_trace_chromatography_v1.tsv"
MZML_FILE = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)

OUT_CENTROIDS = REPORTS / "feature10275997325121874585_all_trace_scan_centroids_v1.tsv"
OUT_MATCHES = REPORTS / "feature10275997325121874585_trace_centroid_matches_v1.tsv"

FEATURE_ID = "10275997325121874585"
MZ_PADDING = 0.20
MATCH_TOLERANCE_DA = 0.015

# mzML CV accessions
MS_LEVEL = "MS:1000511"
SCAN_START_TIME = "MS:1000016"
MZ_ARRAY = "MS:1000514"
INTENSITY_ARRAY = "MS:1000515"
FLOAT32 = "MS:1000521"
FLOAT64 = "MS:1000523"
ZLIB_COMPRESSION = "MS:1000574"
NO_COMPRESSION = "MS:1000576"
LITTLE_ENDIAN = "MS:1000141"
BIG_ENDIAN = "MS:1000140"
MINUTES = "UO:0000031"
SECONDS = "UO:0000010"


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


def child_elements(element, name):
    return [child for child in element.iter() if local_name(child.tag) == name]


def cv_params(element):
    result = {}
    for child in element.iter():
        if local_name(child.tag) == "cvParam":
            result[child.attrib.get("accession", "")] = child.attrib
    return result


def parse_scan_number(spectrum_id):
    match = re.search(r"(?:^|\s)scan=(\d+)", spectrum_id)
    return int(match.group(1)) if match else None


def decode_binary_array(binary_array):
    params = cv_params(binary_array)

    if MZ_ARRAY in params:
        array_type = "mz"
    elif INTENSITY_ARRAY in params:
        array_type = "intensity"
    else:
        return None, None

    is_float32 = FLOAT32 in params
    is_float64 = FLOAT64 in params
    if is_float32 == is_float64:
        raise ValueError(f"Cannot determine binary numeric type for {array_type} array")

    endian = ">" if BIG_ENDIAN in params else "<"
    code = "f" if is_float32 else "d"
    width = 4 if is_float32 else 8

    binary_text = ""
    for child in binary_array.iter():
        if local_name(child.tag) == "binary":
            binary_text = (child.text or "").strip()
            break

    if not binary_text:
        return array_type, []

    raw = base64.b64decode(binary_text)
    if ZLIB_COMPRESSION in params:
        raw = zlib.decompress(raw)
    elif NO_COMPRESSION not in params:
        # Many mzML files omit an explicit no-compression CV term.
        # In that case, treat the payload as uncompressed if its size is valid.
        if len(raw) % width:
            raise ValueError(
                f"Unrecognised compression for {array_type} array; "
                f"decoded byte count={len(raw)}"
            )

    if len(raw) % width:
        raise ValueError(
            f"Binary array length {len(raw)} is not divisible by {width}"
        )

    count = len(raw) // width
    return array_type, list(struct.unpack(endian + code * count, raw))


def spectrum_arrays(spectrum):
    arrays = {}
    for binary_array in child_elements(spectrum, "binaryDataArray"):
        array_type, values = decode_binary_array(binary_array)
        if array_type:
            arrays[array_type] = values

    mzs = arrays.get("mz", [])
    intensities = arrays.get("intensity", [])
    if len(mzs) != len(intensities):
        raise ValueError(
            f"m/z and intensity array lengths differ: {len(mzs)} vs {len(intensities)}"
        )
    return mzs, intensities


def get_scan_rt_seconds(spectrum):
    for scan in child_elements(spectrum, "scan"):
        params = cv_params(scan)
        if SCAN_START_TIME in params:
            value = float(params[SCAN_START_TIME]["value"])
            unit = params[SCAN_START_TIME].get("unitAccession", "")
            if unit == MINUTES:
                return value * 60.0
            if unit == SECONDS:
                return value
            # mzML scan-start-time is usually minutes; preserve unknown units
            # as recorded rather than silently converting.
            return value
    return math.nan


def choose_column(fieldnames, candidates, required=True):
    normalised = {name.strip().lower(): name for name in fieldnames}
    for candidate in candidates:
        if candidate.lower() in normalised:
            return normalised[candidate.lower()]
    if required:
        raise KeyError(
            f"Could not find any of {candidates}. Available columns: {fieldnames}"
        )
    return None


def main():
    if not TRACE_FILE.exists():
        raise FileNotFoundError(TRACE_FILE)
    if not MZML_FILE.exists():
        raise FileNotFoundError(MZML_FILE)

    with TRACE_FILE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    if not rows:
        raise RuntimeError(f"No rows in {TRACE_FILE}")

    feature_col = choose_column(
        fieldnames, ["feature_id", "featureid", "feature"]
    )
    scan_col = choose_column(
        fieldnames,
        ["scan_index_zero_based", "spectrum_ordinal_zero_based",
         "scan_index", "spectrum_index"],
    )
    trace_mz_col = choose_column(
        fieldnames,
        ["observed_mz", "mz", "trace_mz", "centroid_mz"],
    )
    trace_index_col = choose_column(
        fieldnames, ["trace_index", "isotope_index", "peak_index"], required=False
    )
    trace_rt_col = choose_column(
        fieldnames,
        ["rt_seconds", "retention_time_seconds", "rt_sec", "rt"],
        required=False,
    )

    feature_rows = [
        row for row in rows
        if row[feature_col].strip() == FEATURE_ID
    ]
    if not feature_rows:
        raise RuntimeError(
            f"No rows for feature {FEATURE_ID} in {TRACE_FILE}"
        )

    trace_scans = defaultdict(list)
    for row in feature_rows:
        ordinal = int(float(row[scan_col]))
        trace_mz = float(row[trace_mz_col])
        trace_index = (
            row[trace_index_col].strip() if trace_index_col else ""
        )
        trace_rt = (
            float(row[trace_rt_col])
            if trace_rt_col and row[trace_rt_col].strip()
            else math.nan
        )
        trace_scans[ordinal].append({
            "trace_mz": trace_mz,
            "trace_index": trace_index,
            "trace_rt": trace_rt,
        })

    min_mz = min(
        item["trace_mz"]
        for scan_rows in trace_scans.values()
        for item in scan_rows
    ) - MZ_PADDING
    max_mz = max(
        item["trace_mz"]
        for scan_rows in trace_scans.values()
        for item in scan_rows
    ) + MZ_PADDING

    print(f"Feature: {FEATURE_ID}")
    print(f"Trace observations: {len(feature_rows)}")
    print(f"Unique trace scans: {len(trace_scans)}")
    print(f"Centroid extraction range: {min_mz:.6f}–{max_mz:.6f} m/z")
    print("Reading original mzML; this may take a little while...")

    centroid_rows = []
    trace_match_rows = []
    ordinal = -1
    target_ordinals = set(trace_scans)

    # Stream spectrum elements to avoid loading the complete mzML tree.
    for event, spectrum in ET.iterparse(MZML_FILE, events=("end",)):
        if local_name(spectrum.tag) != "spectrum":
            continue

        ordinal += 1
        if ordinal not in target_ordinals:
            spectrum.clear()
            continue

        params = cv_params(spectrum)
        level = int(float(params[MS_LEVEL]["value"])) if MS_LEVEL in params else None
        spectrum_id = spectrum.attrib.get("id", "")
        native_scan = parse_scan_number(spectrum_id)
        rt_seconds = get_scan_rt_seconds(spectrum)

        if level != 1:
            raise RuntimeError(
                f"Trace ordinal {ordinal} maps to MS level {level}, not MS1"
            )

        mzs, intensities = spectrum_arrays(spectrum)
        local_centroids = []
        for peak_index, (mz, intensity) in enumerate(zip(mzs, intensities)):
            if min_mz <= mz <= max_mz:
                local_centroids.append({
                    "ordinal": ordinal,
                    "native_scan": native_scan,
                    "spectrum_id": spectrum_id,
                    "rt_seconds": rt_seconds,
                    "peak_index": peak_index,
                    "mz": float(mz),
                    "intensity": float(intensity),
                })

        centroid_rows.extend(local_centroids)

        for trace in trace_scans[ordinal]:
            if local_centroids:
                best = min(
                    local_centroids,
                    key=lambda peak: abs(peak["mz"] - trace["trace_mz"])
                )
                error_mda = (best["mz"] - trace["trace_mz"]) * 1000.0
                within_tolerance = abs(best["mz"] - trace["trace_mz"]) <= MATCH_TOLERANCE_DA
            else:
                best = None
                error_mda = math.nan
                within_tolerance = False

            trace_match_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index_zero_based": ordinal,
                "native_scan_number": native_scan,
                "spectrum_id": spectrum_id,
                "ms_level": level,
                "spectrum_rt_seconds": rt_seconds,
                "trace_rt_seconds": trace["trace_rt"],
                "trace_index": trace["trace_index"],
                "trace_mz": trace["trace_mz"],
                "nearest_centroid_mz": best["mz"] if best else math.nan,
                "nearest_centroid_intensity": best["intensity"] if best else math.nan,
                "nearest_centroid_error_mDa": error_mda,
                "within_0_015_Da": within_tolerance,
                "local_centroid_count": len(local_centroids),
            })

        spectrum.clear()

    missing = target_ordinals - {
        row["scan_index_zero_based"] for row in trace_match_rows
    }
    if missing:
        raise RuntimeError(
            f"Failed to process trace scan ordinals: {sorted(missing)}"
        )

    with OUT_CENTROIDS.open("w", newline="") as handle:
        columns = [
            "feature_id", "scan_index_zero_based", "native_scan_number",
            "spectrum_id", "rt_seconds", "peak_index", "centroid_mz",
            "centroid_intensity",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for peak in centroid_rows:
            writer.writerow({
                "feature_id": FEATURE_ID,
                "scan_index_zero_based": peak["ordinal"],
                "native_scan_number": peak["native_scan"],
                "spectrum_id": peak["spectrum_id"],
                "rt_seconds": f'{peak["rt_seconds"]:.6f}',
                "peak_index": peak["peak_index"],
                "centroid_mz": f'{peak["mz"]:.9f}',
                "centroid_intensity": f'{peak["intensity"]:.6f}',
            })

    with OUT_MATCHES.open("w", newline="") as handle:
        columns = [
            "feature_id", "scan_index_zero_based", "native_scan_number",
            "spectrum_id", "ms_level", "spectrum_rt_seconds",
            "trace_rt_seconds", "trace_index", "trace_mz",
            "nearest_centroid_mz", "nearest_centroid_intensity",
            "nearest_centroid_error_mDa", "within_0_015_Da",
            "local_centroid_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in trace_match_rows:
            writer.writerow({
                key: (
                    f"{value:.9f}"
                    if isinstance(value, float) and math.isfinite(value)
                    and key in {"trace_mz", "nearest_centroid_mz"}
                    else f"{value:.6f}"
                    if isinstance(value, float) and math.isfinite(value)
                    else value
                )
                for key, value in row.items()
            })

    matched = sum(row["within_0_015_Da"] for row in trace_match_rows)
    print(f"Native centroid export: {OUT_CENTROIDS}")
    print(f"Trace-to-centroid match report: {OUT_MATCHES}")
    print(f"Centroid rows exported: {len(centroid_rows)}")
    print(f"Trace rows matched within ±0.015 Da: {matched}/{len(trace_match_rows)}")
    print(
        "Unique scans with centroid data:",
        len({row["ordinal"] for row in centroid_rows}),
        "/",
        len(trace_scans),
    )
    print("Note: centroid proximity validates peak correspondence, not isotope identity.")


if __name__ == "__main__":
    main()
