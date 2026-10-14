#!/usr/bin/env python3
"""Inspect native MS1 centroids around a candidate charge-2 envelope."""

from pathlib import Path
import base64
import csv
import math
import statistics
import struct
import xml.etree.ElementTree as ET
import zlib

ROOT = Path.cwd()
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
MZML_FILE = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
REPORTS = ROOT / "reports"
REPORTS.mkdir(exist_ok=True)

FEATURE_ID = "10275997325121874585"
SCANS = [5852, 5865, 5878, 5891, 5904, 5917]
CORE_TRACES = [1, 2, 3, 4]
CANDIDATE_TRACES = [5, 6]
SPACING = 1.00335483507 / 2.0

# Inspect a local region around the envelope. Keep broad enough to see
# nearby competing peaks, but small enough to remain interpretable.
MZ_MARGIN = 1.2
RT_MARGIN_SEC = 2.0
CORE_MATCH_TOL_MDA = 3.0


def lname(tag):
    return tag.rsplit("}", 1)[-1]


def iter_named(element, name):
    for child in element.iter():
        if lname(child.tag) == name:
            yield child


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0].keys()), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_traces():
    traces = {}
    with TRACE_FILE.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["feature_id"] != FEATURE_ID:
                continue
            scan = int(row["scan_index"])
            idx = int(row["trace_index"])
            if scan not in SCANS:
                continue
            traces[(scan, idx)] = {
                "mz": float(row["observed_mz"]),
                "rt": float(row["rt_sec"]),
                "intensity": float(row["peak_intensity"]),
            }
    return traces


def decode_array(array):
    params = {
        p.attrib.get("accession"): p.attrib
        for p in iter_named(array, "cvParam")
    }
    encoded = "".join(
        (node.text or "").strip()
        for node in iter_named(array, "binary")
    )
    if not encoded:
        return []

    raw = base64.b64decode(encoded)
    if "MS:1000574" in params:
        raw = zlib.decompress(raw)

    if "MS:1000523" in params:
        width, fmt = 8, "d"
    elif "MS:1000521" in params:
        width, fmt = 4, "f"
    else:
        raise ValueError("Unknown mzML numeric-array precision")

    if len(raw) % width:
        raise ValueError("Invalid binary array length")

    n = len(raw) // width
    return struct.unpack("<" + fmt * n, raw)


def parse_spectrum(spectrum):
    params = {
        p.attrib.get("accession"): p.attrib
        for p in iter_named(spectrum, "cvParam")
    }
    ms_level = int(params.get("MS:1000511", {}).get("value", "0"))
    rt_value = params.get("MS:1000016", {}).get("value")
    rt_unit = params.get("MS:1000016", {}).get("unitAccession", "")
    rt_sec = float(rt_value) if rt_value is not None else math.nan
    if rt_unit == "UO:0000031":
        rt_sec *= 60.0

    mzs, intensities = [], []
    for array in iter_named(spectrum, "binaryDataArray"):
        cv = {
            p.attrib.get("accession"): p.attrib
            for p in iter_named(array, "cvParam")
        }
        values = decode_array(array)
        if "MS:1000514" in cv:
            mzs = values
        elif "MS:1000515" in cv:
            intensities = values

    if len(mzs) != len(intensities):
        raise ValueError("m/z and intensity array lengths differ")

    return ms_level, rt_sec, list(zip(mzs, intensities))


def main():
    if not TRACE_FILE.exists():
        raise FileNotFoundError(TRACE_FILE)
    if not MZML_FILE.exists():
        raise FileNotFoundError(MZML_FILE)

    traces = read_traces()
    native = {}
    ordinal = -1

    # Trace scan_index is a zero-based spectrum ordinal, not the native
    # scan number. Read through the last requested ordinal.
    for _, element in ET.iterparse(MZML_FILE, events=("end",)):
        if lname(element.tag) != "spectrum":
            continue
        ordinal += 1
        if ordinal in SCANS:
            ms_level, rt_sec, peaks = parse_spectrum(element)
            if ms_level != 1:
                raise RuntimeError(
                    f"Expected MS1 at ordinal {ordinal}, got MS{ms_level}"
                )
            native[ordinal] = {
                "rt": rt_sec,
                "native_scan": ordinal + 1,
                "spectrum_id": element.attrib.get("id", ""),
                "peaks": peaks,
            }
        element.clear()
        if ordinal >= max(SCANS):
            break

    if len(native) != len(SCANS):
        raise RuntimeError(f"Only found {len(native)} of {len(SCANS)} spectra")

    all_rows = []
    summary_rows = []

    for scan in SCANS:
        spectrum = native[scan]
        core = [
            traces[(scan, idx)]
            for idx in CORE_TRACES
            if (scan, idx) in traces
        ]
        if len(core) != 4:
            raise RuntimeError(f"Missing one or more core traces at scan {scan}")

        intercept = statistics.mean(
            traces[(scan, idx)]["mz"] - (idx - 1) * SPACING
            for idx in CORE_TRACES
        )
        predicted = {
            idx: intercept + (idx - 1) * SPACING
            for idx in range(1, 7)
        }

        rt_values = [
            traces[(scan, idx)]["rt"]
            for idx in CORE_TRACES
            if (scan, idx) in traces
        ]
        center_rt = statistics.mean(rt_values)
        n_local = 0

        for peak_no, (mz, intensity) in enumerate(spectrum["peaks"]):
            if not math.isfinite(mz) or not math.isfinite(intensity):
                continue
            if intensity <= 0:
                continue
            if not (intercept - MZ_MARGIN <= mz <= intercept + 6 * SPACING + MZ_MARGIN):
                continue

            # Candidate trace labels are nearest-m/z diagnostics, not
            # definitive isotope identities.
            nearest_idx = min(predicted, key=lambda i: abs(mz - predicted[i]))
            nearest_error = (mz - predicted[nearest_idx]) * 1000.0
            trace_near = [
                idx for idx in (5, 6)
                if (scan, idx) in traces
                and abs(mz - traces[(scan, idx)]["mz"]) <= 0.01
            ]

            row = {
                "feature_id": FEATURE_ID,
                "scan_index_zero_based": scan,
                "native_scan_number": spectrum["native_scan"],
                "spectrum_id": spectrum["spectrum_id"],
                "ms_level": 1,
                "rt_sec": f"{spectrum['rt']:.5f}",
                "centroid_ordinal_in_spectrum": peak_no,
                "centroid_mz": f"{mz:.9f}",
                "centroid_intensity": f"{intensity:.3f}",
                "core_intercept_mz": f"{intercept:.9f}",
                "nearest_sequential_isotope_index": nearest_idx - 1,
                "nearest_model_mz": f"{predicted[nearest_idx]:.9f}",
                "error_to_nearest_model_mDa": f"{nearest_error:.6f}",
                "within_3mDa_of_nearest_model": abs(nearest_error) <= CORE_MATCH_TOL_MDA,
                "near_trace5_or_trace6_observation": ",".join(map(str, trace_near)),
                "rt_distance_from_core_sec": f"{spectrum['rt'] - center_rt:.5f}",
            }
            all_rows.append(row)
            n_local += 1

        summary_rows.append({
            "feature_id": FEATURE_ID,
            "scan_index_zero_based": scan,
            "native_scan_number": spectrum["native_scan"],
            "rt_sec": f"{spectrum['rt']:.5f}",
            "core_intercept_mz": f"{intercept:.9f}",
            "native_centroids_in_local_mz_region": n_local,
            "trace5_mz": (
                f"{traces[(scan, 5)]['mz']:.9f}"
                if (scan, 5) in traces else ""
            ),
            "trace6_mz": (
                f"{traces[(scan, 6)]['mz']:.9f}"
                if (scan, 6) in traces else ""
            ),
        })

    out_all = REPORTS / "local_ms1_centroids_feature10275997325121874585_v1.tsv"
    out_summary = REPORTS / "local_ms1_centroids_feature10275997325121874585_v1_summary.tsv"
    write_tsv(out_all, all_rows)
    write_tsv(out_summary, summary_rows)

    print(f"Feature: {FEATURE_ID}")
    print(f"Selected MS1 spectra: {len(native)}")
    print(f"Local m/z margin: ±{MZ_MARGIN:.2f} Da around the six-isotope prediction range")
    print(f"Centroid detail report: {out_all}")
    print(f"Per-scan summary: {out_summary}")
    print("\nLocal centroid counts:")
    for row in summary_rows:
        print(
            f"scan {row['scan_index_zero_based']} "
            f"(native {row['native_scan_number']}): "
            f"{row['native_centroids_in_local_mz_region']} centroids"
        )
    print(
        "\nCaution: proximity to a theoretical isotope position is not proof "
        "of isotope-group membership. Co-elution, isotope pattern, charge, "
        "and consistency across scans must be assessed together."
    )


if __name__ == "__main__":
    main()
