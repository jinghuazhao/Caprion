#!/usr/bin/env python3
"""Inspect traces 5 and 6 in native MS1 spectra by validated spectrum ordinal.

The trace scan_index is treated as a zero-based mzML spectrum ordinal, not
the native scan number. No original XLSX isotope-group assignments are used.
"""

from pathlib import Path
import base64
import re
import zlib
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

ROOT = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
MZML_FILE = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
OUT = ROOT / "reports"

FEATURE_ID = "10275997325121874585"
COMMON_SCANS = [5852, 5865, 5878, 5891, 5904, 5917]
TRACE_INDICES = [5, 6]
SPACING = 1.00335483507 / 2.0
WINDOW_DA = 0.15

NS = "{http://psi.hupo.org/ms/mzml}"

PER_SCAN = OUT / "native_centroids_traces5_6_v3_per_scan.tsv"
CANDIDATES = OUT / "native_centroids_traces5_6_v3_candidates.tsv"
SUMMARY = OUT / "native_centroids_traces5_6_v3_summary.tsv"


def decode_array(array_element):
    accessions = {
        cv.attrib.get("accession")
        for cv in array_element.findall(NS + "cvParam")
    }
    binary = array_element.find(NS + "binary")
    if binary is None or not binary.text:
        return np.array([], dtype=float)

    data = base64.b64decode(binary.text)
    if "MS:1000574" in accessions:
        data = zlib.decompress(data)

    if "MS:1000523" in accessions:
        dtype = "<f8"
    elif "MS:1000521" in accessions:
        dtype = "<f4"
    else:
        raise ValueError("Unsupported mzML binary-array precision")

    return np.frombuffer(data, dtype=dtype).astype(float)


def parse_spectrum(element, ordinal):
    ms_level = None
    for cv in element.findall(NS + "cvParam"):
        if cv.attrib.get("accession") == "MS:1000511":
            ms_level = int(float(cv.attrib["value"]))
            break

    rt_sec = np.nan
    for cv in element.findall(".//" + NS + "scan/" + NS + "cvParam"):
        if cv.attrib.get("accession") == "MS:1000016":
            rt = float(cv.attrib["value"])
            rt_sec = (
                rt * 60.0
                if cv.attrib.get("unitAccession") == "UO:0000031"
                else rt
            )
            break

    mz = None
    intensity = None
    for array in element.findall(
        ".//" + NS + "binaryDataArray"
    ):
        accessions = {
            cv.attrib.get("accession")
            for cv in array.findall(NS + "cvParam")
        }
        if "MS:1000514" in accessions:
            mz = decode_array(array)
        elif "MS:1000515" in accessions:
            intensity = decode_array(array)

    sid = element.attrib.get("id", "")
    match = re.search(r"scan=(\d+)", sid)
    native_scan = int(match.group(1)) if match else None

    return {
        "ordinal": ordinal,
        "spectrum_id": sid,
        "native_scan": native_scan,
        "ms_level": ms_level,
        "rt_sec": rt_sec,
        "mz": mz,
        "intensity": intensity,
    }


def nearest_centroid(mz, intensity, target, window):
    indices = np.flatnonzero(np.abs(mz - target) <= window)
    if len(indices) == 0:
        return None

    idx = indices[np.argmin(np.abs(mz[indices] - target))]
    return {
        "mz": float(mz[idx]),
        "intensity": float(intensity[idx]),
        "error_mDa": float((mz[idx] - target) * 1000.0),
        "n_candidates": int(len(indices)),
        "max_intensity": float(np.max(intensity[indices])),
        "sum_intensity": float(np.sum(intensity[indices])),
    }


def main():
    if not TRACE_FILE.exists():
        raise FileNotFoundError(TRACE_FILE)
    if not MZML_FILE.exists():
        raise FileNotFoundError(MZML_FILE)

    OUT.mkdir(parents=True, exist_ok=True)
    traces = pd.read_csv(TRACE_FILE, sep="\t", low_memory=False)

    required = {
        "feature_id", "trace_index", "scan_index",
        "observed_mz", "peak_intensity",
    }
    missing = required - set(traces.columns)
    if missing:
        raise ValueError(f"Missing trace columns: {sorted(missing)}")

    traces = traces.loc[
        traces["feature_id"].astype(str).eq(FEATURE_ID)
        & pd.to_numeric(traces["scan_index"], errors="coerce").isin(
            COMMON_SCANS
        )
    ].copy()

    for col in ("trace_index", "scan_index", "observed_mz", "peak_intensity"):
        traces[col] = pd.to_numeric(traces[col], errors="coerce")
    traces = traces.dropna(
        subset=["trace_index", "scan_index", "observed_mz", "peak_intensity"]
    )
    traces["trace_index"] = traces["trace_index"].astype(int)
    traces["scan_index"] = traces["scan_index"].astype(int)

    duplicate = traces.duplicated(["scan_index", "trace_index"], keep=False)
    if duplicate.any():
        raise ValueError("Duplicate scan/trace keys found.")

    trace_by_scan = {
        scan: group.set_index("trace_index")
        for scan, group in traces.groupby("scan_index")
    }

    # Capture only requested zero-based ordinals. The ordinal is incremented
    # for every mzML spectrum, including MS2 spectra.
    selected = {}
    wanted = set(COMMON_SCANS)
    ordinal = -1

    with MZML_FILE.open("rb") as handle:
        for _, element in ET.iterparse(handle, events=("end",)):
            if element.tag == NS + "spectrum":
                ordinal += 1
                if ordinal in wanted:
                    selected[ordinal] = parse_spectrum(element, ordinal)
                element.clear()

                if len(selected) == len(wanted):
                    break

    per_scan_rows = []
    candidate_rows = []

    for scan in COMMON_SCANS:
        spectrum = selected.get(scan)
        group = trace_by_scan.get(scan)

        if spectrum is None:
            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index_zero_based": scan,
                "status": "spectrum_ordinal_not_found",
            })
            continue

        if spectrum["ms_level"] != 1:
            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index_zero_based": scan,
                "spectrum_id": spectrum["spectrum_id"],
                "native_scan_number": spectrum["native_scan"],
                "ms_level": spectrum["ms_level"],
                "rt_sec": spectrum["rt_sec"],
                "status": "selected_spectrum_is_not_MS1",
            })
            continue

        if (
            group is None
            or not all(i in group.index for i in (1, 2, 3, 4))
        ):
            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index_zero_based": scan,
                "spectrum_id": spectrum["spectrum_id"],
                "native_scan_number": spectrum["native_scan"],
                "ms_level": spectrum["ms_level"],
                "rt_sec": spectrum["rt_sec"],
                "status": "missing_base_traces_1_to_4",
            })
            continue

        mz = spectrum["mz"]
        intensity = spectrum["intensity"]
        if mz is None or intensity is None:
            raise ValueError(f"Missing binary arrays at ordinal {scan}")

        base_mz = np.array(
            [float(group.loc[i, "observed_mz"]) for i in (1, 2, 3, 4)]
        )
        isotope_index = np.arange(4, dtype=float)
        intercept = float(np.mean(base_mz - isotope_index * SPACING))

        for trace_index in TRACE_INDICES:
            predicted = intercept + (trace_index - 1) * SPACING
            has_trace = trace_index in group.index

            observed = (
                float(group.loc[trace_index, "observed_mz"])
                if has_trace else np.nan
            )
            trace_intensity = (
                float(group.loc[trace_index, "peak_intensity"])
                if has_trace else np.nan
            )
            trace_residual = (
                (observed - predicted) * 1000.0
                if has_trace else np.nan
            )

            centroid = nearest_centroid(mz, intensity, predicted, WINDOW_DA)
            status = "centroid_found" if centroid else "no_centroid_in_window"

            if centroid:
                trace_to_centroid = (
                    (observed - centroid["mz"]) * 1000.0
                    if has_trace else np.nan
                )
                candidate_rows.append({
                    "feature_id": FEATURE_ID,
                    "scan_index_zero_based": scan,
                    "native_scan_number": spectrum["native_scan"],
                    "trace_index": trace_index,
                    "spectrum_id": spectrum["spectrum_id"],
                    "rt_sec": spectrum["rt_sec"],
                    "predicted_mz_from_traces1_4": predicted,
                    "centroid_mz": centroid["mz"],
                    "centroid_intensity": centroid["intensity"],
                    "centroid_error_vs_prediction_mDa": centroid["error_mDa"],
                    "candidate_count_in_window": centroid["n_candidates"],
                    "max_intensity_in_window": centroid["max_intensity"],
                    "sum_intensity_in_window": centroid["sum_intensity"],
                })
            else:
                trace_to_centroid = np.nan

            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index_zero_based": scan,
                "native_scan_number": spectrum["native_scan"],
                "spectrum_id": spectrum["spectrum_id"],
                "ms_level": spectrum["ms_level"],
                "rt_sec": spectrum["rt_sec"],
                "trace_index": trace_index,
                "trace_observation_present": has_trace,
                "trace_observed_mz": observed,
                "trace_peak_intensity": trace_intensity,
                "four_trace_intercept_mz": intercept,
                "expected_spacing_Da": SPACING,
                "predicted_mz_from_traces1_4": predicted,
                "trace_residual_vs_prediction_mDa": trace_residual,
                "native_centroid_status": status,
                "nearest_native_centroid_mz": (
                    centroid["mz"] if centroid else np.nan
                ),
                "nearest_native_centroid_intensity": (
                    centroid["intensity"] if centroid else np.nan
                ),
                "centroid_error_vs_prediction_mDa": (
                    centroid["error_mDa"] if centroid else np.nan
                ),
                "centroid_candidates_in_window": (
                    centroid["n_candidates"] if centroid else 0
                ),
                "trace_to_centroid_delta_mDa": trace_to_centroid,
            })

    per_scan_df = pd.DataFrame(per_scan_rows)
    candidate_df = pd.DataFrame(candidate_rows)

    per_scan_df.to_csv(PER_SCAN, sep="\t", index=False)
    candidate_df.to_csv(CANDIDATES, sep="\t", index=False)

    summary_rows = []
    if not per_scan_df.empty and "trace_index" in per_scan_df:
        for trace_index, part in per_scan_df.groupby("trace_index"):
            valid = part.loc[
                part["native_centroid_status"].eq("centroid_found")
            ]
            trace_res = pd.to_numeric(
                part["trace_residual_vs_prediction_mDa"], errors="coerce"
            ).dropna()
            centroid_res = pd.to_numeric(
                valid["centroid_error_vs_prediction_mDa"], errors="coerce"
            ).dropna()
            delta = pd.to_numeric(
                valid["trace_to_centroid_delta_mDa"], errors="coerce"
            ).dropna()

            summary_rows.append({
                "trace_index": int(trace_index),
                "n_scans_inspected": len(part),
                "n_trace_observations": int(
                    part["trace_observation_present"].sum()
                ),
                "n_native_centroids_found": len(valid),
                "median_abs_trace_residual_mDa": (
                    float(np.median(np.abs(trace_res))
                          ) if len(trace_res) else np.nan
                ),
                "max_abs_trace_residual_mDa": (
                    float(np.max(np.abs(trace_res)))
                    if len(trace_res) else np.nan
                ),
                "median_abs_centroid_error_mDa": (
                    float(np.median(np.abs(centroid_res)))
                    if len(centroid_res) else np.nan
                ),
                "max_abs_centroid_error_mDa": (
                    float(np.max(np.abs(centroid_res)))
                    if len(centroid_res) else np.nan
                ),
                "median_trace_to_centroid_delta_mDa": (
                    float(np.median(delta)) if len(delta) else np.nan
                ),
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SUMMARY, sep="\t", index=False)

    print(f"Feature: {FEATURE_ID}")
    print(f"Expected charge-2 spacing: {SPACING:.12f} Da")
    print("Mapping: zero-based mzML spectrum ordinal")
    print(f"MS1 spectra selected: {sum(s.get('ms_level') == 1 for s in selected.values())}")
    print(f"Per-scan report: {PER_SCAN}")
    print(f"Centroid report: {CANDIDATES}")
    print(f"Summary report: {SUMMARY}")
    print("\nPer-scan inspection:")
    print(per_scan_df.to_string(index=False))
    print("\nSummary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
