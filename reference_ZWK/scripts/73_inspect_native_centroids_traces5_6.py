#!/usr/bin/env python3
"""Inspect native MS1 centroids for candidate traces 5 and 6.

The four-trace envelope provides a theoretical-spacing reference. This script
compares traces 5 and 6 against native centroids near their predicted m/z.

Diagnostic only: it does not assign isotope groups and does not use the
original XLSX isotope-group assignments.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import gzip
import struct
import base64
import re

ROOT = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
MZML_FILE = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
OUT = ROOT / "reports"

FEATURE_ID = "10275997325121874585"
COMMON_SCANS = [5852, 5865, 5878, 5891, 5904, 5917]
EXPECTED_SPACING = 1.00335483507 / 2.0
SEARCH_WINDOW_DA = 0.15
TRACE_INDICES_TO_INSPECT = [5, 6]

PER_SCAN_OUT = OUT / "native_centroids_traces5_6_v2_per_scan.tsv"
CANDIDATES_OUT = OUT / "native_centroids_traces5_6_v2_candidates.tsv"
SUMMARY_OUT = OUT / "native_centroids_traces5_6_v2_summary.tsv"


def decode_binary_array(array_element, default_dtype):
    """Decode a mzML binaryDataArray, supporting zlib and 32/64-bit floats."""
    cv_params = {
        cv.attrib.get("accession"): cv.attrib.get("name", "")
        for cv in array_element.findall(
            "{http://psi.hupo.org/ms/mzml}cvParam"
        )
    }

    binary = array_element.find("{http://psi.hupo.org/ms/mzml}binary")
    if binary is None or not binary.text:
        return np.array([], dtype=default_dtype)

    raw = base64.b64decode(binary.text)
    if "MS:1000574" in cv_params:  # zlib compression
        raw = zlib_decompress(raw)

    if "MS:1000523" in cv_params:  # 64-bit float
        dtype = "<f8"
    elif "MS:1000521" in cv_params:  # 32-bit float
        dtype = "<f4"
    else:
        dtype = default_dtype

    return np.frombuffer(raw, dtype=dtype)


def zlib_decompress(data):
    import zlib
    return zlib.decompress(data)


def parse_scan_number(spectrum_element):
    """Extract scan number from common mzML native spectrum IDs."""
    spectrum_id = spectrum_element.attrib.get("id", "")
    match = re.search(r"(?:scan|spectrum)=(\d+)", spectrum_id)
    if match:
        return int(match.group(1))

    match = re.search(r"scan=(\d+)", spectrum_id)
    return int(match.group(1)) if match else None


def parse_spectrum(spectrum_element, ns):
    """Return m/z array, intensity array and scan metadata."""
    arrays = spectrum_element.findall(
        ".//{http://psi.hupo.org/ms/mzml}binaryDataArray"
    )
    mz = None
    intensity = None

    for array in arrays:
        cv_names = {
            cv.attrib.get("accession")
            for cv in array.findall(
                "{http://psi.hupo.org/ms/mzml}cvParam"
            )
        }
        decoded = decode_binary_array(array, "<f8")
        if "MS:1000514" in cv_names:  # m/z array
            mz = decoded
        elif "MS:1000515" in cv_names:  # intensity array
            intensity = decoded

    if mz is None or intensity is None:
        return None

    params = spectrum_element.findall(
        ".//{http://psi.hupo.org/ms/mzml}scan/"
        "{http://psi.hupo.org/ms/mzml}cvParam"
    )
    rt_sec = np.nan
    for param in params:
        if param.attrib.get("accession") == "MS:1000016":
            rt = float(param.attrib.get("value", "nan"))
            unit = param.attrib.get("unitAccession", "")
            rt_sec = rt * 60.0 if unit == "UO:0000031" else rt
            break

    ms_level = None
    for param in spectrum_element.findall(
        "{http://psi.hupo.org/ms/mzml}cvParam"
    ):
        if param.attrib.get("accession") == "MS:1000511":
            ms_level = int(float(param.attrib.get("value", "0")))
            break

    return {
        "spectrum_id": spectrum_element.attrib.get("id", ""),
        "scan_number": parse_scan_number(spectrum_element),
        "ms_level": ms_level,
        "rt_sec": rt_sec,
        "mz": np.asarray(mz, dtype=float),
        "intensity": np.asarray(intensity, dtype=float),
    }


def nearest_centroid(mz_array, intensity_array, target_mz, window):
    mask = np.abs(mz_array - target_mz) <= window
    indices = np.where(mask)[0]
    if len(indices) == 0:
        return None

    best_index = indices[np.argmin(np.abs(mz_array[indices] - target_mz))]
    return {
        "centroid_mz": float(mz_array[best_index]),
        "centroid_intensity": float(intensity_array[best_index]),
        "centroid_error_mDa": float(
            (mz_array[best_index] - target_mz) * 1000.0
        ),
        "candidate_count": int(len(indices)),
        "max_intensity_in_window": float(np.max(intensity_array[indices])),
        "sum_intensity_in_window": float(np.sum(intensity_array[indices])),
    }


def main():
    if not TRACE_FILE.exists():
        raise FileNotFoundError(TRACE_FILE)
    if not MZML_FILE.exists():
        raise FileNotFoundError(MZML_FILE)

    OUT.mkdir(parents=True, exist_ok=True)
    traces = pd.read_csv(TRACE_FILE, sep="\t", low_memory=False)
    required = {
        "feature_id", "trace_index", "scan_index", "observed_mz",
        "peak_intensity",
    }
    missing = required.difference(traces.columns)
    if missing:
        raise ValueError(f"Trace file missing columns: {sorted(missing)}")

    traces = traces.loc[
        traces["feature_id"].astype(str).eq(FEATURE_ID)
        & traces["scan_index"].isin(COMMON_SCANS)
    ].copy()

    for col in ("trace_index", "scan_index", "observed_mz", "peak_intensity"):
        traces[col] = pd.to_numeric(traces[col], errors="coerce")
    traces = traces.dropna(
        subset=["trace_index", "scan_index", "observed_mz", "peak_intensity"]
    )
    traces["trace_index"] = traces["trace_index"].astype(int)
    traces["scan_index"] = traces["scan_index"].astype(int)

    if traces.empty:
        raise ValueError("No trace observations found for the common scans.")

    # Parse mzML once and retain MS1 spectra with scan numbers matching the
    # trace export's scan_index. Report missing matches rather than guessing.
    ns = {"mzml": "http://psi.hupo.org/ms/mzml"}
    spectra_by_scan = {}

    # iterparse reduces memory use compared with loading the complete XML tree.
    with open(MZML_FILE, "rb") as handle:
        context = ET.iterparse(handle, events=("end",))
        for _, element in context:
            if element.tag.endswith("}spectrum"):
                parsed = parse_spectrum(element, ns)
                if (
                    parsed is not None
                    and parsed["ms_level"] == 1
                    and parsed["scan_number"] in COMMON_SCANS
                ):
                    spectra_by_scan[parsed["scan_number"]] = parsed
                element.clear()

    missing_scans = sorted(set(COMMON_SCANS) - set(spectra_by_scan))
    if missing_scans:
        print(
            "WARNING: no MS1 spectra matched these native scan numbers: "
            + ", ".join(map(str, missing_scans))
        )

    per_scan_rows = []
    candidate_rows = []

    for scan in COMMON_SCANS:
        group = traces.loc[traces["scan_index"].eq(scan)].copy()
        spectrum = spectra_by_scan.get(scan)

        if spectrum is None:
            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index": scan,
                "status": "no_matching_native_ms1_spectrum",
            })
            continue

        trace_map = {
            int(row["trace_index"]): row
            for _, row in group.iterrows()
        }

        # Fit an intercept from traces 1-4, holding isotope spacing fixed.
        if not all(i in trace_map for i in (1, 2, 3, 4)):
            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index": scan,
                "status": "missing_one_or_more_traces_1_to_4",
                "spectrum_id": spectrum["spectrum_id"],
                "rt_sec": spectrum["rt_sec"],
            })
            continue

        base_mz = np.array(
            [float(trace_map[i]["observed_mz"]) for i in (1, 2, 3, 4)]
        )
        base_index = np.arange(4, dtype=float)
        intercept = float(np.mean(base_mz - base_index * EXPECTED_SPACING))

        for trace_index in TRACE_INDICES_TO_INSPECT:
            row = trace_map.get(trace_index)
            assumed_index = trace_index - 1
            predicted_mz = intercept + assumed_index * EXPECTED_SPACING

            if row is None:
                observed_mz = np.nan
                trace_intensity = np.nan
                trace_error_mDa = np.nan
            else:
                observed_mz = float(row["observed_mz"])
                trace_intensity = float(row["peak_intensity"])
                trace_error_mDa = (observed_mz - predicted_mz) * 1000.0

            centroid = nearest_centroid(
                spectrum["mz"], spectrum["intensity"],
                predicted_mz, SEARCH_WINDOW_DA
            )

            if centroid is None:
                centroid_status = "no_centroid_in_window"
                centroid_mz = np.nan
                centroid_intensity = np.nan
                centroid_error_mDa = np.nan
                candidate_count = 0
                max_intensity = np.nan
                sum_intensity = np.nan
            else:
                centroid_status = "centroid_found"
                centroid_mz = centroid["centroid_mz"]
                centroid_intensity = centroid["centroid_intensity"]
                centroid_error_mDa = centroid["centroid_error_mDa"]
                candidate_count = centroid["candidate_count"]
                max_intensity = centroid["max_intensity_in_window"]
                sum_intensity = centroid["sum_intensity_in_window"]

                candidate_rows.append({
                    "feature_id": FEATURE_ID,
                    "scan_index": scan,
                    "trace_index": trace_index,
                    "spectrum_id": spectrum["spectrum_id"],
                    "rt_sec": spectrum["rt_sec"],
                    "predicted_mz_from_traces1_4": predicted_mz,
                    "centroid_mz": centroid_mz,
                    "centroid_intensity": centroid_intensity,
                    "centroid_error_mDa": centroid_error_mDa,
                    "candidate_count_in_window": candidate_count,
                    "max_intensity_in_window": max_intensity,
                    "sum_intensity_in_window": sum_intensity,
                })

            per_scan_rows.append({
                "feature_id": FEATURE_ID,
                "scan_index": scan,
                "trace_index": trace_index,
                "spectrum_id": spectrum["spectrum_id"],
                "rt_sec": spectrum["rt_sec"],
                "trace_observation_present": row is not None,
                "trace_observed_mz": observed_mz,
                "trace_peak_intensity": trace_intensity,
                "four_trace_intercept_mz": intercept,
                "expected_spacing_Da": EXPECTED_SPACING,
                "predicted_mz_from_traces1_4": predicted_mz,
                "trace_residual_vs_prediction_mDa": trace_error_mDa,
                "native_centroid_status": centroid_status,
                "nearest_native_centroid_mz": centroid_mz,
                "nearest_native_centroid_intensity": centroid_intensity,
                "nearest_centroid_error_vs_prediction_mDa": centroid_error_mDa,
                "centroid_candidates_within_window": candidate_count,
                "max_intensity_in_window": max_intensity,
                "sum_intensity_in_window": sum_intensity,
                "trace_to_centroid_delta_mDa": (
                    (observed_mz - centroid_mz) * 1000.0
                    if row is not None and np.isfinite(centroid_mz)
                    else np.nan
                ),
            })

    per_scan_df = pd.DataFrame(per_scan_rows)
    candidate_df = pd.DataFrame(candidate_rows)
    per_scan_df.to_csv(PER_SCAN_OUT, sep="\t", index=False)
    candidate_df.to_csv(CANDIDATES_OUT, sep="\t", index=False)

    summary_rows = []
    if not per_scan_df.empty and "trace_index" in per_scan_df.columns:
        valid = per_scan_df.loc[
            per_scan_df["trace_index"].isin(TRACE_INDICES_TO_INSPECT)
        ]
        for trace_index, group in valid.groupby("trace_index"):
            found = group.loc[
                group["native_centroid_status"].eq("centroid_found")
            ]
            trace_resid = pd.to_numeric(
                group["trace_residual_vs_prediction_mDa"], errors="coerce"
            ).dropna()
            centroid_resid = pd.to_numeric(
                found["nearest_centroid_error_vs_prediction_mDa"],
                errors="coerce"
            ).dropna()

            summary_rows.append({
                "trace_index": int(trace_index),
                "n_common_scans": len(group),
                "n_trace_observations": int(
                    group["trace_observation_present"].sum()
                ),
                "n_scans_with_native_centroid": len(found),
                "median_abs_trace_residual_mDa": (
                    float(np.median(np.abs(trace_resid)))
                    if len(trace_resid) else np.nan
                ),
                "max_abs_trace_residual_mDa": (
                    float(np.max(np.abs(trace_resid)))
                    if len(trace_resid) else np.nan
                ),
                "median_abs_native_centroid_error_mDa": (
                    float(np.median(np.abs(centroid_resid)))
                    if len(centroid_resid) else np.nan
                ),
                "max_abs_native_centroid_error_mDa": (
                    float(np.max(np.abs(centroid_resid)))
                    if len(centroid_resid) else np.nan
                ),
                "median_trace_to_centroid_delta_mDa": (
                    float(np.median(pd.to_numeric(
                        found["trace_to_centroid_delta_mDa"],
                        errors="coerce"
                    ).dropna()))
                    if "trace_to_centroid_delta_mDa" in found.columns
                    and found["trace_to_centroid_delta_mDa"].notna().any()
                    else np.nan
                ),
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(SUMMARY_OUT, sep="\t", index=False)

    print(f"Feature: {FEATURE_ID}")
    print(f"mzML: {MZML_FILE}")
    print(f"Search window: +/- {SEARCH_WINDOW_DA:.3f} Da")
    print(f"Common scans: {COMMON_SCANS}")
    print(f"Native MS1 spectra matched: {len(spectra_by_scan)}")
    print(f"Per-scan results: {PER_SCAN_OUT}")
    print(f"Centroid candidates: {CANDIDATES_OUT}")
    print(f"Summary: {SUMMARY_OUT}")
    print("\nPer-scan observations:")
    print(per_scan_df.to_string(index=False))
    print("\nSummary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
