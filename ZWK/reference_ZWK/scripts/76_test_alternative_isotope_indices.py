#!/usr/bin/env python3
"""Test alternative isotope indices for traces 5 and 6 of a candidate MS1 envelope."""

from pathlib import Path
import csv
import math
import statistics
import xml.etree.ElementTree as ET

ROOT = Path.cwd()
TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
MZML_FILE = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
REPORTS = ROOT / "reports"
REPORTS.mkdir(exist_ok=True)

FEATURE_ID = "10275997325121874585"
SCAN_INDICES = [5852, 5865, 5878, 5891, 5904, 5917]
CORE_TRACES = [1, 2, 3, 4]
TEST_TRACES = [5, 6]
SPACING = 1.00335483507 / 2.0
WINDOW_DA = 0.15
TOLERANCE_MDA = 3.0
# Include plausible indices below and above the original sequential assignment.
ISOTOPE_INDICES = list(range(-2, 10))


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


def children(element, name):
    return [x for x in element.iter() if local_name(x.tag) == name]


def read_trace_table(path):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    selected = []
    for row in rows:
        if row.get("feature_id") != FEATURE_ID:
            continue
        if int(row["scan_index"]) not in SCAN_INDICES:
            continue
        selected.append({
            "feature_id": row["feature_id"],
            "scan_index_zero_based": int(row["scan_index"]),
            "trace_index": int(row["trace_index"]),
            "rt_sec": float(row["rt_sec"]),
            "observed_mz": float(row["observed_mz"]),
            "peak_intensity": float(row["peak_intensity"]),
        })
    return selected


def decode_array(binary_array):
    # mzML binary arrays are base64 encoded and may be compressed.
    import base64
    import struct
    import zlib

    cv_params = {
        p.attrib.get("accession"): p.attrib
        for p in children(binary_array, "cvParam")
    }
    encoded = "".join(
        (x.text or "").strip()
        for x in children(binary_array, "binary")
    )
    if not encoded:
        return []

    raw = base64.b64decode(encoded)
    if "MS:1000574" in cv_params:  # zlib compression
        raw = zlib.decompress(raw)

    if "MS:1000523" in cv_params:  # 64-bit float
        width, fmt = 8, "d"
    elif "MS:1000521" in cv_params:  # 32-bit float
        width, fmt = 4, "f"
    else:
        raise ValueError("Unrecognised mzML binary-array numeric type")

    if len(raw) % width:
        raise ValueError("Binary array length is not a multiple of its data width")

    n = len(raw) // width
    return list(struct.unpack("<" + fmt * n, raw))


def parse_spectrum(spectrum):
    cv = {p.attrib.get("accession"): p.attrib for p in children(spectrum, "cvParam")}
    ms_level = int(cv.get("MS:1000511", {}).get("value", "0"))
    rt_value = cv.get("MS:1000016", {}).get("value")
    rt_unit = cv.get("MS:1000016", {}).get("unitAccession", "")
    rt = float(rt_value) if rt_value is not None else math.nan
    # mzML retention time is often encoded in minutes; convert to seconds.
    if rt_unit == "UO:0000031":
        rt *= 60.0

    mzs, intensities = [], []
    for array in children(spectrum, "binaryDataArray"):
        params = {p.attrib.get("accession"): p.attrib for p in children(array, "cvParam")}
        values = decode_array(array)
        if "MS:1000514" in params:
            mzs = values
        elif "MS:1000515" in params:
            intensities = values

    if len(mzs) != len(intensities):
        raise ValueError(f"m/z and intensity array lengths differ: {spectrum.attrib.get('id')}")

    return ms_level, rt, list(zip(mzs, intensities))


def read_selected_native_spectra(path, wanted_indices):
    found = {}
    max_index = max(wanted_indices)
    ordinal = -1

    # Clear each processed spectrum to limit memory usage.
    for event, element in ET.iterparse(path, events=("end",)):
        if local_name(element.tag) != "spectrum":
            continue
        ordinal += 1

        if ordinal in wanted_indices:
            ms_level, rt, peaks = parse_spectrum(element)
            if ms_level != 1:
                raise ValueError(
                    f"Expected MS1 at zero-based ordinal {ordinal}; got MS{ms_level}"
                )
            found[ordinal] = {
                "native_scan_number": ordinal + 1,
                "spectrum_id": element.attrib.get("id", ""),
                "rt_sec": rt,
                "peaks": peaks,
            }

        element.clear()
        if ordinal >= max_index:
            break

    missing = sorted(set(wanted_indices) - set(found))
    if missing:
        raise RuntimeError(f"Could not find requested spectrum ordinals: {missing}")
    return found


def nearest_peak(peaks, target_mz, window_da):
    candidates = [
        (mz, intensity)
        for mz, intensity in peaks
        if abs(mz - target_mz) <= window_da
    ]
    if not candidates:
        return None, None, 0
    mz, intensity = min(candidates, key=lambda x: abs(x[0] - target_mz))
    return mz, intensity, len(candidates)


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


def main():
    for p in [TRACE_FILE, MZML_FILE]:
        if not p.exists():
            raise FileNotFoundError(p)

    traces = read_trace_table(TRACE_FILE)
    if not traces:
        raise RuntimeError("No trace observations found for the requested feature/scans")

    by_key = {
        (r["scan_index_zero_based"], r["trace_index"]): r for r in traces
    }
    native = read_selected_native_spectra(MZML_FILE, SCAN_INDICES)

    fit_rows = []
    per_scan = []
    for scan in SCAN_INDICES:
        core = [
            by_key[(scan, idx)]
            for idx in CORE_TRACES
            if (scan, idx) in by_key
        ]
        if len(core) != len(CORE_TRACES):
            continue

        # Estimate m0 independently for each scan, using only traces 1-4.
        intercepts = [
            r["observed_mz"] - (r["trace_index"] - 1) * SPACING
            for r in core
        ]
        m0 = statistics.mean(intercepts)
        spectrum = native[scan]

        for trace_idx in TEST_TRACES:
            r = by_key.get((scan, trace_idx))
            if r is None:
                continue

            for isotope_index in ISOTOPE_INDICES:
                predicted = m0 + isotope_index * SPACING
                centroid_mz, centroid_intensity, n_candidates = nearest_peak(
                    spectrum["peaks"], predicted, WINDOW_DA
                )
                trace_residual = (r["observed_mz"] - predicted) * 1000.0
                centroid_residual = (
                    (centroid_mz - predicted) * 1000.0
                    if centroid_mz is not None else math.nan
                )
                fit_rows.append({
                    "feature_id": FEATURE_ID,
                    "scan_index_zero_based": scan,
                    "native_scan_number": spectrum["native_scan_number"],
                    "rt_sec": spectrum["rt_sec"],
                    "trace_index": trace_idx,
                    "tested_isotope_index": isotope_index,
                    "core_intercept_mz_traces1_4": f"{m0:.9f}",
                    "expected_spacing_Da": f"{SPACING:.12f}",
                    "trace_observed_mz": f"{r['observed_mz']:.9f}",
                    "trace_intensity": f"{r['peak_intensity']:.3f}",
                    "predicted_mz": f"{predicted:.9f}",
                    "trace_residual_mDa": f"{trace_residual:.6f}",
                    "nearest_centroid_mz": (
                        f"{centroid_mz:.9f}" if centroid_mz is not None else ""
                    ),
                    "nearest_centroid_intensity": (
                        f"{centroid_intensity:.3f}"
                        if centroid_intensity is not None else ""
                    ),
                    "centroid_error_mDa": (
                        f"{centroid_residual:.6f}"
                        if math.isfinite(centroid_residual) else ""
                    ),
                    "centroids_within_0_15_Da": n_candidates,
                    "trace_matches_nearest_centroid": (
                        abs(r["observed_mz"] - centroid_mz) < 1e-8
                        if centroid_mz is not None else False
                    ),
                })

    # Summarise isotope-index fits from trace m/z values. A fit is only
    # considered eligible when the trace has observations in at least 4 scans.
    for trace_idx in TEST_TRACES:
        for isotope_index in ISOTOPE_INDICES:
            rows = [
                r for r in fit_rows
                if r["trace_index"] == trace_idx
                and r["tested_isotope_index"] == isotope_index
            ]
            residuals = [abs(float(r["trace_residual_mDa"])) for r in rows]
            if not residuals:
                continue
            passing = sum(x <= TOLERANCE_MDA for x in residuals)
            summary = {
                "trace_index": trace_idx,
                "tested_isotope_index": isotope_index,
                "n_scans": len(rows),
                "median_abs_trace_residual_mDa": f"{statistics.median(residuals):.6f}",
                "max_abs_trace_residual_mDa": f"{max(residuals):.6f}",
                "trace_observations_within_3mDa": passing,
                "fraction_within_3mDa": f"{passing / len(rows):.6f}",
                "all_observations_within_3mDa": passing == len(rows),
            }
            per_scan.append(summary)

    per_scan.sort(key=lambda r: (
        r["trace_index"],
        float(r["median_abs_trace_residual_mDa"]),
        int(r["tested_isotope_index"]),
    ))

    best_by_trace = {}
    for trace_idx in TEST_TRACES:
        candidates = [r for r in per_scan if r["trace_index"] == trace_idx]
        if candidates:
            best_by_trace[trace_idx] = candidates[0]

    fit_path = REPORTS / "alternative_isotope_indices_v1_all_fits.tsv"
    summary_path = REPORTS / "alternative_isotope_indices_v1_summary.tsv"
    best_path = REPORTS / "alternative_isotope_indices_v1_best.tsv"

    write_tsv(fit_path, fit_rows)
    write_tsv(summary_path, per_scan)
    write_tsv(best_path, list(best_by_trace.values()))

    print(f"Feature: {FEATURE_ID}")
    print(f"Expected charge-2 spacing: {SPACING:.12f} Da")
    print(f"Selected MS1 spectra: {len(native)}")
    print(f"Trace observations tested: {sum(1 for r in traces if r['trace_index'] in TEST_TRACES)}")
    print(f"Tested isotope indices: {ISOTOPE_INDICES}")
    print(f"Fit report: {fit_path}")
    print(f"Summary report: {summary_path}")
    print(f"Best-index report: {best_path}")
    print("\nBest index by trace (ranked by median absolute residual):")
    for trace_idx in TEST_TRACES:
        candidates = [r for r in per_scan if r["trace_index"] == trace_idx]
        if not candidates:
            print(f"Trace {trace_idx}: no eligible observations")
            continue
        best = candidates[0]
        print(
            f"Trace {trace_idx}: isotope index {best['tested_isotope_index']}; "
            f"median abs residual {best['median_abs_trace_residual_mDa']} mDa; "
            f"max abs residual {best['max_abs_trace_residual_mDa']} mDa; "
            f"{best['trace_observations_within_3mDa']}/{best['n_scans']} "
            "observations within 3 mDa"
        )
    print(
        "\nInterpretation caution: this is a diagnostic index search, not an "
        "automatic isotope-group assignment. Candidate indices must be judged "
        "against spacing, scan consistency, and competing-envelope possibilities."
    )


if __name__ == "__main__":
    main()
