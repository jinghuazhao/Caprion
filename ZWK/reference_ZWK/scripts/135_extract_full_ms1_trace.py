#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from pyteomics import mzml

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)

RT_MIN = 740.0
RT_MAX = 790.0
MZ_MIN = 1340.8
MZ_MAX = 1343.0

OUT_SCANS = Path(
    "reports/feature10275997325121874585_mzml_ms1_scan_inventory_v1.tsv"
)
OUT_PEAKS = Path(
    "reports/feature10275997325121874585_mzml_ms1_centroids_v1.tsv"
)
OUT_SUMMARY = Path(
    "reports/feature10275997325121874585_mzml_ms1_extraction_summary_v1.tsv"
)


def get_rt_seconds(spec):
    """Convert mzML scan start time from minutes to seconds."""
    scans = spec.get("scanList", {}).get("scan", [])
    if not scans:
        return np.nan

    value = scans[0].get("scan start time")
    if value is None:
        return np.nan

    # This mzML stores scan start times in minutes.
    return float(value) * 60.0


def main():
    if not MZML.exists():
        raise SystemExit(f"mzML not found: {MZML}")

    scan_rows = []
    peak_rows = []
    n_ms1_total = 0
    n_in_window = 0

    with mzml.MzML(str(MZML), use_index=True) as reader:
        for spec in reader:
            ms_level = spec.get("ms level")
            if ms_level != 1:
                continue

            n_ms1_total += 1
            rt = get_rt_seconds(spec)
            if not np.isfinite(rt) or not (RT_MIN <= rt <= RT_MAX):
                continue

            n_in_window += 1
            mz = np.asarray(spec.get("m/z array", []), dtype=float)
            intensity = np.asarray(
                spec.get("intensity array", []), dtype=float
            )

            valid = (
                np.isfinite(mz)
                & np.isfinite(intensity)
                & (mz >= MZ_MIN)
                & (mz <= MZ_MAX)
            )
            indices = np.flatnonzero(valid)

            scan_id = spec.get("id", "")
            scan_number = spec.get("index", None)
            if scan_number is not None:
                scan_number = int(scan_number) + 1

            scan_rows.append({
                "spectrum_index_zero_based": spec.get("index", np.nan),
                "spectrum_id": scan_id,
                "scan_number_from_index": scan_number,
                "rt_seconds": rt,
                "ms_level": ms_level,
                "n_peaks_in_mz_window": len(indices),
                "base_peak_intensity_in_window": (
                    float(np.max(intensity[indices])) if len(indices) else 0.0
                ),
            })

            for k in indices:
                peak_rows.append({
                    "spectrum_index_zero_based": spec.get("index", np.nan),
                    "spectrum_id": scan_id,
                    "scan_number_from_index": scan_number,
                    "rt_seconds": rt,
                    "centroid_mz": mz[k],
                    "centroid_intensity": intensity[k],
                })

    scans = pd.DataFrame(scan_rows)
    peaks = pd.DataFrame(peak_rows)

    OUT_SCANS.parent.mkdir(parents=True, exist_ok=True)
    scans.to_csv(OUT_SCANS, sep="\t", index=False, na_rep="NA")
    peaks.to_csv(OUT_PEAKS, sep="\t", index=False, na_rep="NA")

    summary = pd.DataFrame([
        {"metric": "ms1_scans_total", "value": n_ms1_total},
        {"metric": "ms1_scans_in_rt_window", "value": n_in_window},
        {"metric": "rt_min_seconds", "value": RT_MIN},
        {"metric": "rt_max_seconds", "value": RT_MAX},
        {"metric": "mz_min", "value": MZ_MIN},
        {"metric": "mz_max", "value": MZ_MAX},
        {"metric": "scans_with_window_peaks",
         "value": int((scans["n_peaks_in_mz_window"] > 0).sum())
         if len(scans) else 0},
        {"metric": "scans_without_window_peaks",
         "value": int((scans["n_peaks_in_mz_window"] == 0).sum())
         if len(scans) else 0},
        {"metric": "centroid_rows_in_window", "value": len(peaks)},
    ])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"Total MS1 scans in file: {n_ms1_total}")
    print(f"MS1 scans in RT window: {n_in_window}")
    print(f"Scans with peaks in m/z window: "
          f"{(scans['n_peaks_in_mz_window'] > 0).sum() if len(scans) else 0}")
    print(f"Scans without peaks in m/z window: "
          f"{(scans['n_peaks_in_mz_window'] == 0).sum() if len(scans) else 0}")
    print(f"Centroids extracted: {len(peaks)}")
    print("\nOutputs:")
    print(f"  {OUT_SCANS}")
    print(f"  {OUT_PEAKS}")
    print(f"  {OUT_SUMMARY}")

    if len(scans):
        print("\nFirst and last scans:")
        print(scans[[
            "spectrum_index_zero_based", "spectrum_id",
            "rt_seconds", "n_peaks_in_mz_window"
        ]].iloc[[0, -1]].to_string(index=False))


if __name__ == "__main__":
    main()
