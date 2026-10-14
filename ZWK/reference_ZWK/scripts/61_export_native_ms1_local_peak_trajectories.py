#!/usr/bin/env python3
"""Export all native MS1 centroid peaks in a local m/z window.

No nearest-target matching, isotope-group merging, or monoisotopic-anchor
assignment is performed. Peak ranks are local to each scan.
"""
from pathlib import Path
import csv
import math
import statistics
import sys
from collections import defaultdict

try:
    import pyopenms as oms
except ImportError:
    sys.exit("ERROR: pyopenms is unavailable in this Python environment.")

ROOT = Path.cwd()
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
REPORTS = ROOT / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)

SCAN_START = 5800
SCAN_END = 6020
MZ_MIN = 1340.90
MZ_MAX = 1343.70
EXPECTED_Z2_SPACING = 1.00335483507 / 2.0

PEAK_OUT = REPORTS / "native_ms1_all_local_centroids_scan5800_6020_v2.tsv"
GAP_OUT = REPORTS / "native_ms1_all_local_centroid_gaps_scan5800_6020_v2.tsv"
SUMMARY_OUT = REPORTS / "native_ms1_local_peak_trajectory_summary_v2.tsv"


def fmt(x, n=6):
    if x is None:
        return ""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return ""
    return f"{x:.{n}f}" if math.isfinite(x) else ""


def write_tsv(path, rows, fields=None):
    if not rows and fields is None:
        print(f"WARNING: no rows to write: {path}")
        return
    if rows:
        fields = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    if not MZML.is_file():
        sys.exit(f"ERROR: mzML file not found: {MZML}")

    print(f"Loading {MZML}")
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(MZML), exp)
    spectra = exp.getSpectra()
    print(f"Loaded {len(spectra):,} spectra")

    peak_rows = []
    gap_rows = []
    trajectory = defaultdict(list)
    n_ms1 = 0
    n_scans_with_peaks = 0

    for scan_idx in range(max(0, SCAN_START),
                          min(len(spectra), SCAN_END + 1)):
        spec = spectra[scan_idx]
        if int(spec.getMSLevel()) != 1:
            continue

        n_ms1 += 1
        rt = float(spec.getRT())
        mzs, intensities = spec.get_peaks()

        peaks = sorted(
            [
                (float(mz), float(intensity))
                for mz, intensity in zip(mzs, intensities)
                if MZ_MIN <= float(mz) <= MZ_MAX
            ],
            key=lambda p: p[0],
        )

        if peaks:
            n_scans_with_peaks += 1

        for rank, (mz, intensity) in enumerate(peaks, start=1):
            trajectory[rank].append({
                "scan_index": scan_idx,
                "rt_sec": rt,
                "mz": mz,
                "intensity": intensity,
            })

            peak_rows.append({
                "scan_index": scan_idx,
                "rt_sec": fmt(rt, 5),
                "local_peak_rank": rank,
                "mz": fmt(mz, 6),
                "intensity": fmt(intensity, 3),
                "window_mz_min": fmt(MZ_MIN, 3),
                "window_mz_max": fmt(MZ_MAX, 3),
                "rank_note": "scan-local rank only; not an isotope index",
            })

        for rank in range(1, len(peaks)):
            left_mz, left_intensity = peaks[rank - 1]
            right_mz, right_intensity = peaks[rank]
            gap_da = right_mz - left_mz
            nearest_steps = int(round(gap_da / EXPECTED_Z2_SPACING))
            residual_mda = (
                gap_da - nearest_steps * EXPECTED_Z2_SPACING
            ) * 1000.0

            gap_rows.append({
                "scan_index": scan_idx,
                "rt_sec": fmt(rt, 5),
                "left_local_peak_rank": rank,
                "left_mz": fmt(left_mz, 6),
                "left_intensity": fmt(left_intensity, 3),
                "right_local_peak_rank": rank + 1,
                "right_mz": fmt(right_mz, 6),
                "right_intensity": fmt(right_intensity, 3),
                "gap_da": fmt(gap_da, 9),
                "gap_mDa": fmt(gap_da * 1000.0, 3),
                "nearest_integer_z2_steps": nearest_steps,
                "residual_from_nearest_integer_z2_steps_mDa": fmt(
                    residual_mda, 3
                ),
                "gap_is_near_one_z2_step_within_2mDa": int(
                    abs(gap_da - EXPECTED_Z2_SPACING) * 1000.0 <= 2.0
                ),
            })

    summary_rows = []
    print(
        f"\nSelected {n_ms1} MS1 scans; "
        f"{n_scans_with_peaks} contained peaks in the local window."
    )

    print("\nScan-local rank summaries (rank can change identity between scans)")
    print("rank  n_scans  RT_min      RT_max      median_mz    "
          "mz_min       mz_max       intensity_sum")

    for rank in sorted(trajectory):
        obs = trajectory[rank]
        mz_values = [x["mz"] for x in obs]
        rts = [x["rt_sec"] for x in obs]
        intensity_sum = sum(x["intensity"] for x in obs)

        row = {
            "local_peak_rank": rank,
            "n_scans": len(obs),
            "rt_min_sec": fmt(min(rts), 5),
            "rt_max_sec": fmt(max(rts), 5),
            "median_mz": fmt(statistics.median(mz_values), 6),
            "mz_min": fmt(min(mz_values), 6),
            "mz_max": fmt(max(mz_values), 6),
            "intensity_sum": fmt(intensity_sum, 3),
            "interpretation": (
                "scan-local rank; do not assume rank tracks the same ion "
                "across scans"
            ),
        }
        summary_rows.append(row)
        print(
            f"{rank:>4}  {len(obs):>7}  {min(rts):>10.5f}  "
            f"{max(rts):>10.5f}  {statistics.median(mz_values):>11.6f}  "
            f"{min(mz_values):>11.6f}  {max(mz_values):>11.6f}  "
            f"{intensity_sum:>13.3f}"
        )

    write_tsv(PEAK_OUT, peak_rows, fields=[
        "scan_index", "rt_sec", "local_peak_rank", "mz", "intensity",
        "window_mz_min", "window_mz_max", "rank_note",
    ])
    write_tsv(GAP_OUT, gap_rows, fields=[
        "scan_index", "rt_sec", "left_local_peak_rank", "left_mz",
        "left_intensity", "right_local_peak_rank", "right_mz",
        "right_intensity", "gap_da", "gap_mDa",
        "nearest_integer_z2_steps",
        "residual_from_nearest_integer_z2_steps_mDa",
        "gap_is_near_one_z2_step_within_2mDa",
    ])
    write_tsv(SUMMARY_OUT, summary_rows, fields=[
        "local_peak_rank", "n_scans", "rt_min_sec", "rt_max_sec",
        "median_mz", "mz_min", "mz_max", "intensity_sum",
        "interpretation",
    ])

    print("\nFiles written:")
    print(f"  {PEAK_OUT}")
    print(f"  {GAP_OUT}")
    print(f"  {SUMMARY_OUT}")
    print(
        "\nCaution: scan-local ranks are not persistent feature identities. "
        "Use the per-scan m/z values to inspect peak trajectories; do not "
        "infer an isotope envelope solely from adjacent peak gaps."
    )


if __name__ == "__main__":
    main()
