#!/usr/bin/env python3
"""Inspect all native MS1 centroid peaks near a provisional charge-2 ladder.

Diagnostic only: does not assign monoisotopic anchors, merge isotope groups,
or consult the original XLSX isotope-group assignments.
"""
from pathlib import Path
import csv
import math
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "ERROR: pyopenms is unavailable in this Python environment. "
        "Load/use the same environment used for Script 57."
    )

ROOT = Path.cwd()
MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/"
    "mzML/szwk900404i19201xms1.mzML"
)
REPORTS = ROOT / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)

# Native MS1 spectrum indices and RTs previously checked by Script 57.
TARGETS = [
    {"scan_index": 5954, "rt_sec_expected": 769.41804},
    {"scan_index": 5981, "rt_sec_expected": 772.57938},
    {"scan_index": 5988, "rt_sec_expected": 773.44866},
]

MZ_MIN = 1340.9
MZ_MAX = 1343.7
EXPECTED_SPACING_Z2 = 1.00335483507 / 2.0

PEAK_OUT = REPORTS / "native_ms1_local_centroid_peaks_scan5954_5981_5988.tsv"
GAP_OUT = REPORTS / "native_ms1_local_centroid_adjacent_gaps_scan5954_5981_5988.tsv"
SUMMARY_OUT = REPORTS / "native_ms1_local_centroid_peak_summary_scan5954_5981_5988.tsv"


def fmt(value, digits=6):
    if value is None or not math.isfinite(float(value)):
        return ""
    return f"{float(value):.{digits}f}"


def main():
    if not MZML.is_file():
        sys.exit(f"ERROR: mzML file not found: {MZML}")

    print(f"Input mzML: {MZML}")
    print(f"Local centroid m/z window: {MZ_MIN:.4f}–{MZ_MAX:.4f}")
    print("Loading mzML; this may take a little while...")

    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(MZML), exp)
    spectra = exp.getSpectra()
    print(f"Loaded spectra: {len(spectra):,}")

    peak_rows = []
    gap_rows = []
    summary_rows = []

    for target in TARGETS:
        idx = target["scan_index"]
        if idx < 0 or idx >= len(spectra):
            print(f"WARNING: requested spectrum index {idx} is outside the file.")
            continue

        spec = spectra[idx]
        rt = float(spec.getRT())
        ms_level = int(spec.getMSLevel())
        peaks_mz, peaks_int = spec.get_peaks()

        local = []
        for mz, intensity in zip(peaks_mz, peaks_int):
            mz = float(mz)
            intensity = float(intensity)
            if MZ_MIN <= mz <= MZ_MAX:
                local.append((mz, intensity))

        local.sort(key=lambda item: item[0])

        rt_delta = rt - target["rt_sec_expected"]
        summary_rows.append({
            "scan_index": idx,
            "ms_level": ms_level,
            "rt_sec_native": fmt(rt, 5),
            "rt_sec_expected": fmt(target["rt_sec_expected"], 5),
            "rt_delta_sec": fmt(rt_delta, 6),
            "local_peak_count": len(local),
            "window_mz_min": fmt(MZ_MIN, 4),
            "window_mz_max": fmt(MZ_MAX, 4),
            "expected_z2_isotope_spacing_da": fmt(EXPECTED_SPACING_Z2, 9),
        })

        print(
            f"\nScan index {idx}: MS level {ms_level}; "
            f"RT={rt:.5f} s; RT delta={rt_delta:+.6f} s; "
            f"{len(local)} local centroid peaks"
        )
        print("  peak  m/z             intensity       gap from previous (mDa)")

        previous_mz = None
        for rank, (mz, intensity) in enumerate(local, start=1):
            gap_da = None if previous_mz is None else mz - previous_mz
            gap_mda = None if gap_da is None else gap_da * 1000.0
            nearest_z2_steps = (
                None if gap_da is None
                else int(round(gap_da / EXPECTED_SPACING_Z2))
            )
            nearest_z2_resid_mda = (
                None if gap_da is None
                else (
                    gap_da - nearest_z2_steps * EXPECTED_SPACING_Z2
                ) * 1000.0
            )

            peak_rows.append({
                "scan_index": idx,
                "ms_level": ms_level,
                "rt_sec_native": fmt(rt, 5),
                "rt_sec_expected": fmt(target["rt_sec_expected"], 5),
                "rt_delta_sec": fmt(rt_delta, 6),
                "peak_rank_by_mz": rank,
                "mz": fmt(mz, 6),
                "intensity": fmt(intensity, 3),
                "gap_from_previous_da": "" if gap_da is None else fmt(gap_da, 9),
                "gap_from_previous_mDa": "" if gap_mda is None else fmt(gap_mda, 3),
                "nearest_z2_spacing_steps": "" if nearest_z2_steps is None else nearest_z2_steps,
                "residual_from_nearest_z2_steps_mDa": (
                    "" if nearest_z2_resid_mda is None
                    else fmt(nearest_z2_resid_mda, 3)
                ),
            })

            if gap_da is None:
                print(f"  {rank:>4}  {mz:>13.6f}  {intensity:>14.3f}       —")
            else:
                print(
                    f"  {rank:>4}  {mz:>13.6f}  {intensity:>14.3f}  "
                    f"{gap_mda:>12.3f}"
                )

            if previous_mz is not None:
                gap_rows.append({
                    "scan_index": idx,
                    "ms_level": ms_level,
                    "rt_sec_native": fmt(rt, 5),
                    "rt_delta_sec": fmt(rt_delta, 6),
                    "left_peak_rank": rank - 1,
                    "left_mz": fmt(previous_mz, 6),
                    "right_peak_rank": rank,
                    "right_mz": fmt(mz, 6),
                    "gap_da": fmt(gap_da, 9),
                    "gap_mDa": fmt(gap_mda, 3),
                    "nearest_z2_spacing_steps": nearest_z2_steps,
                    "expected_gap_for_nearest_steps_da": fmt(
                        nearest_z2_steps * EXPECTED_SPACING_Z2, 9
                    ),
                    "residual_from_nearest_z2_steps_mDa": fmt(
                        nearest_z2_resid_mda, 3
                    ),
                    "left_intensity": fmt(local[rank - 2][1], 3),
                    "right_intensity": fmt(intensity, 3),
                })

            previous_mz = mz

        if not local:
            print("  No centroid peaks in this m/z window.")

    def write_tsv(path, rows):
        if not rows:
            print(f"WARNING: no rows to write for {path.name}")
            return
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0].keys()), delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

    write_tsv(PEAK_OUT, peak_rows)
    write_tsv(GAP_OUT, gap_rows)
    write_tsv(SUMMARY_OUT, summary_rows)

    print("\nFiles written:")
    print(f"  {PEAK_OUT}")
    print(f"  {GAP_OUT}")
    print(f"  {SUMMARY_OUT}")

    print("\nSummary:")
    for row in summary_rows:
        print(
            f"  scan {row['scan_index']}: MS{row['ms_level']}, "
            f"RT={row['rt_sec_native']} s, "
            f"{row['local_peak_count']} peaks in window"
        )

    print(
        "\nInterpretation caution: adjacent centroid peaks are not automatically "
        "isotopes. A gap near 0.501677 Da is compatible with one isotope step "
        "at charge 2, but does not by itself prove a shared isotope envelope "
        "or identify the monoisotopic peak."
    )


if __name__ == "__main__":
    main()
