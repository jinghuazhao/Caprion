#!/usr/bin/env python3
"""Track native centroid peaks around a provisional charge-2 isotope ladder.

Uses native mzML MS1 spectra only. Does not use the original XLSX IG
assignments, merge features, or assign monoisotopic anchors.
"""
from pathlib import Path
import csv
import math
import statistics
import sys

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

# Provisional trace m/z values from the first pilot feature.
# These are tracking targets, not asserted isotope indices.
TARGETS = [
    {"trace": 1, "mz": 1341.067749, "reference_scan": 5981},
    {"trace": 2, "mz": 1341.570068, "reference_scan": 5981},
    {"trace": 3, "mz": 1342.067993, "reference_scan": 5981},
    {"trace": 4, "mz": 1342.570496, "reference_scan": 5981},
    {"trace": 5, "mz": 1343.074219, "reference_scan": 5954},
    {"trace": 6, "mz": 1343.570129, "reference_scan": 5954},
]

# Include spectra around the observed elution interval.
SCAN_START = 5800
SCAN_END = 6020

# Peak matching settings: nearest native centroid peak within this window.
# Keep this relatively narrow so neighbouring centroid peaks are not silently
# conflated. Results outside the window remain unmatched.
MATCH_TOL_DA = 0.015
EXPECTED_SPACING_Z2 = 1.00335483507 / 2.0

TRACK_OUT = REPORTS / "native_ms1_centroid_tracks_charge2_v2.tsv"
GAP_OUT = REPORTS / "native_ms1_centroid_track_adjacent_gaps_charge2_v2.tsv"
SUMMARY_OUT = REPORTS / "native_ms1_centroid_track_summary_charge2_v2.tsv"


def write_tsv(path, rows):
    if not rows:
        print(f"WARNING: no rows to write: {path}")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0].keys()), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def fmt(x, n=6):
    if x is None:
        return ""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(x):
        return ""
    return f"{x:.{n}f}"


def median_or_none(values):
    return statistics.median(values) if values else None


def main():
    if not MZML.is_file():
        sys.exit(f"ERROR: mzML not found: {MZML}")

    print(f"Loading {MZML}")
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(MZML), exp)
    spectra = exp.getSpectra()
    print(f"Loaded {len(spectra):,} spectra")

    scan_indices = [
        i for i in range(max(0, SCAN_START), min(len(spectra), SCAN_END + 1))
        if int(spectra[i].getMSLevel()) == 1
    ]
    print(
        f"Inspecting MS1 spectra with indices {SCAN_START}–{SCAN_END}; "
        f"{len(scan_indices)} MS1 spectra selected"
    )

    # Store each target's nearest peak in each MS1 spectrum. Matching is
    # independent per target; a centroid can therefore match more than one
    # target if the targets are closer than twice the tolerance.
    track_rows = []
    target_observations = {t["trace"]: [] for t in TARGETS}

    for scan_idx in scan_indices:
        spec = spectra[scan_idx]
        rt = float(spec.getRT())
        mzs, intensities = spec.get_peaks()
        local_peaks = sorted(
            (
                (float(mz), float(intensity))
                for mz, intensity in zip(mzs, intensities)
                if 1340.90 <= float(mz) <= 1343.70
            ),
            key=lambda x: x[0],
        )

        for target in TARGETS:
            target_mz = target["mz"]
            candidates = [
                (mz, intensity)
                for mz, intensity in local_peaks
                if abs(mz - target_mz) <= MATCH_TOL_DA
            ]
            if candidates:
                matched_mz, matched_intensity = min(
                    candidates, key=lambda p: abs(p[0] - target_mz)
                )
                error_mDa = (matched_mz - target_mz) * 1000.0
                found = True
            else:
                matched_mz = None
                matched_intensity = None
                error_mDa = None
                found = False

            row = {
                "scan_index": scan_idx,
                "rt_sec": fmt(rt, 5),
                "trace_target": target["trace"],
                "target_mz": fmt(target_mz, 6),
                "matched": int(found),
                "matched_mz": fmt(matched_mz, 6),
                "intensity": fmt(matched_intensity, 3),
                "mz_error_mDa": fmt(error_mDa, 3),
                "local_centroid_peak_count": len(local_peaks),
                "match_tolerance_da": fmt(MATCH_TOL_DA, 3),
            }
            track_rows.append(row)

            if found:
                target_observations[target["trace"]].append({
                    "scan_index": scan_idx,
                    "rt_sec": rt,
                    "mz": matched_mz,
                    "intensity": matched_intensity,
                    "error_mDa": error_mDa,
                })

    # Adjacent trace targets: compare observed centroids only where both
    # targets have a match in the same MS1 spectrum.
    gap_rows = []
    for left, right in zip(TARGETS[:-1], TARGETS[1:]):
        left_by_scan = {
            x["scan_index"]: x for x in target_observations[left["trace"]]
        }
        right_by_scan = {
            x["scan_index"]: x for x in target_observations[right["trace"]]
        }
        common_scans = sorted(set(left_by_scan) & set(right_by_scan))

        for scan_idx in common_scans:
            a = left_by_scan[scan_idx]
            b = right_by_scan[scan_idx]
            gap_da = b["mz"] - a["mz"]
            residual_mDa = (gap_da - EXPECTED_SPACING_Z2) * 1000.0
            gap_rows.append({
                "scan_index": scan_idx,
                "rt_sec": fmt(a["rt_sec"], 5),
                "left_trace_target": left["trace"],
                "right_trace_target": right["trace"],
                "left_observed_mz": fmt(a["mz"], 6),
                "right_observed_mz": fmt(b["mz"], 6),
                "observed_gap_da": fmt(gap_da, 9),
                "observed_gap_mDa": fmt(gap_da * 1000.0, 3),
                "expected_z2_gap_da": fmt(EXPECTED_SPACING_Z2, 9),
                "spacing_residual_mDa": fmt(residual_mDa, 3),
                "left_intensity": fmt(a["intensity"], 3),
                "right_intensity": fmt(b["intensity"], 3),
            })

    summary_rows = []
    print("\nPer-target trajectory summary")
    print("trace  target_mz    matched_scans  RT_min      RT_max      "
          "median_mz    median_abs_error_mDa  intensity_sum")

    for target in TARGETS:
        obs = target_observations[target["trace"]]
        mz_values = [x["mz"] for x in obs]
        errors = [abs(x["error_mDa"]) for x in obs]
        intensities = [x["intensity"] for x in obs]
        rts = [x["rt_sec"] for x in obs]

        summary = {
            "trace_target": target["trace"],
            "target_mz": fmt(target["mz"], 6),
            "scan_start": SCAN_START,
            "scan_end": SCAN_END,
            "match_tolerance_da": fmt(MATCH_TOL_DA, 3),
            "matched_scans": len(obs),
            "rt_min_sec": fmt(min(rts) if rts else None, 5),
            "rt_max_sec": fmt(max(rts) if rts else None, 5),
            "median_observed_mz": fmt(median_or_none(mz_values), 6),
            "mz_min": fmt(min(mz_values) if mz_values else None, 6),
            "mz_max": fmt(max(mz_values) if mz_values else None, 6),
            "median_abs_mz_error_mDa": fmt(median_or_none(errors), 3),
            "max_abs_mz_error_mDa": fmt(max(errors) if errors else None, 3),
            "intensity_sum": fmt(sum(intensities), 3),
            "note": "tracking target only; not an isotope-index or mono-anchor assignment",
        }
        summary_rows.append(summary)
        print(
            f"{target['trace']:>5}  {target['mz']:.6f}  "
            f"{len(obs):>13}  "
            f"{(min(rts) if rts else float('nan')):>10.5f}  "
            f"{(max(rts) if rts else float('nan')):>10.5f}  "
            f"{(median_or_none(mz_values) if mz_values else float('nan')):>11.6f}  "
            f"{(median_or_none(errors) if errors else float('nan')):>20.3f}  "
            f"{sum(intensities):>13.3f}"
        )

    write_tsv(TRACK_OUT, track_rows)
    write_tsv(GAP_OUT, gap_rows)
    write_tsv(SUMMARY_OUT, summary_rows)

    print("\nAdjacent target-spacing summary")
    for left, right in zip(TARGETS[:-1], TARGETS[1:]):
        rows = [
            r for r in gap_rows
            if r["left_trace_target"] == left["trace"]
            and r["right_trace_target"] == right["trace"]
        ]
        residuals = [abs(float(r["spacing_residual_mDa"])) for r in rows]
        within_2 = (
            sum(x <= 2.0 for x in residuals) / len(residuals)
            if residuals else None
        )
        print(
            f"trace {left['trace']}–{right['trace']}: "
            f"common scans={len(rows)}, "
            f"median |spacing residual|="
            f"{fmt(median_or_none(residuals), 3) if residuals else 'NA'} mDa, "
            f"fraction <=2 mDa="
            f"{fmt(within_2, 3) if within_2 is not None else 'NA'}"
        )

    print("\nFiles written:")
    print(f"  {TRACK_OUT}")
    print(f"  {GAP_OUT}")
    print(f"  {SUMMARY_OUT}")
    print(
        "\nCaution: independent nearest-peak matching can reuse a centroid for "
        "multiple targets and can follow unrelated peaks if m/z trajectories "
        "cross or drift. Inspect the output before drawing isotope-envelope "
        "conclusions."
    )


if __name__ == "__main__":
    main()
