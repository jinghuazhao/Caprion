#!/usr/bin/env python3
"""Conservative m/z-continuity tracking of native MS1 centroid peaks.

This is a local diagnostic, not a validated feature detector. It does not use
the original XLSX assignments or make isotope-group/monoisotopic assignments.
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

# Conservative local tracking settings. A centroid may move by at most this
# amount per missing-scan-adjusted step. These are diagnostic parameters.
MAX_GAP_SCANS = 2
MAX_STEP_MDA = 12.0
MIN_TRACK_OBSERVATIONS = 3

OBS_OUT = REPORTS / "native_ms1_centroid_continuity_observations_v2.tsv"
TRACK_OUT = REPORTS / "native_ms1_centroid_continuity_track_summary_v2.tsv"
AMBIG_OUT = REPORTS / "native_ms1_centroid_continuity_ambiguities_v2.tsv"


def fmt(x, n=6):
    if x is None:
        return ""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return ""
    return f"{x:.{n}f}" if math.isfinite(x) else ""


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    if not MZML.is_file():
        sys.exit(f"ERROR: mzML not found: {MZML}")

    print(f"Loading {MZML}")
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(MZML), exp)
    spectra = exp.getSpectra()
    print(f"Loaded {len(spectra):,} spectra")

    # Retain each MS1 scan as a time-ordered set of centroid peaks.
    scans = []
    for scan_idx in range(max(0, SCAN_START),
                          min(len(spectra), SCAN_END + 1)):
        spec = spectra[scan_idx]
        if int(spec.getMSLevel()) != 1:
            continue

        rt = float(spec.getRT())
        mzs, intensities = spec.get_peaks()
        peaks = sorted(
            [
                (float(mz), float(intensity))
                for mz, intensity in zip(mzs, intensities)
                if MZ_MIN <= float(mz) <= MZ_MAX
            ],
            key=lambda x: x[0],
        )
        scans.append({
            "scan_index": scan_idx,
            "rt_sec": rt,
            "peaks": peaks,
        })

    print(f"MS1 scans selected: {len(scans)}")
    if not scans:
        sys.exit("ERROR: no MS1 scans in the selected index interval.")

    # Greedy continuity tracking:
    # for each new scan, form candidate links from recent active tracks.
    # A link is accepted only when it is mutually unique (the track's closest
    # peak and the peak's closest track). Ties/competing links are recorded
    # rather than arbitrarily forced.
    tracks = []
    ambiguities = []
    next_track_id = 1

    for scan_pos, scan in enumerate(scans):
        scan_idx = scan["scan_index"]
        rt = scan["rt_sec"]
        peaks = scan["peaks"]

        active_tracks = []
        for track in tracks:
            last = track["obs"][-1]
            scan_distance = scan_idx - last["scan_index"]
            if 1 <= scan_distance <= MAX_GAP_SCANS + 1:
                active_tracks.append(track)

        # Candidate links are restricted by m/z continuity. Convert the
        # permitted step to a scan-adjusted allowance for missing scans.
        links_by_track = defaultdict(list)
        links_by_peak = defaultdict(list)

        for track in active_tracks:
            last = track["obs"][-1]
            scan_distance = scan_idx - last["scan_index"]
            allowed_mda = MAX_STEP_MDA * scan_distance

            for peak_i, (mz, intensity) in enumerate(peaks):
                delta_mda = (mz - last["mz"]) * 1000.0
                if abs(delta_mda) <= allowed_mda:
                    candidate = {
                        "track": track,
                        "peak_i": peak_i,
                        "abs_delta_mda": abs(delta_mda),
                        "signed_delta_mda": delta_mda,
                        "scan_distance": scan_distance,
                    }
                    links_by_track[track["id"]].append(candidate)
                    links_by_peak[peak_i].append(candidate)

        accepted = []
        used_track_ids = set()
        used_peak_indices = set()

        # Accept only unique closest links in both directions.
        for track in active_tracks:
            candidates = links_by_track.get(track["id"], [])
            if not candidates:
                continue
            best_delta = min(x["abs_delta_mda"] for x in candidates)
            best = [
                x for x in candidates
                if abs(x["abs_delta_mda"] - best_delta) < 1e-9
            ]
            if len(best) != 1:
                ambiguities.append({
                    "scan_index": scan_idx,
                    "rt_sec": fmt(rt, 5),
                    "track_id": track["id"],
                    "ambiguity_type": "track_has_tied_peak_candidates",
                    "candidate_peak_indices": ",".join(
                        str(x["peak_i"] + 1) for x in best
                    ),
                    "candidate_mz": ",".join(
                        fmt(peaks[x["peak_i"]][0], 6) for x in best
                    ),
                    "detail": "equal closest m/z distances",
                })
                continue

            candidate = best[0]
            peak_candidates = links_by_peak[candidate["peak_i"]]
            peak_best_delta = min(
                x["abs_delta_mda"] for x in peak_candidates
            )
            peak_best = [
                x for x in peak_candidates
                if abs(x["abs_delta_mda"] - peak_best_delta) < 1e-9
            ]
            if len(peak_best) != 1:
                ambiguities.append({
                    "scan_index": scan_idx,
                    "rt_sec": fmt(rt, 5),
                    "track_id": track["id"],
                    "ambiguity_type": "peak_has_tied_track_candidates",
                    "candidate_peak_indices": str(candidate["peak_i"] + 1),
                    "candidate_mz": fmt(peaks[candidate["peak_i"]][0], 6),
                    "detail": "multiple tracks equally close to this peak",
                })
                continue

            if peak_best[0]["track"]["id"] != track["id"]:
                continue

            accepted.append(candidate)
            used_track_ids.add(track["id"])
            used_peak_indices.add(candidate["peak_i"])

        # Record cases where a peak or track had competing links, even if
        # mutual-closest logic did not accept a link.
        for peak_i, candidates in links_by_peak.items():
            distinct_tracks = sorted({
                c["track"]["id"] for c in candidates
            })
            if len(distinct_tracks) > 1:
                ambiguities.append({
                    "scan_index": scan_idx,
                    "rt_sec": fmt(rt, 5),
                    "track_id": "",
                    "ambiguity_type": "peak_competes_for_multiple_tracks",
                    "candidate_peak_indices": str(peak_i + 1),
                    "candidate_mz": fmt(peaks[peak_i][0], 6),
                    "detail": "candidate track IDs=" + ",".join(map(str, distinct_tracks)),
                })

        # Append accepted observations.
        for candidate in accepted:
            track = candidate["track"]
            peak_i = candidate["peak_i"]
            mz, intensity = peaks[peak_i]
            previous = track["obs"][-1]
            track["obs"].append({
                "scan_index": scan_idx,
                "rt_sec": rt,
                "mz": mz,
                "intensity": intensity,
                "scan_local_peak_rank": peak_i + 1,
                "delta_from_previous_mDa": candidate["signed_delta_mda"],
                "scan_distance": candidate["scan_distance"],
            })

        # Start a new track for every peak not linked to an existing track.
        for peak_i, (mz, intensity) in enumerate(peaks):
            if peak_i in used_peak_indices:
                continue
            tracks.append({
                "id": next_track_id,
                "obs": [{
                    "scan_index": scan_idx,
                    "rt_sec": rt,
                    "mz": mz,
                    "intensity": intensity,
                    "scan_local_peak_rank": peak_i + 1,
                    "delta_from_previous_mDa": None,
                    "scan_distance": 0,
                }],
            })
            next_track_id += 1

    # Summarise all tracks, including short tracks for transparency.
    obs_rows = []
    summary_rows = []
    for track in tracks:
        obs = track["obs"]
        mzs = [x["mz"] for x in obs]
        rts = [x["rt_sec"] for x in obs]
        deltas = [
            abs(x["delta_from_previous_mDa"])
            for x in obs
            if x["delta_from_previous_mDa"] is not None
        ]
        intensity_sum = sum(x["intensity"] for x in obs)

        for x in obs:
            obs_rows.append({
                "track_id": track["id"],
                "scan_index": x["scan_index"],
                "rt_sec": fmt(x["rt_sec"], 5),
                "scan_local_peak_rank": x["scan_local_peak_rank"],
                "mz": fmt(x["mz"], 6),
                "intensity": fmt(x["intensity"], 3),
                "delta_from_previous_mDa": fmt(
                    x["delta_from_previous_mDa"], 3
                ),
                "scan_distance_from_previous": x["scan_distance"],
                "track_n_observations": len(obs),
                "track_mz_min": fmt(min(mzs), 6),
                "track_mz_max": fmt(max(mzs), 6),
            })

        summary_rows.append({
            "track_id": track["id"],
            "n_observations": len(obs),
            "first_scan_index": obs[0]["scan_index"],
            "last_scan_index": obs[-1]["scan_index"],
            "rt_min_sec": fmt(min(rts), 5),
            "rt_max_sec": fmt(max(rts), 5),
            "median_mz": fmt(statistics.median(mzs), 6),
            "mz_min": fmt(min(mzs), 6),
            "mz_max": fmt(max(mzs), 6),
            "mz_range_mDa": fmt((max(mzs) - min(mzs)) * 1000.0, 3),
            "median_abs_step_mDa": fmt(
                statistics.median(deltas) if deltas else None, 3
            ),
            "max_abs_step_mDa": fmt(max(deltas) if deltas else None, 3),
            "intensity_sum": fmt(intensity_sum, 3),
            "track_status": (
                "short_provisional_track"
                if len(obs) < MIN_TRACK_OBSERVATIONS
                else "provisional_mz_continuity_track"
            ),
            "note": (
                "greedy continuity track; not a validated feature or isotope index"
            ),
        })

    write_tsv(OBS_OUT, obs_rows, [
        "track_id", "scan_index", "rt_sec", "scan_local_peak_rank",
        "mz", "intensity", "delta_from_previous_mDa",
        "scan_distance_from_previous", "track_n_observations",
        "track_mz_min", "track_mz_max",
    ])
    write_tsv(TRACK_OUT, summary_rows, [
        "track_id", "n_observations", "first_scan_index", "last_scan_index",
        "rt_min_sec", "rt_max_sec", "median_mz", "mz_min", "mz_max",
        "mz_range_mDa", "median_abs_step_mDa", "max_abs_step_mDa",
        "intensity_sum", "track_status", "note",
    ])
    write_tsv(AMBIG_OUT, ambiguities, [
        "scan_index", "rt_sec", "track_id", "ambiguity_type",
        "candidate_peak_indices", "candidate_mz", "detail",
    ])

    long_tracks = [
        x for x in summary_rows
        if int(x["n_observations"]) >= MIN_TRACK_OBSERVATIONS
    ]
    long_tracks.sort(key=lambda x: (
        -int(x["n_observations"]), float(x["median_mz"])
    ))

    print(f"MS1 scans processed: {len(scans)}")
    print(f"Centroid observations exported: {len(obs_rows)}")
    print(f"Provisional tracks: {len(tracks)}")
    print(f"Tracks with >= {MIN_TRACK_OBSERVATIONS} observations: {len(long_tracks)}")
    print(f"Ambiguity records: {len(ambiguities)}")

    print("\nTop provisional tracks by observation count")
    print("track_id  n_obs  RT_min      RT_max      median_mz    "
          "m/z_range_mDa  median_abs_step_mDa")
    for row in long_tracks[:30]:
        print(
            f"{int(row['track_id']):>8}  {int(row['n_observations']):>5}  "
            f"{float(row['rt_min_sec']):>10.5f}  "
            f"{float(row['rt_max_sec']):>10.5f}  "
            f"{float(row['median_mz']):>11.6f}  "
            f"{float(row['mz_range_mDa']):>13.3f}  "
            f"{float(row['median_abs_step_mDa']):>19.3f}"
        )

    print("\nFiles written:")
    print(f"  {OBS_OUT}")
    print(f"  {TRACK_OUT}")
    print(f"  {AMBIG_OUT}")
    print(
        "\nCaution: greedy m/z continuity is not equivalent to OpenMS feature "
        "detection. Closely spaced peaks, missing scans, and competing tracks "
        "can still cause split or incorrect tracks. Inspect these results "
        "before using any track in an isotope-group analysis."
    )


if __name__ == "__main__":
    main()
