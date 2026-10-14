#!/usr/bin/env python3
"""Audit chromatographic co-elution of persistent native-centroid tracks."""

from pathlib import Path
import csv
import math
import statistics

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
INPUT = REPORTS / "local_centroid_tracks_v1_detail.tsv"

FEATURE_ID = "10275997325121874585"
SCANS = [5852, 5865, 5878, 5891, 5904, 5917]
TRACK_IDS = [2, 3, 4, 5, 6, 7]


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def fmt(x):
    return "" if x is None else f"{x:.6f}"


if not INPUT.exists():
    raise FileNotFoundError(INPUT)

source = read_tsv(INPUT)

# One centroid observation per scan and track is expected.
data = {}
for row in source:
    scan = int(row["scan_index_zero_based"])
    track = int(row["track_id"])
    if scan not in SCANS or track not in TRACK_IDS:
        continue

    key = (scan, track)
    if key in data:
        raise ValueError(f"Duplicate scan/track observation: {key}")

    data[key] = {
        "rt": float(row["rt_sec"]),
        "mz": float(row["centroid_mz"]),
        "intensity": float(row["centroid_intensity"]),
    }

# Track-level coverage and retention-time summaries.
track_summary = []
for track in TRACK_IDS:
    observations = [
        (scan, data[(scan, track)])
        for scan in SCANS if (scan, track) in data
    ]
    observations.sort(key=lambda x: x[0])

    scans = [scan for scan, _ in observations]
    rts = [obs["rt"] for _, obs in observations]
    intensities = [obs["intensity"] for _, obs in observations]

    total_intensity = sum(intensities)
    weighted_rt = (
        sum(obs["rt"] * obs["intensity"] for _, obs in observations)
        / total_intensity
        if total_intensity > 0 else None
    )

    gaps = [
        scans[i + 1] - scans[i]
        for i in range(len(scans) - 1)
    ]
    rt_steps = [
        rts[i + 1] - rts[i]
        for i in range(len(rts) - 1)
    ]

    track_summary.append({
        "track_id": track,
        "n_observations": len(observations),
        "n_missing_of_selected_scans": len(SCANS) - len(observations),
        "first_scan_index": scans[0] if scans else "",
        "last_scan_index": scans[-1] if scans else "",
        "first_rt_sec": fmt(rts[0]) if rts else "",
        "last_rt_sec": fmt(rts[-1]) if rts else "",
        "rt_range_sec": fmt(max(rts) - min(rts)) if rts else "",
        "median_rt_sec": fmt(statistics.median(rts)) if rts else "",
        "intensity_weighted_rt_sec": fmt(weighted_rt),
        "median_intensity": fmt(statistics.median(intensities))
        if intensities else "",
        "median_consecutive_rt_step_sec": fmt(statistics.median(rt_steps))
        if rt_steps else "",
        "median_scan_index_gap": fmt(statistics.median(gaps))
        if gaps else "",
    })

write_tsv(
    REPORTS / "track_chromatography_summary_v1.tsv",
    track_summary,
)

# Pairwise overlap is calculated from actual observed scan indices and RT
# intervals. Interval overlap is descriptive and does not imply identity.
pair_summary = []
for i, track_a in enumerate(TRACK_IDS):
    for track_b in TRACK_IDS[i + 1:]:
        obs_a = {
            scan: obs for scan in SCANS
            if (scan, track_a) in data
            for obs in [data[(scan, track_a)]]
        }
        obs_b = {
            scan: obs for scan in SCANS
            if (scan, track_b) in data
            for obs in [data[(scan, track_b)]]
        }

        shared_scans = sorted(set(obs_a) & set(obs_b))
        rt_diffs = [
            abs(obs_a[scan]["rt"] - obs_b[scan]["rt"])
            for scan in shared_scans
        ]

        rts_a = [v["rt"] for v in obs_a.values()]
        rts_b = [v["rt"] for v in obs_b.values()]

        if rts_a and rts_b:
            overlap_sec = max(
                0.0,
                min(max(rts_a), max(rts_b))
                - max(min(rts_a), min(rts_b))
            )
            union_start = min(min(rts_a), min(rts_b))
            union_end = max(max(rts_a), max(rts_b))
            union_sec = union_end - union_start
            overlap_fraction = (
                overlap_sec / union_sec if union_sec > 0 else None
            )
        else:
            overlap_sec = None
            overlap_fraction = None

        pair_summary.append({
            "track_a": track_a,
            "track_b": track_b,
            "n_shared_scans": len(shared_scans),
            "shared_scan_indices": ",".join(map(str, shared_scans)),
            "median_abs_rt_difference_sec": fmt(
                statistics.median(rt_diffs) if rt_diffs else None
            ),
            "max_abs_rt_difference_sec": fmt(
                max(rt_diffs) if rt_diffs else None
            ),
            "rt_interval_overlap_sec": fmt(overlap_sec),
            "rt_interval_overlap_fraction_of_union": fmt(overlap_fraction),
        })

write_tsv(
    REPORTS / "track_chromatography_pairwise_v1.tsv",
    pair_summary,
)

# Scan-level table makes missing observations explicit.
scan_rows = []
for scan in SCANS:
    observed_rts = [
        data[(scan, track)]["rt"]
        for track in TRACK_IDS if (scan, track) in data
    ]
    scan_rows.append({
        "scan_index_zero_based": scan,
        "n_tracks_observed": len(observed_rts),
        "min_rt_sec": fmt(min(observed_rts) if observed_rts else None),
        "max_rt_sec": fmt(max(observed_rts) if observed_rts else None),
        "rt_span_sec": fmt(
            max(observed_rts) - min(observed_rts)
            if observed_rts else None
        ),
        **{
            f"track_{track}_rt_sec": fmt(data[(scan, track)]["rt"])
            if (scan, track) in data else ""
            for track in TRACK_IDS
        },
    })

write_tsv(
    REPORTS / "track_chromatography_by_scan_v1.tsv",
    scan_rows,
)

print(f"Feature: {FEATURE_ID}")
print(f"Input: {INPUT}")
print("\nTrack coverage and retention-time summaries:")
for row in track_summary:
    print(
        f"track {row['track_id']}: "
        f"observations={row['n_observations']}/{len(SCANS)}, "
        f"RT range={row['rt_range_sec']} s, "
        f"weighted RT={row['intensity_weighted_rt_sec']} s, "
        f"median intensity={row['median_intensity']}"
    )

print("\nPairwise co-elution summaries:")
for row in pair_summary:
    print(
        f"tracks {row['track_a']}-{row['track_b']}: "
        f"shared scans={row['n_shared_scans']}, "
        f"median |RT difference|="
        f"{row['median_abs_rt_difference_sec']} s, "
        f"RT interval overlap="
        f"{row['rt_interval_overlap_sec']} s"
    )

print("\nReports:")
for name in [
    "track_chromatography_summary_v1.tsv",
    "track_chromatography_pairwise_v1.tsv",
    "track_chromatography_by_scan_v1.tsv",
]:
    print(REPORTS / name)

print(
    "\nCaution: the selected observations are six scans from a narrow "
    "retention-time window. Co-elution is supportive evidence only; "
    "it does not prove isotope identity or establish an independent IG."
)
