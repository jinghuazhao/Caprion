#!/usr/bin/env python3
"""Track local native MS1 centroids across scans by m/z proximity."""

from pathlib import Path
import csv
import statistics

ROOT = Path.cwd()
INPUT = ROOT / "reports/local_ms1_centroids_feature10275997325121874585_v1.tsv"
REPORTS = ROOT / "reports"
REPORTS.mkdir(exist_ok=True)

MZ_TOLERANCE = 0.015  # Da; diagnostic matching tolerance
MIN_SCANS = 3

with INPUT.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

for r in rows:
    r["scan"] = int(r["scan_index_zero_based"])
    r["mz"] = float(r["centroid_mz"])
    r["intensity"] = float(r["centroid_intensity"])

# Build tracks sequentially. Each centroid can join the nearest existing
# track that has no observation in its scan and is within the m/z tolerance.
tracks = []
for scan in sorted({r["scan"] for r in rows}):
    scan_rows = sorted(
        (r for r in rows if r["scan"] == scan),
        key=lambda r: r["mz"]
    )
    for r in scan_rows:
        candidates = []
        for j, track in enumerate(tracks):
            if any(x["scan"] == scan for x in track):
                continue
            last_mz = statistics.median(x["mz"] for x in track)
            delta = abs(r["mz"] - last_mz)
            if delta <= MZ_TOLERANCE:
                candidates.append((delta, j))
        if candidates:
            _, j = min(candidates)
            tracks[j].append(r)
        else:
            tracks.append([r])

summary = []
detail = []
for track_id, track in enumerate(tracks, start=1):
    scans = sorted({r["scan"] for r in track})
    mzs = [r["mz"] for r in track]
    intensities = [r["intensity"] for r in track]
    summary.append({
        "track_id": track_id,
        "n_scans": len(scans),
        "scan_indices": ",".join(map(str, scans)),
        "median_mz": f"{statistics.median(mzs):.9f}",
        "mz_range_mDa": f"{(max(mzs) - min(mzs)) * 1000:.3f}",
        "median_intensity": f"{statistics.median(intensities):.3f}",
        "max_intensity": f"{max(intensities):.3f}",
        "min_intensity": f"{min(intensities):.3f}",
        "persistent_at_least_min_scans": len(scans) >= MIN_SCANS,
    })
    for r in track:
        detail.append({
            "track_id": track_id,
            "scan_index_zero_based": r["scan"],
            "native_scan_number": r["native_scan_number"],
            "rt_sec": r["rt_sec"],
            "centroid_mz": f"{r['mz']:.9f}",
            "centroid_intensity": f"{r['intensity']:.3f}",
            "error_to_nearest_model_mDa": r["error_to_nearest_model_mDa"],
            "near_trace5_or_trace6_observation":
                r["near_trace5_or_trace6_observation"],
        })

summary.sort(key=lambda r: (
    -r["n_scans"], float(r["median_mz"])
))

def write_tsv(path, data):
    if not data:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(data[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(data)

summary_path = REPORTS / "local_centroid_tracks_v1_summary.tsv"
detail_path = REPORTS / "local_centroid_tracks_v1_detail.tsv"
write_tsv(summary_path, summary)
write_tsv(detail_path, detail)

print(f"Input centroids: {len(rows)}")
print(f"Tracks formed: {len(tracks)}")
print(f"m/z matching tolerance: ±{MZ_TOLERANCE:.3f} Da")
print(f"Summary: {summary_path}")
print(f"Detail: {detail_path}")
print("\nTracks with at least three scan observations:")
for r in summary:
    if r["n_scans"] >= MIN_SCANS:
        print(
            f"track {r['track_id']}: scans={r['n_scans']}, "
            f"median m/z={r['median_mz']}, "
            f"range={r['mz_range_mDa']} mDa, "
            f"median intensity={r['median_intensity']}"
        )

print(
    "\nCaution: sequential nearest-neighbour tracking is exploratory. "
    "Tracks can merge or split, and persistence alone does not establish "
    "an isotope relationship."
)
