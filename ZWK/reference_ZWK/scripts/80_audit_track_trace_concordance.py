#!/usr/bin/env python3
"""Audit concordance between native-centroid tracks and OpenMS trace peaks."""

from pathlib import Path
import csv
import statistics

ROOT = Path.cwd()
REPORTS = ROOT / "reports"

TRACE_FILE = ROOT / "openms_pilot/isotope_trace_peaks.tsv"
TRACK_FILE = REPORTS / "local_centroid_tracks_v1_detail.tsv"

FEATURE_ID = "10275997325121874585"
SCANS = [5852, 5865, 5878, 5891, 5904, 5917]
TRACE_INDICES = [1, 2, 3, 4, 5, 6]

# A match within 5 mDa is a useful diagnostic threshold here, not proof
# of isotope identity.
MATCH_TOLERANCE_MDA = 5.0


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


if not TRACE_FILE.exists():
    raise FileNotFoundError(TRACE_FILE)
if not TRACK_FILE.exists():
    raise FileNotFoundError(TRACK_FILE)

trace_rows = read_tsv(TRACE_FILE)
track_rows = read_tsv(TRACK_FILE)

traces = {}
for row in trace_rows:
    if row["feature_id"] != FEATURE_ID:
        continue
    scan = int(row["scan_index"])
    idx = int(row["trace_index"])
    if scan not in SCANS or idx not in TRACE_INDICES:
        continue
    traces[(scan, idx)] = {
        "mz": float(row["observed_mz"]),
        "intensity": float(row["peak_intensity"]),
        "rt": float(row["rt_sec"]),
    }

tracks = {}
for row in track_rows:
    scan = int(row["scan_index_zero_based"])
    track_id = int(row["track_id"])
    tracks[(scan, track_id)] = {
        "mz": float(row["centroid_mz"]),
        "intensity": float(row["centroid_intensity"]),
        "rt": float(row["rt_sec"]),
    }

# Only include tracks observed in at least three of the selected scans.
track_scan_counts = {}
for scan, track_id in tracks:
    track_scan_counts[track_id] = track_scan_counts.get(track_id, 0) + 1
persistent_ids = sorted(
    track_id for track_id, n in track_scan_counts.items() if n >= 3
)

matches = []
for scan in SCANS:
    for trace_idx in TRACE_INDICES:
        trace = traces.get((scan, trace_idx))
        if trace is None:
            continue

        candidates = []
        for track_id in persistent_ids:
            centroid = tracks.get((scan, track_id))
            if centroid is None:
                continue
            delta_mda = (centroid["mz"] - trace["mz"]) * 1000.0
            candidates.append((abs(delta_mda), track_id, delta_mda, centroid))

        candidates.sort(key=lambda x: x[0])
        best = candidates[0] if candidates else None
        second = candidates[1] if len(candidates) > 1 else None

        if best is None:
            best_id, best_delta, best_intensity, best_abs = "", "", "", ""
            within = False
            ambiguous = False
            second_delta = ""
        else:
            _, best_id, best_delta, centroid = best
            best_delta = f"{best_delta:.6f}"
            best_abs = f"{abs(float(best_delta)):.6f}"
            best_intensity = f"{centroid['intensity']:.3f}"
            within = abs(float(best_delta)) <= MATCH_TOLERANCE_MDA
            second_delta = (
                f"{second[2]:.6f}" if second is not None else ""
            )
            # Flag as ambiguous if another track is also within the
            # tolerance and is nearly as close as the best match.
            ambiguous = (
                second is not None
                and second[0] <= MATCH_TOLERANCE_MDA
                and second[0] - best[0] < 1.0
            )

        matches.append({
            "feature_id": FEATURE_ID,
            "scan_index_zero_based": scan,
            "trace_index": trace_idx,
            "trace_mz": f"{trace['mz']:.9f}",
            "trace_intensity": f"{trace['intensity']:.3f}",
            "nearest_track_id": best_id,
            "nearest_track_mz_delta_mDa": best_delta,
            "nearest_track_abs_delta_mDa": best_abs,
            "nearest_track_intensity": best_intensity,
            "second_nearest_track_delta_mDa": second_delta,
            "within_5mDa": within,
            "match_ambiguous": ambiguous,
        })

summary = []
for trace_idx in TRACE_INDICES:
    rows = [r for r in matches if int(r["trace_index"]) == trace_idx]
    deltas = [
        abs(float(r["nearest_track_mz_delta_mDa"]))
        for r in rows if r["nearest_track_mz_delta_mDa"] != ""
    ]
    passing = sum(bool(r["within_5mDa"]) for r in rows)
    ambiguous = sum(bool(r["match_ambiguous"]) for r in rows)
    summary.append({
        "trace_index": trace_idx,
        "n_observations": len(rows),
        "n_within_5mDa": passing,
        "fraction_within_5mDa": (
            f"{passing / len(rows):.6f}" if rows else ""
        ),
        "median_abs_track_delta_mDa": (
            f"{statistics.median(deltas):.6f}" if deltas else ""
        ),
        "max_abs_track_delta_mDa": (
            f"{max(deltas):.6f}" if deltas else ""
        ),
        "n_ambiguous_matches": ambiguous,
    })

match_path = REPORTS / "track_trace_concordance_v1_detail.tsv"
summary_path = REPORTS / "track_trace_concordance_v1_summary.tsv"
write_tsv(match_path, matches)
write_tsv(summary_path, summary)

print(f"Feature: {FEATURE_ID}")
print(f"Persistent tracks considered: {persistent_ids}")
print(f"Matching tolerance: ±{MATCH_TOLERANCE_MDA:.1f} mDa")
print(f"Detail report: {match_path}")
print(f"Summary report: {summary_path}")
print("\nTrace-to-track concordance:")
for row in summary:
    print(
        f"trace {row['trace_index']}: "
        f"{row['n_within_5mDa']}/{row['n_observations']} within 5 mDa; "
        f"median absolute delta={row['median_abs_track_delta_mDa']} mDa; "
        f"max absolute delta={row['max_abs_track_delta_mDa']} mDa; "
        f"ambiguous={row['n_ambiguous_matches']}"
    )

print(
    "\nCaution: this is a nearest-track concordance audit. "
    "It does not prove isotope identity or establish an independent IG."
)
