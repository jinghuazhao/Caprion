#!/usr/bin/env python3
"""Audit intensity coherence across persistent centroid tracks."""

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


def pearson(x, y):
    if len(x) < 3:
        return None
    mx, my = statistics.mean(x), statistics.mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    denom = math.sqrt(sum(v*v for v in dx) * sum(v*v for v in dy))
    return sum(a*b for a, b in zip(dx, dy)) / denom if denom else None


def ranks(values):
    """Average ranks for ties."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = ((i + 1) + j) / 2.0
        for k in range(i, j):
            result[order[k]] = rank
        i = j
    return result


def fmt(value):
    return "" if value is None else f"{value:.6f}"


if not INPUT.exists():
    raise FileNotFoundError(INPUT)

rows = read_tsv(INPUT)
intensity = {}

for row in rows:
    scan = int(row["scan_index_zero_based"])
    track = int(row["track_id"])
    if scan not in SCANS or track not in TRACK_IDS:
        continue
    intensity[(scan, track)] = float(row["centroid_intensity"])

# Wide table: one row per scan, one intensity column per track.
wide = []
for scan in SCANS:
    record = {"scan_index_zero_based": scan}
    for track in TRACK_IDS:
        value = intensity.get((scan, track))
        record[f"track_{track}_intensity"] = fmt(value)
    wide.append(record)

write_tsv(REPORTS / "track_intensity_matrix_v1.tsv", wide)

# Pairwise correlations across scans.
correlations = []
for i, track_a in enumerate(TRACK_IDS):
    for track_b in TRACK_IDS[i + 1:]:
        paired = [
            (intensity[(scan, track_a)], intensity[(scan, track_b)])
            for scan in SCANS
            if (scan, track_a) in intensity and (scan, track_b) in intensity
        ]
        x = [p[0] for p in paired]
        y = [p[1] for p in paired]
        if len(paired) >= 3 and all(v > 0 for v in x + y):
            log_x = [math.log10(v) for v in x]
            log_y = [math.log10(v) for v in y]
            r_raw = pearson(x, y)
            r_log = pearson(log_x, log_y)
            rho = pearson(ranks(x), ranks(y))
        else:
            r_raw = r_log = rho = None

        correlations.append({
            "track_a": track_a,
            "track_b": track_b,
            "n_paired_scans": len(paired),
            "pearson_raw": fmt(r_raw),
            "spearman": fmt(rho),
            "pearson_log10": fmt(r_log),
        })

write_tsv(
    REPORTS / "track_intensity_correlations_v1.tsv",
    correlations,
)

# Adjacent-track intensity ratios within each scan.
ratios = []
for scan in SCANS:
    for a, b in zip(TRACK_IDS[:-1], TRACK_IDS[1:]):
        ia = intensity.get((scan, a))
        ib = intensity.get((scan, b))
        ratio = ib / ia if ia is not None and ia > 0 and ib is not None else None
        ratios.append({
            "scan_index_zero_based": scan,
            "track_lower": a,
            "track_upper": b,
            "intensity_lower": fmt(ia),
            "intensity_upper": fmt(ib),
            "upper_over_lower_ratio": fmt(ratio),
        })

write_tsv(
    REPORTS / "track_adjacent_intensity_ratios_v1.tsv",
    ratios,
)

print(f"Feature: {FEATURE_ID}")
print(f"Input: {INPUT}")
print(f"Selected scans: {SCANS}")
print(f"Tracks: {TRACK_IDS}")
print("\nPairwise intensity correlations:")
for row in correlations:
    print(
        f"tracks {row['track_a']}-{row['track_b']}: "
        f"n={row['n_paired_scans']}, "
        f"Pearson raw={row['pearson_raw']}, "
        f"Spearman={row['spearman']}, "
        f"Pearson log10={row['pearson_log10']}"
    )

print("\nMedian adjacent-track intensity ratios:")
for a, b in zip(TRACK_IDS[:-1], TRACK_IDS[1:]):
    values = [
        float(row["upper_over_lower_ratio"])
        for row in ratios
        if row["track_lower"] == a
        and row["track_upper"] == b
        and row["upper_over_lower_ratio"] != ""
    ]
    med = statistics.median(values) if values else None
    print(f"track {a}->{b}: n={len(values)}, median ratio={fmt(med)}")

print("\nReports:")
print(REPORTS / "track_intensity_matrix_v1.tsv")
print(REPORTS / "track_intensity_correlations_v1.tsv")
print(REPORTS / "track_adjacent_intensity_ratios_v1.tsv")
print(
    "\nCaution: correlations and intensity ratios are descriptive. "
    "They do not prove isotope identity or establish an independent IG."
)
