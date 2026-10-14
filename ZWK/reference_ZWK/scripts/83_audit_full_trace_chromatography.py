#!/usr/bin/env python3
"""Audit full RT and intensity profiles for one OpenMS feature."""

from pathlib import Path
import csv
import math
import statistics

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
INPUT = ROOT / "openms_pilot/isotope_trace_peaks.tsv"

FEATURE_ID = "10275997325121874585"
TRACE_INDICES = [1, 2, 3, 4, 5, 6]


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


def pearson(x, y):
    if len(x) < 3:
        return None
    mx, my = statistics.mean(x), statistics.mean(y)
    dx = [v - mx for v in x]
    dy = [v - my for v in y]
    denom = math.sqrt(
        sum(v * v for v in dx) * sum(v * v for v in dy)
    )
    return sum(a * b for a, b in zip(dx, dy)) / denom if denom else None


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


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


if not INPUT.exists():
    raise FileNotFoundError(INPUT)

source = read_tsv(INPUT)

# Preserve all observations in the input; use scan index as the observation
# identifier and RT as the chromatographic coordinate.
data = {}
for row in source:
    if row["feature_id"] != FEATURE_ID:
        continue
    trace = int(row["trace_index"])
    if trace not in TRACE_INDICES:
        continue

    scan = int(row["scan_index"])
    key = (trace, scan)
    if key in data:
        raise ValueError(f"Duplicate trace/scan observation: {key}")

    data[key] = {
        "rt": float(row["rt_sec"]),
        "mz": float(row["observed_mz"]),
        "intensity": float(row["peak_intensity"]),
    }

# Full plot-ready table.
profile_rows = []
for trace in TRACE_INDICES:
    observations = sorted(
        [
            (scan, values)
            for (idx, scan), values in data.items()
            if idx == trace
        ],
        key=lambda item: (item[1]["rt"], item[0]),
    )

    for scan, obs in observations:
        profile_rows.append({
            "feature_id": FEATURE_ID,
            "trace_index": trace,
            "scan_index_zero_based": scan,
            "rt_sec": fmt(obs["rt"]),
            "observed_mz": fmt(obs["mz"]),
            "peak_intensity": fmt(obs["intensity"]),
            "log10_intensity": fmt(
                math.log10(obs["intensity"])
                if obs["intensity"] > 0 else None
            ),
        })

write_tsv(
    REPORTS / "full_trace_chromatography_v1.tsv",
    profile_rows,
)

# Trace-level summaries.
summary = []
for trace in TRACE_INDICES:
    observations = sorted(
        [
            (scan, values)
            for (idx, scan), values in data.items()
            if idx == trace
        ],
        key=lambda item: (item[1]["rt"], item[0]),
    )

    if not observations:
        summary.append({
            "trace_index": trace,
            "n_observations": 0,
        })
        continue

    scans = [scan for scan, _ in observations]
    rts = [obs["rt"] for _, obs in observations]
    intensities = [obs["intensity"] for _, obs in observations]
    total = sum(intensities)

    apex_i = max(range(len(observations)), key=lambda i: intensities[i])
    apex_scan = scans[apex_i]
    apex_rt = rts[apex_i]
    apex_intensity = intensities[apex_i]

    weighted_rt = (
        sum(rt * intensity for rt, intensity in zip(rts, intensities)) / total
        if total > 0 else None
    )

    summary.append({
        "trace_index": trace,
        "n_observations": len(observations),
        "first_scan_index": scans[0],
        "last_scan_index": scans[-1],
        "first_rt_sec": fmt(min(rts)),
        "last_rt_sec": fmt(max(rts)),
        "rt_range_sec": fmt(max(rts) - min(rts)),
        "apex_scan_index": apex_scan,
        "apex_rt_sec": fmt(apex_rt),
        "apex_intensity": fmt(apex_intensity),
        "median_rt_sec": fmt(statistics.median(rts)),
        "intensity_weighted_rt_sec": fmt(weighted_rt),
        "sum_intensity": fmt(total),
        "median_intensity": fmt(statistics.median(intensities)),
    })

write_tsv(
    REPORTS / "full_trace_chromatography_summary_v1.tsv",
    summary,
)

# Correlations only on exact shared scan indices. Do not interpolate absent
# observations, because that could manufacture artificial agreement.
correlations = []
for i, trace_a in enumerate(TRACE_INDICES):
    for trace_b in TRACE_INDICES[i + 1:]:
        scans_a = {
            scan for (trace, scan) in data if trace == trace_a
        }
        scans_b = {
            scan for (trace, scan) in data if trace == trace_b
        }
        shared = sorted(scans_a & scans_b)

        xa = [data[(trace_a, scan)]["intensity"] for scan in shared]
        xb = [data[(trace_b, scan)]["intensity"] for scan in shared]

        valid = [
            (a, b) for a, b in zip(xa, xb) if a > 0 and b > 0
        ]
        x = [a for a, _ in valid]
        y = [b for _, b in valid]

        r_raw = pearson(x, y)
        rho = spearman(x, y)
        r_log = (
            pearson([math.log10(v) for v in x],
                    [math.log10(v) for v in y])
            if len(x) >= 3 else None
        )

        correlations.append({
            "trace_a": trace_a,
            "trace_b": trace_b,
            "n_shared_scans": len(shared),
            "n_positive_intensity_pairs": len(valid),
            "pearson_raw": fmt(r_raw),
            "spearman": fmt(rho),
            "pearson_log10": fmt(r_log),
        })

write_tsv(
    REPORTS / "full_trace_intensity_correlations_v1.tsv",
    correlations,
)

# Apex RT offsets relative to trace 1.
apex_by_trace = {
    int(row["trace_index"]): float(row["apex_rt_sec"])
    for row in summary
    if row.get("apex_rt_sec") not in (None, "")
}

apex_offsets = []
if 1 in apex_by_trace:
    for trace in TRACE_INDICES:
        if trace not in apex_by_trace:
            continue
        apex_offsets.append({
            "trace_index": trace,
            "apex_rt_sec": fmt(apex_by_trace[trace]),
            "apex_offset_vs_trace1_sec": fmt(
                apex_by_trace[trace] - apex_by_trace[1]
            ),
        })

write_tsv(
    REPORTS / "full_trace_apex_offsets_v1.tsv",
    apex_offsets,
)

print(f"Feature: {FEATURE_ID}")
print(f"Input: {INPUT}")
print(f"Trace observations retained: {len(profile_rows)}")

print("\nFull-profile summaries:")
for row in summary:
    if row.get("n_observations", 0) == 0:
        print(f"trace {row['trace_index']}: no observations")
        continue
    print(
        f"trace {row['trace_index']}: "
        f"n={row['n_observations']}, "
        f"RT={row['first_rt_sec']}–{row['last_rt_sec']} s, "
        f"apex={row['apex_rt_sec']} s, "
        f"weighted RT={row['intensity_weighted_rt_sec']} s"
    )

print("\nPairwise full-profile intensity correlations:")
for row in correlations:
    print(
        f"traces {row['trace_a']}-{row['trace_b']}: "
        f"n={row['n_positive_intensity_pairs']}, "
        f"Pearson raw={row['pearson_raw']}, "
        f"Spearman={row['spearman']}, "
        f"Pearson log10={row['pearson_log10']}"
    )

print("\nReports:")
for name in [
    "full_trace_chromatography_v1.tsv",
    "full_trace_chromatography_summary_v1.tsv",
    "full_trace_intensity_correlations_v1.tsv",
    "full_trace_apex_offsets_v1.tsv",
]:
    print(REPORTS / name)

print(
    "\nCaution: profile correlations and apex alignment are descriptive. "
    "They do not establish isotope identity or an independent IG."
)
