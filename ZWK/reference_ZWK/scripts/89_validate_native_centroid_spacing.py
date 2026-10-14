#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
FEATURE = "10275997325121874585"
THEORETICAL = 1.00335483507 / 2
TOLERANCE_DA = 0.015
MDA = 1000.0

TRACE = REPORTS / "full_trace_chromatography_v1.tsv"
CENTROIDS = REPORTS / f"local_ms1_centroids_feature{FEATURE}_v1.tsv"
MAPPING = REPORTS / "native_spectrum_index_mapping_v2.tsv"

DETAIL = REPORTS / "native_centroid_spacing_validation_v2.tsv"
SUMMARY = REPORTS / "native_centroid_spacing_validation_summary_v2.tsv"


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        print(f"No rows to write: {path}")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


# Map zero-based spectrum ordinal to native scan and MS level.
mapping_rows = read_tsv(MAPPING)
mapping = {}
for r in mapping_rows:
    idx = int(r["spectrum_ordinal_zero_based"])
    mapping[idx] = {
        "native_scan_number": int(r["native_scan_number"]),
        "ms_level": int(r["ms_level"]),
        "rt_sec": float(r["rt_sec"]),
    }

# Centroids indexed by the shared zero-based spectrum ordinal.
centroid_rows = read_tsv(CENTROIDS)
centroids = defaultdict(list)
for r in centroid_rows:
    if r["feature_id"] != FEATURE:
        continue
    idx = int(r["scan_index_zero_based"])
    centroids[idx].append({
        "native_scan_number": int(r["native_scan_number"]),
        "ms_level": int(r["ms_level"]),
        "mz": float(r["centroid_mz"]),
        "intensity": float(r["centroid_intensity"]),
    })

# Trace observations indexed by scan and trace index.
trace_rows = read_tsv(TRACE)
traces = defaultdict(dict)
for r in trace_rows:
    if r["feature_id"] != FEATURE:
        continue
    idx = int(r["scan_index_zero_based"])
    ti = int(r["trace_index"])
    traces[idx][ti] = {
        "mz": float(r["observed_mz"]),
        "rt_sec": float(r["rt_sec"]),
        "intensity": float(r["peak_intensity"]),
    }

detail = []
for scan_idx in sorted(traces):
    obs = traces[scan_idx]
    if scan_idx not in mapping:
        map_scan = ""
        map_level = ""
        map_status = "NO_MAPPING"
    else:
        map_scan = mapping[scan_idx]["native_scan_number"]
        map_level = mapping[scan_idx]["ms_level"]
        map_status = "OK" if map_level == 1 else "NOT_MS1"

    native_scan = map_scan
    if scan_idx in centroids and centroids[scan_idx]:
        centroid_scan_values = {c["native_scan_number"] for c in centroids[scan_idx]}
        centroid_level_values = {c["ms_level"] for c in centroids[scan_idx]}
        if len(centroid_scan_values) != 1 or len(centroid_level_values) != 1:
            map_status = "CENTROID_METADATA_INCONSISTENT"
        elif map_scan != "" and int(map_scan) not in centroid_scan_values:
            map_status = "NATIVE_SCAN_MISMATCH"
        elif map_level != "" and int(map_level) not in centroid_level_values:
            map_status = "MS_LEVEL_MISMATCH"

    for lower in sorted(obs):
        upper = lower + 1
        if upper not in obs:
            continue

        lo = obs[lower]["mz"]
        hi = obs[upper]["mz"]
        trace_spacing = hi - lo
        trace_resid = (trace_spacing - THEORETICAL) * MDA

        scan_centroids = centroids.get(scan_idx, [])
        lo_candidates = [
            c for c in scan_centroids if abs(c["mz"] - lo) <= TOLERANCE_DA
        ]
        hi_candidates = [
            c for c in scan_centroids if abs(c["mz"] - hi) <= TOLERANCE_DA
        ]

        # Count unique centroid pairs only when both trace peaks have
        # exactly one candidate and those candidates are distinct.
        unique_pair = len(lo_candidates) == 1 and len(hi_candidates) == 1
        if unique_pair:
            clo, chi = lo_candidates[0], hi_candidates[0]
            unique_pair = (
                clo["mz"] != chi["mz"]
                and clo["native_scan_number"] == chi["native_scan_number"]
                and clo["ms_level"] == chi["ms_level"] == 1
            )

        native_spacing = ""
        native_resid = ""
        native_pass = ""
        lower_delta = ""
        upper_delta = ""
        if unique_pair:
            native_spacing = chi["mz"] - clo["mz"]
            native_resid = (native_spacing - THEORETICAL) * MDA
            native_pass = abs(native_resid) <= 3.0
            lower_delta = (clo["mz"] - lo) * MDA
            upper_delta = (chi["mz"] - hi) * MDA

        detail.append({
            "feature_id": FEATURE,
            "scan_index_zero_based": scan_idx,
            "native_scan_number": native_scan,
            "mapping_status": map_status,
            "rt_sec": f'{obs[lower]["rt_sec"]:.6f}',
            "lower_trace_index": lower,
            "upper_trace_index": upper,
            "trace_observed_spacing_Da": f"{trace_spacing:.9f}",
            "trace_spacing_residual_mDa": f"{trace_resid:.6f}",
            "trace_spacing_within_3mDa": abs(trace_resid) <= 3.0,
            "lower_centroid_candidates": len(lo_candidates),
            "upper_centroid_candidates": len(hi_candidates),
            "unique_centroid_pair": unique_pair,
            "lower_trace_centroid_delta_mDa": lower_delta,
            "upper_trace_centroid_delta_mDa": upper_delta,
            "native_centroid_spacing_Da": (
                f"{native_spacing:.9f}" if native_spacing != "" else ""
            ),
            "native_spacing_residual_mDa": (
                f"{native_resid:.6f}" if native_resid != "" else ""
            ),
            "native_spacing_within_3mDa": native_pass,
        })

write_tsv(DETAIL, detail)

by_pair = defaultdict(list)
for r in detail:
    by_pair[(r["lower_trace_index"], r["upper_trace_index"])].append(r)

summary = []
for (lower, upper), rows in sorted(by_pair.items()):
    trace_pass = sum(r["trace_spacing_within_3mDa"] is True
                     or str(r["trace_spacing_within_3mDa"]).lower() == "true"
                     for r in rows)
    unique = [r for r in rows if r["unique_centroid_pair"] is True
              or str(r["unique_centroid_pair"]).lower() == "true"]
    native_pass = sum(
        str(r["native_spacing_within_3mDa"]).lower() == "true"
        for r in unique
    )
    summary.append({
        "feature_id": FEATURE,
        "lower_trace_index": lower,
        "upper_trace_index": upper,
        "n_trace_pairs": len(rows),
        "trace_spacing_within_3mDa": f"{trace_pass}/{len(rows)}",
        "scans_with_centroids": sum(
            int(r["lower_centroid_candidates"]) > 0
            and int(r["upper_centroid_candidates"]) > 0 for r in rows
        ),
        "unique_centroid_match_pairs": f"{len(unique)}/{len(rows)}",
        "native_spacing_within_3mDa": (
            f"{native_pass}/{len(unique)}" if unique else "NA"
        ),
        "mapping_statuses": ";".join(sorted(set(r["mapping_status"] for r in rows))),
    })

write_tsv(SUMMARY, summary)

print(f"Feature: {FEATURE}")
print(f"Detail: {DETAIL}")
print(f"Summary: {SUMMARY}")
print("\nPair summary:")
for r in summary:
    print(
        f'{r["lower_trace_index"]}-{r["upper_trace_index"]}: '
        f'trace ±3 mDa={r["trace_spacing_within_3mDa"]}; '
        f'centroid candidates in both windows={r["scans_with_centroids"]}/'
        f'{r["n_trace_pairs"]}; '
        f'unique centroid pairs={r["unique_centroid_match_pairs"]}; '
        f'native spacing ±3 mDa={r["native_spacing_within_3mDa"]}; '
        f'mapping={r["mapping_statuses"]}'
    )

print("\nCaution: centroid correspondence is not independent proof of isotope identity.")
