#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
FEATURE = "10275997325121874585"

DETAIL_IN = REPORTS / "native_centroid_spacing_validation_v2.tsv"
TRACE_IN = REPORTS / "full_trace_chromatography_v1.tsv"
CENTROID_IN = REPORTS / f"local_ms1_centroids_feature{FEATURE}_v1.tsv"
MAPPING_IN = REPORTS / "native_spectrum_index_mapping_v2.tsv"

DETAIL_OUT = REPORTS / "native_centroid_pair_audit_v1.tsv"
SUMMARY_OUT = REPORTS / "native_centroid_pair_audit_summary_v1.tsv"
MAPPING_GAPS_OUT = REPORTS / "native_centroid_pair_mapping_gaps_v1.tsv"

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows):
    if not rows:
        print(f"No rows to write: {path}")
        return
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t",
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)

detail = read_tsv(DETAIL_IN)
trace_rows = read_tsv(TRACE_IN)
centroid_rows = read_tsv(CENTROID_IN)
mapping_rows = read_tsv(MAPPING_IN)

# Map zero-based spectrum ordinal to native scan number and MS level.
mapping = {}
for r in mapping_rows:
    idx = int(r["spectrum_ordinal_zero_based"])
    mapping[idx] = {
        "native_scan_number": int(r["native_scan_number"]),
        "ms_level": int(r["ms_level"]),
        "rt_sec": float(r["rt_sec"]),
        "native_spectrum_id": r["native_spectrum_id"],
    }

traces = defaultdict(dict)
for r in trace_rows:
    if r["feature_id"] != FEATURE:
        continue
    idx = int(r["scan_index_zero_based"])
    ti = int(r["trace_index"])
    traces[idx][ti] = {
        "mz": float(r["observed_mz"]),
        "rt_sec": float(r["rt_sec"]),
    }

centroids = defaultdict(list)
for r in centroid_rows:
    if r["feature_id"] != FEATURE:
        continue
    idx = int(r["scan_index_zero_based"])
    centroids[idx].append({
        "native_scan_number": int(r["native_scan_number"]),
        "ms_level": int(r["ms_level"]),
        "spectrum_id": r["spectrum_id"],
        "rt_sec": float(r["rt_sec"]),
        "centroid_ordinal": int(r["centroid_ordinal_in_spectrum"]),
        "mz": float(r["centroid_mz"]),
        "intensity": float(r["centroid_intensity"]),
    })

audit = []
gaps = []

for r in detail:
    idx = int(r["scan_index_zero_based"])
    lower = int(r["lower_trace_index"])
    upper = int(r["upper_trace_index"])
    lo_trace = traces.get(idx, {}).get(lower)
    hi_trace = traces.get(idx, {}).get(upper)
    cs = centroids.get(idx, [])

    maprec = mapping.get(idx)
    trace_scan = maprec["native_scan_number"] if maprec else ""
    trace_level = maprec["ms_level"] if maprec else ""
    trace_mapping_status = "OK" if maprec else "NO_MAPPING"

    # Use the native scan and MS1 level from the centroid file itself.
    native_scans = sorted({c["native_scan_number"] for c in cs})
    centroid_levels = sorted({c["ms_level"] for c in cs})
    centroid_scan = native_scans[0] if len(native_scans) == 1 else ""
    centroid_level = centroid_levels[0] if len(centroid_levels) == 1 else ""

    if not maprec:
        gap_reason = "trace_scan_index_absent_from_mapping_table"
    elif maprec["ms_level"] != 1:
        gap_reason = "mapping_row_is_not_MS1"
    elif not cs:
        gap_reason = "no_centroid_rows_for_scan_index"
    elif len(native_scans) != 1:
        gap_reason = "multiple_native_scan_numbers_in_centroid_rows"
    elif len(centroid_levels) != 1:
        gap_reason = "multiple_MS_levels_in_centroid_rows"
    elif centroid_scan != maprec["native_scan_number"]:
        gap_reason = "native_scan_number_disagreement"
    elif centroid_level != maprec["ms_level"]:
        gap_reason = "MS_level_disagreement"
    else:
        gap_reason = ""

    if gap_reason:
        gaps.append({
            "feature_id": FEATURE,
            "scan_index_zero_based": idx,
            "lower_trace_index": lower,
            "upper_trace_index": upper,
            "trace_rt_sec": f'{lo_trace["rt_sec"]:.6f}' if lo_trace else "",
            "mapping_found": bool(maprec),
            "mapped_native_scan_number": trace_scan,
            "mapped_ms_level": trace_level,
            "centroid_rows": len(cs),
            "centroid_native_scan_numbers": ",".join(map(str, native_scans)),
            "centroid_ms_levels": ",".join(map(str, centroid_levels)),
            "gap_reason": gap_reason,
        })

    # Find centroids nearest each trace m/z, sorted by absolute delta.
    lo_candidates = []
    hi_candidates = []
    if lo_trace:
        lo_candidates = sorted(
            [(abs(c["mz"] - lo_trace["mz"]), c) for c in cs],
            key=lambda x: x[0]
        )
    if hi_trace:
        hi_candidates = sorted(
            [(abs(c["mz"] - hi_trace["mz"]), c) for c in cs],
            key=lambda x: x[0]
        )

    lo_match = lo_candidates[0][1] if lo_candidates else None
    hi_match = hi_candidates[0][1] if hi_candidates else None

    lo_delta = (
        (lo_match["mz"] - lo_trace["mz"]) * 1000
        if lo_match and lo_trace else ""
    )
    hi_delta = (
        (hi_match["mz"] - hi_trace["mz"]) * 1000
        if hi_match and hi_trace else ""
    )

    native_spacing = ""
    native_residual = ""
    native_pass = ""
    if lo_match and hi_match:
        native_spacing = hi_match["mz"] - lo_match["mz"]
        native_residual = (native_spacing - 1.00335483507 / 2) * 1000
        native_pass = abs(native_residual) <= 3

    audit.append({
        "feature_id": FEATURE,
        "scan_index_zero_based": idx,
        "native_scan_number_from_mapping": trace_scan,
        "mapping_ms_level": trace_level,
        "mapping_status": trace_mapping_status,
        "trace_rt_sec": f'{lo_trace["rt_sec"]:.6f}' if lo_trace else "",
        "lower_trace_index": lower,
        "upper_trace_index": upper,
        "lower_trace_mz": f'{lo_trace["mz"]:.9f}' if lo_trace else "",
        "upper_trace_mz": f'{hi_trace["mz"]:.9f}' if hi_trace else "",
        "centroid_rows": len(cs),
        "centroid_native_scan_numbers": ",".join(map(str, native_scans)),
        "centroid_ms_levels": ",".join(map(str, centroid_levels)),
        "lower_nearest_centroid_ordinal": (
            lo_match["centroid_ordinal"] if lo_match else ""
        ),
        "lower_nearest_centroid_mz": (
            f'{lo_match["mz"]:.9f}' if lo_match else ""
        ),
        "lower_nearest_centroid_intensity": (
            f'{lo_match["intensity"]:.3f}' if lo_match else ""
        ),
        "lower_centroid_delta_mDa": (
            f'{lo_delta:.6f}' if lo_delta != "" else ""
        ),
        "upper_nearest_centroid_ordinal": (
            hi_match["centroid_ordinal"] if hi_match else ""
        ),
        "upper_nearest_centroid_mz": (
            f'{hi_match["mz"]:.9f}' if hi_match else ""
        ),
        "upper_nearest_centroid_intensity": (
            f'{hi_match["intensity"]:.3f}' if hi_match else ""
        ),
        "upper_centroid_delta_mDa": (
            f'{hi_delta:.6f}' if hi_delta != "" else ""
        ),
        "native_centroid_spacing_Da": (
            f'{native_spacing:.9f}' if native_spacing != "" else ""
        ),
        "native_spacing_residual_mDa": (
            f'{native_residual:.6f}' if native_residual != "" else ""
        ),
        "native_spacing_within_3mDa": native_pass,
        "mapping_or_centroid_issue": gap_reason,
    })

write_tsv(DETAIL_OUT, audit)
write_tsv(MAPPING_GAPS_OUT, gaps)

by_pair = defaultdict(list)
for r in audit:
    by_pair[(r["lower_trace_index"], r["upper_trace_index"])].append(r)

summary = []
for (lower, upper), rows in sorted(by_pair.items()):
    unique_rows = [
        r for r in rows
        if r["lower_nearest_centroid_ordinal"] != ""
        and r["upper_nearest_centroid_ordinal"] != ""
    ]
    passing = sum(
        str(r["native_spacing_within_3mDa"]).lower() == "true"
        for r in unique_rows
    )
    summary.append({
        "feature_id": FEATURE,
        "lower_trace_index": lower,
        "upper_trace_index": upper,
        "n_trace_pairs": len(rows),
        "n_rows_with_centroid_matches": len(unique_rows),
        "n_mapping_or_centroid_issues": sum(bool(r["mapping_or_centroid_issue"]) for r in rows),
        "native_spacing_within_3mDa": (
            f"{passing}/{len(unique_rows)}" if unique_rows else "NA"
        ),
        "native_residuals_mDa": ",".join(
            r["native_spacing_residual_mDa"] for r in unique_rows
        ),
    })

write_tsv(SUMMARY_OUT, summary)

print(f"Detail: {DETAIL_OUT}")
print(f"Summary: {SUMMARY_OUT}")
print(f"Mapping/centroid issues: {MAPPING_GAPS_OUT}")
print(f"Audit rows: {len(audit)}")
print(f"Mapping/centroid issue rows: {len(gaps)}")
print("\nPair summary:")
for r in summary:
    print(
        f'{r["lower_trace_index"]}-{r["upper_trace_index"]}: '
        f'centroid matches={r["n_rows_with_centroid_matches"]}/'
        f'{r["n_trace_pairs"]}; '
        f'native spacing within ±3 mDa={r["native_spacing_within_3mDa"]}; '
        f'issues={r["n_mapping_or_centroid_issues"]}'
    )
