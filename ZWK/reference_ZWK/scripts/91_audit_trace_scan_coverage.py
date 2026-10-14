#!/usr/bin/env python3
import csv
from collections import Counter
from pathlib import Path

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
FEATURE = "10275997325121874585"

TRACE = REPORTS / "full_trace_chromatography_v1.tsv"
MAPPING = REPORTS / "native_spectrum_index_mapping_v2.tsv"
CENTROIDS = REPORTS / f"local_ms1_centroids_feature{FEATURE}_v1.tsv"

DETAIL = REPORTS / "trace_scan_coverage_audit_v1.tsv"
SUMMARY = REPORTS / "trace_scan_coverage_summary_v1.tsv"

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

trace_rows = read_tsv(TRACE)
mapping_rows = read_tsv(MAPPING)
centroid_rows = read_tsv(CENTROIDS)

mapping = {}
for r in mapping_rows:
    idx = int(r["spectrum_ordinal_zero_based"])
    mapping[idx] = {
        "native_scan_number": int(r["native_scan_number"]),
        "ms_level": int(r["ms_level"]),
        "rt_sec": float(r["rt_sec"]),
        "native_spectrum_id": r["native_spectrum_id"],
    }

centroid_scans = {}
for r in centroid_rows:
    if r["feature_id"] != FEATURE:
        continue
    idx = int(r["scan_index_zero_based"])
    centroid_scans.setdefault(idx, set()).add(
        (int(r["native_scan_number"]), int(r["ms_level"]))
    )

# One record per trace observation, retaining the trace index.
detail = []
for r in trace_rows:
    if r["feature_id"] != FEATURE:
        continue

    idx = int(r["scan_index_zero_based"])
    ti = int(r["trace_index"])
    rt = float(r["rt_sec"])
    mz = float(r["observed_mz"])

    m = mapping.get(idx)
    cent = centroid_scans.get(idx, set())

    if m is None:
        mapping_status = "NO_MAPPING_ROW"
        native_scan = ""
        ms_level = ""
        mapping_rt = ""
        spectrum_id = ""
    else:
        native_scan = m["native_scan_number"]
        ms_level = m["ms_level"]
        mapping_rt = f'{m["rt_sec"]:.6f}'
        spectrum_id = m["native_spectrum_id"]
        mapping_status = "MAPPED_MS1" if m["ms_level"] == 1 else "MAPPED_NON_MS1"

    if not cent:
        centroid_status = "NOT_IN_LOCAL_CENTROID_EXPORT"
    elif len(cent) > 1:
        centroid_status = "MULTIPLE_SCAN_OR_LEVEL_VALUES"
    else:
        cscan, clevel = next(iter(cent))
        if m is None:
            centroid_status = "CENTROID_PRESENT_BUT_MAPPING_MISSING"
        elif cscan != native_scan:
            centroid_status = "NATIVE_SCAN_DISAGREEMENT"
        elif clevel != ms_level:
            centroid_status = "MS_LEVEL_DISAGREEMENT"
        else:
            centroid_status = "CENTROID_EXPORT_PRESENT"

    detail.append({
        "feature_id": FEATURE,
        "scan_index_zero_based": idx,
        "trace_index": ti,
        "trace_rt_sec": f"{rt:.6f}",
        "trace_observed_mz": f"{mz:.9f}",
        "mapping_status": mapping_status,
        "mapped_native_scan_number": native_scan,
        "mapped_ms_level": ms_level,
        "mapping_rt_sec": mapping_rt,
        "trace_minus_mapping_rt_sec": (
            f"{rt - m['rt_sec']:.6f}" if m else ""
        ),
        "native_spectrum_id": spectrum_id,
        "local_centroid_status": centroid_status,
        "centroid_native_scan_level_values": ";".join(
            f"{scan}:{level}" for scan, level in sorted(cent)
        ),
    })

write_tsv(DETAIL, detail)

# One row per unique trace scan index.
by_scan = {}
for r in detail:
    idx = int(r["scan_index_zero_based"])
    if idx not in by_scan:
        by_scan[idx] = r

counts = Counter(
    (r["mapping_status"], r["local_centroid_status"])
    for r in by_scan.values()
)
summary = [{
    "mapping_status": k[0],
    "local_centroid_status": k[1],
    "unique_trace_scans": v,
} for k, v in sorted(counts.items())]
write_tsv(SUMMARY, summary)

print(f"Feature: {FEATURE}")
print(f"Trace observations: {len(detail)}")
print(f"Unique trace scans: {len(by_scan)}")
print(f"Detail: {DETAIL}")
print(f"Summary: {SUMMARY}")
print("\nCoverage summary:")
for r in summary:
    print(
        f'{r["mapping_status"]} / {r["local_centroid_status"]}: '
        f'{r["unique_trace_scans"]} scans'
    )

print("\nTrace scans with no mapping row:")
for idx, r in sorted(by_scan.items()):
    if r["mapping_status"] == "NO_MAPPING_ROW":
        print(
            f'index={idx}, RT={r["trace_rt_sec"]}, '
            f'traces at scan={sum(int(x["scan_index_zero_based"]) == idx for x in detail)}'
        )
