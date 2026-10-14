#!/usr/bin/env python3
"""Validate candidate isotope spacing against native MS1 centroid observations."""

from pathlib import Path
import csv
import statistics
from collections import defaultdict

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
FEATURE_ID = "10275997325121874585"
THEORETICAL_SPACING_DA = 1.00335483507 / 2
TOLERANCE_MDA = 3.0

TRACE_FILE = REPORTS / "full_trace_chromatography_v1.tsv"
CENTROID_FILE = REPORTS / "local_ms1_centroids_feature10275997325121874585_v1.tsv"
MAPPING_FILE = REPORTS / "native_spectrum_index_mapping_v2.tsv"

DETAIL_OUT = REPORTS / "native_centroid_spacing_validation_v1.tsv"
SUMMARY_OUT = REPORTS / "native_centroid_spacing_validation_summary_v1.tsv"


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


for path in (TRACE_FILE, CENTROID_FILE, MAPPING_FILE):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing input: {path}\n"
            "Check the existing report filenames before rerunning."
        )

trace_rows = [
    r for r in read_tsv(TRACE_FILE)
    if r["feature_id"] == FEATURE_ID
]
centroid_rows = read_tsv(CENTROID_FILE)
mapping_rows = read_tsv(MAPPING_FILE)

# The trace table uses zero-based mzML spectrum indices. Map them explicitly
# to native scan identifiers rather than assuming they are interchangeable.
mapping = {}
for r in mapping_rows:
    idx_key = next(
        (k for k in r if "index" in k.lower() and "spectrum" in k.lower()),
        None,
    )
    scan_key = next(
        (k for k in r if "scan" in k.lower() and "number" in k.lower()),
        None,
    )
    if idx_key and scan_key and r[idx_key] and r[scan_key]:
        try:
            mapping[int(r[idx_key])] = r[scan_key]
        except ValueError:
            pass

traces_by_scan = defaultdict(dict)
for r in trace_rows:
    scan = int(r["scan_index_zero_based"])
    trace = int(r["trace_index"])
    traces_by_scan[scan][trace] = {
        "mz": float(r["observed_mz"]),
        "intensity": float(r["peak_intensity"]),
        "rt": float(r["rt_sec"]),
    }

# Discover centroid columns by name to accommodate the existing export schema.
if not centroid_rows:
    raise ValueError(f"No rows found in {CENTROID_FILE}")

fields = centroid_rows[0].keys()


def find_col(tokens):
    for col in fields:
        low = col.lower()
        if all(token in low for token in tokens):
            return col
    return None


scan_col = find_col(["scan"])
mz_col = find_col(["mz"])
intensity_col = find_col(["intensity"])

if not scan_col or not mz_col or not intensity_col:
    raise ValueError(
        "Could not identify centroid columns. "
        f"Available columns: {list(fields)}"
    )

centroids_by_scan = defaultdict(list)
for r in centroid_rows:
    try:
        scan = int(float(r[scan_col]))
        mz = float(r[mz_col])
        intensity = float(r[intensity_col])
    except (ValueError, TypeError):
        continue
    centroids_by_scan[scan].append({"mz": mz, "intensity": intensity})

detail = []

for ordinal in sorted(traces_by_scan):
    traces = traces_by_scan[ordinal]
    native_scan = mapping.get(ordinal, "")
    try:
        native_scan_int = int(str(native_scan).split("=")[-1])
    except (ValueError, TypeError):
        native_scan_int = None

    local = (
        centroids_by_scan.get(native_scan_int, [])
        if native_scan_int is not None else []
    )

    # Match each trace m/z to native centroid peaks within 0.015 Da.
    matched = {}
    for trace, obs in traces.items():
        candidates = [
            c for c in local if abs(c["mz"] - obs["mz"]) <= 0.015
        ]
        candidates.sort(key=lambda c: abs(c["mz"] - obs["mz"]))
        matched[trace] = {
            "candidate_count": len(candidates),
            "matched_mz": candidates[0]["mz"] if candidates else None,
            "matched_intensity": candidates[0]["intensity"] if candidates else None,
            "delta_mDa": (
                (candidates[0]["mz"] - obs["mz"]) * 1000
                if candidates else None
            ),
        }

    for lower in range(1, 6):
        upper = lower + 1
        if lower not in traces or upper not in traces:
            continue

        lo = traces[lower]
        hi = traces[upper]
        observed_spacing = hi["mz"] - lo["mz"]
        residual_mDa = (observed_spacing - THEORETICAL_SPACING_DA) * 1000

        lo_match = matched[lower]
        hi_match = matched[upper]

        native_spacing = ""
        native_residual = ""
        if lo_match["matched_mz"] is not None and hi_match["matched_mz"] is not None:
            native_spacing = hi_match["matched_mz"] - lo_match["matched_mz"]
            native_residual = (
                native_spacing - THEORETICAL_SPACING_DA
            ) * 1000

        detail.append({
            "feature_id": FEATURE_ID,
            "scan_index_zero_based": ordinal,
            "native_scan_id": native_scan,
            "rt_sec": f'{lo["rt"]:.6f}',
            "lower_trace_index": lower,
            "upper_trace_index": upper,
            "trace_observed_spacing_Da": f"{observed_spacing:.9f}",
            "trace_spacing_residual_mDa": f"{residual_mDa:.6f}",
            "trace_spacing_within_3mDa": abs(residual_mDa) <= TOLERANCE_MDA,
            "lower_centroid_candidates": lo_match["candidate_count"],
            "upper_centroid_candidates": hi_match["candidate_count"],
            "lower_trace_centroid_delta_mDa": (
                f'{lo_match["delta_mDa"]:.6f}'
                if lo_match["delta_mDa"] is not None else ""
            ),
            "upper_trace_centroid_delta_mDa": (
                f'{hi_match["delta_mDa"]:.6f}'
                if hi_match["delta_mDa"] is not None else ""
            ),
            "native_centroid_spacing_Da": (
                f"{native_spacing:.9f}" if native_spacing != "" else ""
            ),
            "native_spacing_residual_mDa": (
                f"{native_residual:.6f}" if native_residual != "" else ""
            ),
            "native_spacing_within_3mDa": (
                abs(native_residual) <= TOLERANCE_MDA
                if native_residual != "" else ""
            ),
        })

with DETAIL_OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(detail[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(detail)

summary = []
for lower in range(1, 6):
    upper = lower + 1
    subset = [
        r for r in detail
        if int(r["lower_trace_index"]) == lower
        and int(r["upper_trace_index"]) == upper
    ]
    trace_residuals = [float(r["trace_spacing_residual_mDa"]) for r in subset]
    native_residuals = [
        float(r["native_spacing_residual_mDa"])
        for r in subset if r["native_spacing_residual_mDa"] != ""
    ]
    trace_pass = sum(
        str(r["trace_spacing_within_3mDa"]).lower() == "true"
        for r in subset
    )
    native_pass = sum(
        str(r["native_spacing_within_3mDa"]).lower() == "true"
        for r in subset
        if r["native_spacing_within_3mDa"] != ""
    )
    unique_pair = sum(
        int(r["lower_centroid_candidates"]) == 1
        and int(r["upper_centroid_candidates"]) == 1
        for r in subset
    )

    summary.append({
        "lower_trace_index": lower,
        "upper_trace_index": upper,
        "n_matched_scans": len(subset),
        "trace_spacing_median_abs_residual_mDa": (
            f"{statistics.median(abs(x) for x in trace_residuals):.6f}"
            if trace_residuals else ""
        ),
        "trace_spacing_within_3mDa": f"{trace_pass}/{len(subset)}",
        "native_centroid_n": len(native_residuals),
        "native_spacing_median_abs_residual_mDa": (
            f"{statistics.median(abs(x) for x in native_residuals):.6f}"
            if native_residuals else ""
        ),
        "native_spacing_within_3mDa": (
            f"{native_pass}/{len(native_residuals)}"
            if native_residuals else ""
        ),
        "unique_centroid_match_pairs": f"{unique_pair}/{len(subset)}",
    })

with SUMMARY_OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(summary)

print(f"Feature: {FEATURE_ID}")
print(f"Native centroid detail: {DETAIL_OUT}")
print(f"Summary: {SUMMARY_OUT}")
print("\nPair summary:")
for r in summary:
    print(
        f'{r["lower_trace_index"]}-{r["upper_trace_index"]}: '
        f'trace spacing ±3 mDa={r["trace_spacing_within_3mDa"]}; '
        f'native centroid spacing ±3 mDa={r["native_spacing_within_3mDa"]}; '
        f'unique centroid match pairs={r["unique_centroid_match_pairs"]}'
    )

print(
    "\nCaution: matching trace m/z to centroids verifies correspondence to "
    "the exported peaks; it does not independently prove isotope identity."
)
