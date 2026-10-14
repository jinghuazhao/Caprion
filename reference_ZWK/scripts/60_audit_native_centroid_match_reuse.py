#!/usr/bin/env python3
"""Audit duplicate centroid reuse and isotope-spacing residuals from Script 59."""
from pathlib import Path
import csv
import math
import statistics
import sys
from collections import defaultdict

ROOT = Path.cwd()
REPORTS = ROOT / "reports"

TRACK_IN = REPORTS / "native_ms1_centroid_tracks_charge2_v2.tsv"
GAP_IN = REPORTS / "native_ms1_centroid_track_adjacent_gaps_charge2_v2.tsv"

REUSE_OUT = REPORTS / "native_ms1_centroid_duplicate_reuse_audit_v2.tsv"
PAIR_OUT = REPORTS / "native_ms1_centroid_spacing_distinct_peaks_v2.tsv"
SUMMARY_OUT = REPORTS / "native_ms1_centroid_spacing_distinct_peak_summary_v2.tsv"

TRACE_IDS = [1, 2, 3, 4, 5, 6]
TOL_MDA = 2.0


def read_tsv(path):
    if not path.is_file():
        sys.exit(f"ERROR: missing input file: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    if not rows:
        print(f"WARNING: no rows to write: {path}")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0].keys()), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def num(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def fmt(value, digits=3):
    return "" if value is None else f"{value:.{digits}f}"


def median(values):
    return statistics.median(values) if values else None


def main():
    track_rows = read_tsv(TRACK_IN)
    gap_rows = read_tsv(GAP_IN)

    # Key observations by scan and target.
    by_scan = defaultdict(dict)
    for row in track_rows:
        if row["matched"] != "1":
            continue
        scan = int(row["scan_index"])
        trace = int(row["trace_target"])
        mz = num(row["matched_mz"])
        if mz is None:
            continue
        by_scan[scan][trace] = {
            "mz": mz,
            "intensity": num(row["intensity"]),
            "error_mDa": num(row["mz_error_mDa"]),
            "rt_sec": num(row["rt_sec"]),
        }

    reuse_rows = []
    for scan in sorted(by_scan):
        matches = by_scan[scan]
        # Native centroid m/z values are printed to 6 decimals by Script 59.
        # Exact equality here identifies reused selected centroids.
        mz_to_traces = defaultdict(list)
        for trace, obs in matches.items():
            mz_to_traces[round(obs["mz"], 6)].append(trace)

        for mz, traces in sorted(mz_to_traces.items()):
            if len(traces) > 1:
                reuse_rows.append({
                    "scan_index": scan,
                    "rt_sec": fmt(next(iter(matches.values()))["rt_sec"], 5),
                    "reused_mz": f"{mz:.6f}",
                    "trace_targets_reusing_peak": ",".join(map(str, sorted(traces))),
                    "number_of_targets": len(traces),
                    "interpretation": "same selected centroid reused by multiple targets",
                })

    # Recompute adjacent gaps directly from Script 59's matched centroid
    # assignments. A pair is distinct if the two selected centroid m/z values
    # differ. This does not resolve all possible assignment ambiguities.
    pair_rows = []
    pair_summary = defaultdict(list)

    for left, right in zip(TRACE_IDS[:-1], TRACE_IDS[1:]):
        for scan in sorted(by_scan):
            matches = by_scan[scan]
            if left not in matches or right not in matches:
                continue

            a = matches[left]
            b = matches[right]
            gap_da = b["mz"] - a["mz"]
            residual_mDa = (gap_da - 1.00335483507 / 2.0) * 1000.0
            distinct = round(a["mz"], 6) != round(b["mz"], 6)
            within = abs(residual_mDa) <= TOL_MDA

            row = {
                "scan_index": scan,
                "rt_sec": fmt(a["rt_sec"], 5),
                "left_trace": left,
                "right_trace": right,
                "left_mz": fmt(a["mz"], 6),
                "right_mz": fmt(b["mz"], 6),
                "left_intensity": fmt(a["intensity"], 3),
                "right_intensity": fmt(b["intensity"], 3),
                "observed_gap_da": fmt(gap_da, 9),
                "spacing_residual_mDa": fmt(residual_mDa, 3),
                "distinct_selected_centroids": int(distinct),
                "abs_residual_le_2mDa": int(within),
            }
            pair_rows.append(row)

            pair_summary[(left, right)].append({
                "residual": residual_mDa,
                "distinct": distinct,
                "within": within,
                "scan": scan,
                "rt": a["rt_sec"],
            })

    summary_rows = []
    print("Duplicate centroid reuse")
    print(f"  Scans with any matched targets: {len(by_scan)}")
    print(f"  Reuse events: {len(reuse_rows)}")
    if reuse_rows:
        for row in reuse_rows:
            print(
                f"  scan {row['scan_index']} RT {row['rt_sec']}: "
                f"m/z {row['reused_mz']} reused by traces "
                f"{row['trace_targets_reusing_peak']}"
            )
    else:
        print("  No exact centroid reuse detected.")

    print("\nAdjacent pair summary (all matched pairs vs distinct-centroid pairs)")
    for left, right in zip(TRACE_IDS[:-1], TRACE_IDS[1:]):
        values = pair_summary[(left, right)]
        distinct = [v for v in values if v["distinct"]]
        residuals = [v["residual"] for v in distinct]
        abs_res = [abs(x) for x in residuals]
        n_within = sum(x <= TOL_MDA for x in abs_res)
        frac = n_within / len(abs_res) if abs_res else None

        summary_rows.append({
            "left_trace": left,
            "right_trace": right,
            "n_matched_pairs": len(values),
            "n_distinct_centroid_pairs": len(distinct),
            "n_distinct_pairs_with_abs_residual_le_2mDa": n_within,
            "fraction_distinct_pairs_within_2mDa": fmt(frac, 4),
            "median_signed_spacing_residual_mDa": fmt(median(residuals), 3),
            "median_abs_spacing_residual_mDa": fmt(median(abs_res), 3),
            "max_abs_spacing_residual_mDa": fmt(max(abs_res) if abs_res else None, 3),
            "distinct_pair_scans": ",".join(str(v["scan"]) for v in distinct),
        })
        print(
            f"  {left}-{right}: matched={len(values)}, distinct={len(distinct)}, "
            f"median |residual|={fmt(median(abs_res), 3)} mDa, "
            f"within 2 mDa={n_within}/{len(distinct)} "
            f"({fmt(frac, 3) if frac is not None else 'NA'})"
        )

    write_tsv(REUSE_OUT, reuse_rows)
    write_tsv(PAIR_OUT, pair_rows)
    write_tsv(SUMMARY_OUT, summary_rows)

    print("\nFiles written:")
    print(f"  {REUSE_OUT}")
    print(f"  {PAIR_OUT}")
    print(f"  {SUMMARY_OUT}")
    print(
        "\nCaution: different m/z values do not guarantee correct isotope "
        "assignment; distinct-centroid filtering only rules out exact reuse "
        "of the same selected peak."
    )


if __name__ == "__main__":
    main()
