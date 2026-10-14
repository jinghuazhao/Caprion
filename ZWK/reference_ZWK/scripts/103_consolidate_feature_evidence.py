#!/usr/bin/env python3
"""Consolidate evidence for one provisional MS1 isotope-envelope candidate."""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

REPORTS = Path("reports")
FEATURE = "10275997325121874585"
EXPECTED_SPACING = 1.00335483507 / 2.0
TOLERANCE_MDA = 3.0

TRACK_FILE = REPORTS / f"feature{FEATURE}_centroid_tracks_v1.tsv"
PAIR_FILE = REPORTS / f"feature{FEATURE}_shared_scan_track_spacing_v1.tsv"
CANDIDATE_FILE = REPORTS / f"feature{FEATURE}_envelope_candidate_scores_v1.tsv"
FAILURE_FILE = REPORTS / f"feature{FEATURE}_envelope_spacing_failures_v1.tsv"
TRACK67_FILE = REPORTS / f"feature{FEATURE}_tracks_6_7_scan_inspection_v1.tsv"

OUT_TSV = REPORTS / f"feature{FEATURE}_evidence_summary_v1.tsv"
OUT_MD = REPORTS / f"feature{FEATURE}_evidence_report_v1.md"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_files():
    paths = [
        TRACK_FILE, PAIR_FILE, CANDIDATE_FILE,
        FAILURE_FILE, TRACK67_FILE,
    ]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise SystemExit(
            "Missing required input report(s):\n  " + "\n  ".join(missing)
        )


def main():
    require_files()

    tracks = read_tsv(TRACK_FILE)
    pair_rows = read_tsv(PAIR_FILE)
    candidates = read_tsv(CANDIDATE_FILE)
    failures = read_tsv(FAILURE_FILE)
    track67 = read_tsv(TRACK67_FILE)

    track_stats = defaultdict(list)
    for row in tracks:
        tid = int(row["track_id"])
        track_stats[tid].append(row)

    summary = []

    for tid, rows in sorted(track_stats.items()):
        mzs = [float(row["centroid_mz"]) for row in rows]
        intensities = [float(row["centroid_intensity"]) for row in rows]
        rts = [float(row["retention_time_seconds"]) for row in rows]
        summary.append({
            "record_type": "track",
            "candidate": "",
            "track_or_pair": str(tid),
            "n_observations": len(rows),
            "median_mz": f"{median(mzs):.6f}",
            "mz_range_mDa": f"{(max(mzs) - min(mzs)) * 1000:.3f}",
            "median_intensity": f"{median(intensities):.1f}",
            "rt_min_s": f"{min(rts):.5f}",
            "rt_max_s": f"{max(rts):.5f}",
            "median_abs_residual_mDa": "",
            "within_3_mDa": "",
            "notes": "",
        })

    for row in pair_rows:
        summary.append({
            "record_type": "adjacent_pair",
            "candidate": row.get("candidate", ""),
            "track_or_pair": (
                f'{row.get("lower_track", "")}-{row.get("upper_track", "")}'
            ),
            "n_observations": row.get("n_shared_scans", ""),
            "median_mz": "",
            "mz_range_mDa": "",
            "median_intensity": "",
            "rt_min_s": "",
            "rt_max_s": "",
            "median_abs_residual_mDa": row.get(
                "median_absolute_residual_mDa", ""
            ),
            "within_3_mDa": (
                f'{row.get("within_3_mDa", "")}/'
                f'{row.get("n_shared_scans", "")}'
            ),
            "notes": "",
        })

    for row in candidates:
        summary.append({
            "record_type": "candidate_envelope",
            "candidate": row.get("candidate", ""),
            "track_or_pair": row.get("track_sequence", ""),
            "n_observations": row.get("complete_scans", ""),
            "median_mz": "",
            "mz_range_mDa": "",
            "median_intensity": "",
            "rt_min_s": "",
            "rt_max_s": "",
            "median_abs_residual_mDa": row.get(
                "median_scan_median_abs_residual_mDa", ""
            ),
            "within_3_mDa": (
                f'{row.get("all_spacings_within_3_mDa", "")}/'
                f'{row.get("complete_scans", "")}'
            ),
            "notes": "Candidate-level scan diagnostic",
        })

    write_fields = [
        "record_type", "candidate", "track_or_pair", "n_observations",
        "median_mz", "mz_range_mDa", "median_intensity", "rt_min_s",
        "rt_max_s", "median_abs_residual_mDa", "within_3_mDa", "notes",
    ]
    with OUT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=write_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    # Derive the track 6–7 diagnostic directly from its scan-level output.
    residuals67 = [float(row["residual_mDa"]) for row in track67]
    passed67 = sum(abs(value) <= TOLERANCE_MDA for value in residuals67)

    candidate_lines = []
    for row in candidates:
        candidate_lines.append(
            f'- **{row.get("candidate", "unknown")}** '
            f'({row.get("track_sequence", "tracks unspecified")}): '
            f'{row.get("all_spacings_within_3_mDa", "?")}/'
            f'{row.get("complete_scans", "?")} complete scans pass '
            f'all adjacent spacings within ±{TOLERANCE_MDA:.0f} mDa; '
            f'median scan-level median absolute residual '
            f'{row.get("median_scan_median_abs_residual_mDa", "?")} mDa.'
        )

    md = f"""# Feature-level isotope-envelope evidence

- Feature ID: `{FEATURE}`
- Charge considered: 2
- Theoretical isotope spacing: `{EXPECTED_SPACING:.12f} Da`
- Spacing tolerance used in candidate diagnostics: ±{TOLERANCE_MDA:.1f} mDa
- Status: **provisional; common-precursor identity not established**

## Evidence summary

### Candidate-envelope results

{chr(10).join(candidate_lines)}

### Track 6–7 transition

- Shared scans: {len(track67)}
- Scans within ±{TOLERANCE_MDA:.0f} mDa: {passed67}/{len(track67)}
- Median absolute residual: {median(abs(value) for value in residuals67):.3f} mDa
- Residual range: {min(residuals67):+.3f} to {max(residuals67):+.3f} mDa

The track 6–7 transition is poorly supported by the spacing criterion. Track 7 should not be included in the provisional envelope on the current evidence.

## Provisional interpretation

The four-track sequence 2–5 is the conservative candidate core. Track 6 is a possible extension, but its transition from track 5 has weaker spacing support. Track 7 is not supported as a consistent isotope successor to track 6.

These are diagnostics of centroid spacing and track recurrence, not calibrated probabilities that peaks belong to the same precursor. The track data are derived from the same acquisitions, so the observations are not independent validations. Intensity ordering and centroid-to-trace agreement do not by themselves prove isotope identity.

## Independence and limitations

- This report does not use the original XLSX isotope-group assignments to construct the candidate.
- No definitive independent isotope-group ID is assigned here.
- No peptide-to-isotope-group link is established by this report.
- Independent feature detection, isotope grouping and subsequent comparison with the original workbook remain separate stages.
- The current report does not establish a combined-panel false discovery rate.

## Input reports

- `{TRACK_FILE}`
- `{PAIR_FILE}`
- `{CANDIDATE_FILE}`
- `{FAILURE_FILE}`
- `{TRACK67_FILE}`
"""
    OUT_MD.write_text(md)

    print(f"Evidence summary TSV: {OUT_TSV}")
    print(f"Evidence report Markdown: {OUT_MD}")
    print(f"Track rows: {len(track_stats)}")
    print(f"Pair rows: {len(pair_rows)}")
    print(f"Candidate rows: {len(candidates)}")
    print(
        f"Track 6–7 spacing: {passed67}/{len(track67)} within "
        f"±{TOLERANCE_MDA:.0f} mDa"
    )
    print("Provisional interpretation: tracks 2–5 are the conservative core.")


def write_fields_unused():
    # Kept intentionally empty; the output schema is defined in main().
    return None


if __name__ == "__main__":
    main()
