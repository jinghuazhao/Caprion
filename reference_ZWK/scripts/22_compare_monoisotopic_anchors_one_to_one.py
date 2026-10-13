#!/usr/bin/env python3
"""Reassess isotope anchors using one-to-one observed-peak matching.

This is an exploratory diagnostic, not a validated monoisotopic assignment.
Original script 07 and its output files are left unchanged.
"""

from pathlib import Path
from collections import defaultdict
from bisect import bisect_left, bisect_right

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

BASE = Path("openms_pilot")
CANDIDATES = BASE / "psm_feature_offset_candidates_both.tsv"

TRACE_FILES = {
    "first_pilot": BASE / "isotope_trace_peaks.tsv",
    "second_acquisition": BASE / "isotope_trace_peaks_second.tsv",
}

OUT_SCAN = BASE / "anchor_one_to_one_v2_per_scan.tsv"
OUT_SUMMARY = BASE / "anchor_one_to_one_v2_summary.tsv"
OUT_BEST = BASE / "anchor_one_to_one_v2_best.tsv"

ISOTOPE_SPACING = 1.00335483507
PPM_TOL = 10.0
RT_WINDOW = 2.0
ANCHOR_OFFSETS = range(-2, 3)
MAX_ISOTOPE_PEAK = 4

# Assignment objective: maximise the number of valid matches first,
# then minimise total absolute ppm error. The cardinality bonus exceeds
# the maximum possible difference in mass-error penalties.
CARDINALITY_BONUS = (MAX_ISOTOPE_PEAK + 2) * (PPM_TOL + 1)
UNMATCHED_COST = 0.0
INVALID_COST = 1e6


def match_isotopes(mzs, intensities, anchor_mz, isotope_step):
    """Match M0-M4 to distinct observed peaks within the ppm tolerance."""
    n_theoretical = MAX_ISOTOPE_PEAK + 1

    theoretical = anchor_mz + np.arange(n_theoretical) * isotope_step
    observed = np.asarray(mzs, dtype=float)
    intensity = np.asarray(intensities, dtype=float)

    valid = (
        np.isfinite(observed)
        & np.isfinite(intensity)
        & (observed > 0)
        & (intensity > 0)
    )
    observed = observed[valid]
    intensity = intensity[valid]

    if len(observed):
        order = np.argsort(observed)
        observed = observed[order]
        intensity = intensity[order]

    n_observed = len(observed)

    # Columns 0..n_observed-1 are observed peaks; the final five
    # columns are dummy assignments representing unmatched positions.
    cost = np.full(
        (n_theoretical, n_observed + n_theoretical),
        INVALID_COST,
        dtype=float,
    )

    ppm_errors = np.full((n_theoretical, n_observed), np.nan)

    for k, target in enumerate(theoretical):
        if not np.isfinite(target) or target <= 0:
            continue

        errors = (observed - target) / target * 1e6
        ppm_errors[k] = errors

        within = np.abs(errors) <= PPM_TOL
        # Every valid match receives a cardinality bonus. This makes
        # an additional valid match preferable to any lower-error
        # assignment with fewer matches.
        valid_indices = np.flatnonzero(within)
        cost[k, valid_indices] = (
            -CARDINALITY_BONUS
            + np.abs(errors[valid_indices])
            + valid_indices * 1e-10
        )

        cost[k, n_observed:] = UNMATCHED_COST

    row_indices, col_indices = linear_sum_assignment(cost)

    matched = {}
    for k, j in zip(row_indices, col_indices):
        if j >= n_observed:
            continue
        if cost[k, j] >= UNMATCHED_COST:
            continue

        matched[int(k)] = {
            "mz": float(observed[j]),
            "intensity": float(intensity[j]),
            "ppm": float(ppm_errors[k, j]),
        }

    positions = sorted(matched)
    n_found = len(positions)

    consecutive = 0
    for k in range(n_theoretical):
        if k not in matched:
            break
        consecutive += 1

    spacing_errors = []
    for k in range(MAX_ISOTOPE_PEAK):
        if k in matched and (k + 1) in matched:
            gap = matched[k + 1]["mz"] - matched[k]["mz"]
            spacing_errors.append(
                (gap - isotope_step) / isotope_step * 1e6
            )

    median_abs_spacing_error = (
        float(np.median(np.abs(spacing_errors)))
        if spacing_errors else np.nan
    )

    median_abs_mass_error = (
        float(np.median([
            abs(matched[k]["ppm"]) for k in positions
        ]))
        if positions else np.nan
    )

    # Diagnostic components are based solely on one-to-one matches.
    coverage = n_found / n_theoretical
    continuity = consecutive / n_theoretical

    mass_quality = (
        float(np.exp(-0.5 * (median_abs_mass_error / 5.0) ** 2))
        if np.isfinite(median_abs_mass_error) else 0.0
    )

    # Spacing error remains a diagnostic; it is not included in the score
    # until a suitable tolerance is calibrated against centroid precision.
    score = (
        0.60 * coverage
        + 0.25 * continuity
        + 0.15 * mass_quality
    )

    result = {
        "isotope_positions_found": n_found,
        "isotope_positions": ",".join(map(str, positions)),
        "consecutive_isotopes_from_anchor": consecutive,
        "median_abs_adjacent_spacing_error_ppm":
            median_abs_spacing_error,
        "median_abs_individual_peak_error_ppm": median_abs_mass_error,
        "anchor_score_one_to_one_v2": score,
        "anchor_peak_intensity": (
            matched[0]["intensity"] if 0 in matched else np.nan
        ),
        "n_observed_peaks_available": n_observed,
    }

    for k in range(n_theoretical):
        peak = matched.get(k)
        result[f"M{k}_observed_mz"] = (
            peak["mz"] if peak else np.nan
        )
        result[f"M{k}_ppm_error"] = (
            peak["ppm"] if peak else np.nan
        )
        result[f"M{k}_intensity"] = (
            peak["intensity"] if peak else np.nan
        )

    return result


def main():
    candidates = pd.read_csv(CANDIDATES, sep="\t", low_memory=False)

    required = {
        "acquisition", "Spectrum", "feature_id", "charge",
        "isotope_offset", "psm_rt_sec", "psm_mz_used",
    }
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table missing columns: {sorted(missing)}")

    candidates["feature_id"] = candidates["feature_id"].astype(str)

    for col in ["charge", "isotope_offset", "psm_rt_sec", "psm_mz_used"]:
        candidates[col] = pd.to_numeric(candidates[col], errors="coerce")

    candidates = candidates.dropna(
        subset=["charge", "isotope_offset", "psm_rt_sec", "psm_mz_used"]
    ).copy()
    candidates = candidates[candidates["charge"] > 0].copy()
    candidates["charge"] = candidates["charge"].astype(int)
    candidates["isotope_offset"] = candidates["isotope_offset"].astype(int)

    candidate_keys = [
        "acquisition", "Spectrum", "feature_id", "charge",
        "isotope_offset",
    ]
    candidates = candidates.drop_duplicates(candidate_keys).reset_index(drop=True)

    scan_rows = []

    for acquisition, trace_path in TRACE_FILES.items():
        if not trace_path.exists():
            raise FileNotFoundError(trace_path)

        print(f"Loading {trace_path}", flush=True)
        traces = pd.read_csv(trace_path, sep="\t", low_memory=False)

        required_trace = {
            "feature_id", "scan_index", "rt_sec",
            "observed_mz", "peak_intensity",
        }
        missing = required_trace - set(traces.columns)
        if missing:
            raise ValueError(f"{trace_path} missing columns: {sorted(missing)}")

        traces["feature_id"] = traces["feature_id"].astype(str)
        for col in ["scan_index", "rt_sec", "observed_mz", "peak_intensity"]:
            traces[col] = pd.to_numeric(traces[col], errors="coerce")

        traces = traces.dropna(
            subset=[
                "feature_id", "scan_index", "rt_sec",
                "observed_mz", "peak_intensity",
            ]
        )
        traces = traces[
            (traces["peak_intensity"] > 0)
            & (traces["observed_mz"] > 0)
        ].copy()

        by_scan = {}
        for (fid, scan), group in traces.groupby(
            ["feature_id", "scan_index"], sort=False
        ):
            group = group.sort_values("observed_mz")
            by_scan[(fid, int(scan))] = (
                float(group["rt_sec"].median()),
                group["observed_mz"].to_numpy(float),
                group["peak_intensity"].to_numpy(float),
            )

        scans_by_feature = defaultdict(list)
        for (fid, scan), (rt, mzs, intensities) in by_scan.items():
            scans_by_feature[fid].append((rt, scan, mzs, intensities))

        for fid in scans_by_feature:
            scans_by_feature[fid].sort(key=lambda x: x[0])

        rts_by_feature = {
            fid: [item[0] for item in items]
            for fid, items in scans_by_feature.items()
        }

        subset = candidates[candidates["acquisition"] == acquisition]
        print(
            f"{acquisition}: {len(subset):,} candidates; "
            f"{len(by_scan):,} feature-scan groups",
            flush=True,
        )

        for i, row in enumerate(subset.itertuples(index=False), start=1):
            fid = str(row.feature_id)
            charge = int(row.charge)
            step = ISOTOPE_SPACING / charge
            psm_mz = float(row.psm_mz_used)
            psm_rt = float(row.psm_rt_sec)

            feature_scans = scans_by_feature.get(fid, [])
            feature_rts = rts_by_feature.get(fid, [])

            left = bisect_left(feature_rts, psm_rt - RT_WINDOW)
            right = bisect_right(feature_rts, psm_rt + RT_WINDOW)

            for rt, scan, mzs, intensities in feature_scans[left:right]:
                for offset in ANCHOR_OFFSETS:
                    anchor_mz = psm_mz - offset * step
                    matched = match_isotopes(
                        mzs, intensities, anchor_mz, step
                    )

                    scan_rows.append({
                        "acquisition": acquisition,
                        "Spectrum": row.Spectrum,
                        "feature_id": fid,
                        "charge": charge,
                        "nominated_isotope_offset": int(row.isotope_offset),
                        "scan_index": scan,
                        "rt_sec": rt,
                        "psm_rt_sec": psm_rt,
                        "psm_mz": psm_mz,
                        "anchor_offset": offset,
                        "anchor_mz": anchor_mz,
                        **matched,
                    })

            if i % 1000 == 0:
                print(
                    f"  {acquisition}: processed {i:,}/{len(subset):,}",
                    flush=True,
                )

        del traces, by_scan, scans_by_feature, rts_by_feature

    per_scan = pd.DataFrame(scan_rows)
    per_scan.to_csv(OUT_SCAN, sep="\t", index=False)
    print(f"Per-scan results: {OUT_SCAN}", flush=True)

    group_cols = [
        "acquisition", "Spectrum", "feature_id", "charge",
        "nominated_isotope_offset", "anchor_offset",
    ]

    summary = per_scan.groupby(group_cols, dropna=False).agg(
        scans_examined=("anchor_score_one_to_one_v2", "size"),
        mean_anchor_score=("anchor_score_one_to_one_v2", "mean"),
        median_anchor_score=("anchor_score_one_to_one_v2", "median"),
        median_isotope_positions_found=("isotope_positions_found", "median"),
        median_consecutive_isotopes=(
            "consecutive_isotopes_from_anchor", "median"
        ),
        scans_with_3_consecutive_isotopes=(
            "consecutive_isotopes_from_anchor",
            lambda x: int((x >= 3).sum()),
        ),
        median_abs_spacing_error_ppm=(
            "median_abs_adjacent_spacing_error_ppm", "median"
        ),
        median_abs_peak_error_ppm=(
            "median_abs_individual_peak_error_ppm", "median"
        ),
    ).reset_index()

    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    summary["anchor_rank"] = summary.groupby(
        [
            "acquisition", "Spectrum", "feature_id", "charge",
            "nominated_isotope_offset",
        ]
    )["mean_anchor_score"].rank(method="min", ascending=False)

    best = summary[summary["anchor_rank"] == 1].copy()
    best = best.rename(columns={
        "anchor_offset": "best_anchor_offset",
        "mean_anchor_score": "best_anchor_score",
    })

    best.to_csv(OUT_BEST, sep="\t", index=False)

    print(f"Summary: {OUT_SUMMARY}")
    print(f"Best offsets (ties retained): {OUT_BEST}")
    print(f"Input candidates: {len(candidates):,}")
    print(f"Per-scan rows: {len(per_scan):,}")


if __name__ == "__main__":
    main()
