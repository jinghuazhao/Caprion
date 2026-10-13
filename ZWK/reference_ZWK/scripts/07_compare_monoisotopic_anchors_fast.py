from pathlib import Path
from collections import defaultdict
from bisect import bisect_left, bisect_right
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
CANDIDATES = BASE / "psm_feature_offset_candidates_both.tsv"

TRACE_FILES = {
    "first_pilot": BASE / "isotope_trace_peaks.tsv",
    "second_acquisition": BASE / "isotope_trace_peaks_second.tsv",
}

SCAN_OUT = BASE / "monoisotopic_anchor_comparison_per_scan_fast.tsv"
SUMMARY_OUT = BASE / "monoisotopic_anchor_comparison_summary_fast.tsv"
BEST_OUT = BASE / "monoisotopic_anchor_best_per_candidate_fast.tsv"

ISOTOPE_SPACING = 1.00335483507
PPM_TOL = 10.0
RT_WINDOW = 2.0
ANCHOR_OFFSETS = range(-2, 3)
MAX_ISOTOPE_PEAK = 4
BATCH_SIZE = 50_000


def nearest_peak(mzs, intensities, theoretical_mz):
    """Return the nearest peak within the ppm tolerance, if present."""
    if theoretical_mz <= 0 or not np.isfinite(theoretical_mz):
        return None

    lo = theoretical_mz * (1 - PPM_TOL * 1e-6)
    hi = theoretical_mz * (1 + PPM_TOL * 1e-6)

    left = np.searchsorted(mzs, lo, side="left")
    right = np.searchsorted(mzs, hi, side="right")

    if left >= right:
        return None

    indices = np.arange(left, right)
    i = indices[np.argmin(np.abs(mzs[indices] - theoretical_mz))]
    ppm = (mzs[i] - theoretical_mz) / theoretical_mz * 1e6

    return {
        "mz": float(mzs[i]),
        "intensity": float(intensities[i]),
        "ppm": float(ppm),
    }


def score_anchor(mzs, intensities, anchor_mz, isotope_step):
    peaks = [
        nearest_peak(mzs, intensities, anchor_mz + k * isotope_step)
        for k in range(MAX_ISOTOPE_PEAK + 1)
    ]

    present = [k for k, peak in enumerate(peaks) if peak is not None]

    consecutive = 0
    for peak in peaks:
        if peak is None:
            break
        consecutive += 1

    errors = []
    for k in range(MAX_ISOTOPE_PEAK):
        p1, p2 = peaks[k], peaks[k + 1]
        if p1 is not None and p2 is not None:
            observed = p2["mz"] - p1["mz"]
            errors.append(
                (observed - isotope_step) / isotope_step * 1e6
            )

    median_abs_error = (
        float(np.median(np.abs(errors))) if errors else np.nan
    )

    coverage = len(present) / (MAX_ISOTOPE_PEAK + 1)
    continuity = consecutive / (MAX_ISOTOPE_PEAK + 1)
    spacing_quality = (
        max(0.0, 1.0 - median_abs_error / 50.0)
        if np.isfinite(median_abs_error) else 0.0
    )

    # Ranking heuristic, not a probability.
    score = 0.55 * coverage + 0.30 * continuity + 0.15 * spacing_quality

    return {
        "isotope_positions_found": len(present),
        "isotope_positions": ",".join(map(str, present)),
        "consecutive_isotopes_from_anchor": consecutive,
        "median_abs_adjacent_spacing_error_ppm": median_abs_error,
        "anchor_score": score,
        "anchor_peak_intensity": (
            peaks[0]["intensity"] if peaks[0] else np.nan
        ),
        "M1_intensity": peaks[1]["intensity"] if peaks[1] else np.nan,
        "M2_intensity": peaks[2]["intensity"] if peaks[2] else np.nan,
        "M3_intensity": peaks[3]["intensity"] if peaks[3] else np.nan,
        "M4_intensity": peaks[4]["intensity"] if peaks[4] else np.nan,
    }


def main():
    cand = pd.read_csv(CANDIDATES, sep="\t", low_memory=False)

    required_cand = {
        "acquisition", "Spectrum", "feature_id", "charge",
        "isotope_offset", "psm_rt_sec", "psm_mz_used"
    }
    missing = required_cand - set(cand.columns)
    if missing:
        raise ValueError(f"Candidate table missing columns: {sorted(missing)}")

    cand["feature_id"] = cand["feature_id"].astype(str)

    for col in ["charge", "isotope_offset", "psm_rt_sec", "psm_mz_used"]:
        cand[col] = pd.to_numeric(cand[col], errors="coerce")

    cand = cand.dropna(
        subset=["charge", "isotope_offset", "psm_rt_sec", "psm_mz_used"]
    ).copy()

    cand = cand[cand["charge"] > 0].copy()
    cand["charge"] = cand["charge"].astype(int)
    cand["isotope_offset"] = cand["isotope_offset"].astype(int)

    candidate_keys = [
        "acquisition", "Spectrum", "feature_id", "charge", "isotope_offset"
    ]
    cand = cand.drop_duplicates(candidate_keys).reset_index(drop=True)

    scan_summary_parts = []
    best_parts = []
    total_candidates = len(cand)

    with open(SCAN_OUT, "w") as scan_handle:
        header_written = False

        for acquisition, trace_path in TRACE_FILES.items():
            if not trace_path.exists():
                raise FileNotFoundError(trace_path)

            print(f"\nLoading {trace_path}", flush=True)

            traces = pd.read_csv(trace_path, sep="\t", low_memory=False)
            required_trace = {
                "feature_id", "scan_index", "rt_sec",
                "observed_mz", "peak_intensity"
            }
            missing = required_trace - set(traces.columns)
            if missing:
                raise ValueError(
                    f"{trace_path} missing columns: {sorted(missing)}"
                )

            traces["feature_id"] = traces["feature_id"].astype(str)
            for col in [
                "scan_index", "rt_sec", "observed_mz", "peak_intensity"
            ]:
                traces[col] = pd.to_numeric(traces[col], errors="coerce")

            traces = traces.dropna(
                subset=[
                    "feature_id", "scan_index", "rt_sec",
                    "observed_mz", "peak_intensity"
                ]
            )
            traces = traces[traces["peak_intensity"] > 0]

            # One sorted centroid array per feature and scan.
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

            # Index scans by feature. Retention times are sorted for binary search.
            scans_by_feature = defaultdict(list)
            for (fid, scan), (rt, mzs, intensities) in by_scan.items():
                scans_by_feature[fid].append(
                    (rt, scan, mzs, intensities)
                )

            for fid in scans_by_feature:
                scans_by_feature[fid].sort(key=lambda x: x[0])

            rts_by_feature = {
                fid: [item[0] for item in items]
                for fid, items in scans_by_feature.items()
            }

            subset = cand[cand["acquisition"] == acquisition]
            print(
                f"{acquisition}: {len(subset):,} candidates, "
                f"{len(by_scan):,} feature-scan groups",
                flush=True,
            )

            batch = []
            processed = 0

            for row in subset.itertuples(index=False):
                fid = str(row.feature_id)
                z = int(row.charge)
                psm_mz = float(row.psm_mz_used)
                rt0 = float(row.psm_rt_sec)
                step = ISOTOPE_SPACING / z

                feature_scans = scans_by_feature.get(fid, [])
                feature_rts = rts_by_feature.get(fid, [])

                # Binary search retrieves only scans within RT_WINDOW.
                left = bisect_left(feature_rts, rt0 - RT_WINDOW)
                right = bisect_right(feature_rts, rt0 + RT_WINDOW)

                for rt, scan, mzs, intensities in feature_scans[left:right]:
                    for anchor_offset in ANCHOR_OFFSETS:
                        # anchor_offset > 0 means anchor below PSM m/z.
                        anchor_mz = psm_mz - anchor_offset * step
                        scored = score_anchor(mzs, intensities, anchor_mz, step)

                        batch.append({
                            "acquisition": acquisition,
                            "Spectrum": row.Spectrum,
                            "feature_id": fid,
                            "charge": z,
                            "nominated_isotope_offset": int(row.isotope_offset),
                            "scan_index": scan,
                            "rt_sec": rt,
                            "psm_rt_sec": rt0,
                            "psm_mz": psm_mz,
                            "anchor_offset": anchor_offset,
                            "anchor_mz": anchor_mz,
                            **scored,
                        })

                        if len(batch) >= BATCH_SIZE:
                            pd.DataFrame(batch).to_csv(
                                scan_handle, sep="\t", index=False,
                                header=not header_written
                            )
                            header_written = True
                            batch.clear()

                processed += 1
                if processed % 1000 == 0:
                    print(
                        f"  {acquisition}: processed {processed:,}/"
                        f"{len(subset):,} candidates",
                        flush=True,
                    )

            if batch:
                pd.DataFrame(batch).to_csv(
                    scan_handle, sep="\t", index=False,
                    header=not header_written
                )
                header_written = True
                batch.clear()

            del traces, by_scan, scans_by_feature, rts_by_feature

    print(f"\nPer-scan output: {SCAN_OUT}", flush=True)

    # Summarise the output incrementally by reading it in chunks.
    # A candidate/anchor group can span chunks, so aggregate sufficient statistics
    # after concatenating the smaller per-group summaries.
    print("Summarising scan-level output...", flush=True)
    chunks = []
    # The per-scan output renames isotope_offset to
    # nominated_isotope_offset, so use the output column names here.
    output_candidate_keys = [
        "acquisition", "Spectrum", "feature_id", "charge",
        "nominated_isotope_offset"
    ]

    for chunk in pd.read_csv(SCAN_OUT, sep="\t", chunksize=200_000):
        group_cols = output_candidate_keys + ["anchor_offset"]
        grouped = chunk.groupby(group_cols, dropna=False).agg(
            scans_examined=("anchor_score", "size"),
            sum_score=("anchor_score", "sum"),
            median_anchor_score=("anchor_score", "median"),
            median_isotope_positions_found=("isotope_positions_found", "median"),
            median_consecutive_isotopes=("consecutive_isotopes_from_anchor", "median"),
            scans_with_3_consecutive_isotopes=(
                "consecutive_isotopes_from_anchor",
                lambda x: int((x >= 3).sum())
            ),
            median_abs_spacing_error_ppm=(
                "median_abs_adjacent_spacing_error_ppm", "median"
            ),
        ).reset_index()
        chunks.append(grouped)

    if chunks:
        partial = pd.concat(chunks, ignore_index=True)
        group_cols = candidate_keys + ["anchor_offset"]

        # Combine chunk summaries; median fields here are medians of chunk medians,
        # not exact global medians when groups span multiple chunks.
        summary = partial.groupby(group_cols, dropna=False).agg(
            scans_examined=("scans_examined", "sum"),
            sum_score=("sum_score", "sum"),
            median_anchor_score=("median_anchor_score", "mean"),
            median_isotope_positions_found=("median_isotope_positions_found", "mean"),
            median_consecutive_isotopes=("median_consecutive_isotopes", "mean"),
            scans_with_3_consecutive_isotopes=("scans_with_3_consecutive_isotopes", "sum"),
            median_abs_spacing_error_ppm=("median_abs_spacing_error_ppm", "mean"),
        ).reset_index()

        summary["mean_anchor_score"] = (
            summary["sum_score"] / summary["scans_examined"]
        )
        summary.to_csv(SUMMARY_OUT, sep="\t", index=False)

        summary["anchor_rank"] = summary.groupby(output_candidate_keys)[
            "mean_anchor_score"
        ].rank(method="min", ascending=False)

        best = summary[summary["anchor_rank"] == 1].copy()
        best = best.rename(columns={
            "anchor_offset": "best_anchor_offset",
            "mean_anchor_score": "best_anchor_score",
            "median_isotope_positions_found": "best_mean_positions_found",
            "median_consecutive_isotopes": "best_mean_consecutive_isotopes",
        })

        best[output_candidate_keys + [
            "best_anchor_offset", "best_anchor_score",
            "best_mean_positions_found", "best_mean_consecutive_isotopes"
        ]].to_csv(BEST_OUT, sep="\t", index=False)

    else:
        pd.DataFrame().to_csv(SUMMARY_OUT, sep="\t", index=False)
        pd.DataFrame().to_csv(BEST_OUT, sep="\t", index=False)

    print(f"Anchor summary: {SUMMARY_OUT}")
    print(f"Best anchor per candidate: {BEST_OUT}")
    print(f"Input candidates: {total_candidates:,}")


if __name__ == "__main__":
    main()
