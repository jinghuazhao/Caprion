#!/usr/bin/env python3
"""Rescore isotope anchors using individual peak m/z errors from traces."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "monoisotopic_anchor_comparison_per_scan_fast.tsv"
TRACE_FILES = {
    "first_pilot": BASE / "isotope_trace_peaks.tsv",
    "second_acquisition": BASE / "isotope_trace_peaks_second.tsv",
}
OUT_SCAN = BASE / "anchor_mass_accuracy_per_scan.tsv"
OUT_SUMMARY = BASE / "anchor_mass_accuracy_summary.tsv"
OUT_BEST = BASE / "anchor_mass_accuracy_best.tsv"

ISOTOPE_SPACING = 1.00335483507
MAX_ISOTOPE_PEAK = 4
PPM_TOL = 10.0
MASS_ERROR_SCALE_PPM = 5.0


def nearest_peak(mzs, theoretical_mz):
    lo = theoretical_mz * (1 - PPM_TOL * 1e-6)
    hi = theoretical_mz * (1 + PPM_TOL * 1e-6)

    left = np.searchsorted(mzs, lo, side="left")
    right = np.searchsorted(mzs, hi, side="right")

    if left >= right:
        return np.nan, np.nan

    indices = np.arange(left, right)
    idx = indices[np.argmin(np.abs(mzs[indices] - theoretical_mz))]
    observed = float(mzs[idx])
    ppm = (observed - theoretical_mz) / theoretical_mz * 1e6
    return observed, float(ppm)


def load_traces():
    frames = []
    for acquisition, path in TRACE_FILES.items():
        if not path.exists():
            raise SystemExit(f"Missing trace file: {path}")

        df = pd.read_csv(path, sep="\t", low_memory=False)
        df["acquisition"] = acquisition
        df["feature_id"] = df["feature_id"].astype(str)

        for col in ["scan_index", "observed_mz"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        df = df.dropna(subset=["scan_index", "observed_mz"])
        frames.append(df)

    traces = pd.concat(frames, ignore_index=True)
    traces["scan_index"] = traces["scan_index"].astype(int)
    return traces


def main():
    anchors = pd.read_csv(INPUT, sep="\t", low_memory=False)
    anchors["feature_id"] = anchors["feature_id"].astype(str)

    traces = load_traces()

    trace_groups = {}
    for key, group in traces.groupby(
        ["acquisition", "feature_id", "scan_index"],
        sort=False,
    ):
        trace_groups[key] = np.sort(group["observed_mz"].to_numpy(float))

    output = []

    for row in anchors.itertuples(index=False):
        key = (row.acquisition, str(row.feature_id), int(row.scan_index))
        mzs = trace_groups.get(key)

        z = int(row.charge)
        step = ISOTOPE_SPACING / z
        anchor_mz = float(row.anchor_mz)

        errors = []
        observed_values = []
        for k in range(MAX_ISOTOPE_PEAK + 1):
            theoretical = anchor_mz + k * step

            if mzs is None or len(mzs) == 0:
                observed, ppm = np.nan, np.nan
            else:
                observed, ppm = nearest_peak(mzs, theoretical)

            observed_values.append(observed)
            errors.append(ppm)

        valid = np.isfinite(errors)
        n_found = int(np.sum(valid))
        median_abs_peak_ppm = (
            float(np.median(np.abs(np.asarray(errors)[valid])))
            if n_found else np.nan
        )

        positions = pd.to_numeric(
            pd.Series([row.isotope_positions_found]), errors="coerce"
        ).iloc[0]
        consecutive = int(row.consecutive_isotopes_from_anchor)

        coverage = float(positions) / 5 if pd.notna(positions) else 0.0
        continuity = consecutive / 5.0

        # A score component based on individual peak ppm errors.
        mass_quality = (
            float(np.exp(-0.5 * (median_abs_peak_ppm /
                                 MASS_ERROR_SCALE_PPM) ** 2))
            if np.isfinite(median_abs_peak_ppm) else 0.0
        )

        score = (
            0.45 * coverage
            + 0.25 * continuity
            + 0.30 * mass_quality
        )

        result = row._asdict()
        result.update({
            "trace_peaks_available": mzs is not None,
            "n_individual_peak_matches": n_found,
            "median_abs_individual_peak_error_ppm": median_abs_peak_ppm,
            "individual_mass_quality": mass_quality,
            "mass_accuracy_score": score,
        })

        for k, (observed, ppm) in enumerate(zip(observed_values, errors)):
            result[f"M{k}_observed_mz"] = observed
            result[f"M{k}_ppm_error"] = ppm

        output.append(result)

    per_scan = pd.DataFrame(output)
    per_scan.to_csv(OUT_SCAN, sep="\t", index=False)

    group_cols = [
        "acquisition", "Spectrum", "feature_id", "charge",
        "nominated_isotope_offset", "anchor_offset",
    ]

    summary = (
        per_scan.groupby(group_cols, dropna=False)
        .agg(
            mean_score=("mass_accuracy_score", "mean"),
            median_score=("mass_accuracy_score", "median"),
            median_abs_individual_peak_error_ppm=(
                "median_abs_individual_peak_error_ppm", "median"
            ),
            mean_coverage=("isotope_positions_found", "mean"),
            mean_continuity=("consecutive_isotopes_from_anchor", "mean"),
            scans=("scan_index", "size"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    candidate_cols = [
        "acquisition", "Spectrum", "feature_id",
        "charge", "nominated_isotope_offset",
    ]

    best_groups = []
    for _, group in summary.groupby(candidate_cols, dropna=False):
        maximum = group["mean_score"].max()
        winners = group[np.isclose(group["mean_score"], maximum)].copy()
        winners["n_tied_best_offsets"] = len(winners)
        best_groups.append(winners)

    best = (
        pd.concat(best_groups, ignore_index=True)
        if best_groups else pd.DataFrame()
    )
    best.to_csv(OUT_BEST, sep="\t", index=False)

    print(f"Input anchor rows: {len(anchors):,}")
    print(f"Scored scan rows: {len(per_scan):,}")
    print(f"Summary rows: {len(summary):,}")
    print(f"Best rows, retaining ties: {len(best):,}")
    print("\nTie counts:")
    print(best["n_tied_best_offsets"].value_counts().sort_index())
    print("\nWrote:")
    print(OUT_SCAN)
    print(OUT_SUMMARY)
    print(OUT_BEST)


if __name__ == "__main__":
    main()
