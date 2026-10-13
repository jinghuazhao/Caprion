#!/usr/bin/env python3
"""Inspect observed isotope sequences for selected winning anchors."""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
BEST_FILE = BASE / "anchor_one_to_one_v2_best.tsv"
SCAN_FILE = BASE / "anchor_one_to_one_v2_per_scan.tsv"
TRACE_FILES = {
    "first_pilot": BASE / "isotope_trace_peaks.tsv",
    "second_acquisition": BASE / "isotope_trace_peaks_second.tsv",
}
OUT = BASE / "winning_isotope_sequence_examples_v2.tsv"

ISOTOPE_SPACING = 1.00335483507


def main():
    best = pd.read_csv(BEST_FILE, sep="\t")
    scans = pd.read_csv(SCAN_FILE, sep="\t")

    # Include the known example and one representative candidate for
    # each winning offset, selected deterministically by score.
    examples = []

    known = best[
        (best["Spectrum"] == "szwk900404i19201xms1.00492.00492.2")
        & (best["feature_id"].astype(str) == "4064338817141814541")
    ]
    if not known.empty:
        examples.append(known.iloc[0])

    for offset in sorted(best["best_anchor_offset"].dropna().unique()):
        subset = best[best["best_anchor_offset"] == offset]
        if subset.empty:
            continue
        row = subset.sort_values(
            ["best_anchor_score", "Spectrum"],
            ascending=[False, True],
        ).iloc[0]
        examples.append(row)

    selected = pd.DataFrame(examples).drop_duplicates(
        subset=[
            "acquisition", "Spectrum", "feature_id", "charge",
            "nominated_isotope_offset", "best_anchor_offset",
        ]
    )

    keys = [
        "acquisition", "Spectrum", "feature_id", "charge",
        "nominated_isotope_offset",
    ]
    scan_keys = keys + ["anchor_offset"]

    rows = []
    for row in selected.itertuples(index=False):
        matched = scans[
            (scans["acquisition"] == row.acquisition)
            & (scans["Spectrum"] == row.Spectrum)
            & (scans["feature_id"].astype(str) == str(row.feature_id))
            & (scans["charge"] == row.charge)
            & (
                scans["nominated_isotope_offset"]
                == row.nominated_isotope_offset
            )
            & (scans["anchor_offset"] == row.best_anchor_offset)
        ]
        if matched.empty:
            continue

        scan_row = matched.sort_values("scan_index").iloc[0]
        z = int(row.charge)
        expected = ISOTOPE_SPACING / z

        mzs = [
            float(scan_row[f"M{k}_observed_mz"])
            if pd.notna(scan_row[f"M{k}_observed_mz"]) else np.nan
            for k in range(5)
        ]

        for k, mz in enumerate(mzs):
            if not np.isfinite(mz):
                continue

            next_mz = mzs[k + 1] if k < 4 else np.nan
            gap = next_mz - mz if np.isfinite(next_mz) else np.nan
            gap_error = (
                (gap - expected) / expected * 1e6
                if np.isfinite(gap) else np.nan
            )

            rows.append({
                "acquisition": row.acquisition,
                "Spectrum": row.Spectrum,
                "feature_id": row.feature_id,
                "charge": z,
                "scan_index": int(scan_row["scan_index"]),
                "rt_sec": scan_row["rt_sec"],
                "nominated_isotope_offset": row.nominated_isotope_offset,
                "best_anchor_offset": row.best_anchor_offset,
                "best_anchor_score": row.best_anchor_score,
                "isotope_position": k,
                "observed_mz": mz,
                "individual_peak_error_ppm": scan_row[f"M{k}_ppm_error"],
                "next_observed_mz": next_mz,
                "observed_gap_mz": gap,
                "expected_gap_mz": expected if np.isfinite(gap) else np.nan,
                "gap_error_ppm": gap_error,
            })

    result = pd.DataFrame(rows)
    result.to_csv(OUT, sep="\t", index=False)

    print("Output:", OUT)
    print("Selected candidates:", len(selected))
    print("Peak-position rows:", len(result))

    if not result.empty:
        print("\nSelected candidate metadata:")
        print(
            result.groupby(
                [
                    "acquisition", "Spectrum", "feature_id", "charge",
                    "best_anchor_offset",
                ],
                dropna=False,
            ).agg(
                scan_index=("scan_index", "first"),
                rt_sec=("rt_sec", "first"),
                best_score=("best_anchor_score", "first"),
                positions=("isotope_position", "count"),
                median_abs_peak_error_ppm=(
                    "individual_peak_error_ppm",
                    lambda x: np.nanmedian(np.abs(x)),
                ),
                median_abs_gap_error_ppm=(
                    "gap_error_ppm",
                    lambda x: np.nanmedian(np.abs(x)),
                ),
            ).to_string())


if __name__ == "__main__":
    main()
