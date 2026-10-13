#!/usr/bin/env python3
"""Audit adjacent isotope-gap deviations for all winning v2 scans."""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "winning_anchor_spacing_audit_v2.tsv"
OUT = BASE / "all_winning_gap_deviations_v2.tsv"
SUMMARY = BASE / "all_winning_gap_deviation_summary_v2.tsv"

ISOTOPE_SPACING = 1.00335483507


def main():
    d = pd.read_csv(INPUT, sep="\t")
    records = []

    for row in d.itertuples(index=False):
        z = int(row.charge)
        if z <= 0:
            continue

        expected = ISOTOPE_SPACING / z
        for k in range(4):
            a = getattr(row, f"M{k}_observed_mz")
            b = getattr(row, f"M{k + 1}_observed_mz")
            if not (np.isfinite(a) and np.isfinite(b)):
                continue

            gap = b - a
            deviation = gap - expected
            midpoint = (a + b) / 2

            records.append({
                "acquisition": row.acquisition,
                "Spectrum": row.Spectrum,
                "feature_id": row.feature_id,
                "charge": z,
                "best_anchor_offset": row.best_anchor_offset,
                "scan_index": row.scan_index,
                "gap_index": f"M{k}-M{k+1}",
                "observed_gap_mz": gap,
                "expected_gap_mz": expected,
                "gap_deviation_mDa": deviation * 1000,
                "abs_gap_deviation_mDa": abs(deviation) * 1000,
                "gap_deviation_ppm_at_peak_mz":
                    deviation / midpoint * 1e6,
                "abs_gap_deviation_ppm_at_peak_mz":
                    abs(deviation / midpoint * 1e6),
            })

    gaps = pd.DataFrame(records)
    gaps.to_csv(OUT, sep="\t", index=False)

    print("Gap detail:", OUT)
    print("Gap rows:", len(gaps))

    if gaps.empty:
        return

    print("\n=== Overall absolute gap deviation ===")
    for col in [
        "abs_gap_deviation_mDa",
        "abs_gap_deviation_ppm_at_peak_mz",
    ]:
        print(f"\n{col}")
        print(gaps[col].describe(
            percentiles=[0.5, 0.75, 0.9, 0.95, 0.99]
        ).to_string())

    print("\n=== Fractions within absolute gap-deviation thresholds ===")
    for threshold in [0.25, 0.5, 1.0, 2.0, 5.0]:
        fraction = (gaps["abs_gap_deviation_mDa"] <= threshold).mean()
        print(f"Within {threshold:g} mDa: {fraction:.4f}")

    summary = gaps.groupby("best_anchor_offset").agg(
        gap_rows=("abs_gap_deviation_mDa", "size"),
        median_abs_gap_mDa=("abs_gap_deviation_mDa", "median"),
        p90_abs_gap_mDa=(
            "abs_gap_deviation_mDa", lambda x: x.quantile(0.90)
        ),
        p95_abs_gap_mDa=(
            "abs_gap_deviation_mDa", lambda x: x.quantile(0.95)
        ),
        median_abs_gap_ppm_at_peak_mz=(
            "abs_gap_deviation_ppm_at_peak_mz",
            lambda x: x.abs().median(),
        ),
    ).reset_index()

    summary.to_csv(SUMMARY, sep="\t", index=False)
    print("\nSummary by winning offset:", SUMMARY)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
