#!/usr/bin/env python3
"""Audit charge-dependent m/z spacing between traces within each feature."""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUTS = [
    ("first_pilot", BASE / "trace_structure_audit_v2.tsv"),
    ("second_acquisition", BASE / "trace_structure_audit_v2.tsv"),
]
OUT = BASE / "within_feature_trace_spacing_v2.tsv"
SUMMARY = BASE / "within_feature_trace_spacing_summary_v2.tsv"

ISOTOPE_SPACING = 1.00335483507
MAX_STEPS = 4


def main():
    traces = pd.read_csv(INPUTS[0][1], sep="\t")

    required = {
        "acquisition", "feature_id", "trace_index",
        "feature_charge", "mz_median", "n_scans",
    }
    missing = required - set(traces.columns)
    if missing:
        raise ValueError(f"Trace table lacks columns: {sorted(missing)}")

    # trace_structure_audit_v2.tsv contains both acquisitions.
    rows = []

    for (acquisition, feature_id, charge), g in traces.groupby(
        ["acquisition", "feature_id", "feature_charge"], sort=False
    ):
        charge = int(charge)
        if charge <= 0 or len(g) < 2:
            continue

        g = g.sort_values("mz_median")
        values = g[
            ["trace_index", "mz_median", "n_scans"]
        ].to_dict("records")
        expected_step = ISOTOPE_SPACING / charge

        for i in range(len(values) - 1):
            for j in range(i + 1, len(values)):
                a = values[i]
                b = values[j]
                observed_gap = b["mz_median"] - a["mz_median"]
                if observed_gap <= 0:
                    continue

                steps = int(round(observed_gap / expected_step))
                if steps < 1 or steps > MAX_STEPS:
                    continue

                expected_gap = steps * expected_step
                deviation = observed_gap - expected_gap

                rows.append({
                    "acquisition": acquisition,
                    "feature_id": feature_id,
                    "charge": charge,
                    "lower_trace_index": a["trace_index"],
                    "upper_trace_index": b["trace_index"],
                    "lower_trace_mz": a["mz_median"],
                    "upper_trace_mz": b["mz_median"],
                    "lower_trace_scans": a["n_scans"],
                    "upper_trace_scans": b["n_scans"],
                    "observed_gap_mz": observed_gap,
                    "isotope_steps": steps,
                    "expected_gap_mz": expected_gap,
                    "gap_deviation_mDa": deviation * 1000,
                    "abs_gap_deviation_mDa": abs(deviation) * 1000,
                    "gap_deviation_ppm_at_mz": (
                        deviation / ((a["mz_median"] + b["mz_median"]) / 2)
                        * 1e6
                    ),
                    "both_traces_seen_in_5plus_scans": (
                        a["n_scans"] >= 5 and b["n_scans"] >= 5
                    ),
                })

    result = pd.DataFrame(rows)
    result.to_csv(OUT, sep="\t", index=False)
    print(f"Detail: {OUT}")
    print(f"Trace-pair rows: {len(result)}")

    if result.empty:
        return

    print("\n=== Spacing deviations by isotope step ===")
    summary = result.groupby(
        ["acquisition", "isotope_steps"]
    ).agg(
        trace_pairs=("feature_id", "size"),
        features=("feature_id", "nunique"),
        median_abs_gap_mDa=("abs_gap_deviation_mDa", "median"),
        p90_abs_gap_mDa=(
            "abs_gap_deviation_mDa", lambda x: x.quantile(0.90)
        ),
        p95_abs_gap_mDa=(
            "abs_gap_deviation_mDa", lambda x: x.quantile(0.95)
        ),
        fraction_within_1mDa=(
            "abs_gap_deviation_mDa", lambda x: (x <= 1).mean()
        ),
        fraction_within_2mDa=(
            "abs_gap_deviation_mDa", lambda x: (x <= 2).mean()
        ),
        fraction_both_traces_persistent=(
            "both_traces_seen_in_5plus_scans", "mean"
        ),
    ).reset_index()

    summary.to_csv(SUMMARY, sep="\t", index=False)
    print(f"Summary: {SUMMARY}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
