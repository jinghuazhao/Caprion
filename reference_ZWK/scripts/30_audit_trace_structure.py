#!/usr/bin/env python3
"""Audit trace and feature structure in existing OpenMS peak tables."""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUTS = [
    ("first_pilot", BASE / "isotope_trace_peaks.tsv"),
    ("second_acquisition", BASE / "isotope_trace_peaks_second.tsv"),
]
OUT = BASE / "trace_structure_audit_v2.tsv"
SUMMARY = BASE / "trace_structure_summary_v2.tsv"


def main():
    detail = []
    summary = []

    for acquisition, path in INPUTS:
        print(f"Reading {path}")
        d = pd.read_csv(path, sep="\t")

        required = {
            "feature_id", "trace_index", "feature_charge",
            "scan_index", "observed_mz", "peak_intensity",
        }
        missing = required - set(d.columns)
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")

        d = d.dropna(subset=list(required)).copy()
        d["feature_charge"] = pd.to_numeric(
            d["feature_charge"], errors="coerce"
        )
        d["observed_mz"] = pd.to_numeric(
            d["observed_mz"], errors="coerce"
        )
        d["peak_intensity"] = pd.to_numeric(
            d["peak_intensity"], errors="coerce"
        )
        d = d.dropna(
            subset=["feature_charge", "observed_mz", "peak_intensity"]
        )

        # Summarise trace-level mass and scan coverage.
        trace_keys = ["feature_id", "trace_index", "feature_charge"]
        traces = d.groupby(trace_keys, dropna=False).agg(
            n_rows=("scan_index", "size"),
            n_scans=("scan_index", "nunique"),
            mz_min=("observed_mz", "min"),
            mz_median=("observed_mz", "median"),
            mz_max=("observed_mz", "max"),
            intensity_sum=("peak_intensity", "sum"),
        ).reset_index()
        traces.insert(0, "acquisition", acquisition)
        detail.append(traces)

        # How many traces are associated with each feature?
        per_feature = traces.groupby(
            ["feature_id", "feature_charge"], dropna=False
        ).agg(
            n_traces=("trace_index", "nunique"),
            n_scans_across_traces=("n_scans", "max"),
            mz_span=("mz_max", lambda x: np.nan),
        ).reset_index()

        # Use the observed trace m/z ranges to calculate feature-level span.
        ranges = traces.groupby(
            ["feature_id", "feature_charge"], dropna=False
        ).agg(
            feature_mz_min=("mz_min", "min"),
            feature_mz_max=("mz_max", "max"),
        ).reset_index()
        per_feature = per_feature.drop(columns=["mz_span"]).merge(
            ranges, on=["feature_id", "feature_charge"], how="left"
        )
        per_feature["feature_mz_span"] = (
            per_feature["feature_mz_max"] - per_feature["feature_mz_min"]
        )

        s = {
            "acquisition": acquisition,
            "input_rows": len(d),
            "unique_features": d["feature_id"].nunique(),
            "unique_traces": len(traces),
            "features_with_multiple_traces": int(
                (per_feature["n_traces"] > 1).sum()
            ),
            "fraction_features_multiple_traces": (
                (per_feature["n_traces"] > 1).mean()
                if len(per_feature) else np.nan
            ),
            "median_traces_per_feature": per_feature["n_traces"].median(),
            "max_traces_per_feature": per_feature["n_traces"].max(),
            "median_scans_per_trace": traces["n_scans"].median(),
            "traces_seen_in_2plus_scans": int(
                (traces["n_scans"] >= 2).sum()
            ),
            "traces_seen_in_3plus_scans": int(
                (traces["n_scans"] >= 3).sum()
            ),
            "traces_seen_in_5plus_scans": int(
                (traces["n_scans"] >= 5).sum()
            ),
            "median_feature_mz_span": per_feature["feature_mz_span"].median(),
        }
        summary.append(s)
        print(pd.Series(s).to_string())

    pd.concat(detail, ignore_index=True).to_csv(
        OUT, sep="\t", index=False
    )
    pd.DataFrame(summary).to_csv(
        SUMMARY, sep="\t", index=False
    )
    print(f"\nTrace detail: {OUT}")
    print(f"Summary: {SUMMARY}")


if __name__ == "__main__":
    main()
