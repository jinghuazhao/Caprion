#!/usr/bin/env python3
"""Audit scan-level evidence for selected candidate isotope-group pairs."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "candidate_isotope_group_pair_profiles_v2.tsv"
OUT_TRACE = BASE / "candidate_isotope_group_scan_profile_trace_pairs_v2.tsv"
OUT_PAIR = BASE / "candidate_isotope_group_scan_profile_pair_summary_v2.tsv"

ISOTOPE_MASS = 1.00335483507


def mad(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return np.nan
    med = np.median(x)
    return float(np.median(np.abs(x - med)))


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def load_profiles():
    df = pd.read_csv(
        INPUT,
        sep="\t",
        dtype={
            "pair_id": str,
            "feature_id": str,
            "trace_index": str,
            "candidate_ig_id": str,
        },
    )

    required = {
        "pair_id", "acquisition", "feature_id", "pair_side",
        "trace_index", "charge", "scan_index", "rt_sec",
        "observed_mz", "peak_intensity",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    for col in [
        "charge", "scan_index", "rt_sec", "observed_mz", "peak_intensity"
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(
        subset=["charge", "scan_index", "rt_sec", "observed_mz"]
    ).copy()
    df["charge"] = df["charge"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)
    df["peak_intensity"] = df["peak_intensity"].fillna(0.0)
    return df


def trace_profiles(pg, side):
    result = []
    sg = pg[pg["pair_side"].astype(str).str.upper() == side]

    for trace_index, tg in sg.groupby("trace_index", sort=False):
        # A trace may have multiple peak records for a scan. Aggregate them
        # so the scan contributes only one intensity and one observed m/z.
        scans = (
            tg.groupby("scan_index", as_index=False)
            .agg(
                rt_sec=("rt_sec", "median"),
                observed_mz=("observed_mz", "median"),
                peak_intensity=("peak_intensity", "sum"),
            )
            .sort_values("scan_index")
        )

        if scans.empty:
            continue

        result.append({
            "trace_index": str(trace_index),
            "candidate_ig_id": (
                tg["candidate_ig_id"].iloc[0]
                if "candidate_ig_id" in tg else ""
            ),
            "mz_median": float(scans["observed_mz"].median()),
            "mz_min": float(scans["observed_mz"].min()),
            "mz_max": float(scans["observed_mz"].max()),
            "n_scans": int(scans["scan_index"].nunique()),
            "rt_min_sec": float(scans["rt_sec"].min()),
            "rt_max_sec": float(scans["rt_sec"].max()),
            "profile": scans.set_index("scan_index"),
        })

    return result


def compare_traces(a, b, spacing):
    pa = a["profile"]
    pb = b["profile"]
    common = pa.index.intersection(pb.index)

    delta_mz = b["mz_median"] - a["mz_median"]
    nearest_step = int(round(delta_mz / spacing))
    spacing_error_mDa = (delta_mz - nearest_step * spacing) * 1000.0

    row = {
        "trace_A": a["trace_index"],
        "trace_B": b["trace_index"],
        "candidate_ig_A": a["candidate_ig_id"],
        "candidate_ig_B": b["candidate_ig_id"],
        "trace_A_mz_median": a["mz_median"],
        "trace_B_mz_median": b["mz_median"],
        "delta_mz_B_minus_A": delta_mz,
        "nearest_integer_spacing_steps": nearest_step,
        "nearest_integer_spacing_error_mDa": spacing_error_mDa,
        "trace_A_n_scans": a["n_scans"],
        "trace_B_n_scans": b["n_scans"],
        "common_scan_count": len(common),
        "common_scan_fraction_of_A": (
            len(common) / a["n_scans"] if a["n_scans"] else np.nan
        ),
        "common_scan_fraction_of_B": (
            len(common) / b["n_scans"] if b["n_scans"] else np.nan
        ),
        "trace_A_rt_min_sec": a["rt_min_sec"],
        "trace_A_rt_max_sec": a["rt_max_sec"],
        "trace_B_rt_min_sec": b["rt_min_sec"],
        "trace_B_rt_max_sec": b["rt_max_sec"],
        "rt_interval_overlap_sec": max(
            0.0,
            min(a["rt_max_sec"], b["rt_max_sec"])
            - max(a["rt_min_sec"], b["rt_min_sec"]),
        ),
        "intensity_correlation_common_scans": np.nan,
        "intensity_ratio_median_B_over_A": np.nan,
        "mz_delta_median_common_scans": np.nan,
        "mz_delta_mad_common_scans_mDa": np.nan,
        "mz_delta_range_common_scans_mDa": np.nan,
    }

    if len(common):
        ca = pa.loc[common].sort_index()
        cb = pb.loc[common].sort_index()

        row["intensity_correlation_common_scans"] = safe_corr(
            ca["peak_intensity"], cb["peak_intensity"]
        )

        ia = ca["peak_intensity"].to_numpy(dtype=float)
        ib = cb["peak_intensity"].to_numpy(dtype=float)
        positive = np.isfinite(ia) & np.isfinite(ib) & (ia > 0) & (ib > 0)
        if positive.any():
            row["intensity_ratio_median_B_over_A"] = float(
                np.median(ib[positive] / ia[positive])
            )

        d = cb["observed_mz"].to_numpy() - ca["observed_mz"].to_numpy()
        d = d[np.isfinite(d)]
        if len(d):
            row["mz_delta_median_common_scans"] = float(np.median(d))
            row["mz_delta_mad_common_scans_mDa"] = mad(d) * 1000.0
            row["mz_delta_range_common_scans_mDa"] = (
                float(np.max(d) - np.min(d)) * 1000.0
            )

    return row


def main():
    df = load_profiles()
    trace_rows = []
    pair_rows = []

    for pair_id, pg in df.groupby("pair_id", sort=False):
        charges = pg["charge"].unique()
        if len(charges) != 1:
            print(f"WARNING: skipping {pair_id}: multiple charges")
            continue

        charge = int(charges[0])
        spacing = ISOTOPE_MASS / charge
        a_traces = trace_profiles(pg, "A")
        b_traces = trace_profiles(pg, "B")

        if not a_traces or not b_traces:
            print(f"WARNING: skipping {pair_id}: missing A or B traces")
            continue

        pair_comparisons = []
        for a in a_traces:
            for b in b_traces:
                row = compare_traces(a, b, spacing)
                row.update({
                    "pair_id": pair_id,
                    "acquisition": pg["acquisition"].iloc[0],
                    "feature_id": str(pg["feature_id"].iloc[0]),
                    "charge": charge,
                    "expected_spacing_mz": spacing,
                })
                trace_rows.append(row)
                pair_comparisons.append(row)

        comp = pd.DataFrame(pair_comparisons)
        common = comp["common_scan_count"].to_numpy()
        correlations = comp["intensity_correlation_common_scans"].dropna()
        spacing_errors = comp["nearest_integer_spacing_error_mDa"].abs()

        pair_rows.append({
            "pair_id": pair_id,
            "acquisition": pg["acquisition"].iloc[0],
            "feature_id": str(pg["feature_id"].iloc[0]),
            "charge": charge,
            "n_traces_A": len(a_traces),
            "n_traces_B": len(b_traces),
            "n_cross_group_trace_pairs": len(comp),
            "max_common_scan_count": int(np.max(common)) if len(common) else 0,
            "median_common_scan_count": float(np.median(common)) if len(common) else np.nan,
            "trace_pairs_with_3plus_common_scans": int(np.sum(common >= 3)),
            "median_abs_nearest_spacing_error_mDa": float(spacing_errors.median()),
            "trace_pairs_with_spacing_error_le_2mDa": int(np.sum(spacing_errors <= 2.0)),
            "median_intensity_correlation": (
                float(correlations.median()) if len(correlations) else np.nan
            ),
            "trace_pairs_with_intensity_correlation": len(correlations),
            "interpretation_limit": (
                "Shared scans and intensity similarity are descriptive; "
                "they do not prove that groups should be merged."
            ),
        })

    trace_df = pd.DataFrame(trace_rows)
    pair_df = pd.DataFrame(pair_rows)

    trace_df.to_csv(OUT_TRACE, sep="\t", index=False)
    pair_df.to_csv(OUT_PAIR, sep="\t", index=False)

    print(f"Trace-pair output: {OUT_TRACE}")
    print(f"Pair summary:      {OUT_PAIR}")
    print(f"Pairs summarized:  {len(pair_df)}")
    print(f"Cross-group trace comparisons: {len(trace_df)}")

    if not pair_df.empty:
        print("\n=== Pair-level scan-profile summary ===")
        cols = [
            "acquisition", "feature_id", "charge",
            "n_traces_A", "n_traces_B",
            "max_common_scan_count",
            "trace_pairs_with_3plus_common_scans",
            "median_abs_nearest_spacing_error_mDa",
            "trace_pairs_with_spacing_error_le_2mDa",
            "median_intensity_correlation",
        ]
        print(pair_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
