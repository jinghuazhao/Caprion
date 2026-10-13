#!/usr/bin/env python3
"""Audit scan-level alignment and m/z-delta integrity for selected IG pairs."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "candidate_isotope_group_selected_scan_trajectories_v2.tsv"

OUT_SCAN = BASE / "candidate_isotope_group_scan_integrity_detail_v2.tsv"
OUT_PAIR = BASE / "candidate_isotope_group_scan_integrity_summary_v2.tsv"
OUT_TRACE = BASE / "candidate_isotope_group_scan_integrity_trace_ids_v2.tsv"

SELECTED_FEATURES = {
    "10275997325121874585",
    "7118767039353031980",
    "8339255300597246646",
    "12167939761234023158",
}

REQUIRED = {
    "pair_id", "feature_id", "acquisition", "charge", "scan_index",
    "trace_A", "trace_B", "rt_A_sec", "rt_B_sec",
    "mz_A", "mz_B", "mz_delta_B_minus_A", "intensity_A", "intensity_B",
}


def main():
    df = pd.read_csv(
        INPUT,
        sep="\t",
        dtype={
            "pair_id": str,
            "feature_id": str,
            "trace_A": str,
            "trace_B": str,
        },
    )

    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    numeric = [
        "charge", "scan_index", "rt_A_sec", "rt_B_sec",
        "mz_A", "mz_B", "mz_delta_B_minus_A",
        "intensity_A", "intensity_B",
    ]
    for col in numeric:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["feature_id"] = df["feature_id"].astype(str)
    df = df[df["feature_id"].isin(SELECTED_FEATURES)].copy()
    df = df.dropna(subset=[
        "charge", "scan_index", "mz_A", "mz_B", "mz_delta_B_minus_A"
    ]).copy()

    df["charge"] = df["charge"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)

    if df.empty:
        raise ValueError("No selected feature rows found.")

    # Check duplicated keys before any aggregation.
    key = ["pair_id", "scan_index"]
    duplicate_mask = df.duplicated(key, keep=False)
    df["duplicate_pair_scan_key"] = duplicate_mask

    # Recompute delta directly from the exported m/z values.
    df["mz_delta_recomputed"] = df["mz_B"] - df["mz_A"]
    df["delta_recalculation_error_mDa"] = (
        df["mz_delta_B_minus_A"] - df["mz_delta_recomputed"]
    ) * 1000.0

    df["rt_delta_B_minus_A_sec"] = df["rt_B_sec"] - df["rt_A_sec"]

    detail_frames = []
    summary_rows = []
    trace_rows = []

    for pair_id, pg in df.groupby("pair_id", sort=False):
        pg = pg.sort_values(["scan_index", "rt_A_sec"]).copy()

        feature_id = str(pg["feature_id"].iloc[0])
        acquisition = pg["acquisition"].iloc[0]
        charge_values = pg["charge"].unique()

        trace_a_ids = sorted(pg["trace_A"].dropna().astype(str).unique())
        trace_b_ids = sorted(pg["trace_B"].dropna().astype(str).unique())

        for side, ids in [("A", trace_a_ids), ("B", trace_b_ids)]:
            for trace_id in ids:
                tg = pg[pg[f"trace_{side}"].astype(str) == trace_id]
                trace_rows.append({
                    "pair_id": pair_id,
                    "feature_id": feature_id,
                    "acquisition": acquisition,
                    "side": side,
                    "trace_id": trace_id,
                    "n_rows": len(tg),
                    "n_unique_scan_indices": tg["scan_index"].nunique(),
                    "scan_min": tg["scan_index"].min(),
                    "scan_max": tg["scan_index"].max(),
                    "candidate_ig_ids": (
                        ";".join(sorted(
                            tg.get(f"candidate_ig_{side}", pd.Series(dtype=str))
                            .dropna().astype(str).unique()
                        ))
                    ),
                })

        # Abrupt changes are evaluated in scan order. This is a diagnostic
        # threshold, not an instrument-specific acceptance criterion.
        pg["delta_change_from_previous_mDa"] = (
            pg["mz_delta_recomputed"].diff() * 1000.0
        )
        pg["large_delta_change_gt_5mDa"] = (
            pg["delta_change_from_previous_mDa"].abs() > 5.0
        )
        pg["large_delta_change_gt_20mDa"] = (
            pg["delta_change_from_previous_mDa"].abs() > 20.0
        )

        # Flag scans where the exported delta does not agree with mz_B - mz_A.
        pg["delta_recalculation_mismatch_gt_0.01mDa"] = (
            pg["delta_recalculation_error_mDa"].abs() > 0.01
        )

        detail_frames.append(pg)

        d = pg["mz_delta_recomputed"].to_numpy(dtype=float)
        changes = pg["delta_change_from_previous_mDa"].dropna().abs()

        summary_rows.append({
            "pair_id": pair_id,
            "feature_id": feature_id,
            "acquisition": acquisition,
            "charge_values": ";".join(map(str, charge_values)),
            "n_rows": len(pg),
            "n_unique_scan_indices": pg["scan_index"].nunique(),
            "duplicate_pair_scan_key_rows": int(duplicate_mask.loc[pg.index].sum()),
            "trace_A_ids": ";".join(trace_a_ids),
            "trace_B_ids": ";".join(trace_b_ids),
            "scan_min": int(pg["scan_index"].min()),
            "scan_max": int(pg["scan_index"].max()),
            "median_rt_delta_sec": float(pg["rt_delta_B_minus_A_sec"].median()),
            "max_abs_rt_delta_sec": float(pg["rt_delta_B_minus_A_sec"].abs().max()),
            "median_delta_recomputed_mz": float(np.median(d)),
            "median_delta_recomputed_mDa": float(np.median(d) * 1000.0),
            "delta_min_mDa": float(np.min(d) * 1000.0),
            "delta_max_mDa": float(np.max(d) * 1000.0),
            "delta_range_mDa": float((np.max(d) - np.min(d)) * 1000.0),
            "max_abs_adjacent_delta_change_mDa": (
                float(changes.max()) if len(changes) else np.nan
            ),
            "adjacent_changes_gt_5mDa": int(
                pg["large_delta_change_gt_5mDa"].sum()
            ),
            "adjacent_changes_gt_20mDa": int(
                pg["large_delta_change_gt_20mDa"].sum()
            ),
            "delta_recalculation_mismatches": int(
                pg["delta_recalculation_mismatch_gt_0.01mDa"].sum()
            ),
            "interpretation_limit": (
                "Integrity audit only; abrupt changes can reflect measurement "
                "variation or profile alignment issues."
            ),
        })

    detail = pd.concat(detail_frames, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    trace_ids = pd.DataFrame(trace_rows)

    detail.to_csv(OUT_SCAN, sep="\t", index=False)
    summary.to_csv(OUT_PAIR, sep="\t", index=False)
    trace_ids.to_csv(OUT_TRACE, sep="\t", index=False)

    print(f"Scan detail:  {OUT_SCAN}")
    print(f"Pair summary: {OUT_PAIR}")
    print(f"Trace IDs:    {OUT_TRACE}")
    print(f"Pairs audited: {len(summary)}")
    print(f"Detail rows:   {len(detail)}")

    print("\n=== Pair integrity summary ===")
    cols = [
        "feature_id", "n_rows", "n_unique_scan_indices",
        "duplicate_pair_scan_key_rows", "trace_A_ids", "trace_B_ids",
        "median_delta_recomputed_mDa", "delta_range_mDa",
        "max_abs_adjacent_delta_change_mDa",
        "adjacent_changes_gt_5mDa", "adjacent_changes_gt_20mDa",
        "delta_recalculation_mismatches",
    ]
    print(summary[cols].to_string(index=False))

    print("\n=== First 12 scan rows per pair ===")
    cols = [
        "feature_id", "scan_index", "trace_A", "trace_B",
        "rt_A_sec", "rt_B_sec", "mz_A", "mz_B",
        "mz_delta_recomputed", "delta_change_from_previous_mDa",
        "duplicate_pair_scan_key",
    ]
    for pair_id, pg in detail.groupby("pair_id", sort=False):
        print(f"\nPair: {pair_id}")
        print(pg[cols].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
