#!/usr/bin/env python3
"""Diagnose whether nearby candidate groups resemble shifted isotope ladders."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
PAIR_FILE = BASE / "candidate_isotope_group_overlap_review_v2.tsv"
MEMBER_FILE = BASE / "candidate_isotope_groups_v2.tsv"
OUT = BASE / "candidate_isotope_group_pair_relationships_v2.tsv"
SUMMARY = BASE / "candidate_isotope_group_pair_relationships_summary_v2.tsv"

ISOTOPE_SPACING = 1.00335483507
MAX_SHIFT_STEPS = 4
MZ_TOLERANCE_MDA = 2.0


def main():
    pairs = pd.read_csv(PAIR_FILE, sep="\t")
    members = pd.read_csv(MEMBER_FILE, sep="\t")

    required = {
        "acquisition", "candidate_ig_id_a", "candidate_ig_id_b",
        "feature_id", "charge",
    }
    missing = required - set(pairs.columns)
    if missing:
        raise ValueError(f"Pair file lacks columns: {sorted(missing)}")

    member_required = {
        "acquisition", "candidate_ig_id", "feature_id", "charge",
        "trace_mz_median", "trace_n_scans",
    }
    missing = member_required - set(members.columns)
    if missing:
        raise ValueError(f"Membership file lacks columns: {sorted(missing)}")

    lookup = {}
    for key, g in members.groupby(
        ["acquisition", "feature_id", "charge"], sort=False
    ):
        lookup[key] = {
            ig: x.sort_values("trace_mz_median")
            for ig, x in g.groupby("candidate_ig_id", sort=False)
        }

    rows = []

    for _, p in pairs.iterrows():
        key = (
            p["acquisition"], p["feature_id"], int(p["charge"])
        )
        groups = lookup.get(key, {})
        a = groups.get(p["candidate_ig_id_a"])
        b = groups.get(p["candidate_ig_id_b"])
        if a is None or b is None:
            continue

        za = int(p["charge"])
        step = ISOTOPE_SPACING / za
        ma = a["trace_mz_median"].to_numpy(float)
        mb = b["trace_mz_median"].to_numpy(float)

        # Compare each group to the other after shifting by integer
        # isotope steps. Match each trace in A to its nearest shifted
        # trace in B; this is descriptive, not a definitive assignment.
        best = None
        for shift in range(-MAX_SHIFT_STEPS, MAX_SHIFT_STEPS + 1):
            shifted_b = mb + shift * step
            distances = np.abs(ma[:, None] - shifted_b[None, :]) * 1000
            nearest = distances.min(axis=1)
            matched = nearest <= MZ_TOLERANCE_MDA
            result = {
                "shift_steps": shift,
                "n_a_traces": len(ma),
                "n_a_traces_matched": int(matched.sum()),
                "fraction_a_matched": float(matched.mean()) if len(ma) else np.nan,
                "median_nearest_error_mDa": float(np.median(nearest)),
                "max_nearest_error_mDa": float(np.max(nearest)),
            }
            if best is None or (
                result["n_a_traces_matched"],
                -result["median_nearest_error_mDa"],
            ) > (
                best["n_a_traces_matched"],
                -best["median_nearest_error_mDa"],
            ):
                best = result

        # Use observed trace RT spans if present in the group table.
        rta = (
            float(a["rt_min"].min()) if "rt_min" in a.columns else np.nan
        )
        rtb = (
            float(b["rt_min"].min()) if "rt_min" in b.columns else np.nan
        )
        rta_max = (
            float(a["rt_max"].max()) if "rt_max" in a.columns else np.nan
        )
        rtb_max = (
            float(b["rt_max"].max()) if "rt_max" in b.columns else np.nan
        )

        if np.isfinite(rta) and np.isfinite(rtb_max):
            rt_overlap = max(0.0, min(rta_max, rtb_max) - max(rta, rtb))
            rt_span = max(rta_max, rtb_max) - min(rta, rtb)
            rt_overlap_fraction = rt_overlap / rt_span if rt_span > 0 else np.nan
        else:
            rt_overlap_fraction = np.nan

        if best["fraction_a_matched"] >= 0.8 and best["median_nearest_error_mDa"] <= 1:
            label = "strong_shifted_mz_ladder_similarity"
        elif best["fraction_a_matched"] >= 0.5 and best["median_nearest_error_mDa"] <= 2:
            label = "partial_shifted_mz_ladder_similarity"
        else:
            label = "no_clear_shifted_ladder_similarity"

        rows.append({
            "acquisition": p["acquisition"],
            "feature_id": p["feature_id"],
            "charge": za,
            "candidate_ig_id_a": p["candidate_ig_id_a"],
            "candidate_ig_id_b": p["candidate_ig_id_b"],
            "n_traces_a": len(ma),
            "n_traces_b": len(mb),
            "rt_overlap_fraction": rt_overlap_fraction,
            **{f"best_{k}": v for k, v in best.items()},
            "relationship_diagnostic": label,
            "interpretation": (
                "Review only; do not merge groups or assign monoisotopic anchors."
            ),
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)

    summary = (
        out.groupby(["acquisition", "relationship_diagnostic"])
        .agg(
            pairs=("candidate_ig_id_a", "size"),
            median_fraction_a_matched=("best_fraction_a_matched", "median"),
            median_mz_error_mDa=("best_median_nearest_error_mDa", "median"),
            median_rt_overlap=("rt_overlap_fraction", "median"),
        )
        .reset_index()
    )
    summary.to_csv(SUMMARY, sep="\t", index=False)

    print(f"Pair relationships: {OUT}")
    print(f"Summary: {SUMMARY}")
    print(f"Pairs classified: {len(out)}")
    print("\n=== Relationship diagnostics ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
