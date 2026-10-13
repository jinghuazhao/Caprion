#!/usr/bin/env python3
"""Re-audit candidate IG pairs with RT ranges derived from trace peaks."""

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path("openms_pilot")
PAIR_FILE = BASE / "candidate_isotope_group_overlap_review_v2.tsv"
MEMBER_FILE = BASE / "candidate_isotope_groups_v2.tsv"
OUT = BASE / "candidate_isotope_group_pair_relationships_rt_v2.tsv"
SUMMARY = BASE / "candidate_isotope_group_pair_relationships_rt_summary_v2.tsv"

TRACE_FILES = {
    "first_pilot": BASE / "isotope_trace_peaks.tsv",
    "second_acquisition": BASE / "isotope_trace_peaks_second.tsv",
}

ISOTOPE_SPACING = 1.00335483507
MAX_SHIFT_STEPS = 4
MZ_TOLERANCE_MDA = 2.0


def main():
    pairs = pd.read_csv(PAIR_FILE, sep="\t", dtype={"feature_id": str})
    members = pd.read_csv(
        MEMBER_FILE, sep="\t",
        dtype={"feature_id": str, "trace_index": str},
    )

    # Derive RT bounds from observed trace peaks; add acquisition from filename.
    rt_parts = []
    for acquisition, path in TRACE_FILES.items():
        if not path.exists():
            raise FileNotFoundError(path)

        x = pd.read_csv(
            path, sep="\t",
            usecols=["feature_id", "trace_index", "rt_sec"],
            dtype={"feature_id": str, "trace_index": str},
        )
        x["acquisition"] = acquisition
        x["rt_sec"] = pd.to_numeric(x["rt_sec"], errors="coerce")
        rt_parts.append(x)

    peaks = pd.concat(rt_parts, ignore_index=True)
    rt = (
        peaks.groupby(
            ["acquisition", "feature_id", "trace_index"],
            as_index=False,
        )
        .agg(
            trace_rt_min=("rt_sec", "min"),
            trace_rt_max=("rt_sec", "max"),
            trace_rt_n=("rt_sec", "count"),
        )
    )

    members["acquisition"] = members["acquisition"].astype(str)
    rt["acquisition"] = rt["acquisition"].astype(str)

    members_rt = members.merge(
        rt,
        on=["acquisition", "feature_id", "trace_index"],
        how="left",
        validate="many_to_one",
    )

    missing = int(members_rt["trace_rt_min"].isna().sum())
    print(f"Candidate trace memberships: {len(members_rt)}")
    print(f"Memberships without matched RT: {missing}")

    lookup = {}
    for key, g in members_rt.groupby(
        ["acquisition", "feature_id", "charge"], sort=False
    ):
        lookup[key] = {
            ig: x.sort_values("trace_mz_median")
            for ig, x in g.groupby("candidate_ig_id", sort=False)
        }

    rows = []
    unmatched_pairs = 0

    for _, p in pairs.iterrows():
        key = (
            str(p["acquisition"]),
            str(p["feature_id"]),
            int(p["charge"]),
        )
        groups = lookup.get(key, {})
        a = groups.get(p["candidate_ig_id_a"])
        b = groups.get(p["candidate_ig_id_b"])

        if a is None or b is None:
            unmatched_pairs += 1
            continue

        charge = int(p["charge"])
        isotope_step = ISOTOPE_SPACING / charge
        ma = a["trace_mz_median"].to_numpy(float)
        mb = b["trace_mz_median"].to_numpy(float)

        # Compare the two m/z ladders under integer isotope shifts.
        best = None
        for shift in range(-MAX_SHIFT_STEPS, MAX_SHIFT_STEPS + 1):
            shifted_b = mb + shift * isotope_step
            errors = np.abs(ma[:, None] - shifted_b[None, :]) * 1000
            nearest = errors.min(axis=1)
            matched = nearest <= MZ_TOLERANCE_MDA

            candidate = {
                "shift_steps": shift,
                "fraction_a_matched": float(matched.mean()),
                "n_a_traces_matched": int(matched.sum()),
                "median_nearest_error_mDa": float(np.median(nearest)),
                "max_nearest_error_mDa": float(np.max(nearest)),
            }
            if best is None or (
                candidate["n_a_traces_matched"],
                -candidate["median_nearest_error_mDa"],
            ) > (
                best["n_a_traces_matched"],
                -best["median_nearest_error_mDa"],
            ):
                best = candidate

        amin = pd.to_numeric(a["trace_rt_min"], errors="coerce").min()
        amax = pd.to_numeric(a["trace_rt_max"], errors="coerce").max()
        bmin = pd.to_numeric(b["trace_rt_min"], errors="coerce").min()
        bmax = pd.to_numeric(b["trace_rt_max"], errors="coerce").max()

        overlap = union = fraction = np.nan
        if np.isfinite([amin, amax, bmin, bmax]).all():
            overlap = max(0.0, min(amax, bmax) - max(amin, bmin))
            union = max(amax, bmax) - min(amin, bmin)
            fraction = overlap / union if union > 0 else (
                1.0 if amin == bmin else 0.0
            )

        if best["fraction_a_matched"] >= 0.8 and best["median_nearest_error_mDa"] <= 1:
            mz_label = "strong_shifted_mz_similarity"
        elif best["fraction_a_matched"] >= 0.5 and best["median_nearest_error_mDa"] <= 2:
            mz_label = "partial_shifted_mz_similarity"
        else:
            mz_label = "no_clear_shifted_mz_similarity"

        if not np.isfinite(fraction):
            rt_label = "RT_unavailable"
        elif fraction >= 0.5:
            rt_label = "substantial_RT_overlap"
        elif fraction > 0:
            rt_label = "partial_RT_overlap"
        else:
            rt_label = "no_RT_overlap"

        rows.append({
            "acquisition": p["acquisition"],
            "feature_id": p["feature_id"],
            "charge": charge,
            "candidate_ig_id_a": p["candidate_ig_id_a"],
            "candidate_ig_id_b": p["candidate_ig_id_b"],
            "n_traces_a": len(ma),
            "n_traces_b": len(mb),
            "a_rt_min_sec": amin,
            "a_rt_max_sec": amax,
            "b_rt_min_sec": bmin,
            "b_rt_max_sec": bmax,
            "rt_overlap_sec": overlap,
            "rt_union_span_sec": union,
            "rt_overlap_fraction": fraction,
            **{f"best_{k}": v for k, v in best.items()},
            "mz_relationship": mz_label,
            "rt_relationship": rt_label,
            "interpretation": "Diagnostic only; no merging or anchor assignment.",
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)

    print(f"Detailed audit: {OUT}")
    print(f"Pairs expanded: {len(out)}")
    print(f"Pairs without matching groups: {unmatched_pairs}")

    if out.empty:
        return

    summary = (
        out.groupby(
            ["acquisition", "mz_relationship", "rt_relationship"],
            dropna=False,
        )
        .agg(
            pairs=("candidate_ig_id_a", "size"),
            median_mz_match_fraction=("best_fraction_a_matched", "median"),
            median_mz_error_mDa=("best_median_nearest_error_mDa", "median"),
            median_rt_overlap_fraction=("rt_overlap_fraction", "median"),
        )
        .reset_index()
    )
    summary.to_csv(SUMMARY, sep="\t", index=False)

    print(f"Summary: {SUMMARY}")
    print("\n=== Pair relationship summary ===")
    print(summary.to_string(index=False))
    print("\n=== RT overlap fraction ===")
    print(out["rt_overlap_fraction"].describe().to_string())


if __name__ == "__main__":
    main()
