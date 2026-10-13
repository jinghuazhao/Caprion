#!/usr/bin/env python3
"""Compare one-ladder and two-ladder explanations for candidate IG pairs."""

from pathlib import Path
import itertools
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
PROFILE_FILE = BASE / "candidate_isotope_group_pair_profiles_v2.tsv"
OUT = BASE / "candidate_isotope_group_ladder_alignment_audit_v2.tsv"

ISOTOPE_MASS = 1.00335483507
MAX_ISOTOPE_INDEX = 12
MZ_TOLERANCE_MDA = 2.0
MIN_TRACE_SCANS = 3
MIN_COMMON_SCANS = 3


def fit_ladder(mzs, spacing, max_index=MAX_ISOTOPE_INDEX):
    """Fit a ladder with integer isotope indices, allowing an arbitrary offset."""
    mzs = np.asarray(mzs, dtype=float)
    n = len(mzs)
    if n == 0:
        return None

    # A one-trace ladder is not informative about isotope spacing.
    if n == 1:
        return {
            "indices": [0],
            "residuals_mDa": [0.0],
            "median_abs_residual_mDa": np.nan,
            "fraction_within_2mDa": np.nan,
            "index_span": 0,
        }

    best = None
    best_key = None

    # Fix the lowest-m/z trace at index 0 and enumerate increasing indices.
    # This is deliberately bounded; a successful fit is not proof of identity.
    for tail in itertools.combinations(
        range(1, max_index + 1), n - 1
    ):
        indices = np.array((0,) + tail, dtype=int)
        residuals = (mzs - (mzs[0] + indices * spacing)) * 1000.0
        abs_residuals = np.abs(residuals)
        key = (
            float(np.median(abs_residuals)),
            float(np.quantile(abs_residuals, 0.90)),
            int(indices[-1]),
        )
        if best_key is None or key < best_key:
            best_key = key
            best = {
                "indices": indices.tolist(),
                "residuals_mDa": residuals.tolist(),
                "median_abs_residual_mDa": key[0],
                "p90_abs_residual_mDa": key[1],
                "fraction_within_2mDa": float(
                    np.mean(abs_residuals <= MZ_TOLERANCE_MDA)
                ),
                "index_span": int(indices[-1]),
            }

    return best


def summarize_trace_group(g):
    """Summarize one trace, retaining scan-level observations for co-observation."""
    g = g.sort_values(["rt_sec", "scan_index"])
    return {
        "trace_index": str(g["trace_index"].iloc[0]),
        "mz_median": float(g["observed_mz"].median()),
        "scans": int(g["scan_index"].nunique()),
        "rt_min_sec": float(g["rt_sec"].min()),
        "rt_max_sec": float(g["rt_sec"].max()),
        "scan_set": set(g["scan_index"].astype(int)),
    }


def rt_overlap_fraction(a, b):
    lo = max(a["rt_min_sec"], b["rt_min_sec"])
    hi = min(a["rt_max_sec"], b["rt_max_sec"])
    overlap = max(0.0, hi - lo)
    span = max(a["rt_max_sec"], b["rt_max_sec"]) - min(
        a["rt_min_sec"], b["rt_min_sec"]
    )
    return overlap / span if span > 0 else np.nan


def main():
    if not PROFILE_FILE.exists():
        raise FileNotFoundError(
            f"Missing {PROFILE_FILE}; run script 37 first."
        )

    df = pd.read_csv(
        PROFILE_FILE,
        sep="\t",
        dtype={
            "feature_id": str,
            "trace_index": str,
            "pair_id": str,
            "candidate_ig_id": str,
        },
    )

    required = {
        "acquisition", "feature_id", "charge", "pair_id", "pair_side",
        "candidate_ig_id", "trace_index", "scan_index",
        "rt_sec", "observed_mz", "peak_intensity",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Profile file lacks columns: {sorted(missing)}")

    for c in ["charge", "scan_index", "rt_sec", "observed_mz", "peak_intensity"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(
        subset=["charge", "scan_index", "rt_sec", "observed_mz"]
    ).copy()
    df["charge"] = df["charge"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)

    results = []

    for pair_id, pg in df.groupby("pair_id", sort=False):
        sides = {}
        for side in ["A", "B"]:
            sg = pg[pg["pair_side"] == side]
            traces = []
            for _, tg in sg.groupby("trace_index", sort=False):
                t = summarize_trace_group(tg)
                if t["scans"] >= MIN_TRACE_SCANS:
                    traces.append(t)
            sides[side] = sorted(traces, key=lambda x: x["mz_median"])

        a, b = sides["A"], sides["B"]
        charge_values = pg["charge"].dropna().unique()
        if len(charge_values) != 1 or not a or not b:
            results.append({
                "pair_id": pair_id,
                "audit_status": "insufficient_trace_or_charge_data",
            })
            continue

        charge = int(charge_values[0])
        spacing = ISOTOPE_MASS / charge

        all_traces = sorted(a + b, key=lambda x: x["mz_median"])
        mz_all = [t["mz_median"] for t in all_traces]
        mz_a = [t["mz_median"] for t in a]
        mz_b = [t["mz_median"] for t in b]

        combined_fit = fit_ladder(mz_all, spacing)
        fit_a = fit_ladder(mz_a, spacing)
        fit_b = fit_ladder(mz_b, spacing)

        # Pairwise RT and scan co-observation.
        overlaps = []
        common_scans = []
        for ta in a:
            for tb in b:
                overlaps.append(rt_overlap_fraction(ta, tb))
                common_scans.append(len(ta["scan_set"] & tb["scan_set"]))

        # A simple two-ladder fit score, weighted by the number of traces
        # in each side. Single-trace sides have no informative spacing score.
        separate_components = []
        for traces, fit in [(a, fit_a), (b, fit_b)]:
            if len(traces) >= 2 and fit is not None:
                separate_components.append(
                    fit["median_abs_residual_mDa"] * len(traces)
                )

        separate_score = (
            sum(separate_components) / sum(
                len(x) for x in [a, b] if len(x) >= 2
            )
            if separate_components else np.nan
        )

        combined_score = (
            combined_fit["median_abs_residual_mDa"]
            if combined_fit else np.nan
        )

        # Heuristic only: require low residuals and a meaningful score gap.
        if not np.isfinite(combined_score) or not np.isfinite(separate_score):
            label = "inconclusive"
        elif (
            combined_score <= MZ_TOLERANCE_MDA
            and combined_score + 0.5 < separate_score
            and max(common_scans, default=0) >= MIN_COMMON_SCANS
        ):
            label = "combined_ladder_candidate"
        elif (
            separate_score <= MZ_TOLERANCE_MDA
            and separate_score + 0.5 < combined_score
        ):
            label = "two_ladders_candidate"
        else:
            label = "inconclusive"

        results.append({
            "pair_id": pair_id,
            "acquisition": pg["acquisition"].iloc[0],
            "feature_id": pg["feature_id"].iloc[0],
            "charge": charge,
            "expected_spacing_mz": spacing,
            "n_traces_A": len(a),
            "n_traces_B": len(b),
            "n_traces_combined": len(all_traces),
            "combined_ladder_median_abs_residual_mDa": combined_score,
            "combined_ladder_p90_abs_residual_mDa": (
                combined_fit["p90_abs_residual_mDa"]
                if combined_fit else np.nan
            ),
            "combined_ladder_fraction_within_2mDa": (
                combined_fit["fraction_within_2mDa"]
                if combined_fit else np.nan
            ),
            "combined_ladder_indices": (
                ",".join(map(str, combined_fit["indices"]))
                if combined_fit else ""
            ),
            "A_ladder_median_abs_residual_mDa": (
                fit_a["median_abs_residual_mDa"] if fit_a else np.nan
            ),
            "A_ladder_indices": (
                ",".join(map(str, fit_a["indices"])) if fit_a else ""
            ),
            "B_ladder_median_abs_residual_mDa": (
                fit_b["median_abs_residual_mDa"] if fit_b else np.nan
            ),
            "B_ladder_indices": (
                ",".join(map(str, fit_b["indices"])) if fit_b else ""
            ),
            "separate_ladders_weighted_residual_mDa": separate_score,
            "trace_pair_rt_overlap_median": (
                float(np.nanmedian(overlaps)) if overlaps else np.nan
            ),
            "trace_pair_common_scans_max": (
                max(common_scans) if common_scans else 0
            ),
            "heuristic_relationship": label,
            "interpretation_limit": (
                "Bounded m/z fit only; not a validated IG merge or mono-anchor."
            ),
        })

    out = pd.DataFrame(results)
    out.to_csv(OUT, sep="\t", index=False)

    print(f"Output: {OUT}")
    print(f"Pairs audited: {len(out)}")
    if "heuristic_relationship" in out.columns:
        print("\n=== Heuristic relationship counts ===")
        print(out["heuristic_relationship"].value_counts(dropna=False).to_string())

    cols = [
        "acquisition", "feature_id", "charge",
        "n_traces_A", "n_traces_B",
        "combined_ladder_median_abs_residual_mDa",
        "separate_ladders_weighted_residual_mDa",
        "trace_pair_rt_overlap_median",
        "trace_pair_common_scans_max",
        "heuristic_relationship",
    ]
    print("\n=== Alignment audit preview ===")
    print(out[[c for c in cols if c in out.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
