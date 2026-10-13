#!/usr/bin/env python3
"""Export per-trace isotope-index assignments for selected candidate IG pairs."""

from pathlib import Path
import itertools
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
PROFILE_FILE = BASE / "candidate_isotope_group_pair_profiles_v2.tsv"

OUT = BASE / "candidate_isotope_group_ladder_assignments_v2.tsv"
PAIR_SUMMARY = BASE / "candidate_isotope_group_ladder_assignment_summary_v2.tsv"

ISOTOPE_MASS = 1.00335483507
MAX_ISOTOPE_INDEX = 12


def fit_ladder(mzs, spacing):
    """Fit increasing integer isotope indices, anchoring the lowest m/z at M0."""
    mzs = np.asarray(mzs, dtype=float)
    n = len(mzs)
    if n == 0:
        return None

    if n == 1:
        return {
            "indices": [0],
            "residuals_mDa": [0.0],
            "median_abs_residual_mDa": np.nan,
            "p90_abs_residual_mDa": np.nan,
        }

    best = None
    best_key = None

    for tail in itertools.combinations(range(1, MAX_ISOTOPE_INDEX + 1), n - 1):
        indices = np.array((0,) + tail, dtype=int)
        residuals = (mzs - (mzs[0] + indices * spacing)) * 1000.0
        absolute = np.abs(residuals)
        key = (
            float(np.median(absolute)),
            float(np.quantile(absolute, 0.90)),
            int(indices[-1]),
        )
        if best_key is None or key < best_key:
            best_key = key
            best = {
                "indices": indices.tolist(),
                "residuals_mDa": residuals.tolist(),
                "median_abs_residual_mDa": key[0],
                "p90_abs_residual_mDa": key[1],
            }

    return best


def main():
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

    for c in ["charge", "scan_index", "rt_sec", "observed_mz", "peak_intensity"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(
        subset=["charge", "scan_index", "rt_sec", "observed_mz"]
    ).copy()
    df["charge"] = df["charge"].astype(int)
    df["scan_index"] = df["scan_index"].astype(int)

    assignment_rows = []
    summary_rows = []

    for pair_id, pg in df.groupby("pair_id", sort=False):
        charge_values = pg["charge"].dropna().unique()
        if len(charge_values) != 1:
            continue

        charge = int(charge_values[0])
        spacing = ISOTOPE_MASS / charge
        side_traces = {}

        for side in ["A", "B"]:
            side_traces[side] = []
            sg = pg[pg["pair_side"] == side]

            for _, tg in sg.groupby("trace_index", sort=False):
                tg = tg.sort_values(["rt_sec", "scan_index"])
                if tg["scan_index"].nunique() < 3:
                    continue

                side_traces[side].append({
                    "trace_index": str(tg["trace_index"].iloc[0]),
                    "candidate_ig_id": tg["candidate_ig_id"].iloc[0],
                    "mz_median": float(tg["observed_mz"].median()),
                    "mz_min": float(tg["observed_mz"].min()),
                    "mz_max": float(tg["observed_mz"].max()),
                    "scans": int(tg["scan_index"].nunique()),
                    "rt_min_sec": float(tg["rt_sec"].min()),
                    "rt_max_sec": float(tg["rt_sec"].max()),
                    "intensity_sum": float(tg["peak_intensity"].sum()),
                })

        a = sorted(side_traces["A"], key=lambda x: x["mz_median"])
        b = sorted(side_traces["B"], key=lambda x: x["mz_median"])
        combined = sorted(a + b, key=lambda x: x["mz_median"])

        if not combined:
            continue

        fits = {
            "A": fit_ladder([x["mz_median"] for x in a], spacing),
            "B": fit_ladder([x["mz_median"] for x in b], spacing),
            "combined": fit_ladder(
                [x["mz_median"] for x in combined], spacing
            ),
        }

        for model, traces in [
            ("A_only", a),
            ("B_only", b),
            ("combined", combined),
        ]:
            fit_key = model.replace("_only", "")
            fit = fits[fit_key]
            if fit is None:
                continue

            for i, trace in enumerate(traces):
                assignment_rows.append({
                    "pair_id": pair_id,
                    "acquisition": pg["acquisition"].iloc[0],
                    "feature_id": pg["feature_id"].iloc[0],
                    "charge": charge,
                    "expected_spacing_mz": spacing,
                    "model": model,
                    "pair_side": (
                        "A" if model == "A_only"
                        else "B" if model == "B_only"
                        else (
                            "A" if trace in a else "B"
                        )
                    ),
                    "candidate_ig_id": trace["candidate_ig_id"],
                    "trace_index": trace["trace_index"],
                    "trace_mz_median": trace["mz_median"],
                    "trace_mz_min": trace["mz_min"],
                    "trace_mz_max": trace["mz_max"],
                    "scans": trace["scans"],
                    "rt_min_sec": trace["rt_min_sec"],
                    "rt_max_sec": trace["rt_max_sec"],
                    "intensity_sum": trace["intensity_sum"],
                    "isotope_index_relative_to_lowest_trace": fit["indices"][i],
                    "fitted_mz": (
                        traces[0]["mz_median"]
                        + fit["indices"][i] * spacing
                    ),
                    "residual_mDa": fit["residuals_mDa"][i],
                    "median_abs_residual_mDa": fit["median_abs_residual_mDa"],
                    "p90_abs_residual_mDa": fit["p90_abs_residual_mDa"],
                    "interpretation_limit": (
                        "Relative isotope index only; lowest observed trace "
                        "is not validated as monoisotopic."
                    ),
                })

        summary_rows.append({
            "pair_id": pair_id,
            "acquisition": pg["acquisition"].iloc[0],
            "feature_id": pg["feature_id"].iloc[0],
            "charge": charge,
            "n_traces_A": len(a),
            "n_traces_B": len(b),
            "A_median_abs_residual_mDa": (
                fits["A"]["median_abs_residual_mDa"] if fits["A"] else np.nan
            ),
            "B_median_abs_residual_mDa": (
                fits["B"]["median_abs_residual_mDa"] if fits["B"] else np.nan
            ),
            "combined_median_abs_residual_mDa": (
                fits["combined"]["median_abs_residual_mDa"]
                if fits["combined"] else np.nan
            ),
            "combined_p90_abs_residual_mDa": (
                fits["combined"]["p90_abs_residual_mDa"]
                if fits["combined"] else np.nan
            ),
            "A_relative_indices": (
                ",".join(map(str, fits["A"]["indices"]))
                if fits["A"] else ""
            ),
            "B_relative_indices": (
                ",".join(map(str, fits["B"]["indices"]))
                if fits["B"] else ""
            ),
            "combined_relative_indices": (
                ",".join(map(str, fits["combined"]["indices"]))
                if fits["combined"] else ""
            ),
            "interpretation_limit": (
                "Bounded integer-index fit; diagnostic only."
            ),
        })

    assignments = pd.DataFrame(assignment_rows)
    summary = pd.DataFrame(summary_rows)

    assignments.to_csv(OUT, sep="\t", index=False)
    summary.to_csv(PAIR_SUMMARY, sep="\t", index=False)

    print(f"Trace assignments: {OUT}")
    print(f"Pair summary:      {PAIR_SUMMARY}")
    print(f"Pairs summarized: {len(summary)}")
    print(f"Assignment rows:   {len(assignments)}")

    print("\n=== Pair-level model residuals ===")
    cols = [
        "acquisition", "feature_id", "charge",
        "A_median_abs_residual_mDa",
        "B_median_abs_residual_mDa",
        "combined_median_abs_residual_mDa",
        "combined_p90_abs_residual_mDa",
        "A_relative_indices",
        "B_relative_indices",
        "combined_relative_indices",
    ]
    print(summary[[c for c in cols if c in summary.columns]].to_string(index=False))


if __name__ == "__main__":
    main()
