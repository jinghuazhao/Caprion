#!/usr/bin/env python3
"""Jointly fit a four-peak, z=2 isotope-like envelope across MS1 scans."""

from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

INPUT = Path(
    "reports/feature10275997325121874585_all_trace_scan_centroids_v1.tsv"
)
OUT_SCAN = Path(
    "reports/feature10275997325121874585_joint_fit_scans_v1.tsv"
)
OUT_SUMMARY = Path(
    "reports/feature10275997325121874585_joint_fit_summary_v1.tsv"
)

# Approximate first-isotope m/z and isotope spacing.
FIRST_MZ = 1341.068
SPACING = 1.003354835 / 2
N_PEAKS = 4

# Candidate search range; includes nearby interfering peaks.
MZ_MIN, MZ_MAX = 1340.90, 1342.80

# Scoring scales, in Da. These are fitting tolerances, not instrument specs.
SPACING_SCALE = 0.006
CENTRE_SCALE = 0.025
CONTINUITY_SCALE = 0.008


def candidates_for_scan(g):
    """Enumerate plausible four-peak envelopes in one scan."""
    g = g[
        g["centroid_mz"].between(MZ_MIN, MZ_MAX)
    ].sort_values("centroid_mz")

    mz = g["centroid_mz"].to_numpy(float)
    inten = g["centroid_intensity"].to_numpy(float)

    # Keep a manageable set of the strongest centroids if the region is busy.
    if len(mz) > 24:
        keep = np.argsort(inten)[-24:]
        mz, inten = mz[keep], inten[keep]
        order = np.argsort(mz)
        mz, inten = mz[order], inten[order]

    candidates = []
    for ix in combinations(range(len(mz)), N_PEAKS):
        peaks = mz[list(ix)]
        intensities = inten[list(ix)]
        spacings = np.diff(peaks)

        # Reject combinations with implausible isotope spacing.
        residuals = spacings - SPACING
        if np.any(np.abs(residuals) > 0.020):
            continue

        # First peak should be near the observed trace's first track.
        first_error = peaks[0] - FIRST_MZ
        if abs(first_error) > 0.080:
            continue

        spacing_cost = np.mean((residuals / SPACING_SCALE) ** 2)
        centre_cost = (first_error / CENTRE_SCALE) ** 2

        # Weak intensity prior: avoid selecting near-zero noise when
        # a strong coherent envelope is available.
        intensity_cost = -0.03 * np.log1p(
            np.maximum(intensities, 0).sum()
        )

        cost = spacing_cost + centre_cost + intensity_cost
        candidates.append({
            "mz": peaks,
            "intensity": intensities,
            "spacings": spacings,
            "spacing_residuals": residuals,
            "local_cost": cost,
        })

    candidates.sort(key=lambda x: x["local_cost"])
    return candidates[:30]


def main():
    if not INPUT.exists():
        raise SystemExit(f"Missing input file: {INPUT}")

    df = pd.read_csv(INPUT, sep="\t")
    required = {
        "native_scan_number", "rt_seconds",
        "centroid_mz", "centroid_intensity",
    }
    if not required.issubset(df.columns):
        raise SystemExit(f"Required columns: {sorted(required)}")

    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(required))

    scan_data = []
    for scan, g in df.groupby("native_scan_number", sort=True):
        scan_data.append({
            "scan": int(scan),
            "rt": float(g["rt_seconds"].median()),
            "candidates": candidates_for_scan(g),
        })

    # Dynamic programming: choose a sequence of envelopes that balances
    # local isotope-spacing fit against smooth movement between scans.
    # This helps avoid switching to unrelated peaks from one scan to the next.
    valid = [s for s in scan_data if s["candidates"]]
    if not valid:
        raise SystemExit(
            "No four-peak combinations found. Check the m/z range and input."
        )

    costs = []
    back = []

    for i, scan in enumerate(valid):
        local = np.array(
            [c["local_cost"] for c in scan["candidates"]], dtype=float
        )

        if i == 0:
            costs.append(local)
            back.append(np.full(len(local), -1, dtype=int))
            continue

        prev_scan = valid[i - 1]
        prev_candidates = prev_scan["candidates"]
        curr_candidates = scan["candidates"]
        prev_costs = costs[-1]

        # Account for the number of seconds between observations.
        dt = max(scan["rt"] - prev_scan["rt"], 1.0)
        scale = CONTINUITY_SCALE * np.sqrt(dt / 1.4)
        scale = min(max(scale, 0.004), 0.025)

        new_cost = np.full(len(curr_candidates), np.inf)
        new_back = np.full(len(curr_candidates), -1, dtype=int)

        for j, curr in enumerate(curr_candidates):
            for k, prev in enumerate(prev_candidates):
                movement = curr["mz"] - prev["mz"]
                transition = np.mean((movement / scale) ** 2)
                total = prev_costs[k] + local[j] + transition
                if total < new_cost[j]:
                    new_cost[j] = total
                    new_back[j] = k

        costs.append(new_cost)
        back.append(new_back)

    # Backtrack the globally lowest-cost path.
    chosen_indices = [0] * len(valid)
    chosen_indices[-1] = int(np.argmin(costs[-1]))
    for i in range(len(valid) - 1, 0, -1):
        chosen_indices[i - 1] = int(
            back[i][chosen_indices[i]]
        )

    rows = []
    for scan, idx in zip(valid, chosen_indices):
        c = scan["candidates"][idx]
        row = {
            "native_scan_number": scan["scan"],
            "rt_seconds": scan["rt"],
            "local_fit_cost": c["local_cost"],
        }
        for j in range(N_PEAKS):
            row[f"peak{j+1}_mz"] = c["mz"][j]
            row[f"peak{j+1}_intensity"] = c["intensity"][j]
        for j in range(N_PEAKS - 1):
            row[f"spacing{j+1}_Da"] = c["spacings"][j]
            row[f"spacing{j+1}_residual_mDa"] = (
                c["spacing_residuals"][j] * 1000
            )
        row["median_abs_spacing_residual_mDa"] = (
            np.median(np.abs(c["spacing_residuals"])) * 1000
        )
        row["all_spacings_within_5mDa"] = bool(
            np.all(np.abs(c["spacing_residuals"]) <= 0.005)
        )
        rows.append(row)

    out = pd.DataFrame(rows).sort_values("rt_seconds")
    OUT_SCAN.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_SCAN, sep="\t", index=False, na_rep="NA")

    summary = []
    for j in range(1, N_PEAKS):
        col = f"spacing{j}_residual_mDa"
        x = out[col].dropna()
        summary.append({
            "metric": f"spacing{j}_residual_mDa",
            "n_scans": len(x),
            "median_mDa": x.median(),
            "median_abs_mDa": x.abs().median(),
            "q05_mDa": x.quantile(0.05),
            "q95_mDa": x.quantile(0.95),
            "within_3mDa": int((x.abs() <= 3).sum()),
            "within_5mDa": int((x.abs() <= 5).sum()),
        })

    # Summarize intensity co-movement: correlations alone are descriptive,
    # not proof that the peaks belong to the same precursor.
    intensities = out[
        [f"peak{j}_intensity" for j in range(1, N_PEAKS + 1)]
    ]
    for j in range(1, N_PEAKS + 1):
        x = out[f"peak{j}_mz"]
        summary.append({
            "metric": f"peak{j}_mz_stability",
            "n_scans": x.notna().sum(),
            "median_mDa": np.nanmedian(x),
            "median_abs_mDa": (
                np.nanmedian(np.abs(x - np.nanmedian(x))) * 1000
            ),
            "q05_mDa": x.quantile(0.05),
            "q95_mDa": x.quantile(0.95),
            "within_3mDa": np.nan,
            "within_5mDa": np.nan,
        })

    for a in range(N_PEAKS):
        for b in range(a + 1, N_PEAKS):
            summary.append({
                "metric": f"intensity_correlation_peak{a+1}_peak{b+1}",
                "n_scans": len(out),
                "median_mDa": intensities.iloc[:, a].corr(
                    intensities.iloc[:, b]
                ),
                "median_abs_mDa": np.nan,
                "q05_mDa": np.nan,
                "q95_mDa": np.nan,
                "within_3mDa": np.nan,
                "within_5mDa": np.nan,
            })

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUT_SUMMARY, sep="\t", index=False, na_rep="NA")

    print(f"Scans with fitted envelopes: {len(out)}")
    print(
        "Scans with all spacings within 5 mDa: "
        f"{out['all_spacings_within_5mDa'].sum()}/{len(out)}"
    )
    print("\nJoint-fit spacing residuals (mDa):")
    for j in range(1, N_PEAKS):
        x = out[f"spacing{j}_residual_mDa"].dropna()
        print(
            f"  {j}->{j+1}: median={x.median():+.3f}; "
            f"median absolute={x.abs().median():.3f}; "
            f"within 3 mDa={(x.abs() <= 3).sum()}/{len(x)}; "
            f"within 5 mDa={(x.abs() <= 5).sum()}/{len(x)}"
        )

    print("\nOutputs:")
    print(f"  {OUT_SCAN}")
    print(f"  {OUT_SUMMARY}")
    print(
        "\nCaution: this is a heuristic joint fit. Inspect the output scan "
        "table, particularly missing scans and large residuals. A good "
        "spacing fit does not by itself establish monoisotopic assignment "
        "or peptide identity."
    )


if __name__ == "__main__":
    main()
