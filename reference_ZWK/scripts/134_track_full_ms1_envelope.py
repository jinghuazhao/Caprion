#!/usr/bin/env python3
"""Track a four-peak isotope-like envelope across all centroid-table scans."""

from pathlib import Path
import numpy as np
import pandas as pd

SOURCE = Path(
    "reports/feature10275997325121874585_all_trace_scan_centroids_v1.tsv"
)
PREVIOUS_FIT = Path(
    "reports/feature10275997325121874585_joint_fit_scans_v1.tsv"
)
OUT_SCAN = Path(
    "reports/feature10275997325121874585_full_trace_envelope_v1.tsv"
)
OUT_SUMMARY = Path(
    "reports/feature10275997325121874585_full_trace_envelope_summary_v1.tsv"
)

# Empirical starting positions from the existing joint fit.
START_MZ = np.array([1341.0680, 1341.5697, 1342.0714, 1342.5731])
EXPECTED_SPACING = 1.003354835 / 2
N_TRACKS = 4

# Candidate centroids around each predicted track position.
SEARCH_HALF_WIDTH = 0.075
MAX_TRACK_STEP_DA = 0.035

# Relative scoring weights.
MZ_WEIGHT = 1.0
SPACING_WEIGHT = 1.5
INTENSITY_WEIGHT = 0.25


def robust_scale(x, floor=1.0):
    x = np.asarray(x, dtype=float)
    if len(x) < 3:
        return floor
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return max(1.4826 * mad, floor)


def main():
    if not SOURCE.exists():
        raise SystemExit(f"Missing input: {SOURCE}")

    raw = pd.read_csv(SOURCE, sep="\t")
    required = {
        "native_scan_number", "rt_seconds",
        "centroid_mz", "centroid_intensity",
    }
    if not required.issubset(raw.columns):
        raise SystemExit(f"Required columns: {sorted(required)}")

    for c in required:
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    raw = raw.dropna(subset=list(required)).copy()
    raw = raw.sort_values(["native_scan_number", "centroid_mz"])

    # Preserve every scan represented in the source file.
    scans = []
    for scan, g in raw.groupby("native_scan_number", sort=True):
        scans.append({
            "scan": int(scan),
            "rt": float(g["rt_seconds"].median()),
            "mz": g["centroid_mz"].to_numpy(float),
            "intensity": g["centroid_intensity"].to_numpy(float),
        })

    # Use existing fit, if available, to estimate expected intensity ratios
    # and typical centroid variability. These are soft priors only.
    prior_intensity = np.ones(N_TRACKS)
    mz_sigma = np.full(N_TRACKS, 0.008)
    if PREVIOUS_FIT.exists():
        fit = pd.read_csv(PREVIOUS_FIT, sep="\t")
        for j in range(N_TRACKS):
            col = f"peak{j+1}_intensity"
            if col in fit:
                vals = fit[col].dropna().to_numpy(float)
                if len(vals):
                    prior_intensity[j] = max(np.median(vals), 1.0)
            col_mz = f"peak{j+1}_mz"
            if col_mz in fit:
                vals = fit[col_mz].dropna().to_numpy(float)
                if len(vals) >= 3:
                    mz_sigma[j] = max(
                        robust_scale(vals, floor=0.001) / 1000,
                        0.003,
                    )

    rows = []
    previous_mz = START_MZ.copy()
    previous_intensity = prior_intensity.copy()
    previous_rt = None

    for s in scans:
        mz = s["mz"]
        inten = np.maximum(s["intensity"], 0)
        rt = s["rt"]

        # Predict modest m/z drift from the previous scan.
        predicted = previous_mz.copy()
        dt = 1.4 if previous_rt is None else max(rt - previous_rt, 0.1)

        # Candidate lists include None to represent a missing track.
        candidate_lists = []
        for j in range(N_TRACKS):
            target = predicted[j]
            idx = np.where(np.abs(mz - target) <= SEARCH_HALF_WIDTH)[0]
            # Keep nearest candidates, avoiding a combinatorial explosion.
            idx = sorted(idx, key=lambda k: abs(mz[k] - target))[:8]
            candidates = [(float(mz[k]), float(inten[k])) for k in idx]
            candidates.append((np.nan, np.nan))  # missing-peak option
            candidate_lists.append(candidates)

        # Enumerate ordered track assignments. The four peaks must remain
        # in increasing m/z order; missing peaks are permitted.
        best = None
        best_cost = np.inf

        def recurse(j, chosen, cost):
            nonlocal best, best_cost
            if cost >= best_cost:
                return
            if j == N_TRACKS:
                valid_mz = [p[0] for p in chosen if np.isfinite(p[0])]
                if any(b <= a for a, b in zip(valid_mz, valid_mz[1:])):
                    return

                # Score adjacent spacings when both neighbouring peaks exist.
                for k in range(N_TRACKS - 1):
                    a, b = chosen[k][0], chosen[k + 1][0]
                    if np.isfinite(a) and np.isfinite(b):
                        residual = (b - a) - EXPECTED_SPACING
                        cost += SPACING_WEIGHT * (residual / 0.010) ** 2

                if cost < best_cost:
                    best_cost, best = cost, chosen.copy()
                return

            for candidate in candidate_lists[j]:
                cmz, cint = candidate
                if not np.isfinite(cmz):
                    recurse(j + 1, chosen + [candidate], cost + 4.0)
                    continue

                step = cmz - predicted[j]
                if abs(step) > MAX_TRACK_STEP_DA:
                    continue

                mz_cost = MZ_WEIGHT * (step / max(mz_sigma[j], 0.005)) ** 2

                # Intensity continuity is a soft penalty, since isotope
                # intensities can change differently with interference.
                if cint > 0 and previous_intensity[j] > 0:
                    log_change = np.log(cint) - np.log(previous_intensity[j])
                    intensity_cost = INTENSITY_WEIGHT * min(
                        (log_change / 1.5) ** 2, 4.0
                    )
                else:
                    intensity_cost = 0.0

                recurse(
                    j + 1,
                    chosen + [candidate],
                    cost + mz_cost + intensity_cost,
                )

        recurse(0, [], 0.0)

        if best is None:
            best = [(np.nan, np.nan)] * N_TRACKS
            best_cost = np.nan

        row = {
            "native_scan_number": s["scan"],
            "rt_seconds": rt,
            "assignment_cost": best_cost,
            "n_tracks_found": sum(np.isfinite(x[0]) for x in best),
        }

        chosen_mz = []
        chosen_intensity = []
        for j, (pmz, pint) in enumerate(best, start=1):
            row[f"track{j}_mz"] = pmz
            row[f"track{j}_intensity"] = pint
            chosen_mz.append(pmz)
            chosen_intensity.append(pint)

        for j in range(N_TRACKS - 1):
            a, b = chosen_mz[j], chosen_mz[j + 1]
            if np.isfinite(a) and np.isfinite(b):
                spacing = b - a
                residual = (spacing - EXPECTED_SPACING) * 1000
            else:
                spacing, residual = np.nan, np.nan
            row[f"spacing{j+1}_Da"] = spacing
            row[f"spacing{j+1}_residual_mDa"] = residual

        residuals = [
            row[f"spacing{j}_residual_mDa"]
            for j in range(1, N_TRACKS)
            if np.isfinite(row[f"spacing{j}_residual_mDa"])
        ]
        row["median_abs_spacing_residual_mDa"] = (
            float(np.median(np.abs(residuals))) if residuals else np.nan
        )
        row["all_spacings_within_5mDa"] = (
            len(residuals) == N_TRACKS - 1
            and all(abs(x) <= 5 for x in residuals)
        )
        rows.append(row)

        # Update track priors only for observed peaks.
        for j, (pmz, pint) in enumerate(best):
            if np.isfinite(pmz):
                previous_mz[j] = pmz
                if pint > 0:
                    previous_intensity[j] = pint
        previous_rt = rt

    out = pd.DataFrame(rows)
    OUT_SCAN.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_SCAN, sep="\t", index=False, na_rep="NA")

    summary_rows = []
    for j in range(1, N_TRACKS + 1):
        col = f"track{j}_mz"
        vals = out[col].dropna()
        summary_rows.append({
            "metric": f"track{j}_mz",
            "n_scans": len(vals),
            "median": vals.median() if len(vals) else np.nan,
            "min": vals.min() if len(vals) else np.nan,
            "max": vals.max() if len(vals) else np.nan,
            "median_abs_deviation_mDa": (
                np.median(np.abs(vals - vals.median())) * 1000
                if len(vals) else np.nan
            ),
        })

    summary_rows.append({
        "metric": "n_tracks_found",
        "n_scans": len(out),
        "median": out["n_tracks_found"].median(),
        "min": out["n_tracks_found"].min(),
        "max": out["n_tracks_found"].max(),
        "median_abs_deviation_mDa": np.nan,
    })
    summary_rows.append({
        "metric": "all_spacings_within_5mDa",
        "n_scans": len(out),
        "median": out["all_spacings_within_5mDa"].mean(),
        "min": out["all_spacings_within_5mDa"].sum(),
        "max": len(out),
        "median_abs_deviation_mDa": np.nan,
    })

    pd.DataFrame(summary_rows).to_csv(
        OUT_SUMMARY, sep="\t", index=False, na_rep="NA"
    )

    print(f"Centroid rows: {len(raw):,}")
    print(f"Scans represented: {len(out):,}")
    print(f"Scans with 4 tracks: {(out['n_tracks_found'] == 4).sum()}")
    print(f"Scans with all spacings within 5 mDa: "
          f"{out['all_spacings_within_5mDa'].sum()}")
    print("\nTracks found per scan:")
    print(out["n_tracks_found"].value_counts().sort_index().to_string())
    print("\nOutput files:")
    print(f"  {OUT_SCAN}")
    print(f"  {OUT_SUMMARY}")
    print("\nWarning: this is a continuity-guided heuristic. Review the output; "
          "the algorithm can propagate an incorrect assignment across scans.")


if __name__ == "__main__":
    main()
