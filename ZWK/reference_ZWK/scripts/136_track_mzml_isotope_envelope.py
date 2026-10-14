#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

IN_SCANS = Path(
    "reports/feature10275997325121874585_mzml_ms1_scan_inventory_v1.tsv"
)
IN_PEAKS = Path(
    "reports/feature10275997325121874585_mzml_ms1_centroids_v1.tsv"
)
OUT = Path(
    "reports/feature10275997325121874585_mzml_isotope_tracks_v1.tsv"
)
SUMMARY = Path(
    "reports/feature10275997325121874585_mzml_isotope_tracks_summary_v1.tsv"
)

EXPECTED_SPACING = 1.003354835 / 2
START = 1341.068
HALF_WIDTH = 0.075
MISSING_COST = 3.0
NTRACKS = 4


def main():
    scans = pd.read_csv(IN_SCANS, sep="\t")
    peaks = pd.read_csv(IN_PEAKS, sep="\t")

    for df in (scans, peaks):
        df["spectrum_index_zero_based"] = pd.to_numeric(
            df["spectrum_index_zero_based"], errors="coerce"
        )

    # Include all scans, even those with no peaks in the selected m/z range.
    peaks_by_scan = {
        int(idx): group.sort_values("centroid_mz")
        for idx, group in peaks.groupby("spectrum_index_zero_based")
    }

    out_rows = []
    previous = np.array([
        START,
        START + EXPECTED_SPACING,
        START + 2 * EXPECTED_SPACING,
        START + 3 * EXPECTED_SPACING,
    ])

    for _, scan in scans.sort_values("rt_seconds").iterrows():
        idx = int(scan["spectrum_index_zero_based"])
        group = peaks_by_scan.get(idx, pd.DataFrame())
        if len(group):
            mz = group["centroid_mz"].to_numpy(float)
            inten = group["centroid_intensity"].to_numpy(float)
        else:
            mz, inten = np.array([]), np.array([])

        # Build candidate sets for each expected track; missing is explicit.
        candidates = []
        for j in range(NTRACKS):
            expected = previous[j]
            ix = np.where(np.abs(mz - expected) <= HALF_WIDTH)[0]
            ix = sorted(ix, key=lambda k: abs(mz[k] - expected))[:10]
            candidates.append([(mz[k], inten[k]) for k in ix] + [(np.nan, np.nan)])

        best_set, best_cost = None, np.inf

        def search(j, chosen, cost):
            nonlocal best_set, best_cost
            if cost >= best_cost:
                return
            if j == NTRACKS:
                valid = [x[0] for x in chosen if np.isfinite(x[0])]
                if any(b <= a for a, b in zip(valid, valid[1:])):
                    return
                for k in range(NTRACKS - 1):
                    a, b = chosen[k][0], chosen[k + 1][0]
                    if np.isfinite(a) and np.isfinite(b):
                        cost += 1.5 * (((b - a) - EXPECTED_SPACING) / 0.010) ** 2
                if cost < best_cost:
                    best_set, best_cost = chosen.copy(), cost
                return

            for cmz, cint in candidates[j]:
                if not np.isfinite(cmz):
                    search(j + 1, chosen + [(cmz, cint)], cost + MISSING_COST)
                    continue

                residual = cmz - previous[j]
                c = (residual / 0.020) ** 2
                search(j + 1, chosen + [(cmz, cint)], cost + c)

        search(0, [], 0.0)

        row = {
            "spectrum_index_zero_based": idx,
            "spectrum_id": scan["spectrum_id"],
            "rt_seconds": scan["rt_seconds"],
            "raw_peak_count_in_window": len(mz),
            "n_tracks_found": sum(np.isfinite(x[0]) for x in best_set),
            "assignment_cost": best_cost,
        }

        for j, (pmz, pint) in enumerate(best_set, start=1):
            row[f"track{j}_mz"] = pmz
            row[f"track{j}_intensity"] = pint

        spacing_residuals = []
        for j in range(1, NTRACKS):
            a, b = best_set[j - 1][0], best_set[j][0]
            if np.isfinite(a) and np.isfinite(b):
                spacing = b - a
                residual = (spacing - EXPECTED_SPACING) * 1000
            else:
                spacing, residual = np.nan, np.nan
            row[f"spacing{j}_Da"] = spacing
            row[f"spacing{j}_residual_mDa"] = residual
            if np.isfinite(residual):
                spacing_residuals.append(residual)

        row["all_3_spacings_observed"] = len(spacing_residuals) == 3
        row["all_spacings_within_5mDa"] = (
            len(spacing_residuals) == 3
            and all(abs(x) <= 5 for x in spacing_residuals)
        )
        out_rows.append(row)

        # Update reference positions only when actual peaks were selected.
        for j, (pmz, _) in enumerate(best_set):
            if np.isfinite(pmz):
                previous[j] = pmz

    result = pd.DataFrame(out_rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT, sep="\t", index=False, na_rep="NA")

    summary = pd.DataFrame([
        {"metric": "scans", "value": len(result)},
        {"metric": "scans_with_4_tracks",
         "value": int((result.n_tracks_found == 4).sum())},
        {"metric": "scans_with_3_tracks",
         "value": int((result.n_tracks_found == 3).sum())},
        {"metric": "scans_with_2_or_fewer_tracks",
         "value": int((result.n_tracks_found <= 2).sum())},
        {"metric": "scans_with_all_3_spacings_observed",
         "value": int(result.all_3_spacings_observed.sum())},
        {"metric": "scans_with_all_spacings_within_5mDa",
         "value": int(result.all_spacings_within_5mDa.sum())},
        {"metric": "fraction_scans_with_all_spacings_within_5mDa",
         "value": float(result.all_spacings_within_5mDa.mean())},
    ])
    summary.to_csv(SUMMARY, sep="\t", index=False)

    print(summary.to_string(index=False))
    print("\nFirst and last 10 scans:")
    print(result[[
        "rt_seconds", "n_tracks_found",
        "track1_mz", "track2_mz", "track3_mz", "track4_mz",
        "spacing1_residual_mDa",
        "spacing2_residual_mDa",
        "spacing3_residual_mDa",
    ]].head(10).to_string(index=False))
    print(f"\nFull trace: {OUT}")
    print(f"Summary:    {SUMMARY}")


if __name__ == "__main__":
    main()
