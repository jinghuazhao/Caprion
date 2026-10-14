#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

INPUT = Path(
    "reports/feature10275997325121874585_mzml_ms1_centroids_v1.tsv"
)
OUT_SCAN = Path(
    "reports/feature10275997325121874585_raw_isotope_spacing_v1.tsv"
)
OUT_SUMMARY = Path(
    "reports/feature10275997325121874585_raw_isotope_summary_v1.tsv"
)
OUT_CORR = Path(
    "reports/feature10275997325121874585_raw_isotope_correlations_v1.tsv"
)

EXPECTED = 1.003354835 / 2
CENTRES = np.array([1341.0680, 1341.5697, 1342.0714, 1342.5731])
HALF_WIDTH = 0.025

# Main chromatographic region and wider context.
MAIN_RT = (755.0, 769.0)
WIDE_RT = (740.0, 790.0)


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return np.nan
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def spearman(x, y):
    x, y = pd.Series(x), pd.Series(y)
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return np.nan
    return float(x[ok].corr(y[ok], method="spearman"))


def main():
    d = pd.read_csv(INPUT, sep="\t")
    for col in ["rt_seconds", "centroid_mz", "centroid_intensity"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d = d.dropna(
        subset=["rt_seconds", "centroid_mz", "centroid_intensity"]
    ).copy()

    # For each MS1 scan, choose the nearest raw centroid to each fixed
    # provisional centre. A missing match remains NA; no peak is invented.
    rows = []
    for scan_idx, g in d.groupby("spectrum_index_zero_based", sort=True):
        rt = float(g["rt_seconds"].median())
        row = {
            "spectrum_index_zero_based": scan_idx,
            "spectrum_id": g["spectrum_id"].iloc[0],
            "rt_seconds": rt,
        }
        selected = []
        for j, centre in enumerate(CENTRES, start=1):
            near = g.loc[(g["centroid_mz"] - centre).abs() <= HALF_WIDTH].copy()
            if near.empty:
                mz, intensity = np.nan, np.nan
            else:
                near["distance"] = (near["centroid_mz"] - centre).abs()
                # Nearest centroid, with higher intensity breaking exact ties.
                near = near.sort_values(
                    ["distance", "centroid_intensity"],
                    ascending=[True, False],
                )
                mz = float(near.iloc[0]["centroid_mz"])
                intensity = float(near.iloc[0]["centroid_intensity"])
            row[f"peak{j}_mz"] = mz
            row[f"peak{j}_intensity"] = intensity
            selected.append((mz, intensity))

        spacings = []
        residuals = []
        for j in range(3):
            a, b = selected[j][0], selected[j + 1][0]
            if np.isfinite(a) and np.isfinite(b):
                spacing = b - a
                residual = (spacing - EXPECTED) * 1000
            else:
                spacing, residual = np.nan, np.nan
            row[f"spacing{j+1}_Da"] = spacing
            row[f"spacing{j+1}_residual_mDa"] = residual
            spacings.append(spacing)
            residuals.append(residual)

        row["n_peaks_found"] = sum(
            np.isfinite(mz) for mz, _ in selected
        )
        observed = [x for x in residuals if np.isfinite(x)]
        row["n_spacings_observed"] = len(observed)
        row["max_abs_spacing_residual_mDa"] = (
            max(abs(x) for x in observed) if observed else np.nan
        )
        row["all_3_spacings_observed"] = len(observed) == 3
        row["all_spacings_within_3mDa"] = (
            len(observed) == 3 and all(abs(x) <= 3 for x in observed)
        )
        row["all_spacings_within_5mDa"] = (
            len(observed) == 3 and all(abs(x) <= 5 for x in observed)
        )
        row["mean_intensity_4_peaks"] = (
            np.mean([x[1] for x in selected])
            if all(np.isfinite(x[1]) for x in selected)
            else np.nan
        )
        rows.append(row)

    result = pd.DataFrame(rows).sort_values("rt_seconds")
    OUT_SCAN.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT_SCAN, sep="\t", index=False, na_rep="NA")

    # Correlations in the main chromatographic region, using scans where
    # both selected raw centroids are observed.
    main = result[
        result["rt_seconds"].between(*MAIN_RT)
    ].copy()
    corr_rows = []
    for i in range(1, 5):
        for j in range(i + 1, 5):
            x = main[f"peak{i}_intensity"]
            y = main[f"peak{j}_intensity"]
            ok = x.notna() & y.notna()
            corr_rows.append({
                "region": f"{MAIN_RT[0]}-{MAIN_RT[1]}s",
                "peak_pair": f"{i}-{j}",
                "n_scans": int(ok.sum()),
                "pearson": pearson(x, y),
                "spearman": spearman(x, y),
            })
    pd.DataFrame(corr_rows).to_csv(
        OUT_CORR, sep="\t", index=False, na_rep="NA"
    )

    summary_rows = []
    for name, frame in [
        ("main_region", main),
        ("wide_region", result[
            result["rt_seconds"].between(*WIDE_RT)
        ]),
    ]:
        n = len(frame)
        for metric, mask in [
            ("four_peaks_observed", frame["n_peaks_found"] == 4),
            ("all_three_spacings_observed", frame["all_3_spacings_observed"]),
            ("all_spacings_within_3mDa", frame["all_spacings_within_3mDa"]),
            ("all_spacings_within_5mDa", frame["all_spacings_within_5mDa"]),
        ]:
            count = int(mask.sum())
            summary_rows.append({
                "region": name,
                "rt_start_s": frame["rt_seconds"].min() if n else np.nan,
                "rt_end_s": frame["rt_seconds"].max() if n else np.nan,
                "metric": metric,
                "n_scans": n,
                "count": count,
                "fraction": count / n if n else np.nan,
            })

        for j in range(1, 4):
            vals = frame[f"spacing{j}_residual_mDa"].dropna()
            summary_rows.append({
                "region": name,
                "rt_start_s": frame["rt_seconds"].min() if n else np.nan,
                "rt_end_s": frame["rt_seconds"].max() if n else np.nan,
                "metric": f"spacing{j}_residual_mDa",
                "n_scans": n,
                "count": len(vals),
                "fraction": np.nan,
                "median_residual_mDa": vals.median() if len(vals) else np.nan,
                "median_abs_residual_mDa": (
                    vals.abs().median() if len(vals) else np.nan
                ),
                "max_abs_residual_mDa": (
                    vals.abs().max() if len(vals) else np.nan
                ),
            })

    pd.DataFrame(summary_rows).to_csv(
        OUT_SUMMARY, sep="\t", index=False, na_rep="NA"
    )

    print(f"Scans analysed: {len(result)}")
    print(f"Main region: {MAIN_RT[0]}-{MAIN_RT[1]} s, {len(main)} scans")
    print("\nMain-region counts:")
    for col in [
        "n_peaks_found",
        "all_3_spacings_observed",
        "all_spacings_within_3mDa",
        "all_spacings_within_5mDa",
    ]:
        if col == "n_peaks_found":
            print(result.loc[
                result.rt_seconds.between(*MAIN_RT), col
            ].value_counts().sort_index().to_string())
        else:
            print(f"{col}: {int(main[col].sum())}/{len(main)}")

    print("\nIntensity correlations:")
    print(pd.DataFrame(corr_rows).to_string(index=False))
    print("\nOutputs:")
    print(f"  {OUT_SCAN}")
    print(f"  {OUT_SUMMARY}")
    print(f"  {OUT_CORR}")
    print("\nNote: matching uses fixed m/z windows and nearest centroids; "
          "this is not an isotope-envelope identity proof.")


if __name__ == "__main__":
    main()
