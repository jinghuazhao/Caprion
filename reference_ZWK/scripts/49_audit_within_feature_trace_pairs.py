#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
data = root / "openms_pilot"

trace_files = {
    "first_pilot": data / "isotope_trace_peaks.tsv",
    "second_acquisition": data / "isotope_trace_peaks_second.tsv",
}

out_detail = data / "candidate_isotope_group_within_feature_pair_detail_v2.tsv"
out_summary = data / "candidate_isotope_group_within_feature_pair_summary_v2.tsv"

target_features = {
    ("first_pilot", "10275997325121874585"),
    ("second_acquisition", "12167939761234023158"),
}

C13 = 1.00335483507
MZ_TOLERANCE = 2.0  # mDa, diagnostic threshold only

groups_path = data / "candidate_isotope_groups_v2.tsv"
groups = pd.read_csv(groups_path, sep="\t", low_memory=False)

required_groups = {
    "acquisition", "candidate_ig_id", "feature_id",
    "charge", "trace_index"
}
missing = required_groups - set(groups.columns)
if missing:
    raise ValueError(f"Candidate group table lacks columns: {sorted(missing)}")

groups["feature_id"] = groups["feature_id"].astype(str)
groups["trace_index"] = pd.to_numeric(groups["trace_index"], errors="coerce")
groups["charge"] = pd.to_numeric(groups["charge"], errors="coerce")

membership = {}
for _, row in groups.iterrows():
    key = (
        str(row["acquisition"]),
        str(row["feature_id"]),
        int(row["charge"]) if pd.notna(row["charge"]) else None,
        int(row["trace_index"]) if pd.notna(row["trace_index"]) else None,
    )
    membership.setdefault(key, []).append(str(row["candidate_ig_id"]))

detail_rows = []
summary_rows = []

for acquisition, path in trace_files.items():
    if not path.exists():
        print(f"Skipping missing file: {path}", flush=True)
        continue

    print(f"Reading {path.name}...", flush=True)
    raw = pd.read_csv(path, sep="\t", low_memory=False)

    required_raw = {
        "feature_id", "feature_charge", "trace_index",
        "scan_index", "rt_sec", "observed_mz", "peak_intensity"
    }
    missing = required_raw - set(raw.columns)
    if missing:
        raise ValueError(f"{path.name} lacks columns: {sorted(missing)}")

    raw["feature_id"] = raw["feature_id"].astype(str)
    raw["feature_charge"] = pd.to_numeric(raw["feature_charge"], errors="coerce")
    raw["trace_index"] = pd.to_numeric(raw["trace_index"], errors="coerce")
    raw["scan_index"] = pd.to_numeric(raw["scan_index"], errors="coerce")
    raw["rt_sec"] = pd.to_numeric(raw["rt_sec"], errors="coerce")
    raw["observed_mz"] = pd.to_numeric(raw["observed_mz"], errors="coerce")
    raw["peak_intensity"] = pd.to_numeric(raw["peak_intensity"], errors="coerce")

    target_ids = {
        feature_id for acq, feature_id in target_features
        if acq == acquisition
    }
    raw = raw.loc[raw["feature_id"].isin(target_ids)].copy()

    for feature_id, feature in raw.groupby("feature_id", sort=False):
        charges = feature["feature_charge"].dropna().unique()
        if len(charges) != 1:
            print(
                f"WARNING: {acquisition}/{feature_id}: "
                f"found charges {charges.tolist()}",
                flush=True,
            )
            continue

        charge = int(charges[0])
        expected_spacing = C13 / charge

        # A trace can contain more than one row for a scan. Collapse to one
        # intensity and one median m/z per trace/scan before pairing.
        scan_trace = (
            feature.dropna(subset=["trace_index", "scan_index"])
            .groupby(["trace_index", "scan_index"], as_index=False)
            .agg(
                rt_sec=("rt_sec", "median"),
                observed_mz=("observed_mz", "median"),
                peak_intensity=("peak_intensity", "sum"),
            )
        )

        trace_indices = sorted(
            int(x) for x in scan_trace["trace_index"].dropna().unique()
        )
        print(
            f"  {acquisition}/{feature_id}, z={charge}: "
            f"{len(trace_indices)} traces",
            flush=True,
        )

        traces = {}
        for trace_index in trace_indices:
            t = scan_trace.loc[
                scan_trace["trace_index"] == trace_index
            ].copy()
            t = t.drop_duplicates("scan_index").set_index("scan_index").sort_index()
            traces[trace_index] = t

        for trace_a, trace_b in combinations(trace_indices, 2):
            a = traces[trace_a]
            b = traces[trace_b]
            common = a.index.intersection(b.index)

            key_a = (acquisition, str(feature_id), charge, trace_a)
            key_b = (acquisition, str(feature_id), charge, trace_b)
            ig_a = ";".join(sorted(set(membership.get(key_a, []))))
            ig_b = ";".join(sorted(set(membership.get(key_b, []))))
            same_ig = bool(ig_a and ig_b and set(ig_a.split(";")) & set(ig_b.split(";")))

            row = {
                "acquisition": acquisition,
                "feature_id": str(feature_id),
                "charge": charge,
                "trace_A": trace_a,
                "trace_B": trace_b,
                "candidate_IG_A": ig_a,
                "candidate_IG_B": ig_b,
                "same_candidate_IG": same_ig,
                "n_scans_A": len(a),
                "n_scans_B": len(b),
                "n_common_scans": len(common),
                "expected_isotope_spacing_Da": expected_spacing,
            }

            if len(common) == 0:
                row.update({
                    "median_mz_delta_mDa": np.nan,
                    "median_abs_mz_delta_MAD_mDa": np.nan,
                    "mz_delta_range_mDa": np.nan,
                    "nearest_isotope_step": np.nan,
                    "median_spacing_residual_mDa": np.nan,
                    "fraction_spacing_residual_le_2mDa": np.nan,
                    "intensity_correlation": np.nan,
                    "rt_delta_median_sec": np.nan,
                })
            else:
                aa = a.loc[common]
                bb = b.loc[common]

                # Delta is B minus A. Use the median m/z at each scan.
                delta = (
                    bb["observed_mz"].to_numpy()
                    - aa["observed_mz"].to_numpy()
                )
                valid = np.isfinite(delta)
                delta = delta[valid]

                if len(delta):
                    median_delta = float(np.median(delta))
                    delta_mad = float(np.median(np.abs(delta - median_delta)))
                    delta_range = float(np.max(delta) - np.min(delta))

                    # Isotope spacing is charge-aware. Signed integer steps
                    # allow either trace to have the lower m/z.
                    steps = np.rint(delta / expected_spacing)
                    residuals = delta - steps * expected_spacing

                    median_step = int(np.rint(median_delta / expected_spacing))
                    median_residual = float(
                        np.median(np.abs(residuals))
                    )
                    fraction_residual = float(
                        np.mean(np.abs(residuals) <= MZ_TOLERANCE / 1000.0)
                    )

                    row.update({
                        "median_mz_delta_mDa": median_delta * 1000.0,
                        "median_abs_mz_delta_MAD_mDa": delta_mad * 1000.0,
                        "mz_delta_range_mDa": delta_range * 1000.0,
                        "nearest_isotope_step": median_step,
                        "median_spacing_residual_mDa": median_residual * 1000.0,
                        "fraction_spacing_residual_le_2mDa": fraction_residual,
                    })
                else:
                    row.update({
                        "median_mz_delta_mDa": np.nan,
                        "median_abs_mz_delta_MAD_mDa": np.nan,
                        "mz_delta_range_mDa": np.nan,
                        "nearest_isotope_step": np.nan,
                        "median_spacing_residual_mDa": np.nan,
                        "fraction_spacing_residual_le_2mDa": np.nan,
                    })

                ia = aa["peak_intensity"].to_numpy(dtype=float)
                ib = bb["peak_intensity"].to_numpy(dtype=float)
                good = np.isfinite(ia) & np.isfinite(ib)
                corr = np.nan
                if good.sum() >= 3 and np.std(ia[good]) > 0 and np.std(ib[good]) > 0:
                    corr = float(np.corrcoef(ia[good], ib[good])[0, 1])

                rt_delta = (
                    bb["rt_sec"].to_numpy(dtype=float)
                    - aa["rt_sec"].to_numpy(dtype=float)
                )
                rt_delta = rt_delta[np.isfinite(rt_delta)]

                row.update({
                    "intensity_correlation": corr,
                    "rt_delta_median_sec": (
                        float(np.median(rt_delta)) if len(rt_delta) else np.nan
                    ),
                })

            detail_rows.append(row)

            # Summary for each pair, with all available scan-level metrics.
            summary_rows.append(row.copy())

detail = pd.DataFrame(detail_rows)
summary = pd.DataFrame(summary_rows)

if detail.empty:
    raise RuntimeError("No within-feature trace pairs were generated.")

detail = detail.sort_values(
    ["acquisition", "feature_id", "trace_A", "trace_B"]
)
summary = summary.sort_values(
    ["acquisition", "feature_id", "n_common_scans"],
    ascending=[True, True, False],
)

detail.to_csv(out_detail, sep="\t", index=False)
summary.to_csv(out_summary, sep="\t", index=False)

print(f"\nPair detail: {out_detail}")
print(f"Pair summary: {out_summary}")
print(f"Total trace pairs: {len(detail)}")

display_cols = [
    "acquisition", "feature_id", "charge",
    "trace_A", "trace_B", "candidate_IG_A", "candidate_IG_B",
    "same_candidate_IG", "n_common_scans",
    "median_mz_delta_mDa", "median_abs_mz_delta_MAD_mDa",
    "nearest_isotope_step", "median_spacing_residual_mDa",
    "fraction_spacing_residual_le_2mDa", "intensity_correlation",
]
display_cols = [c for c in display_cols if c in summary.columns]

print("\n=== Within-feature trace-pair summary ===")
print(summary[display_cols].to_string(index=False))
