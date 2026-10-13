#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
data = root / "openms_pilot"

src = data / "candidate_isotope_group_trace_pair_stability_summary_v2.tsv"
detail_src = data / "candidate_isotope_group_trace_pair_stability_detail_v2.tsv"

summary_out = data / "candidate_isotope_group_charge_aware_spacing_summary_v2.tsv"
detail_out = data / "candidate_isotope_group_charge_aware_spacing_detail_v2.tsv"

print(f"Reading summary: {src}", flush=True)
if not src.exists():
    raise FileNotFoundError(src)

summary = pd.read_csv(src, sep="\t", low_memory=False)

# Recover charge from the corresponding pair ID, e.g. "...:z2:...".
def charge_from_pair_id(value):
    import re
    m = re.search(r":z(\d+):", str(value))
    return int(m.group(1)) if m else np.nan

summary["charge"] = summary["pair_id"].map(charge_from_pair_id)
if summary["charge"].isna().any():
    bad = summary.loc[summary["charge"].isna(), "pair_id"].unique()
    raise ValueError(f"Could not infer charge from pair_id: {bad.tolist()}")

summary["charge"] = summary["charge"].astype(int)
summary["expected_adjacent_isotope_spacing_Da"] = (
    1.00335483507 / summary["charge"]
)
summary["expected_adjacent_isotope_spacing_mDa"] = (
    1000.0 * summary["expected_adjacent_isotope_spacing_Da"]
)

delta = summary["median_delta_mDa"].to_numpy(dtype=float)
spacing = summary["expected_adjacent_isotope_spacing_mDa"].to_numpy(dtype=float)

# Test the nearest integer number of charge-adjusted isotope spacings.
summary["nearest_isotope_step_count"] = np.rint(delta / spacing).astype(int)
summary["residual_to_nearest_charge_aware_spacing_mDa"] = (
    delta - summary["nearest_isotope_step_count"] * spacing
)
summary["absolute_spacing_residual_mDa"] = (
    summary["residual_to_nearest_charge_aware_spacing_mDa"].abs()
)

# These flags are descriptive screening criteria, not proof of a common ladder.
summary["spacing_residual_le_2mDa"] = (
    summary["absolute_spacing_residual_mDa"] <= 2.0
)
summary["stable_median_screen"] = (
    (summary["n_scans"] >= 5)
    & (summary["delta_MAD_mDa"] <= 2.0)
    & (summary["fraction_within_2mDa_of_median"] >= 0.70)
)
summary["passes_both_screening_rules"] = (
    summary["spacing_residual_le_2mDa"]
    & summary["stable_median_screen"]
)
summary["interpretation"] = np.select(
    [
        summary["passes_both_screening_rules"],
        summary["spacing_residual_le_2mDa"],
        summary["stable_median_screen"],
    ],
    [
        "charge_spacing_and_stability_screen_pass",
        "spacing_only_pass",
        "stability_only_pass",
    ],
    default="neither_screen_pass",
)

summary.to_csv(summary_out, sep="\t", index=False)

# Add the same charge-aware spacing test to each individual scan row.
if detail_src.exists():
    detail = pd.read_csv(detail_src, sep="\t", low_memory=False)
    detail["charge"] = detail["pair_id"].map(charge_from_pair_id)
    if detail["charge"].isna().any():
        raise ValueError("Could not infer charge for every detail row.")

    detail["charge"] = detail["charge"].astype(int)
    detail["expected_adjacent_isotope_spacing_mDa"] = (
        1000.0 * 1.00335483507 / detail["charge"]
    )
    detail["nearest_isotope_step_count"] = np.rint(
        detail["mz_delta_mDa"]
        / detail["expected_adjacent_isotope_spacing_mDa"]
    ).astype(int)
    detail["residual_to_nearest_charge_aware_spacing_mDa"] = (
        detail["mz_delta_mDa"]
        - detail["nearest_isotope_step_count"]
        * detail["expected_adjacent_isotope_spacing_mDa"]
    )
    detail["absolute_spacing_residual_mDa"] = (
        detail["residual_to_nearest_charge_aware_spacing_mDa"].abs()
    )
    detail.to_csv(detail_out, sep="\t", index=False)
    print(f"Scan-level detail: {detail_out}")

print(f"Summary output: {summary_out}")
print(f"Trace-pair rows: {len(summary)}")
print("Screen passing both criteria:",
      int(summary["passes_both_screening_rules"].sum()))
print("\n=== Charge-aware isotope-spacing audit ===")
cols = [
    "feature_id", "charge", "trace_A", "trace_B", "n_scans",
    "median_delta_mDa", "delta_MAD_mDa",
    "expected_adjacent_isotope_spacing_mDa",
    "nearest_isotope_step_count",
    "residual_to_nearest_charge_aware_spacing_mDa",
    "fraction_within_2mDa_of_median",
    "interpretation",
]
print(summary[cols].sort_values(
    ["feature_id", "trace_A", "trace_B"]
).to_string(index=False))
