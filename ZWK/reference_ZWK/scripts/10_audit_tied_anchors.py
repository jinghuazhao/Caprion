#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path("openms_pilot")

BEST_FILE = BASE / "monoisotopic_anchor_best_per_candidate_fast.tsv"
SUMMARY_FILE = BASE / "monoisotopic_anchor_comparison_summary_fast.tsv"
AUDIT_FILE = BASE / "anchor_agreement_audit.tsv"

OUT_TIES = BASE / "tied_anchor_detail.tsv"
OUT_SUMMARY = BASE / "tied_anchor_summary.tsv"

KEYS = [
    "acquisition",
    "Spectrum",
    "feature_id",
    "charge",
    "nominated_isotope_offset",
]


def require_columns(df, columns, name):
    missing = sorted(set(columns) - set(df.columns))
    if missing:
        raise ValueError(f"{name} is missing columns: {missing}")


best = pd.read_csv(BEST_FILE, sep="\t", low_memory=False)
summary = pd.read_csv(SUMMARY_FILE, sep="\t", low_memory=False)
audit = pd.read_csv(AUDIT_FILE, sep="\t", low_memory=False)

require_columns(
    best,
    KEYS + [
        "best_anchor_offset",
        "best_anchor_score",
        "best_mean_positions_found",
        "best_mean_consecutive_isotopes",
        "scans_examined",
        "scans_with_3_consecutive_isotopes",
        "median_abs_spacing_error_ppm",
    ],
    "Best-anchor table",
)

require_columns(
    summary,
    KEYS + [
        "anchor_offset",
        "mean_anchor_score",
        "median_anchor_score",
        "mean_isotope_positions_found",
        "mean_consecutive_isotopes",
        "scans_with_3_consecutive_isotopes",
        "median_abs_spacing_error_ppm",
    ],
    "Anchor summary",
)

require_columns(
    audit,
    KEYS + [
        "n_best_offsets",
        "best_offsets",
        "best_anchor_score",
        "nominated_mean_score",
    ],
    "Candidate audit",
)

# Restrict to candidates explicitly identified as having tied best offsets.
tied_keys = audit.loc[
    audit["n_best_offsets"] > 1,
    KEYS + ["n_best_offsets", "best_offsets"],
].drop_duplicates(KEYS)

if tied_keys.empty:
    print("No tied candidates found.")
    raise SystemExit(0)

# Keep only rows for the tied candidates and their tied best offsets.
tied = best.merge(
    tied_keys,
    on=KEYS,
    how="inner",
    validate="many_to_one",
)

# Add the per-offset summary metrics. There should be one summary row per
# candidate and anchor offset; fail explicitly if this assumption is false.
if summary.duplicated(KEYS + ["anchor_offset"]).any():
    raise ValueError(
        "Anchor summary contains duplicate candidate/anchor-offset keys. "
        "Inspect the summary before interpreting this audit."
    )

tied = tied.merge(
    summary,
    left_on=KEYS + ["best_anchor_offset"],
    right_on=KEYS + ["anchor_offset"],
    how="left",
    validate="many_to_one",
    suffixes=("_best", "_summary"),
)

# Compare each tied offset's spacing error with the minimum spacing error
# among that candidate's tied offsets. This is a descriptive comparison,
# not a probability or statistical significance test.
err_col = "median_abs_spacing_error_ppm_best"
if err_col not in tied.columns:
    candidates = [
        c for c in tied.columns
        if "median_abs_spacing_error_ppm" in c
    ]
    raise ValueError(
        f"Could not identify spacing-error column {err_col!r}. "
        f"Available candidates: {candidates}"
    )

tied["minimum_tied_spacing_error_ppm"] = tied.groupby(
    KEYS, dropna=False
)[err_col].transform("min")

tied["spacing_error_excess_ppm"] = (
    tied[err_col] - tied["minimum_tied_spacing_error_ppm"]
)

tied["spacing_error_difference_ppm"] = tied.groupby(
    KEYS, dropna=False
)[err_col].transform("max") - tied.groupby(
    KEYS, dropna=False
)[err_col].transform("min")

# Rank the evidence for review: lower spacing error is preferable, followed
# by greater isotope-position and consecutive-isotope support.
tied = tied.sort_values(
    KEYS + [
        err_col,
        "best_mean_consecutive_isotopes",
        "best_mean_positions_found",
        "best_anchor_offset",
    ],
    ascending=[
        True, True, True, True, True,
        True, False, False, True,
    ],
    na_position="last",
)

tied.to_csv(OUT_TIES, sep="\t", index=False)

# Candidate-level summary of how much the tied offsets differ.
candidate_summary = tied.groupby(
    KEYS, dropna=False
).agg(
    n_tied_rows=("best_anchor_offset", "size"),
    n_distinct_tied_offsets=("best_anchor_offset", "nunique"),
    tied_offsets=(
        "best_anchor_offset",
        lambda x: ",".join(map(str, sorted(set(x.dropna())))),
    ),
    best_score_range=(
        "best_anchor_score",
        lambda x: x.max() - x.min(),
    ),
    spacing_error_range_ppm=(
        "spacing_error_difference_ppm",
        "max",
    ),
    minimum_spacing_error_ppm=(
        "minimum_tied_spacing_error_ppm",
        "min",
    ),
    maximum_consecutive_isotopes=(
        "best_mean_consecutive_isotopes",
        "max",
    ),
    maximum_positions_found=(
        "best_mean_positions_found",
        "max",
    ),
).reset_index()

candidate_summary["spacing_error_tie_class"] = pd.cut(
    candidate_summary["spacing_error_range_ppm"],
    bins=[-float("inf"), 1, 10, 100, float("inf")],
    labels=[
        "range_le_1_ppm",
        "range_1_to_10_ppm",
        "range_10_to_100_ppm",
        "range_gt_100_ppm",
    ],
)

candidate_summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

print(f"Tied-offset detail: {OUT_TIES}")
print(f"Tied-offset summary: {OUT_SUMMARY}")
print()
print("Tied candidates:", len(candidate_summary))
print("\nTied offsets per candidate:")
print(
    candidate_summary["n_distinct_tied_offsets"]
    .value_counts()
    .sort_index()
    .to_string()
)
print("\nSpacing-error range classes:")
print(
    candidate_summary["spacing_error_tie_class"]
    .value_counts(dropna=False)
    .to_string()
)
print("\nLargest spacing-error differences:")
print(
    candidate_summary.sort_values(
        "spacing_error_range_ppm", ascending=False
    ).head(15).to_string(index=False)
)
