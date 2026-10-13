from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
BEST_FILE = BASE / "monoisotopic_anchor_best_per_candidate_fast.tsv"
SUMMARY_FILE = BASE / "monoisotopic_anchor_comparison_summary_fast.tsv"
SCAN_FILE = BASE / "monoisotopic_anchor_comparison_per_scan_fast.tsv"

OUT_CANDIDATES = BASE / "anchor_agreement_audit.tsv"
OUT_GROUPS = BASE / "anchor_agreement_group_summary.tsv"

keys = [
    "acquisition", "Spectrum", "feature_id",
    "charge", "nominated_isotope_offset",
]

best = pd.read_csv(BEST_FILE, sep="\t", low_memory=False)
summary = pd.read_csv(SUMMARY_FILE, sep="\t", low_memory=False)

# Build a unique candidate-level table while retaining the number of tied
# best offsets and all tied offset values.
best_offsets = best.groupby(keys, dropna=False).agg(
    n_best_offsets=("best_anchor_offset", "nunique"),
    best_offsets=(
        "best_anchor_offset",
        lambda x: ",".join(map(str, sorted(set(x)))),
    ),
    best_anchor_score=("best_anchor_score", "max"),
    best_mean_positions_found=("best_mean_positions_found", "max"),
    best_mean_consecutive_isotopes=(
        "best_mean_consecutive_isotopes", "max"
    ),
    scans_examined=("scans_examined", "max"),
).reset_index()

# A tie is not silently resolved. "Best offset" is only unique when one
# distinct anchor offset was retained.
best_offsets["best_anchor_offset_unique"] = (
    best_offsets["n_best_offsets"] == 1
)

unique_best = best.drop_duplicates(keys).copy()
unique_best = unique_best[keys + [
    "best_anchor_offset",
]]
unique_best = unique_best.rename(columns={
    "best_anchor_offset": "representative_best_anchor_offset"
})

# Nominated-anchor metrics. Do not rename or duplicate the key column.
nominated = summary[
    summary["anchor_offset"] == summary["nominated_isotope_offset"]
].copy()

nominated = nominated[keys + [
    "mean_anchor_score",
    "median_anchor_score",
    "mean_isotope_positions_found",
    "median_isotope_positions_found",
    "mean_consecutive_isotopes",
    "median_consecutive_isotopes",
    "scans_with_3_consecutive_isotopes",
    "median_abs_spacing_error_ppm",
]].rename(columns={
    "mean_anchor_score": "nominated_mean_score",
    "median_anchor_score": "nominated_median_score",
    "mean_isotope_positions_found": "nominated_mean_positions",
    "median_isotope_positions_found": "nominated_median_positions",
    "mean_consecutive_isotopes": "nominated_mean_consecutive",
    "median_consecutive_isotopes": "nominated_median_consecutive",
    "scans_with_3_consecutive_isotopes": "nominated_scans_with_3_consecutive",
    "median_abs_spacing_error_ppm": "nominated_spacing_error_ppm",
})

audit = best_offsets.merge(
    unique_best, on=keys, how="left", validate="one_to_one"
)
audit = audit.merge(
    nominated, on=keys, how="left", validate="one_to_one"
)

audit["anchor_agrees_unique"] = (
    audit["best_anchor_offset_unique"]
    & (audit["representative_best_anchor_offset"]
       == audit["nominated_isotope_offset"])
)

audit["score_advantage_best_minus_nominated"] = (
    audit["best_anchor_score"] - audit["nominated_mean_score"]
)

audit["best_series_recurrent"] = (
    (audit["best_mean_consecutive_isotopes"] >= 3)
    & (audit["scans_examined"] >= 3)
)
audit["nominated_series_recurrent"] = (
    (audit["nominated_mean_consecutive"] >= 3)
    & (audit["nominated_scans_with_3_consecutive"] >= 3)
)

audit["interpretation_group"] = np.select(
    [
        ~audit["best_anchor_offset_unique"],
        audit["anchor_agrees_unique"] & audit["best_series_recurrent"],
        audit["anchor_agrees_unique"],
        ~audit["anchor_agrees_unique"] & audit["best_series_recurrent"]
            & ~audit["nominated_series_recurrent"],
        ~audit["anchor_agrees_unique"] & audit["best_series_recurrent"]
            & audit["nominated_series_recurrent"],
    ],
    [
        "tied_best_anchor",
        "agreement_recurrent_series",
        "agreement_without_recurrent_series",
        "disagreement_best_series_only",
        "disagreement_both_series_recurrent",
    ],
    default="disagreement_weak_or_ambiguous",
)

audit.to_csv(OUT_CANDIDATES, sep="\t", index=False)

group_summary = audit.groupby(
    ["acquisition", "interpretation_group"], dropna=False
).agg(
    candidates=("Spectrum", "size"),
    median_best_score=("best_anchor_score", "median"),
    median_score_advantage=(
        "score_advantage_best_minus_nominated", "median"
    ),
    fraction_best_recurrent=("best_series_recurrent", "mean"),
    fraction_nominated_recurrent=("nominated_series_recurrent", "mean"),
).reset_index()

group_summary.to_csv(OUT_GROUPS, sep="\t", index=False)

print(f"Candidate audit: {OUT_CANDIDATES}")
print(f"Group summary:   {OUT_GROUPS}")
print("\nCandidate counts:")
print("Best-table rows:", len(best))
print("Unique candidate keys:", len(audit))
print("Candidates with tied best offsets:",
      int((audit["n_best_offsets"] > 1).sum()))
print("\nInterpretation groups:")
print(audit.groupby(
    ["acquisition", "interpretation_group"], dropna=False
).size().to_string())
print("\nAgreement among candidates with a unique best offset:")
unique = audit[audit["best_anchor_offset_unique"]]
print(
    f"{int(unique['anchor_agrees_unique'].sum())} / {len(unique)} "
    f"({unique['anchor_agrees_unique'].mean():.1%})"
    if len(unique) else "No unique best offsets."
)
