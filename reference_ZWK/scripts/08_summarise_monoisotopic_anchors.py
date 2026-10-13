from pathlib import Path
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "monoisotopic_anchor_comparison_per_scan_fast.tsv"
SUMMARY = BASE / "monoisotopic_anchor_comparison_summary_fast.tsv"
BEST = BASE / "monoisotopic_anchor_best_per_candidate_fast.tsv"

if not INPUT.exists():
    raise FileNotFoundError(INPUT)

keys = [
    "acquisition",
    "Spectrum",
    "feature_id",
    "charge",
    "nominated_isotope_offset",
    "anchor_offset",
]

# The per-scan file is already available. Summarise it without rerunning
# the feature/scan matching stage.
df = pd.read_csv(INPUT, sep="\t", low_memory=False)

if df.empty:
    raise ValueError(f"No rows in {INPUT}")

required = set(keys) | {
    "anchor_score",
    "isotope_positions_found",
    "consecutive_isotopes_from_anchor",
    "median_abs_adjacent_spacing_error_ppm",
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in {INPUT}: {sorted(missing)}")

print(f"Loaded {len(df):,} scan/anchor rows", flush=True)

summary = df.groupby(keys, dropna=False).agg(
    scans_examined=("anchor_score", "size"),
    mean_anchor_score=("anchor_score", "mean"),
    median_anchor_score=("anchor_score", "median"),
    mean_isotope_positions_found=("isotope_positions_found", "mean"),
    median_isotope_positions_found=("isotope_positions_found", "median"),
    mean_consecutive_isotopes=(
        "consecutive_isotopes_from_anchor", "mean"
    ),
    median_consecutive_isotopes=(
        "consecutive_isotopes_from_anchor", "median"
    ),
    scans_with_3_consecutive_isotopes=(
        "consecutive_isotopes_from_anchor",
        lambda x: int((x >= 3).sum()),
    ),
    median_abs_spacing_error_ppm=(
        "median_abs_adjacent_spacing_error_ppm", "median"
    ),
).reset_index()

summary.to_csv(SUMMARY, sep="\t", index=False)

# Rank candidate anchors by their mean score across scans.
candidate_keys = [
    "acquisition",
    "Spectrum",
    "feature_id",
    "charge",
    "nominated_isotope_offset",
]

summary["anchor_rank"] = summary.groupby(
    candidate_keys, dropna=False
)["mean_anchor_score"].rank(
    method="min", ascending=False
)

best = summary[summary["anchor_rank"] == 1].copy()
best = best.rename(columns={
    "anchor_offset": "best_anchor_offset",
    "mean_anchor_score": "best_anchor_score",
    "mean_isotope_positions_found": "best_mean_positions_found",
    "mean_consecutive_isotopes": "best_mean_consecutive_isotopes",
})

best[candidate_keys + [
    "best_anchor_offset",
    "best_anchor_score",
    "best_mean_positions_found",
    "best_mean_consecutive_isotopes",
    "scans_examined",
    "scans_with_3_consecutive_isotopes",
    "median_abs_spacing_error_ppm",
]].to_csv(BEST, sep="\t", index=False)

print(f"\nSummary: {SUMMARY}")
print(f"Best anchors: {BEST}")

print("\nBest-anchor offsets:")
print(best["best_anchor_offset"].value_counts().sort_index().to_string())

print("\nAgreement with nominated isotope offset:")
# Convention: target_mz = psm_mz - isotope_offset * isotope_spacing / charge.
# Therefore the corresponding anchor_offset should equal isotope_offset.
# The script's anchor definition is psm_mz - anchor_offset * spacing,
# so agreement requires anchor_offset == nominated_isotope_offset.
print(
    "Matching anchor offset:",
    int((best["best_anchor_offset"] ==
         best["nominated_isotope_offset"]).sum()),
    "/",
    len(best),
)
