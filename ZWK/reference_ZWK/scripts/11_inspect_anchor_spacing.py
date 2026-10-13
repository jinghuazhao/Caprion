#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

base = Path("openms_pilot")
summary_file = base / "monoisotopic_anchor_comparison_summary_fast.tsv"
ties_file = base / "tied_anchor_summary.tsv"
detail_file = base / "tied_anchor_detail.tsv"

summary = pd.read_csv(summary_file, sep="\t", low_memory=False)
ties = pd.read_csv(ties_file, sep="\t", low_memory=False)
detail = pd.read_csv(detail_file, sep="\t", low_memory=False)

print("=== Summary columns ===")
print(summary.columns.tolist())

print("\n=== Spacing-error summary ===")
for col in summary.columns:
    if "spacing" in col.lower() or "ppm" in col.lower():
        print(f"\n{col}")
        print(summary[col].describe().to_string())

print("\n=== Example summary rows ===")
print(summary.head(5).to_string(index=False))

# Inspect candidates with the largest reported spacing-error ranges.
largest = ties.dropna(
    subset=["spacing_error_range_ppm"]
).nlargest(10, "spacing_error_range_ppm")

keys = [
    "acquisition",
    "Spectrum",
    "feature_id",
    "charge",
    "nominated_isotope_offset",
]

print("\n=== Largest tied candidates: offset-level evidence ===")
for _, row in largest.iterrows():
    mask = pd.Series(True, index=detail.index)
    for key in keys:
        mask &= detail[key].eq(row[key])

    subset = detail.loc[mask].copy()
    print("\nCandidate:")
    print({key: row[key] for key in keys})
    print(
        subset.to_string(index=False)
        if len(subset)
        else "No matching detail rows found."
    )

print("\n=== Units reminder ===")
print(
    "Check the exact formula used to calculate the spacing-error column "
    "before treating its values as ppm or using them to rank anchors."
)
