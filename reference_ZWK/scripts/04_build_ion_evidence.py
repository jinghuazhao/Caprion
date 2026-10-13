from pathlib import Path
import pandas as pd

project = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
source = project / "search_98/combined_ion.tsv"
out_dir = project / "reports"
out_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(source, sep="\t", low_memory=False)

identity_cols = [
    "Peptide Sequence",
    "Modified Sequence",
    "M/Z",
    "Charge",
    "Assigned Modifications",
    "Protein",
    "Protein ID",
    "Entry Name",
    "Gene",
    "Protein Description",
]
identity_cols = [c for c in identity_cols if c in df.columns]

# Build a long-format table from the per-run evidence columns.
runs = sorted({
    c.removesuffix(" Match Type").strip()
    for c in df.columns
    if c.endswith(" Match Type")
})

parts = []
for run in runs:
    columns = {
        "Match Type": f"{run} Match Type",
        "Apex Retention Time": f"{run} Apex Retention Time",
        "Intensity": f"{run} Intensity",
        "Spectral Count": f"{run} Spectral Count",
    }
    present = {k: v for k, v in columns.items() if v in df.columns}
    if "Match Type" not in present:
        continue

    part = df[identity_cols].copy()
    part["Run"] = run
    for label, column in present.items():
        part[label] = df[column]
    parts.append(part)

long = pd.concat(parts, ignore_index=True)

# Retain rows with an actual match or quantitative evidence.
has_match = long["Match Type"].notna() & long["Match Type"].ne("unmatched")
has_intensity = (
    pd.to_numeric(long["Intensity"], errors="coerce").notna()
    if "Intensity" in long.columns
    else False
)
has_rt = (
    pd.to_numeric(long["Apex Retention Time"], errors="coerce").notna()
    if "Apex Retention Time" in long.columns
    else False
)
long = long.loc[has_match | has_intensity | has_rt].copy()

evidence_file = out_dir / "ion_evidence_long.tsv"
long.to_csv(evidence_file, sep="\t", index=False)

# One row per ion, summarising the distinction between direct and transferred evidence.
summary = (
    long.assign(
        is_msms=long["Match Type"].eq("MS/MS"),
        is_mbr=long["Match Type"].eq("MBR"),
    )
    .groupby(identity_cols, dropna=False)
    .agg(
        Runs_with_MS_MS=("is_msms", "sum"),
        Runs_with_MBR=("is_mbr", "sum"),
    )
    .reset_index()
)

summary_file = out_dir / "ion_evidence_summary.tsv"
summary.to_csv(summary_file, sep="\t", index=False)

print(f"Runs represented: {len(runs)}")
print(f"Long-format evidence rows: {len(long):,}")
print(f"Unique ion entries: {len(summary):,}")
print(f"Wrote: {evidence_file}")
print(f"Wrote: {summary_file}")
