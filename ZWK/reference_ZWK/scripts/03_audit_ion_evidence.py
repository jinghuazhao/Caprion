from pathlib import Path
import pandas as pd

project = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
file = project / "search_98/combined_ion.tsv"
out = project / "reports/ion_evidence_audit.txt"

df = pd.read_csv(file, sep="\t", low_memory=False)

# Run-specific evidence columns
match_cols = [c for c in df.columns if c.endswith(" Match Type")]
rt_cols = [c for c in df.columns if c.endswith(" Apex Retention Time")]
intensity_cols = [c for c in df.columns if c.endswith(" Intensity")]
spectral_cols = [c for c in df.columns if c.endswith(" Spectral Count")]

lines = [
    f"Input: {file}",
    f"Rows (ion-level entries): {len(df):,}",
    f"Columns: {len(df.columns)}",
    f"Runs with Match Type columns: {len(match_cols)}",
    f"Runs with Apex Retention Time columns: {len(rt_cols)}",
    f"Runs with Intensity columns: {len(intensity_cols)}",
    f"Runs with Spectral Count columns: {len(spectral_cols)}",
    "",
    "Match Type values across all runs:",
]

if match_cols:
    counts = (
        df[match_cols]
        .melt(value_name="Match Type")["Match Type"]
        .dropna()
        .astype(str)
        .value_counts()
    )
    lines.extend(f"  {k}: {v:,}" for k, v in counts.items())
else:
    lines.append("  No Match Type columns found")

lines += [
    "",
    "Ion identity fields:",
]
for col in [
    "Peptide Sequence", "Modified Sequence", "M/Z", "Charge",
    "Assigned Modifications", "Protein", "Protein ID",
]:
    if col in df.columns:
        lines.append(f"  {col}: present")

lines += [
    "",
    "Rows with per-run evidence:",
]
for label, cols in [
    ("Intensity", intensity_cols),
    ("Apex Retention Time", rt_cols),
    ("Spectral Count", spectral_cols),
]:
    if cols:
        values = df[cols].apply(pd.to_numeric, errors="coerce")
        lines.append(
            f"  {label}: {int(values.notna().any(axis=1).sum()):,} ions "
            f"have at least one non-missing value"
        )

out.write_text("\n".join(lines) + "\n")
print(out.read_text())
