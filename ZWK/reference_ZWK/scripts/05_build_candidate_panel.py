from pathlib import Path
import pandas as pd

project = Path.home() / "Caprion/analysis/ZWK/reference_ZWK"
source = project / "search_98/combined_ion.tsv"
out_dir = project / "reports"

df = pd.read_csv(source, sep="\t", low_memory=False)

identity = [
    "Peptide Sequence", "Modified Sequence", "M/Z", "Charge",
    "Assigned Modifications", "Protein", "Protein ID",
    "Entry Name", "Gene", "Protein Description",
]
identity = [c for c in identity if c in df.columns]

match_cols = [c for c in df.columns if c.endswith(" Match Type")]
match = df[match_cols].astype("string")

df["Runs with direct MS/MS"] = match.eq("MS/MS").sum(axis=1)
df["Runs with MBR"] = match.eq("MBR").sum(axis=1)
df["Runs with either match"] = match.isin(["MS/MS", "MBR"]).sum(axis=1)

# Candidate panel: direct MS/MS evidence in at least two runs.
panel = df.loc[df["Runs with direct MS/MS"] >= 2].copy()

panel_cols = identity + [
    "Runs with direct MS/MS",
    "Runs with MBR",
    "Runs with either match",
]
panel_cols = [c for c in panel_cols if c in panel.columns]

output = out_dir / "candidate_reference_ion_panel.tsv"
panel[panel_cols].to_csv(output, sep="\t", index=False)

print(f"All ion entries: {len(df):,}")
print(f"Candidate panel (direct MS/MS in >=2 runs): {len(panel):,}")
print(f"Written: {output}")
