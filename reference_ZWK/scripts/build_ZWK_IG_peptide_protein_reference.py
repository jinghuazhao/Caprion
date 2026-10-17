
#!/usr/bin/env python3
"""Build an evidence-labelled IG–peptide–protein reference from existing data."""

from pathlib import Path
import json
import re
import pandas as pd

ROOT = Path.cwd()
WORKBOOK = Path.home() / "Caprion/pre_qc_data/pilot/ZWK_EDR_20191002.xlsx"
CANDIDATES = ROOT / "reports/ZWK_stage_B/ZWK_stage_B_IG_PSM_candidates.tsv"
OUT = ROOT / "reports/ZWK_reference"

def sequence_key(value):
    """Conservative sequence normalisation for matching modified peptide strings."""
    if pd.isna(value):
        return ""
    s = str(value).strip().upper()

    # Remove common modification annotations before extracting amino-acid letters.
    s = re.sub(r"\[[^\]]*\]|\([^)]*\)", "", s)

    # Convert conventional flanking notation K.PEPTIDE.R to PEPTIDE.
    parts = s.split(".")
    if len(parts) == 3 and len(parts[0]) <= 1 and len(parts[2]) <= 1:
        s = parts[1]

    return "".join(re.findall(r"[A-Z]", s))

def clean_list(values):
    if hasattr(values, "dropna"):
        values = values.dropna()

    return sorted({
        str(x).strip()
        for x in values
        if pd.notna(x)
        and str(x).strip()
        and str(x).strip().lower() != "nan"
    })

def main():
    OUT.mkdir(parents=True, exist_ok=True)

    for p in (WORKBOOK, CANDIDATES):
        if not p.is_file():
            raise FileNotFoundError(f"Required input not found: {p}")

    print("Reading Stage B candidates...")
    candidates = pd.read_csv(CANDIDATES, sep="\t", low_memory=False)

    required = {
        "sample_id", "reconstructed_ig_id", "charge", "ig_mz",
        "ig_rt_median_seconds", "psm_record_id", "peptide",
        "mz_error_ppm", "rt_error_seconds", "link_status"
    }
    missing = required - set(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table lacks required columns: {sorted(missing)}")

    candidates["peptide_sequence_key"] = candidates["peptide"].map(sequence_key)

    print("Reading workbook Normalized Peptides...")
    workbook = pd.read_excel(
        WORKBOOK, sheet_name="Normalized Peptides", header=4
    )

    required_wb = {
        "Isotope Group ID", "Protein", "Modified Peptide Sequence",
        "Monoisotopic m/z", "Max Isotope Time Centroid", "Charge"
    }
    missing = required_wb - set(workbook.columns)
    if missing:
        raise ValueError(f"Workbook lacks required columns: {sorted(missing)}")

    workbook["peptide_sequence_key"] = (
        workbook["Modified Peptide Sequence"].map(sequence_key)
    )
    workbook = workbook.loc[
        workbook["peptide_sequence_key"].ne("")
        & workbook["Protein"].notna()
    ].copy()

    # One mapping row per normalised sequence and protein.
    mapping = (
        workbook.groupby(
            ["peptide_sequence_key", "Protein"], dropna=False
        )
        .agg(
            workbook_peptide_forms=(
                "Modified Peptide Sequence",
                lambda x: ";".join(clean_list(x))
            ),
            workbook_ig_ids=(
                "Isotope Group ID",
                lambda x: ";".join(clean_list(x))
            ),
            workbook_row_count=("Protein", "size"),
            workbook_mz_values=(
                "Monoisotopic m/z",
                lambda x: ";".join(clean_list(x))
            ),
            workbook_charge_values=(
                "Charge", lambda x: ";".join(clean_list(x))
            ),
            workbook_rt_values=(
                "Max Isotope Time Centroid",
                lambda x: ";".join(clean_list(x))
            ),
        )
        .reset_index()
        .rename(columns={"Protein": "workbook_protein"})
    )

    mapping_path = OUT / "ZWK_peptide_protein_mapping.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False)

    # Summarise all workbook protein mappings per sequence.
    protein_summary = (
        mapping.groupby("peptide_sequence_key")
        .agg(
            workbook_proteins=(
                "workbook_protein",
                lambda x: ";".join(clean_list(x))
            ),
            n_workbook_proteins=("workbook_protein", "nunique"),
            workbook_peptide_forms=(
                "workbook_peptide_forms",
                lambda x: ";".join(clean_list(
                    item for value in x for item in str(value).split(";")
                ))
            ),
            workbook_ig_ids=(
                "workbook_ig_ids",
                lambda x: ";".join(clean_list(
                    item for value in x for item in str(value).split(";")
                ))
            ),
        )
        .reset_index()
    )

    # Retain all Stage B candidate rows. Do not join reconstructed IG IDs
    # directly to workbook IG IDs: their identity has not been established.
    reference = candidates.merge(
        protein_summary,
        on="peptide_sequence_key",
        how="left",
        validate="many_to_one"
    )

    reference["n_workbook_proteins"] = (
        reference["n_workbook_proteins"].fillna(0).astype(int)
    )
    reference["peptide_protein_status"] = "not_found_in_workbook"
    reference.loc[
        reference["n_workbook_proteins"] == 1,
        "peptide_protein_status"
    ] = "workbook_sequence_match_single_protein"
    reference.loc[
        reference["n_workbook_proteins"] > 1,
        "peptide_protein_status"
    ] = "workbook_sequence_match_shared_or_ambiguous"

    reference["evidence_interpretation"] = (
        "Stage B IG–PSM link is a candidate, not a confirmed identification; "
        "protein mapping is based on normalised peptide-sequence matching "
        "to the workbook."
    )

    reference_path = OUT / "ZWK_IG_peptide_protein_reference.tsv"
    reference.to_csv(reference_path, sep="\t", index=False)

    # Long-form join: one row per candidate link and mapped workbook protein.
    # Unmatched candidate links are retained with a blank protein field.
    long_reference = candidates.merge(
        mapping,
        on="peptide_sequence_key",
        how="left",
        validate="many_to_many"
    )
    long_reference["peptide_protein_status"] = "workbook_sequence_match"
    long_reference.loc[
        long_reference["workbook_protein"].isna(),
        "peptide_protein_status"
    ] = "not_found_in_workbook"

    long_path = OUT / "ZWK_IG_peptide_protein_reference_long.tsv"
    long_reference.to_csv(long_path, sep="\t", index=False)

    # Summarise the build and its limitations.
    summary = {
        "candidate_input": str(CANDIDATES),
        "workbook_input": str(WORKBOOK),
        "psm_ig_candidate_rows": int(len(candidates)),
        "unique_candidate_ig_ids": int(candidates[
            ["sample_id", "reconstructed_ig_id"]
        ].drop_duplicates().shape[0]),
        "unique_candidate_peptide_strings": int(
            candidates["peptide"].nunique()
        ),
        "workbook_normalized_peptide_rows_used": int(len(workbook)),
        "unique_normalized_sequence_keys_in_workbook": int(
            workbook["peptide_sequence_key"].nunique()
        ),
        "unique_sequence_protein_mappings": int(len(mapping)),
        "candidate_rows_with_single_protein_mapping": int(
            (reference["n_workbook_proteins"] == 1).sum()
        ),
        "candidate_rows_with_multiple_protein_mappings": int(
            (reference["n_workbook_proteins"] > 1).sum()
        ),
        "candidate_rows_without_workbook_mapping": int(
            (reference["n_workbook_proteins"] == 0).sum()
        ),
        "candidate_links_confirmed": False,
        "workbook_ig_ids_joined_directly_to_reconstructed_ig_ids": False,
        "normalization_note": (
            "Modification annotations are removed and alphabetic characters "
            "retained for sequence matching. Review unusual notation and "
            "potential sequence-key collisions before treating matches as exact."
        ),
        "pool_mapping_note": (
            "No mapping between pooled acquisitions and individual ZWK samples "
            "is inferred."
        ),
        "outputs": {
            "reference": str(reference_path),
            "reference_long": str(long_path),
            "peptide_protein_mapping": str(mapping_path),
        },
    }

    summary_path = OUT / "ZWK_reference_build_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print("\n=== BUILD SUMMARY ===")
    print(json.dumps(summary, indent=2))
    print("\nOutputs written under:", OUT)

if __name__ == "__main__":
    main()
