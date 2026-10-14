#!/usr/bin/env python3
"""Audit MSstats export granularity without loading whole files into memory."""

import csv
from pathlib import Path

FILES = [
    Path("search_98/msstats.csv"),
    Path("search_98/msstats_ptm.csv"),
]
OUT = Path("reports/msstats_granularity_audit_v1.tsv")

def audit(path):
    runs = set()
    charges = set()
    peptides = set()
    precursor_keys = set()
    conditions = set()
    intensity_count = 0
    nonzero_intensity_count = 0
    total_rows = 0
    header = []

    with path.open("r", newline="", errors="replace") as handle:
        reader = csv.DictReader(handle)
        header = reader.fieldnames or []

        for row in reader:
            total_rows += 1
            runs.add(row.get("Run", ""))
            charges.add(row.get("PrecursorCharge", ""))
            peptides.add(row.get("PeptideSequence", ""))
            precursor_keys.add((
                row.get("PeptideSequence", ""),
                row.get("PrecursorCharge", ""),
            ))
            conditions.add(row.get("Condition", ""))

            value = row.get("Intensity", "")
            if value not in ("", None):
                intensity_count += 1
                try:
                    if float(value) > 0:
                        nonzero_intensity_count += 1
                except ValueError:
                    pass

    lower_headers = {column.lower() for column in header}
    return {
        "file": str(path),
        "size_bytes": path.stat().st_size,
        "rows": total_rows,
        "unique_runs": len(runs),
        "unique_peptide_sequences": len(peptides),
        "unique_peptide_charge_pairs": len(precursor_keys),
        "unique_precursor_charges": len(charges),
        "unique_conditions": len(conditions),
        "intensity_values_present": intensity_count,
        "positive_intensity_values": nonzero_intensity_count,
        "has_mz_column": any("mz" in column for column in lower_headers),
        "has_retention_time_column": any(
            token in column for column in lower_headers
            for token in ("retention time", "retention_time", "rt_seconds")
        ),
        "has_feature_id_column": any(
            "feature" in column or "isotope group" in column
            for column in lower_headers
        ),
        "headers": " | ".join(header),
    }


def main():
    results = []
    for path in FILES:
        if not path.exists():
            print(f"Missing: {path}")
            continue
        result = audit(path)
        results.append(result)
        print(f"\n{result['file']}")
        for key, value in result.items():
            if key not in {"file", "headers"}:
                print(f"  {key}: {value}")
        print(f"  headers: {result['headers']}")

    if results:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        with OUT.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(results[0]), delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(results)
        print(f"\nAudit TSV: {OUT}")


if __name__ == "__main__":
    main()
