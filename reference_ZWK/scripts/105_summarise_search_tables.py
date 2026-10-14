#!/usr/bin/env python3
"""Summarise distinct tabular outputs in the completed FragPipe search."""

from collections import defaultdict
from pathlib import Path
import csv

ROOT = Path("search_98")
OUT = Path("reports/search_98_table_types_v1.tsv")

groups = defaultdict(list)

for path in ROOT.rglob("*"):
    if not path.is_file():
        continue
    if path.suffix.lower() not in {".tsv", ".txt", ".csv"}:
        continue
    if "MSBooster_plots" in path.parts:
        continue

    try:
        size = path.stat().st_size
        with path.open("r", newline="", errors="replace") as handle:
            delimiter = "\t" if path.suffix.lower() == ".tsv" else (
                "," if path.suffix.lower() == ".csv" else None
            )
            if delimiter:
                header = next(csv.reader(handle, delimiter=delimiter), [])
            else:
                header = [handle.readline().strip()]
        header_text = " | ".join(header[:30])
    except Exception as exc:
        size = path.stat().st_size
        header_text = f"[read error: {type(exc).__name__}: {exc}]"

    # Group by relative parent and basename; repeated per-run outputs
    # will therefore be easy to distinguish from global outputs.
    key = (path.name, header_text)
    groups[key].append((path, size))

rows = []
for (name, header), files in sorted(groups.items()):
    rows.append({
        "filename": name,
        "file_count": len(files),
        "total_size_bytes": sum(size for _, size in files),
        "example_path": str(files[0][0]),
        "header_preview": header,
    })

with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(rows[0]),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Distinct filename/header groups: {len(rows)}")
print(f"Summary: {OUT}")

keywords = (
    "feature", "isotope", "precursor", "chrom", "ms1",
    "quant", "ion.tsv", "trace", "match",
)
for row in rows:
    if any(k in row["filename"].lower() or k in row["header_preview"].lower()
           for k in keywords):
        print(
            f'\n{row["filename"]} | files={row["file_count"]} '
            f'| total={row["total_size_bytes"]:,} bytes'
        )
        print(f'  example: {row["example_path"]}')
        print(f'  header: {row["header_preview"]}')
