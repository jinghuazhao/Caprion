#!/usr/bin/env python3
"""Inventory existing search outputs for MS1 feature-level evidence."""

import csv
from pathlib import Path

ROOT = Path("search_98")
OUT = Path("reports/ionquant_output_inventory_v1.tsv")

KEYWORDS = (
    "ionquant", "feature", "isotope", "ms1", "precursor",
    "chrom", "trace", "quant", "ion.tsv",
)

rows = []

for path in sorted(ROOT.rglob("*")):
    if not path.is_file():
        continue

    name = path.name.lower()
    if not any(keyword in name for keyword in KEYWORDS):
        continue

    row = {
        "relative_path": str(path),
        "size_bytes": path.stat().st_size,
        "extension": path.suffix.lower(),
        "header_preview": "",
        "row_count_estimate": "",
        "inspection_status": "",
    }

    if path.suffix.lower() in {".tsv", ".txt", ".csv"}:
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        try:
            with path.open("r", newline="", errors="replace") as handle:
                reader = csv.reader(handle, delimiter=delimiter)
                header = next(reader, [])
                row["header_preview"] = " | ".join(header[:40])

                # Count data rows without retaining the table in memory.
                row_count = sum(1 for _ in reader)
                row["row_count_estimate"] = row_count
                row["inspection_status"] = "header and rows counted"
        except Exception as exc:
            row["inspection_status"] = f"read error: {type(exc).__name__}: {exc}"
    else:
        row["inspection_status"] = "file inventoried; contents not parsed"

    rows.append(row)

if not rows:
    raise SystemExit(f"No matching files found beneath {ROOT}")

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(rows[0]),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Files inventoried: {len(rows)}")
print(f"Inventory: {OUT}")

for row in rows:
    print(f'\n{row["relative_path"]}')
    print(f'  size: {row["size_bytes"]:,} bytes')
    if row["header_preview"]:
        print(f'  header: {row["header_preview"]}')
        print(f'  data rows: {row["row_count_estimate"]}')
    else:
        print(f'  status: {row["inspection_status"]}')
