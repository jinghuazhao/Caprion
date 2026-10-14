#!/usr/bin/env python3
"""Locate FragPipe PSM tables and identify the exact acquisition."""

from pathlib import Path
import csv

ROOT = Path("search_98")
ACQUISITION = "szwk900404i19201xms1"

print("Search root:", ROOT.resolve())
print("Acquisition:", ACQUISITION)

candidates = sorted(
    p for p in ROOT.rglob("psm.tsv")
    if p.is_file()
)

print("\nPSM tables found:", len(candidates))

exact = []
related = []

for p in candidates:
    try:
        with p.open(newline="", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, [])
            first_row = next(reader, [])
    except (OSError, StopIteration):
        continue

    # Match acquisition names in the file path, header, or first row.
    text = " ".join([str(p), *header, *first_row])
    if ACQUISITION in text:
        related.append((p, header, first_row))

    if ACQUISITION in str(p):
        exact.append((p, header, first_row))

print("\nTables whose paths contain the acquisition:", len(exact))
for p, header, first_row in exact:
    print("\nPATH:", p)
    print("HEADER:", "\t".join(header))
    print("FIRST ROW:", "\t".join(first_row[:8]))

print("\nTables whose path/header/first row mention the acquisition:",
      len(related))
for p, header, first_row in related[:20]:
    print("\nPATH:", p)
    print("HEADER:", "\t".join(header))
    print("FIRST ROW:", "\t".join(first_row[:8]))

if not candidates:
    print("\nNo psm.tsv files found. Inspect other search output formats next.")
