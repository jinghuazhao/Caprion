#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

R = Path("reports")
CLEANUP = "--cleanup" in sys.argv

FILES = {
    "igs": R / "ZWK_merged_reconstructed_igs.tsv",
    "quality": R / "ZWK_merged_ig_quality.tsv",
    "peaks": R / "ZWK_merged_ig_peaks.tsv",
}

# Intermediate data tables only. Keep logs, progress, summaries and manifests.
PREFIXES = (
    ["A", "A_resume"]
    + [f"A_retry{i:02d}" for i in range(4)]
    + [f"A_batch{i:02d}" for i in range(4)]
)
SUFFIXES = (
    "_reconstructed_igs.tsv",
    "_ig_quality.tsv",
    "_ig_peaks.tsv",
)

EXPECTED_IGS = 5_574_224
EXPECTED_QUALITY = 5_574_224
EXPECTED_PEAKS = 218_987_171


def fail(message):
    raise RuntimeError(message)


def scan_table(path, required_columns, callback):
    """Stream a TSV, checking its header and every row's column count."""
    count = 0
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        if header is None:
            fail(f"Empty file: {path}")

        if len(header) != len(set(header)):
            fail(f"Duplicate column names in {path}: {header}")

        missing = set(required_columns) - set(header)
        if missing:
            fail(f"{path}: missing columns {sorted(missing)}")

        positions = {name: header.index(name) for name in required_columns}

        for line_no, row in enumerate(reader, start=2):
            if len(row) != len(header):
                fail(
                    f"{path}: line {line_no} has {len(row)} fields; "
                    f"expected {len(header)}"
                )
            callback(row, positions, line_no)
            count += 1

    print(f"  {path.name}: {count:,} rows OK")
    return count


def validate():
    for path in FILES.values():
        if not path.is_file():
            fail(f"Missing merged file: {path}")

    # 1. Reconstructed IGs: establish the authoritative ID set.
    ig_ids = set()

    def check_ig(row, pos, line):
        gid = row[pos["reconstructed_ig_id"]].strip()
        if not gid:
            fail(f"Blank IG ID at reconstructed-IG line {line}")
        if gid in ig_ids:
            fail(f"Duplicate reconstructed IG ID: {gid}")
        ig_ids.add(gid)

    print("\n1. Checking reconstructed IG table")
    n_igs = scan_table(
        FILES["igs"],
        ["reconstructed_ig_id", "sample_id", "source_file"],
        check_ig,
    )
    if n_igs != EXPECTED_IGS:
        fail(f"Expected {EXPECTED_IGS:,} IG rows; found {n_igs:,}")

    # 2. Quality: exactly one record per reconstructed IG.
    quality_ids = set()

    def check_quality(row, pos, line):
        gid = row[pos["reconstructed_ig_id"]].strip()
        if not gid:
            fail(f"Blank quality IG ID at line {line}")
        if gid not in ig_ids:
            fail(f"Quality record references unknown IG: {gid}")
        if gid in quality_ids:
            fail(f"Duplicate quality record for IG: {gid}")
        quality_ids.add(gid)

    print("\n2. Checking quality table")
    n_quality = scan_table(
        FILES["quality"],
        ["reconstructed_ig_id", "sample_id", "source_file"],
        check_quality,
    )
    if n_quality != EXPECTED_QUALITY:
        fail(f"Expected {EXPECTED_QUALITY:,} quality rows; found {n_quality:,}")
    if quality_ids != ig_ids:
        fail("Quality IDs do not exactly match reconstructed IG IDs")

    del quality_ids

    # 3. Peaks: every peak must refer to a known reconstructed IG.
    peak_count = 0

    def check_peak(row, pos, line):
        nonlocal peak_count
        gid = row[pos["reconstructed_ig_id"]].strip()
        if not gid:
            fail(f"Blank peak IG ID at line {line}")
        if gid not in ig_ids:
            fail(f"Peak references unknown IG: {gid} (line {line})")
        peak_count += 1

    print("\n3. Checking peak table (this may take a while)")
    n_peaks = scan_table(
        FILES["peaks"],
        ["reconstructed_ig_id", "sample_id", "source_file",
         "spectrum_index_zero_based", "isotope_index", "centroid_mz"],
        check_peak,
    )
    if n_peaks != EXPECTED_PEAKS:
        fail(f"Expected {EXPECTED_PEAKS:,} peak rows; found {n_peaks:,}")

    print("\nVALIDATION PASSED")
    print(f"Reconstructed IGs: {n_igs:,}")
    print(f"Quality records:   {n_quality:,} (one per IG)")
    print(f"Peak records:      {n_peaks:,}")
    print("All quality and peak IG references are valid.")

    if not CLEANUP:
        print("\nNo files deleted.")
        print("After reviewing the results, run with --cleanup to remove")
        print("intermediate data tables while retaining audit files.")
        return

    # Delete only known intermediate data tables, never merged outputs.
    candidates = [
        R / f"{prefix}{suffix}"
        for prefix in PREFIXES
        for suffix in SUFFIXES
    ]
    existing = [p for p in candidates if p.is_file()]

    print(f"\nCleanup will remove {len(existing)} intermediate TSV files.")
    for p in existing:
        print(f"  {p} ({p.stat().st_size:,} bytes)")

    # Only proceed after every validation above has passed.
    for p in existing:
        p.unlink()

    print(f"\nCleanup complete: removed {len(existing)} intermediate TSV files.")
    print("Merged files and audit logs were preserved.")


if __name__ == "__main__":
    try:
        validate()
    except Exception as exc:
        print(f"\nVALIDATION FAILED: {exc}", file=sys.stderr)
        print("No cleanup was performed.", file=sys.stderr)
        sys.exit(1)
