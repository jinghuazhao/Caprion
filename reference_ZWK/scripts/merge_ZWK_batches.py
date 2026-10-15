#!/usr/bin/env python3
import csv
import re
import sys
from pathlib import Path
from collections import Counter

REPORTS = Path("reports")
MANIFEST = Path("config/ZWK-reference-panel.fp-manifest")
PREFIXES = ["A", "A_resume"] + [f"A_retry{i:02d}" for i in range(4)]

TABLES = {
    "reconstructed_igs": "reconstructed_igs.tsv",
    "ig_quality": "ig_quality.tsv",
    "ig_peaks": "ig_peaks.tsv",
}

OUTPUTS = {
    name: REPORTS / f"ZWK_merged_{name}.tsv"
    for name in TABLES
}


def basename(value):
    return Path(value.strip()).name


def read_tsv(path):
    f = path.open(newline="", encoding="utf-8-sig")
    reader = csv.DictReader(f, delimiter="\t")
    if not reader.fieldnames:
        f.close()
        raise ValueError(f"Missing TSV header: {path}")
    return f, reader


def main():
    # Do not overwrite any existing merged output.
    existing = [str(p) for p in OUTPUTS.values() if p.exists()]
    if existing:
        sys.exit(
            "ERROR: output already exists; nothing overwritten:\n  "
            + "\n  ".join(existing)
            + "\nMove or rename these files before rerunning."
        )

    if not MANIFEST.exists():
        sys.exit(f"Missing manifest: {MANIFEST}")

    expected = set()
    for line in MANIFEST.read_text().splitlines():
        for field in line.split():
            if re.search(r"xms1\.mzML$", field):
                expected.add(basename(field))

    if len(expected) != 94:
        sys.exit(f"Expected 94 distinct MS1 files, found {len(expected)}.")

    # COMPLETE files are tracked separately for each run.
    done_by_prefix = {}
    completed_union = set()

    for prefix in PREFIXES:
        path = REPORTS / f"{prefix}_progress.tsv"
        if not path.exists():
            sys.exit(f"Missing progress file: {path}")

        with path.open(newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if not reader.fieldnames or not {
                "source_file", "status"
            }.issubset(reader.fieldnames):
                sys.exit(f"Unexpected progress-file columns: {path}")

            done = {
                basename(row["source_file"])
                for row in reader
                if row.get("status", "").strip().upper() == "COMPLETE"
            }

        done_by_prefix[prefix] = done
        completed_union.update(done)

    missing = expected - completed_union
    unexpected = completed_union - expected
    if missing or unexpected:
        sys.exit(
            f"Source coverage failed: missing={len(missing)}, "
            f"unexpected={len(unexpected)}. No outputs created."
        )

    print(f"Manifest coverage: {len(completed_union)}/94 MS1 files")

    # Write temporary outputs; publish only after validation.
    temp_outputs = {
        name: path.with_suffix(path.suffix + ".tmp")
        for name, path in OUTPUTS.items()
    }
    for path in temp_outputs.values():
        if path.exists():
            sys.exit(f"Temporary output already exists: {path}")

    total_rows = Counter()
    retained_ids = set()
    duplicate_ids = Counter()
    source_rows = Counter()
    schema_by_table = {}

    try:
        for table_name, suffix in TABLES.items():
            output_tmp = temp_outputs[table_name]
            reference_header = None

            with output_tmp.open("w", newline="", encoding="utf-8") as out:
                writer = None

                for prefix in PREFIXES:
                    path = REPORTS / f"{prefix}_{suffix}"
                    if not path.exists():
                        raise FileNotFoundError(f"Missing table: {path}")

                    f, reader = read_tsv(path)
                    try:
                        header = reader.fieldnames

                        if "source_file" not in header:
                            raise ValueError(
                                f"{path}: no source_file column; "
                                f"columns={header}"
                            )

                        if reference_header is None:
                            reference_header = header
                            writer = csv.DictWriter(
                                out, fieldnames=header,
                                delimiter="\t", lineterminator="\n",
                                extrasaction="raise"
                            )
                            writer.writeheader()
                        elif header != reference_header:
                            raise ValueError(
                                f"Column mismatch in {path}\n"
                                f"Expected: {reference_header}\n"
                                f"Found:    {header}"
                            )

                        for row in reader:
                            source = basename(row["source_file"])
                            if source not in done_by_prefix[prefix]:
                                continue

                            writer.writerow(row)
                            total_rows[table_name] += 1
                            source_rows[(table_name, source)] += 1

                            if table_name == "reconstructed_igs":
                                if "reconstructed_ig_id" not in row:
                                    raise ValueError(
                                        f"{path}: missing reconstructed_ig_id"
                                    )
                                gid = row["reconstructed_ig_id"].strip()
                                if not gid:
                                    raise ValueError(
                                        f"{path}: blank reconstructed_ig_id"
                                    )
                                if gid in retained_ids:
                                    duplicate_ids[gid] += 1
                                else:
                                    retained_ids.add(gid)
                    finally:
                        f.close()

            schema_by_table[table_name] = reference_header

        if duplicate_ids:
            examples = list(duplicate_ids.items())[:20]
            raise ValueError(
                f"Duplicate reconstructed IG IDs: {len(duplicate_ids)}; "
                f"examples: {examples}"
            )

        if total_rows["reconstructed_igs"] != len(retained_ids):
            raise ValueError("Reconstructed IG row count does not match unique IDs.")

        # All validation passed; publish outputs.
        for name, final_path in OUTPUTS.items():
            temp_outputs[name].replace(final_path)

    except Exception:
        for path in temp_outputs.values():
            path.unlink(missing_ok=True)
        raise

    print("\nMERGE PASSED")
    print(f"Retained reconstructed IG rows: {total_rows['reconstructed_igs']:,}")
    print(f"Unique reconstructed IG IDs:    {len(retained_ids):,}")
    print(f"Duplicate reconstructed IDs:    {len(duplicate_ids):,}")
    print(f"Retained quality rows:          {total_rows['ig_quality']:,}")
    print(f"Retained peak rows:             {total_rows['ig_peaks']:,}")

    print("\nCreated files:")
    for name, path in OUTPUTS.items():
        print(f"  {path} ({path.stat().st_size:,} bytes)")

    print("\nNote: quality and peak row counts need not equal IG row counts.")


if __name__ == "__main__":
    main()
