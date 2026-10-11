#!/usr/bin/env bash
# 06_build_manifest.sh
# Build and validate a production FragPipe manifest from RAW metadata and mzML files.
# Use the four-column FragPipe manifest format; preserve existing outputs and pilot data.

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
RAW_DIR="${HOME}/Caprion/analysis/ZWK/raw"
MZML_DIR="${HOME}/Caprion/analysis/ZWK/mzML"
REPORT_DIR="${PIPELINE}/reports"

# Versioned outputs: never overwrite the original three-column outputs.
MANIFEST="${PIPELINE}/production_v2.fp-manifest"
SAMPLE_METADATA="${PIPELINE}/production_sample_metadata_v2.tsv"
REPORT="${REPORT_DIR}/manifest_validation_v2.txt"

mkdir -p "$REPORT_DIR"

for f in "$MANIFEST" "$SAMPLE_METADATA" "$REPORT"; do
    if [[ -e "$f" ]]; then
        echo "ERROR: Refusing to overwrite existing file: $f" >&2
        exit 1
    fi
done

if [[ ! -d "$RAW_DIR" || ! -d "$MZML_DIR" ]]; then
    echo "ERROR: RAW or mzML directory not found." >&2
    exit 1
fi

python3 - "$RAW_DIR" "$MZML_DIR" "$MANIFEST" "$SAMPLE_METADATA" "$REPORT" <<'PY'
from pathlib import Path
from collections import Counter
import re
import sys

raw_dir, mzml_dir, manifest_path, metadata_path, report_path = map(Path, sys.argv[1:])

def read_fields(path):
    fields = {}
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line and not line.startswith("#"):
            key, value = line.split("=", 1)
            fields[key.strip()] = value.strip()
    return fields

metadata_files = sorted(raw_dir.glob("*-metadata.txt"))
mzml_files = sorted(mzml_dir.glob("*.mzML"))

if not metadata_files:
    sys.exit("ERROR: No *-metadata.txt files found.")
if not mzml_files:
    sys.exit("ERROR: No mzML files found.")

mzml_by_stem = {}
for p in mzml_files:
    key = p.stem.casefold()
    if key in mzml_by_stem:
        sys.exit(f"ERROR: Duplicate mzML stem: {p.stem}")
    mzml_by_stem[key] = p

records = []
problems = []
seen_raw = set()

for p in metadata_files:
    fields = read_fields(p)
    raw_name = fields.get("RAW file path", "")
    sample_id = fields.get("Sample id", "")
    created = fields.get("Creation date", "")
    row_number = fields.get("Sample row number", "")

    if not raw_name or not sample_id:
        problems.append(f"Missing RAW file path or Sample id in {p.name}")
        continue

    stem = Path(raw_name).stem.casefold()
    if stem in seen_raw:
        problems.append(f"Duplicate RAW metadata stem: {stem}")
        continue
    seen_raw.add(stem)

    mzml = mzml_by_stem.get(stem)
    if mzml is None:
        problems.append(f"No matching mzML for {raw_name} ({sample_id})")
        continue

    m = re.fullmatch(r"(IgY14|MARS14)_Pool(\d+)_Fraction\s+(\d+)", sample_id)
    if m:
        method, pool, fraction = m.groups()
        experiment = f"{method}_Pool{pool}"
        category = "fraction"
    elif sample_id.startswith("PQC-3"):
        experiment = "PQC-3_QC"
        method, pool, fraction = "PQC-3", "", ""
        category = "qc"
    else:
        problems.append(f"Unrecognised Sample id format: {sample_id} ({p.name})")
        continue

    records.append({
        "mzml": str(mzml.resolve()),
        "experiment": experiment,
        "sample_id": sample_id,
        "raw_name": raw_name,
        "creation_date": created,
        "row_number": row_number,
        "method": method,
        "pool": pool,
        "fraction": fraction,
        "category": category,
    })

metadata_stems = {
    Path(read_fields(p).get("RAW file path", "")).stem.casefold()
    for p in metadata_files
    if read_fields(p).get("RAW file path", "")
}
mzml_stems = set(mzml_by_stem)
mzml_without_metadata = sorted(mzml_stems - metadata_stems)

if mzml_without_metadata:
    problems.extend(
        f"mzML has no matching metadata: {mzml_by_stem[s].name}"
        for s in mzml_without_metadata
    )

paths = [r["mzml"] for r in records]
if len(paths) != len(set(paths)):
    problems.append("Duplicate mzML paths in generated mapping.")

if problems:
    print("ERROR: Manifest validation failed; no output files written.", file=sys.stderr)
    for problem in problems:
        print(f" - {problem}", file=sys.stderr)
    sys.exit(1)

if len(records) != len(metadata_files) or len(records) != len(mzml_files):
    sys.exit(
        f"ERROR: Expected one-to-one mapping, got {len(metadata_files)} metadata, "
        f"{len(mzml_files)} mzML and {len(records)} matched records."
    )

# FragPipe format: path, experiment, empty third field, acquisition type.
manifest_text = "".join(
    f'{r["mzml"]}\t{r["experiment"]}\t\tDDA\n'
    for r in sorted(records, key=lambda x: x["mzml"])
)

metadata_header = [
    "mzML_path", "RAW_file", "sample_id", "experiment",
    "method_or_control", "pool", "fraction", "creation_date",
    "sample_row_number", "category"
]
metadata_lines = ["\t".join(metadata_header)]

for r in sorted(records, key=lambda x: x["mzml"]):
    vals = [
        r["mzml"], r["raw_name"], r["sample_id"], r["experiment"],
        r["method"], r["pool"], r["fraction"], r["creation_date"],
        r["row_number"], r["category"]
    ]
    metadata_lines.append("\t".join(v.replace("\t", " ") for v in vals))

counts = Counter(r["experiment"] for r in records)
sample_counts = Counter(r["sample_id"] for r in records)

report_lines = [
    "Production manifest validation",
    "==============================",
    f"Metadata files: {len(metadata_files)}",
    f"mzML files: {len(mzml_files)}",
    f"Matched records: {len(records)}",
    f"Manifest rows: {len(records)}",
    f"Unique mzML paths: {len(set(paths))}",
    "Manifest format: four tab-separated fields",
    "Third field: empty",
    "Acquisition type: DDA",
    "",
    "Runs by experiment:",
]
report_lines.extend(f"  {name}: {n}" for name, n in sorted(counts.items()))
report_lines.extend(["", "Sample IDs occurring more than once:"])

duplicates = [(name, n) for name, n in sample_counts.items() if n > 1]
report_lines.extend(
    [f"  {name}: {n}" for name, n in sorted(duplicates)]
    if duplicates else ["  None"]
)

report_lines.extend([
    "",
    "Known incomplete fraction coverage:",
    "  IgY14 Pool5 is missing fractions 1, 10 and 11 in the metadata inventory.",
    "",
    "Grouping note:",
    "  Experiments are assigned by depletion method and pool; PQC-3 is separate.",
    "  These labels preserve provenance and are not biological-replicate claims.",
    "",
    "Status: PASS",
])

created = []
try:
    for path, content in [
        (manifest_path, manifest_text),
        (metadata_path, "\n".join(metadata_lines) + "\n"),
        (report_path, "\n".join(report_lines) + "\n"),
    ]:
        with path.open("x") as out:
            out.write(content)
        created.append(path)
except Exception:
    for path in created:
        path.unlink(missing_ok=True)
    raise

print("\n".join(report_lines))
print("\nCreated:")
print(f"  {manifest_path}")
print(f"  {metadata_path}")
print(f"  {report_path}")
PY
