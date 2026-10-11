#!/usr/bin/env bash
# 03_validate_mzml.sh
# Validate existing mzML files without modifying them.

set -eo pipefail

BASE="${HOME}/Caprion/analysis/ZWK"
MZML_DIR="${BASE}/mzML"
REPORT_DIR="${HOME}/Caprion/analysis/pipeline/reports"
REPORT="${REPORT_DIR}/mzml_validation.tsv"

mkdir -p "$REPORT_DIR"

if [[ ! -d "$MZML_DIR" ]]; then
    echo "ERROR: mzML directory not found: $MZML_DIR" >&2
    exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required." >&2
    exit 1
fi

python3 - "$MZML_DIR" "$REPORT" <<'PY'
import os
import sys
import gzip
import xml.etree.ElementTree as ET

mzml_dir, report = sys.argv[1:3]

files = sorted(
    f for f in os.listdir(mzml_dir)
    if f.lower().endswith(".mzml")
    and os.path.isfile(os.path.join(mzml_dir, f))
)

if not files:
    print(f"ERROR: No mzML files found in {mzml_dir}", file=sys.stderr)
    sys.exit(1)

def local_name(tag):
    return tag.rsplit("}", 1)[-1]

def inspect(path):
    size = os.path.getsize(path)
    if size == 0:
        return size, "FAIL", "empty file", 0, 0

    opener = gzip.open if path.lower().endswith(".gz") else open
    ms2 = 0
    spectra = 0
    root_seen = False

    try:
        # iterparse avoids loading the entire XML document into memory.
        with opener(path, "rb") as handle:
            for event, elem in ET.iterparse(handle, events=("start", "end")):
                tag = local_name(elem.tag)

                if event == "start" and tag == "mzML":
                    root_seen = True

                if event == "end" and tag == "spectrum":
                    spectra += 1
                    for child in elem:
                        if local_name(child.tag) != "cvParam":
                            continue
                        if child.attrib.get("accession") == "MS:1000511":
                            if child.attrib.get("value") == "2":
                                ms2 += 1
                            break
                    elem.clear()

        if not root_seen:
            return size, "FAIL", "mzML root element not found", spectra, ms2
        if spectra == 0:
            return size, "FAIL", "no spectra found", spectra, ms2
        if ms2 == 0:
            return size, "WARN", "no MS2 spectra found", spectra, ms2

        return size, "PASS", "XML parsed; spectra and MS2 found", spectra, ms2

    except Exception as exc:
        return size, "FAIL", f"XML parsing error: {exc}", spectra, ms2

results = []

for filename in files:
    path = os.path.join(mzml_dir, filename)
    size, status, message, spectra, ms2 = inspect(path)
    results.append((filename, size, status, spectra, ms2, message))

with open(report, "w", encoding="utf-8") as out:
    out.write("filename\tsize_bytes\tstatus\ttotal_spectra\tms2_spectra\tmessage\n")
    for row in results:
        filename, size, status, spectra, ms2, message = row
        out.write(
            f"{filename}\t{size}\t{status}\t{spectra}\t{ms2}\t"
            f"{message.replace(chr(9), ' ').replace(chr(10), ' ')}\n"
        )

passed = sum(row[2] == "PASS" for row in results)
warned = sum(row[2] == "WARN" for row in results)
failed = sum(row[2] == "FAIL" for row in results)

print(f"mzML files checked: {len(results)}")
print(f"PASS: {passed}")
print(f"WARN: {warned}")
print(f"FAIL: {failed}")
print(f"Report: {report}")

if failed:
    print("\nFiles requiring attention:")
    for filename, size, status, spectra, ms2, message in results:
        if status == "FAIL":
            print(f"{filename}: {message}")
    sys.exit(1)
PY
