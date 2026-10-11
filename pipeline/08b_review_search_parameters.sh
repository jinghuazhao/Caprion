#!/usr/bin/env bash
# =============================================================================
# 08b_review_search_parameters.sh
#
# Purpose:
#   Extract FragPipe/MSFragger search and quantification settings for review
#   before the full 98-file production run.
#
# Inputs:
#   LFQ-MBR-Sage-comparison.workflow
#   dry_run_98_v2/fragpipe.workflow
#
# Output:
#   Terminal report only; no files are modified.
#
# Safety:
#   Read-only. Does not launch FragPipe or overwrite existing results.
#
# Execution:
#   bash 08b_review_search_parameters.sh
# =============================================================================

set -eo pipefail

PIPELINE="${HOME}/Caprion/analysis/pipeline"
WORKFLOW="${PIPELINE}/LFQ-MBR-Sage-comparison.workflow"
GENERATED="${PIPELINE}/dry_run_98_v2/fragpipe.workflow"

for file in "$WORKFLOW" "$GENERATED"; do
    if [[ ! -f "$file" ]]; then
        printf 'ERROR: required file missing: %s\n' "$file" >&2
        exit 1
    fi
done

python - "$WORKFLOW" "$GENERATED" <<'PY'
import sys

def read_properties(path):
    properties = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            properties.setdefault(key.strip(), []).append(value.strip())
    return properties

prepared_path, generated_path = sys.argv[1:]
prepared = read_properties(prepared_path)
generated = read_properties(generated_path)

terms = (
    "database.",
    "workflow.input.",
    "msfragger.",
    "ionquant.",
    "msbooster.",
    "percolator.",
    "proteinprophet.",
    "philosopher.",
    "enzyme",
    "digest",
    "tolerance",
    "mass",
    "modification",
    "isotope",
    "calibration",
    "precursor",
    "fragment",
    "glyco",
    "labile",
    "open-search",
    "open_search",
)

def relevant(key):
    lower = key.lower()
    return any(term in lower for term in terms)

keys = sorted(
    key for key in set(prepared) | set(generated)
    if relevant(key)
)

print("FRAGPIPE SEARCH-PARAMETER AUDIT")
print("===============================")

for key in keys:
    left = prepared.get(key, ["<not present>"])
    right = generated.get(key, ["<not present>"])

    if left == right:
        print(f"\n{key} = {'; '.join(left)}")
    else:
        print(f"\n{key}")
        print(f"  Prepared:  {'; '.join(left)}")
        print(f"  Generated: {'; '.join(right)}")

print("\nAudit only: no settings were changed.")
print("Review the output before starting production.")
PY
