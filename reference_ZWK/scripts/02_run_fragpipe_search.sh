#!/usr/bin/env bash
set -Eeuo pipefail

module load ceuadmin/FragPipe/22.0

PROJECT="$HOME/Caprion/analysis/ZWK/reference_ZWK"
WORKFLOW="$PROJECT/config/ZWK-reference-panel-candidate.workflow"
MANIFEST="$PROJECT/config/ZWK-reference-panel.fp-manifest"
WORKDIR="$PROJECT/search_98"
LOG="$PROJECT/reports/fragpipe_search_98.log"

mkdir -p "$PROJECT/reports"

# Do not overwrite a previous search.
if [[ -e "$WORKDIR" ]]; then
    echo "ERROR: Work directory already exists: $WORKDIR" >&2
    exit 1
fi

# Validate required files.
[[ -f "$WORKFLOW" ]] || { echo "ERROR: Missing workflow: $WORKFLOW" >&2; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: Missing manifest: $MANIFEST" >&2; exit 1; }

awk -F '\t' '
    NF != 4 || $1 == "" || $2 == "" || $4 != "DDA" { bad=1 }
    END {
        if (NR != 98 || bad) {
            print "ERROR: Manifest must contain 98 valid four-column DDA rows."
            exit 1
        }
    }
' "$MANIFEST"

while IFS=$'\t' read -r mzml experiment unused acquisition; do
    [[ -f "$mzml" ]] || { echo "ERROR: Missing mzML: $mzml" >&2; exit 1; }
done < "$MANIFEST"

DB=$(sed -n 's/^database\.db-path=//p' "$WORKFLOW" | head -n 1)
[[ -n "$DB" && -f "$DB" ]] || { echo "ERROR: Missing database: $DB" >&2; exit 1; }

FRAGPIPE=$(command -v fragpipe || true)
if [[ -z "$FRAGPIPE" ]]; then
    echo "ERROR: fragpipe executable not found after module load." >&2
    echo "Check the module's PATH and launcher location." >&2
    exit 1
fi

mkdir "$WORKDIR"

echo "FragPipe executable: $FRAGPIPE"
echo "Workflow: $WORKFLOW"
echo "Manifest: $MANIFEST"
echo "Database: $DB"
echo "Workdir: $WORKDIR"
echo "Started: $(date --iso-8601=seconds)"

set +e
"$FRAGPIPE" \
    --headless \
    --workflow "$WORKFLOW" \
    --manifest "$MANIFEST" \
    --workdir "$WORKDIR" \
    --threads 8 \
    --ram 24 \
    2>&1 | tee "$LOG"
status=${PIPESTATUS[0]}
set -e

echo "Finished: $(date --iso-8601=seconds)" echo "FragPipe exit status: $status" exit "$status"
