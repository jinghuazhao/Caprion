#!/usr/bin/env bash
# 05_prepare_workflow.sh
# Create a separate workflow copy; preserve source files unchanged.

set -eo pipefail

PIPELINE_DIR="${HOME}/Caprion/analysis/pipeline"
SOURCE="${HOME}/Caprion/analysis/ZWK/discovery_fresh/LFQ-MBR-Sage-comparison.workflow"
DB="${HOME}/Caprion/analysis/ZWK/discovery_fresh/database_canonical/2026-10-09-decoys-UP000005640-canonical-proteins-20420.fasta.fas"

DEST="${PIPELINE_DIR}/LFQ-MBR-Sage-comparison.workflow"
CONFIG="${PIPELINE_DIR}/workflow_config.txt"

if [[ ! -f "$SOURCE" ]]; then
    echo "ERROR: Source workflow not found: $SOURCE" >&2
    exit 1
fi

if [[ ! -s "$DB" ]]; then
    echo "ERROR: Validated FASTA not found: $DB" >&2
    exit 1
fi

if [[ -e "$DEST" || -e "$CONFIG" ]]; then
    echo "ERROR: Destination or configuration already exists; refusing to overwrite." >&2
    echo "Workflow: $DEST" >&2
    echo "Config:   $CONFIG" >&2
    exit 1
fi

cp -p "$SOURCE" "$DEST"

# Strip CR only for validation; preserve the copied file's original line endings.
if ! tr -d '\r' < "$DEST" | grep -Fqx 'database.decoy-tag=rev_'; then
    echo "ERROR: Expected decoy-tag setting missing." >&2
    rm -f "$DEST"
    exit 1
fi

for setting in \
    'ionquant.maxlfq=1' \
    'ionquant.mbr=1' \
    'ionquant.run-ionquant=true' \
    'ionquant.use-lfq=true' \
    'workflow.input.data-type.regular-ms=true'
do
    if ! tr -d '\r' < "$DEST" | grep -Fqx "$setting"; then
        echo "ERROR: Expected setting missing: $setting" >&2
        rm -f "$DEST"
        exit 1
    fi
done

cat > "$CONFIG" <<EOF
# Pipeline workflow configuration
workflow_source=$SOURCE
workflow_copy=$DEST
validated_database=$DB
fragpipe_install=/usr/local/Cluster-Apps/ceuadmin/FragPipe/22.0
decoy_prefix=rev_
threads=8
ram_gb=24
EOF

echo "Workflow copy created: $DEST"
echo "Configuration written: $CONFIG"
echo
echo "Validated settings:"
tr -d '\r' < "$DEST" |
    grep -nE 'database.decoy-tag=|ionquant.maxlfq=|ionquant.mbr=|ionquant.run-ionquant=|ionquant.use-lfq=|workflow.input.data-type.regular-ms='

echo
echo "Source workflow remains unchanged."
echo "No FragPipe search was run."
