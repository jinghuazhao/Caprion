#!/usr/bin/env bash
# 02_convert_missing.sh
# Convert only RAW files without a matching mzML file.
# Existing mzML files are never overwritten.

set -eo pipefail

BASE="${HOME}/Caprion/analysis/ZWK"
RAW_DIR="${BASE}/raw"
MZML_DIR="${BASE}/mzML"
REPORT_DIR="${HOME}/Caprion/analysis/pipeline/reports"
REPORT="${REPORT_DIR}/conversion_log.tsv"

MODULE_NAME="ceuadmin/ThermoRawFileParser/1.4.4"

mkdir -p "$REPORT_DIR" "$MZML_DIR"

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: RAW directory not found: $RAW_DIR" >&2
    exit 1
fi

# Load the module only if conversion is needed.
find_missing_raw() {
    local raw_file raw_stem
    while IFS= read -r -d '' raw_file; do
        raw_stem="${raw_file##*/}"
        raw_stem="${raw_stem%.*}"

        if [[ ! -f "$MZML_DIR/${raw_stem}.mzML" ]] &&
           [[ ! -f "$MZML_DIR/${raw_stem}.mzml" ]]; then
            printf '%s\0' "$raw_file"
        fi
    done < <(
        find "$RAW_DIR" -maxdepth 1 -type f \
            -iname '*.raw' -print0 | sort -z
    )
}

mapfile -d '' -t MISSING_RAW < <(find_missing_raw)

if (( ${#MISSING_RAW[@]} == 0 )); then
    echo "No missing mzML conversions."
    echo "RAW files: all have matching mzML files."
    exit 0
fi

echo "RAW files requiring conversion: ${#MISSING_RAW[@]}"
printf 'Loading module: %s\n' "$MODULE_NAME"

if ! type module >/dev/null 2>&1; then
    echo "ERROR: Environment Modules is unavailable in this shell." >&2
    echo "Start a login shell or initialise the module environment." >&2
    exit 1
fi

module load "$MODULE_NAME"

PARSER="${THERMORAWFILEPARSER:-/usr/local/Cluster-Apps/ceuadmin/ThermoRawFileParser/1.4.4/ThermoRawFileParser.exe}"

if [[ ! -f "$PARSER" ]]; then
    echo "ERROR: ThermoRawFileParser not found: $PARSER" >&2
    exit 1
fi

if ! command -v mono >/dev/null 2>&1; then
    echo "ERROR: mono is unavailable after loading the parser module." >&2
    exit 1
fi

if [[ ! -f "$REPORT" ]]; then
    printf 'raw_file\tmzml_file\tstatus\tdetail\n' > "$REPORT"
fi

for RAW_FILE in "${MISSING_RAW[@]}"; do
    STEM="${RAW_FILE##*/}"
    STEM="${STEM%.*}"
    OUTPUT="${MZML_DIR}/${STEM}.mzML"

    # Recheck before conversion to avoid overwriting an existing file.
    if [[ -e "$OUTPUT" ]]; then
        printf '%s\t%s\tSKIPPED\tOutput appeared before conversion\n' \
            "$RAW_FILE" "$OUTPUT" >> "$REPORT"
        continue
    fi

    echo "Converting: ${RAW_FILE##*/}"

    # Convert to a temporary directory, then move the result into place.
    TMP_DIR="$(mktemp -d "${MZML_DIR}/.convert_${STEM}.XXXXXX")"

    if mono "$PARSER" -i="$RAW_FILE" -o="$TMP_DIR"; then
        GENERATED="${TMP_DIR}/${STEM}.mzML"

        if [[ ! -s "$GENERATED" ]]; then
            GENERATED="$(find "$TMP_DIR" -maxdepth 1 -type f \
                -iname '*.mzML' -print -quit)"
        fi

        if [[ -n "${GENERATED:-}" && -s "$GENERATED" ]]; then
            if [[ -e "$OUTPUT" ]]; then
                printf '%s\t%s\tSKIPPED\tOutput already exists\n' \
                    "$RAW_FILE" "$OUTPUT" >> "$REPORT"
            else
                mv "$GENERATED" "$OUTPUT"
                printf '%s\t%s\tPASS\tConverted successfully\n' \
                    "$RAW_FILE" "$OUTPUT" >> "$REPORT"
            fi
        else
            printf '%s\t%s\tFAIL\tNo non-empty mzML produced\n' \
                "$RAW_FILE" "$OUTPUT" >> "$REPORT"
            rm -rf "$TMP_DIR"
            echo "ERROR: No non-empty mzML produced for $RAW_FILE" >&2
            exit 1
        fi
    else
        printf '%s\t%s\tFAIL\tThermoRawFileParser failed\n' \
            "$RAW_FILE" "$OUTPUT" >> "$REPORT"
        rm -rf "$TMP_DIR"
        echo "ERROR: Conversion failed for $RAW_FILE" >&2
        exit 1
    fi

    rm -rf "$TMP_DIR"
done

echo "Conversion log: $REPORT"
echo "Run 01_inventory.sh again to verify the final inventory."
