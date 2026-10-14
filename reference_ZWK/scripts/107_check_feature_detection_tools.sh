#!/usr/bin/env bash
set -u

OUT="reports/feature_detection_tool_inventory_v1.txt"
mkdir -p reports

{
    echo "=== Timestamp ==="
    date

    echo
    echo "=== Host and OS ==="
    hostname
    uname -a
    cat /etc/os-release 2>/dev/null || true

    echo
    echo "=== Environment ==="
    printf 'CONDA_PREFIX=%s\n' "${CONDA_PREFIX:-}"
    printf 'VIRTUAL_ENV=%s\n' "${VIRTUAL_ENV:-}"
    printf 'PATH=%s\n' "$PATH"
    printf 'LD_LIBRARY_PATH=%s\n' "${LD_LIBRARY_PATH:-}"

    echo
    echo "=== Relevant executables ==="
    for cmd in OpenMS TOPPFeatureFinderMetabo FeatureFinderMetabo \
               FeatureFinderCentroided FeatureFinderIdentification \
               FileConverter PeakPickerHiRes MSFragger IonQuant R python \
               python3 apptainer singularity; do
        printf '%-32s' "$cmd"
        command -v "$cmd" 2>/dev/null || echo "not found on PATH"
    done

    echo
    echo "=== OpenMS-related environment modules ==="
    if command -v module >/dev/null 2>&1; then
        module avail OpenMS 2>&1 || true
        module avail openms 2>&1 || true
    else
        echo "module command not available in this shell"
    fi

    echo
    echo "=== Python package availability ==="
    python - <<'PY'
import importlib.util
for package in ("pyopenms", "numpy", "pandas", "pymzml", "psims"):
    print(f"{package}: {'available' if importlib.util.find_spec(package) else 'not installed'}")
PY

    echo
    echo "=== Common installation locations ==="
    for dir in \
        "$HOME/.local/bin" \
        "$HOME/miniconda3/bin" \
        "$HOME/mambaforge/bin" \
        "$CEUADMIN/openms" \
        /usr/local/Cluster-Apps/ceuadmin/openms \
        /usr/local/bin; do
        [ -d "$dir" ] || continue
        echo "-- $dir"
        find "$dir" -maxdepth 3 \
            \( -iname '*FeatureFinderMetabo*' \
            -o -iname '*OpenMS*' \
            -o -iname 'FileConverter' \
            -o -iname 'PeakPickerHiRes' \) \
            -print 2>/dev/null | head -50
    done

    echo
    echo "=== Existing mzML sample ==="
    ls -lh \
      /rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/ZWK/mzML/szwk900404i19201xms1.mzML

} 2>&1 | tee "$OUT"

echo
echo "Inventory written to $OUT"
