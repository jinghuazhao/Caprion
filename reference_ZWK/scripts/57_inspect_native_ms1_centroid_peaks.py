#!/usr/bin/env python3
"""Inspect native mzML MS1 centroid peaks near selected charge-2 trace m/z values.

Diagnostic only. Does not assign monoisotopic anchors or alter isotope groups.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

FEATURE_ID = "10275997325121874585"
CHARGE = 2

# Retention times in seconds, taken from the trace table.
TARGETS = {
    5954: 769.41804,
    5981: 772.57938,
    5988: 773.44866,
}

TARGET_MZ = [
    1341.064209,
    1341.067871,
    1341.560547,
    1341.571289,
    1342.063232,
    1342.067505,
    1342.568115,
    1342.569946,
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mzml",
        required=True,
        help="Original mzML acquisition used to generate isotope_trace_peaks.tsv",
    )
    parser.add_argument(
        "--mz-tolerance",
        type=float,
        default=0.03,
        help="Half-width of the centroid peak search window in Da (default: 0.03)",
    )
    parser.add_argument(
        "--rt-tolerance",
        type=float,
        default=2.0,
        help="Retention-time tolerance in seconds (default: 2.0)",
    )
    args = parser.parse_args()

    mzml = Path(args.mzml).expanduser()
    if not mzml.is_file():
        sys.exit(f"ERROR: mzML file not found: {mzml}")

    try:
        import pyopenms as oms
    except ImportError:
        sys.exit(
            "ERROR: pyopenms is unavailable. Load the OpenMS environment "
            "used for the pilot and retry."
        )

    print(f"Input mzML: {mzml}")
    print(f"Feature ID: {FEATURE_ID}; charge: {CHARGE}")
    print(f"m/z window: +/- {args.mz_tolerance:.4f} Da")
    print(f"RT window: +/- {args.rt_tolerance:.2f} seconds")

    experiment = oms.MSExperiment()
    oms.MzMLFile().load(str(mzml), experiment)
    spectra = list(experiment)

    print(f"Spectra loaded: {len(spectra)}")

    results = []

    for trace_scan_index, target_rt in TARGETS.items():
        nearby = []

        for spectrum_index, spectrum in enumerate(spectra):
            if spectrum.getMSLevel() != 1:
                continue

            rt = float(spectrum.getRT())

            if abs(rt - target_rt) <= args.rt_tolerance:
                nearby.append((spectrum_index, spectrum, rt))

        print(
            f"\nTrace-table scan {trace_scan_index}; "
            f"target RT={target_rt:.5f} s; "
            f"nearby MS1 spectra={len(nearby)}"
        )

        if not nearby:
            print("  No MS1 spectrum found in this RT window.")
            continue

        for spectrum_index, spectrum, rt in nearby:
            mz_array, intensity_array = spectrum.get_peaks()

            print(
                f"\n  mzML spectrum index={spectrum_index}; "
                f"RT={rt:.5f} s; centroid peaks={len(mz_array)}"
            )

            for target_mz in TARGET_MZ:
                matches = []

                for mz, intensity in zip(mz_array, intensity_array):
                    mz = float(mz)
                    intensity = float(intensity)

                    if abs(mz - target_mz) <= args.mz_tolerance:
                        matches.append((mz, intensity))

                matches.sort(key=lambda item: abs(item[0] - target_mz))

                if not matches:
                    print(f"    target {target_mz:.6f}: no peak in window")
                    results.append({
                        "trace_table_scan_index": trace_scan_index,
                        "target_rt_sec": target_rt,
                        "mzml_spectrum_index": spectrum_index,
                        "mzml_rt_sec": rt,
                        "target_mz": target_mz,
                        "observed_mz": None,
                        "intensity": None,
                        "mass_error_mDa": None,
                    })
                    continue

                closest_mz, closest_intensity = matches[0]
                error_mDa = (closest_mz - target_mz) * 1000

                print(
                    f"    target {target_mz:.6f}: "
                    f"{len(matches)} peak(s); "
                    f"closest m/z={closest_mz:.6f}; "
                    f"intensity={closest_intensity:.1f}; "
                    f"error={error_mDa:+.3f} mDa"
                )

                for observed_mz, intensity in matches:
                    results.append({
                        "trace_table_scan_index": trace_scan_index,
                        "target_rt_sec": target_rt,
                        "mzml_spectrum_index": spectrum_index,
                        "mzml_rt_sec": rt,
                        "target_mz": target_mz,
                        "observed_mz": observed_mz,
                        "intensity": intensity,
                        "mass_error_mDa": (
                            observed_mz - target_mz
                        ) * 1000,
                    })

    output = Path("reports/native_ms1_centroid_peak_inspection_v2.tsv")
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(output, sep="\t", index=False)

    print(f"\nWrote: {output.resolve()}")
    print("Diagnostic only: no IG assignments or anchors were changed.")


if __name__ == "__main__":
    main()
