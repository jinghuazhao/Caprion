#!/usr/bin/env python3
"""Diagnose isotope peak and adjacent-spacing errors for one candidate."""

from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
TRACE_FILE = BASE / "isotope_trace_peaks_second.tsv"
OUT = BASE / "isotope_spacing_diagnostic.tsv"

SPECTRUM = "szwk021704i19101xms1.20054.20054.4"
FEATURE_ID = "10796541564565835142"
SCAN_INDEX = 20060

ISOTOPE_SPACING = 1.00335483507
CHARGE = 4
PSM_MZ = 969.4884
NOMINATED_OFFSET = 2
MAX_ISOTOPE_PEAK = 4
PPM_TOL = 10.0


def main():
    trace = pd.read_csv(TRACE_FILE, sep="\t", low_memory=False)
    trace["feature_id"] = trace["feature_id"].astype(str)

    scan = trace[
        (trace["feature_id"] == FEATURE_ID)
        & (pd.to_numeric(trace["scan_index"], errors="coerce") == SCAN_INDEX)
    ].copy()

    if scan.empty:
        raise SystemExit("No trace rows found for the requested feature and scan.")

    scan["observed_mz"] = pd.to_numeric(scan["observed_mz"], errors="coerce")
    scan = scan.dropna(subset=["observed_mz"]).sort_values("observed_mz")
    observed_mzs = scan["observed_mz"].to_numpy(float)

    step = ISOTOPE_SPACING / CHARGE
    anchor_mz = PSM_MZ - NOMINATED_OFFSET * step

    rows = []
    matched = []

    for k in range(MAX_ISOTOPE_PEAK + 1):
        theoretical_mz = anchor_mz + k * step
        ppm_errors = (
            (observed_mzs - theoretical_mz) / theoretical_mz * 1e6
        )
        idx = int(np.argmin(np.abs(ppm_errors)))
        observed_mz = float(observed_mzs[idx])
        ppm_error = float(ppm_errors[idx])
        accepted = abs(ppm_error) <= PPM_TOL

        matched.append(observed_mz if accepted else np.nan)

        rows.append({
            "spectrum": SPECTRUM,
            "feature_id": FEATURE_ID,
            "scan_index": SCAN_INDEX,
            "charge": CHARGE,
            "nominated_offset": NOMINATED_OFFSET,
            "isotope_index": k,
            "theoretical_mz": theoretical_mz,
            "nearest_observed_mz": observed_mz,
            "individual_peak_error_ppm": ppm_error,
            "within_tolerance": accepted,
        })

    for k in range(MAX_ISOTOPE_PEAK):
        m1, m2 = matched[k], matched[k + 1]
        if np.isfinite(m1) and np.isfinite(m2):
            observed_spacing = m2 - m1
            spacing_error_ppm = (
                (observed_spacing - step) / step * 1e6
            )
        else:
            observed_spacing = np.nan
            spacing_error_ppm = np.nan

        rows[k]["next_observed_spacing_mz"] = observed_spacing
        rows[k]["adjacent_spacing_error_ppm"] = spacing_error_ppm

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)

    print(out.to_string(index=False))
    print(f"\nExpected isotope spacing: {step:.9f} m/z")
    print(f"Diagnostic written to: {OUT}")


if __name__ == "__main__":
    main()
