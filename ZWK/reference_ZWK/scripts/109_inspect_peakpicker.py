#!/usr/bin/env python3
"""Inspect PyOpenMS high-resolution peak picker without processing mzML files."""

import pyopenms as oms

print("PyOpenMS:", oms.__version__)

picker = oms.PeakPickerHiRes()

print("\n=== PeakPickerHiRes parameters ===")
params = picker.getParameters()
print(params)

print("\n=== Public methods ===")
for name in dir(picker):
    if not name.startswith("_"):
        print(name)

print("\n=== Method documentation ===")
for name in ("pick", "pickExperiment", "pickSpectrum", "getParameters",
             "setParameters", "getDefaults"):
    method = getattr(picker, name, None)
    if method is not None:
        print(f"\n--- {name} ---")
        print(getattr(method, "__doc__", None) or "(no docstring)")

print("\n=== Relevant data classes ===")
for name in ("MSSpectrum", "MSExperiment", "MzMLFile",
             "Peak1D", "Feature", "FeatureMap", "MassTrace"):
    cls = getattr(oms, name, None)
    print(f"{name}: {'available' if cls is not None else 'not exposed'}")
    if cls is not None and name == "MSSpectrum":
        for method in dir(cls):
            if any(k in method.lower() for k in
                   ("peak", "precursor", "rt", "mslevel", "mz", "intensity")):
                if not method.startswith("_"):
                    print(" ", method)

print("\nInspection only; no data were loaded or modified.")
