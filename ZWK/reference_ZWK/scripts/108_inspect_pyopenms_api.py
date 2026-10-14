#!/usr/bin/env python3
"""Inspect installed PyOpenMS feature-detection capabilities without processing data."""

import inspect
import sys

import pyopenms as oms

print("=== Python environment ===")
print("Python:", sys.version.replace("\n", " "))
print("PyOpenMS version:", getattr(oms, "__version__", "not exposed"))
print("OpenMS version:", getattr(oms, "VersionInfo", None))

if hasattr(oms, "VersionInfo"):
    try:
        print("OpenMS version string:", oms.VersionInfo.getVersion())
    except Exception as exc:
        print("Version query unavailable:", repr(exc))

print("\n=== Relevant classes ===")
class_names = (
    "FeatureFinder",
    "FeatureFinderAlgorithmPicked",
    "FeatureFinderAlgorithmMetabo",
    "FeatureMap",
    "Feature",
    "MSExperiment",
    "MzMLFile",
    "PeakPickerHiRes",
    "MassTrace",
    "FeatureFinderMetabo",
)

for name in class_names:
    obj = getattr(oms, name, None)
    print(f"{name}: {'available' if obj is not None else 'not exposed'}")
    if obj is not None and name in (
        "FeatureFinder",
        "FeatureFinderAlgorithmMetabo",
        "PeakPickerHiRes",
    ):
        try:
            print("  doc:", (inspect.getdoc(obj) or "")[:500].replace("\n", " "))
        except Exception:
            pass

print("\n=== FeatureFinder methods ===")
if hasattr(oms, "FeatureFinder"):
    for name in dir(oms.FeatureFinder):
        if not name.startswith("_"):
            print(" ", name)

print("\n=== Algorithm names and default parameters ===")
if hasattr(oms, "FeatureFinder"):
    try:
        finder = oms.FeatureFinder()
        names = finder.getParameters()
        print("FeatureFinder parameters:")
        print(names)
    except Exception as exc:
        print("Could not query FeatureFinder parameters:", repr(exc))

    for algorithm in ("FeatureFinderMetabo", "FeatureFinderCentroided"):
        try:
            params = finder.getParameters(algorithm)
            print(f"\n{algorithm} parameters:")
            print(params)
        except Exception as exc:
            print(f"\n{algorithm}: parameter query failed: {exc!r}")

print("\n=== Feature and metadata access ===")
for cls_name in ("Feature", "FeatureMap", "MSExperiment", "MzMLFile"):
    cls = getattr(oms, cls_name, None)
    if cls is None:
        continue
    print(f"\n{cls_name}:")
    for name in dir(cls):
        if any(token in name.lower() for token in (
            "mz", "intensity", "charge", "quality", "rt",
            "load", "store", "feature", "meta", "size",
        )):
            if not name.startswith("_"):
                print(" ", name)

print("\nInspection only: no mzML files were loaded or modified.")
