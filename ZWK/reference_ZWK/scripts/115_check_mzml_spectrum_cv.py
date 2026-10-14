#!/usr/bin/env python3
"""Inspect mzML CV terms for one spectrum using PyOpenMS metadata."""

from pathlib import Path
import pyopenms as oms

MZML = Path(
    "/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/"
    "analysis/ZWK/mzML/szwk900404i19201xms1.mzML"
)
RT_TARGET = 761.13018
RT_TOL = 0.01

exp = oms.MSExperiment()
oms.MzMLFile().load(str(MZML), exp)

matches = [
    (i, s) for i, s in enumerate(exp)
    if s.getMSLevel() == 1 and abs(s.getRT() - RT_TARGET) <= RT_TOL
]
if len(matches) != 1:
    raise RuntimeError(f"Expected one target spectrum; found {len(matches)}")

index, spectrum = matches[0]

print("File:", MZML)
print("Spectrum index:", index)
print("RT:", spectrum.getRT())
print("MS level:", spectrum.getMSLevel())
print("PyOpenMS spectrum type:", spectrum.getType())

print("\n=== Available spectrum metadata ===")
try:
    print(spectrum.getMetaInfo())
except Exception as exc:
    print("getMetaInfo unavailable:", repr(exc))

print("\n=== Data-array metadata ===")
try:
    print("Float data arrays:", spectrum.getFloatDataArrays())
except Exception as exc:
    print("Float data array query failed:", repr(exc))
try:
    print("Integer data arrays:", spectrum.getIntegerDataArrays())
except Exception as exc:
    print("Integer data array query failed:", repr(exc))
try:
    print("String data arrays:", spectrum.getStringDataArrays())
except Exception as exc:
    print("String data array query failed:", repr(exc))

print("\n=== Spectrum type enum members ===")
enum_cls = getattr(oms.SpectrumSettings, "SpectrumType", None)
if enum_cls is not None:
    for name in dir(enum_cls):
        if not name.startswith("_"):
            try:
                print(name, "=", getattr(enum_cls, name))
            except Exception:
                pass

print("\nNo files were modified.")
