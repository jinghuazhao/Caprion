#!/usr/bin/env python3
"""Read-only validation of ZWK reference-panel discovery inputs."""

from collections import Counter
from pathlib import Path
import sys

BASE = Path.home() / "Caprion" / "analysis"
ZWK = BASE / "ZWK"
MZML_DIR = ZWK / "mzML"
MANIFEST = BASE / "pipeline" / "production_v2.fp-manifest"
HISTORICAL = ZWK / "default-rev"


def main():
    errors = []

    print("ZWK reference-panel input check")
    print("=" * 40)
    print(f"Batch directory: {ZWK}")
    print(f"Source directory: {MZML_DIR}")
    print(f"Manifest: {MANIFEST}")
    print(f"Historical results: {HISTORICAL}")

    for label, path in (
        ("Batch directory", ZWK),
        ("Source directory", MZML_DIR),
        ("Pipeline directory", MANIFEST.parent),
        ("Manifest", MANIFEST),
        ("Historical results", HISTORICAL),
    ):
        print(f"{label} exists: {path.exists()}")
        if not path.exists():
            errors.append(f"{label} not found: {path}")

    if errors:
        print("\nFAIL: Required paths are missing.")
        for error in errors:
            print(f"  - {error}")
        return 1

    source_files = sorted(MZML_DIR.glob("*.mzML"))
    source = {p.resolve() for p in source_files}
    entries = []

    for lineno, line in enumerate(
        MANIFEST.read_text(errors="replace").splitlines(), start=1
    ):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 3:
            errors.append(f"Malformed manifest line {lineno}: {line}")
            continue

        entries.append(
            (Path(fields[0]).resolve(), fields[1], fields[2])
        )

    paths = [path for path, _, _ in entries]
    counts = Counter(paths)
    manifest_paths = set(paths)

    missing = sorted(path for path in manifest_paths if not path.is_file())
    duplicates = sorted(path for path, n in counts.items() if n > 1)
    absent = sorted(source - manifest_paths)
    outside = sorted(manifest_paths - source)

    print("\nCoverage")
    print("-" * 40)
    print(f"Source mzML files: {len(source)}")
    print(f"Manifest entries: {len(entries)}")
    print(f"Unique manifest paths: {len(manifest_paths)}")
    print(f"Missing paths on disk: {len(missing)}")
    print(f"Duplicate manifest paths: {len(duplicates)}")
    print(f"Source files absent from manifest: {len(absent)}")
    print(f"Manifest paths outside source set: {len(outside)}")

    print("\nManifest categories")
    print("-" * 40)
    for category, count in sorted(
        Counter(category for _, category, _ in entries).items()
    ):
        print(f"{category}: {count}")

    for label, paths_to_report in (
        ("Missing paths", missing),
        ("Duplicate paths", duplicates),
        ("Source files absent from manifest", absent),
        ("Manifest paths outside source set", outside),
    ):
        if paths_to_report:
            errors.append(f"{label}: {len(paths_to_report)}")
            for path in paths_to_report[:10]:
                print(f"  {path}")

    if len(source) != 98 or len(entries) != 98:
        errors.append("Expected exactly 98 source files and manifest entries.")

    print("\nHistorical results")
    print("-" * 40)
    print(f"Historical directory: {HISTORICAL}")
    print(
        "Historical interact pepXML files:",
        len(list(HISTORICAL.glob("interact-*.pep.xml"))),
    )
    print(
        "Historical original pepXML files:",
        len(list(HISTORICAL.glob("*.pepXML"))),
    )

    if errors:
        print("\nFAIL")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("\nPASS: all 98 source files match the manifest exactly once.")
    print("NOTE: this checks input coverage, not search-result validity.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
PY
