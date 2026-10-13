#!/usr/bin/env python3
"""Audit historical and proposed ZWK search configurations without modifying them."""

from pathlib import Path
import hashlib
import re

BASE = Path.home() / "Caprion" / "analysis"
ZWK = BASE / "ZWK"
PIPELINE = BASE / "pipeline"
OUT = ZWK / "reference_ZWK" / "reports" / "01_search_config_audit.txt"

CONFIGS = {
    "historical": ZWK / "default-rev",
    "proposed": PIPELINE / "production_runs" / "production_98_20261011_152549",
}

FILES = {
    "workflow": ("fragpipe.workflow", "workflow_used.workflow"),
    "fragger parameters": ("fragger.params", "fragger.params"),
    "MSBooster parameters": ("msbooster_params.txt", "msbooster_params.txt"),
    "manifest": ("fragpipe-files.fp-manifest", "manifest_used.fp-manifest"),
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_text(path):
    return path.read_text(errors="replace").splitlines()


def relevant_lines(path):
    patterns = re.compile(
        r"database|fasta|manifest|workflow|ionquant\.mbr|"
        r"isotope_error|num_threads|protein|speclibgen|"
        r"rt-cal|iontype_[cz]|pinPepXMLDirectory|mzmlDirectory",
        re.IGNORECASE,
    )
    return [
        line.strip()
        for line in read_text(path)
        if patterns.search(line) and line.strip()
    ]


def main():
    report = [
        "ZWK search configuration audit",
        "=" * 72,
        "Read-only audit; original configurations are not modified.",
        "",
    ]

    resolved_files = {}

    for label, root in CONFIGS.items():
        report += [f"[{label.upper()}]", f"Directory: {root}", ""]

        if not root.is_dir():
            report += ["ERROR: directory not found", ""]
            continue

        for description, names in FILES.items():
            filename = names[0] if label == "historical" else names[1]
            path = root / filename

            report.append(f"{description}: {path}")
            if not path.is_file():
                report += ["  MISSING", ""]
                continue

            digest = sha256(path)
            report.append(f"  SHA256: {digest}")
            resolved_files[(label, description)] = path

            if description == "manifest":
                lines = [
                    line for line in read_text(path)
                    if line.strip() and not line.lstrip().startswith("#")
                ]
                report.append(f"  Non-comment lines: {len(lines)}")

            if description in (
                "workflow", "fragger parameters", "MSBooster parameters"
            ):
                report.append("  Relevant settings:")
                for line in relevant_lines(path)[:80]:
                    report.append(f"    {line}")

            report.append("")

        report.append("Database references:")
        candidates = set()

        for key in (
            "workflow", "fragger parameters", "MSBooster parameters"
        ):
            path = resolved_files.get((label, key))
            if path is None:
                continue
            for line in read_text(path):
                if re.search(r"database|fasta", line, re.IGNORECASE):
                    for token in re.findall(r"(?<!\S)/[^ \t\"']+", line):
                        token = token.rstrip(",;)")
                        if token.lower().endswith((".fasta", ".fas", ".fa")):
                            candidates.add(Path(token))

        if not candidates:
            report.append("  No absolute FASTA paths detected.")
        for db in sorted(candidates, key=str):
            report.append(f"  {db} | exists={db.is_file()}")

        report.append("")

    report += [
        "INTERPRETATION NOTES",
        "-" * 72,
        "1. A matching file count does not establish equivalent search settings.",
        "2. Database existence does not establish biological suitability.",
        "3. Compare settings and database provenance before selecting a workflow.",
        "4. No search has been launched by this script.",
        "",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(report) + "\n")

    print(f"Audit report written to: {OUT}")
    print(f"Report lines: {len(report)}")
    print("No search was launched and no existing configuration was changed.")


if __name__ == "__main__":
    main()
