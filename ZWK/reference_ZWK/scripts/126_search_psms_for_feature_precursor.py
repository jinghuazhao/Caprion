#!/usr/bin/env python3
"""Search FragPipe PSM tables for charge-2 precursors near m/z 1341.068."""

import csv
from pathlib import Path

TARGET_MZ = 1341.068
TARGET_CHARGE = 2
TARGET_RT = 764.23206

STRICT_TOLERANCE = 0.020
BROAD_TOLERANCE = 0.050
RT_WINDOW = 10.0

ROOT = Path("search_98")
OUT_ALL = Path(
    "reports/feature10275997325121874585_cross_acquisition_psm_candidates_v1.tsv"
)
OUT_NEAR = Path(
    "reports/feature10275997325121874585_cross_acquisition_psm_candidates_near_rt_v1.tsv"
)


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


tables = sorted(ROOT.rglob("psm.tsv"))
if not tables:
    raise FileNotFoundError(f"No psm.tsv tables found under {ROOT}")

candidates = []
files_read = 0
rows_read = 0
missing_headers = []

for path in tables:
    try:
        handle = path.open(newline="", errors="replace")
    except OSError as exc:
        print(f"WARNING: cannot open {path}: {exc}")
        continue

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Spectrum", "Charge", "Retention", "Observed M/Z"}

        if not required.issubset(set(reader.fieldnames or [])):
            missing_headers.append(str(path))
            continue

        files_read += 1

        for row in reader:
            rows_read += 1
            charge = num(row.get("Charge"))
            mz = num(row.get("Observed M/Z"))
            rt = num(row.get("Retention"))

            if charge is None or int(charge) != TARGET_CHARGE:
                continue
            if mz is None or abs(mz - TARGET_MZ) > BROAD_TOLERANCE:
                continue

            rt_delta = rt - TARGET_RT if rt is not None else None
            spectrum = (row.get("Spectrum") or "").strip()

            candidates.append({
                "PSM table": str(path),
                "Spectrum": spectrum,
                "Spectrum File": row.get("Spectrum File", ""),
                "Charge": int(charge),
                "Observed M/Z": mz,
                "m/z error from target (mDa)": round(
                    (mz - TARGET_MZ) * 1000, 3
                ),
                "Within strict ±0.020 Da": (
                    abs(mz - TARGET_MZ) <= STRICT_TOLERANCE
                ),
                "Retention (s)": rt if rt is not None else "",
                "RT difference from target (s)": (
                    round(rt_delta, 5) if rt_delta is not None else ""
                ),
                "Within ±10 s of target RT": (
                    abs(rt_delta) <= RT_WINDOW
                    if rt_delta is not None else False
                ),
                "Peptide": row.get("Peptide", ""),
                "Modified Peptide": row.get("Modified Peptide", ""),
                "Expectation": row.get("Expectation", ""),
                "Probability": row.get("Probability", ""),
                "Protein": row.get("Protein", ""),
                "Gene": row.get("Gene", ""),
            })

candidates.sort(
    key=lambda r: (
        abs(r["m/z error from target (mDa)"]),
        abs(r["RT difference from target (s)"])
        if isinstance(r["RT difference from target (s)"], (int, float))
        else float("inf"),
    )
)

fields = [
    "PSM table", "Spectrum", "Spectrum File", "Charge",
    "Observed M/Z", "m/z error from target (mDa)",
    "Within strict ±0.020 Da", "Retention (s)",
    "RT difference from target (s)", "Within ±10 s of target RT",
    "Peptide", "Modified Peptide", "Expectation", "Probability",
    "Protein", "Gene",
]

OUT_ALL.parent.mkdir(parents=True, exist_ok=True)

for output, rows in (
    (OUT_ALL, candidates),
    (OUT_NEAR, [
        row for row in candidates
        if row["Within ±10 s of target RT"] is True
    ]),
):
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

near_rt = [
    row for row in candidates
    if row["Within ±10 s of target RT"] is True
]
strict = [
    row for row in candidates
    if row["Within strict ±0.020 Da"] is True
]

print("PSM tables found:", len(tables))
print("PSM tables read:", files_read)
print("Rows read:", rows_read)
print("Tables skipped for missing required columns:", len(missing_headers))
print(f"Charge-{TARGET_CHARGE} PSMs within ±{BROAD_TOLERANCE:.3f} Da:",
      len(candidates))
print(f"Of these, within ±{STRICT_TOLERANCE:.3f} Da:", len(strict))
print(f"Of these, within ±{RT_WINDOW:.1f} s of target RT:", len(near_rt))
print("All precursor candidates:", OUT_ALL)
print("RT-near candidates:", OUT_NEAR)

for row in candidates[:30]:
    print(
        f"error={row['m/z error from target (mDa)']:+.3f} mDa; "
        f"RT delta={row['RT difference from target (s)']}; "
        f"spectrum={row['Spectrum']}; "
        f"peptide={row['Peptide']}; "
        f"modified={row['Modified Peptide']}"
    )

if not candidates:
    print(
        "\nNo charge-2 PSM in these tables has an observed precursor "
        "within ±0.050 Da of the target. This does not prove the feature "
        "has no peptide identity; it means these reported PSMs provide "
        "no such candidate under this mass window."
    )
