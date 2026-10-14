#!/usr/bin/env python3
"""Check candidate isotope-track overlap with recorded MS2 isolation windows."""

from pathlib import Path
import csv
from collections import Counter

INPUT = Path(
    "reports/feature10275997325121874585_ms2_precursor_evidence_v1.tsv"
)
OUT = Path(
    "reports/feature10275997325121874585_ms2_isolation_overlap_v1.tsv"
)

TRACKS = {
    2: 1341.067749,
    3: 1341.570496,
    4: 1342.068359,
    5: 1342.571045,
    6: 1343.074463,
    7: 1343.573608,
}

with INPUT.open(newline="") as f:
    records = list(csv.DictReader(f, delimiter="\t"))

out_rows = []
for r in records:
    if r["precursor_mz"] == "":
        continue

    lower = float(r["isolation_lower_mz"])
    upper = float(r["isolation_upper_mz"])
    covered = [
        track for track, mz in TRACKS.items()
        if lower <= mz <= upper
    ]

    out_rows.append({
        "scan_index": r["scan_index"],
        "rt": r["rt"],
        "rt_delta_s": r["rt_delta_s"],
        "precursor_mz": r["precursor_mz"],
        "precursor_charge": r["precursor_charge"],
        "isolation_lower_mz": lower,
        "isolation_upper_mz": upper,
        "isolation_width_da": upper - lower,
        "covered_track_ids": ",".join(map(str, covered)),
        "covers_core_tracks_2_5": int(
            all(t in covered for t in (2, 3, 4, 5))
        ),
        "covers_any_core_track": int(
            any(t in covered for t in (2, 3, 4, 5))
        ),
        "covers_track_6_or_7": int(
            any(t in covered for t in (6, 7))
        ),
    })

fields = list(out_rows[0]) if out_rows else []
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(out_rows)

print("Input records:", len(records))
print("Records with precursor metadata:", len(out_rows))
print("Output:", OUT)

if out_rows:
    print("\n=== Isolation-window overlap ===")
    print(
        "Windows covering any core track (2–5):",
        sum(r["covers_any_core_track"] for r in out_rows),
        "/", len(out_rows)
    )
    print(
        "Windows covering all core tracks (2–5):",
        sum(r["covers_core_tracks_2_5"] for r in out_rows),
        "/", len(out_rows)
    )
    print(
        "Windows covering track 6 or 7:",
        sum(r["covers_track_6_or_7"] for r in out_rows),
        "/", len(out_rows)
    )

    print("\n=== Precursor charge counts ===")
    for charge, count in sorted(Counter(
        r["precursor_charge"] for r in out_rows
    ).items()):
        print(f"Charge {charge or 'unknown'}: {count}")

    print("\n=== Most relevant overlapping windows ===")
    selected = [
        r for r in out_rows if r["covers_any_core_track"]
    ]
    for r in selected[:20]:
        print(
            f"scan={r['scan_index']} RT={float(r['rt']):.5f} "
            f"delta={float(r['rt_delta_s']):+.3f}s "
            f"precursor={float(r['precursor_mz']):.6f} "
            f"z={r['precursor_charge']} "
            f"window=[{r['isolation_lower_mz']:.6f},"
            f"{r['isolation_upper_mz']:.6f}] "
            f"tracks={r['covered_track_ids'] or 'none'}"
        )

print("\nNo mzML files modified; no IG assignments made.")
