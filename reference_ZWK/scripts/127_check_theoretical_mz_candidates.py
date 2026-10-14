#!/usr/bin/env python3
"""Compare target precursor m/z against observed and theoretical PSM m/z."""

import csv
import heapq
from pathlib import Path

TARGET_MZ = 1341.068
TARGET_CHARGE = 2
TARGET_RT = 764.23206
TOLERANCE_DA = 0.050
N_NEAREST = 50

ROOT = Path("search_98")
REPORTS = Path("reports")

OBS_OUT = REPORTS / "feature10275997325121874585_observed_mz_candidates_v2.tsv"
CALC_OUT = REPORTS / "feature10275997325121874585_calculated_mz_candidates_v2.tsv"
NEAREST_OUT = REPORTS / "feature10275997325121874585_nearest_theoretical_mz_psms_v1.tsv"
RT_OUT = REPORTS / "feature10275997325121874585_calculated_mz_rt_ranked_v1.tsv"


def num(value):
    """Convert a value to a finite float, or return None."""
    try:
        value = float(value)
        return value if value == value and abs(value) != float("inf") else None
    except (TypeError, ValueError):
        return None


tables = sorted(ROOT.rglob("psm.tsv"))
observed_matches = []
calculated_matches = []
calculated_rt_matches = []
nearest_calculated = []

rows_read = 0
charge_rows = 0
skipped_tables = 0

required = {
    "Spectrum", "Charge", "Observed M/Z",
    "Calculated M/Z", "Retention", "Peptide",
}

for table in tables:
    try:
        with table.open(newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                skipped_tables += 1
                continue

            for row in reader:
                rows_read += 1

                if num(row.get("Charge")) != TARGET_CHARGE:
                    continue

                charge_rows += 1
                observed = num(row.get("Observed M/Z"))
                calculated = num(row.get("Calculated M/Z"))
                rt = num(row.get("Retention"))

                obs_error = observed - TARGET_MZ if observed is not None else None
                calc_error = calculated - TARGET_MZ if calculated is not None else None
                rt_delta = rt - TARGET_RT if rt is not None else None

                record = {
                    "PSM table": str(table),
                    "Spectrum": row.get("Spectrum", ""),
                    "Spectrum File": row.get("Spectrum File", ""),
                    "Charge": row.get("Charge", ""),
                    "Peptide": row.get("Peptide", ""),
                    "Modified Peptide": row.get("Modified Peptide", ""),
                    "Retention (s)": rt if rt is not None else "",
                    "RT delta from target (s)": round(rt_delta, 5) if rt_delta is not None else "",
                    "Absolute RT delta (s)": round(abs(rt_delta), 5) if rt_delta is not None else "",
                    "Observed M/Z": observed if observed is not None else "",
                    "Observed error (mDa)": round(obs_error * 1000, 3) if obs_error is not None else "",
                    "Observed error (ppm)": round(obs_error / TARGET_MZ * 1e6, 3) if obs_error is not None else "",
                    "Calculated M/Z": calculated if calculated is not None else "",
                    "Calculated error (mDa)": round(calc_error * 1000, 3) if calc_error is not None else "",
                    "Calculated error (ppm)": round(calc_error / TARGET_MZ * 1e6, 3) if calc_error is not None else "",
                    "Expectation": row.get("Expectation", ""),
                    "Probability": row.get("Probability", ""),
                }

                if obs_error is not None and abs(obs_error) <= TOLERANCE_DA:
                    observed_matches.append(record)

                if calc_error is None:
                    continue

                # Retain the 50 smallest absolute theoretical mass errors.
                heapq.heappush(
                    nearest_calculated,
                    (-abs(calc_error), rows_read, record),
                )
                if len(nearest_calculated) > N_NEAREST:
                    heapq.heappop(nearest_calculated)

                if abs(calc_error) <= TOLERANCE_DA:
                    calculated_matches.append(record)
                    if rt_delta is not None:
                        calculated_rt_matches.append(record)

    except OSError as exc:
        print(f"WARNING: could not read {table}: {exc}")
        skipped_tables += 1


REPORTS.mkdir(parents=True, exist_ok=True)

fields = [
    "PSM table", "Spectrum", "Spectrum File", "Charge",
    "Peptide", "Modified Peptide", "Retention (s)",
    "RT delta from target (s)", "Absolute RT delta (s)",
    "Observed M/Z", "Observed error (mDa)", "Observed error (ppm)",
    "Calculated M/Z", "Calculated error (mDa)",
    "Calculated error (ppm)", "Expectation", "Probability",
]


def write_table(path, records):
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(records)


observed_matches.sort(
    key=lambda r: abs(float(r["Observed error (mDa)"]))
)
calculated_matches.sort(
    key=lambda r: abs(float(r["Calculated error (mDa)"]))
)

nearest_records = [
    entry[2]
    for entry in sorted(
        nearest_calculated,
        key=lambda entry: abs(float(entry[2]["Calculated error (mDa)"])),
    )
]

calculated_rt_matches.sort(
    key=lambda r: (
        abs(float(r["RT delta from target (s)"])),
        abs(float(r["Calculated error (mDa)"])),
    )
)

write_table(OBS_OUT, observed_matches)
write_table(CALC_OUT, calculated_matches)
write_table(NEAREST_OUT, nearest_records)
write_table(RT_OUT, calculated_rt_matches)

print("PSM tables found:", len(tables))
print("PSM tables skipped:", skipped_tables)
print("Rows read:", rows_read)
print(f"Charge-{TARGET_CHARGE} rows:", charge_rows)
print(f"Target observed precursor m/z: {TARGET_MZ:.6f}")
print(f"Target charge: {TARGET_CHARGE}")
print(f"Target retention time: {TARGET_RT:.5f} s")
print(f"Tolerance: +/- {TOLERANCE_DA:.3f} Da")
print("Observed-m/z candidates:", len(observed_matches))
print("Calculated-m/z candidates:", len(calculated_matches))
print("Nearest theoretical candidates written:", len(nearest_records))
print("RT-ranked theoretical candidates:", len(calculated_rt_matches))
print("Observed-m/z report:", OBS_OUT)
print("Calculated-m/z report:", CALC_OUT)
print("Nearest-theoretical report:", NEAREST_OUT)
print("RT-ranked report:", RT_OUT)

print("\nNearest theoretical-mass candidates:")
for row in nearest_records[:10]:
    print(
        f"calculated_error={row['Calculated error (mDa)']} mDa; "
        f"observed_error={row['Observed error (mDa)']} mDa; "
        f"RT={row['Retention (s)']}; "
        f"RT_delta={row['RT delta from target (s)']} s; "
        f"peptide={row['Peptide']}; spectrum={row['Spectrum']}"
    )

print("\nBest retention-time matches within theoretical-mass tolerance:")
for row in calculated_rt_matches[:10]:
    print(
        f"RT_delta={row['RT delta from target (s)']} s; "
        f"calculated_error={row['Calculated error (mDa)']} mDa; "
        f"observed_error={row['Observed error (mDa)']} mDa; "
        f"peptide={row['Peptide']}; spectrum={row['Spectrum']}"
    )
