#!/usr/bin/env python3
import csv
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

import numpy as np

REPORTS = Path("reports")
OUT = REPORTS / "ZWK_QC"
OUT.mkdir(exist_ok=True)

IG_FILE = REPORTS / "ZWK_merged_reconstructed_igs.tsv"
QUALITY_FILE = REPORTS / "ZWK_merged_ig_quality.tsv"
PEAK_FILE = REPORTS / "ZWK_merged_ig_peaks.tsv"


def rows(path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def num(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan


def describe(values):
    x = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if not len(x):
        return {"n": 0}
    return {
        "n": len(x),
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "min": float(np.min(x)),
        "q25": float(np.quantile(x, 0.25)),
        "median": float(np.median(x)),
        "q75": float(np.quantile(x, 0.75)),
        "max": float(np.max(x)),
    }


def write_tsv(path, records):
    records = list(records)
    if not records:
        print(f"No records to write: {path}")
        return
    fields = list(dict.fromkeys(k for r in records for k in r))
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(records)
    print(f"Wrote {path} ({len(records):,} rows)")


# ----------------------------------------------------------------------
# 1. Reconstructed IG table: overall, sample and source-file summaries
# ----------------------------------------------------------------------
print("Reading reconstructed IG table...")
igs = []
ig_ids = set()
sample_counts = Counter()
source_counts = Counter()
charge_counts = Counter()
sample_source_counts = Counter()

for r in rows(IG_FILE):
    gid = r["reconstructed_ig_id"]
    if gid in ig_ids:
        raise RuntimeError(f"Duplicate IG ID found: {gid}")
    ig_ids.add(gid)

    sample = r["sample_id"]
    source = Path(r["source_file"]).name
    charge = r["charge"].strip()

    sample_counts[sample] += 1
    source_counts[source] += 1
    charge_counts[charge] += 1
    sample_source_counts[(sample, source)] += 1

    igs.append({
        "reconstructed_ig_id": gid,
        "sample_id": sample,
        "source_file": source,
        "charge": charge,
        "precursor_mz_estimate": num(r["precursor_mz_estimate"]),
        "rt_start_seconds": num(r["rt_start_seconds"]),
        "rt_end_seconds": num(r["rt_end_seconds"]),
        "rt_median_seconds": num(r["rt_median_seconds"]),
        "n_scans": num(r["n_scans"]),
        "n_isotope_peak_observations": num(r["n_isotope_peak_observations"]),
        "median_abs_spacing_residual_mDa": num(r["median_abs_spacing_residual_mDa"]),
        "max_abs_spacing_residual_mDa": num(r["max_abs_spacing_residual_mDa"]),
    })

if not igs:
    raise RuntimeError("No reconstructed IG records found.")

metric_names = [
    "precursor_mz_estimate", "rt_start_seconds", "rt_end_seconds",
    "rt_median_seconds", "n_scans", "n_isotope_peak_observations",
    "median_abs_spacing_residual_mDa", "max_abs_spacing_residual_mDa",
]

overall = []
for metric in metric_names:
    d = describe([r[metric] for r in igs])
    overall.append({"section": "IG", "metric": metric, **d})

overall.extend([
    {"section": "count", "metric": "reconstructed_IGs", "n": len(igs)},
    {"section": "count", "metric": "samples", "n": len(sample_counts)},
    {"section": "count", "metric": "source_files", "n": len(source_counts)},
])
write_tsv(OUT / "overall_summary.tsv", overall)

write_tsv(OUT / "charge_distribution.tsv", [
    {"charge": k, "n_IGs": v, "fraction": v / len(igs)}
    for k, v in sorted(charge_counts.items())
])

write_tsv(OUT / "sample_IG_counts.tsv", [
    {"sample_id": k, "n_IGs": v}
    for k, v in sorted(sample_counts.items())
])

write_tsv(OUT / "source_file_IG_counts.tsv", [
    {"source_file": k, "n_IGs": v}
    for k, v in sorted(source_counts.items())
])

write_tsv(OUT / "sample_source_IG_counts.tsv", [
    {"sample_id": sample, "source_file": source, "n_IGs": n}
    for (sample, source), n in sorted(sample_source_counts.items())
])

# ----------------------------------------------------------------------
# 2. Quality table: validate IDs and summarise quality fields
# ----------------------------------------------------------------------
print("Reading quality table...")
quality_ids = set()
quality_metrics = defaultdict(list)
quality_notes = Counter()

for r in rows(QUALITY_FILE):
    gid = r["reconstructed_ig_id"]
    if gid not in ig_ids:
        raise RuntimeError(f"Quality record references unknown IG: {gid}")
    if gid in quality_ids:
        raise RuntimeError(f"Duplicate quality record for IG: {gid}")
    quality_ids.add(gid)

    for field in (
        "n_scans", "median_scan_gap_seconds", "min_peaks_per_scan",
        "median_peaks_per_scan", "max_peaks_per_scan",
        "median_abs_spacing_residual_mDa", "max_abs_spacing_residual_mDa",
    ):
        quality_metrics[field].append(num(r[field]))

    note = r.get("quality_note", "").strip() or "(blank)"
    quality_notes[note] += 1

if quality_ids != ig_ids:
    raise RuntimeError(
        f"Quality/IG ID mismatch: IGs={len(ig_ids):,}, "
        f"quality={len(quality_ids):,}"
    )

write_tsv(OUT / "quality_summary.tsv", [
    {"metric": field, **describe(values)}
    for field, values in quality_metrics.items()
])

write_tsv(OUT / "quality_notes.tsv", [
    {"quality_note": k, "n_IGs": v}
    for k, v in quality_notes.most_common()
])

# ----------------------------------------------------------------------
# 3. Peak table: stream every row; count peaks and isotope indices
# ----------------------------------------------------------------------
print("Streaming peak table; this can take several minutes...")
peak_counts = Counter()
peak_sample_counts = Counter()
isotope_counts = Counter()
peak_source_counts = Counter()
n_peaks = 0

with PEAK_FILE.open("r", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    required = {
        "reconstructed_ig_id", "sample_id", "source_file", "isotope_index"
    }
    if not reader.fieldnames or not required.issubset(reader.fieldnames):
        raise RuntimeError(f"Unexpected peak-table header: {reader.fieldnames}")

    for row_number, r in enumerate(reader, start=2):
        gid = r["reconstructed_ig_id"]
        if gid not in ig_ids:
            raise RuntimeError(
                f"Peak line {row_number} references unknown IG: {gid}"
            )

        n_peaks += 1
        peak_counts[gid] += 1
        peak_sample_counts[r["sample_id"]] += 1
        peak_source_counts[Path(r["source_file"]).name] += 1
        isotope_counts[r["isotope_index"].strip()] += 1

        if n_peaks % 25_000_000 == 0:
            print(f"  Processed {n_peaks:,} peak rows")

write_tsv(OUT / "peaks_per_IG_summary.tsv", [
    {"metric": "peak_rows_per_IG", **describe(list(peak_counts.values()))},
    {"metric": "IGs_with_at_least_one_peak", "n": len(peak_counts)},
    {"metric": "IGs_without_peak_rows", "n": len(ig_ids - peak_counts.keys())},
])

write_tsv(OUT / "isotope_index_distribution.tsv", [
    {"isotope_index": k, "n_peak_rows": v, "fraction": v / n_peaks}
    for k, v in sorted(
        isotope_counts.items(),
        key=lambda item: (0, int(item[0])) if item[0].lstrip("-").isdigit()
        else (1, item[0])
    )
])

write_tsv(OUT / "sample_peak_counts.tsv", [
    {"sample_id": k, "n_peak_rows": v}
    for k, v in sorted(peak_sample_counts.items())
])

write_tsv(OUT / "source_file_peak_counts.tsv", [
    {"source_file": k, "n_peak_rows": v}
    for k, v in sorted(peak_source_counts.items())
])

print("\nQC summary complete.")
print(f"IGs: {len(igs):,}")
print(f"Samples: {len(sample_counts):,}")
print(f"Source files: {len(source_counts):,}")
print(f"Peak rows streamed: {n_peaks:,}")
print(f"Output directory: {OUT}")
