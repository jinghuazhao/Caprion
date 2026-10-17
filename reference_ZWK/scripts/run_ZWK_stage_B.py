#!/usr/bin/env python3
"""Stage B: link reconstructed ZWK isotope groups to FragPipe PSM evidence.

Candidate links are based on sample, charge, precursor m/z and retention time.
They are not confirmed identifications or guaranteed MS2-to-MS1 scan links.
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

PROTON_MASS = 1.007276466621


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--igs", type=Path,
        default=Path("reports/ZWK_merged_reconstructed_igs.tsv"),
    )
    p.add_argument(
        "--peaks", type=Path,
        default=Path("reports/ZWK_merged_ig_peaks.tsv"),
    )
    p.add_argument("--psm-glob", default="search_98/*/psm.tsv")
    p.add_argument(
        "--outdir", type=Path,
        default=Path("reports/ZWK_stage_B"),
    )
    p.add_argument("--ppm", type=float, default=10.0)
    p.add_argument("--rt-tol", type=float, default=10.0)
    p.add_argument("--chunk-size", type=int, default=250_000)
    return p.parse_args()


def sample_from_path(path):
    return path.parent.name.strip()


def scan_from_spectrum(value):
    # Example: samplexms1.00488.00488.2
    m = re.search(r"\.(\d+)\.\d+\.\d+$", str(value))
    return int(m.group(1)) if m else None


def scan_from_native_id(value):
    m = re.search(r"(?:^|\s)scan=(\d+)", str(value))
    return int(m.group(1)) if m else None


def numeric_column(df, name):
    if name not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[name], errors="coerce")


def load_psms(path):
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df.columns = df.columns.astype(str).str.strip()

    required = {"Spectrum", "Charge", "Retention", "Peptide"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing PSM columns: {sorted(missing)}")

    charge = numeric_column(df, "Charge")

    # Prefer explicit m/z fields if present. Otherwise convert neutral mass
    # to m/z using (M + z * proton_mass) / z.
    mz = numeric_column(df, "Calibrated Observed M/Z")
    mz = mz.fillna(numeric_column(df, "Observed M/Z"))

    calibrated_mass = numeric_column(df, "Calibrated Observed Mass")
    observed_mass = numeric_column(df, "Observed Mass")

    valid_charge = charge.gt(0)
    derived_calibrated = (calibrated_mass + charge * PROTON_MASS) / charge
    derived_observed = (observed_mass + charge * PROTON_MASS) / charge

    mz = mz.fillna(derived_calibrated.where(valid_charge))
    mz = mz.fillna(derived_observed.where(valid_charge))

    result = pd.DataFrame({
        "psm_file": str(path),
        "psm_row_number": np.arange(2, len(df) + 2),
        "spectrum": df["Spectrum"].fillna(""),
        "psm_scan_number": df["Spectrum"].map(scan_from_spectrum),
        "charge": charge,
        "psm_rt_seconds": numeric_column(df, "Retention"),
        "psm_mz": mz,
        "peptide": df["Peptide"].fillna(""),
    })

    for col, target in (
        ("Modified Peptide", "modified_peptide"),
        ("Protein", "protein"),
        ("Gene", "gene"),
        ("Probability", "psm_probability"),
        ("Expectation", "psm_expectation"),
    ):
        if col in df:
            result[target] = df[col]

    result["sample_id"] = sample_from_path(path)
    result = result.dropna(
        subset=["charge", "psm_rt_seconds", "psm_mz"]
    ).copy()
    result = result[result["charge"] > 0].copy()
    result["charge"] = result["charge"].astype(int)
    return result


def partition_igs(path, tmpdir, chunk_size):
    required = {
        "reconstructed_ig_id", "sample_id", "source_file", "charge",
        "precursor_mz_estimate", "rt_start_seconds", "rt_end_seconds",
        "rt_median_seconds", "n_scans",
    }
    counts = Counter()
    written = set()

    for chunk in pd.read_csv(
        path, sep="\t", chunksize=chunk_size, low_memory=False,
        dtype={
            "reconstructed_ig_id": str,
            "sample_id": str,
            "source_file": str,
        },
    ):
        chunk.columns = chunk.columns.astype(str).str.strip()
        missing = required - set(chunk.columns)
        if missing:
            raise ValueError(f"{path}: missing IG columns: {sorted(missing)}")

        chunk = chunk[list(required)].copy()
        chunk["sample_id"] = chunk["sample_id"].str.strip()
        chunk["charge"] = pd.to_numeric(chunk["charge"], errors="coerce")
        chunk["precursor_mz_estimate"] = pd.to_numeric(
            chunk["precursor_mz_estimate"], errors="coerce"
        )
        chunk["rt_median_seconds"] = pd.to_numeric(
            chunk["rt_median_seconds"], errors="coerce"
        )
        chunk = chunk.dropna(
            subset=["charge", "precursor_mz_estimate", "rt_median_seconds"]
        )
        chunk = chunk[chunk["charge"] > 0].copy()
        chunk["charge"] = chunk["charge"].astype(int)

        for sample, group in chunk.groupby("sample_id", sort=False):
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", sample)
            dest = tmpdir / f"{safe}.tsv"
            group.to_csv(
                dest, sep="\t", index=False, mode="a",
                header=sample not in written,
            )
            written.add(sample)
            counts[sample] += len(group)

    return counts


def build_scan_index(path, chunk_size):
    """Index distinct scan numbers represented by reconstructed IG peaks."""
    scans = {}

    for chunk in pd.read_csv(
        path, sep="\t", chunksize=chunk_size, dtype=str,
        usecols=lambda c: c.strip() in {"sample_id", "spectrum_id"},
    ):
        chunk.columns = chunk.columns.astype(str).str.strip()
        for sample, spectrum in chunk[["sample_id", "spectrum_id"]].dropna().itertuples(
            index=False, name=None
        ):
            scan = scan_from_native_id(spectrum)
            if scan is not None:
                scans.setdefault(sample.strip(), set()).add(scan)

    return scans


def main():
    args = arguments()
    args.outdir.mkdir(parents=True, exist_ok=True)

    psm_paths = sorted(Path(p) for p in glob.glob(args.psm_glob))
    if not psm_paths:
        raise SystemExit(f"No PSM files matched: {args.psm_glob}")

    all_psms = []
    psm_counts = Counter()
    for path in psm_paths:
        psms = load_psms(path)
        all_psms.append(psms)
        psm_counts[sample_from_path(path)] += len(psms)

    psms = pd.concat(all_psms, ignore_index=True)
    psms["psm_record_id"] = (
        psms["psm_file"] + "::" + psms["psm_row_number"].astype(str)
    )

    peak_scan_index = build_scan_index(args.peaks, args.chunk_size)

    candidate_rows = []
    sample_psm_stats = {}
    with tempfile.TemporaryDirectory(prefix="zwk_stage_b_") as tmp:
        tmpdir = Path(tmp)
        ig_counts = partition_igs(args.igs, tmpdir, args.chunk_size)

        for sample, sample_psms in psms.groupby("sample_id", sort=False):
            ig_path = tmpdir / (
                re.sub(r"[^A-Za-z0-9_.-]", "_", sample) + ".tsv"
            )
            if not ig_path.exists():
                sample_psm_stats[sample] = {
                    "sample_id": sample,
                    "n_usable_psms": len(sample_psms),
                    "n_psms_with_candidate": 0,
                    "n_psms_without_candidate": len(sample_psms),
                    "n_candidate_links": 0,
                    "n_igs": 0,
                    "n_igs_with_candidate": 0,
                    "n_psms_with_multiple_candidates": 0,
                }
                continue

            igs = pd.read_csv(ig_path, sep="\t", low_memory=False)
            n_igs = len(igs)
            igs["charge"] = igs["charge"].astype(int)

            link_count_by_psm = Counter()
            linked_igs = set()
            links = 0

            # Process each charge independently.
            for charge, charge_psms in sample_psms.groupby("charge", sort=False):
                charge_igs = igs[igs["charge"] == charge]
                if charge_igs.empty:
                    continue

                charge_igs = charge_igs.sort_values("precursor_mz_estimate")
                ig_mz = charge_igs["precursor_mz_estimate"].to_numpy(float)
                ig_rt = charge_igs["rt_median_seconds"].to_numpy(float)
                ig_records = charge_igs.to_dict("records")

                for p in charge_psms.itertuples(index=False):
                    tol = abs(p.psm_mz) * args.ppm * 1e-6
                    left = np.searchsorted(ig_mz, p.psm_mz - tol, side="left")
                    right = np.searchsorted(ig_mz, p.psm_mz + tol, side="right")
                    if right <= left:
                        continue

                    possible = np.arange(left, right)
                    rt_error = ig_rt[possible] - p.psm_rt_seconds
                    keep = np.abs(rt_error) <= args.rt_tol

                    for idx, rte in zip(possible[keep], rt_error[keep]):
                        ig = ig_records[idx]
                        mz_error = ig["precursor_mz_estimate"] - p.psm_mz
                        ppm_error = mz_error / p.psm_mz * 1e6

                        candidate_rows.append({
                            "sample_id": sample,
                            "reconstructed_ig_id": ig["reconstructed_ig_id"],
                            "source_file": ig["source_file"],
                            "charge": charge,
                            "ig_mz": ig["precursor_mz_estimate"],
                            "ig_rt_median_seconds": ig["rt_median_seconds"],
                            "psm_file": p.psm_file,
                            "psm_row_number": p.psm_row_number,
                            "psm_record_id": p.psm_record_id,
                            "spectrum": p.spectrum,
                            "psm_scan_number": p.psm_scan_number,
                            "psm_rt_seconds": p.psm_rt_seconds,
                            "psm_mz": p.psm_mz,
                            "peptide": p.peptide,
                            "modified_peptide": getattr(p, "modified_peptide", ""),
                            "protein": getattr(p, "protein", ""),
                            "gene": getattr(p, "gene", ""),
                            "psm_probability": getattr(p, "psm_probability", ""),
                            "psm_expectation": getattr(p, "psm_expectation", ""),
                            "mz_error_da": mz_error,
                            "mz_error_ppm": ppm_error,
                            "rt_error_seconds": rte,
                            "scan_number_present_in_A_peak_table": (
                                p.psm_scan_number is not None
                                and p.psm_scan_number in peak_scan_index.get(sample, set())
                            ),
                            "link_status": "candidate_not_confirmed",
                        })
                        link_count_by_psm[p.psm_record_id] += 1
                        linked_igs.add(ig["reconstructed_ig_id"])
                        links += 1

            n_with = sum(v > 0 for v in link_count_by_psm.values())
            sample_psm_stats[sample] = {
                "sample_id": sample,
                "n_usable_psms": len(sample_psms),
                "n_psms_with_candidate": n_with,
                "n_psms_without_candidate": len(sample_psms) - n_with,
                "n_candidate_links": links,
                "n_igs": n_igs,
                "n_igs_with_candidate": len(linked_igs),
                "n_psms_with_multiple_candidates": sum(
                    v > 1 for v in link_count_by_psm.values()
                ),
            }

    candidates_path = args.outdir / "ZWK_stage_B_IG_PSM_candidates.tsv"
    summary_path = args.outdir / "ZWK_stage_B_sample_summary.tsv"
    run_summary_path = args.outdir / "ZWK_stage_B_run_summary.json"

    columns = [
        "sample_id", "reconstructed_ig_id", "source_file", "charge",
        "ig_mz", "ig_rt_median_seconds", "psm_file", "psm_row_number",
        "psm_record_id", "spectrum", "psm_scan_number", "psm_rt_seconds",
        "psm_mz", "peptide", "modified_peptide", "protein", "gene",
        "psm_probability", "psm_expectation", "mz_error_da", "mz_error_ppm",
        "rt_error_seconds", "scan_number_present_in_A_peak_table",
        "link_status",
    ]
    pd.DataFrame(candidate_rows, columns=columns).to_csv(
        candidates_path, sep="\t", index=False
    )
    pd.DataFrame(list(sample_psm_stats.values())).sort_values(
        "sample_id"
    ).to_csv(summary_path, sep="\t", index=False)

    summary = {
        "psm_files_found": len(psm_paths),
        "usable_psm_rows": len(psms),
        "usable_ig_rows": int(sum(ig_counts.values())),
        "samples_with_igs": len(ig_counts),
        "ppm_tolerance": args.ppm,
        "rt_tolerance_seconds": args.rt_tol,
        "candidate_links": len(candidate_rows),
        "unique_igs_with_candidate": len({
            r["reconstructed_ig_id"] for r in candidate_rows
        }),
        "unique_psm_records_with_candidate": len({
            r["psm_record_id"] for r in candidate_rows
        }),
        "candidate_links_are_confirmed_identifications": False,
        "scan_number_flag_is_direct_ms2_ms1_link": False,
        "psm_mz_source": (
            "explicit m/z if present; otherwise calibrated observed mass "
            "converted to m/z; otherwise observed mass converted to m/z"
        ),
        "outputs": {
            "candidates": str(candidates_path),
            "sample_summary": str(summary_path),
        },
    }
    run_summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary, indent=2))
    print(f"Sample summary: {summary_path}")


if __name__ == "__main__":
    main()
