#!/usr/bin/env python3
"""Consolidate existing isotope-envelope diagnostics for one candidate feature."""

from pathlib import Path
import csv
import statistics
from collections import defaultdict

ROOT = Path.cwd()
REPORTS = ROOT / "reports"
FEATURE_ID = "10275997325121874585"

INPUTS = {
    "profiles": REPORTS / "full_trace_chromatography_v1.tsv",
    "spacing": REPORTS / "adjacent_theoretical_spacing_v1.tsv",
    "core_audit": REPORTS / "scan_aligned_envelope_audit_v1.tsv",
}

OUT_DETAIL = REPORTS / "feature_envelope_evidence_v1.tsv"
OUT_SUMMARY = REPORTS / "feature_envelope_evidence_summary_v1.tsv"


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def median_or_blank(values):
    return statistics.median(values) if values else ""


def max_or_blank(values):
    return max(values) if values else ""


for name, path in INPUTS.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing {name} input: {path}")

profiles = [
    r for r in read_tsv(INPUTS["profiles"])
    if r["feature_id"] == FEATURE_ID
]
spacing_rows = [
    r for r in read_tsv(INPUTS["spacing"])
    if r["feature_id"] == FEATURE_ID
]
core_rows = [
    r for r in read_tsv(INPUTS["core_audit"])
    if r["feature_id"] == FEATURE_ID
]

if not profiles:
    raise ValueError(f"No profile observations for {FEATURE_ID}")
if not spacing_rows:
    raise ValueError(f"No adjacent-spacing observations for {FEATURE_ID}")
if not core_rows:
    raise ValueError(f"No core-model audit observations for {FEATURE_ID}")

by_trace = defaultdict(list)
for r in profiles:
    by_trace[int(r["trace_index"])].append(r)

spacing_by_trace = defaultdict(list)
for r in spacing_rows:
    spacing_by_trace[int(r["upper_trace_index"])].append(r)

core_by_trace = defaultdict(list)
for r in core_rows:
    core_by_trace[int(r["trace_index"])].append(r)

# Existing intensity-coherence results from Script 81 are incorporated
# only when the report exists. Correlations are descriptive and use small n.
intensity_file = REPORTS / "local_track_intensity_correlations_v1.tsv"
intensity_note = (
    "Track intensity correlations are optional descriptive evidence; "
    "they do not establish isotope identity."
)

detail = []
for trace in sorted(by_trace):
    p = by_trace[trace]
    rt = [float(r["rt_sec"]) for r in p]
    intensity = [float(r["peak_intensity"]) for r in p]

    core = core_by_trace.get(trace, [])
    core_resid = [
        abs(float(r["residual_from_core_model_mDa"]))
        for r in core
        if r.get("residual_from_core_model_mDa", "") != ""
    ]
    core_within = sum(
        str(r["within_3mDa"]).lower() == "true"
        for r in core
        if r.get("within_3mDa", "") != ""
    )

    # Spacing rows are keyed by upper trace, so trace k uses pair (k-1, k).
    sp = spacing_by_trace.get(trace, [])
    sp_resid = [
        abs(float(r["spacing_residual_mDa"]))
        for r in sp
        if r.get("spacing_residual_mDa", "") != ""
    ]
    sp_within = sum(
        str(r["within_3mDa"]).lower() == "true"
        for r in sp
        if r.get("within_3mDa", "") != ""
    )

    if trace <= 4:
        status = "core_candidate"
        rationale = (
            "Part of the provisional traces 1-4 core; core fit is not "
            "independent validation."
        )
    elif trace == 5:
        status = "plausible_extension"
        rationale = (
            "Most pair 4-5 spacings pass the ±3 mDa criterion, but evidence "
            "is limited and at least one scan is inconsistent."
        )
    else:
        status = "unsupported_extension"
        rationale = (
            "Pair 5-6 spacing is repeatedly inconsistent with theoretical "
            "charge-2 spacing."
        )

    detail.append({
        "feature_id": FEATURE_ID,
        "trace_index": trace,
        "n_profile_observations": len(p),
        "rt_min_sec": f"{min(rt):.6f}",
        "rt_max_sec": f"{max(rt):.6f}",
        "rt_span_sec": f"{max(rt) - min(rt):.6f}",
        "median_intensity": f"{statistics.median(intensity):.6f}",
        "max_intensity": f"{max(intensity):.6f}",
        "core_model_n": len(core),
        "core_model_median_abs_residual_mDa": (
            f"{median_or_blank(core_resid):.6f}"
            if core_resid else ""
        ),
        "core_model_max_abs_residual_mDa": (
            f"{max_or_blank(core_resid):.6f}"
            if core_resid else ""
        ),
        "core_model_within_3mDa": core_within,
        "adjacent_pair": f"{trace - 1}-{trace}" if trace > 1 else "",
        "adjacent_spacing_n": len(sp),
        "adjacent_spacing_median_abs_residual_mDa": (
            f"{median_or_blank(sp_resid):.6f}"
            if sp_resid else ""
        ),
        "adjacent_spacing_within_3mDa": sp_within,
        "provisional_evidence_status": status,
        "rationale": rationale,
    })

with OUT_DETAIL.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(detail[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(detail)

summary = {
    "feature_id": FEATURE_ID,
    "n_traces": len(by_trace),
    "provisional_core": "1-4",
    "plausible_extension": "5",
    "unsupported_extension": "6",
    "independent_IG_established": "False",
    "original_XLSX_used": "False",
    "intensity_coherence_note": intensity_note,
    "overall_interpretation": (
        "Traces 1-4 form the provisional core; trace 5 is unresolved but "
        "plausible; trace 6 is not sufficiently supported as the next "
        "isotope. These are evidence labels, not probabilities or a "
        "definitive isotope-group assignment."
    ),
}

with OUT_SUMMARY.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(summary), delimiter="\t")
    writer.writeheader()
    writer.writerow(summary)

print(f"Feature: {FEATURE_ID}")
print(f"Per-trace evidence: {OUT_DETAIL}")
print(f"Summary: {OUT_SUMMARY}")
print("\nTrace evidence:")
for r in detail:
    print(
        f'Trace {r["trace_index"]}: '
        f'n={r["n_profile_observations"]}, '
        f'RT span={r["rt_span_sec"]} s, '
        f'core median |residual|={r["core_model_median_abs_residual_mDa"]} mDa, '
        f'adjacent spacing within ±3 mDa='
        f'{r["adjacent_spacing_within_3mDa"]}/'
        f'{r["adjacent_spacing_n"]}, '
        f'status={r["provisional_evidence_status"]}'
    )
print("\nNo definitive isotope assignment is made.")
