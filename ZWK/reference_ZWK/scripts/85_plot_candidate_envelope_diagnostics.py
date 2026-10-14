#!/usr/bin/env python3
"""Plot mass-spacing residuals and chromatographic trace profiles."""

from pathlib import Path
import csv
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path.cwd()
REPORTS = ROOT / "reports"

RESIDUAL_FILE = REPORTS / "scan_aligned_envelope_audit_v1.tsv"
PROFILE_FILE = REPORTS / "full_trace_chromatography_v1.tsv"

FEATURE_ID = "10275997325121874585"
TRACE_INDICES = [1, 2, 3, 4, 5, 6]
SCANS = [5852, 5865, 5878, 5891, 5904, 5917]


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def as_float(value):
    return float(value) if value not in (None, "") else None


for path in (RESIDUAL_FILE, PROFILE_FILE):
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run the corresponding earlier script first."
        )

residual_rows = read_tsv(RESIDUAL_FILE)
profile_rows = read_tsv(PROFILE_FILE)

# Panel 1: mass residuals against the core-derived charge-2 envelope model.
fig, ax = plt.subplots(figsize=(10, 5.5))

for trace in TRACE_INDICES:
    rows = sorted(
        [
            r for r in residual_rows
            if int(r["trace_index"]) == trace
            and r["residual_from_core_model_mDa"] != ""
        ],
        key=lambda r: int(r["scan_index_zero_based"]),
    )
    if not rows:
        continue

    x = [int(r["scan_index_zero_based"]) for r in rows]
    y = [float(r["residual_from_core_model_mDa"]) for r in rows]
    ax.plot(x, y, marker="o", linewidth=1.2, label=f"Trace {trace}")

ax.axhline(0, linewidth=1, color="black")
ax.axhline(3, linestyle="--", linewidth=1, color="grey")
ax.axhline(-3, linestyle="--", linewidth=1, color="grey")
ax.set_xlabel("Zero-based mzML spectrum index")
ax.set_ylabel("Observed minus expected m/z (mDa)")
ax.set_title("Candidate isotope-envelope mass residuals")
ax.set_xticks(SCANS)
ax.tick_params(axis="x", rotation=30)
ax.grid(True, alpha=0.25)
ax.legend(ncol=3)
fig.tight_layout()

residual_png = REPORTS / "candidate_envelope_mass_residuals_v1.png"
fig.savefig(residual_png, dpi=180)
plt.close(fig)

# Panel 2: log10 intensity versus retention time.
fig, ax = plt.subplots(figsize=(10, 5.5))

for trace in TRACE_INDICES:
    rows = sorted(
        [
            r for r in profile_rows
            if int(r["trace_index"]) == trace
            and as_float(r["peak_intensity"]) is not None
            and float(r["peak_intensity"]) > 0
        ],
        key=lambda r: float(r["rt_sec"]),
    )
    if not rows:
        continue

    x = [float(r["rt_sec"]) for r in rows]
    y = [math.log10(float(r["peak_intensity"])) for r in rows]
    ax.plot(x, y, marker="o", markersize=3, linewidth=1.2,
            label=f"Trace {trace}")

ax.set_xlabel("Retention time (s)")
ax.set_ylabel("log10 peak intensity")
ax.set_title("Full chromatographic profiles of candidate traces")
ax.grid(True, alpha=0.25)
ax.legend(ncol=3)
fig.tight_layout()

profile_png = REPORTS / "candidate_envelope_chromatography_v1.png"
fig.savefig(profile_png, dpi=180)
plt.close(fig)

print(f"Feature: {FEATURE_ID}")
print(f"Mass-residual plot: {residual_png}")
print(f"Chromatographic-profile plot: {profile_png}")
print("\nThe mass-residual plot uses traces 1–4 to define the core model.")
print(
    "Caution: these plots are diagnostic. They do not establish isotope "
    "identity or an independent isotope group."
)
