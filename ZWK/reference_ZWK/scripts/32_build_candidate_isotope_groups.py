#!/usr/bin/env python3
"""Build provisional isotope-group candidates from OpenMS traces.

This is a candidate-generation procedure, not validated monoisotopic assignment.
The original XLSX isotope-group assignments are not used.
"""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("openms_pilot")
INPUT = BASE / "trace_structure_audit_v2.tsv"
DETAIL_OUT = BASE / "candidate_isotope_groups_v2.tsv"
SUMMARY_OUT = BASE / "candidate_isotope_groups_summary_v2.tsv"

ISOTOPE_SPACING = 1.00335483507
MAX_ISOTOPE_STEPS = 4
GAP_TOLERANCE_MDA = 2.0
MIN_TRACE_SCANS = 5


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.rank[a] < self.rank[b]:
            a, b = b, a
        self.parent[b] = a
        if self.rank[a] == self.rank[b]:
            self.rank[a] += 1


def main():
    d = pd.read_csv(INPUT, sep="\t")
    required = {
        "acquisition", "feature_id", "trace_index", "feature_charge",
        "mz_median", "n_scans", "intensity_sum",
    }
    missing = required - set(d.columns)
    if missing:
        raise ValueError(f"Input lacks columns: {sorted(missing)}")

    detail = []
    summary = []

    group_keys = ["acquisition", "feature_id", "feature_charge"]

    for (acquisition, feature_id, charge), g in d.groupby(
        group_keys, sort=False
    ):
        charge = int(charge)
        if charge <= 0:
            continue

        # Keep only traces persistent across several scans.
        g = g[g["n_scans"] >= MIN_TRACE_SCANS].copy()
        if len(g) < 2:
            continue

        g = g.sort_values("mz_median").reset_index(drop=True)
        mz = g["mz_median"].to_numpy(dtype=float)
        n = len(g)
        expected_step = ISOTOPE_SPACING / charge

        uf = UnionFind(n)
        accepted_edges = []

        for i in range(n - 1):
            for j in range(i + 1, n):
                gap = mz[j] - mz[i]
                if gap <= 0:
                    continue

                steps = int(round(gap / expected_step))
                if not 1 <= steps <= MAX_ISOTOPE_STEPS:
                    continue

                deviation_mDa = (
                    gap - steps * expected_step
                ) * 1000.0

                if abs(deviation_mDa) <= GAP_TOLERANCE_MDA:
                    uf.union(i, j)
                    accepted_edges.append({
                        "i": i,
                        "j": j,
                        "steps": steps,
                        "gap_mDa": deviation_mDa,
                    })

        components = {}
        for i in range(n):
            components.setdefault(uf.find(i), []).append(i)

        group_number = 0
        for indices in components.values():
            if len(indices) < 2:
                continue

            group_number += 1
            group_id = (
                f"{acquisition}:{feature_id}:{charge}:C{group_number}"
            )
            component_edges = [
                e for e in accepted_edges
                if e["i"] in indices and e["j"] in indices
            ]

            # A connected component can contain inconsistent or branching
            # links. Flag it instead of silently treating it as a clean ladder.
            trace_mz = mz[indices]
            adjacent_gaps = np.diff(trace_mz)
            adjacent_steps = np.rint(
                adjacent_gaps / expected_step
            ).astype(int)

            component_edge_errors = [
                abs(e["gap_mDa"]) for e in component_edges
            ]

            clean_ladder = (
                len(indices) >= 2
                and len(adjacent_steps) >= 1
                and np.all(adjacent_steps >= 1)
                and np.all(adjacent_steps <= MAX_ISOTOPE_STEPS)
                and np.all(np.abs(
                    adjacent_gaps
                    - adjacent_steps * expected_step
                ) * 1000 <= GAP_TOLERANCE_MDA)
            )

            if clean_ladder:
                status = "candidate_regular_ladder"
            else:
                status = "candidate_ambiguous_component"

            for i in indices:
                r = g.iloc[i]
                detail.append({
                    "acquisition": acquisition,
                    "candidate_ig_id": group_id,
                    "feature_id": feature_id,
                    "charge": charge,
                    "trace_index": r["trace_index"],
                    "trace_mz_median": r["mz_median"],
                    "trace_mz_min": r["mz_min"],
                    "trace_mz_max": r["mz_max"],
                    "trace_n_scans": r["n_scans"],
                    "trace_intensity_sum": r["intensity_sum"],
                    "component_trace_count": len(indices),
                    "component_edge_count": len(component_edges),
                    "median_abs_component_edge_deviation_mDa": (
                        np.median(component_edge_errors)
                        if component_edge_errors else np.nan
                    ),
                    "candidate_status": status,
                    "monoisotopic_anchor_assigned": False,
                })

            summary.append({
                "acquisition": acquisition,
                "candidate_ig_id": group_id,
                "feature_id": feature_id,
                "charge": charge,
                "n_traces": len(indices),
                "n_accepted_edges": len(component_edges),
                "mz_min": float(np.min(trace_mz)),
                "mz_max": float(np.max(trace_mz)),
                "mz_span": float(np.max(trace_mz) - np.min(trace_mz)),
                "median_abs_edge_deviation_mDa": (
                    np.median(component_edge_errors)
                    if component_edge_errors else np.nan
                ),
                "candidate_status": status,
                "monoisotopic_anchor_assigned": False,
            })

    detail_df = pd.DataFrame(detail)
    summary_df = pd.DataFrame(summary)

    detail_df.to_csv(DETAIL_OUT, sep="\t", index=False)
    summary_df.to_csv(SUMMARY_OUT, sep="\t", index=False)

    print(f"Candidate trace table: {DETAIL_OUT}")
    print(f"Candidate group summary: {SUMMARY_OUT}")
    print(f"Candidate groups: {len(summary_df)}")
    print(f"Trace memberships: {len(detail_df)}")

    if not summary_df.empty:
        print("\n=== Candidate group status ===")
        print(summary_df["candidate_status"].value_counts().to_string())

        print("\n=== Number of traces per candidate group ===")
        print(summary_df["n_traces"].describe().to_string())

        print("\n=== Groups by acquisition ===")
        print(
            summary_df.groupby(
                ["acquisition", "candidate_status"]
            ).size().rename("groups").to_string()
        )


if __name__ == "__main__":
    main()
