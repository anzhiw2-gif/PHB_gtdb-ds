#!/usr/bin/env python3
"""Apply frozen, conservative decision rules to the four-panel tree summary."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def support(row: dict[str, str]) -> float | None:
    value = row.get("candidate_reference_mrca_support", "")
    try:
        return float(value)
    except ValueError:
        return None


def decide(row: dict[str, str]) -> tuple[str, str]:
    panel = row["panel"]
    value = support(row)
    distance = float(row["nearest_reference_distance"])
    if panel == "extracellular_explicit":
        if value is not None and value >= 90:
            return "retain_extracellular_dPHASCL2_like_candidate", "nearest dPHASCL2 reference placement has support >=90; candidate-only"
        return "unresolved_extracellular_placement", "nearest reference is plausible but MRCA support is absent or weak"
    if panel == "phaC_vs_phaZ":
        return "PhaC_vs_PhaZ_unresolved", "challenge panel contains PhaC-like architecture; tree proximity to nPHAMCL references cannot establish depolymerase function"
    if panel == "structure_anomaly":
        return "cross_family_architecture_unresolved", "candidate assignment and nearest reference family disagree; support does not resolve the architecture conflict"
    if panel == "ephaz_competition":
        if value is not None and value >= 90 and distance <= 1.5:
            return "retain_dPHASCL1_like_candidate", "near dPHASCL1 reference placement with support >=90 and distance <=1.5; candidate-only"
        return "ePhaZ_competition_unresolved", "placement support or distance is insufficient for a resolved family call"
    raise ValueError(f"unknown panel: {panel}")


def derive(summary: Path, output: Path) -> dict[str, object]:
    with summary.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["phylo_decision"], row["phylo_reason"] = decide(row)
        row["phenotype_boundary"] = "candidate_only_not_phenotype_validation"
    fields = sorted({key for row in rows for key in row})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    counts = {}
    for row in rows:
        counts[row["phylo_decision"]] = counts.get(row["phylo_decision"], 0) + 1
    manifest = {"status": "completed_candidate_only", "decision_counts": counts, "input": str(summary.resolve()), "output": str(output.resolve()), "phenotype_boundary": "Tree decisions are candidate-only and do not prove PHB/PHA degradation.", "rules": {"strong_support": ">=90", "ephaz_distance": "<=1.5 branch-length units", "phac_panel": "always unresolved because PhaC-like challenge sequences are not depolymerase positives"}}
    output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(derive(args.summary, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
