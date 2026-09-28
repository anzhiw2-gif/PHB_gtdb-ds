#!/usr/bin/env python3
"""Reconcile candidate, score, and PhaDED assignment tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit(candidates_path: Path, scores_path: Path, assignments_path: Path, output_path: Path) -> dict[str, object]:
    with candidates_path.open(encoding="utf-8-sig", newline="") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    with scores_path.open(encoding="utf-8-sig", newline="") as handle:
        scores = list(csv.DictReader(handle, delimiter="\t"))
    with assignments_path.open(encoding="utf-8-sig", newline="") as handle:
        assignments = list(csv.DictReader(handle, delimiter="\t"))
    candidate_ids = [row.get("accession", "") for row in candidates]
    assignment_ids = [row.get("accession", "") for row in assignments]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError("candidate accession table contains duplicates")
    if len(assignment_ids) != len(set(assignment_ids)):
        raise ValueError("assignment table contains duplicates")
    if candidate_ids != assignment_ids:
        raise ValueError("assignment accession order does not reconcile to candidate table")
    score_keys = [(row.get("accession", ""), row.get("profile_id", "")) for row in scores]
    if len(score_keys) != len(set(score_keys)):
        raise ValueError("score table contains duplicate accession/profile rows")
    candidate_by_accession = {row["accession"]: row for row in candidates}
    status_counts = Counter(row["assignment_status"] for row in assignments)
    by_layer: dict[str, Counter[str]] = defaultdict(Counter)
    for row in assignments:
        by_layer[candidate_by_accession[row["accession"]].get("layer", "")][row["assignment_status"]] += 1
    report = {
        "schema_version": 1,
        "candidate_count": len(candidates),
        "assignment_count": len(assignments),
        "score_row_count": len(scores),
        "candidate_assignment_reconciled": True,
        "status_counts": dict(sorted(status_counts.items())),
        "layer_status_counts": {layer: dict(sorted(counts.items())) for layer, counts in sorted(by_layer.items())},
        "inputs": {
            "candidate_source_map": {"path": str(candidates_path.resolve()), "size": candidates_path.stat().st_size, "sha256": sha256(candidates_path)},
            "score_table": {"path": str(scores_path.resolve()), "size": scores_path.stat().st_size, "sha256": sha256(scores_path)},
            "assignment_table": {"path": str(assignments_path.resolve()), "size": assignments_path.stat().st_size, "sha256": sha256(assignments_path)},
        },
        "phenotype_boundary": "All labels are PhaDED candidate architecture/homology evidence, not validated PHB degradation phenotype.",
    }
    output_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.candidates, args.scores, args.assignments, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
