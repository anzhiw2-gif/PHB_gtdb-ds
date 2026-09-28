#!/usr/bin/env python3
"""Merge feature evidence into a separate PhaDED assignment evidence table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge(assignments: Path, features: Path, output: Path) -> dict[str, object]:
    with assignments.open(encoding="utf-8-sig", newline="") as handle:
        assignment_rows = list(csv.DictReader(handle, delimiter="\t"))
    with features.open(encoding="utf-8-sig", newline="") as handle:
        feature_rows = list(csv.DictReader(handle, delimiter="\t"))
    feature_by_accession = {row["accession"]: row for row in feature_rows}
    if len(feature_by_accession) != len(feature_rows):
        raise ValueError("feature table has duplicate accessions")
    if set(feature_by_accession) != {row["accession"] for row in assignment_rows}:
        raise ValueError("feature and assignment accession sets differ")
    fields = list(assignment_rows[0]) + ["architecture_consistency", "feature_evidence_status", "assignment_review"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in assignment_rows:
            feature = feature_by_accession[row["accession"]]
            merged = dict(row)
            merged["architecture_consistency"] = feature.get("architecture_consistency", "pending")
            merged["feature_evidence_status"] = feature.get("feature_evidence_status", "pending")
            merged["assignment_review"] = "review_required" if merged["architecture_consistency"] == "conflicting" else "no_conflict_recorded"
            writer.writerow(merged)
    report = {
        "candidate_count": len(assignment_rows),
        "inputs": {"assignments": {"path": str(assignments.resolve()), "size": assignments.stat().st_size, "sha256": sha256(assignments)}, "features": {"path": str(features.resolve()), "size": features.stat().st_size, "sha256": sha256(features)}},
        "output": {"path": str(output.resolve()), "size": output.stat().st_size, "sha256": sha256(output)},
        "status": "candidate-only",
        "phenotype_boundary": "Feature and assignment evidence remain candidate-only and do not prove PHB degradation.",
    }
    (output.parent / "assignment_evidence_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(merge(args.assignments, args.features, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
