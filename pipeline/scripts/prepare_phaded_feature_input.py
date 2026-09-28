#!/usr/bin/env python3
"""Bind layer SignalP labels and PhaDED mapping labels for feature extraction."""

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


def prepare(source_map: Path, assignments: Path, output: Path) -> dict[str, object]:
    with source_map.open(encoding="utf-8-sig", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    with assignments.open(encoding="utf-8-sig", newline="") as handle:
        assignment_rows = list(csv.DictReader(handle, delimiter="\t"))
    by_accession = {row["accession"]: row for row in assignment_rows}
    if len(by_accession) != len(assignment_rows) or len(source_rows) != len(assignment_rows):
        raise ValueError("source and assignment tables must have identical unique accessions")
    fields = list(source_rows[0]) + ["phaded_superfamily", "signalp_class", "signalp_evidence"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for source in source_rows:
            accession = source["accession"]
            assignment = by_accession.get(accession)
            if assignment is None:
                raise ValueError(f"assignment missing for {accession}")
            row = dict(source)
            row["phaded_superfamily"] = assignment.get("phaded_superfamily_best", "") or "pending_unassigned"
            row["signalp_class"] = source.get("signal_type", "") or "pending"
            row["signalp_evidence"] = "accepted_layer_signalp" if source.get("layer", "").startswith("ePhaZ_curated_core") else "not_bound_in_layer_contract"
            writer.writerow(row)
    report = {
        "candidate_count": len(source_rows),
        "inputs": {
            "source_map": {"path": str(source_map.resolve()), "size": source_map.stat().st_size, "sha256": sha256(source_map)},
            "assignments": {"path": str(assignments.resolve()), "size": assignments.stat().st_size, "sha256": sha256(assignments)},
        },
        "output": {"path": str(output.resolve()), "size": output.stat().st_size, "sha256": sha256(output)},
        "status": "candidate-only",
    }
    (output.parent / "feature_input_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-map", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source_map, args.assignments, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
