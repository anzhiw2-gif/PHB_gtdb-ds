#!/usr/bin/env python3
"""Merge baseline motif/SignalP features with independent Pfam architecture evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Mapping


PFAM_COLUMNS = ["pfam_accessions", "pfam_coordinates", "pfam_interpro_state", "architecture_consistency", "assignment_review"]


def merged_fields(base_fields: list[str]) -> list[str]:
    """Return the baseline schema plus Pfam fields not already present."""
    fields = list(base_fields)
    fields.extend(field for field in PFAM_COLUMNS if field not in fields)
    return fields


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge(feature_rows: list[Mapping[str, str]], pfam_rows: list[Mapping[str, str]]) -> list[dict[str, str]]:
    by_accession = {row.get("accession", ""): row for row in pfam_rows}
    if len(by_accession) != len(pfam_rows):
        raise ValueError("Pfam evidence has duplicate or missing accessions")
    output = []
    for feature in feature_rows:
        accession = feature.get("accession", "")
        pfam = by_accession.pop(accession, None)
        if pfam is None:
            raise ValueError(f"baseline feature has no Pfam evidence: {accession}")
        row = dict(feature)
        for field in PFAM_COLUMNS:
            row[field] = str(pfam.get(field, ""))
        output.append(row)
    if by_accession:
        raise ValueError("Pfam evidence has accessions absent from baseline feature table")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--pfam-evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.features.open(encoding="utf-8-sig", newline="") as handle:
        features = list(csv.DictReader(handle, delimiter="\t"))
    with args.pfam_evidence.open(encoding="utf-8-sig", newline="") as handle:
        pfam = list(csv.DictReader(handle, delimiter="\t"))
    rows = merge(features, pfam)
    fields = merged_fields(list(features[0]))
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    manifest = {
        "status": "completed_candidate_only",
        "candidate_count": len(rows),
        "inputs": {"features": {"path": str(args.features.resolve()), "sha256": sha256(args.features)}, "pfam_evidence": {"path": str(args.pfam_evidence.resolve()), "sha256": sha256(args.pfam_evidence)}},
        "output": {"path": str(args.output.resolve()), "sha256": sha256(args.output)},
        "interpretation": "Pfam-derived architecture consistency supplements motif/SignalP fields and is not phenotype validation.",
    }
    (args.output.parent / "phaded_feature_consistency_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"candidate_count": len(rows), "status": "completed_candidate_only"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
