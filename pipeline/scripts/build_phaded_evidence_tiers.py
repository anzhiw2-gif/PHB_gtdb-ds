#!/usr/bin/env python3
"""Build a complete, candidate-only PhaDED evidence ledger with explicit grades."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Mapping


GRADE_FIELD = "evidence_grade"
LITERATURE_FIELD = "literature_evidence_tier"
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_sequence_integrity(declared: str, sequence: str) -> str:
    """Re-evaluate stale integrity labels against the bound FASTA sequence."""
    seq = sequence.strip().upper()
    if seq.endswith("*") and seq.count("*") == 1:
        seq = seq[:-1]
        if all(residue in AMINO_ACIDS for residue in seq):
            return "terminal_stop_only" if seq.startswith("M") else "possible_N_truncation"
    if any(residue not in AMINO_ACIDS for residue in seq):
        return "invalid_internal_character"
    if not seq.startswith("M"):
        return "possible_N_truncation"
    return "valid"


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = None
    chunks: list[str] = []
    for line in path.read_text(encoding="ascii").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                if current in records or not chunks:
                    raise ValueError(f"duplicate or empty FASTA record: {current}")
                records[current] = "".join(chunks)
            current = line[1:].split()[0]
            chunks = []
        else:
            if current is None:
                raise ValueError("sequence appears before FASTA header")
            chunks.append(line)
    if current is not None:
        if current in records or not chunks:
            raise ValueError(f"duplicate or empty FASTA record: {current}")
        records[current] = "".join(chunks)
    if not records:
        raise ValueError("FASTA is empty")
    return records


def grade_row(row: Mapping[str, str]) -> dict[str, str]:
    """Assign a transparent sequence-evidence grade without making a phenotype call."""
    assignment = str(row.get("assignment_status", "")).strip()
    architecture = str(row.get("architecture_consistency", "pending")).strip()
    integrity = str(row.get("sequence_integrity", "pending")).strip()
    if architecture == "conflicting":
        grade = "HOLD_architecture_conflict"
        basis = "PhaDED profile present but Pfam/InterPro architecture is conflicting; retain for review."
    elif integrity in {"invalid_internal_character", "invalid_character", "possible_N_truncation"}:
        grade = "HOLD_sequence_integrity"
        basis = "Sequence contains an integrity warning; this is a quality hold, not a biological negative."
    elif assignment in {"unassigned_PhaDED_like", "", "pending"}:
        grade = "PENDING_profile_unassigned"
        basis = "No assigned PhaDED profile at the current threshold; retain as unresolved candidate."
    elif architecture == "consistent":
        grade = "L2_profile_plus_architecture_consistent"
        basis = "Assigned PhaDED profile with architecture marked consistent."
    elif architecture == "partial":
        grade = "L1_profile_plus_partial_architecture"
        basis = "Assigned PhaDED profile with partial targeted architecture evidence."
    else:
        grade = "L0_profile_only"
        basis = "Assigned PhaDED profile; independent architecture evidence is pending/not tested."
    result = dict(row)
    result[GRADE_FIELD] = grade
    result["evidence_basis"] = basis
    result["structure_evidence_status"] = str(row.get("structure_evidence_status", "not_tested_full_library"))
    result["phylogeny_evidence_status"] = str(row.get("phylogeny_evidence_status", "not_tested_full_library"))
    result["interpro_evidence_status"] = str(row.get("interpro_evidence_status", row.get("pfam_interpro_state", "pending")))
    result[LITERATURE_FIELD] = str(row.get(LITERATURE_FIELD, "none_exact_candidate_match"))
    return result


def literature_map(bridge_path: Path) -> dict[str, str]:
    """Return candidate accession -> literature tier for exact bridge rows only."""
    mapping: dict[str, str] = {}
    with bridge_path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tier = row.get("evidence_tier", "")
            if tier not in {"direct_literature_supported", "experimental_positive"}:
                continue
            for accession in row.get("candidate_accessions", "").split(";"):
                accession = accession.strip()
                if accession:
                    mapping[accession] = "Tier_A_exact_accession_literature"
    return mapping


def build(assignments: Path, features: Path, bridge: Path, output: Path, sequences: Path | None = None) -> dict[str, object]:
    with assignments.open(encoding="utf-8-sig", newline="") as handle:
        assignment_rows = list(csv.DictReader(handle, delimiter="\t"))
    with features.open(encoding="utf-8-sig", newline="") as handle:
        feature_rows = list(csv.DictReader(handle, delimiter="\t"))
    assignment_by_id = {row.get("accession", ""): row for row in assignment_rows}
    feature_by_id = {row.get("accession", ""): row for row in feature_rows}
    if not assignment_by_id or len(assignment_by_id) != len(assignment_rows):
        raise ValueError("assignments contain duplicate or missing accessions")
    if set(assignment_by_id) != set(feature_by_id):
        raise ValueError("assignment and feature accession sets differ")
    exact_literature = literature_map(bridge)
    sequence_by_id = read_fasta(sequences) if sequences is not None else {}
    if sequence_by_id and set(sequence_by_id) != set(assignment_by_id):
        raise ValueError("sequence and assignment accession sets differ")
    rows: list[dict[str, str]] = []
    for accession, assignment in assignment_by_id.items():
        merged = dict(assignment)
        merged.update(feature_by_id[accession])
        merged["assignment_status"] = assignment.get("assignment_status", "")
        merged["architecture_consistency"] = feature_by_id[accession].get("architecture_consistency", "pending")
        merged["assignment_review"] = feature_by_id[accession].get("assignment_review", "no_conflict_recorded")
        if sequence_by_id:
            merged["sequence_integrity"] = normalized_sequence_integrity(
                merged.get("sequence_integrity", "pending"), sequence_by_id[accession]
            )
        merged[LITERATURE_FIELD] = exact_literature.get(accession, "none_exact_candidate_match")
        rows.append(grade_row(merged))
    fields = list(assignment_rows[0])
    for field in feature_rows[0]:
        if field not in fields:
            fields.append(field)
    for field in (GRADE_FIELD, "evidence_basis", "structure_evidence_status", "phylogeny_evidence_status", "interpro_evidence_status", LITERATURE_FIELD):
        if field not in fields:
            fields.append(field)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "status": "completed_candidate_only",
        "candidate_count": len(rows),
        "grade_counts": {},
        "literature_exact_count": sum(row[LITERATURE_FIELD].startswith("Tier_A") for row in rows),
        "phenotype_boundary": "Grades represent sequence/domain/literature evidence and do not prove PHB degradation phenotype.",
        "inputs": {name: {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)} for name, path in (("assignments", assignments), ("features", features), ("literature_bridge", bridge))},
        "output": {"path": str(output.resolve()), "size": output.stat().st_size, "sha256": sha256(output)},
    }
    from collections import Counter
    report["grade_counts"] = dict(Counter(row[GRADE_FIELD] for row in rows))
    (output.parent / "evidence_tier_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--literature-bridge", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sequences", type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.assignments, args.features, args.literature_bridge, args.output, args.sequences), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
