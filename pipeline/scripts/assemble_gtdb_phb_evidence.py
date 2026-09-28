"""Assemble genome-level PHB evidence from frozen candidate and literature tables."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


OUTPUT_FIELDS = [
    "genome",
    "evidence_tier",
    "evidence_label",
    "literature_status",
    "protein_count",
    "candidate_count",
    "superfamily_count",
    "family_count",
    "assignment_statuses",
    "architecture_states",
    "hold_count",
    "source_runs",
    "interpretation_boundary",
]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _index_features(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = row.get("accession", "")
        if accession:
            indexed[accession] = row
    return indexed


def _index_literature(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    indexed: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        genome = row.get("genome_accession", "")
        if genome:
            indexed[genome].append(row)
    return indexed


def _index_holds(rows: list[dict[str, str]]) -> set[str]:
    hold_decisions = {
        "hold_architecture_conflict",
        "hold_gene_model_or_structure",
    }
    return {
        row.get("candidate_id", "")
        for row in rows
        if row.get("manual_review_decision", "") in hold_decisions
    }


def _tier_for_genome(has_exact_positive: bool, has_hold: bool, has_candidate: bool) -> tuple[str, str]:
    if has_exact_positive:
        return "A", "direct_literature_supported"
    if has_hold:
        return "H", "hold_unresolved"
    if has_candidate:
        return "B", "genome_supported_candidate"
    return "C", "sequence_only_candidate"


def assemble_genome_evidence(
    assignment_path: Path,
    feature_path: Path,
    literature_path: Path,
    manual_review_path: Path,
    output_path: Path,
) -> dict[str, dict[str, str]]:
    assignments = _read_tsv(Path(assignment_path))
    features = _index_features(_read_tsv(Path(feature_path)))
    literature = _index_literature(_read_tsv(Path(literature_path)))
    holds = _index_holds(_read_tsv(Path(manual_review_path)))

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in assignments:
        genome = row.get("genome", "")
        if genome:
            grouped[genome].append(row)

    output: dict[str, dict[str, str]] = {}
    for genome in sorted(grouped):
        proteins = grouped[genome]
        literature_rows = literature.get(genome, [])
        exact_positive = any(
            row.get("evidence_status") == "experimental_positive"
            and row.get("exact_identifier_scope") in {"exact_genome", "exact_protein_genome"}
            for row in literature_rows
        )
        hold_rows = [row for row in proteins if row.get("accession", "") in holds]
        assigned = [
            row
            for row in proteins
            if row.get("assignment_status") in {"assigned", "ambiguous_family", "ambiguous_superfamily", "reference_only_nearest"}
            and row.get("phaded_superfamily_best", "")
        ]
        has_candidate = bool(assigned)
        tier, label = _tier_for_genome(exact_positive, bool(hold_rows), has_candidate)
        feature_rows = [features.get(row.get("accession", ""), {}) for row in proteins]
        literature_status = "direct_positive" if exact_positive else ("reviewed_records" if literature_rows else "not_found")
        output[genome] = {
            "genome": genome,
            "evidence_tier": tier,
            "evidence_label": label,
            "literature_status": literature_status,
            "protein_count": str(len(proteins)),
            "candidate_count": str(len(assigned)),
            "superfamily_count": str(len({row.get("phaded_superfamily_best", "") for row in proteins if row.get("phaded_superfamily_best", "")})),
            "family_count": str(len({row.get("phaded_family_best", "") for row in proteins if row.get("phaded_family_best", "")})),
            "assignment_statuses": ";".join(sorted({row.get("assignment_status", "") for row in proteins if row.get("assignment_status", "")})),
            "architecture_states": ";".join(sorted({row.get("architecture_consistency", "") for row in feature_rows if row.get("architecture_consistency", "")})),
            "hold_count": str(len(hold_rows)),
            "source_runs": "20260909_phaded_candidate_mapping_01;20260913_phaded_19_candidate_manual_review_01",
            "interpretation_boundary": "candidate-only; Tier A requires an exact accession-bound experimental record and does not promote homologs",
        }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output.values())
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignment", type=Path, required=True)
    parser.add_argument("--feature", type=Path, required=True)
    parser.add_argument("--literature", type=Path, required=True)
    parser.add_argument("--manual-review", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assemble_genome_evidence(args.assignment, args.feature, args.literature, args.manual_review, args.output)


if __name__ == "__main__":
    main()
