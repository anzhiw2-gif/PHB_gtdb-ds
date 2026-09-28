#!/usr/bin/env python3
"""Reconcile PhaDED assignments, manual architecture holds, and literature bridges.

This report is evidence bookkeeping only. Exact literature accessions do not
transfer phenotype claims to candidates unless the candidate accession joins
the bridge at the same genome/protein/gene identifier.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def value(row: dict[str, str], key: str, default: str = "UNSPECIFIED") -> str:
    text = (row.get(key) or "").strip()
    return text if text else default


def norm_accession(text: str) -> str:
    text = (text or "").strip()
    for prefix in ("RS_", "GB_"):
        if text.startswith(prefix):
            text = text[len(prefix):]
    return text


def write_rows(path: str | Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def stratify(assignments: list[dict[str, str]], evidence: list[dict[str, str]], output: str | Path) -> list[dict[str, str]]:
    evidence_by_id = {value(row, "accession", ""): row for row in evidence}
    dimensions = ("phaded_superfamily_best", "phaded_family_best", "assignment_status", "layer", "architecture_consistency", "feature_evidence_status")
    counts: Counter[tuple[str, str]] = Counter()
    total = len(assignments)
    for row in assignments:
        ev = evidence_by_id.get(value(row, "accession", ""), {})
        for dimension in dimensions:
            source = ev if dimension in ("architecture_consistency", "feature_evidence_status") else row
            counts[(dimension, value(source, dimension))] += 1
    rows = []
    for (dimension, category), count in sorted(counts.items()):
        rows.append({"dimension": dimension, "category": category, "count": str(count), "fraction": f"{count / total:.8f}" if total else "0"})
    write_rows(output, rows, ["dimension", "category", "count", "fraction"])
    return rows


def architecture_holds(manual: list[dict[str, str]], output: str | Path) -> list[dict[str, str]]:
    columns = ["candidate_id", "manual_review_decision", "superfamily", "nearest_reference_family", "sequence_length", "sequence_integrity", "structure_quality", "pfam_architecture", "interpro_state", "phylogeny_support", "architecture_disposition", "exclude_from_retained_candidates", "exclusion_reason", "phenotype_boundary"]
    rows = []
    for source in manual:
        decision = value(source, "manual_review_decision", "pending")
        excluded = decision.startswith("hold_")
        if decision == "hold_architecture_conflict":
            reason = "结构域指纹与目标 superfamily 不协调，保留为 hold"
        elif decision == "hold_gene_model_or_structure":
            reason = "基因模型不完整或结构质量不足，保留为 hold"
        else:
            reason = "未发现需排除的人工复核状态"
        rows.append({
            "candidate_id": value(source, "candidate_id", ""),
            "manual_review_decision": decision,
            "superfamily": value(source, "superfamily"),
            "nearest_reference_family": value(source, "nearest_reference_family"),
            "sequence_length": value(source, "sequence_length"),
            "sequence_integrity": value(source, "sequence_integrity"),
            "structure_quality": value(source, "structure_quality"),
            "pfam_architecture": value(source, "pfam_architecture"),
            "interpro_state": value(source, "interpro_state"),
            "phylogeny_support": value(source, "phylogeny_support"),
            "architecture_disposition": "excluded_hold" if excluded else "retained_candidate_only",
            "exclude_from_retained_candidates": "yes" if excluded else "no",
            "exclusion_reason": reason,
            "phenotype_boundary": value(source, "phenotype_boundary", "candidate_only"),
        })
    write_rows(output, rows, columns)
    return rows


def literature_bridge(assignments: list[dict[str, str]], bridge: list[dict[str, str]], output: str | Path) -> list[dict[str, str]]:
    by_genome: defaultdict[str, list[str]] = defaultdict(list)
    for row in assignments:
        genome = norm_accession(value(row, "genome", ""))
        if genome and genome != "UNSPECIFIED":
            by_genome[genome].append(value(row, "accession", ""))
    grouped: dict[tuple[str, str, str], dict[str, object]] = {}
    for row in bridge:
        key = (value(row, "record_id", ""), value(row, "evidence_tier"), value(row, "match_class"))
        item = grouped.setdefault(key, {"source": row, "genomes": set(), "pmids": set()})
        g = norm_accession(value(row, "gtdb_genome_accession", ""))
        if g and g != "UNSPECIFIED":
            item["genomes"].add(g)
        pmid = value(row, "pmid", "")
        if pmid != "UNSPECIFIED":
            item["pmids"].add(pmid)
    rows = []
    columns = ["record_id", "literature_protein_accession", "literature_genome_accession", "evidence_tier", "match_class", "pmid", "bridge_genome_count", "candidate_exact_genome_match_count", "candidate_taxonomy_context_genome_count", "candidate_accession_count", "candidate_accessions", "candidate_bridge_status", "interpretation"]
    for (record_id, tier, match_class), item in sorted(grouped.items()):
        source = item["source"]
        exact = tier == "direct_literature_supported" and match_class.startswith("exact_")
        matched = []
        exact_genomes = []
        for genome in sorted(item["genomes"]):
            candidate_ids = by_genome.get(genome, [])
            matched.extend(candidate_ids)
            if candidate_ids and exact:
                exact_genomes.append(genome)
        matched = sorted(set(matched))
        status = "exact_candidate_match" if matched and exact else ("external_exact_positive_no_candidate_overlap" if exact else "no_exact_candidate_bridge")
        rows.append({
            "record_id": record_id,
            "literature_protein_accession": value(source, "literature_protein_accession", ""),
            "literature_genome_accession": value(source, "literature_genome_accession", ""),
            "evidence_tier": tier,
            "match_class": match_class,
            "pmid": ";".join(sorted(item["pmids"])),
            "bridge_genome_count": str(len(item["genomes"])),
            "candidate_exact_genome_match_count": str(len(exact_genomes)),
            "candidate_taxonomy_context_genome_count": str(len({norm_accession(g) for g in item["genomes"] if norm_accession(g) in by_genome}) if match_class == "taxonomy_only_match" else 0),
            "candidate_accession_count": str(len(matched)),
            "candidate_accessions": ";".join(matched),
            "candidate_bridge_status": status,
            "interpretation": "精确文献支持仅适用于 accession 相同的记录；无候选重叠不构成候选阴性。" if exact else "taxonomy-only 或 unresolved 记录不提供候选的精确文献支持。",
        })
    write_rows(output, rows, columns)
    return rows


def file_binding(path: str | Path) -> dict[str, object]:
    p = Path(path)
    if not p.is_file():
        return {"path": str(p.resolve()), "size": None, "sha256": None, "status": "pending"}
    return {"path": str(p.resolve()), "size": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "status": "verified"}


def run_reconciliation(assignment_path: str | Path, evidence_path: str | Path, manual_path: str | Path, bridge_path: str | Path, stats_path: str | Path, holds_path: str | Path, matches_path: str | Path) -> dict[str, int]:
    assignments = read_tsv(assignment_path)
    evidence = read_tsv(evidence_path)
    manual = read_tsv(manual_path)
    bridge = read_tsv(bridge_path)
    stratify(assignments, evidence, stats_path)
    architecture_holds(manual, holds_path)
    bridge_rows = literature_bridge(assignments, bridge, matches_path)
    return {
        "candidate_count": len(assignments),
        "manual_review_count": len(manual),
        "exact_candidate_genome_count": sum(int(r["candidate_exact_genome_match_count"]) for r in bridge_rows if r["evidence_tier"] == "direct_literature_supported"),
        "exact_bridge_record_count": sum(1 for r in bridge_rows if r["evidence_tier"] == "direct_literature_supported"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignment", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--manual", required=True)
    parser.add_argument("--bridge", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--holds", required=True)
    parser.add_argument("--matches", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()
    summary = run_reconciliation(args.assignment, args.evidence, args.manual, args.bridge, args.stats, args.holds, args.matches)
    manifest = {"schema_version": "1.0", "run_id": Path(args.manifest).parent.parent.name, "generated_at": date.today().isoformat(), "status": "completed", "summary": summary, "phenotype_boundary": "All PhaDED, domain, structure, phylogeny, and literature joins remain computational evidence; only exact accession literature rows are direct literature support.", "inputs": {name: file_binding(path) for name, path in (("input_contract", Path(args.manifest).parent.parent / "input_contract.json"), ("assignment", args.assignment), ("evidence", args.evidence), ("manual_review_19", args.manual), ("gtdb_bridge", args.bridge))}, "outputs": {name: file_binding(path) for name, path in (("stratified_stats", args.stats), ("architecture_holds", args.holds), ("literature_matches", args.matches))}, "source_code": file_binding(Path(__file__)), "environment": {"python_version": sys.version, "dependencies": "Python standard library only"}}
    Path(args.manifest).write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
