#!/usr/bin/env python3
"""Rank and summarize PhaDED candidates without making phenotype claims."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


PRIORITY_RULE_VERSION = "candidate_priority_v1"
SUMMARY_FIELDS = ["metric", "value", "unit", "interpretation_boundary"]

def _float(row: dict[str, object], key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def rank_candidates(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    ranked: list[dict[str, object]] = []
    for row in rows:
        status = row.get("assignment_status", "")
        architecture = row.get("architecture_consistency", "pending")
        integrity = row.get("sequence_integrity", "valid")
        signal = row.get("signalp_evidence", "no_record")
        score = 0.0
        score += {"assigned": 4.0, "ambiguous_family": 2.0, "ambiguous_superfamily": 1.0}.get(status, 0.0)
        score += {"consistent": 3.0, "partial": 1.0, "pending": 0.0, "conflicting": -3.0}.get(architecture, 0.0)
        score += {
            "valid": 1.0,
            "terminal_stop_only": 1.0,
            "possible_N_truncation": -1.0,
            "invalid_internal_character": -4.0,
            "invalid_character": -4.0,
        }.get(integrity, 0.0)
        score += {"accepted": 1.0, "no_record": 0.0, "OTHER": 0.0}.get(signal, 0.0)
        score += min(max(_float(row, "superfamily_score_gap"), 0.0), 5.0) / 5.0
        if row.get("reference_novelty", "") == "high":
            score += 0.5
        output = dict(row)
        output["priority_score"] = round(score, 4)
        output["review_bucket"] = "ambiguous" if status.startswith("ambiguous") or status == "unassigned_PhaDED_like" else ("priority" if score >= 5 else "review")
        ranked.append(output)
    return sorted(ranked, key=lambda row: (-float(row["priority_score"]), str(row.get("accession", ""))))


def summarize_candidates(rows: list[dict[str, str]]) -> dict[str, int]:
    proteins = {row.get("accession", "") for row in rows if row.get("accession", "")}
    genome_superfamily = {(row.get("genome", ""), row.get("phaded_superfamily_best", "")) for row in rows if row.get("genome", "") and row.get("phaded_superfamily_best", "")}
    genome_family = {(row.get("genome", ""), row.get("phaded_family_best", "")) for row in rows if row.get("genome", "") and row.get("phaded_family_best", "")}
    return {"protein_count": len(proteins), "genome_superfamily_count": len(genome_superfamily), "genome_family_count": len(genome_family)}


def select_representatives(rows: list[dict[str, str]], max_per_class: int = 10) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    seen_sequences: set[str] = set()
    per_class: defaultdict[str, int] = defaultdict(int)
    ordered = sorted(rows, key=lambda row: (0 if row.get("assignment_status") == "assigned" else 1, row.get("accession", "")))
    for row in ordered:
        digest = row.get("sequence_sha256", "")
        if not digest:
            continue
        phaded_class = row.get("phaded_superfamily_best", "") or "unassigned"
        if digest and digest in seen_sequences:
            continue
        if per_class[phaded_class] >= max_per_class:
            continue
        selected.append(dict(row, selection_reason="validated_unique_sequence"))
        per_class[phaded_class] += 1
        if digest:
            seen_sequences.add(digest)
    return selected


def write_downstream_outputs(
    ranked: list[dict[str, object]],
    summary: dict[str, int],
    representatives: list[dict[str, str]],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    ranked_path = output_dir / "priority_candidate_set.tsv"
    ranked_fields = sorted({key for row in ranked for key in row})
    with ranked_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ranked_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(ranked)
    summary_path = output_dir / "taxonomy_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        units = {
            "protein_count": "proteins",
            "genome_superfamily_count": "genome x PhaDED_superfamily",
            "genome_family_count": "genome x PhaDED_family",
        }
        writer.writerows({
            "metric": metric, "value": value, "unit": units[metric],
            "interpretation_boundary": "candidate-only; not phenotype validation",
        } for metric, value in sorted(summary.items()))
    representative_path = output_dir / "representative_selection.tsv"
    representative_fields = sorted({key for row in representatives for key in row}) if representatives else ["accession"]
    with representative_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=representative_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(representatives)
    boundary = output_dir / "interpretation_boundary.md"
    boundary.write_text(
        "# Candidate-only interpretation boundary\n\n"
        "These tables rank computational candidates and preserve uncertainty. "
        "They do not establish PHB/PHA degradation activity or organism phenotype.\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    with args.input.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ranked = rank_candidates(rows)
    summary = summarize_candidates(rows)
    representatives = select_representatives(rows)
    write_downstream_outputs(ranked, summary, representatives, args.output.parent)
    if args.output != args.output.parent / "priority_candidate_set.tsv":
        args.output.write_text((args.output.parent / "priority_candidate_set.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    print(json.dumps({"status": "candidate-only", "priority_rule_version": PRIORITY_RULE_VERSION, **summary, "representative_count": len(representatives)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
