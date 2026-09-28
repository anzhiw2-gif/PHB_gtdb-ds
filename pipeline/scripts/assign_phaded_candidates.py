#!/usr/bin/env python3
"""Assign existing ePhaZ/iPhaZ candidates to PhaDED profiles, preserving ambiguity."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Mapping


ALLOWED_STATUSES = {
    "assigned", "ambiguous_superfamily", "ambiguous_family", "reference_only_nearest",
    "unassigned_PhaDED_like", "not_tested_missing_sequence", "excluded_input_integrity",
}
OUTPUT_FIELDS = [
    "accession", "genome", "family", "layer", "signal_type", "pha_subtype_evidence",
    "phaded_superfamily_best", "phaded_superfamily_second", "superfamily_score_gap",
    "phaded_family_best", "phaded_family_second", "family_score_gap", "assignment_status",
    "evidence_status", "evidence_model", "evidence_threshold",
]


def _number(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"profile score is not numeric: {value!r}")


def _rank(scores: Mapping[str, object]) -> list[tuple[str, float]]:
    return sorted(((name, _number(score)) for name, score in scores.items()), key=lambda item: (-item[1], item[0]))


def _best_two(scores: Mapping[str, object]) -> tuple[str, float, str, float, float]:
    ranked = _rank(scores)
    if not ranked:
        return "", 0.0, "", 0.0, 0.0
    best, best_score = ranked[0]
    second, second_score = ranked[1] if len(ranked) > 1 else ("", 0.0)
    return best, best_score, second, second_score, best_score - second_score


def assign_candidate(
    accession: str,
    scores: Mapping[str, object],
    profile_status: Mapping[str, str],
    *,
    min_score_gap: float = 1.0,
    profile_metadata: Mapping[str, Mapping[str, str]] | None = None,
    evidence_model: str = "PhaDED_profile_score",
    evidence_threshold: str = "declared_in_run_contract",
) -> dict[str, object]:
    if not accession.strip():
        raise ValueError("candidate accession is required")
    if min_score_gap < 0:
        raise ValueError("min_score_gap must be non-negative")
    unknown = sorted(set(scores) - set(profile_status))
    if unknown:
        raise ValueError("profile scores missing status: " + ",".join(unknown))
    if not scores:
        return {
            "accession": accession, "phaded_superfamily_best": "", "phaded_superfamily_second": "",
            "superfamily_score_gap": 0.0, "phaded_family_best": "", "phaded_family_second": "",
            "family_score_gap": 0.0, "assignment_status": "unassigned_PhaDED_like",
            "evidence_status": "no_profile_score", "evidence_model": evidence_model,
            "evidence_threshold": evidence_threshold,
        }
    metadata = profile_metadata or {}
    superfamily_scores: dict[str, float] = defaultdict(float)
    superfamily_profiles: dict[str, list[str]] = defaultdict(list)
    for profile, score in scores.items():
        value = _number(score)
        details = metadata.get(profile, {})
        superfamily = details.get("superfamily", profile)
        superfamily_scores[superfamily] = max(superfamily_scores[superfamily], value)
        superfamily_profiles[superfamily].append(profile)
    sf_best, sf_score, sf_second, sf_second_score, sf_gap = _best_two(superfamily_scores)
    family_best, family_score, family_second, family_second_score, family_gap = _best_two(scores)
    status = "assigned"
    if not sf_best:
        status = "unassigned_PhaDED_like"
    elif len(superfamily_scores) > 1 and sf_gap < min_score_gap:
        status = "ambiguous_superfamily"
    else:
        same_sf = {profile: scores[profile] for profile in superfamily_profiles[sf_best]}
        family_best, family_score, family_second, family_second_score, family_gap = _best_two(same_sf)
        if len(same_sf) > 1 and family_gap < min_score_gap:
            status = "ambiguous_family"
    if len(superfamily_scores) == 1 and len(scores) > 1 and family_gap < min_score_gap:
        status = "ambiguous_family"
    if len(superfamily_scores) > 1 and sf_gap == 0:
        status = "ambiguous_superfamily"
        sf_best = ""
    if len(scores) > 1 and family_gap == 0 and status == "assigned":
        status = "ambiguous_family"
    if family_best and profile_status.get(family_best) == "reference_only":
        status = "reference_only_nearest"
    return {
        "accession": accession, "phaded_superfamily_best": sf_best,
        "phaded_superfamily_second": sf_second, "superfamily_score_gap": sf_gap,
        "phaded_family_best": family_best, "phaded_family_second": family_second,
        "family_score_gap": family_gap, "assignment_status": status,
        "evidence_status": "profile_score", "evidence_model": evidence_model,
        "evidence_threshold": evidence_threshold,
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def map_candidates(
    candidate_rows: list[dict[str, str]],
    score_rows: list[dict[str, str]],
    profile_status: Mapping[str, str],
    profile_metadata: Mapping[str, Mapping[str, str]],
    *,
    min_score_gap: float = 1.0,
    evidence_threshold: str = "declared_in_run_contract",
) -> list[dict[str, object]]:
    scores_by_accession: dict[str, dict[str, float]] = defaultdict(dict)
    for row in score_rows:
        accession = row.get("accession", "").strip()
        profile = row.get("profile_id", "").strip()
        if not accession or not profile:
            raise ValueError("score table requires accession and profile_id")
        if profile not in profile_status:
            raise ValueError(f"score profile missing from profile manifest: {profile}")
        if profile in scores_by_accession[accession]:
            raise ValueError(f"duplicate profile score: {accession}, {profile}")
        scores_by_accession[accession][profile] = _number(row.get("score"))
    output: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in candidate_rows:
        accession = row.get("accession", "").strip()
        if not accession or accession in seen:
            raise ValueError(f"duplicate or missing candidate accession: {accession}")
        seen.add(accession)
        result = assign_candidate(
            accession, scores_by_accession.get(accession, {}), profile_status,
            min_score_gap=min_score_gap, profile_metadata=profile_metadata,
            evidence_threshold=evidence_threshold,
        )
        result.update({
            "genome": row.get("genome", ""), "family": row.get("family", ""),
            "layer": row.get("layer", row.get("source_layer", "")),
            "signal_type": row.get("signal_type", ""),
            "pha_subtype_evidence": row.get("pha_subtype_evidence", row.get("assignment", "")),
        })
        output.append(result)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--profile-status", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-score-gap", type=float, default=1.0)
    parser.add_argument("--evidence-threshold", default="declared_in_run_contract")
    args = parser.parse_args(argv)
    manifest = _read_tsv(args.profile_status)
    status = {row["profile_id"]: row["model_status"] for row in manifest}
    metadata = {
        row["profile_id"]: {"superfamily": row.get("phaded_superfamily", ""), "family": row.get("phaded_family_id", "")}
        for row in manifest
    }
    rows = map_candidates(
        _read_tsv(args.candidates), _read_tsv(args.scores), status, metadata,
        min_score_gap=args.min_score_gap, evidence_threshold=args.evidence_threshold,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    summary = {"candidate_count": len(rows), "status_counts": {status: sum(row["assignment_status"] == status for row in rows) for status in ALLOWED_STATUSES}}
    (args.output.parent / "mapping_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
