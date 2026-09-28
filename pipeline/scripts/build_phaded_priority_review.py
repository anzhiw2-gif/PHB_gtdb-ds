#!/usr/bin/env python3
"""Build ordered PhaDED review and bounded InterProScan candidate tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable


STAGE_ORDER = {
    "conflicting": 0,
    "ePhaZ_competition_review": 1,
    "ePhaZ_curated_core_secreted": 2,
    "ePhaZ_curated_core_nonsecreted": 2,
    "iPhaZ_tier1": 3,
}
ALLOWED_LAYERS = set(STAGE_ORDER) - {"conflicting"}


def _score(row: dict[str, str]) -> float:
    try:
        return float(row.get("priority_score", 0.0))
    except (TypeError, ValueError):
        return 0.0


def _stage(row: dict[str, str]) -> str:
    if row.get("architecture_consistency") == "conflicting":
        return "conflicting"
    return row.get("layer", "")


def build_review_queue(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    """Keep requested strata, deduplicate accessions, and apply review order."""
    selected: dict[str, dict[str, str]] = {}
    for source in rows:
        accession = source.get("accession", "").strip()
        stage = _stage(source)
        if not accession or (stage != "conflicting" and stage not in ALLOWED_LAYERS):
            continue
        candidate = dict(source)
        candidate["review_stage"] = stage
        current = selected.get(accession)
        if current is None or (STAGE_ORDER[stage], -_score(candidate), accession) < (
            STAGE_ORDER[current["review_stage"]], -_score(current), accession
        ):
            selected[accession] = candidate
    queued = sorted(
        selected.values(),
        key=lambda row: (STAGE_ORDER[row["review_stage"]], -_score(row), row["accession"]),
    )
    for index, row in enumerate(queued, start=1):
        row["review_order"] = str(index)
    return queued


def select_interpro_candidates(rows: Iterable[dict[str, str]], cap_per_stratum: int = 10) -> list[dict[str, str]]:
    """Select all conflicts/competition and capped family-aware representatives."""
    if cap_per_stratum < 1:
        raise ValueError("cap_per_stratum must be positive")
    queue = build_review_queue(rows)
    output: list[dict[str, str]] = []
    counts: dict[tuple[str, str, str], int] = {}
    for row in queue:
        stage = row["review_stage"]
        if stage == "conflicting":
            reason = "all_conflicting_records"
        elif stage == "ePhaZ_competition_review":
            reason = "all_ePhaZ_competition_records"
        else:
            key = (stage, row.get("phaded_superfamily_best", ""), row.get("phaded_family_best", ""))
            if counts.get(key, 0) >= cap_per_stratum:
                continue
            counts[key] = counts.get(key, 0) + 1
            reason = "top_priority_per_layer_superfamily_family"
        selected = dict(row)
        selected["interpro_selection_reason"] = reason
        output.append(selected)
    for index, row in enumerate(output, start=1):
        row["interpro_order"] = str(index)
    return output


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    fields = sorted({key for row in rows for key in row}) if rows else ["accession"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cap-per-stratum", type=int, default=10)
    args = parser.parse_args(argv)
    with args.input.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    queue = build_review_queue(rows)
    interpro = select_interpro_candidates(queue, args.cap_per_stratum)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    queue_path = args.output_dir / "priority_review_queue.tsv"
    interpro_path = args.output_dir / "interpro_priority_candidates.tsv"
    _write(queue_path, queue)
    _write(interpro_path, interpro)
    summary = {
        "status": "candidate-only",
        "input": {"path": str(args.input.resolve()), "size": args.input.stat().st_size, "sha256": _sha256(args.input)},
        "review_queue": {"path": str(queue_path.resolve()), "count": len(queue), "sha256": _sha256(queue_path)},
        "interpro_candidates": {"path": str(interpro_path.resolve()), "count": len(interpro), "sha256": _sha256(interpro_path)},
        "selection_rule": {
            "review_order": ["conflicting", "ePhaZ_competition_review", "ePhaZ_curated_core", "iPhaZ_tier1"],
            "conflicting": "all records with architecture_consistency=conflicting",
            "competition": "all ePhaZ_competition_review records",
            "capped_strata": "top priority_score per layer + PhaDED superfamily + family",
            "cap_per_stratum": args.cap_per_stratum,
        },
        "phenotype_boundary": "Review and InterProScan domain evidence are candidate-only and do not prove PHB/PHA degradation.",
    }
    (args.output_dir / "priority_review_manifest.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "candidate-only", "review_count": len(queue), "interpro_count": len(interpro)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
