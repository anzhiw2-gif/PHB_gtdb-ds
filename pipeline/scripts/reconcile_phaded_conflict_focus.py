#!/usr/bin/env python3
"""Reconcile architecture conflicts and superfamily ambiguities without promotion."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def select_focus_rows(rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    seen: set[str] = set()
    for source in rows:
        accession = str(source.get("accession", "")).strip()
        if not accession or accession in seen:
            continue
        architecture = source.get("domain_evidence_status") == "domain_conflict"
        superfamily = source.get("profile_evidence_status") == "profile_ambiguous_superfamily"
        if not architecture and not superfamily:
            continue
        row = dict(source)
        if architecture:
            decision = "hold_architecture_conflict"
            action = "accession-bound Pfam/structure/phylogeny review; no automatic release"
            reason = "architecture conflict has priority over profile, InterPro, motif, structure, or phylogeny support"
        else:
            decision = "ambiguous_superfamily_unresolved"
            action = "accession-bound independent profile/domain/structure/phylogeny review"
            reason = "superfamily ambiguity is unresolved competition, not a biological negative"
        row.update({
            "focus_decision": decision,
            "focus_recommended_action": action,
            "focus_reason": reason,
            "focus_phenotype_boundary": "candidate-only; no computational layer validates PHB/PHA degradation",
        })
        output.append(row)
        seen.add(accession)
    return output


def _write_tsv(path: Path, rows: list[Mapping[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["status"], delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def run_audit(matrix_path: Path, output_dir: Path, run_id: str) -> dict[str, object]:
    rows = select_focus_rows(read_tsv(matrix_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "conflict_superfamily_focus.tsv", rows)
    counts = dict(sorted(Counter(row["focus_decision"] for row in rows).items()))
    motif_counts = dict(sorted(Counter(row.get("motif_evidence_status", "") for row in rows).items()))
    report = {
        "schema_version": "1.0",
        "run_id": run_id,
        "status": "completed_accession_bound_candidate_only",
        "focus_row_count": len(rows),
        "focus_decision_counts": counts,
        "motif_status_counts": motif_counts,
        "architecture_conflict_count": counts.get("hold_architecture_conflict", 0),
        "superfamily_ambiguity_count": counts.get("ambiguous_superfamily_unresolved", 0),
        "priority_rule": "architecture conflict overrides every other evidence layer; superfamily ambiguity remains unresolved",
        "input": {"path": str(matrix_path.resolve()), "size": matrix_path.stat().st_size, "sha256": sha256(matrix_path)},
        "output": {"path": str((output_dir / "conflict_superfamily_focus.tsv").resolve()), "size": (output_dir / "conflict_superfamily_focus.tsv").stat().st_size, "sha256": sha256(output_dir / "conflict_superfamily_focus.tsv")},
        "phenotype_boundary": "PhaDED profile, Pfam, InterPro, motif, structure, and phylogeny evidence remain candidate-only.",
    }
    (output_dir / "conflict_superfamily_focus_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args(argv)
    print(json.dumps(run_audit(args.matrix, args.output_dir, args.run_id), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
