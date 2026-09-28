#!/usr/bin/env python3
"""Create a deterministic, stratified review sample for family ambiguity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
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


def gap_bin(value: str) -> str:
    try:
        gap = float(value)
    except (TypeError, ValueError):
        return "missing"
    if gap < 1:
        return "0-1"
    if gap < 5:
        return "1-5"
    if gap < 20:
        return "5-20"
    return "20+"


def stratum_key(row: Mapping[str, str]) -> str:
    return "|".join((
        row.get("phaded_superfamily_best", "") or "unresolved_superfamily",
        gap_bin(row.get("family_score_gap", "")),
        row.get("domain_evidence_status", "") or "domain_missing",
        row.get("motif_evidence_status", "") or "motif_missing",
    ))


def sample_rows(
    rows: Iterable[Mapping[str, str]], *, fraction: float = 0.01, max_per_stratum: int = 25,
) -> list[dict[str, str]]:
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    strata: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("assignment_status") != "ambiguous_family":
            continue
        strata[stratum_key(row)].append(row)
    output: list[dict[str, str]] = []
    for key, members in sorted(strata.items()):
        target = min(max_per_stratum, max(1, math.ceil(len(members) * fraction)))
        ranked = sorted(members, key=lambda row: hashlib.sha256(row.get("accession", "").encode("utf-8")).hexdigest())
        for rank, source in enumerate(ranked[:target], start=1):
            copied = dict(source)
            copied.update({
                "review_stratum": key,
                "stratum_size": str(len(members)),
                "sample_rank": str(rank),
                "sample_target": str(target),
                "sampling_method": "sha256(accession) ascending within deterministic evidence stratum",
                "review_interpretation": "candidate_only_unresolved",
                "recommended_action": "inspect accession-bound profile/domain/motif conflict; do not lower threshold or promote family",
                "phenotype_boundary": "candidate-only sequence evidence; ambiguity is not a biological negative",
            })
            output.append(copied)
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


def run_audit(matrix_path: Path, output_dir: Path, run_id: str, fraction: float = 0.01, max_per_stratum: int = 25) -> dict[str, object]:
    rows = read_tsv(matrix_path)
    ambiguous = [row for row in rows if row.get("assignment_status") == "ambiguous_family"]
    sample = sample_rows(ambiguous, fraction=fraction, max_per_stratum=max_per_stratum)
    stratum_counts = []
    grouped: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in ambiguous:
        grouped[stratum_key(row)].append(row)
    sampled_keys = Counter(row["review_stratum"] for row in sample)
    for key, members in sorted(grouped.items()):
        stratum_counts.append({
            "review_stratum": key,
            "stratum_size": len(members),
            "sample_size": sampled_keys.get(key, 0),
            "sample_fraction": f"{sampled_keys.get(key, 0) / len(members):.6f}",
            "interpretation": "candidate_only_unresolved; no threshold change",
        })
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "family_ambiguity_stratified_sample.tsv", sample)
    _write_tsv(output_dir / "family_ambiguity_stratum_counts.tsv", stratum_counts)
    report = {
        "schema_version": "1.0",
        "run_id": run_id,
        "status": "completed_read_only_stratified_review",
        "candidate_matrix_rows": len(rows),
        "family_ambiguity_rows": len(ambiguous),
        "stratum_count": len(stratum_counts),
        "sample_count": len(sample),
        "sampling_fraction": fraction,
        "max_per_stratum": max_per_stratum,
        "sampling_method": "deterministic SHA-256 accession ranking within superfamily/gap/domain/motif strata",
        "assignment_boundary": "ambiguous_family remains unresolved candidate-only; no threshold lowering or family promotion",
        "input": {"path": str(matrix_path.resolve()), "size": matrix_path.stat().st_size, "sha256": sha256(matrix_path)},
        "outputs": {
            name: {"path": str((output_dir / name).resolve()), "size": (output_dir / name).stat().st_size, "sha256": sha256(output_dir / name)}
            for name in ("family_ambiguity_stratified_sample.tsv", "family_ambiguity_stratum_counts.tsv")
        },
        "phenotype_boundary": "Profile/domain/motif review is candidate evidence and does not prove PHB/PHA degradation.",
    }
    (output_dir / "family_ambiguity_sampling_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--fraction", type=float, default=0.01)
    parser.add_argument("--max-per-stratum", type=int, default=25)
    args = parser.parse_args(argv)
    print(json.dumps(run_audit(args.matrix, args.output_dir, args.run_id, args.fraction, args.max_per_stratum), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
