#!/usr/bin/env python3
"""Freeze candidate grades and evidence-coverage statistics for the PhaDED amendment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping


LAYER_FIELDS = (
    "profile_evidence_status", "domain_evidence_status", "motif_evidence_status",
    "localization_evidence_status", "structure_evidence_status_v2",
    "phylogeny_evidence_status_v2", "interpro_status", "subtype_call",
    "subtype_confidence", "evidence_grade",
)
PHENOTYPE_BOUNDARY = "All profile/domain/motif/localization/structure/phylogeny/InterPro labels are candidate-only and do not prove PHB/PHA degradation phenotype."


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _counts(rows: Iterable[Mapping[str, str]], field: str) -> dict[str, int]:
    return dict(sorted(Counter(row.get(field, "") for row in rows).items()))


def summarize_rows(rows: list[Mapping[str, str]]) -> dict[str, object]:
    accessions = [str(row.get("accession", "")) for row in rows]
    if any(not accession for accession in accessions) or len(set(accessions)) != len(accessions):
        raise ValueError("matrix must contain unique, non-empty accessions")
    return {
        "candidate_count": len(rows),
        "subtype_call_counts": _counts(rows, "subtype_call"),
        "subtype_confidence_counts": _counts(rows, "subtype_confidence"),
        "evidence_grade_counts": _counts(rows, "evidence_grade"),
        "evidence_layers": {field: _counts(rows, field) for field in LAYER_FIELDS},
        "architecture_conflict_count": sum(row.get("domain_evidence_status") == "domain_conflict" for row in rows),
        "superfamily_ambiguity_count": sum(row.get("profile_evidence_status") == "profile_ambiguous_superfamily" for row in rows),
        "family_ambiguity_count": sum(row.get("profile_evidence_status") == "profile_ambiguous_family" for row in rows),
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
    }


def _write_tsv(path: Path, rows: list[Mapping[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def _record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}


def write_completion(
    matrix_path: Path,
    motif_manifest: Path,
    conflict_manifest: Path,
    ambiguity_manifest: Path,
    profile_gap_path: Path,
    manual_path: Path,
    output_dir: Path,
    run_id: str,
) -> dict[str, object]:
    matrix_rows = read_tsv(matrix_path)
    summary = summarize_rows(matrix_rows)
    profile_rows = read_tsv(profile_gap_path)
    manual_rows = read_tsv(manual_path)
    profile_coverage = dict(sorted(Counter(row.get("coverage_status", "") for row in profile_rows).items()))
    manual_decisions = dict(sorted(Counter(row.get("manual_review_decision", "") for row in manual_rows).items()))
    output_dir.mkdir(parents=True, exist_ok=True)
    coverage_rows: list[dict[str, object]] = []
    for layer, counts in summary["evidence_layers"].items():
        for status, count in counts.items():
            coverage_rows.append({"evidence_layer": layer, "status": status, "count": count, "interpretation": "candidate-only evidence state"})
    _write_tsv(output_dir / "phaded_evidence_coverage.tsv", coverage_rows, ["evidence_layer", "status", "count", "interpretation"])
    grade_rows = []
    for field in ("subtype_call", "subtype_confidence", "evidence_grade"):
        for value, count in summary["evidence_layers" if field not in {"subtype_call", "subtype_confidence", "evidence_grade"} else field + "_counts"].items():
            grade_rows.append({"metric": field, "value": value, "count": count, "interpretation": "candidate-only; no phenotype validation"})
    _write_tsv(output_dir / "phaded_candidate_grade_summary.tsv", grade_rows, ["metric", "value", "count", "interpretation"])
    metrics = [
        {"metric": "candidate_count", "value": summary["candidate_count"], "interpretation": "accession-level candidate records"},
        {"metric": "architecture_conflict_count", "value": summary["architecture_conflict_count"], "interpretation": "hold; no automatic release"},
        {"metric": "superfamily_ambiguity_count", "value": summary["superfamily_ambiguity_count"], "interpretation": "unresolved; not a biological negative"},
        {"metric": "family_ambiguity_count", "value": summary["family_ambiguity_count"], "interpretation": "stratified sampling only; no threshold lowering"},
        {"metric": "reference_only_profile_count", "value": profile_coverage.get("reference_only_family", 0) + profile_coverage.get("reference_only_superfamily", 0), "interpretation": "planned_not_run pending independent positives, negatives, challenge set, and held-out calibration"},
        {"metric": "manual_retain_count", "value": sum(value for key, value in manual_decisions.items() if key.startswith("retain_")), "interpretation": "candidate-only retention"},
        {"metric": "manual_hold_count", "value": sum(value for key, value in manual_decisions.items() if key.startswith("hold_")), "interpretation": "hold from promotion"},
    ]
    _write_tsv(output_dir / "phaded_evidence_completion_summary.tsv", metrics, ["metric", "value", "interpretation"])
    report = [
        f"# PhaDED evidence completion amendment: {run_id}",
        "",
        PHENOTYPE_BOUNDARY,
        "",
        f"- Accession records: **{summary['candidate_count']}**",
        f"- Subtype calls: `{json.dumps(summary['subtype_call_counts'], ensure_ascii=False, sort_keys=True)}`",
        f"- Motif layer: `{json.dumps(summary['evidence_layers']['motif_evidence_status'], ensure_ascii=False, sort_keys=True)}`",
        f"- Architecture conflicts: **{summary['architecture_conflict_count']}**, superfamily ambiguity: **{summary['superfamily_ambiguity_count']}**, family ambiguity: **{summary['family_ambiguity_count']}**",
        f"- Reference-only profiles: **{profile_coverage.get('reference_only_family', 0) + profile_coverage.get('reference_only_superfamily', 0)}**, all remain `planned_not_run`.",
        f"- Manual 19 decisions: `{json.dumps(manual_decisions, ensure_ascii=False, sort_keys=True)}`",
        "",
        "## Evidence boundaries",
        "",
        "1. Motif pattern matches are sequence-pattern evidence only. `GxSxG` is not a complete catalytic proof; His/Asp/oxyanion-hole and SBD/linker/lid remain separately represented.",
        "2. Architecture conflicts remain hold even when InterPro, motif, structure, or phylogeny layers appear supportive.",
        "3. Ambiguous family/superfamily records remain unresolved candidate-only; the stratified sample is a review queue, not a relabeling operation.",
        "4. Reference-only profiles require independent positive references, family-resolved negatives, challenge sequences, and held-out calibration before training or full-library reassignment.",
    ]
    (output_dir / "phaded_evidence_completion_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    inputs = {name: _record(path) for name, path in {
        "matrix": matrix_path, "motif_manifest": motif_manifest, "conflict_manifest": conflict_manifest,
        "ambiguity_manifest": ambiguity_manifest, "profile_gap": profile_gap_path, "manual_review": manual_path,
    }.items()}
    outputs = {name: _record(output_dir / name) for name in (
        "phaded_evidence_coverage.tsv", "phaded_candidate_grade_summary.tsv",
        "phaded_evidence_completion_summary.tsv", "phaded_evidence_completion_report.md",
    )}
    manifest = {
        "schema_version": "1.0", "run_id": run_id, "status": "completed_candidate_only",
        "summary": summary, "profile_coverage": profile_coverage, "manual_decisions": manual_decisions,
        "inputs": inputs, "outputs": outputs,
        "profile_gap_boundary": "reference-only profiles are planned_not_run until independent positives, family-resolved negatives, challenge set, and held-out calibration exist",
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
    }
    results_manifest_path = output_dir / "phaded_evidence_completion_manifest.json"
    results_manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    run_root = output_dir.parent
    input_contract = run_root / "input_contract.json"
    run_manifest = {
        "schema_version": "1.0",
        "run_id": run_id,
        "status": "completed_candidate_only",
        "input_contract": _record(input_contract) if input_contract.exists() else {"status": "pending", "path": str(input_contract.resolve())},
        "results_manifest": _record(results_manifest_path),
        "source_script": _record(Path(__file__).resolve()),
        "execution": {"mode": "local_read_only", "server": "not_used", "cpu_threads": 0, "gpu": "not_used"},
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
    }
    (run_root / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--motif-manifest", type=Path, required=True)
    parser.add_argument("--conflict-manifest", type=Path, required=True)
    parser.add_argument("--ambiguity-manifest", type=Path, required=True)
    parser.add_argument("--profile-gap", type=Path, required=True)
    parser.add_argument("--manual", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args(argv)
    print(json.dumps(write_completion(args.matrix, args.motif_manifest, args.conflict_manifest, args.ambiguity_manifest, args.profile_gap, args.manual, args.output_dir, args.run_id), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
