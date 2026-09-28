#!/usr/bin/env python3
"""Finalize a dated, candidate-only PhaDED subtype reconciliation run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Iterable, Mapping


EVIDENCE_LAYERS = (
    "profile_evidence_status",
    "domain_evidence_status",
    "motif_evidence_status",
    "localization_evidence_status",
    "structure_evidence_status_v2",
    "phylogeny_evidence_status_v2",
    "interpro_status",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _count(rows: Iterable[Mapping[str, str]], field: str) -> dict[str, int]:
    return dict(sorted(Counter(row.get(field, "") for row in rows).items()))


def summarize_rows(rows: list[Mapping[str, str]]) -> dict[str, object]:
    """Return deterministic counts used by the amendment and its QA checks."""
    profile_ambiguity = _count(rows, "profile_evidence_status")
    domain_counts = _count(rows, "domain_evidence_status")
    manual_counts = _count(rows, "manual_review_status")
    return {
        "candidate_count": len(rows),
        "architecture_conflict": domain_counts.get("domain_conflict", 0),
        "superfamily_ambiguity": profile_ambiguity.get("profile_ambiguous_superfamily", 0),
        "family_ambiguity": profile_ambiguity.get("profile_ambiguous_family", 0),
        "manual_focus_total": sum(
            count for state, count in manual_counts.items()
            if state and state != "not_in_19_candidate_review"
        ),
        "manual_hold": sum(
            count for state, count in manual_counts.items() if state.startswith("hold_")
        ),
        "manual_retain": sum(
            count for state, count in manual_counts.items() if state == "retain_candidate_only"
        ),
        "profile_evidence_status_counts": profile_ambiguity,
        "subtype_call_counts": _count(rows, "subtype_call"),
        "subtype_confidence_counts": _count(rows, "subtype_confidence"),
        "interpro_status_counts": _count(rows, "interpro_status"),
        "evidence_layer_counts": {
            status: count
            for field in EVIDENCE_LAYERS
            for status, count in _count(rows, field).items()
        },
        "evidence_layers": {field: _count(rows, field) for field in EVIDENCE_LAYERS},
    }


def _write_tsv(path: Path, rows: list[Mapping[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _input_record(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "size": path.stat().st_size,
        "sha256": sha256(path),
    }


def _base_ledger_check(base_ledger: Path | None, base_manifest: Path | None) -> dict[str, object]:
    if not base_ledger:
        return {"status": "not_provided"}
    rows = read_tsv(base_ledger)
    actual = _count(rows, "evidence_grade")
    result: dict[str, object] = {
        "status": "checked",
        "ledger": _input_record(base_ledger),
        "actual_evidence_grade_counts": actual,
    }
    if base_manifest and base_manifest.exists():
        manifest = json.loads(base_manifest.read_text(encoding="utf-8"))
        expected = manifest.get("grade_counts", {})
        result["manifest"] = _input_record(base_manifest)
        result["manifest_evidence_grade_counts"] = expected
        result["manifest_matches_ledger"] = expected == actual
    return result


def write_reports(
    matrix_path: Path,
    profile_gap_path: Path,
    manual_path: Path,
    output_dir: Path,
    run_id: str,
    base_ledger_path: Path | None = None,
    base_manifest_path: Path | None = None,
    extra_inputs: Mapping[str, Path] | None = None,
) -> dict[str, object]:
    """Write all machine-readable summaries and a candidate-only run manifest."""
    matrix_rows = read_tsv(matrix_path)
    profile_gap_rows = read_tsv(profile_gap_path)
    manual_rows = read_tsv(manual_path)
    matrix_by = {row.get("accession", ""): row for row in matrix_rows if row.get("accession", "")}
    if len(matrix_by) != len(matrix_rows):
        raise ValueError("matrix contains an empty or duplicate accession")
    manual_ids = [row.get("candidate_id", "") for row in manual_rows]
    if any(not value for value in manual_ids) or len(set(manual_ids)) != len(manual_ids):
        raise ValueError("manual review table contains an empty or duplicate candidate_id")

    summary = summarize_rows(matrix_rows)
    output_dir.mkdir(parents=True, exist_ok=True)

    conflict_rows = [
        {"category": "architecture_conflict", "count": summary["architecture_conflict"], "interpretation": "Pfam/architecture conflict; keep hold_architecture_conflict"},
        {"category": "superfamily_ambiguity", "count": summary["superfamily_ambiguity"], "interpretation": "Profile superfamily ambiguity; unresolved is not a negative"},
        {"category": "family_ambiguity", "count": summary["family_ambiguity"], "interpretation": "Family ambiguity; no family call without calibrated profile"},
        {"category": "manual_focus_total", "count": summary["manual_focus_total"], "interpretation": "19 accession-bound manual review records"},
        {"category": "manual_retain_candidate_only", "count": summary["manual_retain"], "interpretation": "Retain as candidate-only"},
        {"category": "manual_hold", "count": summary["manual_hold"], "interpretation": "Hold from promotion"},
        {"category": "profile_unassigned", "count": summary["profile_evidence_status_counts"].get("profile_unassigned", 0) if "profile_evidence_status_counts" in summary else 0, "interpretation": "Reference/profile gap; planned_not_run for new profiles"},
    ]
    _write_tsv(output_dir / "conflict_ambiguity_summary.tsv", conflict_rows, ["category", "count", "interpretation"])

    layer_rows: list[dict[str, object]] = []
    for layer, counts in summary["evidence_layers"].items():
        for status, count in counts.items():
            layer_rows.append({"evidence_layer": layer, "status": status, "count": count})
    _write_tsv(output_dir / "subtype_evidence_layer_counts.tsv", layer_rows, ["evidence_layer", "status", "count"])

    focus_fields = [
        "candidate_id", "review_tier", "manual_review_decision", "superfamily", "nearest_reference_family",
        "subtype_call", "subtype_confidence", "profile_evidence_status", "domain_evidence_status",
        "motif_evidence_status", "motif_completeness", "localization_evidence_status",
        "structure_evidence_status_v2", "phylogeny_evidence_status_v2", "interpro_status",
        "evidence_grade", "literature_evidence_tier", "phenotype_boundary",
    ]
    focus_output: list[dict[str, object]] = []
    for manual in manual_rows:
        accession = manual.get("candidate_id", "")
        row = dict(manual)
        matrix_row = matrix_by.get(accession, {})
        row.update({
            field: matrix_row.get(field, "") or manual.get(field, "")
            for field in focus_fields if field != "candidate_id"
        })
        focus_output.append(row)
    _write_tsv(output_dir / "priority_19_candidate_summary.tsv", focus_output, focus_fields)

    gap_fields = list(profile_gap_rows[0]) if profile_gap_rows else ["profile_id", "coverage_status"]
    gap_output = []
    for row in profile_gap_rows:
        copied = dict(row)
        copied["planned_action"] = (
            "planned_not_run; collect independent positives, family-resolved negatives, and challenge set"
            if row.get("coverage_status", "").startswith("reference_only")
            else "retain_current_calibrated_profile"
        )
        gap_output.append(copied)
    _write_tsv(output_dir / "profile_gap_plan.tsv", gap_output, gap_fields + ["planned_action"])

    base_check = _base_ledger_check(base_ledger_path, base_manifest_path)
    repo_root = Path(__file__).resolve().parents[2]
    source_snapshot = {
        "finalizer": _input_record(Path(__file__).resolve()),
        "interpro_merger": _input_record(repo_root / "pipeline" / "scripts" / "merge_phaded_interpro_evidence.py"),
        "subtype_matrix_builder": _input_record(repo_root / "pipeline" / "scripts" / "build_phaded_subtype_matrix.py"),
        "agents_rules": _input_record(repo_root / "AGENTS.md"),
    }
    report_lines = [
        f"# PhaDED 亚型证据 reconciliation amendment ({run_id})",
        "",
        "本报告为 candidate-only 序列层面结果，不等同于已验证 PHB/PHA 降解表型。",
        "",
        f"- 候选 accession：{summary['candidate_count']}",
        f"- InterPro 状态：{summary['interpro_status_counts']}",
        f"- 架构冲突：{summary['architecture_conflict']}；superfamily ambiguity：{summary['superfamily_ambiguity']}；family ambiguity：{summary['family_ambiguity']}",
        f"- 19 条人工重点候选：{summary['manual_retain']} retain candidate-only，{summary['manual_hold']} hold",
        "- motif 层将 lipase box 与 catalytic Ser/Cys、His、Asp、oxyanion-hole、SBD、linker、lid 分开记录；仅 `GxSxG`/lipase box 不升级为完整催化证据。",
        "- 结构和系统发育只在 accession-bound 审查集合中有状态；全库未测试记录保留 `not_tested_full_library`。",
        "- 37 个 reference-only profile 保留为 `planned_not_run`，不得从父级 superfamily 候选池推导 family call。",
        "",
        "## 历史台账一致性",
        "",
    ]
    if base_check.get("status") == "not_provided":
        report_lines.append("未提供旧台账，未执行历史 grade 一致性核对。")
    elif base_check.get("manifest_matches_ledger") is True:
        report_lines.append("旧台账与其 manifest 的 evidence_grade 计数一致。")
    else:
        report_lines.append("旧台账 TSV 与旧 manifest 的 evidence_grade 计数存在不一致；本 amendment 以旧 TSV 逐行合并，并保留该差异。")
        report_lines.append(f"旧 TSV 实际计数：{base_check.get('actual_evidence_grade_counts', {})}")
        report_lines.append(f"旧 manifest 计数：{base_check.get('manifest_evidence_grade_counts', {})}")
    report_lines.extend([
        "",
        "## 边界",
        "",
        "profile、domain、motif、localization、structure、phylogeny 和 InterPro 证据表示同源性或功能潜力；缺失/未测试不是生物学阴性，也不是实验阳性。",
        "",
    ])
    (output_dir / "reconciliation_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    inputs: dict[str, object] = {
        "matrix": _input_record(matrix_path),
        "profile_gap": _input_record(profile_gap_path),
        "manual_review": _input_record(manual_path),
    }
    if base_ledger_path:
        inputs["base_ledger"] = _input_record(base_ledger_path)
    if base_manifest_path:
        inputs["base_manifest"] = _input_record(base_manifest_path)
    for name, path in (extra_inputs or {}).items():
        inputs[name] = _input_record(path)
    inputs_dir = output_dir.parent / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    source_manifest_path = inputs_dir / "source_manifest.tsv"
    _write_tsv(
        source_manifest_path,
        [
            {"name": name, "path": record["path"], "size": record["size"], "sha256": record["sha256"]}
            for name, record in inputs.items()
        ],
        ["name", "path", "size", "sha256"],
    )
    outputs = {
        name: _input_record(output_dir / name)
        for name in (
            "conflict_ambiguity_summary.tsv", "subtype_evidence_layer_counts.tsv",
            "priority_19_candidate_summary.tsv", "profile_gap_plan.tsv", "reconciliation_report.md",
        )
    }
    manifest = {
        "schema_version": "1.0",
        "run_id": run_id,
        "status": "completed_candidate_only",
        "execution": {"mode": "local_read_only", "server": "not_used", "cpu_threads": 0, "gpu": "not_used"},
        "inputs": inputs,
        "outputs": outputs,
        "summary": summary,
        "source_snapshot": source_snapshot,
        "historical_base_ledger_check": base_check,
        "profile_gap_boundary": "reference_only profiles are planned_not_run until independent positive, family-resolved negative, challenge, and held-out calibration evidence exists",
        "phenotype_boundary": "No computational evidence field proves PHB/PHA degradation phenotype.",
    }
    (output_dir / "reconciliation_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    run_root = output_dir.parent
    contract = {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": date.today().isoformat(),
        "status": "completed_candidate_only",
        "purpose": "Merge validated full-library InterPro accession states and reconcile PhaDED subtype evidence layers",
        "authorization": {
            "candidate_only_execution": True,
            "formal_scan_authorized": False,
            "formal_registry_modified": False,
            "server_execution_started": False,
        },
        "inputs": inputs,
        "source_snapshot": source_snapshot,
        "scope": {
            "candidate_count": summary["candidate_count"],
            "manual_focus_count": summary["manual_focus_total"],
            "reference_only_profile_training": "planned_not_run",
        },
        "phenotype_boundary": "Profile, domain, motif, localization, structure, phylogeny, and InterPro evidence do not prove PHB/PHA degradation phenotype.",
    }
    contract_path = run_root / "input_contract.json"
    contract_path.write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    root_manifest = {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": date.today().isoformat(),
        "status": "completed_candidate_only",
        "input_contract": _input_record(contract_path),
        "results_manifest": _input_record(output_dir / "reconciliation_manifest.json"),
        "inputs_manifest": _input_record(source_manifest_path),
        "results": outputs,
        "source_snapshot": source_snapshot,
        "execution": {"mode": "local_read_only", "server": "not_used", "cpu_threads": 0, "gpu": "not_used"},
        "phenotype_boundary": "No computational evidence field proves PHB/PHA degradation phenotype.",
    }
    (run_root / "run_manifest.json").write_text(json.dumps(root_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--profile-gap", type=Path, required=True)
    parser.add_argument("--manual", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--base-ledger", type=Path)
    parser.add_argument("--base-manifest", type=Path)
    parser.add_argument("--extra-input", action="append", default=[], metavar="NAME=PATH")
    args = parser.parse_args(argv)
    extra: dict[str, Path] = {}
    for item in args.extra_input:
        name, separator, value = item.partition("=")
        if not separator or not name or not value:
            raise SystemExit(f"invalid --extra-input: {item!r}; expected NAME=PATH")
        extra[name] = Path(value)
    print(json.dumps(write_reports(args.matrix, args.profile_gap, args.manual, args.output_dir, args.run_id, args.base_ledger, args.base_manifest, extra), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
