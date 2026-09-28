#!/usr/bin/env python3
"""Finalize the candidate-only four-seed ePhaZ calibration audit package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import date
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return fields, list(reader)


def _record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path), "status": "verified"}


def finalize(run_dir: str | Path) -> dict[str, object]:
    run_dir = Path(run_dir)
    inputs = run_dir / "inputs"
    results = run_dir / "results" / "calibration"
    panel_path = inputs / "subtype_training_panels.tsv"
    metrics_path = results / "leave_one_out_metrics.tsv"
    decision_path = results / "subtype_model_decision.json"
    contract_path = run_dir / "input_contract.json"
    required = [panel_path, metrics_path, decision_path, contract_path]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing calibration inputs: " + ", ".join(missing))

    panel_fields, panel_rows = _read_tsv(panel_path)
    metric_fields, metric_rows = _read_tsv(metrics_path)
    if "subtype" not in panel_fields or "accession" not in panel_fields:
        raise ValueError("panel manifest must include subtype and accession")
    if "subtype" not in metric_fields or "decision" not in metric_fields:
        raise ValueError("metrics must include subtype and decision")
    decision_payload = json.loads(decision_path.read_text(encoding="utf-8"))
    decisions = decision_payload.get("decisions", {})
    if not isinstance(decisions, dict):
        raise ValueError("subtype_model_decision.json decisions must be an object")

    audit_fields = panel_fields + ["subtype_decision", "audit_class", "strict_training_eligible"]
    audit_rows: list[dict[str, str]] = []
    for row in panel_rows:
        subtype = row["subtype"].strip()
        role = row.get("role", "").strip()
        evidence = row.get("evidence_status", "").strip()
        subtype_decision = str(decisions.get(subtype, "missing_decision"))
        if subtype_decision == "missing_decision":
            raise ValueError(f"panel subtype missing decision: {subtype}")
        if role == "challenge_control":
            audit_class = "challenge_control"
        elif role == "formal_negative":
            audit_class = "formal_negative"
        elif evidence == "experimental_positive" and role == "train_positive":
            audit_class = "experimental_train_positive"
        elif evidence == "experimental_positive" and role == "heldout_positive":
            audit_class = "experimental_heldout_positive"
        else:
            audit_class = "unclassified_panel_record"
        eligible = str(
            role == "train_positive"
            and evidence == "experimental_positive"
            and subtype_decision == "calibrated_candidate_model"
        ).lower()
        audit_rows.append({**row, "subtype_decision": subtype_decision, "audit_class": audit_class, "strict_training_eligible": eligible})

    audit_path = results / "ephaz_seed_audit.tsv"
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)

    report_path = results / "calibration_status.md"
    summary_lines = [
        "# 四种子 ePhaZ 校准状态",
        "",
        f"- run_id: `{run_dir.name}`",
        f"- 完成日期: `{date.today().isoformat()}`",
        "- 状态: `completed_candidate_only`",
        "- 正式扫描: 未启动",
        "- `formal_scan_models.tsv`: 未修改",
        "- 历史 `ePhaZ_curated_core`: 未覆盖",
        "- 模型升级授权: `false`",
        "",
        "## 分型结果",
        "",
        "| subtype | train | held-out | genera | held-out recovered | negative hits | decision |",
        "|---|---:|---:|---:|---|---:|---|",
    ]
    for metric in metric_rows:
        summary_lines.append(
            "| {subtype} | {train_positive} | {heldout_positive} | {genera} | {heldout_recovered} | {negative_hits} | `{decision}` |".format(**metric)
        )
    summary_lines.extend([
        "",
        "## 解释",
        "",
        "三个 subtype 均未满足预先声明的严格条件（至少 3 个独立训练阳性、" 
        "有独立 held-out 阳性、至少 3 个属、held-out 全部回收且阴性/挑战集零命中）。",
        "因此当前四种子模型只能作为历史 candidate evidence 或 reference-only 证据，"
        "不能升级为广谱分类器，也不能据此修改正式模型注册表或启动新的 GTDB 全量扫描。",
        "`e_dPHAMCL` 的 `Q84C08` challenge control 被标记为命中，进一步说明该面板不能作为零交叉反应分类器。",
        "所有序列、结构域、HMM 或校准证据仍然是候选同源/功能潜力证据，不是已验证 PHB/PHA 降解表型。",
        "",
        "## 输出",
        "",
        "- `ephaz_seed_audit.tsv`: 逐条 panel 记录及审计分类",
        "- `leave_one_out_metrics.tsv`: 分型指标",
        "- `subtype_model_decision.json`: 不升级模型的机器可读决策",
        "- `calibration_sha256_manifest.json`: 输入和输出 SHA-256",
    ])
    report_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    # Do not include the manifest itself in its own hash inventory.
    files: dict[str, dict[str, object]] = {}
    for path in sorted(inputs.glob("*")):
        if path.is_file():
            files[f"inputs/{path.name}"] = _record(path)
    for path in sorted(results.glob("*")):
        if path.is_file() and path.name != "calibration_sha256_manifest.json":
            files[f"results/calibration/{path.name}"] = _record(path)
    hash_manifest = {
        "run_id": run_dir.name,
        "status": "completed_candidate_only",
        "model_promotion_authorized": False,
        "formal_registry_modified": False,
        "files": files,
    }
    hash_path = results / "calibration_sha256_manifest.json"
    hash_path.write_text(json.dumps(hash_manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")

    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract["status"] = "completed_candidate_only"
    contract.setdefault("authorization", {})
    contract["authorization"].update({
        "formal_scan_authorized": False,
        "formal_registry_modified": False,
        "formal_scan_started": False,
        "server_execution_started": False,
    })
    contract["calibration_outputs"] = {
        "seed_audit": _record(audit_path),
        "metrics": _record(metrics_path),
        "decision": _record(decision_path),
        "status_report": _record(report_path),
        "sha256_manifest": _record(hash_path),
    }
    contract_path.write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {
        "status": "completed_candidate_only",
        "model_promotion_authorized": False,
        "formal_registry_modified": False,
        "outputs": {key: value["path"] for key, value in contract["calibration_outputs"].items()},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(finalize(args.run_dir), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
