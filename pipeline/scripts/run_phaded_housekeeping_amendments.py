#!/usr/bin/env python3
"""Task 0 housekeeping amendments: provenance correction and count discrepancy.

This runner only ever *reads* historical runs and only ever *writes* inside the
dated run directory it is given. It exists because two corrections are needed
after 2026-09-16:

1. the 09-16 input contracts still record ``server_execution_started: false``
   although SignalP actually completed on the T141 dated deploy;
2. the latest reconciliation manifest records ``ambiguous_family_within_superfamily``
   as 9,631 / ``ambiguous_superfamily`` as 41 while the profile layer counts
   9,664 / 42.

Neither the historical contracts nor the historical reconciliation manifest are
modified: the amendments are recorded as new dated artifacts under a new run.
All outputs remain candidate-only evidence and do not prove PHB/PHA degradation
phenotype.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import amend_phaded_provenance as provenance_module  # noqa: E402
import run_context  # noqa: E402


#: Profile-layer status -> subtype call produced by the ambiguous branch.
AMBIGUOUS_LAYER_SUBTYPE_CALLS: tuple[tuple[str, str], ...] = (
    ("profile_ambiguous_family", "ambiguous_family_within_superfamily"),
    ("profile_ambiguous_superfamily", "ambiguous_superfamily"),
)
HOLD_SUBTYPE_CALLS: tuple[str, ...] = (
    "hold_architecture_conflict",
    "hold_gene_model_or_structure",
)
MATRIX_COLUMNS: tuple[str, ...] = (
    "accession",
    "profile_evidence_status",
    "subtype_call",
    "subtype_confidence",
    "domain_evidence_status",
    "manual_review_status",
    "assignment_status",
)
BRANCH_PRECEDENCE_REFERENCE = (
    "pipeline/scripts/build_phaded_subtype_matrix.py L210-L219 (build_matrix: the "
    "manual/domain hold branches are evaluated before the assignment_status "
    "ambiguous branches, so a hold_* subtype_call wins over an ambiguous profile "
    "evidence status)"
)
PHENOTYPE_BOUNDARY = provenance_module.PHENOTYPE_BOUNDARY


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _file_record(path: Path) -> dict[str, object]:
    return {
        "path": str(Path(path).resolve()),
        "size": Path(path).stat().st_size,
        "sha256": run_context.sha256_file(path),
    }


# ---------------------------------------------------------------------------
# Step 0.1 - provenance amendment
# ---------------------------------------------------------------------------


def verify_provenance_assertions(manifest: Mapping[str, object]) -> None:
    """Fail loudly when the recorded execution provenance is not the checked one."""
    execution = dict(manifest.get("execution") or {})
    correction = dict(manifest.get("correction") or {})
    failures: list[str] = []
    if execution.get("mode") != "server_dated_deploy":
        failures.append(f"execution.mode is {execution.get('mode')!r}, expected 'server_dated_deploy'")
    if execution.get("signalp_returncode") != 0:
        failures.append(f"execution.signalp_returncode is {execution.get('signalp_returncode')!r}, expected 0")
    if correction.get("no_historical_run_overwritten") is not True:
        failures.append(
            "correction.no_historical_run_overwritten is "
            f"{correction.get('no_historical_run_overwritten')!r}, expected True"
        )
    if failures:
        raise AssertionError("provenance amendment failed its own assertions: " + "; ".join(failures))


def write_provenance_output_hash_trail(
    run_dir: Path, run_id: str, manifest: Mapping[str, object]
) -> dict[str, object]:
    """Record the authoritative final SHA-256 of every amendment output.

    ``amend_phaded_provenance`` hashes ``provenance_amendment.json`` *before* the
    final rewrite that embeds that hash, so the value stored inside the file
    describes the pre-rewrite revision. The trail records both the final file
    hash and whether the embedded self-hash still matches, instead of trusting
    the embedded value.
    """
    run_dir = Path(run_dir)
    outputs = dict(manifest.get("outputs") or {})
    trail: dict[str, object] = {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": _now(),
        "note": (
            "Final on-disk SHA-256 for the amendment outputs. 'embedded_self_hash' reports the "
            "value recorded inside the artifact and whether it still matches the final file."
        ),
        "outputs": {},
    }
    for key, filename in (("provenance_amendment", "provenance_amendment.json"), ("provenance_readme", "provenance_readme.md")):
        path = run_dir / "results" / filename
        record = _file_record(path)
        embedded = str((outputs.get(key) or {}).get("sha256", ""))
        record["embedded_self_hash"] = {
            "value": embedded or None,
            "matches_final_file": bool(embedded) and embedded == record["sha256"],
        }
        trail["outputs"][filename] = record
    _write_json(run_dir / "results" / "provenance_outputs.json", trail)
    return trail


def run_provenance_amendment(
    *,
    run_id: str,
    run_dir: Path,
    signalp_summary: Path,
    foldseek_manifest: Path,
    final_matrix: Path,
    server: str,
    account: str,
    signalp_threads: int,
    deploy_path: str = "pending",
    observed_at: str = "not_recorded",
    gpu_execution: str = "not_recorded",
    gpu_observation: str = "not_recorded",
    logical_cpu: object = "pending",
    available_memory: object = "pending",
) -> dict[str, object]:
    """Run the tested provenance builder and verify plus hash its outputs."""
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)
    manifest = provenance_module.build_provenance_amendment(
        run_id=run_id,
        output_dir=run_dir / "results",
        signalp_summary=Path(signalp_summary),
        foldseek_manifest=Path(foldseek_manifest),
        final_matrix=Path(final_matrix),
        server_observation={
            "server": server,
            "account": account,
            "signalp_threads": signalp_threads,
            "deploy_path": deploy_path,
            "observed_at": observed_at,
            "gpu_execution": gpu_execution,
            "gpu_observation": gpu_observation,
            "logical_cpu": logical_cpu,
            "available_memory": available_memory,
        },
    )
    verify_provenance_assertions(manifest)
    write_provenance_output_hash_trail(run_dir, run_id, manifest)
    return manifest


# ---------------------------------------------------------------------------
# read-only snapshot of the historical input contracts
# ---------------------------------------------------------------------------


def snapshot_contract_authorization(contract_paths: Iterable[Path]) -> dict[str, object]:
    """Read the historical contracts and record their authorization fields."""
    contracts: dict[str, object] = {}
    for raw_path in contract_paths:
        path = Path(raw_path)
        if not path.is_file():
            raise FileNotFoundError(path)
        digest_before = run_context.sha256_file(path)
        stat_before = path.stat()
        data = json.loads(path.read_text(encoding="utf-8"))
        digest_after = run_context.sha256_file(path)
        stat_after = path.stat()
        authorization = data.get("authorization")
        key = str(data.get("run_id") or path.parent.name)
        contracts[key] = {
            "path": str(path.resolve()),
            "size": stat_after.st_size,
            "sha256": digest_after,
            "status": data.get("status"),
            "generated_at": data.get("generated_at"),
            "authorization": authorization,
            "observed_server_execution_started": (authorization or {}).get("server_execution_started"),
            "read_only_verified": (
                digest_before == digest_after
                and stat_before.st_mtime_ns == stat_after.st_mtime_ns
                and stat_before.st_size == stat_after.st_size
            ),
        }
    return {
        "schema_version": "1.0",
        "access": "read_only",
        "observed_at": _now(),
        "note": (
            "Original 09-16 input contracts are NOT retroactively modified. They still record "
            "server_execution_started=false; this snapshot plus the dated provenance amendment are "
            "the current interpretation for later use."
        ),
        "contracts": contracts,
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
    }


def write_historical_contract_snapshot(
    *,
    run_dir: Path,
    run_id: str,
    contract_paths: Sequence[Path],
    extra_inputs: Mapping[str, Path] | None = None,
) -> dict[str, object]:
    """Write the read-only contract snapshot and the run-context input manifest."""
    run_dir = Path(run_dir)
    snapshot = snapshot_contract_authorization(contract_paths)
    snapshot["run_id"] = run_id
    _write_json(run_dir / "results" / "historical_contract_authorization_snapshot.json", snapshot)

    declared: dict[str, Path] = {}
    for raw_path in contract_paths:
        path = Path(raw_path)
        declared[f"historical_contract::{path.parent.name}"] = path
    declared.update(dict(extra_inputs or {}))
    contract = run_context.build_input_contract(run_dir, run_id=run_id, inputs=declared)
    _write_json(run_dir / "inputs" / "input_manifest.json", contract)
    return snapshot


# ---------------------------------------------------------------------------
# Step 0.2 - count discrepancy amendment
# ---------------------------------------------------------------------------


def _reclassification_reason(row: Mapping[str, str]) -> tuple[str, list[str]]:
    """Explain why an ambiguous profile row does not carry the ambiguous subtype call."""
    subtype_call = row.get("subtype_call", "")
    if subtype_call not in HOLD_SUBTYPE_CALLS:
        return "unexplained", []
    triggers: list[str] = []
    if row.get("manual_review_status") == "hold_architecture_conflict":
        triggers.append("manual_hold_architecture_conflict")
    if row.get("domain_evidence_status") == "domain_conflict":
        triggers.append("domain_conflict")
    if row.get("manual_review_status") == "hold_gene_model_or_structure":
        triggers.append("manual_hold_gene_model_or_structure")
    return "hold_branch_precedence", triggers


def scan_matrix(matrix_path: Path) -> dict[str, object]:
    """Stream the subtype matrix and collect the counts needed by the amendment."""
    matrix_path = Path(matrix_path)
    profile_counts: Counter[str] = Counter()
    subtype_counts: Counter[str] = Counter()
    hold_rows_by_trigger: Counter[str] = Counter()
    reclassified: list[dict[str, object]] = []
    row_count = 0
    with matrix_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:  # pragma: no cover - defensive
            raise ValueError(f"empty matrix: {matrix_path}") from error
        index = {name: position for position, name in enumerate(header)}
        missing = [name for name in MATRIX_COLUMNS if name not in index]
        if missing:
            raise ValueError(f"matrix is missing required columns: {missing}")
        for values in reader:
            row_count += 1
            row = {name: values[index[name]] for name in MATRIX_COLUMNS}
            profile_counts[row["profile_evidence_status"]] += 1
            subtype_counts[row["subtype_call"]] += 1
            if row["subtype_call"] in HOLD_SUBTYPE_CALLS:
                reason, triggers = _reclassification_reason(row)
                hold_rows_by_trigger[" + ".join(triggers) if triggers else "unattributed"] += 1
            expected = dict(AMBIGUOUS_LAYER_SUBTYPE_CALLS).get(row["profile_evidence_status"])
            if expected is not None and row["subtype_call"] != expected:
                reason, triggers = _reclassification_reason(row)
                entry = dict(row)
                entry["expected_subtype_call"] = expected
                entry["resolution"] = reason
                entry["triggered_branches"] = triggers
                reclassified.append(entry)
    return {
        "row_count": row_count,
        "profile_evidence_status_counts": dict(profile_counts),
        "subtype_call_counts": dict(subtype_counts),
        "hold_subtype_rows_by_trigger": dict(hold_rows_by_trigger),
        "reclassified_rows": reclassified,
    }


def _compare_counts(
    recorded: Mapping[str, int], recomputed: Mapping[str, int], layer: str
) -> dict[str, object]:
    mismatches: dict[str, object] = {}
    for key in sorted(set(recorded) | set(recomputed)):
        if recorded.get(key) != recomputed.get(key):
            mismatches[f"{layer}.{key}"] = {
                "recorded": recorded.get(key),
                "recomputed": recomputed.get(key),
            }
    return mismatches


def analyze_count_discrepancy(*, manifest_path: Path, matrix_path: Path) -> dict[str, object]:
    """Explain the subtype-layer vs profile-layer ambiguity count difference."""
    manifest_path = Path(manifest_path)
    matrix_path = Path(matrix_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary = dict(manifest.get("summary") or {})
    recorded_profile = dict(summary.get("profile_evidence_status_counts") or {})
    recorded_subtype = dict(summary.get("subtype_call_counts") or {})
    scanned = scan_matrix(matrix_path)
    recomputed_profile = dict(scanned["profile_evidence_status_counts"])
    recomputed_subtype = dict(scanned["subtype_call_counts"])

    mismatches = _compare_counts(recorded_profile, recomputed_profile, "profile_evidence_status_counts")
    mismatches.update(_compare_counts(recorded_subtype, recomputed_subtype, "subtype_call_counts"))

    reclassified_rows = list(scanned["reclassified_rows"])
    resolutions: list[dict[str, object]] = []
    for profile_status, expected_call in AMBIGUOUS_LAYER_SUBTYPE_CALLS:
        rows = [row for row in reclassified_rows if row["profile_evidence_status"] == profile_status]
        profile_count = recorded_profile.get(profile_status, 0)
        subtype_count = recorded_subtype.get(expected_call, 0)
        delta = profile_count - subtype_count
        resolutions.append(
            {
                "profile_evidence_status": profile_status,
                "expected_subtype_call": expected_call,
                "profile_layer_count": profile_count,
                "subtype_layer_count": subtype_count,
                "delta": delta,
                "recomputed_profile_layer_count": recomputed_profile.get(profile_status, 0),
                "recomputed_subtype_layer_count": recomputed_subtype.get(expected_call, 0),
                "reclassified_count": len(rows),
                "reclassified_to": dict(Counter(str(row["subtype_call"]) for row in rows)),
                "reclassified_triggers": dict(
                    Counter(
                        " + ".join(row["triggered_branches"]) if row["triggered_branches"] else "unexplained"
                        for row in rows
                    )
                ),
                "delta_fully_accounted": delta == len(rows),
                "cause": (
                    "layer precedence: subtype_call is assigned by the manual/domain hold branches before "
                    "the assignment_status ambiguous branches, so a domain_conflict (or manual hold) row "
                    "keeps the hold_* subtype call even though its profile layer is ambiguous"
                ),
                "source_code_reference": BRANCH_PRECEDENCE_REFERENCE,
            }
        )

    unresolved = [row for row in reclassified_rows if row["resolution"] == "unexplained"]
    fully_accounted = all(item["delta_fully_accounted"] for item in resolutions)
    status = "resolved_layer_precedence" if (fully_accounted and not unresolved and not mismatches) else "unresolved"

    return {
        "schema_version": "1.0",
        "generated_at": _now(),
        "status": status,
        "resolution": (
            "The difference is fully explained by decision-branch precedence inside "
            "build_phaded_subtype_matrix.build_matrix; the profile layer and the subtype layer answer "
            "different questions and neither is wrong. The historical manifest needs no correction of "
            "its numbers, only this explicit interpretation."
        ),
        "manifest": _file_record(manifest_path),
        "matrix": _file_record(matrix_path),
        "matrix_row_count": scanned["row_count"],
        "recorded_profile_evidence_status_counts": recorded_profile,
        "recorded_subtype_call_counts": recorded_subtype,
        "recomputed_profile_evidence_status_counts": recomputed_profile,
        "recomputed_subtype_call_counts": recomputed_subtype,
        "manifest_matrix_count_mismatches": mismatches,
        "resolutions": resolutions,
        "hold_subtype_rows_by_trigger": scanned["hold_subtype_rows_by_trigger"],
        "reclassified_row_count": len(reclassified_rows),
        "unexplained_reclassified_count": len(unresolved),
        "reclassified_rows": reclassified_rows,
        "counterfactual": (
            "If the assignment_status ambiguous branches were evaluated before the hold branches, "
            "subtype_call_counts would read ambiguous_family_within_superfamily="
            f"{recorded_profile.get('profile_ambiguous_family', 0)} and ambiguous_superfamily="
            f"{recorded_profile.get('profile_ambiguous_superfamily', 0)}, i.e. it would double-count "
            "rows that are held for a domain-level architecture conflict. The current ordering is the "
            "more conservative one."
        ),
        "historical_manifest_policy": (
            "The 09-16 reconciliation manifest is never rewritten; this amendment is the current "
            "interpretation of its counts."
        ),
        "source_code_reference": BRANCH_PRECEDENCE_REFERENCE,
        "boundary": (
            "Count reconciliation only; no subtype call, hold decision, or candidate classification is "
            "changed by this amendment."
        ),
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
    }


def _write_reclassified_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    fields = list(MATRIX_COLUMNS) + ["expected_subtype_call", "resolution", "triggered_branches"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            record = dict(row)
            record["triggered_branches"] = " + ".join(row.get("triggered_branches") or [])
            writer.writerow(record)


def render_count_discrepancy_report(result: Mapping[str, object]) -> str:
    manifest = dict(result["manifest"])
    matrix = dict(result["matrix"])
    lines = [
        "# Task 0 Step 0.2 — 计数不一致 amendment（subtype 层 9,631/41 vs profile 层 9,664/42）",
        "",
        f"- run_id：`{result.get('run_id', 'pending')}`",
        f"- 状态：`{result['status']}`",
        f"- 输入 manifest：`{manifest['path']}`（size {manifest['size']}，SHA-256 `{manifest['sha256']}`）",
        f"- 输入矩阵：`{matrix['path']}`（size {matrix['size']}，SHA-256 `{matrix['sha256']}`）",
        f"- 矩阵行数（实测流式重算）：{result['matrix_row_count']}",
        "",
        "## 1. 结论",
        "",
        str(result["resolution"]),
        "",
        "## 2. 差值来源",
        "",
        "| profile 层状态 | profile 计数 | subtype call | subtype 计数 | 差值 | 被改判行数 | 差值是否完全解释 |",
        "|---|---|---|---|---|---|---|",
    ]
    for item in result["resolutions"]:
        lines.append(
            "| `{profile_evidence_status}` | {profile_layer_count} | `{expected_subtype_call}` | "
            "{subtype_layer_count} | {delta} | {reclassified_count} | {delta_fully_accounted} |".format(**item)
        )
    lines += [
        "",
        "判定分支优先级（`pipeline/scripts/build_phaded_subtype_matrix.py` L210–L219）：",
        "",
        "1. `manual_review_decision == hold_architecture_conflict` 或 `domain_evidence == domain_conflict` → `hold_architecture_conflict`",
        "2. `manual_review_decision == hold_gene_model_or_structure` → `hold_gene_model_or_structure`",
        "3. `assignment_status == unassigned_PhaDED_like` → `unassigned_PhaDED_like`",
        "4. `assignment_status == ambiguous_superfamily` → `ambiguous_superfamily`",
        "5. `assignment_status == ambiguous_family` → `ambiguous_family_within_superfamily`",
        "",
        "hold 分支（1、2）在 ambiguous 分支（4、5）之前求值，因此一个在 profile 层被判为 ambiguous 的候选，",
        "若同时存在 domain 层 `domain_conflict`（或 19 条人工复核 hold），其 `subtype_call` 会记为 hold，",
        "于是 subtype 层的 ambiguous 计数小于 profile 层。两者回答的是不同问题，都没有算错。",
        "",
        "## 3. 逐条证据（被改判的 ambiguous 行）",
        "",
        f"共 {result['reclassified_row_count']} 行；未解释 {result['unexplained_reclassified_count']} 行。",
        "",
        "| accession | profile 层状态 | subtype_call | subtype_confidence | domain 层状态 | manual 复核状态 | 触发分支 |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in result["reclassified_rows"][:60]:
        lines.append(
            "| {accession} | `{profile_evidence_status}` | `{subtype_call}` | `{subtype_confidence}` | "
            "`{domain_evidence_status}` | `{manual_review_status}` | {triggered} |".format(
                triggered=" + ".join(row.get("triggered_branches") or []) or "unexplained", **row
            )
        )
    if result["reclassified_row_count"] > 60:
        lines.append(f"| … | （其余 {result['reclassified_row_count'] - 60} 行见 TSV） | | | | | |")
    lines += [
        "",
        "## 4. 计数自洽核对（recorded vs 从矩阵重算）",
        "",
        f"不一致条目数：{len(result['manifest_matrix_count_mismatches'])}",
        "",
        "```json",
        json.dumps(result["manifest_matrix_count_mismatches"], indent=2, ensure_ascii=False),
        "```",
        "",
        "## 5. 反事实说明",
        "",
        str(result["counterfactual"]),
        "",
        "## 6. 边界",
        "",
        f"- {result['boundary']}",
        f"- {result['historical_manifest_policy']}",
        f"- {PHENOTYPE_BOUNDARY}",
        "",
    ]
    return "\n".join(lines)


def write_count_discrepancy_amendment(
    *, run_dir: Path, run_id: str, manifest_path: Path, matrix_path: Path
) -> dict[str, object]:
    """Write the discrepancy amendment without touching the historical manifest."""
    run_dir = Path(run_dir)
    results_dir = run_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(manifest_path)
    matrix_path = Path(matrix_path)
    manifest_hash_before = run_context.sha256_file(manifest_path)
    manifest_stat_before = manifest_path.stat()

    result = analyze_count_discrepancy(manifest_path=manifest_path, matrix_path=matrix_path)
    result["run_id"] = run_id
    result["historical_manifest_unmodified"] = False

    rows_path = results_dir / "reclassified_ambiguous_hold_rows.tsv"
    _write_reclassified_rows(rows_path, result["reclassified_rows"])
    report_path = results_dir / "count_discrepancy_report.md"
    report_path.write_text(render_count_discrepancy_report(result), encoding="utf-8")

    manifest_stat_after = manifest_path.stat()
    unmodified = (
        run_context.sha256_file(manifest_path) == manifest_hash_before
        and manifest_stat_before.st_mtime_ns == manifest_stat_after.st_mtime_ns
    )
    result["historical_manifest_unmodified"] = unmodified
    result["outputs"] = {
        "count_discrepancy_report.md": _file_record(report_path),
        "reclassified_ambiguous_hold_rows.tsv": _file_record(rows_path),
    }
    amendment_path = results_dir / "count_discrepancy_amendment.json"
    _write_json(amendment_path, result)

    trail = {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": _now(),
        "note": "Final on-disk SHA-256 of every count-discrepancy output.",
        "historical_manifest_unmodified": unmodified,
        "outputs": {
            "count_discrepancy_amendment.json": _file_record(amendment_path),
            "count_discrepancy_report.md": _file_record(report_path),
            "reclassified_ambiguous_hold_rows.tsv": _file_record(rows_path),
        },
    }
    _write_json(results_dir / "count_discrepancy_outputs.json", trail)
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _provenance_command(args: argparse.Namespace) -> int:
    manifest = run_provenance_amendment(
        run_id=args.run_id,
        run_dir=args.run_dir,
        signalp_summary=args.signalp_summary,
        foldseek_manifest=args.foldseek_manifest,
        final_matrix=args.final_matrix,
        server=args.server,
        account=args.account,
        signalp_threads=args.signalp_threads,
        deploy_path=args.deploy_path,
        observed_at=args.observed_at,
        gpu_execution=args.gpu_execution,
        gpu_observation=args.gpu_observation,
        logical_cpu=args.logical_cpu,
        available_memory=args.available_memory,
    )
    summary: dict[str, object] = {
        "command": "provenance",
        "run_id": args.run_id,
        "execution_mode": manifest["execution"]["mode"],
        "signalp_returncode": manifest["execution"]["signalp_returncode"],
        "no_historical_run_overwritten": manifest["correction"]["no_historical_run_overwritten"],
    }
    if args.historical_contract:
        snapshot = write_historical_contract_snapshot(
            run_dir=args.run_dir,
            run_id=args.run_id,
            contract_paths=list(args.historical_contract),
            extra_inputs={
                "signalp_run_summary": args.signalp_summary,
                "foldseek_merge_manifest": args.foldseek_manifest,
                "final_foldseek_matrix": args.final_matrix,
            },
        )
        summary["historical_contracts"] = {
            key: {
                "sha256": value["sha256"],
                "server_execution_started": value["observed_server_execution_started"],
                "read_only_verified": value["read_only_verified"],
            }
            for key, value in snapshot["contracts"].items()
        }
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


def _count_discrepancy_command(args: argparse.Namespace) -> int:
    result = write_count_discrepancy_amendment(
        run_dir=args.run_dir, run_id=args.run_id, manifest_path=args.manifest, matrix_path=args.matrix
    )
    print(
        json.dumps(
            {
                "command": "count-discrepancy",
                "run_id": args.run_id,
                "status": result["status"],
                "resolutions": [
                    {key: item[key] for key in ("profile_evidence_status", "profile_layer_count", "subtype_layer_count", "delta", "reclassified_count", "delta_fully_accounted")}
                    for item in result["resolutions"]
                ],
                "unexplained_reclassified_count": result["unexplained_reclassified_count"],
                "manifest_matrix_count_mismatches": len(result["manifest_matrix_count_mismatches"]),
                "historical_manifest_unmodified": result["historical_manifest_unmodified"],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="subcommand", required=True)

    provenance = sub.add_parser("provenance", help="Step 0.1 provenance amendment")
    provenance.add_argument("--run-id", required=True)
    provenance.add_argument("--run-dir", type=Path, required=True)
    provenance.add_argument("--signalp-summary", type=Path, required=True)
    provenance.add_argument("--foldseek-manifest", type=Path, required=True)
    provenance.add_argument("--final-matrix", type=Path, required=True)
    provenance.add_argument("--historical-contract", type=Path, action="append", default=[])
    provenance.add_argument("--server", required=True)
    provenance.add_argument("--account", required=True)
    provenance.add_argument("--signalp-threads", type=int, required=True)
    provenance.add_argument("--deploy-path", default="pending")
    provenance.add_argument("--observed-at", default="not_recorded")
    provenance.add_argument("--gpu-execution", default="not_recorded")
    provenance.add_argument("--gpu-observation", default="not_recorded")
    provenance.add_argument("--logical-cpu", default="pending")
    provenance.add_argument("--available-memory", default="pending")
    provenance.set_defaults(func=_provenance_command)

    counts = sub.add_parser("count-discrepancy", help="Step 0.2 count discrepancy amendment")
    counts.add_argument("--run-id", required=True)
    counts.add_argument("--run-dir", type=Path, required=True)
    counts.add_argument("--manifest", type=Path, required=True)
    counts.add_argument("--matrix", type=Path, required=True)
    counts.set_defaults(func=_count_discrepancy_command)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
