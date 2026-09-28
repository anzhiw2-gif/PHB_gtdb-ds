#!/usr/bin/env python3
"""Run the **authorized** PhaDED calibrated-candidate promotion for a scoped run.

This runner *consumes* ``finalize_phaded_reference_panel_calibration.py`` unmodified.  It
adds exactly one thing the module cannot express: a run whose authorization is **scoped**.
The module has a single global ``promote_calibrated_model`` switch, so "the gate passes but
the promotion decision for this profile is deliberately deferred" has no representation in
its decision table -- and the family gate must never be falsified to fake that outcome.
The runner therefore splits the readiness by its explicit ``promotion_scope`` column and
invokes the module once per scope, then merges the module's own rows:

* ``authorized_promotion_frame`` -> ``finalize(..., promote_calibrated_model=<flag>)``.
  With the flag on, exactly the profiles whose gate returns true become
  ``calibrated_candidate_model``; with the flag off nothing is promoted.
* ``deferred_design_spec_special_family`` -> ``finalize(..., promote_calibrated_model=False)``
  in its own invocation, which is the module's own ``candidate_gate_passed_not_promoted``
  path.  A deferred profile that passes its gate is therefore **visibly** reported as a
  passing-but-unpromoted candidate: the deferral is neither a hidden suppression nor a
  silent promotion.
* the **whole frame** with the flag off is run as a control, proving that no profile can
  reach ``calibrated_candidate_model`` without the explicit authorization flag.

The governing ``leaveout_calibration_decisions.tsv`` and
``leaveout_calibration_gate_detail.tsv`` written into the run's ``results/`` are the union
of the scope invocations in readiness order.  Every merged row is checked row-for-row
against the module output it came from, so the merge cannot rewrite a decision; the
per-scope module outputs stay in the run as ``scopes/<scope>/``.  ``promotion_verification.json``
records every invariant, and the runner refuses to finish if any of them is false.

Boundary: a promotion is a candidate-gate outcome under an explicit authorization.  It is
not an experimental phenotype claim, and it carries no statement about any individual hit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
FINALIZE_SCRIPT = Path(__file__).resolve().parent / "finalize_phaded_reference_panel_calibration.py"

SCOPE_AUTHORIZED = "authorized_promotion_frame"
SCOPE_DEFERRED = "deferred_design_spec_special_family"
SCOPES = (SCOPE_AUTHORIZED, SCOPE_DEFERRED)
SCOPE_FIELD = "promotion_scope"

PROMOTED = "calibrated_candidate_model"
NOT_PROMOTED = "candidate_gate_passed_not_promoted"
NOT_RUN = "not_run"

SCOPES_DIR = "scopes"
CONTROL_DIR = "control_flag_off"
DECISIONS_NAME = "leaveout_calibration_decisions.tsv"
GATE_DETAIL_NAME = "leaveout_calibration_gate_detail.tsv"
REPORT_NAME = "leaveout_calibration_report.json"
STATUS_NAME = "calibration_status.md"

OUTPUT_GATE_TABLE = "per_family_gate_table.tsv"
OUTPUT_SCOPE_MAP = "promotion_scope_map.tsv"
OUTPUT_VERIFICATION = "promotion_verification.json"

GATE_TABLE_FIELDS = (
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "qualified_e2_e3_positive_count", "independent_genus_count",
    "required_qualified_e2_e3_positive_count", "required_independent_genus_count",
    "positive_count", "required_positive_count", "positive_basis",
    "heldout_positive", "negative_count", "challenge_count",
    "family_resolved_negative_status", "challenge_status",
    "gate_kind", "gate_passed", "blocking_slots", "promotion_authorized",
    "promotion_scope", "readiness_row_origin", "final_decision", "calibration_status", "promoted",
)
SCOPE_MAP_FIELDS = (
    "profile_id", "phaded_family_id", "profile_kind", "promotion_scope", "source_invocation",
    "source_decisions_file", "source_gate_detail_file",
)

READINESS_REQUIRED = ("profile_id", "profile_kind", "phaded_family_id")


class PromotionError(RuntimeError):
    """Raised when the promotion cannot be executed or verified; always fails closed."""


def load_finalize_module():
    """Load the frozen finalize script as a module (never modified, only consumed)."""
    spec = importlib.util.spec_from_file_location("finalize_phaded_reference_panel_calibration",
                                                  FINALIZE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _display_path(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO)).replace("\\", "/")
    except ValueError:
        return str(Path(path).resolve()).replace("\\", "/")


def _describe(path: Path) -> dict[str, object]:
    return {"path": _display_path(path), "size": path.stat().st_size, "sha256": _sha256(path)}


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = [{key: (value or "") for key, value in row.items() if key is not None}
                for row in reader]
    return fieldnames, rows


def _write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def promote(
    readiness: str | Path,
    output_dir: str | Path,
    *,
    run_id: str,
    promote_calibrated_model: bool = False,
    defer_profiles: tuple[str, ...] | list[str] = (),
) -> dict[str, object]:
    """Execute the scoped promotion and verify it; any violated invariant raises."""
    readiness, output_dir = Path(readiness), Path(output_dir)
    module = load_finalize_module()
    fieldnames, rows = _read_rows(readiness)
    if SCOPE_FIELD not in fieldnames:
        raise PromotionError(f"readiness has no {SCOPE_FIELD} column; promotion scope unknown")
    missing = [field for field in READINESS_REQUIRED if field not in fieldnames]
    if missing:
        raise PromotionError("readiness missing columns: " + ",".join(missing))
    if not rows:
        raise PromotionError("readiness is empty")

    seen: set[str] = set()
    for row in rows:
        profile_id = row["profile_id"].strip()
        if not profile_id or profile_id in seen:
            raise PromotionError(f"empty or duplicate profile_id: {profile_id!r}")
        seen.add(profile_id)
        scope = row[SCOPE_FIELD].strip()
        if scope not in SCOPES:
            raise PromotionError(f"{profile_id}: unknown {SCOPE_FIELD} {scope!r}")

    in_file_deferred = [row["profile_id"].strip() for row in rows
                        if row[SCOPE_FIELD].strip() == SCOPE_DEFERRED]
    acknowledged = [str(profile_id) for profile_id in defer_profiles]
    if sorted(in_file_deferred) != sorted(acknowledged):
        raise PromotionError(
            "the run must acknowledge exactly the deferred profiles of the readiness: "
            f"readiness={sorted(in_file_deferred)} acknowledged={sorted(acknowledged)}")

    authorized_rows = [row for row in rows if row[SCOPE_FIELD].strip() == SCOPE_AUTHORIZED]
    deferred_rows = [row for row in rows if row[SCOPE_FIELD].strip() == SCOPE_DEFERRED]
    if not authorized_rows:
        raise PromotionError(
            f"the {SCOPE_AUTHORIZED} scope is empty; an authorized promotion run has no frame")

    outputs = (output_dir / DECISIONS_NAME, output_dir / GATE_DETAIL_NAME,
               output_dir / REPORT_NAME, output_dir / STATUS_NAME,
               output_dir / OUTPUT_GATE_TABLE, output_dir / OUTPUT_SCOPE_MAP,
               output_dir / OUTPUT_VERIFICATION)
    for path in outputs:
        if path.exists():
            raise FileExistsError(path)
    for path in (output_dir / SCOPES_DIR, output_dir / CONTROL_DIR):
        if path.exists():
            raise FileExistsError(path)

    scopes_dir = output_dir / SCOPES_DIR
    invocations: dict[str, dict[str, object]] = {}
    for scope, scope_rows, flag in ((SCOPE_AUTHORIZED, authorized_rows, promote_calibrated_model),
                                    (SCOPE_DEFERRED, deferred_rows, False)):
        if not scope_rows:
            continue
        scope_dir = scopes_dir / scope
        scope_dir.mkdir(parents=True, exist_ok=True)
        scope_readiness = scope_dir / "readiness.tsv"
        _write_rows(scope_readiness, fieldnames, scope_rows)
        module.finalize(scope_readiness, scope_dir, run_id=run_id,
                        promote_calibrated_model=flag)
        invocations[scope] = {
            "promotion_authorized": flag,
            "profile_count": len(scope_rows),
            "readiness": _describe(scope_readiness),
            "decisions": _describe(scope_dir / DECISIONS_NAME),
            "gate_detail": _describe(scope_dir / GATE_DETAIL_NAME),
            "report": _describe(scope_dir / REPORT_NAME),
        }

    control_dir = output_dir / CONTROL_DIR
    control_dir.mkdir(parents=True, exist_ok=True)
    control_readiness = control_dir / "readiness.tsv"
    _write_rows(control_readiness, fieldnames, rows)
    module.finalize(control_readiness, control_dir, run_id=run_id,
                    promote_calibrated_model=False)

    # ---- merge the module's own rows, in readiness order --------------------------
    merged: dict[str, list[dict[str, str]]] = {"decisions": [], "gate_detail": []}
    source_rows: dict[str, dict[str, dict[str, str]]] = {}
    for scope in SCOPES:
        if scope not in invocations:
            continue
        scope_dir = scopes_dir / scope
        _, scope_decisions = _read_rows(scope_dir / DECISIONS_NAME)
        _, scope_detail = _read_rows(scope_dir / GATE_DETAIL_NAME)
        source_rows.setdefault("decisions", {}).update(
            {row["profile_id"]: row for row in scope_decisions})
        source_rows.setdefault("gate_detail", {}).update(
            {row["profile_id"]: row for row in scope_detail})
    for row in rows:
        profile_id = row["profile_id"].strip()
        for key in ("decisions", "gate_detail"):
            if profile_id not in source_rows[key]:
                raise PromotionError(f"{profile_id}: no {key} row in any scope invocation")
            merged[key].append(source_rows[key][profile_id])

    _, control_decisions = _read_rows(control_dir / DECISIONS_NAME)
    _, control_detail = _read_rows(control_dir / GATE_DETAIL_NAME)

    # ---- invariants --------------------------------------------------------------
    scope_of = {row["profile_id"].strip(): row[SCOPE_FIELD].strip() for row in rows}
    checks: dict[str, bool] = {}
    checks["merged_row_count_equals_readiness"] = (
        len(merged["decisions"]) == len(rows) and len(merged["gate_detail"]) == len(rows))
    checks["merged_profile_order_equals_readiness"] = (
        [row["profile_id"] for row in merged["decisions"]] == [row["profile_id"].strip() for row in rows]
        and [row["profile_id"] for row in merged["gate_detail"]] == [row["profile_id"].strip() for row in rows])
    checks["merged_decisions_equal_scope_rows"] = all(
        row == source_rows["decisions"][row["profile_id"]] for row in merged["decisions"])
    checks["merged_gate_detail_equal_scope_rows"] = all(
        row == source_rows["gate_detail"][row["profile_id"]] for row in merged["gate_detail"])
    checks["promotion_authorized_true_only_in_authorized_frame"] = all(
        (row["promotion_authorized"] == "true")
        == (promote_calibrated_model and scope_of[row["profile_id"]] == SCOPE_AUTHORIZED)
        for row in merged["gate_detail"])
    checks["promoted_implies_gate_passed_and_authorized"] = all(
        row["gate_passed"] == "true" and row["promotion_authorized"] == "true"
        for row in merged["gate_detail"] if row["decision"] == PROMOTED)
    checks["no_deferred_profile_promoted"] = all(
        row["decision"] != PROMOTED for row in merged["decisions"]
        if scope_of[row["profile_id"]] == SCOPE_DEFERRED)
    checks["flag_off_control_has_no_promotion"] = (
        all(row["decision"] != PROMOTED for row in control_decisions)
        and all(row["calibration_status"] != "passed" for row in control_detail))
    if not promote_calibrated_model:
        checks["no_promotion_without_the_flag"] = all(
            row["decision"] != PROMOTED for row in merged["decisions"])
    failed = sorted(name for name, value in checks.items() if not value)
    if failed:
        raise PromotionError("promotion verification failed: " + ",".join(failed))

    # ---- governing outputs ---------------------------------------------------------
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_rows(output_dir / DECISIONS_NAME, list(module.FIELDS), merged["decisions"])
    _write_rows(output_dir / GATE_DETAIL_NAME, list(module.GATE_FIELDS), merged["gate_detail"])

    decision_counts = Counter(row["decision"] for row in merged["decisions"])
    gate_counts = Counter(row["gate_kind"] for row in merged["gate_detail"])
    promoted = [row["profile_id"] for row in merged["decisions"] if row["decision"] == PROMOTED]
    deferred = [{
        "profile_id": row["profile_id"],
        "gate_passed": source_rows["gate_detail"][row["profile_id"]]["gate_passed"],
        "blocking_slots": source_rows["gate_detail"][row["profile_id"]]["blocking_slots"],
        "decision": row["decision"],
        "calibration_status": row["calibration_status"],
        "promotion_authorized": source_rows["gate_detail"][row["profile_id"]]["promotion_authorized"],
    } for row in merged["decisions"] if scope_of[row["profile_id"]] == SCOPE_DEFERRED]

    summary = {
        "profile_count": len(merged["decisions"]),
        "calibrated_candidate_model": decision_counts.get(PROMOTED, 0),
        "candidate_gate_passed_not_promoted": decision_counts.get(NOT_PROMOTED, 0),
        "reference_only_insufficient_panel": decision_counts.get("reference_only_insufficient_panel", 0),
        "planned_not_run_superfamily_requires_family_resolution": decision_counts.get(
            "planned_not_run_superfamily_requires_family_resolution", 0),
        "new_family_call_created": False,
        "promotion_authorized": promote_calibrated_model,
        "gate_kind_counts": dict(sorted(gate_counts.items())),
        "family_gate_evaluated": sum(row["profile_kind"] == "family" for row in merged["gate_detail"]),
        "superfamily_gate_evaluated": sum(row["profile_kind"] == "superfamily" for row in merged["gate_detail"]),
        "profiles_passing_family_gate": gate_counts.get("family", 0),
        "profiles_passing_superfamily_gate": gate_counts.get("superfamily", 0),
        "authorized_scope_profile_count": len(authorized_rows),
        "deferred_scope_profile_count": len(deferred_rows),
    }

    readiness_lookup = {row["profile_id"].strip(): row for row in rows}
    detail_lookup = {row["profile_id"]: row for row in merged["gate_detail"]}
    decision_lookup = {row["profile_id"]: row for row in merged["decisions"]}
    gate_table = []
    for row in rows:
        profile_id = row["profile_id"].strip()
        detail, decision = detail_lookup[profile_id], decision_lookup[profile_id]
        gate_table.append({
            "profile_id": profile_id,
            "profile_kind": row["profile_kind"],
            "phaded_superfamily": row["phaded_superfamily"],
            "phaded_family_id": row["phaded_family_id"],
            "qualified_e2_e3_positive_count": detail["qualified_e2_e3_positive_count"],
            "independent_genus_count": detail["independent_genus_count"],
            "required_qualified_e2_e3_positive_count": str(module.MINIMUM_QUALIFIED_E2_E3_POSITIVE_COUNT),
            "required_independent_genus_count": detail["required_independent_genus_count"],
            "positive_count": detail["positive_count"],
            "required_positive_count": detail["required_positive_count"],
            "positive_basis": detail["positive_basis"],
            "heldout_positive": detail["heldout_positive"],
            "negative_count": detail["negative_count"],
            "challenge_count": detail["challenge_count"],
            "family_resolved_negative_status": decision["family_resolved_negative_status"],
            "challenge_status": decision["challenge_status"],
            "gate_kind": detail["gate_kind"],
            "gate_passed": detail["gate_passed"],
            "blocking_slots": detail["blocking_slots"],
            "promotion_authorized": detail["promotion_authorized"],
            "promotion_scope": scope_of[profile_id],
            "readiness_row_origin": row.get("readiness_row_origin", ""),
            "final_decision": decision["decision"],
            "calibration_status": decision["calibration_status"],
            "promoted": "true" if decision["decision"] == PROMOTED else "false",
        })
    _write_rows(output_dir / OUTPUT_GATE_TABLE, list(GATE_TABLE_FIELDS), gate_table)

    scope_map = [{
        "profile_id": row["profile_id"].strip(),
        "phaded_family_id": row["phaded_family_id"],
        "profile_kind": row["profile_kind"],
        "promotion_scope": scope_of[row["profile_id"].strip()],
        "source_invocation": scope_of[row["profile_id"].strip()],
        "source_decisions_file": _display_path(scopes_dir / scope_of[row["profile_id"].strip()] / DECISIONS_NAME),
        "source_gate_detail_file": _display_path(scopes_dir / scope_of[row["profile_id"].strip()] / GATE_DETAIL_NAME),
    } for row in rows]
    _write_rows(output_dir / OUTPUT_SCOPE_MAP, list(SCOPE_MAP_FIELDS), scope_map)

    verification = {
        "schema_version": 1,
        "run_id": run_id,
        "task": "promote_phaded_calibrated_candidate",
        "promotion_authorized": promote_calibrated_model,
        "checks": checks,
        "promoted": promoted,
        "deferred": deferred,
        "authorized_scope": SCOPE_AUTHORIZED,
        "deferred_scope": SCOPE_DEFERRED,
        "control_invocation": {
            "path": _display_path(control_dir),
            "promotion_authorized": False,
            "promoted": [row["profile_id"] for row in control_decisions
                         if row["decision"] == PROMOTED],
        },
        "finalize_module": _describe(FINALIZE_SCRIPT),
        "boundary": (
            "A promotion is a candidate-gate outcome under an explicit authorization. It is "
            "not an experimental phenotype claim and says nothing about any individual hit."),
    }
    (output_dir / OUTPUT_VERIFICATION).write_text(
        json.dumps(verification, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    decisions_path = output_dir / DECISIONS_NAME
    gate_detail_path = output_dir / GATE_DETAIL_NAME
    authorized_gate_block = json.loads(
        (scopes_dir / SCOPE_AUTHORIZED / REPORT_NAME).read_text(encoding="utf-8"))["gate"]
    report = {
        "schema_version": 1,
        "run_id": run_id,
        "task": "promote_phaded_calibrated_candidate",
        "status": "completed_candidate_only",
        "consumes": _describe(FINALIZE_SCRIPT),
        "input": _describe(readiness),
        "output": _describe(decisions_path),
        "gate_detail": _describe(gate_detail_path),
        "scope_invocations": invocations,
        "control_invocation": {
            "path": _display_path(control_dir),
            "promotion_authorized": False,
            "profile_count": len(rows),
        },
        "merge": {
            "rule": "union of the per-scope module outputs, emitted in readiness order",
            "source": "scopes/<promotion_scope>/leaveout_calibration_{decisions,gate_detail}.tsv",
            "verified_row_for_row": True,
            "rewritten_rows": 0,
            "row_count": len(merged["decisions"]),
        },
        "summary": summary,
        "gate": authorized_gate_block,
        "promotion": {
            "authorized": promote_calibrated_model,
            "flag": "--promote-calibrated-model",
            "default": False,
            "authorized_scope": SCOPE_AUTHORIZED,
            "deferred_scope": SCOPE_DEFERRED,
            "promoted": promoted,
            "deferred": deferred,
            "scope_expression": (
                "the module has one global promotion switch and no per-profile deferral, so the "
                "run splits the readiness by promotion_scope and consumes the module once per scope"),
        },
        "verification": _display_path(output_dir / OUTPUT_VERIFICATION),
        "boundary": (
            "Candidate-only sequence evidence. A promotion records that the project's candidate "
            "gate passed and that an authorized promotion was executed; it is not an experimental "
            "phenotype claim and does not mean any individual hit degrades PHB/PHA."),
    }
    (output_dir / REPORT_NAME).write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    (output_dir / STATUS_NAME).write_text(
        f"# PhaDED scoped calibration promotion ({run_id})\n\n"
        f"- Profiles evaluated: {summary['profile_count']} "
        f"({summary['authorized_scope_profile_count']} authorized frame + "
        f"{summary['deferred_scope_profile_count']} deferred)\n"
        f"- `calibrated_candidate_model`: {summary['calibrated_candidate_model']}\n"
        f"- `candidate_gate_passed_not_promoted`: {summary['candidate_gate_passed_not_promoted']}\n"
        f"- `reference_only_insufficient_panel`: {summary['reference_only_insufficient_panel']}\n"
        f"- `planned_not_run_superfamily_requires_family_resolution`: "
        f"{summary['planned_not_run_superfamily_requires_family_resolution']}\n"
        f"- family gate evaluated for: {summary['family_gate_evaluated']} profile(s); "
        f"passed: {summary['profiles_passing_family_gate']}\n"
        f"- superfamily gate evaluated for: {summary['superfamily_gate_evaluated']} profile(s); "
        f"passed: {summary['profiles_passing_superfamily_gate']}\n"
        f"- promotion to `calibrated_candidate_model` authorized: "
        f"{'yes' if promote_calibrated_model else 'no (--promote-calibrated-model not given)'}\n"
        f"- promoted: {', '.join(promoted) if promoted else 'none'}\n"
        f"- deferred by design (gate result printed for transparency): "
        f"{', '.join(item['profile_id'] for item in deferred) if deferred else 'none'}\n"
        "- The governing decision and gate-detail files are the union of the per-scope module "
        "outputs in readiness order; every row was verified against its source row.\n"
        "- A promotion is a candidate-gate outcome, not a phenotype claim, and it does not "
        "mean any individual hit degrades PHB/PHA.\n",
        encoding="utf-8",
    )
    return {"summary": summary, "promoted": promoted, "deferred": deferred,
            "verification": verification, "outputs": {
                "decisions": _describe(decisions_path),
                "gate_detail": _describe(gate_detail_path),
                "report": _describe(output_dir / REPORT_NAME),
                "gate_table": _describe(output_dir / OUTPUT_GATE_TABLE),
                "verification": _describe(output_dir / OUTPUT_VERIFICATION),
            }}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readiness", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--promote-calibrated-model", action="store_true",
                        help="explicit authorization for the authorized_promotion_frame scope")
    parser.add_argument("--defer-profile", action="append", default=[],
                        help="must equal the readiness's deferred_design_spec_special_family set")
    args = parser.parse_args(argv)
    result = promote(args.readiness, args.output_dir, run_id=args.run_id,
                     promote_calibrated_model=args.promote_calibrated_model,
                     defer_profiles=args.defer_profile)
    print(json.dumps({"summary": result["summary"], "promoted": result["promoted"],
                      "deferred": result["deferred"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
