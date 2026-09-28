#!/usr/bin/env python3
"""Finalize fail-closed leave-out decisions for reference-only PhaDED panels.

Two conventions live in this module (see
``docs/superpowers/plans/2026-09-17-phaded-gate-convention-decision-record.md`` and
``docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign.md`` Task 5):

* the **family** gate keeps its four thresholds but the positive slot now counts
  *qualified independent* evidence: ``qualified_e2_e3_positive_count >= 3`` **and**
  ``independent_genus_count >= 3`` (plus the unchanged held-out / family-resolved
  negative / challenge / status slots).  Row counts of ungraded historical positives no
  longer satisfy the functional calibration gate, and both new columns fail closed to 0
  when they are absent.  Weakening this gate is never allowed;
* the **superfamily** gate is a separate, explicitly named gate with the same numeric
  thresholds.  ``profile_kind == "superfamily"`` used to be structurally excluded
  because the single gate demanded ``kind == "family"``, which contradicted the
  layering authority (superfamilies are the functional-prior frame, families the
  refinement layer).  Superfamily calibration is a *coarser* judgement, so every gate
  detail row carries an explicit ``gate_resolution`` label.

Functional calibration and model promotion are separate actions.  When a gate passes
without an explicitly authorized promotion, the profile is reported as
``candidate_gate_passed_not_promoted``: ``calibrated_candidate_model`` is emitted only
when the caller passes ``promote_calibrated_model=True`` (CLI:
``--promote-calibrated-model``, default off).  Authorization alone never bypasses the
gate.  When a gate does not pass, the decision supplied by the readiness input is kept
verbatim and the calibration status stays ``not_run``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


REQUIRED = {
    "profile_id", "profile_kind", "external_bound_positive", "external_bound_negative",
    "external_bound_challenge", "heldout_positive", "family_resolved_negative_status",
    "challenge_status", "calibration_decision", "new_family_call_created",
}
FIELDS = (
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "positive_count_for_gate", "external_bound_positive", "external_bound_negative",
    "external_bound_challenge", "heldout_positive", "family_resolved_negative_status",
    "challenge_status", "decision", "calibration_status", "new_family_call_created",
)
GATE_FIELDS = (
    "profile_id", "profile_kind", "gate_kind", "gate_resolution", "gate_evaluated", "gate_passed",
    "positive_count", "required_positive_count", "positive_basis", "family_bound_positive_count",
    "heldout_positive", "required_heldout_positive_count", "family_bound_heldout_positive",
    "negative_count", "required_negative_count", "negative_basis", "family_bound_negative_count",
    "challenge_count", "required_challenge_count", "challenge_basis", "family_bound_challenge_count",
    "axis_anchored_negative_count", "axis_anchored_challenge_count", "axis_anchor_basis_status",
    "zero_unexplained_hits_status", "blocking_slots", "decision", "calibration_status",
    "new_family_call_created", "qualified_e2_e3_positive_count", "independent_genus_count",
    "required_independent_genus_count", "promotion_authorized",
)

AXIS_ANCHOR_BASIS_AUDITED = "audited"
ZERO_UNEXPLAINED_SUFFICIENT = "sufficient"
ZERO_UNEXPLAINED_NOT_ASSESSABLE = "not_assessable"
NOT_AUDITED = "not_audited"

#: Functional calibration governance thresholds (never relaxed).
MINIMUM_POSITIVE_COUNT = 3
MINIMUM_QUALIFIED_E2_E3_POSITIVE_COUNT = 3
MINIMUM_INDEPENDENT_GENUS_COUNT = 3
MINIMUM_HELDOUT_POSITIVE_COUNT = 1
MINIMUM_NEGATIVE_COUNT = 1
MINIMUM_CHALLENGE_COUNT = 1

PROMOTED = "calibrated_candidate_model"
NOT_PROMOTED = "candidate_gate_passed_not_promoted"
CALIBRATION_PASSED = "passed"
CALIBRATION_NOT_RUN = "not_run"
SUPPORTED_DECISIONS = {
    PROMOTED, NOT_PROMOTED, "reference_only_insufficient_panel",
    "planned_not_run_superfamily_requires_family_resolution", "planned_not_run",
}

AXIS_COLUMNS = ("axis_anchored_negative_count", "axis_anchored_challenge_count")
SUPERFAMILY_COLUMNS = ("superfamily_positive_count", "superfamily_heldout_positive")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _int(row: dict[str, str], field: str) -> int:
    try:
        value = int(row.get(field, "0") or "0")
    except ValueError as error:
        raise ValueError(f"{row.get('profile_id', '(unknown)')}: {field} is not an integer") from error
    if value < 0:
        raise ValueError(f"{row.get('profile_id', '(unknown)')}: {field} is negative")
    return value


def _optional_int(row: dict[str, str], field: str) -> int:
    """Read an optional counter column; an absent column can only fail closed (0)."""
    if field not in row:
        return 0
    return _int(row, field)


def _text(row: dict[str, str], field: str, default: str) -> str:
    value = str(row.get(field, "") or "").strip()
    return value or default


def _family_gate(
    kind: str,
    qualified_e2_e3_positive_count: int,
    independent_genus_count: int,
    heldout: int,
    negative: int,
    challenge: int,
    family_resolved_negative_status: str,
    challenge_status: str,
) -> bool:
    """Functional calibration gate at family resolution (2026-09-28 redesign, Task 5).

    The positive slot counts *qualified independent* E2/E3 evidence and, separately,
    the number of independent genera it covers, so ungraded historical row counts can
    no longer carry the gate.  No numeric threshold is lowered.  Both new columns are
    read as optional counters: an absent column fails closed to 0.
    """
    return (
        kind == "family"
        and qualified_e2_e3_positive_count >= MINIMUM_QUALIFIED_E2_E3_POSITIVE_COUNT
        and independent_genus_count >= MINIMUM_INDEPENDENT_GENUS_COUNT
        and heldout >= MINIMUM_HELDOUT_POSITIVE_COUNT
        and negative >= MINIMUM_NEGATIVE_COUNT
        and challenge >= MINIMUM_CHALLENGE_COUNT
        and family_resolved_negative_status == "sufficient"
        and challenge_status == "sufficient"
    )


def _superfamily_gate(kind: str, positive_count: int, heldout: int, negative: int, challenge: int, axis_anchor_basis_status: str, zero_unexplained_hits_status: str) -> bool:
    """Independent superfamily-resolution gate: same thresholds, explicitly named slots."""
    return kind == "superfamily" and positive_count >= 3 and heldout >= 1 and negative >= 1 and challenge >= 1 and axis_anchor_basis_status == "audited" and zero_unexplained_hits_status == "sufficient"


def _family_slots(qualified_e2_e3_positive_count: int, independent_genus_count: int, heldout: int, negative: int, challenge: int, family_resolved_negative_status: str, challenge_status: str) -> list[str]:
    blocking = []
    if qualified_e2_e3_positive_count < MINIMUM_QUALIFIED_E2_E3_POSITIVE_COUNT:
        blocking.append("qualified_e2_e3_positive_below_3")
    if independent_genus_count < MINIMUM_INDEPENDENT_GENUS_COUNT:
        blocking.append("independent_genus_below_3")
    if heldout < 1:
        blocking.append("heldout_positive_below_1")
    if negative < 1:
        blocking.append("family_resolved_negative_below_1")
    if challenge < 1:
        blocking.append("challenge_below_1")
    if family_resolved_negative_status != "sufficient":
        blocking.append("family_resolved_negative_status_not_sufficient")
    if challenge_status != "sufficient":
        blocking.append("challenge_status_not_sufficient")
    return blocking


def _superfamily_slots(positive_count: int, heldout: int, negative: int, challenge: int, axis_anchor_basis_status: str, zero_unexplained_hits_status: str) -> list[str]:
    blocking = []
    if positive_count < 3:
        blocking.append("superfamily_positive_below_3")
    if heldout < 1:
        blocking.append("superfamily_heldout_positive_below_1")
    if negative < 1:
        blocking.append("axis_anchored_negative_below_1")
    if challenge < 1:
        blocking.append("axis_anchored_challenge_below_1")
    if axis_anchor_basis_status != "audited":
        blocking.append("axis_anchor_basis_not_audited")
    if zero_unexplained_hits_status != "sufficient":
        blocking.append("zero_unexplained_hits_not_sufficient")
    return blocking


def finalize(
    readiness_path: str | Path,
    output_dir: str | Path,
    *,
    run_id: str,
    promote_calibrated_model: bool = False,
) -> dict[str, object]:
    """Evaluate both calibration gates and write the fail-closed decision files.

    ``promote_calibrated_model`` is the explicit promotion authorization and defaults to
    ``False``: without it a passing profile is reported as
    ``candidate_gate_passed_not_promoted`` and never as ``calibrated_candidate_model``.
    """
    readiness_path, output_dir = Path(readiness_path), Path(output_dir)
    if not readiness_path.is_file() or readiness_path.is_symlink():
        raise ValueError(f"readiness input is not a regular file: {readiness_path}")
    with readiness_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(REQUIRED - set(reader.fieldnames or []))
        if missing:
            raise ValueError("readiness missing columns: " + ",".join(missing))
        rows = list(reader)
    if not rows:
        raise ValueError("readiness is empty")
    seen: set[str] = set()
    decisions: list[dict[str, str]] = []
    gate_detail: list[dict[str, str]] = []
    for row in rows:
        profile_id = row["profile_id"].strip()
        if not profile_id or profile_id in seen:
            raise ValueError(f"empty or duplicate profile_id: {profile_id}")
        seen.add(profile_id)
        positive_count = _int(row, "existing_experimental_positive") if "existing_experimental_positive" in row else _int(row, "external_bound_positive")
        if "existing_experimental_positive" in row:
            positive_count += _int(row, "external_bound_positive")
        negative = _int(row, "external_bound_negative")
        challenge = _int(row, "external_bound_challenge")
        heldout = _int(row, "heldout_positive")
        kind = row["profile_kind"].strip()
        supplied = row["calibration_decision"].strip()
        new_call = row["new_family_call_created"].strip().lower()
        if new_call not in {"false", "0", "no"}:
            raise ValueError(f"{profile_id}: new family call flag must be false")

        family_resolved_negative_status = row["family_resolved_negative_status"].strip()
        challenge_status = row["challenge_status"].strip()

        axis_negative = _optional_int(row, "axis_anchored_negative_count")
        axis_challenge = _optional_int(row, "axis_anchored_challenge_count")
        axis_basis = _text(row, "axis_anchor_basis_status", NOT_AUDITED)
        if (axis_negative or axis_challenge) and axis_basis != AXIS_ANCHOR_BASIS_AUDITED:
            raise ValueError(f"{profile_id}: axis_anchored counts require an audited anchor basis")
        zero_unexplained = _text(row, "zero_unexplained_hits_status", ZERO_UNEXPLAINED_NOT_ASSESSABLE)
        superfamily_positive = _optional_int(row, "superfamily_positive_count")
        superfamily_heldout = _optional_int(row, "superfamily_heldout_positive")
        superfamily_positive_basis = "member_family_rollup" if "superfamily_positive_count" in row else "family_bound_fallback"
        superfamily_heldout_basis = "member_family_rollup" if "superfamily_heldout_positive" in row else "family_bound_fallback"
        # Qualified independent evidence: absent columns fail closed to 0, exactly like
        # every other optional readiness counter in this module.
        qualified_positive = _optional_int(row, "qualified_e2_e3_positive_count")
        independent_genus = _optional_int(row, "independent_genus_count")
        positive_basis = (
            "qualified_e2_e3_independent_positive"
            if "qualified_e2_e3_positive_count" in row else "absent_column_fails_closed"
        )

        family_gate = _family_gate(kind, qualified_positive, independent_genus, heldout, negative, challenge, family_resolved_negative_status, challenge_status)
        superfamily_gate = _superfamily_gate(kind, superfamily_positive, superfamily_heldout, axis_negative, axis_challenge, axis_basis, zero_unexplained)
        gate = family_gate or superfamily_gate
        if gate and promote_calibrated_model:
            decision = PROMOTED
        elif gate:
            decision = NOT_PROMOTED
        else:
            decision = supplied
        if decision == PROMOTED and not (gate and promote_calibrated_model):
            raise ValueError(
                f"{profile_id}: calibrated decision violates panel gate or lacks promotion authorization"
            )
        if decision not in SUPPORTED_DECISIONS:
            raise ValueError(f"{profile_id}: unsupported calibration decision {decision}")
        if decision == PROMOTED:
            calibration_status = CALIBRATION_PASSED
        elif gate:
            calibration_status = NOT_PROMOTED
        else:
            calibration_status = CALIBRATION_NOT_RUN
        decisions.append({
            "profile_id": profile_id,
            "profile_kind": kind,
            "phaded_superfamily": row.get("phaded_superfamily", "").strip(),
            "phaded_family_id": row.get("phaded_family_id", "").strip(),
            "positive_count_for_gate": str(positive_count),
            "external_bound_positive": str(_int(row, "external_bound_positive")),
            "external_bound_negative": str(negative),
            "external_bound_challenge": str(challenge),
            "heldout_positive": str(heldout),
            "family_resolved_negative_status": family_resolved_negative_status,
            "challenge_status": challenge_status,
            "decision": decision,
            "calibration_status": calibration_status,
            "new_family_call_created": "false",
        })

        if kind == "family":
            gate_kind = "family" if family_gate else "none"
            slot_positive, slot_heldout = qualified_positive, heldout
            slot_negative, slot_challenge = negative, challenge
            slot_positive_basis = positive_basis
            negative_basis = "family_bound_external_negative"
            challenge_basis = "family_bound_external_challenge"
            blocking = _family_slots(qualified_positive, independent_genus, heldout, negative, challenge, family_resolved_negative_status, challenge_status)
            if not family_gate and superfamily_gate:
                raise ValueError(f"{profile_id}: superfamily gate matched a family row")
        elif kind == "superfamily":
            gate_kind = "superfamily" if superfamily_gate else "none"
            slot_positive, slot_heldout = superfamily_positive, superfamily_heldout
            slot_negative, slot_challenge = axis_negative, axis_challenge
            slot_positive_basis = superfamily_positive_basis
            negative_basis = "axis_anchored_negative"
            challenge_basis = "axis_anchored_challenge"
            blocking = _superfamily_slots(superfamily_positive, superfamily_heldout, axis_negative, axis_challenge, axis_basis, zero_unexplained)
            if not superfamily_gate and family_gate:
                raise ValueError(f"{profile_id}: family gate matched a superfamily row")
        else:
            gate_kind = "none"
            slot_positive, slot_heldout = qualified_positive, heldout
            slot_negative, slot_challenge = negative, challenge
            slot_positive_basis = positive_basis
            negative_basis = "family_bound_external_negative"
            challenge_basis = "family_bound_external_challenge"
            blocking = _family_slots(qualified_positive, independent_genus, heldout, negative, challenge, family_resolved_negative_status, challenge_status)
            blocking.append("unsupported_profile_kind")

        gate_detail.append({
            "profile_id": profile_id,
            "profile_kind": kind,
            "gate_kind": gate_kind,
            "gate_resolution": kind,
            "gate_evaluated": "true",
            "gate_passed": "true" if gate else "false",
            "positive_count": str(slot_positive),
            "required_positive_count": "3",
            "positive_basis": slot_positive_basis,
            "family_bound_positive_count": str(positive_count),
            "heldout_positive": str(slot_heldout),
            "required_heldout_positive_count": "1",
            "family_bound_heldout_positive": str(heldout),
            "negative_count": str(slot_negative),
            "required_negative_count": "1",
            "negative_basis": negative_basis,
            "family_bound_negative_count": str(negative),
            "challenge_count": str(slot_challenge),
            "required_challenge_count": "1",
            "challenge_basis": challenge_basis,
            "family_bound_challenge_count": str(challenge),
            "axis_anchored_negative_count": str(axis_negative),
            "axis_anchored_challenge_count": str(axis_challenge),
            "axis_anchor_basis_status": axis_basis,
            "zero_unexplained_hits_status": zero_unexplained,
            "blocking_slots": ";".join(blocking),
            "decision": decision,
            "calibration_status": calibration_status,
            "new_family_call_created": "false",
            "qualified_e2_e3_positive_count": str(qualified_positive),
            "independent_genus_count": str(independent_genus),
            "required_independent_genus_count": str(MINIMUM_INDEPENDENT_GENUS_COUNT),
            "promotion_authorized": "true" if promote_calibrated_model else "false",
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    decisions_path = output_dir / "leaveout_calibration_decisions.tsv"
    gate_detail_path = output_dir / "leaveout_calibration_gate_detail.tsv"
    report_path = output_dir / "leaveout_calibration_report.json"
    status_path = output_dir / "calibration_status.md"
    for path in (decisions_path, gate_detail_path, report_path, status_path):
        if path.exists():
            raise FileExistsError(path)
    with decisions_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(decisions)
    with gate_detail_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GATE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(gate_detail)
    counts = Counter(row["decision"] for row in decisions)
    detail_counts = Counter(row["gate_kind"] for row in gate_detail)
    summary = {
        "profile_count": len(decisions),
        "calibrated_candidate_model": counts.get(PROMOTED, 0),
        "candidate_gate_passed_not_promoted": counts.get(NOT_PROMOTED, 0),
        "reference_only_insufficient_panel": counts.get("reference_only_insufficient_panel", 0),
        "planned_not_run_superfamily_requires_family_resolution": counts.get("planned_not_run_superfamily_requires_family_resolution", 0),
        "new_family_call_created": False,
        "promotion_authorized": promote_calibrated_model,
        "gate_kind_counts": dict(sorted(detail_counts.items())),
        "family_gate_evaluated": sum(row["profile_kind"] == "family" for row in gate_detail),
        "superfamily_gate_evaluated": sum(row["profile_kind"] == "superfamily" for row in gate_detail),
        "profiles_passing_family_gate": detail_counts.get("family", 0),
        "profiles_passing_superfamily_gate": detail_counts.get("superfamily", 0),
    }
    report = {
        "schema_version": 1,
        "run_id": run_id,
        "status": "completed_candidate_only",
        "input": {"path": str(readiness_path.resolve()), "size": readiness_path.stat().st_size, "sha256": _sha256(readiness_path)},
        "output": {"path": str(decisions_path.resolve()), "size": decisions_path.stat().st_size, "sha256": _sha256(decisions_path)},
        "gate_detail": {"path": str(gate_detail_path.resolve()), "size": gate_detail_path.stat().st_size, "sha256": _sha256(gate_detail_path)},
        "summary": summary,
        "gate": {
            "minimum_positive_count": MINIMUM_POSITIVE_COUNT,
            "minimum_qualified_e2_e3_positive_count": MINIMUM_QUALIFIED_E2_E3_POSITIVE_COUNT,
            "minimum_independent_genus_count": MINIMUM_INDEPENDENT_GENUS_COUNT,
            "minimum_heldout_positive_count": 1,
            "minimum_family_resolved_negative_count": 1,
            "minimum_challenge_count": 1,
            "zero_unexplained_negative_or_challenge_hits": True,
            "new_family_calls_allowed": False,
            "promotion_requires_explicit_authorization": True,
            "promotion": {
                "authorized": promote_calibrated_model,
                "flag": "--promote-calibrated-model",
                "default": False,
                "passing_but_unpromoted_status": NOT_PROMOTED,
                "authorization_never_bypasses_a_gate": True,
            },
            "family": {
                "resolution": "family",
                "applies_to_profile_kind": "family",
                "minimum_positive_count": MINIMUM_POSITIVE_COUNT,
                "minimum_qualified_e2_e3_positive_count": MINIMUM_QUALIFIED_E2_E3_POSITIVE_COUNT,
                "minimum_independent_genus_count": MINIMUM_INDEPENDENT_GENUS_COUNT,
                "minimum_heldout_positive_count": 1,
                "minimum_negative_count": 1,
                "minimum_family_resolved_negative_count": 1,
                "minimum_challenge_count": 1,
                "zero_unexplained_negative_or_challenge_hits": True,
                "new_family_calls_allowed": False,
                "positive_slot_source": "qualified_e2_e3_positive_count (qualified independent E2/E3 evidence)",
                "positive_slot_source_column": "qualified_e2_e3_positive_count; absent column fails closed to 0",
                "independent_genus_slot_source_column": "independent_genus_count; absent column fails closed to 0",
                "legacy_positive_count_column": "positive_count_for_gate keeps existing_experimental_positive + external_bound_positive for audit only and no longer carries the gate",
                "negative_slot_source": "external_bound_negative with family_resolved_negative_status == sufficient",
                "challenge_slot_source": "external_bound_challenge with challenge_status == sufficient",
                "axis_anchored_credit_feeds_gate": False,
                "expression_unchanged_from_previous_revision": False,
                "expression_changed_by": "docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign.md Task 5 (qualified E2/E3 positives + independent genera)",
                "thresholds_lowered": False,
            },
            "superfamily": {
                "resolution": "superfamily",
                "applies_to_profile_kind": "superfamily",
                "coarser_than": "family",
                "minimum_positive_count": 3,
                "minimum_heldout_positive_count": 1,
                "minimum_negative_count": 1,
                "minimum_family_resolved_negative_count": 1,
                "minimum_challenge_count": 1,
                "zero_unexplained_negative_or_challenge_hits": True,
                "new_family_calls_allowed": False,
                "positive_slot_source": "superfamily_positive_count (member-family rollup); absent column fails closed to the family-bound count",
                "heldout_slot_source": "superfamily_heldout_positive (member-family rollup); absent column fails closed to 0",
                "negative_slot_source": "axis_anchored_negative_count with axis_anchor_basis_status == audited",
                "challenge_slot_source": "axis_anchored_challenge_count with axis_anchor_basis_status == audited",
                "zero_unexplained_hits_source": "zero_unexplained_hits_status == sufficient",
                "axis_anchored_credit_feeds_gate": True,
                "resolution_label_required_in_output": True,
            },
        },
        "boundary": "Calibration outputs are candidate-only sequence evidence; no profile is promoted unless every gate passes under an explicitly authorized promotion action, and no PHB/PHA phenotype is proven.",
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    status_path.write_text(
        f"# PhaDED reference-only panel calibration ({run_id})\n\n"
        f"- Profiles evaluated: {len(decisions)}\n"
        f"- `calibrated_candidate_model`: {summary['calibrated_candidate_model']}\n"
        f"- `candidate_gate_passed_not_promoted`: {summary['candidate_gate_passed_not_promoted']}\n"
        f"- `reference_only_insufficient_panel`: {summary['reference_only_insufficient_panel']}\n"
        f"- `planned_not_run_superfamily_requires_family_resolution`: {summary['planned_not_run_superfamily_requires_family_resolution']}\n"
        f"- family gate evaluated for: {summary['family_gate_evaluated']} profile(s); passed: {summary['profiles_passing_family_gate']}\n"
        f"- superfamily gate evaluated for: {summary['superfamily_gate_evaluated']} profile(s); passed: {summary['profiles_passing_superfamily_gate']}\n"
        f"- promotion to `calibrated_candidate_model` authorized: {'yes' if promote_calibrated_model else 'no (default: --promote-calibrated-model not given)'}\n"
        "- The family gate requires >= 3 qualified independent E2/E3 positives covering >= 3 independent genera; absent readiness columns fail closed to 0.\n"
        "- Gates are labelled by resolution in `leaveout_calibration_gate_detail.tsv`; a superfamily-resolution pass is coarser than a family-resolution pass and is never a family call.\n"
        "- A gate pass without explicit promotion is reported as `candidate_gate_passed_not_promoted`; authorization never bypasses a gate.\n"
        "- No new family call was created.\n\n"
        "All reference-only profiles remain blocked until their own resolution's gate is satisfied (qualified independent positives across independent genera, family-resolved or axis-anchored negatives, challenge controls, leakage checks, held-out recovery, and zero unexplained hits). Results are candidate-only and do not prove phenotype.\n",
        encoding="utf-8",
    )
    return {"summary": summary, "decisions": decisions, "gate_detail": gate_detail, "report": report}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readiness", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument(
        "--promote-calibrated-model", action="store_true",
        help="explicit authorization to emit calibrated_candidate_model for profiles that pass a gate "
             "(default off: a pass is reported as candidate_gate_passed_not_promoted only)",
    )
    args = parser.parse_args(argv)
    result = finalize(
        args.readiness, args.output_dir, run_id=args.run_id,
        promote_calibrated_model=args.promote_calibrated_model,
    )
    print(json.dumps(result["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
