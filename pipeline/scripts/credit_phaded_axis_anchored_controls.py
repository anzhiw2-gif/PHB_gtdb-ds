#!/usr/bin/env python3
"""Axis-anchored credit and superfamily-scope conventions for the PhaDED gate.

Two counting-convention defects are corrected here; **no numeric threshold is moved**:

* **Rule #3 — axis-anchored credit.**  The per-profile gate counters in
  ``reconcile_phaded_reference_panel.py`` and
  ``acquire_phaded_reference_panel_amendment.py`` only sum records bound to the
  profile's own family, so every control that the discrimination panel anchors to a
  *discrimination axis* contributes a constant zero and the 960-record panel is
  invisible to calibration.  This module adds an explicit, separately named credit
  channel (``axis_anchored_negative_count`` / ``axis_anchored_challenge_count``) that
  never feeds the legacy counters.
* **Rule #4 — superfamily scope.**  ``profile_kind == "superfamily"`` is structurally
  excluded by the family gate, so this module also provides the member-family roll-up
  that a superfamily-resolution judgement needs.

Credit is fail-closed.  A record only counts when all of the following hold:

1. ``axis_id`` is present *and* declared in the frozen axis definition set;
2. ``anchor_justification`` is present (auditable anchor statement);
3. ``source_doi_or_pmid`` is present (literature provenance);
4. ``anchor_evidence_kind`` is in that axis/end's pre-registered qualifying list --
   in particular a ``candidate_polar_x1_pattern_state`` can never open the negative
   channel, because a polar x1 is a confounder flag and not a biological negative;
5. the credited profile holds a ``resolved_`` family anchor on the **opposite** end of
   the same axis, i.e. the record and the profile actually form a discrimination pair.

Outputs are candidate-only counting evidence.  Nothing here proves a PHB/PHA
degradation phenotype, and nothing here may be used to filter or delete candidates.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path


AUDIT_FIELDS = ("axis_id", "anchor_justification", "source_doi_or_pmid")
AUDITED = "audited"
NOT_AUDITED = "not_audited"

CHANNEL_NEGATIVE = "axis_anchored_negative"
CHANNEL_CHALLENGE = "axis_anchored_challenge"
CHANNEL_NONE = "no_axis_credit"

REASON_NO_OPPOSING_ANCHOR = "no_opposing_end_resolved_anchor"

RESOLVED_PREFIX = "resolved"
ENDS = ("positive", "negative")
OPPOSITE_END = {"positive": "negative", "negative": "positive"}


def _split_kinds(value: object) -> frozenset[str]:
    """Split a ``;``-separated qualifying-kind declaration into a set of tokens."""
    if value is None:
        return frozenset()
    return frozenset(token.strip() for token in str(value).split(";") if token.strip())


def load_axis_definitions(path: str | Path) -> dict[str, dict[str, object]]:
    """Load the frozen ``axis_definitions`` from an ``axis_reachability.json`` payload."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    definitions = payload.get("axis_definitions")
    if not isinstance(definitions, list) or not definitions:
        raise ValueError(f"axis definitions are missing or empty: {path}")
    by_id: dict[str, dict[str, object]] = {}
    for definition in definitions:
        axis_id = str(definition.get("axis_id", "")).strip()
        if not axis_id:
            raise ValueError(f"axis definition without axis_id: {path}")
        if axis_id in by_id:
            raise ValueError(f"duplicate axis definition: {axis_id}")
        by_id[axis_id] = definition
    return by_id


def load_family_superfamily(path: str | Path) -> dict[str, str]:
    """Return ``family_id -> superfamily`` from a frozen profile manifest."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows or "phaded_family_id" not in (rows[0].keys() if rows else ()):
            raise ValueError(f"profile manifest has no phaded_family_id column: {path}")
    mapping: dict[str, str] = {}
    for row in rows:
        family = (row.get("phaded_family_id") or "").strip()
        superfamily = (row.get("phaded_superfamily") or "").strip()
        if not family or not superfamily:
            continue
        existing = mapping.get(family)
        if existing is not None and existing != superfamily:
            raise ValueError(f"{family}: conflicting superfamily {existing!r} vs {superfamily!r}")
        mapping[family] = superfamily
    if not mapping:
        raise ValueError(f"profile manifest declares no family/superfamily pairs: {path}")
    return mapping


def audit_anchor_basis(record: dict) -> str:
    """Return ``"audited"`` or raise when the axis anchor is not auditable (fail-closed)."""
    missing = [field for field in AUDIT_FIELDS if not str(record.get(field, "")).strip()]
    if missing:
        accession = str(record.get("accession", "(unknown)")).strip() or "(unknown)"
        raise ValueError(
            f"{accession}: axis anchor is not auditable, missing " + ",".join(missing)
        )
    return AUDITED


def classify_record(record: dict, axis_definitions: dict[str, dict[str, object]]) -> tuple[str, str]:
    """Return ``(credit_channel, reason)`` for one panel record, fail-closed."""
    axis_id = str(record.get("axis_id", "")).strip()
    if not axis_id:
        raise ValueError("panel record without axis_id cannot be credited")
    if axis_id not in axis_definitions:
        raise ValueError(f"unknown axis_id: {axis_id}")
    audit_anchor_basis(record)
    definition = axis_definitions[axis_id]

    end = str(record.get("axis_end", "")).strip()
    if end not in ENDS:
        raise ValueError(f"{record.get('accession', '(unknown)')}: unsupported axis_end {end!r}")
    role = str(record.get("role", "")).strip()
    kind = str(record.get("anchor_evidence_kind", "")).strip()
    evidence = str(record.get("evidence_status", "")).strip()

    if evidence == "experimental_negative":
        if role != "negative_anchor" or end != "negative":
            return CHANNEL_NONE, f"negative_evidence_outside_the_negative_end_anchor_slot:{role}:{end}"
        qualifying = _split_kinds(definition.get("negative_end_qualifying_evidence_kinds"))
        if kind not in qualifying:
            return CHANNEL_NONE, f"anchor_evidence_kind_not_qualifying:{kind}"
        return CHANNEL_NEGATIVE, f"negative_end_qualifying_evidence_kind:{kind}"

    if role != "challenge":
        return CHANNEL_NONE, f"role_not_a_discrimination_control:{role}"
    key = "positive_end_challenge_evidence_kinds" if end == "positive" else "negative_end_challenge_evidence_kinds"
    qualifying = _split_kinds(definition.get(key))
    if kind not in qualifying:
        return CHANNEL_NONE, f"anchor_evidence_kind_not_qualifying:{kind}"
    return CHANNEL_CHALLENGE, f"challenge_anchor:{key}:{kind}"


def _target_scopes(
    target_profiles: list[dict],
    family_superfamily: dict[str, str],
) -> tuple[dict[str, str], dict[str, dict[str, object]]]:
    """Split targets into family scopes and superfamily scopes (member families)."""
    family_targets: dict[str, str] = {}
    superfamily_targets: dict[str, dict[str, object]] = {}
    members_by_superfamily: dict[str, list[str]] = defaultdict(list)
    for family, superfamily in family_superfamily.items():
        members_by_superfamily[superfamily].append(family)
    for row in target_profiles:
        profile_id = str(row.get("profile_id", "")).strip()
        if not profile_id:
            raise ValueError("target profile without profile_id")
        kind = str(row.get("profile_kind", "")).strip()
        if kind == "family":
            family = str(row.get("phaded_family_id", "")).strip()
            if not family:
                raise ValueError(f"{profile_id}: family target without phaded_family_id")
            family_targets[profile_id] = family
        elif kind == "superfamily":
            superfamily = str(row.get("phaded_superfamily", "")).strip()
            if not superfamily:
                raise ValueError(f"{profile_id}: superfamily target without phaded_superfamily")
            superfamily_targets[profile_id] = {
                "phaded_superfamily": superfamily,
                "member_family_ids": sorted(members_by_superfamily.get(superfamily, [])),
            }
        else:
            raise ValueError(f"{profile_id}: unsupported target profile_kind {kind!r}")
    return family_targets, superfamily_targets


def attribute_axis_credit(
    panel_records: list[dict],
    axis_definitions: dict[str, dict[str, object]],
    target_profiles: list[dict],
    family_superfamily: dict[str, str],
) -> dict[str, object]:
    """Attribute auditable axis records to the profiles they genuinely discriminate.

    A record at end *e* of axis *a* credits a target profile only when that profile
    holds a ``resolved_`` family anchor on the opposite end of *a*.  Records are never
    double counted, and the legacy family-bound counters are not touched at all.
    """
    family_targets, superfamily_targets = _target_scopes(list(target_profiles), family_superfamily)

    anchor_families: dict[tuple[str, str], set[str]] = defaultdict(set)
    audited: list[dict] = []
    for record in panel_records:
        channel, reason = classify_record(record, axis_definitions)
        binding = str(record.get("family_binding_status", "")).strip()
        family = str(record.get("bound_family", "")).strip()
        if binding.startswith(RESOLVED_PREFIX) and family:
            anchor_families[(str(record.get("axis_id", "")).strip(), str(record.get("axis_end", "")).strip())].add(family)
        audited.append(
            {
                "record": record,
                "channel": channel,
                "reason": reason,
                "binding": binding,
                "family": family,
            }
        )

    profiles: dict[str, dict[str, object]] = {}
    credit_ids: dict[tuple[str, str], list[str]] = defaultdict(list)
    record_rows: list[dict] = []

    for entry in audited:
        record = entry["record"]
        channel = entry["channel"]
        binding = entry["binding"]
        family = entry["family"]
        axis_id = str(record.get("axis_id", "")).strip()
        end = str(record.get("axis_end", "")).strip()
        opposite = OPPOSITE_END[end]
        opposite_families = anchor_families.get((axis_id, opposite), set())

        credited: list[str] = []
        if channel != CHANNEL_NONE:
            for profile_id, target_family in family_targets.items():
                if target_family in opposite_families:
                    credited.append(profile_id)
            for profile_id, scope in superfamily_targets.items():
                if any(member in opposite_families for member in scope["member_family_ids"]):
                    credited.append(profile_id)
            credited = sorted(dict.fromkeys(credited))

        if channel == CHANNEL_NONE:
            credit_reason = entry["reason"]
        elif credited:
            credit_reason = "credited_against_opposite_end_resolved_anchor:" + opposite
        else:
            credit_reason = REASON_NO_OPPOSING_ANCHOR

        if channel == CHANNEL_NEGATIVE and credited:
            for profile_id in credited:
                credit_ids[(profile_id, CHANNEL_NEGATIVE)].append(str(record.get("accession", "")).strip())
        if channel == CHANNEL_CHALLENGE and credited:
            for profile_id in credited:
                credit_ids[(profile_id, CHANNEL_CHALLENGE)].append(str(record.get("accession", "")).strip())

        family_resolved_negative = (
            str(record.get("evidence_status", "")).strip() == "experimental_negative"
            and binding.startswith(RESOLVED_PREFIX)
            and bool(str(record.get("bound_reference_only_profile_id", "")).strip())
        )
        family_resolved_challenge = (
            str(record.get("role", "")).strip() == "challenge"
            and binding.startswith(RESOLVED_PREFIX)
            and bool(str(record.get("bound_reference_only_profile_id", "")).strip())
        )
        record_rows.append(
            {
                "accession": str(record.get("accession", "")).strip(),
                "axis_id": axis_id,
                "axis_end": end,
                "opposite_end": opposite,
                "role": str(record.get("role", "")).strip(),
                "evidence_status": str(record.get("evidence_status", "")).strip(),
                "anchor_evidence_kind": str(record.get("anchor_evidence_kind", "")).strip(),
                "family_binding_status": binding,
                "bound_family": family,
                "bound_reference_only_profile_id": str(record.get("bound_reference_only_profile_id", "")).strip(),
                "anchor_basis_status": AUDITED,
                "credit_channel": channel,
                "credit_reason": credit_reason,
                "credited_profile_ids": credited,
                "credited_profile_count": len(credited),
                "family_resolved_negative_slot_credited": family_resolved_negative,
                "family_resolved_challenge_slot_credited": family_resolved_challenge,
                "axis_anchored_slot_credited": bool(credited),
                "source_doi_or_pmid": str(record.get("source_doi_or_pmid", "")).strip(),
            }
        )

    for profile_id in list(family_targets) + list(superfamily_targets):
        negatives = sorted(credit_ids.get((profile_id, CHANNEL_NEGATIVE), []))
        challenges = sorted(credit_ids.get((profile_id, CHANNEL_CHALLENGE), []))
        profiles[profile_id] = {
            "profile_id": profile_id,
            "profile_kind": "family" if profile_id in family_targets else "superfamily",
            "axis_anchored_negative_count": len(negatives),
            "axis_anchored_challenge_count": len(challenges),
            "axis_anchored_negative_credit_ids": negatives,
            "axis_anchored_challenge_credit_ids": challenges,
            "axis_anchor_basis_status": AUDITED if (negatives or challenges) else NOT_AUDITED,
        }

    for profile_id in superfamily_targets:
        profiles[profile_id]["member_family_ids"] = superfamily_targets[profile_id]["member_family_ids"]

    reason_counts: dict[str, int] = defaultdict(int)
    for row in record_rows:
        reason_counts[row["credit_reason"]] += 1
    return {
        "records": record_rows,
        "profiles": profiles,
        "counts": {
            "panel_records": len(record_rows),
            "records_with_axis_credit": sum(bool(row["axis_anchored_slot_credited"]) for row in record_rows),
            "records_without_axis_credit": sum(not row["axis_anchored_slot_credited"] for row in record_rows),
            "negative_channel_records": sum(row["credit_channel"] == CHANNEL_NEGATIVE for row in record_rows),
            "challenge_channel_records": sum(row["credit_channel"] == CHANNEL_CHALLENGE for row in record_rows),
            "no_channel_records": sum(row["credit_channel"] == CHANNEL_NONE for row in record_rows),
            "family_resolved_negative_slot_records": sum(
                row["family_resolved_negative_slot_credited"] for row in record_rows
            ),
            "credit_reason_counts": dict(sorted(reason_counts.items())),
            "target_profiles": len(profiles),
        },
    }


def rollup_superfamily_evidence(
    readiness_rows: list[dict],
    family_superfamily: dict[str, str],
) -> dict[str, dict[str, object]]:
    """Sum member-family evidence for every superfamily-resolution target.

    A superfamily *is* the union of its member families in the frozen layering
    authority, so a member family's independent positive is a superfamily-resolution
    positive by construction.  The roll-up is reported side by side with the raw
    family-bound columns, never instead of them, and it never invents evidence:
    member families without a readiness row contribute 0 and are listed explicitly.
    """
    by_family: dict[str, dict] = {}
    for row in readiness_rows:
        family = str(row.get("phaded_family_id", "")).strip()
        if family and str(row.get("profile_kind", "")).strip() == "family":
            by_family[family] = row
    members_by_superfamily: dict[str, list[str]] = defaultdict(list)
    for family, superfamily in family_superfamily.items():
        members_by_superfamily[superfamily].append(family)

    rollup: dict[str, dict[str, object]] = {}
    for row in readiness_rows:
        if str(row.get("profile_kind", "")).strip() != "superfamily":
            continue
        profile_id = str(row.get("profile_id", "")).strip()
        superfamily = str(row.get("phaded_superfamily", "")).strip()
        members = sorted(members_by_superfamily.get(superfamily, []))
        positives = 0
        heldout = 0
        missing: list[str] = []
        for member in members:
            member_row = by_family.get(member)
            if member_row is None:
                missing.append(member)
                continue
            positives += _as_int(member_row, "existing_experimental_positive") + _as_int(member_row, "external_bound_positive")
            heldout += _as_int(member_row, "heldout_positive")
        rollup[profile_id] = {
            "profile_id": profile_id,
            "phaded_superfamily": superfamily,
            "member_family_ids": members,
            "member_families_without_readiness_row": missing,
            "superfamily_positive_count": positives,
            "superfamily_heldout_positive": heldout,
        }
    return rollup


def _as_int(row: dict, field: str) -> int:
    raw = str(row.get(field, "0") or "0").strip() or "0"
    try:
        value = int(raw)
    except ValueError as error:
        raise ValueError(f"{row.get('profile_id', '(unknown)')}: {field} is not an integer") from error
    if value < 0:
        raise ValueError(f"{row.get('profile_id', '(unknown)')}: {field} is negative")
    return value


__all__ = [
    "AUDITED",
    "AUDIT_FIELDS",
    "CHANNEL_CHALLENGE",
    "CHANNEL_NEGATIVE",
    "CHANNEL_NONE",
    "NOT_AUDITED",
    "REASON_NO_OPPOSING_ANCHOR",
    "attribute_axis_credit",
    "audit_anchor_basis",
    "classify_record",
    "load_axis_definitions",
    "load_family_superfamily",
    "rollup_superfamily_evidence",
]
