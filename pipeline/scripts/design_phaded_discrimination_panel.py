#!/usr/bin/env python3
"""Design an axis-anchored PhaDED challenge panel (Task 5).

The panel gives the 37 reference-only PhaDED profiles negative controls that are
anchored on one of four declared discrimination axes, instead of the previous
unanchored challenge set.  This stage only *consumes* frozen evidence: it runs no
HMM, retrieves no sequence, and creates no family call.

Hard invariants (fail-closed, enforced by ``validate_panel_invariants``):

* every record binds to exactly one declared axis and one axis end, and records
  its anchoring basis (accession, evidence layer, why it sits at that end,
  source DOI/PMID);
* ``annotation_only`` records are never written as formal negatives;
* a formal family-resolved negative requires formal experimental-negative
  evidence *and* an accession binding to an existing family that is represented
  by a reference-only profile;
* records that cannot be bound to an existing family stay ``unresolved``;
* the calibration gate (>=3 independent positives, >=1 held-out positive, >=1
  family-resolved negative, >=1 challenge, zero unexplained hits) is reproduced
  verbatim from ``finalize_phaded_reference_panel_calibration.py`` and no
  threshold is modified anywhere in this stage.

Known evidence-layer defect recorded, not fixed, by this stage: in
``runs/20260916_phaded_reference_panel_acquisition_03/inputs/external_panel_candidates.tsv``
the data rows carry one field more than the header, so every field from
``bound_family`` onward is shifted (``pmid`` empty, ``doi`` holds a PMID,
``source_database`` holds a DOI).  PMID/DOI citations in this run therefore come
from the unshifted ``external_panel_evidence_prior.tsv`` or from the acquisition
report, never from the shifted columns.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


RUN_ID = "20260917_phaded_discrimination_panel_01"

AXIS_IDS = (
    "axis_1_substrate_chain_length",
    "axis_2_particle_state",
    "axis_3_localization",
    "axis_4_fold_neighborhood",
)
AXIS_ENDS = ("positive", "negative")
ROLES = ("positive_anchor", "negative_anchor", "challenge")
EVIDENCE_STATUSES = ("experimental_positive", "experimental_negative", "annotation_only")
RECORD_ORIGINS = (
    "reference_ledger",
    "external_acquisition",
    "candidate_universe",
)
ANCHOR_EVIDENCE_KINDS = (
    "experimental_positive_assay",
    "experimental_negative_assay",
    "supported_domain_state",
    "supported_pattern_state",
    "reported_label",
    "candidate_polar_x1_pattern_state",
)
BINDING_STATUSES = (
    "resolved_existing_family",
    "resolved_existing_family_versionless",
    "ambiguous_existing_family",
    "unresolved",
)
RESOLVED_BINDING_STATUSES = ("resolved_existing_family", "resolved_existing_family_versionless")

ASSIGNMENT_RULES = (
    "R1_lid_detected",
    "R2_mcl_substrate",
    "R3_phaz7_like_superfamily",
    "R4_denatured_experimental_positive",
    "R5_scl_experimental_positive",
    "R6_pf10503_domain_state",
    "R7_lipase_fold_domain_state",
    "R8_candidate_polar_x1",
    "R9_external_declared_anchor",
)

ALPHA_BETA_HYDROLASE_PFAMS = ("PF00326", "PF00756", "PF00561", "PF12146", "PF12697")
PF10503 = "PF10503"
PF06850 = "PF06850"

LID_POSITIVE_STATES = ("supported", "partial_lid_loop1_only")
PHAZ7_LIKE_SUPERFAMILY = "extracellular native-SCL/PhaZ7-like"

# Reproduced verbatim from finalize_phaded_reference_panel_calibration.py L75.
GATE = {
    "minimum_positive_count": 3,
    "minimum_heldout_positive_count": 1,
    "minimum_family_resolved_negative_count": 1,
    "minimum_challenge_count": 1,
    "zero_unexplained_negative_or_challenge_hits": True,
    "source_file": "pipeline/scripts/finalize_phaded_reference_panel_calibration.py",
    "source_symbol": "gate",
    "gate_relaxation_applied": False,
    "gate_relaxation_supported": False,
}

PANEL_FIELDS = (
    "axis_id",
    "axis_name",
    "axis_end",
    "role",
    "accession",
    "record_origin",
    "evidence_status",
    "anchor_evidence_kind",
    "anchor_evidence_basis",
    "family_binding_status",
    "bound_family",
    "bound_reference_only_profile_id",
    "formal_negative_eligible",
    "gate_credit_slot",
    "anchor_justification",
    "source_doi_or_pmid",
    "evidence_source_path",
    "secondary_axis_capability",
)

AXIS_DEFINITION_FIELDS = (
    "axis_id",
    "axis_name",
    "positive_end_state",
    "negative_end_state",
    "positive_end_qualifying_evidence_kinds",
    "negative_end_qualifying_evidence_kinds",
    "positive_end_challenge_evidence_kinds",
    "negative_end_challenge_evidence_kinds",
    "positive_end_anchor_status",
    "negative_end_anchor_status",
    "axis_boundary_note",
    "source_doi_or_pmid",
)
ASSIGNMENT_RULE_FIELDS = (
    "rule_id",
    "priority",
    "source_layer",
    "axis_id",
    "axis_end",
    "qualifying_evidence_kind",
    "condition_field",
    "condition_operator",
    "condition_value",
    "rule_description",
)
DECLARATION_FIELDS = (
    "accession",
    "include",
    "axis_id",
    "axis_end",
    "anchor_evidence_kind",
    "evidence_status",
    "family_binding_status_expected",
    "anchor_justification",
    "source_doi_or_pmid",
    "evidence_source_path",
    "exclusion_reason",
)

BOUNDARY = (
    "Axis-anchored controls are candidate-only sequence/domain/pattern evidence. "
    "They do not establish PHB/PHA degradation phenotype, and no threshold of the "
    "reference-only calibration gate is modified by this stage."
)


def builtin_axis_definitions() -> dict[str, dict[str, str]]:
    """Return the canonical in-code axis definitions used to validate the input table.

    For every axis end the definition declares which evidence kinds anchor that end
    (they produce ``positive_anchor`` / ``negative_anchor``) and which kinds may sit
    there only as ``challenge`` controls.  A kind that is not declared for an end is
    rejected outright by :func:`check_declared_kind_for_end`.
    """
    return {
        "axis_1_substrate_chain_length": {
            "axis_id": "axis_1_substrate_chain_length",
            "axis_name": "substrate chain length (SCL/PHB versus MCL)",
            "positive_end_state": "SCL_PHB_substrate",
            "negative_end_state": "MCL_non_PHB_substrate",
            "positive_end_qualifying_evidence_kinds": "experimental_positive_assay",
            "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
            "positive_end_challenge_evidence_kinds": "",
            "negative_end_challenge_evidence_kinds": "reported_label;experimental_positive_assay",
            "positive_end_anchor_status": "anchored_experimental_family_unresolved",
            "negative_end_anchor_status": "anchored_experimental_family_unresolved",
            "axis_boundary_note": "label-only MCL records are challenges, not formal negatives",
            "source_doi_or_pmid": "10.1186/1471-2105-10-89 (PMID 19296857)",
        },
        "axis_2_particle_state": {
            "axis_id": "axis_2_particle_state",
            "axis_name": "particle state (denatured-only extracellular versus native-active PhaZ7 type)",
            "positive_end_state": "denatured_PHA_only_extracellular",
            "negative_end_state": "native_active_PhaZ7_type",
            "positive_end_qualifying_evidence_kinds": "experimental_positive_assay",
            "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
            "positive_end_challenge_evidence_kinds": "",
            "negative_end_challenge_evidence_kinds": "reported_label;experimental_positive_assay",
            "positive_end_anchor_status": "anchored_experimental_family_resolved",
            "negative_end_anchor_status": "no_formal_negative_record_available",
            "axis_boundary_note": "the native-active end has no formal negative record in the frozen evidence",
            "source_doi_or_pmid": "10.1186/1471-2105-10-89 (PMID 19296857)",
        },
        "axis_3_localization": {
            "axis_id": "axis_3_localization",
            "axis_name": "localization and architecture (extracellular signal peptide plus SBD versus intracellular lid)",
            "positive_end_state": "extracellular_signal_peptide_and_SBD",
            "negative_end_state": "intracellular_no_signal_peptide_lid_detected",
            "positive_end_qualifying_evidence_kinds": "experimental_positive_assay",
            "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
            "positive_end_challenge_evidence_kinds": "reported_label",
            "negative_end_challenge_evidence_kinds": "supported_pattern_state;reported_label",
            "positive_end_anchor_status": "anchored_experimental_external_records_only",
            "negative_end_anchor_status": "anchored_experimental_family_unresolved",
            "axis_boundary_note": "PF06850 is not an extracellular SBD marker in this reference set",
            "source_doi_or_pmid": "10.1186/1471-2105-10-89 (PMID 19296857)",
        },
        "axis_4_fold_neighborhood": {
            "axis_id": "axis_4_fold_neighborhood",
            "axis_name": "fold neighbourhood (PF10503 PHA depolymerase domain versus same-fold lipase/esterase with polar x1)",
            "positive_end_state": "PF10503_domain_present",
            "negative_end_state": "alpha_beta_hydrolase_lipase_esterase_polar_x1",
            "positive_end_qualifying_evidence_kinds": "experimental_positive_assay;supported_domain_state",
            "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
            "positive_end_challenge_evidence_kinds": "",
            "negative_end_challenge_evidence_kinds": (
                "supported_domain_state;candidate_polar_x1_pattern_state;reported_label"
            ),
            "positive_end_anchor_status": "anchored_domain_state_without_phenotype_evidence",
            "negative_end_anchor_status": "no_formal_negative_record_available",
            "axis_boundary_note": "polar x1 is not a biological negative; PF06850 is not fold evidence",
            "source_doi_or_pmid": "10.1186/1471-2105-10-89 (PMID 19296857)",
        },
    }


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    """Read a TSV into a list of plain string dictionaries (BOM tolerant)."""
    path = Path(path)
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"input is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _split_kinds(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in (value or "").split(";") if item.strip())


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _tsv_bool(value: bool) -> str:
    return "true" if value else "false"


def load_axis_definitions(path: str | Path) -> dict[str, dict[str, str]]:
    """Load and validate the declared discrimination axes."""
    rows = read_tsv(path)
    if not rows:
        raise ValueError("axis definitions are empty")
    definitions: dict[str, dict[str, str]] = {}
    for row in rows:
        missing = [field for field in AXIS_DEFINITION_FIELDS if field not in row]
        if missing:
            raise ValueError("axis definition missing columns: " + ",".join(missing))
        axis_id = (row.get("axis_id") or "").strip()
        if axis_id not in AXIS_IDS:
            raise ValueError(f"undeclared axis_id: {axis_id!r}")
        if axis_id in definitions:
            raise ValueError(f"duplicate axis_id: {axis_id}")
        for field in ("axis_name", "positive_end_state", "negative_end_state"):
            if not (row.get(field) or "").strip():
                raise ValueError(f"{axis_id}: {field} must not be empty")
        for field in (
            "positive_end_qualifying_evidence_kinds",
            "negative_end_qualifying_evidence_kinds",
            "source_doi_or_pmid",
        ):
            if not (row.get(field) or "").strip():
                raise ValueError(f"{axis_id}: {field} must not be empty")
        for field in (
            "positive_end_qualifying_evidence_kinds",
            "negative_end_qualifying_evidence_kinds",
            "positive_end_challenge_evidence_kinds",
            "negative_end_challenge_evidence_kinds",
        ):
            if field not in row:
                raise ValueError(f"{axis_id}: axis definition missing column {field}")
            kinds = _split_kinds(row[field])
            if field.endswith("qualifying_evidence_kinds") and not kinds:
                raise ValueError(f"{axis_id}: {field} must declare at least one evidence kind")
            unknown = [kind for kind in kinds if kind not in ANCHOR_EVIDENCE_KINDS]
            if unknown:
                raise ValueError(f"{axis_id}: unknown evidence kind(s): {','.join(unknown)}")
        definitions[axis_id] = {field: (row.get(field) or "").strip() for field in AXIS_DEFINITION_FIELDS}
    return definitions


def load_assignment_rules(path: str | Path) -> list[dict[str, str]]:
    """Load and validate the declared, priority-ordered axis assignment rules."""
    rows = read_tsv(path)
    if not rows:
        raise ValueError("assignment rules are empty")
    rules: list[dict[str, str]] = []
    seen: set[str] = set()
    priorities: set[int] = set()
    for row in rows:
        missing = [field for field in ASSIGNMENT_RULE_FIELDS if field not in row]
        if missing:
            raise ValueError("assignment rule missing columns: " + ",".join(missing))
        rule_id = (row.get("rule_id") or "").strip()
        if rule_id not in ASSIGNMENT_RULES:
            raise ValueError(f"undeclared assignment rule_id: {rule_id!r}")
        if rule_id in seen:
            raise ValueError(f"duplicate rule_id: {rule_id}")
        seen.add(rule_id)
        try:
            priority = int((row.get("priority") or "").strip())
        except ValueError as error:
            raise ValueError(f"{rule_id}: priority must be an integer") from error
        if priority in priorities:
            raise ValueError(f"{rule_id}: duplicate priority {priority}")
        priorities.add(priority)
        if not (row.get("qualifying_evidence_kind") or "").strip():
            raise ValueError(f"{rule_id}: qualifying_evidence_kind must not be empty")
        for field in ("axis_id", "axis_end", "source_layer", "rule_description"):
            if not (row.get(field) or "").strip():
                raise ValueError(f"{rule_id}: {field} must not be empty")
        rules.append({field: (row.get(field) or "").strip() for field in ASSIGNMENT_RULE_FIELDS})
    rules.sort(key=lambda item: int(item["priority"]))
    return rules


def prior_panel_evidence_status(panel_name: str) -> str:
    """Map a frozen acquisition panel name to its documented evidence status."""
    value = (panel_name or "").strip()
    if value in {"independent_experimental_positive", "mcl_pha_experimental_positive"}:
        return "experimental_positive"
    if value in {
        "mcl_pha_non_phb_negative",
        "intracellular_non_ephaz_negative",
        "fragment_or_incomplete_negative",
    }:
        return "experimental_negative"
    if value == "annotation_only_near_neighbor_negative":
        return "annotation_only"
    raise ValueError(f"unknown acquisition panel name: {panel_name!r}")


def normalize_declaration(row: dict[str, str]) -> dict[str, object]:
    """Normalize one declared external axis anchor."""
    missing = [field for field in DECLARATION_FIELDS if field not in row]
    if missing:
        raise ValueError("declaration missing columns: " + ",".join(missing))
    accession = (row.get("accession") or "").strip()
    if not accession:
        raise ValueError("declaration accession must not be empty")
    include_raw = (row.get("include") or "").strip().lower()
    if include_raw not in {"true", "false"}:
        raise ValueError(f"{accession}: include must be true or false")
    declaration: dict[str, object] = {
        "accession": accession,
        "include": include_raw == "true",
        "axis_id": (row.get("axis_id") or "").strip(),
        "axis_end": (row.get("axis_end") or "").strip(),
        "anchor_evidence_kind": (row.get("anchor_evidence_kind") or "").strip(),
        "evidence_status": (row.get("evidence_status") or "").strip(),
        "family_binding_status_expected": (row.get("family_binding_status_expected") or "").strip(),
        "anchor_justification": (row.get("anchor_justification") or "").strip(),
        "source_doi_or_pmid": (row.get("source_doi_or_pmid") or "").strip(),
        "evidence_source_path": (row.get("evidence_source_path") or "").strip(),
        "exclusion_reason": (row.get("exclusion_reason") or "").strip(),
    }
    return declaration


def check_declared_kind_for_end(
    axis_id: str,
    axis_end: str,
    kind: str,
    axis_definitions: dict[str, dict[str, str]],
) -> str:
    """Validate a declared/assigned evidence kind against an axis end and return the role.

    A kind is allowed at an end only when the axis definition declares it there, either
    as an anchoring kind or as a challenge-only kind.  Everything else fails closed, so
    an evidence kind can never be silently placed at an axis end that does not support it.
    """
    if axis_id not in axis_definitions:
        raise ValueError(f"undeclared axis_id: {axis_id!r}")
    if axis_end not in AXIS_ENDS:
        raise ValueError(f"axis_end must be one of {AXIS_ENDS}: {axis_end!r}")
    if kind not in ANCHOR_EVIDENCE_KINDS:
        raise ValueError(f"unknown anchor_evidence_kind: {kind!r}")
    axis = axis_definitions[axis_id]
    allowed = set(_split_kinds(axis[f"{axis_end}_end_qualifying_evidence_kinds"]))
    allowed.update(_split_kinds(axis[f"{axis_end}_end_challenge_evidence_kinds"]))
    if kind not in allowed:
        raise ValueError(
            f"{axis_id}:{axis_end}: evidence kind {kind!r} is not declared for this axis end "
            f"(declared: {','.join(sorted(allowed)) or 'none'})"
        )
    return classify_role(axis_end, kind, axis)


def classify_role(axis_end: str, evidence_kind: str, axis_definition: dict[str, str]) -> str:
    """Derive the role of a record: only declared anchoring kinds anchor an axis end."""
    if axis_end not in AXIS_ENDS:
        raise ValueError(f"axis_end must be one of {AXIS_ENDS}: {axis_end!r}")
    if evidence_kind not in ANCHOR_EVIDENCE_KINDS:
        raise ValueError(f"unknown anchor_evidence_kind: {evidence_kind!r}")
    qualifying = _split_kinds(axis_definition[f"{axis_end}_end_qualifying_evidence_kinds"])
    if evidence_kind in qualifying:
        return "positive_anchor" if axis_end == "positive" else "negative_anchor"
    return "challenge"


def is_formal_negative_eligible(row: dict[str, str]) -> bool:
    """A formal family-resolved negative needs formal negative evidence and a bound family."""
    if (row.get("role") or "").strip() != "negative_anchor":
        return False
    if (row.get("evidence_status") or "").strip() != "experimental_negative":
        return False
    if (row.get("family_binding_status") or "").strip() not in RESOLVED_BINDING_STATUSES:
        return False
    return bool((row.get("bound_reference_only_profile_id") or "").strip())


def build_ledger_index(ledger_rows: list[dict[str, str]]) -> dict[str, object]:
    """Index ledger accessions exactly and by versionless base accession."""
    exact: dict[str, str] = {}
    by_base: dict[str, set[str]] = defaultdict(set)
    for row in ledger_rows:
        accession = (row.get("accession") or "").strip()
        family = (row.get("phaded_family_id") or "").strip()
        if not accession or not family:
            continue
        exact.setdefault(accession, family)
        by_base[accession.split(".", 1)[0]].add(family)
    return {"exact": exact, "by_base": {key: sorted(value) for key, value in by_base.items()}}


def bind_family(accession: str, index: dict[str, object]) -> tuple[str, str]:
    """Bind an accession to an existing ledger family without ever upgrading a non-match."""
    accession = (accession or "").strip()
    exact = index["exact"]  # type: ignore[index]
    by_base = index["by_base"]  # type: ignore[index]
    if accession in exact:
        return "resolved_existing_family", exact[accession]
    families = by_base.get(accession.split(".", 1)[0], [])
    if len(families) == 1:
        return "resolved_existing_family_versionless", families[0]
    if len(families) > 1:
        return "ambiguous_existing_family", ""
    return "unresolved", ""


def validate_panel_invariants(rows: list[dict[str, str]]) -> None:
    """Fail closed when any axis-binding, evidence-status, or gate invariant is violated."""
    if not rows:
        raise ValueError("challenge panel is empty")
    seen: set[str] = set()
    definitions = builtin_axis_definitions()
    for row in rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("panel row has an empty accession")
        if accession in seen:
            raise ValueError(f"duplicate panel accession: {accession}")
        seen.add(accession)
        axis_id = (row.get("axis_id") or "").strip()
        axis_end = (row.get("axis_end") or "").strip()
        role = (row.get("role") or "").strip()
        status = (row.get("evidence_status") or "").strip()
        kind = (row.get("anchor_evidence_kind") or "").strip()
        origin = (row.get("record_origin") or "").strip()
        binding = (row.get("family_binding_status") or "").strip()
        if role not in ROLES:
            raise ValueError(f"{accession}: unsupported role {role!r}")
        if status not in EVIDENCE_STATUSES:
            raise ValueError(f"{accession}: unsupported evidence_status {status!r}")
        if origin not in RECORD_ORIGINS:
            raise ValueError(f"{accession}: unsupported record_origin {origin!r}")
        if binding not in BINDING_STATUSES:
            raise ValueError(f"{accession}: unsupported family_binding_status {binding!r}")
        if axis_id not in AXIS_IDS:
            raise ValueError(f"{accession}: record must bind exactly one declared axis, got {axis_id!r}")
        check_declared_kind_for_end(axis_id, axis_end, kind, definitions)
        expected_role = classify_role(axis_end, kind, definitions[axis_id])
        if role != expected_role:
            raise ValueError(f"{accession}: role {role!r} contradicts axis end/evidence kind ({expected_role})")
        for field in ("anchor_justification", "source_doi_or_pmid", "anchor_evidence_basis", "axis_name"):
            if not (row.get(field) or "").strip():
                raise ValueError(f"{accession}: {field} must not be empty")
        basis = (row.get("anchor_evidence_basis") or "").lower()
        if PF06850.lower() in basis:
            raise ValueError(
                f"{accession}: PF06850 must not be used as an axis anchor (DEFECT_D: it marks the "
                "intracellular Cys-type superfamily, not the extracellular SBD)"
            )
        if origin == "candidate_universe" and kind != "candidate_polar_x1_pattern_state":
            raise ValueError(
                f"{accession}: candidate-universe rows may only carry the polar x1 pattern kind "
                "(this stage does not consume a candidate-level motif completion layer, so no other "
                "evidence kind is available for candidate-universe accessions)"
            )
        if status == "annotation_only" and role == "negative_anchor":
            raise ValueError(f"{accession}: annotation_only must never be written as a formal negative")
        if status == "annotation_only" and _as_bool(row.get("formal_negative_eligible")):
            raise ValueError(f"{accession}: annotation_only must never be gate-eligible as a negative")
        eligible = is_formal_negative_eligible(row)
        if _as_bool(row.get("formal_negative_eligible")) != eligible:
            raise ValueError(f"{accession}: formal_negative_eligible contradicts evidence and binding")
        expected_slot = "family_resolved_negative" if eligible else "none"
        if (row.get("gate_credit_slot") or "").strip() != expected_slot:
            raise ValueError(f"{accession}: gate_credit_slot must be {expected_slot!r}")
        if binding not in RESOLVED_BINDING_STATUSES:
            if (row.get("bound_family") or "").strip() or (row.get("bound_reference_only_profile_id") or "").strip():
                raise ValueError(f"{accession}: unbound record must stay unresolved without a bound family")
        elif not (row.get("bound_family") or "").strip():
            raise ValueError(f"{accession}: resolved binding requires a bound family")
        if (row.get("bound_reference_only_profile_id") or "").strip() and not eligible and role == "negative_anchor":
            raise ValueError(f"{accession}: a bound reference-only negative must be gate eligible")


def summarize_axis_reachability(
    rows: list[dict[str, str]],
    axis_definitions: dict[str, dict[str, str]],
) -> dict[str, dict[str, object]]:
    """Per-axis counts of anchored controls and the remaining gap to one formal negative."""
    summary: dict[str, dict[str, object]] = {}
    for axis_id in AXIS_IDS:
        axis = axis_definitions[axis_id]
        axis_rows = [row for row in rows if row.get("axis_id") == axis_id]
        negative_rows = [row for row in axis_rows if row.get("axis_end") == "negative"]
        negative_anchors = [row for row in negative_rows if row.get("role") == "negative_anchor"]
        summary[axis_id] = {
            "axis_name": axis["axis_name"],
            "positive_end_state": axis["positive_end_state"],
            "negative_end_state": axis["negative_end_state"],
            "positive_end_anchor_status": axis["positive_end_anchor_status"],
            "negative_end_anchor_status": axis["negative_end_anchor_status"],
            "records_total": len(axis_rows),
            "positive_anchors": sum(row.get("role") == "positive_anchor" for row in axis_rows),
            "challenges": sum(row.get("role") == "challenge" for row in axis_rows),
            "anchored_negative_records": len(negative_rows),
            "formal_experimental_negatives": len(negative_anchors),
            "family_resolved_negatives": sum(
                _as_bool(row.get("formal_negative_eligible")) for row in negative_anchors
            ),
            "gate_credit_negatives": sum(
                (row.get("gate_credit_slot") or "") == "family_resolved_negative" for row in negative_anchors
            ),
            "negative_records_with_unresolved_family": sum(
                (row.get("family_binding_status") or "") not in RESOLVED_BINDING_STATUSES for row in negative_rows
            ),
            "negative_records_annotation_only": sum(
                row.get("evidence_status") == "annotation_only" for row in negative_rows
            ),
            "gap_to_one_family_resolved_negative": max(
                0, GATE["minimum_family_resolved_negative_count"]
                - sum(_as_bool(row.get("formal_negative_eligible")) for row in negative_anchors)
            ),
        }
    return summary


def summarize_profile_gaps(
    readiness_rows: list[dict[str, str]],
    hmm_built_by_profile: dict[str, bool] | None = None,
) -> list[dict[str, object]]:
    """Per-profile gate slot status and the remaining gap for every reference-only profile."""
    hmm_built_by_profile = hmm_built_by_profile or {}
    profiles: list[dict[str, object]] = []
    for row in sorted(readiness_rows, key=lambda item: (item.get("profile_id") or "").strip()):
        profile_id = (row.get("profile_id") or "").strip()
        if not profile_id:
            raise ValueError("readiness row without profile_id")

        def _int(field: str) -> int:
            raw = (row.get(field) or "0").strip() or "0"
            try:
                return int(raw)
            except ValueError as error:
                raise ValueError(f"{profile_id}: {field} is not an integer") from error

        positive = _int("existing_experimental_positive") + _int("external_bound_positive")
        heldout = _int("heldout_positive")
        negative = _int("external_bound_negative")
        challenge = _int("external_bound_challenge")
        hmm_built = bool(hmm_built_by_profile.get(profile_id, False))
        blocking: list[str] = []
        positive_gap = max(0, GATE["minimum_positive_count"] - positive)
        heldout_gap = max(0, GATE["minimum_heldout_positive_count"] - heldout)
        negative_gap = max(0, GATE["minimum_family_resolved_negative_count"] - negative)
        challenge_gap = max(0, GATE["minimum_challenge_count"] - challenge)
        if positive_gap:
            blocking.append("positive_gap")
        if heldout_gap:
            blocking.append("heldout_gap")
        if negative_gap:
            blocking.append("negative_gap")
        if challenge_gap:
            blocking.append("challenge_gap")
        blocking.append("zero_unexplained_hits_not_assessable")
        if not hmm_built:
            blocking.append("profile_hmm_not_built")
        gate_passed = not any(
            (positive_gap, heldout_gap, negative_gap, challenge_gap)
        ) and hmm_built
        profiles.append(
            {
                "profile_id": profile_id,
                "profile_kind": (row.get("profile_kind") or "").strip(),
                "phaded_superfamily": (row.get("phaded_superfamily") or "").strip(),
                "phaded_family_id": (row.get("phaded_family_id") or "").strip(),
                "existing_experimental_positive": _int("existing_experimental_positive"),
                "external_bound_positive": _int("external_bound_positive"),
                "current_positive_count_for_gate": positive,
                "required_positive_count": GATE["minimum_positive_count"],
                "positive_gap": positive_gap,
                "current_heldout_positive": heldout,
                "required_heldout_positive_count": GATE["minimum_heldout_positive_count"],
                "heldout_gap": heldout_gap,
                "current_family_resolved_negative": negative,
                "required_family_resolved_negative_count": GATE["minimum_family_resolved_negative_count"],
                "negative_gap": negative_gap,
                "current_challenge": challenge,
                "required_challenge_count": GATE["minimum_challenge_count"],
                "challenge_gap": challenge_gap,
                "axis_anchored_negatives_available_to_this_stage": 0,
                "gate_credit_negatives_from_this_stage": 0,
                "zero_unexplained_hits_status": "not_assessable_no_profile_hmm",
                "profile_hmms_built": hmm_built,
                "gate_passed": gate_passed,
                "blocking_reasons": blocking,
            }
        )
    return profiles


def build_reachability_document(
    panel_rows: list[dict[str, str]],
    axis_definitions: dict[str, dict[str, str]],
    profiles: list[dict[str, object]],
    external_resolutions: list[dict[str, str]],
    unanchored_records: list[dict[str, str]],
    *,
    run_id: str,
) -> dict[str, object]:
    """Assemble the machine-readable per-axis reachability and per-profile gap report."""
    per_axis = summarize_axis_reachability(panel_rows, axis_definitions)
    role_counts = Counter(row["role"] for row in panel_rows)
    status_counts = Counter(row["evidence_status"] for row in panel_rows)
    binding_counts = Counter(row["family_binding_status"] for row in panel_rows)
    totals = {
        "panel_records": len(panel_rows),
        "positive_anchor_records": role_counts.get("positive_anchor", 0),
        "negative_anchor_records": role_counts.get("negative_anchor", 0),
        "challenge_records": role_counts.get("challenge", 0),
        "evidence_status_counts": dict(sorted(status_counts.items())),
        "family_binding_status_counts": dict(sorted(binding_counts.items())),
        "formal_negative_eligible_records": sum(
            _as_bool(row.get("formal_negative_eligible")) for row in panel_rows
        ),
        "gate_credit_family_resolved_negatives": sum(
            (row.get("gate_credit_slot") or "") == "family_resolved_negative" for row in panel_rows
        ),
        "unanchored_records": len(unanchored_records),
        "external_declared_records": len(external_resolutions),
    }
    reference_only_profiles = len(profiles)
    summary = {
        "reference_only_profiles": reference_only_profiles,
        "profiles_meeting_full_gate": sum(bool(profile["gate_passed"]) for profile in profiles),
        "profiles_blocked": sum(not bool(profile["gate_passed"]) for profile in profiles),
        "profiles_with_negative_gap": sum(profile["negative_gap"] > 0 for profile in profiles),
        "profiles_with_heldout_gap": sum(profile["heldout_gap"] > 0 for profile in profiles),
        "profiles_with_positive_gap": sum(profile["positive_gap"] > 0 for profile in profiles),
        "profiles_with_challenge_gap": sum(profile["challenge_gap"] > 0 for profile in profiles),
        "total_positive_gap": sum(int(profile["positive_gap"]) for profile in profiles),
        "total_heldout_gap": sum(int(profile["heldout_gap"]) for profile in profiles),
        "total_family_resolved_negative_gap": sum(int(profile["negative_gap"]) for profile in profiles),
        "total_challenge_gap": sum(int(profile["challenge_gap"]) for profile in profiles),
        "profiles_with_built_hmm": sum(bool(profile["profile_hmms_built"]) for profile in profiles),
    }
    blockers = [
        "No record in the frozen evidence layers carries formal experimental-negative evidence "
        "that is also bound to a reference-only family, so the family-resolved negative slot is 0 for all axes.",
        "The axis-anchored negative controls that do exist are either annotation_only, bound to "
        "trained families (DED_hfam_4 / DED_hfam_8), or unresolved, so they cannot be credited.",
        "The 37 reference-only profiles have no built HMM (hmm_sha256 pending, no hmm_path), so the "
        "panel cannot be scored against them and the zero-unexplained-hit condition cannot be evaluated.",
        "Cross-family axis controls are not counted by the current per-profile gate counters "
        "(reconcile_phaded_reference_panel.py L135-147, acquire_phaded_reference_panel_amendment.py L184-191), "
        "which only credit records bound to the profile's own family.",
    ]
    return {
        "schema_version": 1,
        "run_id": run_id,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "status": "completed_candidate_only",
        "boundary": BOUNDARY,
        "gate": GATE,
        "axis_definitions": [axis_definitions[axis_id] for axis_id in AXIS_IDS],
        "per_axis": per_axis,
        "totals": totals,
        "summary": summary,
        "per_profile": profiles,
        "external_records_resolution": external_resolutions,
        "unanchored_records": unanchored_records,
        "blockers": blockers,
    }


def write_panel_tsv(path: str | Path, rows: list[dict[str, str]]) -> None:
    """Write the panel TSV, refusing to overwrite an existing result file."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing panel output: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(PANEL_FIELDS), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in PANEL_FIELDS})


def write_reachability_json(path: str | Path, document: dict[str, object]) -> None:
    """Write the reachability JSON, refusing to overwrite an existing result file."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing reachability output: {path}")
    path.write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _pfams(row: dict[str, str]) -> set[str]:
    return {item.strip() for item in (row.get("pfam_accessions") or "").split(";") if item.strip()}


def assign_ledger_row(
    ledger_row: dict[str, str],
    manifest_row: dict[str, str] | None,
    axis_definitions: dict[str, dict[str, str]],
) -> dict[str, str] | None:
    """Bind one reference-ledger row to exactly one axis end using the declared priority order."""
    manifest_row = manifest_row or {}
    lid_state = (manifest_row.get("lid_state") or "").strip()
    pfams = _pfams(manifest_row)
    substrate_class = (ledger_row.get("substrate_class") or "").strip()
    superfamily = (ledger_row.get("phaded_superfamily") or "").strip()
    substrate_detail = (ledger_row.get("substrate_detail") or "").strip()
    evidence_status = (ledger_row.get("evidence_status") or "").strip()

    axis_id: str | None = None
    axis_end = ""
    kind = ""
    basis = ""
    if lid_state in LID_POSITIVE_STATES:
        axis_id, axis_end, kind = "axis_3_localization", "negative", "supported_pattern_state"
        basis = f"lid_loop_detected_{lid_state}"
    elif substrate_class == "MCL":
        axis_id, axis_end, kind = "axis_1_substrate_chain_length", "negative", "reported_label"
        basis = "ledger_substrate_class_MCL"
    elif superfamily == PHAZ7_LIKE_SUPERFAMILY:
        axis_id, axis_end, kind = "axis_2_particle_state", "negative", "reported_label"
        basis = "ledger_superfamily_native_active_PhaZ7_type"
    elif substrate_detail == "denatured PHA" and evidence_status == "experimental_positive":
        axis_id, axis_end, kind = "axis_2_particle_state", "positive", "experimental_positive_assay"
        basis = "ledger_denatured_PHA_experimental_positive"
    elif substrate_class == "SCL" and evidence_status == "experimental_positive":
        axis_id, axis_end, kind = "axis_1_substrate_chain_length", "positive", "experimental_positive_assay"
        basis = "ledger_SCL_experimental_positive"
    elif PF10503 in pfams:
        axis_id, axis_end, kind = "axis_4_fold_neighborhood", "positive", "supported_domain_state"
        basis = "pfam_PF10503_binding"
    elif pfams & set(ALPHA_BETA_HYDROLASE_PFAMS):
        axis_id, axis_end, kind = "axis_4_fold_neighborhood", "negative", "supported_domain_state"
        basis = "pfam_alpha_beta_hydrolase_without_PF10503"
    if axis_id is None:
        return None
    check_declared_kind_for_end(axis_id, axis_end, kind, axis_definitions)
    return {
        "axis_id": axis_id,
        "axis_end": axis_end,
        "anchor_evidence_kind": kind,
        "anchor_evidence_basis": basis,
    }


def assign_confounder_row(
    confounder_row: dict[str, str],
    evidence_row: dict[str, str] | None,
    axis_definitions: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Bind one candidate-universe polar-x1 confounder to the fold-neighbourhood negative end."""
    pfams = _pfams(evidence_row or {})
    if pfams & set(ALPHA_BETA_HYDROLASE_PFAMS):
        basis = "x1_non_hydrophobic_task3_and_alpha_beta_hydrolase_pfam"
    elif pfams:
        basis = "x1_non_hydrophobic_task3_only_no_fold_domain"
    else:
        basis = "x1_non_hydrophobic_task3_no_tested_pfam"
    axis_id, axis_end, kind = "axis_4_fold_neighborhood", "negative", "candidate_polar_x1_pattern_state"
    check_declared_kind_for_end(axis_id, axis_end, kind, axis_definitions)
    return {
        "axis_id": axis_id,
        "axis_end": axis_end,
        "anchor_evidence_kind": kind,
        "anchor_evidence_basis": basis,
    }


def secondary_axis_capability(
    ledger_row: dict[str, str],
    manifest_row: dict[str, str] | None,
    assigned_axis: str,
) -> str:
    """List other axis ends whose end-state this ledger row also satisfies (informational only)."""
    manifest_row = manifest_row or {}
    capabilities: list[str] = []
    lid_state = (manifest_row.get("lid_state") or "").strip()
    pfams = _pfams(manifest_row)
    if (ledger_row.get("substrate_class") or "").strip() == "MCL" and assigned_axis != "axis_1_substrate_chain_length":
        capabilities.append("axis_1_substrate_chain_length:negative")
    if (ledger_row.get("phaded_superfamily") or "").strip() == PHAZ7_LIKE_SUPERFAMILY and assigned_axis != "axis_2_particle_state":
        capabilities.append("axis_2_particle_state:negative")
    if lid_state in LID_POSITIVE_STATES and assigned_axis != "axis_3_localization":
        capabilities.append("axis_3_localization:negative")
    if PF10503 in pfams and assigned_axis != "axis_4_fold_neighborhood":
        capabilities.append("axis_4_fold_neighborhood:positive")
    if (ledger_row.get("reported_localization") or "").strip() == "extracellular" and assigned_axis != "axis_3_localization":
        capabilities.append("axis_3_localization:positive_label_only")
    return ";".join(capabilities)


def load_external_declarations(
    path: str | Path,
    frozen_status: dict[str, dict[str, str]] | None = None,
) -> list[dict[str, object]]:
    """Load the declared external axis anchors and validate them fail-closed."""
    declarations = [normalize_declaration(row) for row in read_tsv(path)]
    if not declarations:
        raise ValueError("external axis declarations are empty")
    validate_external_declarations(declarations, frozen_status=frozen_status)
    return declarations


def validate_external_declarations(
    declarations: list[dict[str, object]],
    frozen_status: dict[str, dict[str, str]] | None = None,
) -> None:
    """Validate declared external anchors against the axis rules and frozen acquisition evidence."""
    definitions = builtin_axis_definitions()
    seen: set[str] = set()
    for declaration in declarations:
        accession = str(declaration["accession"])
        if accession in seen:
            raise ValueError(f"duplicate declared accession: {accession}")
        seen.add(accession)
        if not declaration["include"]:
            if not str(declaration["exclusion_reason"]).strip():
                raise ValueError(f"{accession}: excluded declaration requires an exclusion_reason")
            continue
        check_declared_kind_for_end(
            str(declaration["axis_id"]),
            str(declaration["axis_end"]),
            str(declaration["anchor_evidence_kind"]),
            definitions,
        )
        if not str(declaration["anchor_justification"]).strip():
            raise ValueError(f"{accession}: included declaration requires an anchor_justification")
        if not str(declaration["source_doi_or_pmid"]).strip():
            raise ValueError(f"{accession}: included declaration requires source_doi_or_pmid")
        if not str(declaration["evidence_source_path"]).strip():
            raise ValueError(f"{accession}: included declaration requires evidence_source_path")
        status = str(declaration["evidence_status"])
        if status not in EVIDENCE_STATUSES:
            raise ValueError(f"{accession}: unsupported evidence_status {status!r}")
        binding = str(declaration["family_binding_status_expected"])
        if binding not in BINDING_STATUSES:
            raise ValueError(f"{accession}: unsupported family_binding_status_expected {binding!r}")
        if frozen_status and accession in frozen_status:
            frozen = frozen_status[accession]
            frozen_evidence = (frozen.get("evidence_status") or "").strip()
            if frozen_evidence and frozen_evidence != status:
                raise ValueError(
                    f"{accession}: declared evidence_status {status!r} contradicts frozen acquisition "
                    f"evidence {frozen_evidence!r}"
                )
            frozen_binding = (frozen.get("family_binding_status") or "").strip()
            if frozen_binding and frozen_binding != binding:
                raise ValueError(
                    f"{accession}: declared binding {binding!r} contradicts frozen acquisition "
                    f"binding {frozen_binding!r}"
                )


def build_external_panel_rows(
    declarations: list[dict[str, object]],
    axis_definitions: dict[str, dict[str, str]],
    ledger_index: dict[str, object],
    reference_only_profile_by_family: dict[str, str],
    prior_by_accession: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Turn included declarations into panel rows and resolved declarations into an audit list."""
    rows: list[dict[str, str]] = []
    resolutions: list[dict[str, str]] = []
    for declaration in declarations:
        accession = str(declaration["accession"])
        if not declaration["include"]:
            resolutions.append(
                {
                    "accession": accession,
                    "resolution": "excluded",
                    "exclusion_reason": str(declaration["exclusion_reason"]),
                    "evidence_source_path": str(declaration["evidence_source_path"]),
                }
            )
            continue
        axis_id = str(declaration["axis_id"])
        axis_end = str(declaration["axis_end"])
        kind = str(declaration["anchor_evidence_kind"])
        status = str(declaration["evidence_status"])
        binding, family = bind_family(accession, ledger_index)
        expected_binding = str(declaration["family_binding_status_expected"])
        if binding != expected_binding:
            raise ValueError(
                f"{accession}: computed family binding {binding!r} differs from declared {expected_binding!r}"
            )
        bound_profile = reference_only_profile_by_family.get(family, "") if family else ""
        role = classify_role(axis_end, kind, axis_definitions[axis_id])
        prior = prior_by_accession.get(accession, {})
        evidence_path = str(declaration["evidence_source_path"])
        if prior.get("panel"):
            evidence_path = f"{evidence_path} [panel={prior['panel']}]"
        row = {field: "" for field in PANEL_FIELDS}
        row.update(
            axis_id=axis_id,
            axis_name=axis_definitions[axis_id]["axis_name"],
            axis_end=axis_end,
            role=role,
            accession=accession,
            record_origin="external_acquisition",
            evidence_status=status,
            anchor_evidence_kind=kind,
            anchor_evidence_basis=f"external_declared:{prior.get('panel', 'acquisition_03_panel')}",
            family_binding_status=binding,
            bound_family=family,
            bound_reference_only_profile_id=bound_profile,
            anchor_justification=str(declaration["anchor_justification"]),
            source_doi_or_pmid=str(declaration["source_doi_or_pmid"]),
            evidence_source_path=evidence_path,
            secondary_axis_capability="",
        )
        row["formal_negative_eligible"] = _tsv_bool(is_formal_negative_eligible(row))
        row["gate_credit_slot"] = (
            "family_resolved_negative" if _as_bool(row["formal_negative_eligible"]) else "none"
        )
        rows.append(row)
        resolutions.append(
            {
                "accession": accession,
                "resolution": "included_as_axis_anchor",
                "exclusion_reason": "",
                "evidence_source_path": evidence_path,
            }
        )
    return rows, resolutions


def build_ledger_panel_rows(
    ledger_rows: list[dict[str, str]],
    manifest_by_accession: dict[str, dict[str, str]],
    axis_definitions: dict[str, dict[str, str]],
    reference_only_profile_by_family: dict[str, str],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Assign every reference-ledger row to at most one axis end, or record it as unanchored."""
    rows: list[dict[str, str]] = []
    unanchored: list[dict[str, str]] = []
    for ledger_row in ledger_rows:
        accession = (ledger_row.get("accession") or "").strip()
        manifest_row = manifest_by_accession.get(accession, {})
        assignment = assign_ledger_row(ledger_row, manifest_row, axis_definitions)
        if assignment is None:
            unanchored.append(
                {
                    "accession": accession,
                    "record_origin": "reference_ledger",
                    "phaded_superfamily": (ledger_row.get("phaded_superfamily") or "").strip(),
                    "phaded_family_id": (ledger_row.get("phaded_family_id") or "").strip(),
                    "evidence_status": (ledger_row.get("evidence_status") or "").strip(),
                    "reason": "no_supported_state_on_any_declared_axis_label_only_record",
                }
            )
            continue
        axis_id = assignment["axis_id"]
        axis_end = assignment["axis_end"]
        kind = assignment["anchor_evidence_kind"]
        family = (ledger_row.get("phaded_family_id") or "").strip()
        role = classify_role(axis_end, kind, axis_definitions[axis_id])
        row = {field: "" for field in PANEL_FIELDS}
        row.update(
            axis_id=axis_id,
            axis_name=axis_definitions[axis_id]["axis_name"],
            axis_end=axis_end,
            role=role,
            accession=accession,
            record_origin="reference_ledger",
            evidence_status=(ledger_row.get("evidence_status") or "").strip(),
            anchor_evidence_kind=kind,
            anchor_evidence_basis=assignment["anchor_evidence_basis"],
            family_binding_status="resolved_existing_family",
            bound_family=family,
            bound_reference_only_profile_id=reference_only_profile_by_family.get(family, ""),
            anchor_justification=_ledger_justification(ledger_row, manifest_row, assignment, axis_definitions),
            source_doi_or_pmid=_ledger_source(ledger_row),
            evidence_source_path=(
                "runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv; "
                "runs/20260917_phaded_reference_residue_mapping_01/results/reference_mapping_manifest.tsv"
            ),
            secondary_axis_capability=secondary_axis_capability(ledger_row, manifest_row, axis_id),
        )
        row["formal_negative_eligible"] = _tsv_bool(is_formal_negative_eligible(row))
        row["gate_credit_slot"] = (
            "family_resolved_negative" if _as_bool(row["formal_negative_eligible"]) else "none"
        )
        rows.append(row)
    return rows, unanchored


def _ledger_justification(
    ledger_row: dict[str, str],
    manifest_row: dict[str, str],
    assignment: dict[str, str],
    axis_definitions: dict[str, dict[str, str]],
) -> str:
    end_state = axis_definitions[assignment["axis_id"]][
        "positive_end_state" if assignment["axis_end"] == "positive" else "negative_end_state"
    ]
    return (
        f"{assignment['anchor_evidence_basis']} -> {assignment['axis_id']}:{assignment['axis_end']}"
        f" ({end_state}); frozen ledger labels: substrate_class={ledger_row.get('substrate_class', '')}, "
        f"substrate_detail={ledger_row.get('substrate_detail', '')}, "
        f"superfamily={ledger_row.get('phaded_superfamily', '')}, "
        f"evidence_status={ledger_row.get('evidence_status', '')}; "
        f"mapping states: lid_state={manifest_row.get('lid_state', '')}, "
        f"pfam_accessions={manifest_row.get('pfam_accessions', '')}. "
        "Label/domain/pattern evidence only; no PHB/PHA phenotype is claimed."
    )


def _ledger_source(ledger_row: dict[str, str]) -> str:
    doi = (ledger_row.get("primary_doi") or "").strip()
    pmid = (ledger_row.get("pmid") or "").strip()
    if not doi and not pmid:
        raise ValueError(f"{ledger_row.get('accession')}: ledger row without DOI/PMID")
    return f"{doi or 'pending'} (PMID {pmid or 'pending'})"


def build_confounder_panel_rows(
    confounder_rows: list[dict[str, str]],
    evidence_by_accession: dict[str, dict[str, str]],
    axis_definitions: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    """Bind the 563 Task-3 polar-x1 candidate confounders to the fold-neighbourhood negative end."""
    rows: list[dict[str, str]] = []
    for confounder in confounder_rows:
        accession = (confounder.get("accession") or "").strip()
        if not accession:
            raise ValueError("confounder row without accession")
        evidence_row = evidence_by_accession.get(accession, {})
        assignment = assign_confounder_row(confounder, evidence_row, axis_definitions)
        axis_id = assignment["axis_id"]
        axis_end = assignment["axis_end"]
        kind = assignment["anchor_evidence_kind"]
        role = classify_role(axis_end, kind, axis_definitions[axis_id])
        row = {field: "" for field in PANEL_FIELDS}
        row.update(
            axis_id=axis_id,
            axis_name=axis_definitions[axis_id]["axis_name"],
            axis_end=axis_end,
            role=role,
            accession=accession,
            record_origin="candidate_universe",
            evidence_status="annotation_only",
            anchor_evidence_kind=kind,
            anchor_evidence_basis=assignment["anchor_evidence_basis"],
            family_binding_status="unresolved",
            bound_family="",
            bound_reference_only_profile_id="",
            anchor_justification=(
                f"{assignment['anchor_evidence_basis']} -> {axis_id}:{axis_end}; "
                f"Task 3 non-hydrophobic x1={confounder.get('lipase_box_x1', '')} "
                f"(group={confounder.get('x1_group', '')}), confounder_flag="
                f"{confounder.get('confounder_flag', '')}; candidate-universe accession has no family-level "
                "binding evidence, so it stays unresolved and is a challenge control only "
                "(polar x1 is not a biological negative)."
            ),
            source_doi_or_pmid="10.1186/1471-2105-10-89 (PMID 19296857)",
            evidence_source_path=(
                "runs/20260917_phaded_x1_hydrophobicity_01/results/confounder_candidates.tsv; "
                "runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv"
            ),
            secondary_axis_capability="",
        )
        row["formal_negative_eligible"] = _tsv_bool(is_formal_negative_eligible(row))
        row["gate_credit_slot"] = (
            "family_resolved_negative" if _as_bool(row["formal_negative_eligible"]) else "none"
        )
        rows.append(row)
    return rows


def _reference_only_profiles(profile_manifest_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in profile_manifest_rows if (row.get("model_status") or "").strip() == "reference_only"]


def design(
    *,
    axis_definitions_path: str | Path,
    assignment_rules_path: str | Path,
    declarations_path: str | Path,
    ledger_path: str | Path,
    mapping_manifest_path: str | Path,
    confounder_path: str | Path,
    candidate_evidence_path: str | Path,
    readiness_path: str | Path,
    profile_manifest_path: str | Path,
    external_prior_path: str | Path,
    external_bound_path: str | Path,
    output_dir: str | Path,
    expected_reference_only_count: int = 37,
) -> dict[str, object]:
    """Build the axis-anchored challenge panel and the reachability report."""
    axis_definitions = load_axis_definitions(axis_definitions_path)
    if set(axis_definitions) != set(AXIS_IDS):
        raise ValueError(
            "axis definitions must declare exactly the four planned axes: "
            + ",".join(sorted(set(AXIS_IDS) - set(axis_definitions)))
        )
    rules = load_assignment_rules(assignment_rules_path)
    if [rule["rule_id"] for rule in rules] != list(ASSIGNMENT_RULES):
        raise ValueError("assignment rules must declare every implemented rule exactly once")

    ledger_rows = read_tsv(ledger_path)
    manifest_by_accession = {row["accession"]: row for row in read_tsv(mapping_manifest_path)}
    profile_manifest_rows = read_tsv(profile_manifest_path)
    reference_only = _reference_only_profiles(profile_manifest_rows)
    if len(reference_only) != expected_reference_only_count:
        raise ValueError(
            f"expected {expected_reference_only_count} reference-only profiles, observed {len(reference_only)}"
        )
    reference_only_profile_by_family = {
        (row.get("phaded_family_id") or "").strip(): (row.get("profile_id") or "").strip()
        for row in reference_only
        if (row.get("profile_kind") or "").strip() == "family" and (row.get("phaded_family_id") or "").strip()
    }
    hmm_built_by_profile = {
        (row.get("profile_id") or "").strip(): bool((row.get("hmm_sha256") or "").strip())
        for row in reference_only
    }

    ledger_index = build_ledger_index(ledger_rows)
    frozen_status = {
        row["accession"]: row
        for row in read_tsv(external_bound_path)
        if (row.get("accession") or "").strip()
    }
    declarations = load_external_declarations(declarations_path, frozen_status=frozen_status)
    prior_by_accession = {row["accession"]: row for row in read_tsv(external_prior_path) if row.get("accession")}

    ledger_panel_rows, unanchored = build_ledger_panel_rows(
        ledger_rows, manifest_by_accession, axis_definitions, reference_only_profile_by_family
    )
    confounder_rows = read_tsv(confounder_path)
    evidence_by_accession = {row["accession"]: row for row in read_tsv(candidate_evidence_path)}
    confounder_panel_rows = build_confounder_panel_rows(
        confounder_rows, evidence_by_accession, axis_definitions
    )
    external_panel_rows, external_resolutions = build_external_panel_rows(
        declarations, axis_definitions, ledger_index, reference_only_profile_by_family, prior_by_accession
    )

    panel_rows = [*ledger_panel_rows, *external_panel_rows, *confounder_panel_rows]
    validate_panel_invariants(panel_rows)

    readiness_rows = read_tsv(readiness_path)
    if len(readiness_rows) != expected_reference_only_count:
        raise ValueError(
            f"expected {expected_reference_only_count} readiness rows, observed {len(readiness_rows)}"
        )
    profiles = summarize_profile_gaps(readiness_rows, hmm_built_by_profile)
    document = build_reachability_document(
        panel_rows, axis_definitions, profiles, external_resolutions, unanchored, run_id=RUN_ID
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_panel_tsv(output_dir / "challenge_panel_by_axis.tsv", panel_rows)
    write_reachability_json(output_dir / "axis_reachability.json", document)
    return {"panel_rows": panel_rows, "report": document}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=Path("runs") / RUN_ID)
    parser.add_argument("--expected-reference-only-count", type=int, default=37)
    args = parser.parse_args()
    run_dir = args.run_dir
    result = design(
        axis_definitions_path=run_dir / "inputs" / "discrimination_axis_definitions.tsv",
        assignment_rules_path=run_dir / "inputs" / "discrimination_axis_assignment_rules.tsv",
        declarations_path=run_dir / "inputs" / "external_axis_anchor_declarations.tsv",
        ledger_path=Path("runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv"),
        mapping_manifest_path=Path(
            "runs/20260917_phaded_reference_residue_mapping_01/results/reference_mapping_manifest.tsv"
        ),
        confounder_path=Path("runs/20260917_phaded_x1_hydrophobicity_01/results/confounder_candidates.tsv"),
        candidate_evidence_path=Path(
            "runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv"
        ),
        readiness_path=Path(
            "runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/profile_calibration_readiness.tsv"
        ),
        profile_manifest_path=Path(
            "runs/20260915_phaded_reference_panel_acquisition_02/inputs/profile_manifest.tsv"
        ),
        external_prior_path=Path(
            "runs/20260915_phaded_reference_panel_acquisition_02/inputs/external_panel_evidence_prior.tsv"
        ),
        external_bound_path=Path(
            "runs/20260916_phaded_reference_panel_acquisition_03/results/panel_candidates_bound.tsv"
        ),
        output_dir=run_dir / "results",
        expected_reference_only_count=args.expected_reference_only_count,
    )
    print(
        json.dumps(
            {
                "panel_records": result["report"]["totals"]["panel_records"],
                "role_counts": {
                    "positive_anchor": result["report"]["totals"]["positive_anchor_records"],
                    "negative_anchor": result["report"]["totals"]["negative_anchor_records"],
                    "challenge": result["report"]["totals"]["challenge_records"],
                },
                "gate_credit_family_resolved_negatives": result["report"]["totals"][
                    "gate_credit_family_resolved_negatives"
                ],
                "profiles_blocked": result["report"]["summary"]["profiles_blocked"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
