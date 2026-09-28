#!/usr/bin/env python3
"""Compare the frozen v1 PhaDED high-confidence table with the v2 catalog, and
report the rule-sensitivity of the candidate composition under three
preregistered scenarios.

Boundary
--------
Candidate-only and read-only.  This module never reads project data by itself: the
caller passes a v1 label map, a v2 label map and a table of evidence rows, and the
CLI only reads the paths it is given.  It introduces no scoring, no family call,
no disposition change and no deletion; it only *reports* differences.  The frozen
v1 high-confidence table is a **comparison input only** (design doc
``docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md``,
"候选目录": the old high-confidence table is never a v2 authority).

Public API
----------
``compare(v1, v2)``            accession-level removed / added / relabelled sets.
``impact_summary(result)``     JSON-safe counts, label transitions and deltas.
``SCENARIOS`` / ``CRITERIA``   the preregistered scenario set, defined as data.
``resolve_scenario_names(...)``  validate a requested scenario subset.
``evaluate_scenario(...)``     evaluate one scenario over rows of mappings.
``analyze(...)`` / ``build_report(...)``  the full JSON report + TSV rows.
``composition_blocks(...)`` / ``validate_report_guard(...)``  the 76:24 guard.
``main(argv=None) -> int``     CLI: JSON report + accession-level TSV.

The v1 rules encoded by ``current_v1_rules``
--------------------------------------------
Source (READ-ONLY, frozen): ``pipeline/scripts/filter_phaded_high_confidence.py``
(module header: "DEPRECATED / FROZEN v1 SEMANTICS -- COMPARISON INPUT ONLY") and
its preregistered decision table ``inputs/filter_criteria.md`` of run
``runs/20260919_phaded_high_confidence_filter_01``.  The scenario reproduces, in
the frozen evaluation order:

1. **superfamily scope** -- only the eight superfamilies in
   ``filter_phaded_high_confidence.SUPERFAMILY_CRITERIA`` are ever counted; any
   other row is out of v1 scope (the frozen loop ``continue``s, i.e. it is not
   even in the v1 denominator).  Reported as ``not_in_v1_scope``.
2. **common criteria** (frozen ``ok_common``, hold reason ``common_criteria``) --
   ``interpro_status == "interpro_supported"``;
   ``architecture_consistency != "conflicting"``;
   a strong profile, i.e. ``profile_evidence_status`` in
   {``profile_trained_hit``, ``profile_ambiguous_family``} **or**
   ``profile_evidence_status == "discovery_hit"`` with
   ``best_evalue < 1e-10``;
   and the accession is not in the frozen confounder set.
3. **localization hold rule** (hold reason ``localization_conflict``) -- for
   *export-requiring* superfamilies ``signalp_class`` must be in
   {SP, LIPO, TAT, TATLIPO}; for *no-export* superfamilies it must be in
   {OTHER, pending}.  This is the v1 behaviour that held **48,038 of 67,720
   type-1 candidates (71%)** on the missing-secretion-signal rule alone
   (``filter_summary.json`` -> ``per_superfamily["extracellular dPHASCL type 1"]
   .fail_localization = 48038``; see
   ``docs/T141_20260919_phaded_high_confidence_filter_status.md`` and
   ``docs/T141_20260921_phaded_comprehensive_review.md`` section 7.3).
4. **per-superfamily architecture criteria** (frozen ``arch_criteria``;
   hold reason ``architecture_criteria``, or ``nucleophile_conflict`` when the
   family-level nucleophile typing is a conflict) -- evaluated with the frozen
   ``evaluate`` helper, including ``x1_hydrophobic=1`` reading
   ``lipase_box_x1`` in the hydrophobic set and ``pf06850=detected`` reading
   ``sbd_pf06850_binding_state``.

Fail-closed: a column that is present but blank never passes a criterion, exactly
as in the frozen module.  A column that is **absent** from the row mapping cannot
be reconstructed, so that row's count is reported as ``pending`` (an explicit
``pending_proteins`` / ``pending_criteria`` count) instead of being guessed
either way.  The one v1 rule that lives outside the row -- confounder membership,
which comes from a frozen companion file -- is ``pending`` unless the caller
supplies ``--confounders`` (or the row carries ``is_confounder``).

Reconstruction assumptions (documented, not silent)
---------------------------------------------------
* v1 read ``catalytic_domain_type`` from the domain-type table and the other
  architecture fields from the motif table; the flat row here carries one
  ``catalytic_domain_type`` column, which is the frozen output column of the same
  name.
* v1 drew the superfamily either from the trained-profile matrix
  (``phaded_superfamily_best``) or from the discovery classification
  (``superfamily_claim`` with ``superfamily_confidence == "unique"``); here the
  row's ``historical_superfamily`` (then ``superfamily``) is that claim.  When
  ``superfamily_confidence`` is present it must be ``unique`` for the discovery
  branch.
* a row may carry either ``signalp_class`` (v1 column) or the v2
  ``transport_signal_prediction``, which uses the same SignalP classes; the first
  existing column wins.
* ``SUPERFAMILY_CRITERIA``, ``EXPORT``, ``NO_EXPORT``, ``HYDROPHOBIC`` and
  ``nucleophile_type`` are re-encoded here and asserted equal to the frozen
  source by the test module, so the two cannot drift.

The 76:24 guard
---------------
The intracellular:extracellular split of the v1 candidate set (76.1% : 23.9%) is
a *candidate composition under a named rule*, never a biological ratio: the
criteria are asymmetric (type-1 needs secretion signal + geometry + architecture)
and the underlying literature is biased toward experimentally tractable
extracellular enzymes
(``docs/T141_20260921_phaded_comprehensive_review.md`` section 7.3).  Every
composition block therefore carries ``scenario`` + ``rule`` + the fixed
``COMPOSITION_CAVEAT`` string, and ``validate_report_guard`` refuses to write a
report whose keys or labels could present the split as a biological finding
(forbidden underscore-joined fragments such as ``biological_ratio`` /
``true_ratio`` / ``in_vivo_ratio``, and any composition block outside the
canonical ``scenarios[*].compositions[*]`` location).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Mapping, Sequence

__all__ = [
    # comparison
    "compare",
    "impact_summary",
    # preregistered scenarios
    "SCENARIOS",
    "CRITERIA",
    "Scenario",
    "ScenarioContext",
    "ScenarioResult",
    "resolve_scenario_names",
    "evaluate_row",
    "evaluate_scenario",
    # report assembly
    "Analysis",
    "analyze",
    "build_report",
    "composition_blocks",
    "validate_report_guard",
    # row accessors
    "row_superfamily",
    "row_family",
    "row_transport_signal",
    "candidate_localization_class",
    # frozen v1 rules re-encoded for the comparison
    "V1_SUPERFAMILY_CRITERIA",
    "V1_EXPORT_CLASSES",
    "V1_NO_EXPORT_CLASSES",
    "V1_HYDROPHOBIC",
    "V1_EXPORT_REQUIRING",
    "V1_NO_EXPORT_REQUIRED",
    "V1_RULE_CRITERIA",
    "V1_HOLD_REASONS",
    "V1_HOLD_REASON_BY_CRITERION",
    "V1_SOURCE",
    "v1_nucleophile_type",
    # composition guard
    "COMPOSITION_BASIS",
    "COMPOSITION_CAVEAT",
    "FORBIDDEN_LABEL_FRAGMENTS",
    # table IO and output guards
    "read_rows",
    "read_label_map",
    "read_accession_set",
    "bound_file",
    "require_file",
    "prepare_output_dir",
    "output_path",
    "write_accession_tsv",
    "write_outputs",
    # CLI
    "build_parser",
    "run",
    "main",
    "REPORT_JSON",
    "ACCESSIONS_TSV",
    "OUTPUT_FILENAMES",
    "ACCESSION_TSV_COLUMNS",
    "PENDING",
    "SCHEMA_VERSION",
]

SCHEMA_VERSION = "1.0"
GENERATOR = "compare_phaded_catalog_versions.py"

PENDING = "pending"
"""Verdict for a rule that cannot be evaluated from the supplied columns."""

REPORT_JSON = "catalog_version_comparison.json"
ACCESSIONS_TSV = "accession_differences.tsv"
OUTPUT_FILENAMES = (REPORT_JSON, ACCESSIONS_TSV)

ACCESSION_TSV_COLUMNS = (
    "section",
    "accession",
    "change_kind",
    "scenario",
    "v1_label",
    "v2_label",
)

TRUTHY_VALUES = frozenset({"1", "true", "yes", "y", "confounder"})
FALSY_VALUES = frozenset({"0", "false", "no", "n", ""})

SEQUENCE_FAMILY_VALIDATED = "sequence_family_hmm_validated"
HIGH_CONFIDENCE_LABEL = "high_confidence"

CLAIM_BOUNDARY = (
    "Candidate-only comparison report. v2 dispositions describe sequence and "
    "evidence state, not PHB/PHA degradation phenotypes; the frozen v1 "
    "high-confidence table is a comparison input only. Every count here is a "
    "count of candidates under a named rule, not a measured biological "
    "proportion."
)

COMPOSITION_BASIS = "candidate_composition_under_named_rule"

COMPOSITION_CAVEAT = (
    "Candidate composition under the named rule only. This is not a biological "
    "ratio and must not be reported as one: the intracellular:extracellular "
    "split is produced by asymmetric criteria (type-1 requires a predicted "
    "secretion signal plus geometry plus architecture, so 48,038 of 67,720 "
    "type-1 candidates = 71% are held on the secretion-signal rule alone) and by "
    "research bias toward experimentally tractable extracellular enzymes. Never "
    "cite it as the true in vivo, in-genome or environmental proportion "
    "(docs/T141_20260921_phaded_comprehensive_review.md section 7.3)."
)

FORBIDDEN_LABEL_FRAGMENTS = (
    "biological_ratio",
    "true_ratio",
    "in_vivo_ratio",
    "invivo_ratio",
    "real_ratio",
    "genomic_ratio",
    "environmental_ratio",
    "absolute_ratio",
)


# --------------------------------------------------------------------------- #
# frozen v1 rule constants (re-encoded; parity asserted by the test module)
# --------------------------------------------------------------------------- #

V1_SOURCE = (
    "pipeline/scripts/filter_phaded_high_confidence.py (frozen v1 filter, "
    "'FROZEN v1 SEMANTICS -- COMPARISON INPUT ONLY'); preregistered decision "
    "table inputs/filter_criteria.md of runs/20260919_phaded_high_confidence_filter_01"
)

V1_EXPORT_CLASSES = frozenset({"SP", "LIPO", "TAT", "TATLIPO"})
V1_NO_EXPORT_CLASSES = frozenset({"OTHER", "pending"})
V1_HYDROPHOBIC = frozenset({"A", "C", "F", "I", "L", "M", "V", "W", "Y"})
V1_STRONG_PROFILE_STATUSES = frozenset({"profile_trained_hit", "profile_ambiguous_family"})
V1_DISCOVERY_STATUS = "discovery_hit"
V1_DISCOVERY_EVALUE_THRESHOLD = 1e-10
V1_UNIQUE_CONFIDENCE = "unique"

V1_SUPERFAMILY_CRITERIA: dict[str, tuple[str, list[str]]] = {
    "extracellular dPHASCL type 1": (
        "export",
        ["lipase_box_state=supported", "x1_hydrophobic=1", "catalytic_domain_type=type1_verified"],
    ),
    "extracellular dPHASCL type 2": (
        "export",
        ["lipase_box_state=supported", "x1_hydrophobic=1", "catalytic_domain_type=type2_verified"],
    ),
    "intracellular nPHASCL without lipase box": (
        "no_export",
        ["lipase_box_state=not_detected_pattern", "pf06850=detected"],
    ),
    "intracellular nPHASCL with lipase box": (
        "no_export",
        ["lipase_box_state=supported", "x1_hydrophobic=1"],
    ),
    "intracellular nPHAMCL": (
        "no_export",
        ["lipase_box_state=supported", "x1_hydrophobic=1", "lid_state=supported"],
    ),
    "extracellular dPHAMCL": (
        "export",
        ["lipase_box_state=supported", "x1_hydrophobic=1"],
    ),
    "periplasmic PHA depolymerases": (
        "export",
        ["catalytic_domain_type=type2_verified"],
    ),
    "extracellular native-SCL/PhaZ7-like": (
        "export",
        ["ahsmg_state=supported"],
    ),
}

V1_EXPORT_REQUIRING = frozenset(
    name for name, (requirement, _) in V1_SUPERFAMILY_CRITERIA.items() if requirement == "export"
)
V1_NO_EXPORT_REQUIRED = frozenset(
    name for name, (requirement, _) in V1_SUPERFAMILY_CRITERIA.items() if requirement == "no_export"
)

# frozen nucleophile-type labels (v1 defect kept verbatim: AHSMG as a third type)
V1_CYS_SUPERFAMILY = "intracellular nPHASCL without lipase box"
V1_PHAZ7_SUPERFAMILY = "extracellular native-SCL/PhaZ7-like"
V1_NT_SER = "ser"
V1_NT_CYS = "cys"
V1_NT_AHSMG = "ahsmg"
V1_NT_CONFLICT = "conflict"
V1_NT_UNDETERMINED = "undetermined"

V1_HOLD_REASONS = (
    "not_in_v1_scope",
    "common_criteria",
    "localization_conflict",
    "nucleophile_conflict",
    "architecture_criteria",
)
"""Frozen v1 hold/defer reason labels (``hold_reason`` in the frozen v1 output)."""

V1_HOLD_REASON_BY_CRITERION = {
    "v1_superfamily_recognized": "not_in_v1_scope",
    "v1_interpro_supported": "common_criteria",
    "v1_no_architecture_conflict": "common_criteria",
    "v1_strong_profile": "common_criteria",
    "v1_not_confounder": "common_criteria",
    "v1_localization_export_requirement": "localization_conflict",
    "v1_localization_no_export_requirement": "localization_conflict",
    "v1_architecture_criteria": "architecture_criteria",
}

EXCLUDED_REASON_BY_CRITERION = dict(V1_HOLD_REASON_BY_CRITERION)
EXCLUDED_REASON_BY_CRITERION["sequence_family_hmm_validated"] = (
    "model_layer_not_sequence_validated"
)

SUPERFAMILY_COLUMNS = ("historical_superfamily", "superfamily")
FAMILY_COLUMNS = (
    "sequence_family_call",
    "nucleophile_family",
    "phaded_family_best",
    "phaded_superfamily_best",
)
TRANSPORT_COLUMNS = ("transport_signal_prediction", "signalp_class")
ACCESSION_COLUMNS = ("accession", "protein_id")
LOCALIZATION_CLASSES = (
    ("intracellular", "intracellular"),
    ("extracellular", "extracellular"),
    ("periplasmic", "periplasmic"),
)


def v1_nucleophile_type(superfamily: str, lipase_box_state: str, ahsmg_state: str) -> str:
    """Frozen v1 family-level nucleophile typing (verbatim re-encoding)."""
    if superfamily == V1_CYS_SUPERFAMILY:
        return V1_NT_CONFLICT if lipase_box_state == "supported" else V1_NT_CYS
    if superfamily == V1_PHAZ7_SUPERFAMILY:
        return V1_NT_AHSMG if ahsmg_state == "supported" else V1_NT_CONFLICT
    if lipase_box_state == "supported":
        return V1_NT_SER
    if lipase_box_state == "not_detected_pattern":
        return V1_NT_CONFLICT
    return V1_NT_UNDETERMINED


# --------------------------------------------------------------------------- #
# accession-level comparison (v1 vs v2)
# --------------------------------------------------------------------------- #


def compare(
    v1: Mapping[str, str] | None = None,
    v2: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Compare two accession -> label maps.

    Returns a dict whose ``removed`` / ``added`` / ``relabelled`` values are
    **sets**, so the answer is independent of input insertion order:

    * ``removed``    -- in v1 only;
    * ``added``      -- in v2 only;
    * ``relabelled`` -- in both with a different label;
    * ``retained_unchanged`` -- in both with the same label.

    The input label maps are echoed back under ``v1_labels`` / ``v2_labels`` so
    that label transitions can be summarised without re-reading the tables.
    """
    first = dict(v1 or {})
    second = dict(v2 or {})
    first_accessions = set(first)
    second_accessions = set(second)
    shared = first_accessions & second_accessions
    return {
        "removed": first_accessions - second_accessions,
        "added": second_accessions - first_accessions,
        "relabelled": {
            accession for accession in shared if first[accession] != second[accession]
        },
        "retained_unchanged": {
            accession for accession in shared if first[accession] == second[accession]
        },
        "v1_labels": first,
        "v2_labels": second,
    }


def impact_summary(result: Mapping[str, Any]) -> dict[str, Any]:
    """JSON-safe summary of a :func:`compare` result (counts, transitions, deltas)."""
    v1_labels: Mapping[str, str] = result.get("v1_labels", {})
    v2_labels: Mapping[str, str] = result.get("v2_labels", {})
    v1_counts = Counter(v1_labels.values())
    v2_counts = Counter(v2_labels.values())
    transitions = Counter(
        f"{v1_labels[accession]}->{v2_labels[accession]}"
        for accession in sorted(result.get("relabelled", ()))
    )
    labels = sorted(set(v1_counts) | set(v2_counts))
    return {
        "v1_count": len(v1_labels),
        "v2_count": len(v2_labels),
        "removed_count": len(result.get("removed", ())),
        "added_count": len(result.get("added", ())),
        "relabelled_count": len(result.get("relabelled", ())),
        "retained_unchanged_count": len(result.get("retained_unchanged", ())),
        "label_transitions": dict(sorted(transitions.items())),
        "label_counts_v1": {label: v1_counts.get(label, 0) for label in labels},
        "label_counts_v2": {label: v2_counts.get(label, 0) for label in labels},
        "label_deltas": {
            label: v2_counts.get(label, 0) - v1_counts.get(label, 0) for label in labels
        },
    }


# --------------------------------------------------------------------------- #
# row accessors
# --------------------------------------------------------------------------- #


def _field(row: Mapping[str, Any], *keys: str) -> tuple[bool, str]:
    """Return ``(present, text)`` for the first key that exists in the row.

    ``present`` is False only when none of the keys exists in the mapping, which
    is what makes the difference between "cannot be reconstructed" (pending) and
    "present but blank" (fail-closed under the frozen v1 semantics).
    """
    for key in keys:
        if key in row:
            value = row.get(key)
            return True, "" if value is None else str(value).strip()
    return False, ""


def _accession(row: Mapping[str, Any], index: int) -> str:
    for key in ACCESSION_COLUMNS:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    raise ValueError(
        f"row {index} has no {'/'.join(ACCESSION_COLUMNS)} value; "
        "accession-level reporting requires one"
    )


def row_superfamily(row: Mapping[str, Any]) -> str:
    """The row's (historical) superfamily claim, or ``unassigned``."""
    present, value = _field(row, *SUPERFAMILY_COLUMNS)
    return value if present and value else "unassigned"


def row_family(row: Mapping[str, Any]) -> str:
    """The row's sequence-family call, or ``unassigned``."""
    present, value = _field(row, *FAMILY_COLUMNS)
    return value if present and value else "unassigned"


def row_transport_signal(row: Mapping[str, Any]) -> str:
    """The row's transport-signal prediction, or ``not_reported``."""
    present, value = _field(row, *TRANSPORT_COLUMNS)
    return value if present and value else "not_reported"


def candidate_localization_class(row: Mapping[str, Any]) -> str:
    """intracellular / extracellular / periplasmic / unassigned.

    Derived from the (historical) superfamily claim only, because that is the
    basis on which the frozen v1 76:24 candidate composition was computed.  It is
    a *candidate composition* partition, never a biological localization call.
    """
    superfamily = row_superfamily(row).lower()
    for prefix, label in LOCALIZATION_CLASSES:
        if superfamily.startswith(prefix):
            return label
    return "unassigned"


# --------------------------------------------------------------------------- #
# preregistered scenario rules (defined once, as data)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ScenarioContext:
    """Inputs a scenario rule needs that do not live inside a single row."""

    confounders: frozenset[str] | None = None


@dataclass(frozen=True)
class Scenario:
    """One preregistered scenario: a named predicate over a candidate row."""

    name: str
    description: str
    rule: str
    criteria: tuple[str, ...]
    removed_criteria: tuple[str, ...]
    source: str


@dataclass(frozen=True)
class ScenarioResult:
    """Outcome of one scenario over one table of rows."""

    scenario: str
    included: tuple[str, ...]
    excluded: tuple[str, ...]
    pending: tuple[str, ...]
    pending_criteria: Mapping[str, int]
    excluded_reasons: Mapping[str, int]
    first_failing_criteria: Mapping[str, str]
    included_rows: tuple[Mapping[str, Any], ...] = ()


CRITERIA: dict[str, Callable[[Mapping[str, Any], ScenarioContext], bool | str]] = {}


def _criterion(name: str):
    def register(function):
        CRITERIA[name] = function
        return function

    return register


@_criterion("v1_superfamily_recognized")
def _criterion_v1_superfamily_recognized(row, context):
    present, value = _field(row, *SUPERFAMILY_COLUMNS)
    if not present:
        return PENDING
    return value in V1_SUPERFAMILY_CRITERIA


@_criterion("v1_interpro_supported")
def _criterion_v1_interpro_supported(row, context):
    present, value = _field(row, "interpro_status")
    if not present:
        return PENDING
    return value == "interpro_supported"


@_criterion("v1_no_architecture_conflict")
def _criterion_v1_no_architecture_conflict(row, context):
    present, value = _field(row, "architecture_consistency")
    if not present:
        return PENDING
    return value != "conflicting"


@_criterion("v1_strong_profile")
def _criterion_v1_strong_profile(row, context):
    present, status = _field(row, "profile_evidence_status")
    if not present:
        return PENDING
    if status in V1_STRONG_PROFILE_STATUSES:
        return True
    if status != V1_DISCOVERY_STATUS:
        return False
    evalue_present, raw_evalue = _field(row, "profile_best_evalue", "best_evalue")
    if not evalue_present:
        return PENDING
    try:
        evalue = float(raw_evalue)
    except ValueError:
        return PENDING
    if not evalue < V1_DISCOVERY_EVALUE_THRESHOLD:
        return False
    confidence_present, confidence = _field(row, "superfamily_confidence")
    if confidence_present and confidence != V1_UNIQUE_CONFIDENCE:
        return False
    return True


@_criterion("v1_not_confounder")
def _criterion_v1_not_confounder(row, context):
    """Frozen confounder rule, evaluated over every available confounder source.

    A row-level flag (``is_confounder``) and the caller's ``--confounders`` set are
    two statements of the same claim, so a "yes" from either one excludes the row
    (fail-closed); the row passes only when every available source says "no"; and
    with no definite answer the criterion is pending.
    """
    verdicts: list[bool | None] = []
    present, value = _field(row, "is_confounder")
    if present:
        lowered = value.lower()
        if lowered in TRUTHY_VALUES:
            return False
        if lowered in FALSY_VALUES:
            verdicts.append(False)
        else:
            verdicts.append(None)
    if context.confounders is not None:
        accession_present, accession = _field(row, *ACCESSION_COLUMNS)
        if accession_present:
            verdicts.append(accession in context.confounders)
        else:
            verdicts.append(None)
    if any(verdict is True for verdict in verdicts):
        return False
    if verdicts and all(verdict is False for verdict in verdicts):
        return True
    return PENDING


@_criterion("v1_localization_export_requirement")
def _criterion_v1_localization_export_requirement(row, context):
    present, superfamily = _field(row, *SUPERFAMILY_COLUMNS)
    if not present:
        return PENDING
    if superfamily not in V1_EXPORT_REQUIRING:
        return True
    signal_present, signal = _field(row, *TRANSPORT_COLUMNS)
    if not signal_present:
        return PENDING
    return signal in V1_EXPORT_CLASSES


@_criterion("v1_localization_no_export_requirement")
def _criterion_v1_localization_no_export_requirement(row, context):
    present, superfamily = _field(row, *SUPERFAMILY_COLUMNS)
    if not present:
        return PENDING
    if superfamily not in V1_NO_EXPORT_REQUIRED:
        return True
    signal_present, signal = _field(row, *TRANSPORT_COLUMNS)
    if not signal_present:
        return PENDING
    return signal in V1_NO_EXPORT_CLASSES


def _v1_evaluate_architecture(field: str, expected: str, row: Mapping[str, Any]):
    """Frozen ``evaluate`` helper: True / False / PENDING."""
    if field == "x1_hydrophobic":
        present, value = _field(row, "lipase_box_x1")
        if not present:
            return PENDING
        return value in V1_HYDROPHOBIC
    if field == "pf06850":
        present, value = _field(row, "sbd_pf06850_binding_state")
        if not present:
            return PENDING
        return value == "detected"
    present, value = _field(row, field)
    if not present:
        return PENDING
    return value == expected


@_criterion("v1_architecture_criteria")
def _criterion_v1_architecture_criteria(row, context):
    present, superfamily = _field(row, *SUPERFAMILY_COLUMNS)
    if not present:
        return PENDING
    criteria = V1_SUPERFAMILY_CRITERIA.get(superfamily)
    if criteria is None:
        return False
    _requirement, required = criteria
    outcomes = [
        _v1_evaluate_architecture(field, expected, row)
        for field, expected in (item.split("=", 1) for item in required)
    ]
    if any(outcome is False for outcome in outcomes):
        return False
    if any(outcome is not True for outcome in outcomes):
        return PENDING
    return True


@_criterion("sequence_family_hmm_validated")
def _criterion_sequence_family_hmm_validated(row, context):
    present, value = _field(row, "model_layer")
    if not present:
        return PENDING
    return value == SEQUENCE_FAMILY_VALIDATED


V1_RULE_CRITERIA = (
    "v1_superfamily_recognized",
    "v1_interpro_supported",
    "v1_no_architecture_conflict",
    "v1_strong_profile",
    "v1_not_confounder",
    "v1_localization_export_requirement",
    "v1_localization_no_export_requirement",
    "v1_architecture_criteria",
)

SCENARIOS: dict[str, Scenario] = {
    "current_v1_rules": Scenario(
        name="current_v1_rules",
        description=(
            "The frozen v1 high-confidence rule set, reproduced from the frozen v1 "
            "filter module; this is the rule that produced the 40,690 table."
        ),
        rule=(
            "include a row iff its superfamily is one of the eight frozen v1 "
            "superfamilies AND interpro_status == 'interpro_supported' AND "
            "architecture_consistency != 'conflicting' AND the profile evidence is "
            "strong (profile_trained_hit / profile_ambiguous_family, or discovery_hit "
            "with best_evalue < 1e-10) AND the accession is not a frozen confounder "
            "AND the localization hold rule holds (export-requiring superfamilies need "
            "signalp_class in {SP, LIPO, TAT, TATLIPO}; no-export superfamilies need "
            "signalp_class in {OTHER, pending}) AND the frozen per-superfamily "
            "architecture criteria hold"
        ),
        criteria=V1_RULE_CRITERIA,
        removed_criteria=(),
        source=V1_SOURCE,
    ),
    "localization_uncertain_not_excluded": Scenario(
        name="localization_uncertain_not_excluded",
        description=(
            "The v1 rules with the missing-secretion-signal exclusion removed: a row "
            "is no longer held merely because an export-requiring superfamily has no "
            "predicted secretion signal."
        ),
        rule=(
            "identical to current_v1_rules with exactly one criterion removed: "
            "v1_localization_export_requirement. A row whose only v1 failure is "
            "'export-requiring superfamily without a predicted secretion signal' is "
            "included, while every other v1 criterion - including the no-export side "
            "of the same gate - still applies. This is the v1 rule that held 48,038 of "
            "67,720 extracellular dPHASCL type-1 candidates (71%) on the missing "
            "secretion signal alone"
        ),
        criteria=tuple(
            name for name in V1_RULE_CRITERIA if name != "v1_localization_export_requirement"
        ),
        removed_criteria=("v1_localization_export_requirement",),
        source=V1_SOURCE,
    ),
    "sequence_validated_models_only": Scenario(
        name="sequence_validated_models_only",
        description=(
            "Restrict the candidate set to rows whose model layer is the independently "
            "validated sequence-family layer."
        ),
        rule=(
            "include a row iff model_layer == 'sequence_family_hmm_validated'. The "
            "frozen v1 gates are deliberately NOT applied, because v2 does not use a "
            "binary high-confidence column; discovery_hmm_uncalibrated and "
            "reference_query_only rows are excluded. A row whose model_layer column is "
            "absent is counted as pending, never silently excluded or included"
        ),
        criteria=("sequence_family_hmm_validated",),
        removed_criteria=(),
        source=(
            "docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md "
            "('四层模型体系': sequence_family_hmm_validated is an independent layer)"
        ),
    ),
}


def resolve_scenario_names(names: Iterable[str] | str | None = None) -> tuple[str, ...]:
    """Validate a requested scenario subset; unknown names raise, never ignore.

    Accepts a single comma-separated string or an iterable of such strings and
    returns the requested names in preregistration order.
    """
    if names is None:
        return tuple(SCENARIOS)
    if isinstance(names, str):
        chunks: Sequence[str] = (names,)
    else:
        chunks = tuple(names)
    requested: list[str] = []
    for chunk in chunks:
        for name in str(chunk).split(","):
            name = name.strip()
            if name and name not in requested:
                requested.append(name)
    unknown = [name for name in requested if name not in SCENARIOS]
    if unknown:
        raise ValueError(
            f"unknown scenario name(s): {', '.join(sorted(unknown))}; "
            f"valid names: {', '.join(SCENARIOS)}"
        )
    wanted = set(requested)
    return tuple(name for name in SCENARIOS if name in wanted)


def excluded_reason(criterion: str, row: Mapping[str, Any]) -> str:
    """Frozen v1 hold reason for a failing criterion (or a scenario-specific one)."""
    if criterion == "v1_architecture_criteria":
        superfamily_present, superfamily = _field(row, *SUPERFAMILY_COLUMNS)
        lipase_present, lipase_box = _field(row, "lipase_box_state")
        ahsmg_present, ahsmg = _field(row, "ahsmg_state")
        if superfamily_present and lipase_present and ahsmg_present:
            if v1_nucleophile_type(superfamily, lipase_box, ahsmg) == V1_NT_CONFLICT:
                return "nucleophile_conflict"
    return EXCLUDED_REASON_BY_CRITERION.get(criterion, "not_in_v1_scope")


def evaluate_row(
    scenario: Scenario,
    row: Mapping[str, Any],
    context: ScenarioContext,
) -> tuple[bool, tuple[str, ...], str, str]:
    """Evaluate one scenario against one row.

    Returns ``(included, pending_criteria, first_failing_criterion, reason)``.  A
    row is included only when every applied criterion is True; a row with no
    failing criterion but at least one unevaluable criterion is *pending*, never
    silently included or excluded.
    """
    pending: list[str] = []
    for name in scenario.criteria:
        outcome = CRITERIA[name](row, context)
        if outcome is False:
            return False, tuple(pending), name, excluded_reason(name, row)
        if outcome is not True:
            pending.append(name)
    return (not pending), tuple(pending), "", ""


def evaluate_scenario(
    rows: Iterable[Mapping[str, Any]],
    scenario_name: str,
    context: ScenarioContext | None = None,
) -> ScenarioResult:
    """Evaluate one preregistered scenario over rows of mappings."""
    scenario = SCENARIOS.get(scenario_name)
    if scenario is None:
        raise ValueError(
            f"unknown scenario name: {scenario_name}; valid names: {', '.join(SCENARIOS)}"
        )
    active = context if context is not None else ScenarioContext()
    included: list[str] = []
    excluded: list[str] = []
    pending: list[str] = []
    pending_criteria: Counter = Counter()
    excluded_reasons: Counter = Counter()
    first_failing: dict[str, str] = {}
    included_rows: list[Mapping[str, Any]] = []
    for index, row in enumerate(rows):
        accession = _accession(row, index)
        is_included, row_pending, failing, reason = evaluate_row(scenario, row, active)
        for criterion in row_pending:
            pending_criteria[criterion] += 1
        if is_included:
            included.append(accession)
            included_rows.append(row)
        elif failing:
            excluded.append(accession)
            excluded_reasons[reason] += 1
            first_failing[accession] = failing
        else:
            pending.append(accession)
    return ScenarioResult(
        scenario=scenario_name,
        included=tuple(sorted(included)),
        excluded=tuple(sorted(excluded)),
        pending=tuple(sorted(pending)),
        pending_criteria=dict(sorted(pending_criteria.items())),
        excluded_reasons=dict(sorted(excluded_reasons.items())),
        first_failing_criteria=dict(sorted(first_failing.items())),
        included_rows=tuple(included_rows),
    )


# --------------------------------------------------------------------------- #
# composition blocks and the 76:24 guard
# --------------------------------------------------------------------------- #


def _composition(
    scenario: str,
    rule: str,
    field: str,
    counts: Mapping[str, int],
    total: int,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    ordered = {key: counts[key] for key in sorted(counts)}
    block: dict[str, Any] = {
        "basis": COMPOSITION_BASIS,
        "scenario": scenario,
        "rule": rule,
        "field": field,
        "denominator": f"{total} proteins selected by scenario {scenario}",
        "caveat": COMPOSITION_CAVEAT,
        "counts": ordered,
        "fractions": (
            {key: round(value / total, 6) for key, value in ordered.items()} if total else {}
        ),
    }
    if extra:
        block.update(extra)
    return block


def _localization_ratio_label(scenario: str, counts: Mapping[str, int], total: int) -> str:
    intracellular = counts.get("intracellular", 0)
    extracellular = counts.get("extracellular", 0)
    if total <= 0:
        return (
            "intracellular:extracellular candidate composition not computable under "
            f"{scenario}: no protein is selected by this rule"
        )
    return (
        "intracellular:extracellular candidate composition "
        f"{intracellular / total * 100:.1f}:{extracellular / total * 100:.1f} "
        f"under {scenario} (n={total})"
    )


def _walk(node: Any, path: tuple[str, ...] = ()) -> Iterator[tuple[tuple[str, ...], Any]]:
    yield path, node
    if isinstance(node, Mapping):
        for key, value in node.items():
            yield from _walk(value, path + (str(key),))
    elif isinstance(node, (list, tuple)):
        for index, value in enumerate(node):
            yield from _walk(value, path + (str(index),))


def _is_composition_path(path: tuple[str, ...]) -> bool:
    return len(path) == 4 and path[0] == "scenarios" and path[2] == "compositions"


def _looks_like_composition(node: Any) -> bool:
    if not isinstance(node, Mapping):
        return False
    if "counts" not in node:
        return False
    return any(key in node for key in ("basis", "caveat", "scenario", "denominator"))


def composition_blocks(report: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Every canonical composition block of a report, in deterministic order."""
    return [
        node
        for path, node in _walk(report)
        if _is_composition_path(path) and _looks_like_composition(node)
    ]


def _check_label_fragment(text: str, where: str) -> None:
    lowered = text.lower()
    for fragment in FORBIDDEN_LABEL_FRAGMENTS:
        if fragment in lowered:
            raise ValueError(
                f"forbidden label fragment {fragment!r} at {where}: the 76:24 "
                "candidate composition may never be labelled as a biological ratio"
            )


def validate_report_guard(report: Mapping[str, Any]) -> None:
    """Refuse a report that could present a candidate composition as biology.

    Raises ``ValueError`` when a composition block is missing its scenario name,
    rule text or caveat, when it lives outside
    ``scenarios[*].compositions[*]``, or when any report key/label carries a
    forbidden ratio fragment.
    """
    for path, node in _walk(report):
        where = "/".join(path) or "<root>"
        if isinstance(node, Mapping):
            for key, value in node.items():
                _check_label_fragment(str(key), f"{where}.{key}")
                if isinstance(value, str) and value != COMPOSITION_CAVEAT:
                    _check_label_fragment(value, f"{where}.{key}")
        if _looks_like_composition(node):
            if not _is_composition_path(path):
                raise ValueError(
                    f"composition block outside scenarios[*].compositions[*] at {where}: "
                    "a composition must always be reported together with its scenario "
                    "and the fixed caveat"
                )
            block = node
            if block.get("basis") != COMPOSITION_BASIS:
                raise ValueError(
                    f"composition block at {where} has basis {block.get('basis')!r}, "
                    f"expected {COMPOSITION_BASIS!r}"
                )
            if block.get("scenario") not in SCENARIOS:
                raise ValueError(
                    f"composition block at {where} must carry one of the preregistered "
                    f"scenario names {sorted(SCENARIOS)}, got {block.get('scenario')!r}"
                )
            if block.get("scenario") != path[1]:
                raise ValueError(
                    f"composition block at {where} names scenario "
                    f"{block.get('scenario')!r} but is filed under {path[1]!r}"
                )
            if not str(block.get("rule", "")).strip():
                raise ValueError(f"composition block at {where} is missing its rule text")
            if block.get("caveat") != COMPOSITION_CAVEAT:
                raise ValueError(
                    f"composition block at {where} is missing the fixed composition caveat"
                )


# --------------------------------------------------------------------------- #
# report assembly
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Analysis:
    """The JSON report plus the rows of the accession-level TSV."""

    report: dict[str, Any]
    accession_rows: tuple[dict[str, str], ...]


def _scenario_block(
    result: ScenarioResult,
    scenario: Scenario,
    others: set[str],
    compared_against: Sequence[str],
) -> dict[str, Any]:
    rule = scenario.rule
    total = len(result.included)
    superfamily_counts = Counter(row_superfamily(row) for row in result.included_rows)
    family_counts = Counter(row_family(row) for row in result.included_rows)
    transport_counts = Counter(row_transport_signal(row) for row in result.included_rows)
    localization_counts = Counter(
        candidate_localization_class(row) for row in result.included_rows
    )
    genomes = {
        str(row.get("genome")).strip()
        for row in result.included_rows
        if str(row.get("genome") or "").strip()
    }
    genomes_missing = sum(
        1 for row in result.included_rows if not str(row.get("genome") or "").strip()
    )
    selected = set(result.included)
    only_here = sorted(selected - others)
    not_here = sorted(others - selected)
    counts_status = (
        "final"
        if not result.pending
        else f"lower_bound_{len(result.pending)}_rows_pending_not_decidable"
    )
    return {
        "scenario": result.scenario,
        "description": scenario.description,
        "rule": rule,
        "criteria": list(scenario.criteria),
        "removed_criteria": list(scenario.removed_criteria),
        "source": scenario.source,
        "total_proteins": total,
        "total_proteins_status": counts_status,
        "genomes": len(genomes),
        "genomes_status": counts_status,
        "genome_values_missing": genomes_missing,
        "pending_proteins": len(result.pending),
        "pending_criteria": dict(result.pending_criteria),
        "excluded_proteins": len(result.excluded),
        "excluded_reason_counts": dict(result.excluded_reasons),
        "first_failing_criteria": dict(result.first_failing_criteria),
        "compositions": {
            "superfamily": _composition(
                result.scenario, rule, "superfamily", superfamily_counts, total
            ),
            "family": _composition(result.scenario, rule, "family", family_counts, total),
            "transport_signal": _composition(
                result.scenario, rule, "transport_signal", transport_counts, total
            ),
            "candidate_localization": _composition(
                result.scenario,
                rule,
                "candidate_localization",
                localization_counts,
                total,
                extra={
                    "ratio_label": _localization_ratio_label(
                        result.scenario, localization_counts, total
                    ),
                    "partition_basis": (
                        "(historical) superfamily claim only; a candidate composition "
                        "partition, not a biological localization call"
                    ),
                },
            ),
        },
        "differences_vs_other_scenarios": {
            "compared_against": list(compared_against),
            "only_in_scenario_count": len(only_here),
            "not_in_scenario_count": len(not_here),
            "accession_list_location": ACCESSIONS_TSV,
        },
    }


def analyze(
    *,
    rows: Iterable[Mapping[str, Any]],
    v1_labels: Mapping[str, str] | None = None,
    v2_labels: Mapping[str, str] | None = None,
    scenario_names: Iterable[str] | str | None = None,
    context: ScenarioContext | None = None,
) -> Analysis:
    """Build the comparison report and the accession-level TSV rows.

    ``rows`` are plain mappings (dicts read from any table): this module never
    imports or inspects another pipeline module, so the v2 catalog builder stays
    free to change its internals.
    """
    row_list = list(rows)
    names = resolve_scenario_names(scenario_names)
    first = dict(v1_labels or {})
    second = dict(v2_labels or {})
    impact = compare(v1=first, v2=second)

    results = {name: evaluate_scenario(row_list, name, context) for name in names}
    selections = {name: set(result.included) for name, result in results.items()}

    scenarios: dict[str, Any] = {}
    for name in names:
        others: set[str] = set()
        for other in names:
            if other != name:
                others |= selections[other]
        compared_against = [other for other in names if other != name]
        scenarios[name] = _scenario_block(
            results[name], SCENARIOS[name], others, compared_against
        )

    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generator": GENERATOR,
        "claim_boundary": CLAIM_BOUNDARY,
        "preregistered_scenarios": list(SCENARIOS),
        "selected_scenarios": list(names),
        "scenario_definitions": {
            name: {
                "description": scenario.description,
                "rule": scenario.rule,
                "criteria": list(scenario.criteria),
                "removed_criteria": list(scenario.removed_criteria),
                "source": scenario.source,
            }
            for name, scenario in SCENARIOS.items()
        },
        "v1_v2_impact": impact_summary(impact),
        "scenarios": scenarios,
        "row_count": len(row_list),
        "accession_list_location": ACCESSIONS_TSV,
        "guard": {
            "composition_basis": COMPOSITION_BASIS,
            "caveat": COMPOSITION_CAVEAT,
            "rule": (
                "every composition is reported together with the scenario name that "
                "produced it and with the fixed caveat; no key or label may present "
                "the candidate composition as a biological ratio"
            ),
        },
    }

    accession_rows: list[dict[str, str]] = []
    for change_kind in ("removed", "added", "relabelled"):
        for accession in sorted(impact[change_kind]):
            accession_rows.append(
                {
                    "section": "v1_v2",
                    "accession": accession,
                    "change_kind": change_kind,
                    "scenario": "",
                    "v1_label": first.get(accession, ""),
                    "v2_label": second.get(accession, ""),
                }
            )
    for name in names:
        selected = selections[name]
        others = set()
        for other in names:
            if other != name:
                others |= selections[other]
        for accession in sorted(selected - others):
            accession_rows.append(
                {
                    "section": "scenario",
                    "accession": accession,
                    "change_kind": "only_in_scenario",
                    "scenario": name,
                    "v1_label": "",
                    "v2_label": "",
                }
            )
        for accession in sorted(others - selected):
            accession_rows.append(
                {
                    "section": "scenario",
                    "accession": accession,
                    "change_kind": "not_in_scenario",
                    "scenario": name,
                    "v1_label": "",
                    "v2_label": "",
                }
            )
    return Analysis(report=report, accession_rows=tuple(accession_rows))


def build_report(
    rows: Iterable[Mapping[str, Any]],
    v1_labels: Mapping[str, str] | None = None,
    v2_labels: Mapping[str, str] | None = None,
    scenario_names: Iterable[str] | str | None = None,
    context: ScenarioContext | None = None,
) -> dict[str, Any]:
    """The JSON report alone (see :func:`analyze` for the accompanying TSV rows)."""
    return analyze(
        rows=rows,
        v1_labels=v1_labels,
        v2_labels=v2_labels,
        scenario_names=scenario_names,
        context=context,
    ).report


# --------------------------------------------------------------------------- #
# table IO and output guards
# --------------------------------------------------------------------------- #


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV into plain dicts (``utf-8-sig`` tolerant)."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def read_label_map(
    path: Path,
    column: str,
    truthy_label: str | None = None,
) -> dict[str, str]:
    """Read ``accession -> label`` from a TSV, optionally mapping truthy values.

    Missing label column, missing accession values and duplicate accessions all
    raise: the comparison must never silently drop or double-count a row.
    """
    labels: dict[str, str] = {}
    rows = read_rows(path)
    if not rows:
        raise ValueError(f"{path}: no data rows")
    if column not in rows[0]:
        raise ValueError(
            f"{path}: label column {column!r} is absent; available columns: "
            f"{', '.join(sorted(rows[0]))}"
        )
    for index, row in enumerate(rows):
        try:
            accession = _accession(row, index + 2)
        except ValueError as error:
            raise ValueError(f"{path}: {error}") from error
        if accession in labels:
            raise ValueError(f"{path}: duplicate accession {accession!r}")
        value = row.get(column)
        value = "" if value is None else str(value).strip()
        if truthy_label is not None and value.lower() in TRUTHY_VALUES:
            value = truthy_label
        labels[accession] = value
    return labels


def read_accession_set(path: Path) -> frozenset[str]:
    """Read a one-accession-per-row TSV into a set."""
    accessions: set[str] = set()
    for index, row in enumerate(read_rows(path)):
        try:
            accessions.add(_accession(row, index + 2))
        except ValueError as error:
            raise ValueError(f"{path}: {error}") from error
    return frozenset(accessions)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bound_file(path: Path | None) -> dict[str, Any]:
    """Provenance record for an input path (``pending``/``not_supplied`` if absent)."""
    if path is None:
        return {"path": None, "status": "not_supplied", "sha256": None, "size": None}
    candidate = Path(path)
    if not candidate.is_file():
        return {"path": str(candidate), "status": "pending", "sha256": None, "size": None}
    return {
        "path": str(candidate),
        "status": "verified",
        "sha256": sha256_file(candidate),
        "size": candidate.stat().st_size,
    }


def require_file(path: Path | None, flag: str) -> Path:
    if path is None:
        raise ValueError(f"{flag} is required")
    candidate = Path(path)
    if not candidate.is_file():
        raise ValueError(f"{flag}: input path does not exist: {candidate}")
    return candidate


def prepare_output_dir(path: Path) -> Path:
    """Create the output dir, or refuse an existing non-empty one."""
    candidate = Path(path)
    if candidate.exists():
        if not candidate.is_dir():
            raise ValueError(f"--output-dir exists and is not a directory: {candidate}")
        existing = sorted(item.name for item in candidate.iterdir())
        if existing:
            raise ValueError(
                f"--output-dir is not empty; refusing to overwrite: {candidate} "
                f"(contains {existing[:5]})"
            )
    else:
        candidate.mkdir(parents=True, exist_ok=True)
    return candidate


def output_path(output_dir: Path, filename: str) -> Path:
    """Resolve a report filename inside ``output_dir`` (never outside it)."""
    if filename not in OUTPUT_FILENAMES:
        raise ValueError(f"refusing to write unexpected filename: {filename}")
    base = Path(output_dir).resolve()
    target = (base / filename).resolve()
    if target.parent != base:
        raise ValueError(f"refusing to write outside --output-dir: {target}")
    return target


def write_accession_tsv(path: Path, accession_rows: Iterable[Mapping[str, str]]) -> None:
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(ACCESSION_TSV_COLUMNS),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in accession_rows:
            writer.writerow({column: row.get(column, "") for column in ACCESSION_TSV_COLUMNS})


def write_outputs(
    output_dir: Path,
    report: Mapping[str, Any],
    accession_rows: Iterable[Mapping[str, str]],
) -> dict[str, str]:
    """Write the JSON report and the accession-level TSV into ``output_dir``."""
    report_path = output_path(output_dir, REPORT_JSON)
    tsv_path = output_path(output_dir, ACCESSIONS_TSV)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_accession_tsv(tsv_path, accession_rows)
    return {"report": str(report_path), "accessions": str(tsv_path)}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--v1", type=Path, default=None, help="frozen v1 label table (TSV)")
    parser.add_argument("--v2", type=Path, default=None, help="v2 catalog (TSV)")
    parser.add_argument(
        "--evidence",
        type=Path,
        default=None,
        help="rows the scenarios run over (default: the --v2 table)",
    )
    parser.add_argument("--confounders", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--scenarios",
        action="append",
        default=None,
        help="comma-separated scenario subset; repeatable; unknown names raise",
    )
    parser.add_argument("--v1-label-column", default="high_confidence")
    parser.add_argument("--v2-label-column", default="primary_disposition")
    parser.add_argument("--list-scenarios", action="store_true")
    return parser


def run(args: argparse.Namespace) -> dict[str, Any]:
    v1_path = require_file(args.v1, "--v1")
    v2_path = require_file(args.v2, "--v2")
    evidence_path = require_file(args.evidence or args.v2, "--evidence")
    confounders_path = (
        require_file(args.confounders, "--confounders") if args.confounders else None
    )
    names = resolve_scenario_names(args.scenarios)
    if args.output_dir is None:
        raise ValueError("--output-dir is required")
    output_dir = prepare_output_dir(args.output_dir)

    v1_labels = read_label_map(v1_path, args.v1_label_column, truthy_label=HIGH_CONFIDENCE_LABEL)
    v2_labels = read_label_map(v2_path, args.v2_label_column)
    rows = read_rows(evidence_path)
    context = ScenarioContext(
        confounders=read_accession_set(confounders_path) if confounders_path else None
    )

    analysis = analyze(
        rows=rows,
        v1_labels=v1_labels,
        v2_labels=v2_labels,
        scenario_names=names,
        context=context,
    )
    report = dict(analysis.report)
    report["inputs"] = {
        "v1": bound_file(v1_path),
        "v2": bound_file(v2_path),
        "evidence": bound_file(evidence_path),
        "confounders": bound_file(confounders_path),
        "v1_label_column": args.v1_label_column,
        "v2_label_column": args.v2_label_column,
    }
    report["outputs"] = {
        "report": str(output_path(output_dir, REPORT_JSON)),
        "accessions": str(output_path(output_dir, ACCESSIONS_TSV)),
    }
    validate_report_guard(report)
    write_outputs(output_dir, report, analysis.accession_rows)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.list_scenarios:
        print(
            json.dumps(
                {name: scenario.description for name, scenario in SCENARIOS.items()},
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    try:
        report = run(args)
    except (ValueError, FileNotFoundError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
