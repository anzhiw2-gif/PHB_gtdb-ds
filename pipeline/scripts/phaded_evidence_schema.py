#!/usr/bin/env python3
"""Shared PhaDED evidence ontology: controlled vocabularies and row validators.

This module is the single place that decides what a PhaDED evidence field may
mean.  It is standard-library only and never touches project data.

Boundaries
----------
* Rows are plain ``dict[str, str]`` mappings exactly as produced by
  ``csv.DictReader``.  A column that is *absent* from a row is not validated and
  is reported with its documented weakest default; a column that is *present*
  must hold a controlled value, otherwise :class:`ValueError` names the
  offending field and value.
* Row set completeness (for example "every catalog protein carries exactly one
  ``primary_disposition``") belongs to the catalog builder (Task 7) and the flow
  reconciler (Task 8), not to this module.

AHSMG is a Ser motif class
--------------------------
``AHSMG`` is a motif class, never a nucleophile identity:
``NUCLEOPHILE_IDENTITIES`` holds only ``ser``/``cys``/``unresolved``, both
validators reject ``nucleophile_identity=ahsmg``, and
:func:`normalize_nucleophile` maps an AHSMG motif to ``ser`` plus
``motif_class=AHSMG``.

SignalP is a transport prediction, not localization truth
---------------------------------------------------------
:func:`normalize_localization` never infers ``localization_evidence`` from a
SignalP class.  A row may claim ``localization_evidence=experimental_confirmed``
only when it also names a non-SignalP ``localization_evidence_source``.

The discovery layer only scores
-------------------------------
``discovery_hmm_uncalibrated`` may recall and score existing intermediate
products; it never makes a family call and never performs an exclusion action
(AGENTS.md).  :func:`validate_candidate_evidence` rejects a discovery-layer row
that carries a ``sequence_family_call`` or ``family_call``, and rejects an
``excluded_input_quality`` disposition unless the row documents an
input-quality reason.

Sequence-model gates
--------------------
:func:`load_sequence_model_gates` validates
``pipeline/config/phaded_sequence_model_gates.tsv``.  The preregistered
``sequence-gate-v1`` row sets ``min_heldout_recall`` to ``0.8``, but **for
held-out panels smaller than five sequences every held-out sequence must be
recovered**; that rule is encoded in the evaluator
(``evaluate_phaded_sequence_models``, Task 6), not inferred from rounding of
``0.8``.  ``SMALL_HELDOUT_PANEL_THRESHOLD`` records the threshold so Task 6 does
not have to redefine it.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Mapping, Sequence


# ---------------------------------------------------------------------------
# Controlled vocabularies
# ---------------------------------------------------------------------------

MODEL_LAYERS = {
    "reference_query_only",
    "discovery_hmm_uncalibrated",
    "sequence_family_hmm_validated",
    "calibrated_candidate_model",
}
FUNCTIONAL_CALIBRATION_STATUSES = {
    "not_function_calibrated",
    "candidate_gate_passed_not_promoted",
    "calibrated_candidate_model",
    # A contradictory direct experimental negative blocks functional
    # calibration while the sequence itself stays in the reference audit
    # (build_phaded_profiles.py).  It is a *blocked* state, not a calibrated
    # one, and it must remain visible rather than being collapsed into
    # ``not_function_calibrated``.
    "blocked_contradictory_experimental_negative",
}
#: Sequence layer and functional calibration status are independent fields, but
#: they may not contradict each other.  A functional calibration status of
#: ``calibrated_candidate_model`` exists only on the calibrated layer, so the
#: discovery and reference layers can never look functionally calibrated.
ALLOWED_MODEL_STATES = {
    "reference_query_only": frozenset({
        "not_function_calibrated",
        "blocked_contradictory_experimental_negative",
    }),
    "discovery_hmm_uncalibrated": frozenset({
        "not_function_calibrated",
        "candidate_gate_passed_not_promoted",
        "blocked_contradictory_experimental_negative",
    }),
    "sequence_family_hmm_validated": frozenset({
        "not_function_calibrated",
        "candidate_gate_passed_not_promoted",
        "blocked_contradictory_experimental_negative",
    }),
    "calibrated_candidate_model": frozenset({"calibrated_candidate_model"}),
}
NUCLEOPHILE_IDENTITIES = {"ser", "cys", "unresolved"}
MOTIF_CLASSES = {"GxSxG", "AHSMG", "Cys-associated", "other", "unresolved"}
EVIDENCE_GRADES = {"E3", "E2", "E1", "A"}
PRIMARY_DISPOSITIONS = {
    "core_sequence_homolog",
    "probable_sequence_homolog",
    "remote_homolog_candidate",
    "function_unresolved",
    "deferred_structure_review",
    "excluded_input_quality",
}
TRANSPORT_SIGNAL_PREDICTIONS = {"SP", "LIPO", "TAT", "TATLIPO", "OTHER", "not_tested", "invalid_input"}
LOCALIZATION_EVIDENCE = {
    "experimental_confirmed",
    "other_evidence",
    "historical_family_label_only",
    "unknown",
}
CATALYTIC_DOMAIN_TYPES = {"type1_verified", "type2_verified", "other", "unresolved"}

#: ``excluded_input_quality`` is limited to missing, truncated or malformed
#: input.  Scientific uncertainty is never an input error.
INPUT_QUALITY_REASONS = {
    "sequence_missing",
    "sequence_truncated",
    "sequence_too_short",
    "empty_sequence",
    "invalid_residues",
    "duplicate_input",
    "malformed_record",
    "input_error",
}
#: Columns that may carry the input-quality reason of an exclusion.
INPUT_QUALITY_REASON_FIELDS = ("input_quality_reason", "excluded_input_quality_reason", "exclusion_reason")
#: Field aliases that must never be used to encode an exclusion.
FAMILY_CALL_FIELDS = ("sequence_family_call", "family_call")

#: Sequence-model gates (Task 3 config, consumed by Task 6).
SEQUENCE_MODEL_GATE_VERSION = "sequence-gate-v1"
SEQUENCE_MODEL_GATE_COLUMNS = (
    "gate_version",
    "min_training_sequences",
    "min_heldout_sequences",
    "min_heldout_recall",
    "max_unexplained_confounder_hits",
    "require_unique_best_model",
    "require_alignment_hash",
)
SEQUENCE_MODEL_GATE_INTEGER_COLUMNS = (
    "min_training_sequences",
    "min_heldout_sequences",
    "max_unexplained_confounder_hits",
)
SEQUENCE_MODEL_GATE_POSITIVE_COLUMNS = ("min_training_sequences", "min_heldout_sequences")
SEQUENCE_MODEL_GATE_BOOLEAN_COLUMNS = ("require_unique_best_model", "require_alignment_hash")
#: Held-out panels smaller than this must recover every held-out sequence; the
#: evaluator implements the rule and must not infer it from 0.8 rounding.
SMALL_HELDOUT_PANEL_THRESHOLD = 5
SEQUENCE_MODEL_GATES_PATH = Path(__file__).resolve().parents[1] / "config" / "phaded_sequence_model_gates.tsv"


# ---------------------------------------------------------------------------
# Vocabulary normalization tables
# ---------------------------------------------------------------------------

def _key(value: object) -> str:
    """Return a case/separator-insensitive lookup key for a controlled value."""
    return re.sub(r"[^a-z0-9]+", "_", str("" if value is None else value).strip().lower()).strip("_")


NUCLEOPHILE_IDENTITY_ALIASES = {
    "unresolved": "unresolved",
    "unknown": "unresolved",
    "none": "unresolved",
    "pending": "unresolved",
    "ser": "ser",
    "s": "ser",
    "serine": "ser",
    "cys": "cys",
    "c": "cys",
    "cysteine": "cys",
}
MOTIF_CLASS_ALIASES = {
    "unresolved": "unresolved",
    "unknown": "unresolved",
    "none": "unresolved",
    "absent": "unresolved",
    "not_found": "unresolved",
    "not_supported": "unresolved",
    "gxsxg": "GxSxG",
    "gxsxg_like": "GxSxG",
    "gxsxg_motif": "GxSxG",
    "lipase_box": "GxSxG",
    "canonical_ser_motif": "GxSxG",
    "ahsmg": "AHSMG",
    "ahsmg_like": "AHSMG",
    "ahsmg_motif": "AHSMG",
    "cys_associated": "Cys-associated",
    "cys_associated_motif": "Cys-associated",
    "cys_motif": "Cys-associated",
    "cys": "Cys-associated",
    "other": "other",
    "noncanonical": "other",
    "non_canonical": "other",
}
TRANSPORT_SIGNAL_ALIASES = {
    "sp": "SP",
    "lipo": "LIPO",
    "lipoprotein": "LIPO",
    "spii": "LIPO",
    "tat": "TAT",
    "tatlipo": "TATLIPO",
    "tat_lipo": "TATLIPO",
    "other": "OTHER",
    "no_signal": "OTHER",
    "none": "OTHER",
    # SignalP 6 reports a type-IV pilin class that is not part of this project
    # vocabulary; it is recorded as the non-Sec/non-Tat/non-Lipo class.
    "pilin": "OTHER",
    "not_tested": "not_tested",
    "untested": "not_tested",
    "pending": "not_tested",
    "invalid_input": "invalid_input",
    "invalid": "invalid_input",
    "error": "invalid_input",
}
LOCALIZATION_EVIDENCE_ALIASES = {
    "unknown": "unknown",
    "unresolved": "unknown",
    "none": "unknown",
    "pending": "unknown",
    "experimental_confirmed": "experimental_confirmed",
    "experimental": "experimental_confirmed",
    "confirmed": "experimental_confirmed",
    "other_evidence": "other_evidence",
    "other": "other_evidence",
    "historical_family_label_only": "historical_family_label_only",
    "historical_family_label": "historical_family_label_only",
    "historical_label": "historical_family_label_only",
    "historical": "historical_family_label_only",
}
CATALYTIC_DOMAIN_TYPE_ALIASES = {
    "unresolved": "unresolved",
    "unknown": "unresolved",
    "none": "unresolved",
    "type1": "type1_verified",
    "type_1": "type1_verified",
    "typei": "type1_verified",
    "type_i": "type1_verified",
    "type1_verified": "type1_verified",
    "type2": "type2_verified",
    "type_2": "type2_verified",
    "typeii": "type2_verified",
    "type_ii": "type2_verified",
    "type2_verified": "type2_verified",
    "other": "other",
}
INPUT_QUALITY_REASON_ALIASES = {
    "missing": "sequence_missing",
    "missing_sequence": "sequence_missing",
    "truncated": "sequence_truncated",
    "truncation": "sequence_truncated",
    "short": "sequence_too_short",
    "too_short": "sequence_too_short",
    "empty": "empty_sequence",
    "empty_sequence_record": "empty_sequence",
    "bad_residues": "invalid_residues",
    "nonstandard_residues": "invalid_residues",
    "duplicate": "duplicate_input",
    "duplicated_input": "duplicate_input",
    "malformed": "malformed_record",
    "malformed_input": "malformed_record",
    "wrong_input": "input_error",
}
#: Reasons that describe scientific uncertainty rather than broken input.
SCIENTIFIC_UNCERTAINTY_REASON_KEYS = {
    "function_unresolved",
    "unknown_function",
    "scientific_uncertainty",
    "low_confidence",
    "ambiguous_family",
    "remote_homolog",
    "no_family_call",
    "deferred_structure_review",
}
_TRUE_TOKENS = {"true", "yes", "y", "1", "supported", "present", "detected", "eligible"}
_FALSE_TOKENS = {
    "false", "no", "n", "0", "unsupported", "not_supported", "absent",
    "not_detected", "not_eligible", "ineligible", "pending", "unknown", "none",
}

_MODEL_LAYER_BY_KEY = {_key(item): item for item in MODEL_LAYERS}
_FUNCTIONAL_STATUS_BY_KEY = {_key(item): item for item in FUNCTIONAL_CALIBRATION_STATUSES}
_NUCLEOPHILE_IDENTITY_BY_KEY = {_key(item): item for item in NUCLEOPHILE_IDENTITIES}
_MOTIF_CLASS_BY_KEY = {_key(item): item for item in MOTIF_CLASSES}
_EVIDENCE_GRADE_BY_KEY = {_key(item): item for item in EVIDENCE_GRADES}
_PRIMARY_DISPOSITION_BY_KEY = {_key(item): item for item in PRIMARY_DISPOSITIONS}
_TRANSPORT_SIGNAL_BY_KEY = {_key(item): item for item in TRANSPORT_SIGNAL_PREDICTIONS}
_LOCALIZATION_EVIDENCE_BY_KEY = {_key(item): item for item in LOCALIZATION_EVIDENCE}
_CATALYTIC_DOMAIN_TYPE_BY_KEY = {_key(item): item for item in CATALYTIC_DOMAIN_TYPES}
_INPUT_QUALITY_REASON_BY_KEY = {_key(item): item for item in INPUT_QUALITY_REASONS}

#: A source made only of SignalP tokens proves a transport prediction, never
#: experimental localization.  Mixed sources are left for human review.
_SIGNALP_ONLY_SOURCE = re.compile(
    r"signalp[0-9a-z]*(?:_(?:only|prediction|predictions|predicted|class|output|call|based|signal))*"
)


# ---------------------------------------------------------------------------
# Row helpers
# ---------------------------------------------------------------------------

def _text(row: Mapping[str, object], field: str) -> str:
    """Read a TSV field defensively; absent and ``None`` columns read as ``""``."""
    if field not in row:
        return ""
    return _text_value(row[field])


def _text_value(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _resolve(
    value: object,
    field: str,
    canonical: Mapping[str, str],
    allowed: "set[str] | frozenset[str]",
    *,
    default: str | None = None,
    aliases: Mapping[str, str] | None = None,
    hint: str = "",
) -> str:
    raw = _text_value(value)
    if not raw and default is not None:
        return default
    key = _key(raw)
    if aliases is not None and key in aliases:
        return aliases[key]
    if key in canonical:
        return canonical[key]
    raise ValueError(
        f"{field}={raw!r} is not a recognized value; expected one of {sorted(allowed)}{hint}"
    )


def _resolve_nucleophile_identity(value: object, field: str, *, default: str | None = None) -> str:
    raw = _text_value(value)
    if not raw and default is not None:
        return default
    if _key(raw) == "ahsmg":
        raise ValueError(
            f"{field}='ahsmg' is invalid: AHSMG is a Ser motif class (see MOTIF_CLASSES), "
            f"never a third nucleophile identity; use one of {sorted(NUCLEOPHILE_IDENTITIES)}"
        )
    return _resolve(
        raw, field, _NUCLEOPHILE_IDENTITY_BY_KEY, NUCLEOPHILE_IDENTITIES,
        aliases=NUCLEOPHILE_IDENTITY_ALIASES,
    )


def _resolve_optional_flag(row: Mapping[str, object], field: str) -> str:
    raw = _text(row, field)
    if not raw:
        return ""
    key = _key(raw)
    if key in _TRUE_TOKENS:
        return "true"
    if key in _FALSE_TOKENS:
        return "false"
    raise ValueError(f"{field}={raw!r} is not a recognized flag; expected 'true' or 'false'")


def _resolve_input_quality_reason(row: Mapping[str, object]) -> str:
    for field in INPUT_QUALITY_REASON_FIELDS:
        raw = _text(row, field)
        if not raw:
            continue
        key = _key(raw)
        if key in SCIENTIFIC_UNCERTAINTY_REASON_KEYS:
            raise ValueError(
                f"{field}={raw!r} must not be recorded as an input-quality exclusion: "
                "primary_disposition='excluded_input_quality' is limited to missing, truncated or "
                "malformed input, and scientific uncertainty is not an input error"
            )
        if key in INPUT_QUALITY_REASON_ALIASES:
            return INPUT_QUALITY_REASON_ALIASES[key]
        if key in _INPUT_QUALITY_REASON_BY_KEY:
            return _INPUT_QUALITY_REASON_BY_KEY[key]
        raise ValueError(
            f"{field}={raw!r} is not a recognized input-quality reason; "
            f"expected one of {sorted(INPUT_QUALITY_REASONS)}"
        )
    return ""


def _is_signalp_only_source(source: str) -> bool:
    return bool(_SIGNALP_ONLY_SOURCE.fullmatch(_key(source)))


def _resolve_localization(row: Mapping[str, object]) -> "tuple[str, str]":
    localization = _resolve(
        _text(row, "localization_evidence"), "localization_evidence",
        _LOCALIZATION_EVIDENCE_BY_KEY, LOCALIZATION_EVIDENCE,
        default="unknown", aliases=LOCALIZATION_EVIDENCE_ALIASES,
    )
    source = _text(row, "localization_evidence_source")
    if localization == "experimental_confirmed":
        if not source:
            raise ValueError(
                "localization_evidence_source is required when "
                "localization_evidence='experimental_confirmed': a SignalP transport prediction "
                "never establishes experimental localization"
            )
        if _is_signalp_only_source(source):
            raise ValueError(
                f"localization_evidence_source={source!r} cannot support "
                "localization_evidence='experimental_confirmed': SignalP predicts a transport "
                "signal and must not be promoted to experimental localization evidence"
            )
    return localization, source


def _resolve_family_call(row: Mapping[str, object]) -> "tuple[str, str]":
    for field in FAMILY_CALL_FIELDS:
        raw = _text(row, field)
        if raw:
            return field, raw
    return "", ""


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def validate_model_state(model_layer: str, functional_calibration_status: str) -> None:
    """Reject unknown or contradictory sequence-layer / functional-status pairs.

    Sequence discovery, sequence-family validation and functional calibration
    are separate claims.  A row may not carry a functional calibration status on
    a layer that cannot support it, and ``calibrated_candidate_model`` is only
    ever emitted on the calibrated layer by an explicitly authorized action.
    """
    layer = _text_value(model_layer)
    status = _text_value(functional_calibration_status)
    if layer not in MODEL_LAYERS:
        raise ValueError(
            f"model_layer={layer!r} is not a recognized model layer; expected one of {sorted(MODEL_LAYERS)}"
        )
    if status not in FUNCTIONAL_CALIBRATION_STATUSES:
        raise ValueError(
            "functional_calibration_status="
            f"{status!r} is not a recognized functional calibration status; "
            f"expected one of {sorted(FUNCTIONAL_CALIBRATION_STATUSES)}"
        )
    if status not in ALLOWED_MODEL_STATES[layer]:
        raise ValueError(
            f"functional_calibration_status={status!r} is not allowed with model_layer={layer!r}; "
            f"allowed statuses for this layer: {sorted(ALLOWED_MODEL_STATES[layer])}"
        )


def normalize_nucleophile(superfamily: str, motif_state: str, ahsmg_state: str) -> dict:
    """Normalize nucleophile identity and motif class without merging the two.

    ``AHSMG`` is a Ser motif class: an AHSMG record yields
    ``{"nucleophile_identity": "ser", "motif_class": "AHSMG"}`` and never a
    third nucleophile identity.  When no motif evidence is informative the
    superfamily supplies the prior: SCL depolymerase families are Ser
    hydrolases, while MCL families keep an unresolved nucleophile.
    """
    ahsmg = False
    raw_ahsmg = _text_value(ahsmg_state)
    if raw_ahsmg:
        key = _key(raw_ahsmg)
        if key in _TRUE_TOKENS:
            ahsmg = True
        elif key in _FALSE_TOKENS:
            ahsmg = False
        else:
            raise ValueError(
                f"ahsmg_state={raw_ahsmg!r} is not a recognized state; "
                "expected a supported/unsupported style flag such as 'supported', 'true' or 'unsupported'"
            )
    motif_class = _resolve(
        motif_state, "motif_state", _MOTIF_CLASS_BY_KEY, MOTIF_CLASSES,
        default="unresolved", aliases=MOTIF_CLASS_ALIASES,
    )
    if ahsmg:
        # AHSMG is the more specific observation; it always wins over the
        # generic lipase-box/motif state.
        motif_class = "AHSMG"

    if motif_class == "AHSMG":
        identity = "ser"
    elif motif_class == "GxSxG":
        identity = "ser"
    elif motif_class == "Cys-associated":
        identity = "cys"
    else:
        superfamily_key = _key(superfamily)
        if "scl" in superfamily_key:
            identity = "ser"
        else:
            identity = "unresolved"
    return {"nucleophile_identity": identity, "motif_class": motif_class}


def normalize_localization(signalp_class: str, localization_evidence: str) -> dict:
    """Normalize a SignalP class and the independent localization evidence.

    The SignalP class is recorded as ``transport_signal_prediction`` only.  This
    function never derives ``localization_evidence`` from it: the caller must
    supply the localization evidence explicitly, and
    :func:`validate_candidate_evidence` / :func:`validate_reference_evidence`
    require a non-SignalP source before ``experimental_confirmed`` is accepted.
    """
    prediction = _resolve(
        signalp_class, "signalp_class", _TRANSPORT_SIGNAL_BY_KEY, TRANSPORT_SIGNAL_PREDICTIONS,
        default="not_tested", aliases=TRANSPORT_SIGNAL_ALIASES,
    )
    evidence = _resolve(
        localization_evidence, "localization_evidence",
        _LOCALIZATION_EVIDENCE_BY_KEY, LOCALIZATION_EVIDENCE,
        default="unknown", aliases=LOCALIZATION_EVIDENCE_ALIASES,
    )
    return {"transport_signal_prediction": prediction, "localization_evidence": evidence}


def validate_reference_evidence(row: dict) -> dict:
    """Validate the evidence fields of one reference-ledger row.

    Columns that the ledger does not carry yet are reported as ``""`` (or their
    documented weakest default) instead of being invented, so this validator is
    usable both before and after the Task 4 ledger extension.
    """
    layer = _resolve(
        _text(row, "model_layer"), "model_layer", _MODEL_LAYER_BY_KEY, MODEL_LAYERS,
        default="reference_query_only",
    )
    status = _resolve(
        _text(row, "functional_calibration_status"), "functional_calibration_status",
        _FUNCTIONAL_STATUS_BY_KEY, FUNCTIONAL_CALIBRATION_STATUSES,
        default="not_function_calibrated",
    )
    validate_model_state(layer, status)

    identity = _resolve_nucleophile_identity(
        _text(row, "nucleophile_identity"), "nucleophile_identity", default="unresolved",
    )
    motif_class = _resolve(
        _text(row, "motif_class"), "motif_class", _MOTIF_CLASS_BY_KEY, MOTIF_CLASSES,
        default="unresolved", aliases=MOTIF_CLASS_ALIASES,
    )
    if motif_class == "AHSMG" and identity == "cys":
        raise ValueError(
            "nucleophile_identity='cys' contradicts motif_class='AHSMG': "
            "AHSMG is a Ser motif class, not a Cys nucleophile motif"
        )
    lipase_box = _text(row, "lipase_box_or_ahsmg")
    if identity == "cys" and "ahsmg" in _key(lipase_box):
        raise ValueError(
            f"nucleophile_identity='cys' contradicts lipase_box_or_ahsmg={lipase_box!r}: "
            "AHSMG is a Ser motif class, so an AHSMG reference cannot be recorded as Cys"
        )

    localization, localization_source = _resolve_localization(row)
    grade = _resolve(
        _text(row, "experimental_evidence_grade"), "experimental_evidence_grade",
        _EVIDENCE_GRADE_BY_KEY, EVIDENCE_GRADES, default="",
    )
    discovery_eligible = _resolve_optional_flag(row, "discovery_training_eligible")
    functional_eligible = _resolve_optional_flag(row, "functional_calibration_eligible")
    if functional_eligible == "true" and grade not in {"E2", "E3"}:
        raise ValueError(
            "functional_calibration_eligible='true' requires experimental_evidence_grade in "
            f"['E2', 'E3']; found experimental_evidence_grade={grade!r} "
            "(E1 and A cannot be counted as direct functional positives)"
        )
    catalytic_domain_type = _resolve(
        _text(row, "catalytic_domain_type"), "catalytic_domain_type",
        _CATALYTIC_DOMAIN_TYPE_BY_KEY, CATALYTIC_DOMAIN_TYPES,
        default="unresolved", aliases=CATALYTIC_DOMAIN_TYPE_ALIASES,
    )
    return {
        "reference_id": _text(row, "reference_id"),
        "model_layer": layer,
        "functional_calibration_status": status,
        "experimental_evidence_grade": grade,
        "discovery_training_eligible": discovery_eligible,
        "functional_calibration_eligible": functional_eligible,
        "catalytic_domain_type": catalytic_domain_type,
        "localization_evidence": localization,
        "localization_evidence_source": localization_source,
        "nucleophile_identity": identity,
        "motif_class": motif_class,
    }


def validate_candidate_evidence(row: dict) -> dict:
    """Validate the evidence fields of one candidate row.

    Rejections enforced here (each ``ValueError`` names the offending field and
    value):

    * ``nucleophile_identity=ahsmg`` and any other value outside
      ``NUCLEOPHILE_IDENTITIES``; a Cys identity combined with an AHSMG motif
      class is contradictory;
    * ``localization_evidence=experimental_confirmed`` without a stated
      non-SignalP ``localization_evidence_source``;
    * a ``discovery_hmm_uncalibrated`` row that carries a family call
      (``sequence_family_call``/``family_call``) or claims a
      ``core_sequence_homolog`` disposition;
    * an ``excluded_input_quality`` disposition without a documented
      input-quality reason (the discovery layer records problems but never
      filters, deletes or demotes).
    """
    identity = _resolve_nucleophile_identity(
        _text(row, "nucleophile_identity"), "nucleophile_identity", default="unresolved",
    )
    motif_class = _resolve(
        _text(row, "motif_class"), "motif_class", _MOTIF_CLASS_BY_KEY, MOTIF_CLASSES,
        default="unresolved", aliases=MOTIF_CLASS_ALIASES,
    )
    if motif_class == "AHSMG" and identity == "cys":
        raise ValueError(
            "nucleophile_identity='cys' contradicts motif_class='AHSMG': "
            "AHSMG is a Ser motif class, not a Cys nucleophile motif"
        )

    layer = _resolve(
        _text(row, "model_layer"), "model_layer", _MODEL_LAYER_BY_KEY, MODEL_LAYERS,
        default="reference_query_only",
    )
    status = _resolve(
        _text(row, "functional_calibration_status"), "functional_calibration_status",
        _FUNCTIONAL_STATUS_BY_KEY, FUNCTIONAL_CALIBRATION_STATUSES,
        default="not_function_calibrated",
    )
    validate_model_state(layer, status)

    transport_signal_prediction = _resolve(
        _text(row, "transport_signal_prediction"), "transport_signal_prediction",
        _TRANSPORT_SIGNAL_BY_KEY, TRANSPORT_SIGNAL_PREDICTIONS,
        default="not_tested", aliases=TRANSPORT_SIGNAL_ALIASES,
    )
    localization, localization_source = _resolve_localization(row)
    grade = _resolve(
        _text(row, "experimental_evidence_grade"), "experimental_evidence_grade",
        _EVIDENCE_GRADE_BY_KEY, EVIDENCE_GRADES, default="",
    )
    catalytic_domain_type = _resolve(
        _text(row, "catalytic_domain_type"), "catalytic_domain_type",
        _CATALYTIC_DOMAIN_TYPE_BY_KEY, CATALYTIC_DOMAIN_TYPES,
        default="unresolved", aliases=CATALYTIC_DOMAIN_TYPE_ALIASES,
    )
    disposition = _resolve(
        _text(row, "primary_disposition"), "primary_disposition",
        _PRIMARY_DISPOSITION_BY_KEY, PRIMARY_DISPOSITIONS, default="",
    )
    family_call_field, family_call = _resolve_family_call(row)

    input_quality_reason = ""
    if disposition == "excluded_input_quality":
        input_quality_reason = _resolve_input_quality_reason(row)
        if not input_quality_reason:
            detail = ""
            if layer == "discovery_hmm_uncalibrated":
                detail = (
                    "; the discovery layer only scores existing intermediate products and never "
                    "filters, deletes or demotes candidates, so an unexplained exclusion is invalid"
                )
            raise ValueError(
                "primary_disposition='excluded_input_quality' requires a documented input-quality "
                f"reason in one of {list(INPUT_QUALITY_REASON_FIELDS)}{detail}"
            )

    if layer == "discovery_hmm_uncalibrated":
        if family_call:
            raise ValueError(
                f"{family_call_field}={family_call!r} is forbidden with "
                "model_layer='discovery_hmm_uncalibrated': the discovery layer only scores existing "
                "intermediate products and never makes a family call"
            )
        if disposition == "core_sequence_homolog":
            raise ValueError(
                "primary_disposition='core_sequence_homolog' is forbidden with "
                "model_layer='discovery_hmm_uncalibrated': a core sequence homolog requires an "
                "independently validated sequence family model"
            )

    return {
        "model_layer": layer,
        "functional_calibration_status": status,
        "primary_disposition": disposition,
        "sequence_family_call": family_call,
        "experimental_evidence_grade": grade,
        "catalytic_domain_type": catalytic_domain_type,
        "transport_signal_prediction": transport_signal_prediction,
        "localization_evidence": localization,
        "localization_evidence_source": localization_source,
        "nucleophile_identity": identity,
        "motif_class": motif_class,
        "evidence_flags": _text(row, "evidence_flags"),
        "input_quality_reason": input_quality_reason,
    }


def load_sequence_model_gates(path: "str | Path") -> "dict[str, str]":
    """Load and validate the preregistered sequence-model gate row.

    Raises :class:`ValueError` naming the offending column or value when the
    column set, the single gate row, the gate version, a numeric bound or a
    boolean flag is malformed.
    """
    gates_path = Path(path)
    if not gates_path.is_file() or gates_path.is_symlink():
        raise ValueError(f"sequence model gates is not a regular file: {gates_path}")
    with gates_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    if len(fields) != len(set(fields)):
        raise ValueError(f"sequence model gates has duplicate columns: {fields}")
    missing = sorted(set(SEQUENCE_MODEL_GATE_COLUMNS) - set(fields))
    extra = sorted(set(fields) - set(SEQUENCE_MODEL_GATE_COLUMNS))
    if missing or extra:
        raise ValueError(
            "sequence model gates columns do not match the preregistered set; "
            f"missing columns={missing}, unexpected columns={extra}"
        )
    if len(rows) != 1:
        raise ValueError(
            f"sequence model gates must contain exactly one gate row; found {len(rows)}"
        )
    row = {column: _text(rows[0], column) for column in SEQUENCE_MODEL_GATE_COLUMNS}
    if row["gate_version"] != SEQUENCE_MODEL_GATE_VERSION:
        raise ValueError(
            f"gate_version={row['gate_version']!r} is not the preregistered gate version "
            f"{SEQUENCE_MODEL_GATE_VERSION!r}"
        )
    for column in SEQUENCE_MODEL_GATE_INTEGER_COLUMNS:
        value = row[column]
        try:
            number = int(value)
        except ValueError:
            raise ValueError(f"{column}={value!r} must be an integer") from None
        if number < 0:
            raise ValueError(f"{column}={value!r} must not be negative")
        if column in SEQUENCE_MODEL_GATE_POSITIVE_COLUMNS and number < 1:
            raise ValueError(f"{column}={value!r} must be at least 1")
    recall_text = row["min_heldout_recall"]
    try:
        recall = float(recall_text)
    except ValueError:
        raise ValueError(f"min_heldout_recall={recall_text!r} must be a number") from None
    if not 0.0 < recall <= 1.0:
        raise ValueError(f"min_heldout_recall={recall_text!r} must be greater than 0 and at most 1")
    for column in SEQUENCE_MODEL_GATE_BOOLEAN_COLUMNS:
        if row[column] not in {"true", "false"}:
            raise ValueError(f"{column}={row[column]!r} must be 'true' or 'false'")
    return row


def main(argv: Sequence[str] | None = None) -> int:
    """Validate a TSV of candidate or reference rows against the shared ontology."""
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("tsv", help="TSV file of candidate or reference rows")
    parser.add_argument(
        "--kind", choices=("candidate", "reference"), default="candidate",
        help="row family to validate (default: candidate)",
    )
    parser.add_argument(
        "--keep-going", action="store_true",
        help="report every violation instead of stopping at the first one",
    )
    parser.add_argument(
        "--gates", help="also validate a phaded_sequence_model_gates.tsv file",
    )
    args = parser.parse_args(argv)

    if args.gates:
        try:
            load_sequence_model_gates(args.gates)
        except ValueError as error:
            print(f"invalid sequence model gates: {error}", file=sys.stderr)
            return 1

    path = Path(args.tsv)
    if not path.is_file() or path.is_symlink():
        print(f"input is not a regular file: {path}", file=sys.stderr)
        return 2
    validator = validate_candidate_evidence if args.kind == "candidate" else validate_reference_evidence
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            print(f"input has no TSV header: {path}", file=sys.stderr)
            return 2
        rows = list(reader)

    invalid = 0
    checked = 0
    for index, row in enumerate(rows, start=2):
        checked += 1
        try:
            validator(row)
        except ValueError as error:
            invalid += 1
            print(f"{path}:{index}: {error}", file=sys.stderr)
            if not args.keep_going:
                break
    print(json.dumps({"rows": len(rows), "valid": checked - invalid, "invalid": invalid}))
    return 1 if invalid else 0


if __name__ == "__main__":
    raise SystemExit(main())
