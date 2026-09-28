#!/usr/bin/env python3
"""Independently validate PhaDED sequence-family models (Task 6).

This module answers exactly one question: *did this sequence family model earn
the ``sequence_family_hmm_validated`` layer?*  It reads the profile manifest
emitted by ``build_phaded_profiles`` plus the held-out / nearest-family /
competing-fold evaluation inputs, applies the preregistered gates in
``pipeline/config/phaded_sequence_model_gates.tsv``, and writes

* ``sequence_model_evaluation.tsv`` -- one row per evaluated model,
* ``sequence_model_gate_detail.tsv`` -- one row per evaluated model with the
  per-leg gate detail and the stable blocking reason codes,
* ``sequence_model_evaluation_manifest.json`` -- provenance for both.

Boundaries (each one is load-bearing)
-------------------------------------
* **An experimental positive is not a precondition.**  Experimental counts may
  never create, replace or rescue the sequence layer: a family with
  ``experimental_anchor_count=0`` that passes every sequence leg is validated,
  and a family with three real anchors that fails a sequence leg is not
  (approved design assumption 1).
* **This module never calibrates function.**  ``functional_calibration_status``
  is copied from the input row verbatim -- including an unrecognized value --
  and is never upgraded, downgraded or inferred.  No code path returns the
  ``calibrated_candidate_model`` layer; that status exists only on a separate,
  explicitly authorized action.  Every output row additionally carries the
  fixed boundary pair ``sequence_claim`` / ``phenotype_boundary``.
* **Absent values are ``pending`` and fail closed.**  A gate leg whose input is
  missing, empty or unparseable reports ``pending`` and blocks; nothing is
  assumed to have succeeded.
* **``EFFN`` is a diagnostic.**  ``hmmbuild`` EFFN below one is common for small
  families (for example the 7-sequence hfam_70 profile reports EFFN 0.748535)
  and is reported here as a diagnostic field only; it never blocks.
* **Held-out panels smaller than five must recover every sequence.**  For
  panels of five or more the configured ``min_heldout_recall`` applies.  The
  rule is encoded explicitly and reported as ``heldout_rule``; it is never
  inferred from rounding of ``0.8``.
* **Bit reproducibility is not a gate.**  MAFFT 7.525 is not bit-reproducible
  (the same input can yield different alignments), so ``bit_reproducible=false``
  is recorded and the concrete ``alignment_sha256`` is required instead.

This module never writes outside its output directory, never touches
``results/``, ``runs/``, ``deploy/`` or ``pipeline/config/``, and reads no
network or third-party dependency (standard library only).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


def _load_schema() -> Any:
    """Import the shared evidence ontology as a module.

    ``phaded_evidence_schema`` lives next to this file but the scripts
    directory is not a package, so a plain import only works when it is already
    on ``sys.path``.  Every other import path is a hard error.
    """
    try:
        import phaded_evidence_schema as schema  # type: ignore[import-not-found]

        return schema
    except ImportError:
        pass
    schema_path = Path(__file__).resolve().parent / "phaded_evidence_schema.py"
    if not schema_path.is_file():
        raise ImportError(f"cannot locate the shared evidence ontology next to {__file__}")
    spec = importlib.util.spec_from_file_location("phaded_evidence_schema", schema_path)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise ImportError(f"cannot load the shared evidence ontology: {schema_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("phaded_evidence_schema", module)
    spec.loader.exec_module(module)
    return module


SCHEMA = _load_schema()

MODEL_LAYER = "sequence_family_hmm_validated"
#: The only layer this evaluator may emit, besides passing an input row through.
#: The layer that exists solely on an explicitly authorized action is
#: deliberately absent and is pinned by a test.
PROMOTABLE_MODEL_LAYERS = frozenset({"discovery_hmm_uncalibrated", "sequence_family_hmm_validated"})
#: Layers that are read, recorded and blocked -- never promoted.
UNPROMOTABLE_MODEL_LAYERS = frozenset({"reference_query_only"})
MODEL_LAYERS = PROMOTABLE_MODEL_LAYERS | UNPROMOTABLE_MODEL_LAYERS

#: Fixed overclaiming boundary carried by every evaluated model.
SEQUENCE_CLAIM = "sequence-defined family homology"
PHENOTYPE_BOUNDARY = "does not establish PHB/PHA degradation"
#: File naming convention for an absent value (AGENTS.md: never invent a value).
PENDING = "pending"
GATE_VERSION = SCHEMA.SEQUENCE_MODEL_GATE_VERSION
SMALL_HELDOUT_PANEL_THRESHOLD = SCHEMA.SMALL_HELDOUT_PANEL_THRESHOLD
MIN_TRAINING_SEQUENCES = 3
UNIQUE_BEST_MODEL_COLUMN = "unique_best_model"
DEFAULT_GATES_PATH = Path(__file__).resolve().parents[1] / "config" / "phaded_sequence_model_gates.tsv"
VALIDATED_LAYER_MARGIN = "nearest_family_score_exceeds_competing_fold"

#: Reason codes in the order the gate legs are evaluated.  Every code is stable
#: and machine-readable; renaming one is a breaking change for downstream
#: consumers of ``sequence_model_gate_detail.tsv``.
REASON_MODEL_LAYER_UNPROMOTABLE = "model_layer_unpromotable"
REASON_MODEL_LAYER_UNKNOWN = "model_layer_unknown"
REASON_MODEL_LAYER_REFERENCE_ONLY = "model_layer_reference_query_only_unpromotable"
REASON_INPUT_STATUS_NOT_ALLOWED = "input_functional_calibration_status_not_allowed"
REASON_DISCOVERY_FAMILY_CALL = "discovery_family_call_present"
REASON_TRAINING_BELOW_MINIMUM = "training_sequences_below_minimum"
REASON_TRAINING_COUNT_PENDING = "training_sequence_count_pending"
REASON_TRAINING_SET_HASH_INVALID = "training_set_sha256_invalid"
REASON_HELDOUT_PANEL_BELOW_MINIMUM = "heldout_panel_below_minimum"
REASON_HELDOUT_COUNT_PENDING = "heldout_count_pending"
REASON_HELDOUT_RECOVERED_PENDING = "heldout_recovered_pending"
REASON_HELDOUT_RECALL_BELOW_MINIMUM = "heldout_recall_below_minimum"
REASON_HELDOUT_RECOVERED_EXCEEDS_COUNT = "heldout_recovered_exceeds_count"
REASON_CONFOUNDER_HITS_PRESENT = "unexplained_confounder_hits_present"
REASON_CONFOUNDER_HITS_PENDING = "unexplained_confounder_hits_pending"
REASON_NOT_UNIQUE_BEST_MODEL = "not_unique_best_model"
REASON_UNIQUE_BEST_MODEL_PENDING = "unique_best_model_pending"
REASON_ALIGNMENT_HASH_PENDING = "alignment_sha256_pending"
REASON_ALIGNMENT_HASH_INVALID = "alignment_sha256_invalid"
REASON_MARGIN_UNRESOLVED = VALIDATED_LAYER_MARGIN

BLOCKING_REASON_CODES = (
    REASON_MODEL_LAYER_UNPROMOTABLE,
    REASON_MODEL_LAYER_UNKNOWN,
    REASON_MODEL_LAYER_REFERENCE_ONLY,
    REASON_INPUT_STATUS_NOT_ALLOWED,
    REASON_DISCOVERY_FAMILY_CALL,
    REASON_TRAINING_BELOW_MINIMUM,
    REASON_TRAINING_COUNT_PENDING,
    REASON_TRAINING_SET_HASH_INVALID,
    REASON_HELDOUT_PANEL_BELOW_MINIMUM,
    REASON_HELDOUT_COUNT_PENDING,
    REASON_HELDOUT_RECOVERED_PENDING,
    REASON_HELDOUT_RECALL_BELOW_MINIMUM,
    REASON_HELDOUT_RECOVERED_EXCEEDS_COUNT,
    REASON_CONFOUNDER_HITS_PRESENT,
    REASON_CONFOUNDER_HITS_PENDING,
    REASON_NOT_UNIQUE_BEST_MODEL,
    REASON_UNIQUE_BEST_MODEL_PENDING,
    REASON_ALIGNMENT_HASH_PENDING,
    REASON_ALIGNMENT_HASH_INVALID,
    REASON_MARGIN_UNRESOLVED,
)

HELDOUT_RULE_ALL_RECOVERED = "all_recovered_below_five"
HELDOUT_RULE_MIN_RECALL = "min_recall_at_or_above_five"

#: Evaluation inputs the CLI accepts next to (or instead of) the manifest row.
METRIC_COLUMNS = (
    "heldout_count",
    "heldout_recovered",
    "unexplained_confounder_hits",
    UNIQUE_BEST_MODEL_COLUMN,
    "eff_nseq",
    "coverage_fraction",
    "family_assignment_stable",
    "nearest_family_score",
    "competing_fold_score",
)
#: Aliases accepted for the evaluation inputs so a hand-built input TSV may use
#: either spelling.  ``unique_best_model`` is the pinned name used in outputs.
METRIC_ALIASES = {
    "n_heldout": "heldout_count",
    "heldout_n": "heldout_count",
    "heldout_total": "heldout_count",
    "n_heldout_recovered": "heldout_recovered",
    "heldout_hits": "heldout_recovered",
    "heldout_recovered_count": "heldout_recovered",
    "confounder_hits": "unexplained_confounder_hits",
    "unexplained_confounder_count": "unexplained_confounder_hits",
    "unique_best": UNIQUE_BEST_MODEL_COLUMN,
    "eff_nseq_diagnostic": "eff_nseq",
    "coverage": "coverage_fraction",
    "assignment_stability": "family_assignment_stable",
    "nearest_family_best_score": "nearest_family_score",
    "competing_fold_best_score": "competing_fold_score",
}
#: Manifest identity: an evaluation-input TSV may never overwrite these.
IDENTITY_COLUMNS = (
    "profile_id",
    "profile_kind",
    "model_layer",
    "functional_calibration_status",
)

EVALUATION_FIELDS = [
    "profile_id",
    "profile_kind",
    "input_model_layer",
    "model_layer",
    "promoted",
    "sequence_claim",
    "phenotype_boundary",
    "functional_calibration_status",
    "gate_version",
    "training_sequence_count",
    "experimental_anchor_count",
    "annotation_only_count",
    "heldout_count",
    "heldout_recovered",
    "heldout_recall",
    "heldout_rule",
    "unexplained_confounder_hits",
    "unique_best_model",
    "score_margin",
    "coverage_fraction",
    "family_assignment_stable",
    "eff_nseq",
    "eff_nseq_diagnostic",
    "alignment_sha256",
    "hmm_sha256",
    "training_set_sha256",
    "bit_reproducible",
    "blocking_reason_count",
    "blocking_reasons",
]
GATE_DETAIL_FIELDS = [
    "profile_id",
    "model_layer",
    "promoted",
    "gate_version",
    "min_training_sequences",
    "min_heldout_sequences",
    "min_heldout_recall",
    "max_unexplained_confounder_hits",
    "require_unique_best_model",
    "require_alignment_hash",
    "heldout_rule",
    "small_heldout_panel_threshold",
    "heldout_recall",
    "heldout_required_recovered",
    "score_margin",
    "unexplained_confounder_hits",
    "unique_best_model",
    "alignment_sha256_status",
    "eff_nseq",
    "blocking_reasons",
]

_HEX64 = re.compile(r"\A[0-9a-fA-F]{64}\Z")
#: The functional calibration statuses this evaluator may ever report: exactly
#: the statuses the shared ontology allows on a sequence-family layer.  The
#: separately authorized calibrated status is therefore excluded by construction
#: and is never named here.  A status outside this set is copied verbatim and
#: blocks promotion (see :func:`_normalize_status`).
_ALLOWED_FUNCTIONAL_STATUSES = frozenset(
    SCHEMA.ALLOWED_MODEL_STATES["sequence_family_hmm_validated"]
)
_STATUS_BY_KEY = {
    SCHEMA._key(status): status for status in sorted(_ALLOWED_FUNCTIONAL_STATUSES)
}
_LAYER_BY_KEY = {SCHEMA._key(layer): layer for layer in sorted(MODEL_LAYERS)}


# ---------------------------------------------------------------------------
# Row helpers
# ---------------------------------------------------------------------------

def sha256_file(path: "str | Path") -> str:
    """Return the SHA-256 of a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text(row: Mapping[str, object], field: str) -> str:
    """Read a field defensively; an absent or ``None`` column reads as ``""``."""
    if field not in row:
        return ""
    return _text_value(row[field])


def _text_value(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _raw_field(row: Mapping[str, object], names: Sequence[str]) -> "tuple[str, str]":
    """Return ``(column_name, value)`` for the first present, non-empty field."""
    for name in names:
        if name in row:
            text = _text(row, name)
            if text:
                return name, text
    return "", ""


def _is_pending(text: str) -> bool:
    return not text or SCHEMA._key(text) in {"pending", "not_available", "unknown", "na"}


def _as_int(text: str) -> "int | None":
    """Strict integer reader; anything but a plain int is ``None`` (fail closed)."""
    text = text.strip()
    if not text or not re.fullmatch(r"[+-]?\d+", text):
        return None
    return int(text)


def _as_flag(text: str) -> "bool | None":
    """Strict boolean reader for gate flags; unknown values are ``None``."""
    key = SCHEMA._key(text)
    if key in SCHEMA._TRUE_TOKENS:
        return True
    if key in SCHEMA._FALSE_TOKENS:
        return False
    return None


def _as_float(text: str) -> "float | None":
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _format_number(value: "float | None", digits: int = 6) -> str:
    """Render a numeric metric compactly; absent values stay ``pending``."""
    if value is None:
        return PENDING
    return (("%." + str(digits) + "f") % value).rstrip("0").rstrip(".")


def _cell(value: object) -> str:
    """Render a metric cell; absent values stay ``pending``."""
    if value is None:
        return PENDING
    return str(value)


def _flag_cell(value: "bool | None") -> str:
    if value is None:
        return PENDING
    return "true" if value else "false"


def _normalize_model_layer(row: Mapping[str, object]) -> "tuple[str, bool]":
    """Return ``(layer, recognized)`` without ever inventing a layer."""
    text = _text(row, "model_layer")
    return _LAYER_BY_KEY.get(SCHEMA._key(text), ""), SCHEMA._key(text) in _LAYER_BY_KEY


def _normalize_status(row: Mapping[str, object]) -> "tuple[str, bool]":
    """Copy the functional calibration status verbatim; never upgrade it.

    The mapping contains only the non-calibrated statuses, so this function has
    no expression that could produce a calibrated value.  An absent column
    reads as ``not_function_calibrated`` and anything unrecognized is returned
    unchanged and flagged, which blocks promotion instead of being reinterpreted.
    """
    text = _text(row, "functional_calibration_status")
    if not text:
        return "not_function_calibrated", True
    allowed = _STATUS_BY_KEY.get(SCHEMA._key(text))
    if allowed is not None:
        return allowed, True
    return text, False


def _heldout_counts(row: Mapping[str, object]) -> "tuple[int | None, int | None]":
    _, raw_count = _raw_field(row, ("heldout_count", "n_heldout", "heldout_total"))
    _, raw_recovered = _raw_field(
        row, ("heldout_recovered", "n_heldout_recovered", "heldout_hits")
    )
    count = None if _is_pending(raw_count) else _as_int(raw_count)
    recovered = None if _is_pending(raw_recovered) else _as_int(raw_recovered)
    return count, recovered


def _confounder_count(row: Mapping[str, object]) -> "int | None":
    _, raw = _raw_field(
        row,
        (
            "unexplained_confounder_hits",
            "confounder_hits",
            "unexplained_confounder_count",
        ),
    )
    if _is_pending(raw):
        return None
    return _as_int(raw)


def _unique_best_model(row: Mapping[str, object]) -> "bool | None":
    _, raw = _raw_field(row, (UNIQUE_BEST_MODEL_COLUMN, "unique_best"))
    if _is_pending(raw):
        return None
    return _as_flag(raw)


def _alignment_hash(row: Mapping[str, object]) -> "tuple[str, str]":
    """Return ``(status, value)`` where status is verified/pending/invalid."""
    text = _text(row, "alignment_sha256")
    if _is_pending(text):
        return "pending", ""
    if not _HEX64.match(text):
        return "invalid", ""
    return "verified", text.lower()


# ---------------------------------------------------------------------------
# Gates
# ---------------------------------------------------------------------------

def _default_gates() -> dict[str, str]:
    return {
        "gate_version": GATE_VERSION,
        "min_training_sequences": str(MIN_TRAINING_SEQUENCES),
        "min_heldout_sequences": "1",
        "min_heldout_recall": "0.8",
        "max_unexplained_confounder_hits": "0",
        "require_unique_best_model": "true",
        "require_alignment_hash": "true",
    }


def _typed_gates(gates: "Mapping[str, object] | None") -> dict[str, Any]:
    """Read the seven preregistered gate columns defensively.

    Every column is optional and falls back to the preregistered default, so
    the evaluator stays usable with a partial gate mapping while validation of
    the on-disk file remains :func:`phaded_evidence_schema.load_sequence_model_gates`.
    """
    source: dict[str, object] = {}
    for key, value in (gates or {}).items():
        if isinstance(key, str):
            source[key] = value
    defaults = _default_gates()
    version = _text_value(source.get("gate_version")) or defaults["gate_version"]

    def _int_column(column: str) -> int:
        parsed = _as_int(_text_value(source.get(column)))
        if parsed is None or parsed < 0:
            return int(defaults[column])
        return parsed

    recall = _as_float(_text_value(source.get("min_heldout_recall")))
    if recall is None or not 0.0 < recall <= 1.0:
        recall = float(defaults["min_heldout_recall"])
    return {
        "gate_version": version,
        "min_training_sequences": _int_column("min_training_sequences"),
        "min_heldout_sequences": _int_column("min_heldout_sequences"),
        "min_heldout_recall": recall,
        "max_unexplained_confounder_hits": _int_column("max_unexplained_confounder_hits"),
        "require_unique_best_model": _as_flag(
            _text_value(source.get("require_unique_best_model"))
        )
        is not False,
        "require_alignment_hash": _as_flag(
            _text_value(source.get("require_alignment_hash"))
        )
        is not False,
    }


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def _heldout_metrics(
    count: "int | None", recovered: "int | None", gates: Mapping[str, Any]
) -> dict[str, Any]:
    """Held-out recall, required recoveries and the explicit panel rule."""
    recall = float(recovered) / float(count) if count else None
    small_panel = count is not None and count < SMALL_HELDOUT_PANEL_THRESHOLD
    rule = HELDOUT_RULE_ALL_RECOVERED if small_panel else HELDOUT_RULE_MIN_RECALL
    if count is None:
        required: "int | None" = None
    elif small_panel:
        # "For held-out panels smaller than five, all held-out sequences must be
        # recovered": encoded here, never inferred from rounding of the recall.
        required = count
    else:
        required = int(-(-count * gates["min_heldout_recall"] // 1))
    return {
        "heldout_recall": recall,
        "heldout_rule": rule,
        "small_panel": small_panel,
        "required_recovered": required,
    }


def compute_metrics(row: Mapping[str, object], gates: "Mapping[str, object] | None" = None) -> dict[str, Any]:
    """Every reported metric for one model, with ``pending`` for absent values.

    Reported: held-out recall and the explicit held-out rule, the exact
    unexplained confounder hit count, unique-best-model status, score margin
    (nearest-family minus competing-fold, so positive means the family wins),
    coverage, family assignment stability and ``EFFN`` as a **diagnostic**.
    """
    gate_values = _typed_gates(gates)
    count, recovered = _heldout_counts(row)
    heldout = _heldout_metrics(count, recovered, gate_values)

    _, margin_text = _raw_field(row, ("score_margin",))
    margin = None if _is_pending(margin_text) else _as_float(margin_text)
    _, nearest_text = _raw_field(row, ("nearest_family_score", "nearest_family_best_score"))
    _, competing_text = _raw_field(row, ("competing_fold_score", "competing_fold_best_score"))
    nearest = None if _is_pending(nearest_text) else _as_float(nearest_text)
    competing = None if _is_pending(competing_text) else _as_float(competing_text)
    margin_source = "reported"
    if margin is None and nearest is not None and competing is not None:
        margin = nearest - competing
        margin_source = "computed_from_nearest_and_competing"

    training = _as_int(_text(row, "training_sequence_count"))
    experimental = _as_int(_text(row, "experimental_anchor_count"))
    annotation_only = _as_int(_text(row, "annotation_only_count"))
    eff_nseq = _as_float(_text(row, "eff_nseq"))
    coverage = _as_float(_text(row, "coverage_fraction"))
    stability = _text(row, "family_assignment_stable")

    return {
        # Counts stay numeric (``None`` when absent); the TSV/JSON writers
        # render ``None`` as ``pending`` so a missing count is never a ``0``.
        "training_sequence_count": training,
        "experimental_anchor_count": experimental,
        "annotation_only_count": annotation_only,
        "heldout_count": count,
        "heldout_recovered": recovered,
        "heldout_recall": heldout["heldout_recall"],
        "heldout_rule": heldout["heldout_rule"],
        "heldout_required_recovered": heldout["required_recovered"],
        "small_heldout_panel_threshold": SMALL_HELDOUT_PANEL_THRESHOLD,
        "unexplained_confounder_hits": _confounder_count(row),
        "unique_best_model": _unique_best_model(row),
        # Rendered here so the metric is identical everywhere it is reported.
        "score_margin": _format_number(margin),
        "score_margin_source": margin_source,
        "nearest_family_score": nearest,
        "competing_fold_score": competing,
        "coverage_fraction": coverage,
        "family_assignment_stable": stability or PENDING,
        # Diagnostics: reported, never a gate leg.
        "eff_nseq": _format_number(eff_nseq),
        "eff_nseq_diagnostic": _format_number(eff_nseq),
        "bit_reproducible": _text(row, "bit_reproducible") or PENDING,
    }


# ---------------------------------------------------------------------------
# Gate evaluation
# ---------------------------------------------------------------------------

def blocking_reasons(
    row: Mapping[str, object],
    gates: "Mapping[str, object] | None" = None,
    detail: "dict[str, Any] | None" = None,
) -> list[str]:
    """Return the stable blocking reason codes for one manifest row.

    One code per failed gate leg, in evaluation order; an empty list means the
    row earned ``sequence_family_hmm_validated``.  Experimental and annotation
    counts and the input functional calibration status never appear here as a
    reason *for* promotion: an unrecognized calibration status is reported (it
    forbids calling the model validated) but is still copied verbatim.
    """
    gate_values = _typed_gates(gates)
    metrics = compute_metrics(row, gates)
    layer, recognized = _normalize_model_layer(row)
    status, status_allowed = _normalize_status(row)
    if detail is not None:
        detail.update(
            {
                "gate_version": gate_values["gate_version"],
                "min_training_sequences": gate_values["min_training_sequences"],
                "min_heldout_sequences": gate_values["min_heldout_sequences"],
                "min_heldout_recall": gate_values["min_heldout_recall"],
                "max_unexplained_confounder_hits": gate_values[
                    "max_unexplained_confounder_hits"
                ],
                "require_unique_best_model": gate_values["require_unique_best_model"],
                "require_alignment_hash": gate_values["require_alignment_hash"],
                "heldout_rule": metrics["heldout_rule"],
                "small_heldout_panel_threshold": SMALL_HELDOUT_PANEL_THRESHOLD,
                "heldout_recall": _format_number(metrics["heldout_recall"]),
                "heldout_required_recovered": metrics["heldout_required_recovered"],
                "score_margin": metrics["score_margin"],
                "unexplained_confounder_hits": metrics["unexplained_confounder_hits"],
                "unique_best_model": metrics["unique_best_model"],
                "eff_nseq": metrics["eff_nseq_diagnostic"],
                "input_model_layer": layer,
                "functional_calibration_status": status,
            }
        )

    reasons: list[str] = []

    # -- layer bookkeeping -------------------------------------------------
    if not recognized:
        reasons.append(REASON_MODEL_LAYER_UNKNOWN)
    if layer not in PROMOTABLE_MODEL_LAYERS:
        reasons.append(REASON_MODEL_LAYER_UNPROMOTABLE)
    if layer == "reference_query_only":
        reasons.append(REASON_MODEL_LAYER_REFERENCE_ONLY)
    if not status_allowed:
        reasons.append(REASON_INPUT_STATUS_NOT_ALLOWED)
    if layer == "discovery_hmm_uncalibrated":
        call_field, call = _raw_field(row, SCHEMA.FAMILY_CALL_FIELDS)
        if call:
            reasons.append(REASON_DISCOVERY_FAMILY_CALL)

    # -- training set ------------------------------------------------------
    training = metrics["training_sequence_count"]
    if training is None:
        reasons.append(REASON_TRAINING_COUNT_PENDING)
    elif training < gate_values["min_training_sequences"]:
        reasons.append(REASON_TRAINING_BELOW_MINIMUM)
    training_hash = _text(row, "training_set_sha256")
    if training_hash and not _is_pending(training_hash) and not _HEX64.match(training_hash):
        reasons.append(REASON_TRAINING_SET_HASH_INVALID)

    # -- held-out panel ----------------------------------------------------
    count = metrics["heldout_count"]
    recovered = metrics["heldout_recovered"]
    if count is None:
        reasons.append(REASON_HELDOUT_COUNT_PENDING)
    elif count < gate_values["min_heldout_sequences"]:
        reasons.append(REASON_HELDOUT_PANEL_BELOW_MINIMUM)
    if recovered is None:
        reasons.append(REASON_HELDOUT_RECOVERED_PENDING)
    if count is not None and recovered is not None:
        if recovered > count:
            reasons.append(REASON_HELDOUT_RECOVERED_EXCEEDS_COUNT)
        elif metrics["heldout_rule"] == HELDOUT_RULE_ALL_RECOVERED:
            # Fewer than five held-out sequences: every one must be recovered.
            if recovered != count:
                reasons.append(REASON_HELDOUT_RECALL_BELOW_MINIMUM)
        elif metrics["heldout_recall"] < gate_values["min_heldout_recall"]:
            reasons.append(REASON_HELDOUT_RECALL_BELOW_MINIMUM)

    # -- confounder confusion ---------------------------------------------
    confounders = metrics["unexplained_confounder_hits"]
    if confounders is None:
        reasons.append(REASON_CONFOUNDER_HITS_PENDING)
    elif confounders > gate_values["max_unexplained_confounder_hits"]:
        reasons.append(REASON_CONFOUNDER_HITS_PRESENT)

    # -- unique best model -------------------------------------------------
    unique = metrics["unique_best_model"]
    if gate_values["require_unique_best_model"]:
        if unique is None:
            reasons.append(REASON_UNIQUE_BEST_MODEL_PENDING)
        elif unique is False:
            reasons.append(REASON_NOT_UNIQUE_BEST_MODEL)

    # -- alignment provenance ---------------------------------------------
    # A present-but-malformed hash is always an error; the gate selects whether
    # an *absent* hash blocks.  ``bit_reproducible`` is never a gate leg.
    hash_status, _ = _alignment_hash(row)
    if hash_status == "invalid":
        reasons.append(REASON_ALIGNMENT_HASH_INVALID)
    elif hash_status == "pending" and gate_values["require_alignment_hash"]:
        reasons.append(REASON_ALIGNMENT_HASH_PENDING)

    # -- nearest-family / competing-fold margin (reported when available) ---
    margin = _as_float(metrics["score_margin"])
    if margin is not None and margin <= 0:
        reasons.append(REASON_MARGIN_UNRESOLVED)

    return reasons


def evaluate_model(
    row: Mapping[str, object], gates: "Mapping[str, object] | None" = None
) -> dict[str, Any]:
    """Evaluate one profile-manifest row against the preregistered sequence gates.

    Returns the layer the model earned (``sequence_family_hmm_validated`` only
    when :func:`blocking_reasons` is empty), the input layer as received, every
    metric of :func:`compute_metrics`, the blocking reason codes, and the fixed
    overclaiming boundary.  ``functional_calibration_status`` is copied from
    ``row`` verbatim; this function never returns the calibrated layer.
    """
    detail: dict[str, Any] = {}
    reasons = blocking_reasons(row, gates, detail=detail)
    metrics = compute_metrics(row, gates)
    input_layer, recognized = _normalize_model_layer(row)
    status, status_allowed = _normalize_status(row)
    earned = MODEL_LAYER if (recognized and status_allowed and not reasons) else input_layer
    return {
        "profile_id": _text(row, "profile_id") or PENDING,
        "profile_kind": _text(row, "profile_kind") or PENDING,
        "input_model_layer": input_layer or PENDING,
        # An unrecognized layer is echoed as received (it is not invented), so an
        # auditable cell never comes back empty.
        "model_layer": earned or _text(row, "model_layer") or PENDING,
        "promoted": not reasons,
        "blocking_reasons": reasons,
        "metrics": metrics,
        "detail": detail,
        # Overclaiming boundary: the sequence claim, the phenotype boundary and
        # the *input* functional calibration status, verbatim.
        "sequence_claim": SEQUENCE_CLAIM,
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
        "functional_calibration_status": status,
        "functional_calibration_changed": False,
    }


def evaluate_record(row: Mapping[str, object], gates: "Mapping[str, object] | None" = None) -> dict[str, Any]:
    """Evaluate one row and return both the ``str`` evaluation and gate-detail cells.

    Missing values stay ``pending``; no cell is inferred.  This is the shape the
    CLI writes to ``sequence_model_evaluation.tsv`` and
    ``sequence_model_gate_detail.tsv``.
    """
    result = evaluate_model(row, gates)
    metrics = result["metrics"]
    detail = result["detail"]
    hash_status, alignment_hash = _alignment_hash(row)
    evaluation = {
        "profile_id": result["profile_id"],
        "profile_kind": result["profile_kind"],
        "input_model_layer": result["input_model_layer"],
        "model_layer": result["model_layer"],
        "promoted": "true" if result["promoted"] else "false",
        "sequence_claim": result["sequence_claim"],
        "phenotype_boundary": result["phenotype_boundary"],
        "functional_calibration_status": result["functional_calibration_status"],
        "gate_version": detail["gate_version"],
        "training_sequence_count": _cell(metrics["training_sequence_count"]),
        "experimental_anchor_count": _cell(metrics["experimental_anchor_count"]),
        "annotation_only_count": _cell(metrics["annotation_only_count"]),
        "heldout_count": _cell(metrics["heldout_count"]),
        "heldout_recovered": _cell(metrics["heldout_recovered"]),
        "heldout_recall": _format_number(metrics["heldout_recall"]),
        "heldout_rule": metrics["heldout_rule"],
        "unexplained_confounder_hits": _cell(metrics["unexplained_confounder_hits"]),
        "unique_best_model": _flag_cell(metrics["unique_best_model"]),
        "score_margin": metrics["score_margin"],
        "coverage_fraction": _format_number(metrics["coverage_fraction"]),
        "family_assignment_stable": _cell(metrics["family_assignment_stable"]),
        "eff_nseq": metrics["eff_nseq_diagnostic"],
        "eff_nseq_diagnostic": metrics["eff_nseq_diagnostic"],
        "alignment_sha256": alignment_hash or PENDING,
        "hmm_sha256": _text(row, "hmm_sha256") or PENDING,
        "training_set_sha256": _text(row, "training_set_sha256") or PENDING,
        "bit_reproducible": _flag_cell(_as_flag(_text(row, "bit_reproducible"))),
        "blocking_reason_count": str(len(result["blocking_reasons"])),
        "blocking_reasons": ";".join(result["blocking_reasons"]),
    }
    gate_detail = {
        "profile_id": result["profile_id"],
        "model_layer": result["model_layer"],
        "promoted": evaluation["promoted"],
        "gate_version": detail["gate_version"],
        "min_training_sequences": str(detail["min_training_sequences"]),
        "min_heldout_sequences": str(detail["min_heldout_sequences"]),
        "min_heldout_recall": _format_number(detail["min_heldout_recall"]),
        "max_unexplained_confounder_hits": str(detail["max_unexplained_confounder_hits"]),
        "require_unique_best_model": _flag_cell(detail["require_unique_best_model"]),
        "require_alignment_hash": _flag_cell(detail["require_alignment_hash"]),
        "heldout_rule": detail["heldout_rule"],
        "small_heldout_panel_threshold": str(detail["small_heldout_panel_threshold"]),
        "heldout_recall": detail["heldout_recall"],
        "heldout_required_recovered": _cell(detail["heldout_required_recovered"]),
        "score_margin": detail["score_margin"],
        "unexplained_confounder_hits": _cell(detail["unexplained_confounder_hits"]),
        "unique_best_model": _flag_cell(detail["unique_best_model"]),
        "alignment_sha256_status": hash_status,
        "eff_nseq": detail["eff_nseq"],
        "blocking_reasons": evaluation["blocking_reasons"],
    }
    for field in GATE_DETAIL_FIELDS:
        gate_detail.setdefault(field, PENDING)
    return {"evaluation": evaluation, "gate_detail": gate_detail, "result": result}


# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------

def _read_tsv(path: "str | Path") -> "tuple[list[str], list[dict[str, str]]]":
    tsv = Path(path)
    if not tsv.is_file() or tsv.is_symlink():
        raise FileNotFoundError(f"input is not a regular file: {tsv}")
    with tsv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f"input has no TSV header: {tsv}")
        if len(fields) != len(set(fields)):
            raise ValueError(f"input has duplicate columns: {tsv}")
        return fields, list(reader)


def _read_evaluation_inputs(path: "str | Path") -> dict[str, dict[str, str]]:
    """Read the evaluation-input TSV, keyed by ``profile_id``.

    Only metric columns are accepted; the manifest identity columns
    (``profile_id``, ``profile_kind``, ``model_layer``,
    ``functional_calibration_status``) are refused here so an input file can
    never overwrite the layer or the functional calibration status.
    """
    _, rows = _read_tsv(path)
    metrics: dict[str, dict[str, str]] = {}
    for row in rows:
        profile_id = _text(row, "profile_id")
        if not profile_id:
            raise ValueError(f"evaluation input row has no profile_id: {path}")
        if profile_id in metrics:
            raise ValueError(f"evaluation input repeats profile_id={profile_id!r}: {path}")
        entry: dict[str, str] = {}
        for column in METRIC_COLUMNS:
            if column in row:
                entry[column] = _text(row, column)
        metrics[profile_id] = entry
    return metrics


# ---------------------------------------------------------------------------
# Evaluation run
# ---------------------------------------------------------------------------

def build_arg_namespace(argv: "list[str] | None" = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("profiles", help="profile manifest TSV produced by build_phaded_profiles")
    parser.add_argument("output_dir", help="new, empty output directory for this evaluation")
    parser.add_argument("--evaluation", help="optional TSV with held-out/confounder/challenge inputs")
    parser.add_argument(
        "--gates", default=str(DEFAULT_GATES_PATH),
        help="preregistered sequence-model gates TSV (default: %(default)s)",
    )
    parser.add_argument("--run-id", default="", help="run id recorded in the manifest (default: output dir name)")
    parser.add_argument("--tool", default=PENDING, help="tool/command provenance string")
    parser.add_argument("--command", default=PENDING, help="exact command provenance string")
    return parser.parse_args(argv)


def _bound_file(path: "str | Path") -> dict[str, Any]:
    if not Path(path).is_file():
        return {"path": str(path), "status": PENDING, "sha256": None, "size": None}
    return {
        "path": str(path),
        "status": "verified",
        "sha256": sha256_file(path),
        "size": Path(path).stat().st_size,
    }


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, PENDING)) for field in fields) + "\n")


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Evaluate every profile in the manifest and write the three outputs.

    Refuses to touch an existing output directory, validates that every input
    path exists, and keeps every absent value ``pending``.
    """
    profiles_path = Path(args.profiles)
    gates_path = Path(args.gates)
    output_dir = Path(args.output_dir)
    evaluation_path = Path(args.evaluation) if getattr(args, "evaluation", None) else None

    if not profiles_path.is_file() or profiles_path.is_symlink():
        raise FileNotFoundError(f"profile manifest is not a regular file: {profiles_path}")
    if not gates_path.is_file() or gates_path.is_symlink():
        raise FileNotFoundError(f"sequence model gates is not a regular file: {gates_path}")
    if evaluation_path is not None and (not evaluation_path.is_file() or evaluation_path.is_symlink()):
        raise FileNotFoundError(f"evaluation input is not a regular file: {evaluation_path}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(output_dir)

    # Validated by the shared ontology; a malformed gate row is a hard error.
    gates_row = SCHEMA.load_sequence_model_gates(gates_path)
    evaluation_inputs = (
        _read_evaluation_inputs(evaluation_path) if evaluation_path is not None else {}
    )
    _, manifest_rows = _read_tsv(profiles_path)

    evaluated: list[dict[str, Any]] = []
    evaluation_rows: list[dict[str, str]] = []
    detail_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for index, manifest_row in enumerate(manifest_rows):
        merged: dict[str, str] = {key: _text_value(value) for key, value in manifest_row.items()}
        profile_id = _text(merged, "profile_id") or ("row_%d" % (index + 1))
        merged["profile_id"] = profile_id
        if profile_id in seen:
            raise ValueError(f"profile manifest repeats profile_id={profile_id!r}")
        seen.add(profile_id)
        supplied = evaluation_inputs.get(profile_id, {})
        for column in METRIC_COLUMNS:
            if column in supplied:
                merged[column] = supplied[column]
        payload = evaluate_record(merged, gates_row)
        evaluated.append(payload["result"])
        evaluation_rows.append(payload["evaluation"])
        detail_rows.append(payload["gate_detail"])

    unknown_inputs = sorted(set(evaluation_inputs) - seen)
    if unknown_inputs:
        raise ValueError(
            "evaluation input names profiles that the manifest does not contain: "
            + ",".join(unknown_inputs)
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "sequence_model_evaluation.tsv", EVALUATION_FIELDS, evaluation_rows)
    _write_tsv(output_dir / "sequence_model_gate_detail.tsv", GATE_DETAIL_FIELDS, detail_rows)

    validated = [result for result in evaluated if result["promoted"]]
    counts = {
        "models_evaluated": len(evaluated),
        "models_validated": len(validated),
        "models_blocked": len(evaluated) - len(validated),
        "models_with_pending_inputs": sum(
            1 for row in evaluation_rows if PENDING in row.values()
        ),
    }
    manifest = {
        "schema_version": "1.0",
        "run_id": getattr(args, "run_id", "") or output_dir.name,
        "gate_version": gates_row["gate_version"],
        "model_layer": MODEL_LAYER,
        "unpromotable_layers": sorted(UNPROMOTABLE_MODEL_LAYERS),
        "sequence_claim": SEQUENCE_CLAIM,
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
        "functional_calibration_status_changed": False,
        "functional_calibration_status_upgraded": False,
        "family_call_made": False,
        "candidate_rows_changed": 0,
        "blocking_reason_codes": list(BLOCKING_REASON_CODES),
        "small_heldout_panel_threshold": SMALL_HELDOUT_PANEL_THRESHOLD,
        "heldout_rules": {
            "below_threshold": HELDOUT_RULE_ALL_RECOVERED,
            "at_or_above_threshold": HELDOUT_RULE_MIN_RECALL,
        },
        "eff_nseq_role": "diagnostic_only",
        "tool": getattr(args, "tool", PENDING),
        "command": getattr(args, "command", PENDING),
        "threads": {"max_single_task_threads_allowed": 40},
        "counts": counts,
        "inputs": {
            "profile_manifest": _bound_file(profiles_path),
            "sequence_model_gates": _bound_file(gates_path),
            "evaluation_inputs": (
                _bound_file(evaluation_path) if evaluation_path is not None else PENDING
            ),
            "gtdb": {
                "taxonomy": {"path": None, "status": PENDING, "sha256": None, "size": None},
                "metadata": {"path": None, "status": PENDING, "sha256": None, "size": None},
                "tree": {"path": None, "status": PENDING, "sha256": None, "size": None},
            },
        },
        "outputs": {
            "sequence_model_evaluation": str(output_dir / "sequence_model_evaluation.tsv"),
            "sequence_model_gate_detail": str(output_dir / "sequence_model_gate_detail.tsv"),
        },
    }
    (output_dir / "sequence_model_evaluation_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    return {
        "run_id": manifest["run_id"],
        "model_layer": MODEL_LAYER,
        "gate_version": gates_row["gate_version"],
        "counts": counts,
        "gates": {
            column: _text_value(gates_row.get(column))
            for column in SCHEMA.SEQUENCE_MODEL_GATE_COLUMNS
        },
        "rows": evaluated,
        "evaluation_rows": evaluation_rows,
        "gate_detail_rows": detail_rows,
        "manifest": manifest,
    }


def main(argv: "list[str] | None" = None) -> int:
    """Evaluate the profile manifest against the preregistered sequence gates."""
    args = build_arg_namespace(argv)
    try:
        payload = run(args)
    except (FileNotFoundError, FileExistsError, ValueError) as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
