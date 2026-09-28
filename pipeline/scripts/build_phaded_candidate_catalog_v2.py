#!/usr/bin/env python3
"""Build the v2 PhaDED candidate catalog with evidence-separated field semantics.

The frozen v1 pool (``filter_phaded_high_confidence.py``) collapses nucleophile
type, transport prediction, localization evidence, domain architecture and
confidence into one binary "high-confidence" label.  v2 keeps those claims in
separate columns and gives every protein exactly one ``primary_disposition``.

Semantic corrections implemented here (spec
``docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md``)
-----------------------------------------------------------------------
1. **AHSMG is a Ser motif class, never a third nucleophile identity.**  A row
   with ``ahsmg_state=supported`` yields
   ``nucleophile_identity=ser`` and ``motif_class=AHSMG``.  Identity is
   normalized by the shared ontology
   (:func:`phaded_evidence_schema.normalize_nucleophile`), so the two columns
   can never be collapsed back into one label.
2. **A SignalP result describes a predicted transport state only.**  It is
   recorded in ``transport_signal_prediction`` and is never promoted to
   ``localization_evidence``; a SignalP class alone can never produce
   ``localization_evidence=experimental_confirmed`` (the shared validator also
   rejects a SignalP-only ``localization_evidence_source``).
3. **SBD does not participate in the type 1 / type 2 definition.**
   ``catalytic_domain_type`` is emitted exactly as the domain-geometry evidence
   states it; SBD, lid, linker and the with-lipase variant are recorded
   independently in ``accessory_domain_architecture``.

Candidate catalog (six primary dispositions)
--------------------------------------------
``core_sequence_homolog`` / ``probable_sequence_homolog`` /
``remote_homolog_candidate`` / ``function_unresolved`` /
``deferred_structure_review`` / ``excluded_input_quality``.

``excluded_input_quality`` is reachable *only* through a documented
input-quality reason (missing sequence, obvious truncation, wrong input).  A row
whose only problem is scientific uncertainty -- ambiguous family, unresolved
localization, competing function, low confidence -- never receives it, and
scientific-uncertainty wording written into an input-quality field is rejected
outright (:data:`RULE_SCIENTIFIC_UNCERTAINTY_NOT_INPUT_ERROR`).  The single
adjudicated precedence case, "a row carries *both* a documented input-quality
reason and scientific-uncertainty cues", resolves to ``excluded_input_quality``
because a broken input cannot be classified at all -- but the uncertainty is
still recorded in ``evidence_flags`` and ``scientific_uncertainty_present``, so
the exclusion stays auditable, and the *reverse* direction (the same uncertainty
without the documented input-quality reason) never reaches exclusion.  Both
directions are pinned by tests.

Localization is a compatibility annotation
------------------------------------------
A transport prediction that contradicts the direction implied by the historical
superfamily adds ``localization_discordant`` to ``evidence_flags``.  It never by
itself deletes or demotes a sequence homolog.  The stricter "functional view" is
derived only (:func:`functional_view_eligible` / :func:`build_functional_view`)
and can never change :func:`primary_disposition`; the master classifier does not
read the view at all.

The 52 Cys-export demotions
---------------------------
The explicit demotion input adds ``localization_discordant`` and prevents those
rows from entering the functional view.  Its row count is checked against a
*declared, overridable* expectation
(:data:`EXPECTED_CYS_EXPORT_DEMOTIONS`, default 52, documented by
:data:`EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS`) and every input identifier must be
accounted for exactly once; a mismatch raises instead of proceeding silently.

Boundaries
----------
Synthetic inputs only: this module reads no project data and recomputes no
catalog (plan Authorization Checkpoint 5 is not authorized).  It writes nothing
but the three artifacts below, inside ``--output-dir`` only.  It performs no
family call, no functional calibration, no deletion and no promotion.
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

    ``phaded_evidence_schema`` lives next to this file but the scripts directory
    is not a package, so a plain import only works when it is already on
    ``sys.path``.  Every other import path is a hard error.
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

#: Absent values are written as ``pending`` and are never invented (AGENTS.md).
PENDING = "pending"

# --- controlled vocabularies reused from the shared ontology ---------------
# These are aliases of the ontology sets, not local copies: redefining them here
# would let the catalog and the validators drift apart.
PRIMARY_DISPOSITIONS = SCHEMA.PRIMARY_DISPOSITIONS
MODEL_LAYERS = SCHEMA.MODEL_LAYERS
NUCLEOPHILE_IDENTITIES = SCHEMA.NUCLEOPHILE_IDENTITIES
MOTIF_CLASSES = SCHEMA.MOTIF_CLASSES
TRANSPORT_SIGNAL_PREDICTIONS = SCHEMA.TRANSPORT_SIGNAL_PREDICTIONS
LOCALIZATION_EVIDENCE = SCHEMA.LOCALIZATION_EVIDENCE
CATALYTIC_DOMAIN_TYPES = SCHEMA.CATALYTIC_DOMAIN_TYPES
EVIDENCE_GRADES = SCHEMA.EVIDENCE_GRADES

# --- historical superfamilies with a dedicated primary disposition ---------
#: The lipase-box intracellular family needs a structural competition panel
#: before its function can be stated at all (spec: with-lipase hits stay
#: deferred and are never deleted).
SUPERFAMILY_NPHASCL_WITH_LIPASE_BOX = "intracellular nPHASCL with lipase box"
#: The MCL family is homology-plausible but cannot exclude a competing function,
#: so it is never presented as a functional call.
SUPERFAMILY_NPHAMCL = "intracellular nPHAMCL"
#: Historical Cys export family: transport prediction contradicts the
#: intracellular label.
SUPERFAMILY_NPHASCL_WITHOUT_LIPASE_BOX = "intracellular nPHASCL without lipase box"

#: The only layers that can carry a ``core_sequence_homolog`` verdict.
CORE_ELIGIBLE_LAYERS = frozenset({"sequence_family_hmm_validated", "calibrated_candidate_model"})
#: Layers whose hits are recall-only or reference-only: probable at best.
PROBABLE_LAYERS = frozenset({
    "discovery_hmm_uncalibrated",
    "sequence_family_hmm_validated",
    "calibrated_candidate_model",
})

#: Layer precedence when the input does not state ``model_layer``.  A row only
#: reaches a validated layer on an explicit statement; the default for an
#: unstated layer is the recall-only discovery layer.
_MODEL_LAYER_PRECEDENCE = (
    "calibrated_candidate_model",
    "sequence_family_hmm_validated",
    "reference_query_only",
    "discovery_hmm_uncalibrated",
)
#: Fields whose content implies a validated sequence-family layer.
_VALIDATED_LAYER_FIELDS = ("validated_family_call", "family_model_validated", "sequence_model_validated")

#: Sequence-family confidence reported alongside ``sequence_family_call``.  This
#: is separate from ``experimental_evidence_grade`` (functional evidence) and
#: from ``model_layer`` (which model produced the call).
SEQUENCE_FAMILY_CONFIDENCE_VALUES = ("high", "medium", "low", PENDING)

# --- accessory-domain architecture -----------------------------------------
#: Accessory domains are annotation-only extras.  They are deliberately *not*
#: an input to ``catalytic_domain_type``: SBD in particular never defines
#: type 1 versus type 2.
ACCESSORY_DOMAIN_ORDER = ("SBD", "lid", "linker", "with_lipase")
ACCESSORY_DOMAIN_ABSENT = "not_detected"
ACCESSORY_DOMAIN_INPUT_FIELDS = {
    "SBD": ("sbd_state", "sbd_pf06850_binding_state", "sbd"),
    "lid": ("lid_state", "lid"),
    "linker": ("linker_state", "linker"),
    "with_lipase": ("with_lipase_state", "with_lipase"),
}
ACCESSORY_DOMAIN_SEPARATOR = "|"
ACCESSORY_DOMAIN_JOINER = ":"
_SBD_POSITIVE_PFAM_TOKENS = ("pf06850", "sbd")

# --- evidence flags (multi-valued, non-exclusive) --------------------------
EVIDENCE_FLAG_SEPARATOR = "|"
#: Canonical emission order.  Stable and machine-readable; renaming a flag is a
#: breaking change for consumers of ``phaded_candidate_catalog_v2.tsv``.
EVIDENCE_FLAGS = (
    "ser_motif_lipase_box",
    "ser_motif_ahsmg",
    "cys_associated_motif",
    "sbd_accessory_supported",
    "transport_prediction_supports_export",
    "localization_discordant",
    "localization_unresolved",
    "family_assignment_ambiguous",
    "discovery_layer_only",
    "sequence_family_hmm_validated",
    "functional_calibration_pending",
    "functional_gate_passed_not_promoted",
    "function_competing",
    "structure_review_deferred",
    "input_quality_problem",
    "homology_below_validated_layer",
)
FLAG_SER_LIPASE_BOX = "ser_motif_lipase_box"
FLAG_SER_AHSMG = "ser_motif_ahsmg"
FLAG_CYS_ASSOCIATED = "cys_associated_motif"
FLAG_SBD_ACCESSORY = "sbd_accessory_supported"
FLAG_TRANSPORT_EXPORT = "transport_prediction_supports_export"
FLAG_LOCALIZATION_DISCORDANT = "localization_discordant"
FLAG_LOCALIZATION_UNRESOLVED = "localization_unresolved"
FLAG_ASSIGNMENT_AMBIGUOUS = "family_assignment_ambiguous"
FLAG_DISCOVERY_LAYER = "discovery_layer_only"
FLAG_VALIDATED_LAYER = "sequence_family_hmm_validated"
FLAG_CALIBRATION_PENDING = "functional_calibration_pending"
FLAG_GATE_PASSED_NOT_PROMOTED = "functional_gate_passed_not_promoted"
FLAG_FUNCTION_COMPETING = "function_competing"
FLAG_STRUCTURE_DEFERRED = "structure_review_deferred"
FLAG_INPUT_QUALITY_PROBLEM = "input_quality_problem"
FLAG_BELOW_VALIDATED = "homology_below_validated_layer"

#: Flags that describe broken input rather than a scientific reading.
INPUT_QUALITY_FLAGS = (FLAG_INPUT_QUALITY_PROBLEM,)
#: Flags that describe scientific uncertainty.  None of these may ever appear in
#: :data:`INPUT_QUALITY_FLAGS`, and none of them can pull a row into
#: ``excluded_input_quality``.
SCIENTIFIC_UNCERTAINTY_FLAGS = (
    FLAG_ASSIGNMENT_AMBIGUOUS,
    FLAG_FUNCTION_COMPETING,
    FLAG_STRUCTURE_DEFERRED,
    FLAG_LOCALIZATION_UNRESOLVED,
    FLAG_LOCALIZATION_DISCORDANT,
    FLAG_BELOW_VALIDATED,
    FLAG_DISCOVERY_LAYER,
)

#: Transport predictions that assert an export route.
EXPORT_PREDICTIONS = frozenset({"SP", "LIPO", "TAT", "TATLIPO"})
#: ``no_export`` means the family's historical localization label is
#: intracellular; ``any`` means the family is genuinely periplasmic/ambiguous
#: and a prediction is not a contradiction.
_SUPERFAMILY_EXPORT_REQUIREMENTS = {
    SUPERFAMILY_NPHASCL_WITHOUT_LIPASE_BOX: "no_export",
    "intracellular nPHASCL with lipase box": "no_export",
    "intracellular nPHAMCL": "no_export",
    "extracellular dPHASCL type 1": "export",
    "extracellular dPHASCL type 2": "export",
    "extracellular dPHAMCL": "export",
    "extracellular native-SCL/PhaZ7-like": "export",
    "periplasmic PHA depolymerases": "any",
}

# --- stable, machine-readable rejection rules ------------------------------
RULE_SCIENTIFIC_UNCERTAINTY_NOT_INPUT_ERROR = "scientific_uncertainty_is_not_an_input_error"
RULE_DEMOTION_COUNT_MISMATCH = "demotion_input_count_mismatch"
RULE_DEMOTION_ACCOUNTING = "demotion_accounting_not_exactly_once"
RULE_DUPLICATE_CANDIDATE_ACCESSION = "duplicate_candidate_accession"
RULE_UNKNOWN_EVIDENCE_FLAG = "unknown_evidence_flag"
RULE_SCHEMA_VALIDATION_FAILED = "shared_ontology_validation_failed"
RULE_CLASSIFICATION_FAILED = "classification_failed"
RULE_DISPOSITION_NOT_IN_VOCABULARY = "primary_disposition_not_in_vocabulary"

#: Declared expectation for the 2026-09-20 Cys-export demotion set.  The real
#: value is 52 and is *not* a magic number: it is declared here, documented by
#: :data:`EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS`, overridable via
#: ``--expect-demoted``, and every run records which basis it used.
EXPECTED_CYS_EXPORT_DEMOTIONS = 52
EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS = (
    "frozen v1 evidence, read-only: "
    "runs/20260920_phaded_three_gaps_01/results/cys_demotion/cys_demotion_summary.json "
    "declares moved_to_hold=52 and the authoritative accession list is "
    "runs/20260920_phaded_three_gaps_01/results/cys_demotion/hold_candidates_addendum_52.tsv "
    "(52 data rows); basis: SignalP 6.0h predicts an export signal for 52 "
    "'intracellular' Cys candidates, a PREDICTED conflict, reversible and not an "
    "experimental refutation; the demotion was operator-authorized 2026-09-20"
)
BASIS_DECLARED_DEFAULT = "declared_default"
BASIS_EXPLICIT_EXPECTATION = "explicit_expectation"

#: Candidate-input columns.  Only the first is required; every other column is
#: optional and absent values resolve to their documented weakest default.
CANDIDATE_COLUMNS = (
    "accession",
    "genome",
    "superfamily",
    "sequence_family_call",
    "sequence_family_confidence",
    "assignment_unique",
    "model_layer",
    "functional_calibration_status",
    "catalytic_domain_type",
    "sbd_state",
    "lid_state",
    "linker_state",
    "with_lipase_state",
    "motif_state",
    "lipase_box_state",
    "ahsmg_state",
    "signalp_class",
    "localization_evidence",
    "localization_evidence_source",
    "experimental_evidence_grade",
    "input_quality_reason",
    "evidence_flags",
)

#: v2 catalog columns, in order.  One row per protein, exactly one
#: ``primary_disposition``; ``evidence_flags`` is multi-valued and
#: non-exclusive, and ``functional_view_eligible`` is a derived view column
#: that is never an input to the master disposition.
CATALOG_COLUMNS = (
    "accession",
    "genome",
    "historical_superfamily",
    "sequence_family_call",
    "sequence_family_confidence",
    "catalytic_domain_type",
    "accessory_domain_architecture",
    "transport_signal_prediction",
    "localization_evidence",
    "localization_evidence_source",
    "nucleophile_identity",
    "motif_class",
    "experimental_evidence_grade",
    "model_layer",
    "functional_calibration_status",
    "primary_disposition",
    "evidence_flags",
    "input_quality_reason",
    "scientific_uncertainty_present",
    "functional_view_eligible",
)

#: Input columns carried into ``evidence_flags`` when the input already states
#: them (an input flag is validated against the vocabulary, never passed
#: through blindly).
_INPUT_FLAG_ALIASES = {
    "ahsmg": FLAG_SER_AHSMG,
    "ahsmg_ser_motif": FLAG_SER_AHSMG,
    "ser_motif": FLAG_SER_LIPASE_BOX,
    "function_competing": FLAG_FUNCTION_COMPETING,
    "competing_function": FLAG_FUNCTION_COMPETING,
    "localization_discordant": FLAG_LOCALIZATION_DISCORDANT,
    "cys_export_demotion": FLAG_LOCALIZATION_DISCORDANT,
}
_FLAG_SPLIT = re.compile(r"[|,;]+")
_TRUE_TOKENS = frozenset(SCHEMA._TRUE_TOKENS)
_FALSE_TOKENS = frozenset(SCHEMA._FALSE_TOKENS)


class EvidenceError(ValueError):
    """A catalog row or input violates a documented evidence rule."""

    def __init__(self, rule: str, message: str) -> None:
        super().__init__(f"{rule}: {message}")
        self.rule = rule
        self.message = message


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _text(value: object) -> str:
    """Read an evidence field defensively; absent and ``None`` read as ``""``."""
    if value is None:
        return ""
    return str(value).strip()


def _text_value(row: Mapping[str, object], field: str) -> str:
    if field not in row:
        return ""
    return _text(row[field])


def _first_text(row: Mapping[str, object], fields: Sequence[str]) -> str:
    for field in fields:
        value = _text_value(row, field)
        if value:
            return value
    return ""


def _flag_state(value: object) -> "bool | None":
    """Normalize a supported/unsupported style state.

    ``None`` means "not stated" (which is *not* the same as "not detected").
    A "not detected" spelling such as ``not_detected_pattern`` is a *stated
    absence*, i.e. ``False``, not an unknown state.
    """
    raw = _text(value)
    if not raw:
        return None
    key = re.sub(r"[^a-z0-9]+", "_", raw.lower()).strip("_")
    if key in _TRUE_TOKENS:
        return True
    if key in _FALSE_TOKENS or key.startswith("not_detected") or key.startswith("no_"):
        return False
    raise EvidenceError(
        RULE_CLASSIFICATION_FAILED,
        f"{raw!r} is not a recognized supported/unsupported state",
    )


def _canonical(value: object, field: str, allowed: Iterable[str], *, default: str = "") -> str:
    """Resolve ``value`` against a closed vocabulary from the shared ontology."""
    raw = _text(value)
    if not raw:
        return default
    key = re.sub(r"[^a-z0-9]+", "_", raw.lower()).strip("_")
    for candidate in allowed:
        if re.sub(r"[^a-z0-9]+", "_", candidate.lower()).strip("_") == key:
            return candidate
    raise EvidenceError(
        RULE_CLASSIFICATION_FAILED,
        f"{field}={raw!r} is not a recognized value; expected one of {sorted(allowed)}",
    )


def _is_truthy(value: object) -> bool:
    key = re.sub(r"[^a-z0-9]+", "_", _text(value).lower()).strip("_")
    return key in _TRUE_TOKENS


def _is_falsey(value: object) -> bool:
    key = re.sub(r"[^a-z0-9]+", "_", _text(value).lower()).strip("_")
    return key in _FALSE_TOKENS


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: "str | Path") -> "tuple[list[str], list[dict[str, str]]]":
    tsv_path = Path(path)
    if not tsv_path.is_file() or tsv_path.is_symlink():
        raise FileNotFoundError(f"input is not a regular file: {tsv_path}")
    with tsv_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or ())
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t",
            lineterminator="\n", extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _bound_file(path: "str | Path") -> dict[str, Any]:
    """Describe an input file for the manifest (path, size, SHA-256)."""
    file_path = Path(path)
    return {
        "path": str(file_path),
        "size": file_path.stat().st_size,
        "sha256": _sha256(file_path),
    }


# ---------------------------------------------------------------------------
# Evidence-flag vocabulary
# ---------------------------------------------------------------------------

def normalize_evidence_flags(raw: object) -> str:
    """Normalize a multi-valued flag list into canonical order.

    Accepts ``|``, ``,`` or ``;`` separated input, maps documented aliases,
    de-duplicates, and rejects anything outside :data:`EVIDENCE_FLAGS` so the
    emitted vocabulary stays closed and auditable.
    """
    text = _text(raw)
    if not text:
        return ""
    seen: set[str] = set()
    for token in _FLAG_SPLIT.split(text):
        token = token.strip()
        if not token:
            continue
        key = re.sub(r"[^a-z0-9]+", "_", token.lower()).strip("_")
        resolved = _INPUT_FLAG_ALIASES.get(key, token if token in EVIDENCE_FLAGS else "")
        if not resolved:
            for candidate in EVIDENCE_FLAGS:
                if re.sub(r"[^a-z0-9]+", "_", candidate.lower()).strip("_") == key:
                    resolved = candidate
                    break
        if not resolved:
            raise EvidenceError(
                RULE_UNKNOWN_EVIDENCE_FLAG,
                f"evidence flag {token!r} is not in the declared vocabulary {list(EVIDENCE_FLAGS)}",
            )
        seen.add(resolved)
    return EVIDENCE_FLAG_SEPARATOR.join(flag for flag in EVIDENCE_FLAGS if flag in seen)


def _ordered_flags(flags: Iterable[str]) -> str:
    wanted = set(flags)
    return EVIDENCE_FLAG_SEPARATOR.join(flag for flag in EVIDENCE_FLAGS if flag in wanted)


# ---------------------------------------------------------------------------
# Input quality: the guard around ``excluded_input_quality``
# ---------------------------------------------------------------------------

def input_quality_reason(row: Mapping[str, object]) -> str:
    """Return the documented input-quality reason, or ``""`` when there is none.

    Delegates to :func:`phaded_evidence_schema._resolve_input_quality_reason`, so
    a row that writes scientific-uncertainty wording into an input-quality field
    is rejected here as well as by the shared validator.  The rejection is
    re-raised as a named rule so the failing direction is explicit in tests and
    logs.
    """
    try:
        return SCHEMA._resolve_input_quality_reason(row)
    except ValueError as error:
        raise EvidenceError(
            RULE_SCIENTIFIC_UNCERTAINTY_NOT_INPUT_ERROR,
            f"{error}; primary_disposition='excluded_input_quality' is limited to missing, "
            "truncated or malformed input and scientific uncertainty is not an input error",
        ) from None


def _has_documented_input_quality(row: Mapping[str, object]) -> bool:
    return bool(input_quality_reason(row))


# ---------------------------------------------------------------------------
# Field-level semantics
# ---------------------------------------------------------------------------

def resolve_model_layer(row: Mapping[str, object]) -> str:
    """Resolve the sequence/functional model layer of a row.

    An explicit valid ``model_layer`` always wins.  Otherwise an explicit
    functional calibration status names its layer, then any validated-family
    evidence, then the field precedence in :data:`_MODEL_LAYER_PRECEDENCE`.  A
    row that states nothing is treated as recall-only
    (``discovery_hmm_uncalibrated``), never as validated.
    """
    stated = _text_value(row, "model_layer")
    if stated:
        return _canonical(stated, "model_layer", MODEL_LAYERS)
    status = _text_value(row, "functional_calibration_status")
    if status == "calibrated_candidate_model":
        return "calibrated_candidate_model"
    for field in _VALIDATED_LAYER_FIELDS:
        if _is_truthy(_text_value(row, field)):
            return "sequence_family_hmm_validated"
    for field in ("trained_profile_hit", "profile_evidence_status"):
        if _is_truthy(_text_value(row, field)):
            return "sequence_family_hmm_validated"
    for layer in _MODEL_LAYER_PRECEDENCE:
        if _is_truthy(_text_value(row, layer)):
            return layer
    return "discovery_hmm_uncalibrated"


def resolve_functional_calibration_status(row: Mapping[str, object]) -> str:
    """Resolve the functional calibration status independently of the layer.

    Absent status is ``not_function_calibrated``: sequence evidence never
    implies a functional calibration, and this module performs none.
    """
    stated = _text_value(row, "functional_calibration_status")
    if stated:
        return _canonical(
            stated, "functional_calibration_status",
            SCHEMA.FUNCTIONAL_CALIBRATION_STATUSES,
        )
    return "not_function_calibrated"


def motif_state_of(row: Mapping[str, object]) -> "tuple[bool, str]":
    """Resolve ``(ahsmg_supported, explicit_motif_state)`` for the ontology.

    The three historical inputs are kept separate so the shared
    :func:`phaded_evidence_schema.normalize_nucleophile` -- not this module --
    decides the identity/motif pair.
    """
    ahsmg_supported = _flag_state(_text_value(row, "ahsmg_state")) is True
    stated = _text_value(row, "motif_state")
    return ahsmg_supported, stated


def _is_superfamily_key(superfamily: str, needle: str) -> bool:
    key = re.sub(r"[^a-z0-9]+", "", _text(superfamily).lower())
    return needle in key


def _superfamily_is_scl(superfamily: str) -> bool:
    """Whether the historical superfamily is an SCL (Ser hydrolase) family."""
    return "scl" in re.sub(r"[^a-z0-9]+", "", _text(superfamily).lower())


def resolve_nucleophile(row: Mapping[str, object]) -> dict:
    """Return ``nucleophile_identity`` and ``motif_class`` for one row.

    Delegates entirely to :func:`phaded_evidence_schema.normalize_nucleophile`;
    AHSMG is therefore always a Ser motif class and never a third nucleophile
    identity.  Absent GxSxG in an SCL family is the Cys-associated case; when no
    motif evidence is informative and no superfamily prior applies, both stay
    ``unresolved`` instead of guessing.
    """
    ahsmg_supported, stated = motif_state_of(row)
    lipase_box = _flag_state(_text_value(row, "lipase_box_state"))
    superfamily = _text_value(row, "superfamily")

    motif_state = stated
    if not motif_state:
        if lipase_box is True:
            motif_state = "GxSxG"
        elif ahsmg_supported:
            motif_state = "AHSMG"
        elif lipase_box is False and _superfamily_is_scl(superfamily):
            motif_state = "Cys-associated"
        else:
            motif_state = ""
    out = SCHEMA.normalize_nucleophile(
        superfamily, motif_state, "supported" if ahsmg_supported else "",
    )
    if out["motif_class"] == "unresolved" and _is_superfamily_key(superfamily, "nphamcl"):
        # An MCL family is not assumed to carry a canonical Ser lipase box.
        out["motif_class"] = "other"
    return out


def resolve_catalytic_domain_type(row: Mapping[str, object]) -> str:
    """Resolve the type 1 / type 2 call from domain geometry *only*.

    SBD, lid, linker and with-lipase annotations are deliberately not consulted:
    an accessory domain never defines type 1 versus type 2.  An unrecognized
    domain-type token (for example a v1 geometry verdict such as
    ``undetermined_no_oxyanion``) resolves to ``unresolved`` rather than being
    silently relabelled as a verified type.
    """
    stated = _text_value(row, "catalytic_domain_type")
    if not stated:
        return "unresolved"
    try:
        return _canonical(stated, "catalytic_domain_type", CATALYTIC_DOMAIN_TYPES)
    except EvidenceError:
        return "unresolved"


def _accessory_state_is_absent(state: str) -> bool:
    """Whether an accessory-domain state is a plain "not there" reading.

    A stated absence carries no architecture information, so it is omitted from
    ``accessory_domain_architecture``; a supporting or otherwise informative
    state (including a lookalike state such as an undetermined geometry) is kept
    verbatim so nothing is silently rewritten.
    """
    key = re.sub(r"[^a-z0-9]+", "_", state.lower()).strip("_")
    if key in _FALSE_TOKENS or key in {"not_detected_pattern", "undetermined"}:
        return True
    return key.startswith("not_detected") or key.startswith("no_")


def accessory_domain_architecture(row: Mapping[str, object]) -> str:
    """Render the accessory-domain architecture independently of the type call."""
    parts: list[str] = []
    for domain in ACCESSORY_DOMAIN_ORDER:
        state = _first_text(row, ACCESSORY_DOMAIN_INPUT_FIELDS[domain])
        if not state or _accessory_state_is_absent(state):
            continue
        parts.append(f"{domain}{ACCESSORY_DOMAIN_JOINER}{state}")
    return ACCESSORY_DOMAIN_SEPARATOR.join(parts)


def _sbd_supported(row: Mapping[str, object]) -> bool:
    for field in ACCESSORY_DOMAIN_INPUT_FIELDS["SBD"]:
        raw = _text_value(row, field)
        if not raw:
            continue
        if _flag_state(raw) is True:
            return True
        key = re.sub(r"[^a-z0-9]+", "_", raw.lower()).strip("_")
        if any(token in key for token in _SBD_POSITIVE_PFAM_TOKENS):
            return True
    return False


def resolve_transport_signal(row: Mapping[str, object]) -> str:
    """Return the predicted transport state only.

    The demotion input separates the frozen ``signalp_class`` from the
    authoritative ``signalp_class_prediction``; the latter wins when both are
    present.  An unrecognized class resolves to ``invalid_input`` (a non-signal),
    never to an export-predicting class.
    """
    stated = _first_text(row, ("signalp_class_prediction", "signalp_class"))
    if not stated:
        return "not_tested"
    try:
        return _canonical(stated, "transport_signal_prediction", TRANSPORT_SIGNAL_PREDICTIONS)
    except EvidenceError:
        return "invalid_input"


def resolve_localization_evidence(row: Mapping[str, object]) -> "tuple[str, str]":
    """Return ``(localization_evidence, localization_evidence_source)``.

    The SignalP class is never consulted: a transport prediction is not
    localization evidence.  Absent evidence is ``unknown``.
    """
    stated = _text_value(row, "localization_evidence")
    try:
        evidence = _canonical(
            stated, "localization_evidence", LOCALIZATION_EVIDENCE, default="unknown"
        )
    except EvidenceError:
        evidence = "unknown"
    return evidence, _text_value(row, "localization_evidence_source")


def resolve_sequence_family_confidence(row: Mapping[str, object], layer: str) -> str:
    """Resolve the sequence-family confidence, separately from experimental grade.

    When nothing is stated *and* the layer itself had to be defaulted, the
    confidence is emitted as the empty string: a defaulted layer must not
    manufacture a confidence level.  ``pending`` is reserved for a row that
    states a validated layer but no assignment uniqueness.
    """
    stated = _text_value(row, "sequence_family_confidence")
    if stated:
        return _canonical(
            stated, "sequence_family_confidence", SEQUENCE_FAMILY_CONFIDENCE_VALUES
        )
    if not _text_value(row, "model_layer"):
        # The layer itself had to be defaulted, so no confidence was stated by
        # any input: emit the absence, not a confidence level.
        return ""
    unique = _text_value(row, "assignment_unique")
    if layer == "sequence_family_hmm_validated":
        if _is_truthy(unique):
            return "high"
        if _is_falsey(unique):
            return "medium"
        return PENDING
    if layer == "calibrated_candidate_model":
        return "high"
    if layer in {"discovery_hmm_uncalibrated", "reference_query_only"}:
        return "low"
    return PENDING


def localization_discordant(row: Mapping[str, object], prediction: str) -> bool:
    """Whether the prediction contradicts the direction the family implies.

    A prediction is only a contradiction when the historical family implies a
    direction: ``export`` families need an export-predicting class, ``no_export``
    (intracellular) families must not get one.  Peri-ambiguous families and
    ``not_tested``/``invalid_input`` predictions are never contradictions.
    """
    requirement = _SUPERFAMILY_EXPORT_REQUIREMENTS.get(_text_value(row, "superfamily"))
    if requirement is None:
        requirement = "any"
    if prediction in {"not_tested", "invalid_input"}:
        return False
    exports = prediction in EXPORT_PREDICTIONS
    if requirement == "export":
        return not exports
    if requirement == "no_export":
        return exports
    return False


def resolve_family_call(row: Mapping[str, object], layer: str) -> str:
    """Return the sequence-family call for a row, honouring the discovery boundary.

    Two documented conventions meet here:

    * an absent value is the empty string, never an invented token -- so
      "no family call was made" is not confused with a family literally named
      ``pending`` (the ontology's own weakest default for this column is ``""``);
    * the discovery layer only scores existing intermediate products and never
      makes a family call, so a discovery-layer row is emitted empty even when
      the input carries a candidate value.
    """
    if layer == "discovery_hmm_uncalibrated":
        return ""
    return _text_value(row, "sequence_family_call")


def certificate_grade(row: Mapping[str, object]) -> str:
    """Return the experimental evidence grade as the shared ontology defines it.

    ``EVIDENCE_GRADES`` has no ``pending`` member and the ontology's documented
    weakest default for this column is the empty string, so an unavailable grade
    is emitted as ``""`` -- never as an invented grade, and never as a token the
    shared validator would reject.  An explicit ``pending``/``unknown`` is
    treated as "not yet graded", not as a value.
    """
    stated = _text_value(row, "experimental_evidence_grade")
    if not stated:
        return ""
    if re.sub(r"[^a-z0-9]+", "_", stated.lower()).strip("_") in {"pending", "unknown", "none", "na"}:
        return ""
    return _canonical(stated, "experimental_evidence_grade", EVIDENCE_GRADES)


# ---------------------------------------------------------------------------
# Primary disposition
# ---------------------------------------------------------------------------

def primary_disposition(row: Mapping[str, object]) -> str:
    """Return the single consolidated disposition of one candidate row.

    Precedence, in order:

    1. a documented input-quality reason -> ``excluded_input_quality`` (broken
       input cannot be classified at all);
    2. ``intracellular nPHASCL with lipase box`` -> ``deferred_structure_review``;
    3. ``intracellular nPHAMCL`` -> ``function_unresolved``;
    4. a validated (or calibrated) layer with a unique assignment ->
       ``core_sequence_homolog``;
    5. a validated or discovery layer -> ``probable_sequence_homolog``;
    6. anything else -> ``remote_homolog_candidate``.

    Steps 2-3 precede steps 4-5 on purpose: "special family before model layer"
    is a scientific-integrity rule, so a validated model layer can never promote
    a competing-function family (nPHAMCL) or a structure-deferred family.

    ``excluded_input_quality`` is reachable *only* through step 1, which requires
    a documented input-quality reason.  Ambiguous family, unresolved
    localization, competing function and low confidence therefore never route
    there.
    """
    superfamily = _text_value(row, "superfamily")

    # Rule 1: input quality.  Only a documented input-quality reason qualifies;
    # `input_quality_reason()` raises rather than returning a value when the
    # "reason" is scientific-uncertainty wording.
    if input_quality_reason(row):
        return "excluded_input_quality"

    # Rule 2: with-lipase hits need a structural competition panel.
    if superfamily == SUPERFAMILY_NPHASCL_WITH_LIPASE_BOX:
        return "deferred_structure_review"

    # Rule 3: the MCL family cannot exclude a competing function.
    if superfamily == SUPERFAMILY_NPHAMCL:
        return "function_unresolved"

    # Rule 4: a core homolog requires an independently validated sequence family
    # model *and* a unique assignment.
    layer = resolve_model_layer(row)
    if layer in CORE_ELIGIBLE_LAYERS and _is_truthy(_text_value(row, "assignment_unique")):
        return "core_sequence_homolog"

    # Rule 5: any validated or discovery-layer hit is at least probable.
    if layer in PROBABLE_LAYERS:
        return "probable_sequence_homolog"

    # Rule 6: reference-only or unrecognized evidence is remote.
    return "remote_homolog_candidate"


def _disposition_is_valid(disposition: str) -> bool:
    return disposition in PRIMARY_DISPOSITIONS


# ---------------------------------------------------------------------------
# Row classification
# ---------------------------------------------------------------------------

def classify_row(row: Mapping[str, object]) -> dict[str, str]:
    """Classify one candidate row into the v2 evidence columns.

    Every emitted string is one of: a documented default (``pending`` /
    ``unknown`` / ``unresolved``), a value copied from the input, or a value
    derived by the shared ontology.  Nothing is invented and no value is
    silently promoted -- in particular a SignalP class never becomes
    ``localization_evidence`` and an accessory domain never becomes a
    ``catalytic_domain_type``.
    """
    raw = dict(row)
    superfamily = _text_value(raw, "superfamily")
    layer = resolve_model_layer(raw)
    status = resolve_functional_calibration_status(raw)
    try:
        SCHEMA.validate_model_state(layer, status)
    except ValueError as error:
        raise EvidenceError(RULE_SCHEMA_VALIDATION_FAILED, str(error)) from None

    nucleophile = resolve_nucleophile(raw)
    prediction = resolve_transport_signal(raw)
    localization, localization_source = resolve_localization_evidence(raw)
    domain_type = resolve_catalytic_domain_type(raw)
    architecture = accessory_domain_architecture(raw)
    discordant = localization_discordant(raw, prediction)

    # Input quality is resolved unconditionally so that a row which writes
    # scientific-uncertainty wording into an input-quality field is rejected
    # even when its disposition would not have been an exclusion.
    reason = input_quality_reason(raw)
    disposition = primary_disposition(raw)
    if not _disposition_is_valid(disposition):  # pragma: no cover - defensive
        raise EvidenceError(
            RULE_DISPOSITION_NOT_IN_VOCABULARY,
            f"{disposition!r} is not a declared primary disposition",
        )
    unique = _is_truthy(_text_value(raw, "assignment_unique"))

    flags: list[str] = []
    motif_class = nucleophile["motif_class"]
    if motif_class == "AHSMG":
        flags.append(FLAG_SER_AHSMG)
    elif motif_class == "GxSxG":
        flags.append(FLAG_SER_LIPASE_BOX)
    elif motif_class == "Cys-associated":
        flags.append(FLAG_CYS_ASSOCIATED)
    if _sbd_supported(raw):
        flags.append(FLAG_SBD_ACCESSORY)
    if prediction in EXPORT_PREDICTIONS:
        flags.append(FLAG_TRANSPORT_EXPORT)
    if discordant:
        flags.append(FLAG_LOCALIZATION_DISCORDANT)
    if localization == "unknown" and not localization_source:
        flags.append(FLAG_LOCALIZATION_UNRESOLVED)

    if layer == "sequence_family_hmm_validated":
        flags.append(FLAG_VALIDATED_LAYER)
        if not unique:
            flags.append(FLAG_ASSIGNMENT_AMBIGUOUS)
    elif layer == "discovery_hmm_uncalibrated":
        flags.append(FLAG_DISCOVERY_LAYER)
    elif layer == "reference_query_only":
        flags.append(FLAG_BELOW_VALIDATED)

    if status == "candidate_gate_passed_not_promoted":
        flags.append(FLAG_GATE_PASSED_NOT_PROMOTED)
    else:
        flags.append(FLAG_CALIBRATION_PENDING)
    # Scientific-uncertainty flags are properties of the evidence, not rewards
    # for a particular disposition: a competing-function family stays flagged
    # even when the row is separately excluded for broken input.
    if superfamily == SUPERFAMILY_NPHAMCL:
        flags.append(FLAG_FUNCTION_COMPETING)
    if superfamily == SUPERFAMILY_NPHASCL_WITH_LIPASE_BOX:
        flags.append(FLAG_STRUCTURE_DEFERRED)
    if disposition == "excluded_input_quality":
        flags.append(FLAG_INPUT_QUALITY_PROBLEM)

    flags.extend(
        _FLAG_SPLIT.split(normalize_evidence_flags(_text_value(raw, "evidence_flags")))
    )

    scientific = bool(set(flags) & set(SCIENTIFIC_UNCERTAINTY_FLAGS)) or (
        disposition in {"function_unresolved", "deferred_structure_review"}
    )

    return {
        "accession": _text_value(raw, "accession"),
        "genome": _text_value(raw, "genome"),
        "historical_superfamily": superfamily if superfamily else PENDING,
        "sequence_family_call": resolve_family_call(raw, layer),
        "sequence_family_confidence": resolve_sequence_family_confidence(raw, layer),
        "catalytic_domain_type": domain_type,
        "accessory_domain_architecture": architecture,
        "transport_signal_prediction": prediction,
        "localization_evidence": localization,
        "localization_evidence_source": localization_source,
        "nucleophile_identity": nucleophile["nucleophile_identity"],
        "motif_class": motif_class,
        "experimental_evidence_grade": certificate_grade(raw),
        "model_layer": layer,
        "functional_calibration_status": status,
        "primary_disposition": disposition,
        "evidence_flags": _ordered_flags(flags),
        "input_quality_reason": reason,
        "scientific_uncertainty_present": "true" if scientific else "false",
    }


def classify_rows(rows: Iterable[Mapping[str, object]]) -> "list[dict[str, str]]":
    """Classify many rows, failing closed on the first bad row."""
    return [classify_row(row) for row in rows]


# ---------------------------------------------------------------------------
# Derived functional view (never a master disposition)
# ---------------------------------------------------------------------------

def functional_view_eligible(classified_row: Mapping[str, object]) -> bool:
    """Whether a *classified* row may appear in the functional-support view.

    This is a derived annotation.  It reads ``primary_disposition`` but never
    writes it, and :func:`primary_disposition` never reads this function's
    output: the master disposition is the master catalog's, and the functional
    view is a strict subset of it.

    Requirements: a core sequence homolog, a validated sequence-family model
    with a unique assignment, a real experimental grade (never ``A`` or
    ``pending``), localization that is not discordant, unresolved or merely a
    historical label, and no disqualifying flag (ambiguity, competing function,
    deferred structure, broken input).

    ``functional_calibration_status`` is deliberately *not* a requirement: this
    column answers "may this sequence support a functional statement?", whereas
    the calibration status answers "has a functional model been accepted?".
    Treating a sequence-qualified row as ineligible would re-collapse the two
    claims that v2 keeps apart.
    """
    if _text(classified_row.get("primary_disposition")) != "core_sequence_homolog":
        return False
    if _text(classified_row.get("model_layer")) != "sequence_family_hmm_validated":
        return False
    grade = _text(classified_row.get("experimental_evidence_grade"))
    if grade not in {"E3", "E2"}:
        return False
    if _text(classified_row.get("localization_evidence")) not in {
        "experimental_confirmed", "other_evidence",
    }:
        return False
    flags = set(
        token
        for token in _FLAG_SPLIT.split(_text(classified_row.get("evidence_flags")))
        if token
    )
    disqualifying = {
        FLAG_LOCALIZATION_DISCORDANT,
        FLAG_LOCALIZATION_UNRESOLVED,
        FLAG_ASSIGNMENT_AMBIGUOUS,
        FLAG_FUNCTION_COMPETING,
        FLAG_STRUCTURE_DEFERRED,
        FLAG_INPUT_QUALITY_PROBLEM,
        FLAG_BELOW_VALIDATED,
        FLAG_DISCOVERY_LAYER,
    }
    return not (flags & disqualifying)


def build_functional_view(classified_rows: Iterable[Mapping[str, object]]) -> "list[dict[str, str]]":
    """Return the derived functional-support view as a new list.

    The input rows are not mutated: the view is a copy of the eligible subset,
    so computing it cannot change any master catalog disposition.
    """
    return [dict(row) for row in classified_rows if functional_view_eligible(row)]


# ---------------------------------------------------------------------------
# Candidate universe
# ---------------------------------------------------------------------------

def unify_universe(
    candidate_ids: Iterable[object], demotion_ids: Iterable[object]
) -> dict[str, Any]:
    """Union the candidate universe with the explicit demotion identifiers.

    Returns the sorted identifier list (which defines "exactly once") plus the
    per-source counts and the identifiers that exist *only* in the demotion
    input.  Duplicate identifiers inside one source are a hard error: an
    ambiguous key cannot be accounted for exactly once.
    """
    candidates = [(_text(value), "candidates") for value in candidate_ids]
    demotions = [(_text(value), "demotion") for value in demotion_ids]
    for accession, source in candidates + demotions:
        if not accession:
            raise EvidenceError(
                RULE_DEMOTION_ACCOUNTING, f"empty accession in the {source} input"
            )

    per_source: dict[str, list[str]] = {"candidates": [], "demotion": []}
    duplicates: list[str] = []
    for source, values in (("candidates", candidates), ("demotion", demotions)):
        counts: dict[str, int] = {}
        for accession, _ in values:
            counts[accession] = counts.get(accession, 0) + 1
        duplicates.extend(
            f"{source}:{accession}" for accession, count in counts.items() if count > 1
        )
        per_source[source] = sorted(counts)

    return {
        "accessions": sorted(set(per_source["candidates"]) | set(per_source["demotion"])),
        "candidate_count": len(per_source["candidates"]),
        "demotion_count": len(per_source["demotion"]),
        "demotion_only": sorted(set(per_source["demotion"]) - set(per_source["candidates"])),
        "in_both": sorted(set(per_source["demotion"]) & set(per_source["candidates"])),
        "duplicates": sorted(duplicates),
    }


def account_demotions(
    demotion_rows: Sequence[Mapping[str, object]],
    *,
    expected: int,
    classified_accessions: Iterable[object] = (),
    candidate_accessions: Iterable[object] = (),
) -> dict[str, Any]:
    """Account for every explicit demotion identifier exactly once.

    Raises :class:`EvidenceError` when the row count contradicts the declared
    expectation (``RULE_DEMOTION_COUNT_MISMATCH``) or when an identifier is
    repeated, missing from the catalog, or otherwise not exactly once
    (``RULE_DEMOTION_ACCOUNTING``).  Never proceeds silently on a mismatch.
    """
    ids_in_order = [_text(row.get("accession")) for row in demotion_rows]
    ids = [value for value in ids_in_order if value]
    blank = len(ids_in_order) - len(ids)
    counts: dict[str, int] = {}
    for accession in ids:
        counts[accession] = counts.get(accession, 0) + 1
    duplicates = sorted(accession for accession, count in counts.items() if count > 1)
    distinct = sorted(counts)
    catalog = {_text(value) for value in classified_accessions}
    candidates = {_text(value) for value in candidate_accessions}
    missing = sorted(accession for accession in distinct if accession not in catalog)

    if blank:
        raise EvidenceError(
            RULE_DEMOTION_ACCOUNTING,
            f"the demotion input contains {blank} row(s) without an accession",
        )
    if duplicates:
        raise EvidenceError(
            RULE_DEMOTION_ACCOUNTING,
            f"the demotion input repeats {len(duplicates)} accession(s), so they cannot be "
            f"accounted for exactly once: {duplicates[:5]}",
        )
    if len(ids) != expected:
        raise EvidenceError(
            RULE_DEMOTION_COUNT_MISMATCH,
            f"the demotion input has {len(ids)} row(s) but the declared expectation is "
            f"{expected}; refusing to proceed silently (override --expect-demoted to state "
            "a different, auditable basis)",
        )
    if missing:
        raise EvidenceError(
            RULE_DEMOTION_ACCOUNTING,
            f"{len(missing)} demotion accession(s) are absent from the catalog: {missing[:5]}",
        )

    return {
        "declared_expected": expected,
        "input_rows": len(ids_in_order),
        "input_ids": len(distinct),
        "accounted_ids": len(distinct),
        "missing_ids": [],
        "duplicate_ids": [],
        "duplicates_within_universe": [],
        "demotion_only_ids": len([value for value in distinct if value not in candidates]),
        "matched_candidate_ids": len([value for value in distinct if value in candidates]),
        "exactly_once": len(distinct) == len(ids) == expected,
        "counts_match_declared": len(ids) == expected,
    }


# ---------------------------------------------------------------------------
# Catalog assembly
# ---------------------------------------------------------------------------

def build_catalog_rows(
    candidate_rows: Sequence[Mapping[str, object]],
    demotion_rows: Sequence[Mapping[str, object]] = (),
    *,
    expect_demoted: int = EXPECTED_CYS_EXPORT_DEMOTIONS,
    basis: str = "",
    basis_detail: str = "",
) -> dict[str, Any]:
    """Classify the unified universe and account for the demotion input.

    Returns ``{"rows": [...], "summary": {...}, "manifest": {...},
    "functional_view": [...]}``.  Every identifier in the universe appears
    exactly once in ``rows``, and every row is validated by the shared ontology
    before it is returned.
    """
    if expect_demoted < 0:
        raise EvidenceError(
            RULE_DEMOTION_COUNT_MISMATCH,
            f"--expect-demoted must not be negative; got {expect_demoted}",
        )
    if not basis:
        basis = (
            BASIS_DECLARED_DEFAULT
            if expect_demoted == EXPECTED_CYS_EXPORT_DEMOTIONS
            else BASIS_EXPLICIT_EXPECTATION
        )
    if not basis_detail:
        # A caller-supplied basis statement is its own detail; only the declared
        # default borrows the frozen evidence basis.
        basis_detail = (
            EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS
            if basis == BASIS_DECLARED_DEFAULT
            else basis
        )

    candidates = {_text(row.get("accession")): dict(row) for row in candidate_rows}
    demotions = {_text(row.get("accession")): dict(row) for row in demotion_rows}
    universe = unify_universe(candidates, demotions)
    if universe["duplicates"]:
        raise EvidenceError(
            RULE_DEMOTION_ACCOUNTING,
            "an accession is repeated inside one input, so it cannot be accounted for "
            f"exactly once: {universe['duplicates'][:5]}",
        )

    rows: list[dict[str, str]] = []
    for accession in universe["accessions"]:
        # The explicit demotion input carries the authoritative transport
        # prediction for a demoted row; the candidate row keeps its own values
        # for every other field.
        source: dict[str, object] = dict(candidates.get(accession, {}))
        demotion = demotions.get(accession, {})
        if demotion:
            source["accession"] = accession
            for field in ("signalp_class_prediction", "signalp_class_prediction_note"):
                if _text(demotion.get(field)):
                    source[field] = demotion[field]
            if not _text(source.get("signalp_class")):
                source["signalp_class"] = _text(demotion.get("signalp_class"))
            for field in CANDIDATE_COLUMNS:
                if not _text(source.get(field)) and _text(demotion.get(field)):
                    source[field] = demotion[field]
            source["explicit_demotion"] = "true"
        try:
            classified = classify_row(source)
        except EvidenceError as error:
            raise EvidenceError(error.rule, f"{accession}: {error.message}") from None
        rows.append(classified)

    # The functional view is derived, never a master verdict.  It is computed
    # only after every master disposition is final.
    for row in rows:
        row["functional_view_eligible"] = (
            "true" if functional_view_eligible(row) else "false"
        )
    view = build_functional_view(rows)

    accounting = account_demotions(
        demotion_rows, expected=expect_demoted,
        classified_accessions=[row["accession"] for row in rows],
        candidate_accessions=candidates,
    )
    accounting["basis"] = basis
    accounting["basis_detail"] = basis_detail
    accounting["non_demotion_rows"] = len(rows) - accounting["demotion_only_ids"]

    # Fail closed on any row the shared ontology rejects.
    for row in rows:
        try:
            SCHEMA.validate_candidate_evidence(row)
        except ValueError as error:
            raise EvidenceError(
                RULE_SCHEMA_VALIDATION_FAILED, f"{row['accession']}: {error}"
            ) from None

    counts: dict[str, int] = {}
    for row in rows:
        counts[row["primary_disposition"]] = counts.get(row["primary_disposition"], 0) + 1
    flags: dict[str, int] = {}
    for row in rows:
        for flag in _FLAG_SPLIT.split(row["evidence_flags"]):
            if flag:
                flags[flag] = flags.get(flag, 0) + 1

    accessions = [row["accession"] for row in rows]
    summary = {
        "catalog_rows": len(rows),
        "primary_disposition_counts": dict(sorted(counts.items())),
        "primary_disposition_values_seen": sorted(counts),
        "one_disposition_per_protein": len(accessions) == len(set(accessions)),
        "evidence_flag_counts": dict(sorted(flags.items())),
        "evidence_flag_vocabulary": list(EVIDENCE_FLAGS),
        "functional_view_eligible_rows": len(view),
        "functional_view_cannot_change_disposition": True,
        "localization_discordant_rows": sum(
            1 for row in rows if FLAG_LOCALIZATION_DISCORDANT in row["evidence_flags"]
        ),
        "excluded_input_quality_rows": counts.get("excluded_input_quality", 0),
        "scientific_uncertainty_rows": sum(
            1 for row in rows if row["scientific_uncertainty_present"] == "true"
        ),
        "universe": {
            "catalog_rows": len(rows),
            "candidate_count": universe["candidate_count"],
            "demotion_only_count": len(universe["demotion_only"]),
            "in_both_count": len(universe["in_both"]),
        },
        "demotion_accounting": accounting,
    }
    manifest = {
        "schema_version": "1.0",
        "catalog": "phaded_candidate_catalog_v2",
        "semantics": "evidence-separated (spec 2026-09-28)",
        "columns": list(CATALOG_COLUMNS),
        "primary_dispositions": sorted(PRIMARY_DISPOSITIONS),
        "evidence_flag_vocabulary": list(EVIDENCE_FLAGS),
        "counts": {
            "catalog_rows": len(rows),
            "candidates_input_rows": universe["candidate_count"],
            "demotion_input_rows": accounting["input_rows"],
            "functional_view_eligible_rows": len(view),
        },
        "demotion_basis": {
            "declared_expected": expect_demoted,
            "basis": basis,
            "basis_detail": basis_detail,
            "override_flag": "--expect-demoted",
        },
        "boundaries": {
            "phenotype_claimed": False,
            "functional_calibration_performed": False,
            "family_call_made": False,
            "candidates_deleted": False,
            "v1_frozen_results_untouched": True,
            "one_disposition_per_protein": summary["one_disposition_per_protein"],
            "localization_is_annotation_only": True,
            "sbd_defines_catalytic_domain_type": False,
            "ahsmg_is_a_ser_motif_class": True,
        },
        "model_layers": sorted(MODEL_LAYERS),
        "nucleophile_identities": sorted(NUCLEOPHILE_IDENTITIES),
        "motif_classes": sorted(MOTIF_CLASSES),
    }
    return {
        "rows": rows,
        "summary": summary,
        "manifest": manifest,
        "functional_view": view,
    }


# ---------------------------------------------------------------------------
# Artifact writing
# ---------------------------------------------------------------------------

def write_catalog(
    output_dir: Path,
    rows: Sequence[Mapping[str, str]],
    summary: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> dict[str, str]:
    """Write the three declared artifacts inside ``output_dir`` only.

    Refuses an existing non-empty output directory and re-validates every row
    with the shared ontology before a single byte is written.
    """
    out_dir = Path(output_dir)
    if out_dir.exists():
        if not out_dir.is_dir():
            raise FileExistsError(f"output path exists and is not a directory: {out_dir}")
        if any(out_dir.iterdir()):
            raise FileExistsError(
                f"refusing to write into a non-empty output directory: {out_dir}"
            )
    for row in rows:
        try:
            SCHEMA.validate_candidate_evidence(dict(row))
        except ValueError as error:
            raise EvidenceError(
                RULE_SCHEMA_VALIDATION_FAILED,
                f"{row.get('accession', '')}: {error}",
            ) from None
    extra = sorted(set().union(*(set(row) for row in rows)) - set(CATALOG_COLUMNS)) if rows else []
    if extra:
        raise EvidenceError(
            RULE_SCHEMA_VALIDATION_FAILED,
            f"catalog rows carry undeclared columns: {extra}",
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    catalog_path = out_dir / "phaded_candidate_catalog_v2.tsv"
    _write_tsv(catalog_path, CATALOG_COLUMNS, rows)
    summary_path = out_dir / "catalog_summary.json"
    summary_path.write_text(
        json.dumps(dict(summary), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    manifest_path = out_dir / "catalog_manifest.json"
    manifest_path.write_text(
        json.dumps(dict(manifest), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "catalog": str(catalog_path),
        "summary": str(summary_path),
        "manifest": str(manifest_path),
    }


def build_catalog(
    candidates_path: "str | Path",
    demotion_path: "str | Path",
    output_dir: "str | Path",
    *,
    expect_demoted: int = EXPECTED_CYS_EXPORT_DEMOTIONS,
    basis: str = "",
    basis_detail: str = "",
    dry_run: bool = False,
) -> dict[str, Any]:
    """Read the two inputs, build the catalog and write the three artifacts.

    All input paths are validated before any work; ``dry_run`` validates
    everything (including every emitted row) but writes nothing.
    """
    candidates_file = Path(candidates_path)
    demotion_file = Path(demotion_path)
    for path in (candidates_file, demotion_file):
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"input is not a regular file: {path}")
    out_dir = Path(output_dir)
    if out_dir.exists():
        if not out_dir.is_dir():
            raise FileExistsError(f"output path exists and is not a directory: {out_dir}")
        if any(out_dir.iterdir()):
            raise FileExistsError(
                f"refusing to write into a non-empty output directory: {out_dir}"
            )

    candidate_columns, candidate_rows = _read_tsv(candidates_file)
    if "accession" not in candidate_columns:
        raise EvidenceError(
            RULE_CLASSIFICATION_FAILED,
            f"candidates input has no accession column: {candidate_columns}",
        )
    demotion_columns, demotion_rows = _read_tsv(demotion_file)
    if demotion_rows and "accession" not in demotion_columns:
        raise EvidenceError(
            RULE_DEMOTION_ACCOUNTING,
            f"demotion input has no accession column: {demotion_columns}",
        )
    if demotion_rows and not any(_text(row.get("accession")) for row in demotion_rows):
        raise EvidenceError(
            RULE_DEMOTION_ACCOUNTING,
            f"demotion input {demotion_file} has no accession values",
        )

    result = build_catalog_rows(
        candidate_rows, demotion_rows,
        expect_demoted=expect_demoted, basis=basis, basis_detail=basis_detail,
    )
    manifest = dict(result["manifest"])
    manifest["inputs"] = {
        "candidates": _bound_file(candidates_file),
        "demotion": _bound_file(demotion_file),
    }
    manifest["outputs"] = {
        "catalog": str(out_dir / "phaded_candidate_catalog_v2.tsv"),
        "summary": str(out_dir / "catalog_summary.json"),
        "manifest": str(out_dir / "catalog_manifest.json"),
    }
    manifest["run_id"] = out_dir.name
    manifest["dry_run"] = bool(dry_run)

    if not dry_run:
        written = write_catalog(out_dir, result["rows"], result["summary"], manifest)
        manifest["outputs"] = written
    return {
        "run_id": manifest["run_id"],
        "dry_run": bool(dry_run),
        "counts": manifest["counts"],
        "primary_disposition_counts": result["summary"]["primary_disposition_counts"],
        "demotion_accounting": result["summary"]["demotion_accounting"],
        "summary": result["summary"],
        "manifest": manifest,
        "rows": result["rows"],
        "functional_view": result["functional_view"],
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--candidates", type=Path, required=True,
                        help="TSV of the candidate universe (one row per protein)")
    parser.add_argument("--demotion", type=Path, required=True,
                        help="TSV of the explicit demotion/deferred accession set")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="directory for the catalog, summary and manifest")
    parser.add_argument(
        "--expect-demoted", type=int, default=EXPECTED_CYS_EXPORT_DEMOTIONS,
        help=(
            "declared number of demotion rows; defaults to %d. "
            "The default basis: %s" % (
                EXPECTED_CYS_EXPORT_DEMOTIONS, EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS,
            )
        ),
    )
    parser.add_argument("--demotion-basis", default="",
                        help="free-text override of the recorded basis detail")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate every input and row without writing anything")
    return parser


def main(argv: "Sequence[str] | None" = None) -> int:
    """Build the v2 PhaDED candidate catalog (evidence-separated semantics)."""
    args = build_parser().parse_args(argv)
    try:
        result = build_catalog(
            args.candidates, args.demotion, args.output_dir,
            expect_demoted=args.expect_demoted,
            basis_detail=args.demotion_basis,
            dry_run=args.dry_run,
        )
    except (FileNotFoundError, FileExistsError, EvidenceError, ValueError) as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1
    print(json.dumps(
        {
            "run_id": result["run_id"],
            "dry_run": result["dry_run"],
            "counts": result["counts"],
            "primary_disposition_counts": result["primary_disposition_counts"],
            "demotion_accounting": result["demotion_accounting"],
            "outputs": result["manifest"]["outputs"],
        },
        ensure_ascii=False, sort_keys=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
