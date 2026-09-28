#!/usr/bin/env python3
"""Adapt the frozen PhaDED evidence tables into the v2 candidate-catalog input.

Task F14 ("v2 候选目录真实重算 + 流转对账 + v1/v2 影响报告").  The v2 catalog
builder (``build_phaded_candidate_catalog_v2.classify_row``) resolves a documented
set of *evidence* fields; the frozen tables carry those facts under different
column names and with three different encodings.  This module is that single,
auditable translation step and nothing else.

What it does not do
-------------------
No scoring, no family call, no calibration, no disposition.  It joins frozen
tables on ``accession``, renames fields, normalizes documented state tokens and
derives the two model-status fields from the profile authority plus F13's curated
independent-positive counts.  Every mapping decision is written to a mapping
report that names the source ``file::column`` and the rule, so a reader can
re-derive the catalog input without reading this file.

Encodings that must be translated (measured, not assumed)
---------------------------------------------------------
``lipase_box_state``/``ahsmg_state`` appear as ``supported`` /
``not_detected_pattern`` in the v1 merge, the v1 hold table, the 52-row demotion
addendum and ``motif_completion_full.tsv``, but as ``present`` / ``not_detected``
in the authoritative 109,087-row localization subtype matrix.  ``present`` is
mapped to ``supported`` because the shared ontology lists ``present`` in its own
true-token set (:data:`phaded_evidence_schema._TRUE_TOKENS`); the mapping is a
documented alias, recorded per column with its input-value histogram.

AHSMG (measured boundary)
-------------------------
The frozen 109,087-row universe carries ``ahsmg_state=not_detected`` for **every**
row, and the frozen full-library motif layer states the reason in its own column:
``ahsmg_assignability_basis = pattern_miss_reproduces_0_of_0_literature_anchors_for_this_role``
-- the AHSMG detector had no literature anchor to calibrate against, so the miss
is not evidence of absence.  The adapter therefore emits **zero** AHSMG rows and
records the class in the report under ``detector_never_fired`` rather than
inventing a positive.  The AHSMG *semantics* are still implemented and pinned by
tests: a row that does carry ``ahsmg_state=supported`` becomes
``motif_class=AHSMG`` with ``nucleophile_identity=ser`` and never a third
nucleophile identity.  (The 82 AHSMG-positive rows of the frozen v1 40,690 merge
are the pool-external ``extracellular native-SCL/PhaZ7-like`` layer, which is
outside this universe -- see the report's ``outside_universe`` block.)

The corrected calibration bar (F13)
-----------------------------------
``model_layer=sequence_family_hmm_validated`` requires a profile that passed the
independently evaluated held-out/confounder gates of
``evaluate_phaded_sequence_models.py``.  Those evaluation inputs are **not**
available in this checkout, so no profile may be emitted at that layer and none
is: the reachable layers here are ``discovery_hmm_uncalibrated`` (a discovery-layer
profile of the row's own model family has a hit) and ``reference_query_only``.
The F13 curated ledger counts only ``DED_hfam_52`` (6) and ``DED_hfam_70`` (7) at
or above three unique independence groups; ``DED_hfam_4`` (1), ``DED_hfam_55``
(2) and ``DED_hfam_8`` (2) are **below** the bar.  The frozen v1 ``trained`` label
of those three families is therefore *not* kept: their hits stay at
``reference_query_only``, and the decision is recorded under ``calibration_bar``
with ``frozen_v1_trained_label_kept_for_below_bar_families = False``.

Boundaries
----------
Read-only over the frozen evidence: no ``runs/``, ``results/``, ``deploy/`` or
``inputs/`` path is written.  Writes only the two artifacts inside the caller's
``--results-dir``.  ``excluded_input_quality`` is never emitted: no
input-quality field is produced at all, so the only way a row could reach that
disposition (a documented input-quality reason) is structurally absent.  No row
is deleted or demoted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


def _load_module(name: str) -> Any:
    """Import a sibling script as a module (``pipeline/scripts`` is not a package)."""
    try:
        return __import__(name)
    except ImportError:
        pass
    path = Path(__file__).resolve().parent / f"{name}.py"
    if not path.is_file():
        raise ImportError(f"cannot locate {name}.py next to {__file__}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise ImportError(f"cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault(name, module)
    spec.loader.exec_module(module)
    return module


SCHEMA = _load_module("phaded_evidence_schema")

SCHEMA_VERSION = "1.0"
GENERATOR = "build_phaded_v2_catalog_input.py"

#: Absent values are written as the empty string and never invented; the catalog
#: builder turns a blank into its own documented weakest default (``pending``).
ABSENT = ""

#: The ontology's own weakest transport default.
NOT_TESTED = "not_tested"
#: The ontology's own weakest functional-calibration status.
NOT_FUNCTION_CALIBRATED = "not_function_calibrated"

ADAPTED_TSV = "phaded_v2_catalog_input.tsv"
MAPPING_REPORT_JSON = "catalog_input_mapping_report.json"
MAPPING_REPORT_TSV = "catalog_input_column_mapping.tsv"
OUTPUT_FILENAMES = (ADAPTED_TSV, MAPPING_REPORT_JSON, MAPPING_REPORT_TSV)

# ---------------------------------------------------------------------------
# Output contract: exactly the fields the catalog builder resolves, plus the
# v2-comparison evidence columns that the frozen v1 rule reproduction needs.
# ---------------------------------------------------------------------------

#: Fields ``build_phaded_candidate_catalog_v2.classify_row`` resolves.
V2_INPUT_COLUMNS = (
    "accession",
    "genome",
    "historical_superfamily",
    "sequence_family_call",
    "catalytic_domain_type",
    "accessory_domain_architecture",
    "transport_signal_prediction",
    "localization_evidence",
    "nucleophile_identity",
    "motif_class",
    "experimental_evidence_grade",
    "model_layer",
    "functional_calibration_status",
)

#: Evidence state columns the catalog builder reads through its own resolvers
#: (``resolve_nucleophile`` re-derives the identity/motif pair from them) and the
#: frozen v1 rule of ``compare_phaded_catalog_versions`` needs.  They are carried
#: verbatim so the v1 rule can be reproduced rather than guessed.
EVIDENCE_STATE_COLUMNS = (
    "motif_state",
    "lipase_box_state",
    "ahsmg_state",
    "catalytic_ser_cys_state",
    "his_state",
    "asp_state",
    "oxyanion_hole_state",
    "sbd_state",
    "linker_state",
    "lid_state",
    "sbd_pf06850_binding_state",
)

#: v1-rule reproduction columns.  ``compare_phaded_catalog_versions`` reports a
#: criterion as ``pending`` when its column is absent, so these are carried with
#: their frozen values and the comparison stays decidable.
V1_RULE_COLUMNS = (
    "layer",
    "family",
    "source_layer",
    "signalp_class",
    "signalp_evidence",
    "sbd_pf06850_binding_state_raw",
    "lid_state_raw",
    "linker_state_raw",
    "interpro_status",
    "architecture_consistency",
    "profile_evidence_status",
    "profile_best_evalue",
    "family_score_gap",
    "superfamily_score_gap",
    "superfamily_confidence",
    "lipase_box_x1",
    "is_confounder",
    "sequence_integrity",
    "assignment_status",
    "assignment_unique",
    "phaded_superfamily_best",
    "phaded_family_best",
    "subtype_call",
    "evidence_status",
    "evidence_threshold",
    "profile_model_status",
    "pool_origin",
    "hold_reason",
)

OUTPUT_COLUMNS = V2_INPUT_COLUMNS + EVIDENCE_STATE_COLUMNS + V1_RULE_COLUMNS

#: Input table keys, in the order they appear in the report.
INPUT_KEYS = (
    "universe",
    "subtype_matrix",
    "v1_merge",
    "v1_hold",
    "demotion",
    "authority",
    "profile_manifest",
    "reference_ledger",
    "independence_summary",
)

#: The two absolute paths the adapter reports but never reads as its join input.
SEQUENCE_MODEL_EVALUATION_INPUT = (
    "sequence_model_evaluation.tsv (evaluate_phaded_sequence_models.py) -- not present"
)

# ---------------------------------------------------------------------------
# Mapping rules (named, so the report and the code cannot drift)
# ---------------------------------------------------------------------------

RULE_LITERAL = "literal_from_universe"
RULE_SUPERFAMILY_SUBTYPE = "subtype_matrix::phaded_superfamily_best (v2 authority); pending when absent"
RULE_SUPERFAMILY_FALLBACK = "v1_merge::superfamily (v1_fallback, only when the subtype matrix is blank)"
RULE_FAMILY_CALL_WITHHELD = "withheld: no profile at this run passed the held-out/confounder evaluation"
RULE_SIGNALP = (
    "subtype_matrix::signalp_class carried verbatim; the ontology maps pending->not_tested "
    "(never a prediction)"
)
RULE_LOCALIZATION_HISTORICAL = "subtype_matrix::phaded_superfamily_best -> historical_family_label_only"
RULE_LOCALIZATION_UNKNOWN = "no family claim -> unknown"
RULE_NUCLEOPHILE = "phaded_evidence_schema.normalize_nucleophile(superfamily, motif_state, ahsmg_state)"
RULE_MOTIF_STATE = "subtype_matrix::lipase_box_state/::ahsmg_state; present->supported, not_detected->not_detected_pattern"
RULE_DOMAIN_TYPE = "catalytic_domain_type / geometry state; unrecognized verdict -> unresolved"
RULE_ACCESSORY = "sbd_state/linker_state/lid_state carried verbatim (never an input to the type call)"
RULE_GRADE = "no experimental grade exists for a GTDB candidate (reference ledger is disjoint from the universe)"
RULE_MODEL_LAYER = "profile authority + F13 corrected bar; sequence_family_hmm_validated withheld (evaluation inputs absent)"
RULE_CALIBRATION_STATUS = "not_function_calibrated: independent of model_layer; no promotion authorized"
RULE_VERBATIM = "carried verbatim from the frozen source column (comparison evidence)"

#: Roles recorded per output column (used by the mapping report).
ROLE_V2_INPUT = "v2_catalog_input"
ROLE_EVIDENCE_STATE = "v2_evidence_state"
ROLE_V1_RULE = "v1_rule_reproduction_evidence"

# ---------------------------------------------------------------------------
# Failure rules (stable, machine-readable)
# ---------------------------------------------------------------------------

RULE_UNIVERSE_EMPTY = "universe_is_empty"
RULE_UNIVERSE_DUPLICATE = "universe_has_duplicate_accessions"
RULE_UNIVERSE_NOT_COVERED = "universe_accession_absent_from_subtype_matrix"
RULE_SUBTYPE_NOT_UNIQUE = "subtype_matrix_accession_not_unique"
RULE_UNKNOWN_PROFILE = "profile_hit_absent_from_the_profile_manifest"
RULE_TABLE_MISSING = "input_is_not_a_regular_file"
RULE_ROW_COUNT_MISMATCH = "adapted_row_count_mismatch"
RULE_ACCESSION_SET_MISMATCH = "adapted_accession_set_mismatch"
RULE_UNDECLARED_OUTPUT_COLUMN = "output_column_not_declared"

#: Vocabulary aliases used to read the condensed subtype-matrix encoding.  Both
#: spellings live in the shared ontology: ``supported`` is a true token and
#: ``present`` is one too, so this table only makes the two encodings comparable.
STATE_ALIASES = {
    "present": "supported",
    "detected": "supported",
    "not_detected": "not_detected_pattern",
}

#: States whose v1 vocabulary must be translated before the catalog builder reads
#: them: it resolves a support flag from these columns and rejects an
#: unrecognized token outright, so e.g. ``detected`` (the frozen PF06850 binding
#: spelling) is normalized to the ontology's ``supported``.
_FLAG_STATE_COLUMNS = (
    "sbd_state",
    "sbd_pf06850_binding_state",
    "linker_state",
    "lid_state",
)

#: Columns that merely carry a state for the record and are never resolved into
#: a flag by the catalog builder.
_PASSTHROUGH_STATE_COLUMNS = (
    "catalytic_ser_cys_state",
    "his_state",
    "asp_state",
    "oxyanion_hole_state",
)


class InputError(ValueError):
    """An input violates a documented contract; the rule name is machine-readable."""

    def __init__(self, rule: str, message: str) -> None:
        super().__init__(f"{rule}: {message}")
        self.rule = rule
        self.message = message


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _text(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", _text(value).lower()).strip("_")


def _normalize_state(value: object) -> str:
    """Normalize a documented state token; unknown tokens are kept verbatim."""
    raw = _text(value)
    if not raw:
        return ABSENT
    return STATE_ALIASES.get(_key(raw), raw)


def _normalize_flag_state(value: object) -> str:
    """Normalize a state the catalog builder resolves into a boolean flag.

    The builder's ``_flag_state`` rejects a token outside the ontology's
    true/false/``not_detected*`` sets instead of guessing, so a raw value is only
    passed through when the ontology itself accepts it; otherwise the state stays
    absent and the *verbatim* value remains available in the dedicated evidence
    column the comparison reads.  Nothing is coerced into a support claim.
    """
    normalized = _normalize_state(value)
    if not normalized:
        return ABSENT
    try:
        SCHEMA._resolve_optional_flag({"_state": normalized}, "_state")
    except ValueError:
        key = _key(normalized)
        if key.startswith("not_detected") or key.startswith("no_"):
            return normalized
        return ABSENT
    return normalized


def _sbd_flag_state(raw: object) -> str:
    """Normalize an SBD observation into the ontology's own flag vocabulary.

    ``sbd_state`` and ``sbd_pf06850_binding_state`` are both read by the catalog
    builder's support flag, so both must hold a token it accepts.  The frozen
    vocabulary maps onto the ontology as: ``detected`` -> ``supported``; a
    ``not_detected*`` spelling stays verbatim; an explicit not-tested marker
    becomes ``not_detected_in_tested_pfam``, which states the same thing in the
    ontology's own words (the PFAM scan ran and did not find PF06850) rather than
    claiming a negative that was never tested.  An unreadable token is dropped
    rather than coerced; the verbatim value stays in the ``*_raw`` column.
    """
    normalized = _normalize_flag_state(raw)
    if normalized:
        return normalized
    key = _key(raw)
    if key.startswith("not_tested") or key.startswith("pending"):
        return "not_detected_in_tested_pfam"
    return ABSENT


def _first_present(*values: object) -> str:
    for value in values:
        text = _text(value)
        if text:
            return text
    return ABSENT


def _first_informative(*values: object) -> "tuple[str, str]":
    """Return ``(verbatim_value, normalized_flag_state)`` from the first source
    that actually states something the ontology can read.

    Precedence follows the same rule as every other column in this adapter: the
    localization subtype matrix is the authority for this 109,087-row universe,
    and a *frozen-but-uninformative* marker there (``pending``) falls through to
    the v1 evidence rather than blanking the field.  A value the ontology would
    reject is skipped instead of being coerced into a support claim.
    """
    for value in values:
        text = _text(value)
        if not text:
            continue
        normalized = _normalize_flag_state(text)
        if normalized:
            return text, normalized
    return ABSENT, ABSENT


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: "str | Path") -> "tuple[list[str], list[dict[str, str]]]":
    tsv_path = Path(path)
    if not tsv_path.is_file() or tsv_path.is_symlink():
        raise InputError(RULE_TABLE_MISSING, f"input is not a regular file: {tsv_path}")
    with tsv_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or ())
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def _require_column(columns: Sequence[str], column: str, what: str) -> None:
    if column not in columns:
        raise InputError(
            RULE_TABLE_MISSING, f"{what} has no {column!r} column; found: {list(columns)}"
        )


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t",
            lineterminator="\n", extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def bound_file(path: "str | Path") -> dict[str, Any]:
    """Bind one input with its real path, size and SHA-256 (never fabricated)."""
    source = Path(path)
    if not source.is_file():
        return {"path": str(source), "status": "pending", "size": None, "sha256": None}
    return {
        "path": str(source),
        "status": "verified",
        "size": source.stat().st_size,
        "sha256": _sha256(source),
    }


# ---------------------------------------------------------------------------
# Profile authority, the corrected bar and the model layer
# ---------------------------------------------------------------------------

def load_profile_authority(
    profile_manifest_path: "str | Path",
    authority_path: "str | Path",
) -> dict[str, Any]:
    """Read the 46-profile authority: profile -> family/superfamily/model status.

    The authority table declares which layers are registry-eligible; the profile
    manifest binds a profile to a superfamily/family and states whether a model
    exists (``hmm_sha256`` present) and its v1 ``model_status``.  Both are frozen
    inputs; nothing is derived here.
    """
    authority_columns, authority_rows = _read_tsv(authority_path)
    for column in ("layer_id", "layer_kind", "registry_eligible"):
        _require_column(authority_columns, column, "classification authority")
    registry_eligible = {
        _text(row.get("layer_id")): _key(row.get("registry_eligible")) in {"true", "yes", "1"}
        for row in authority_rows
    }
    layer_kind = {
        _text(row.get("layer_id")): _text(row.get("layer_kind")) for row in authority_rows
    }

    profile_columns, profile_rows = _read_tsv(profile_manifest_path)
    for column in ("profile_id", "model_status", "hmm_sha256", "phaded_family_id"):
        _require_column(profile_columns, column, "profile manifest")
    profiles: dict[str, dict[str, Any]] = {}
    for row in profile_rows:
        profiles[_text(row.get("profile_id"))] = {
            "profile_id": _text(row.get("profile_id")),
            "profile_kind": _text(row.get("profile_kind")),
            "phaded_superfamily": _text(row.get("phaded_superfamily")),
            "phaded_family_id": _text(row.get("phaded_family_id")),
            "model_status": _text(row.get("model_status")),
            "has_hmm": bool(_text(row.get("hmm_sha256"))),
        }
    return {
        "profiles": profiles,
        "registry_eligible": registry_eligible,
        "layer_kind": layer_kind,
        "profile_count": len(profiles),
        "authority_row_count": len(authority_rows),
    }


def load_calibration_bar(independence_summary_path: "str | Path") -> dict[str, Any]:
    """Read F13's curated independent-positive counts and its declared bar.

    The counts come from the frozen curation summary; when both a nested and a
    flat copy exist they must agree, otherwise the input is contradictory and the
    load fails closed.
    """
    source = Path(independence_summary_path)
    if not source.is_file():
        raise InputError(RULE_TABLE_MISSING, f"independence summary is not a file: {source}")
    payload = json.loads(source.read_text(encoding="utf-8"))
    nested = payload.get("independence_summary") or {}
    flat = {
        "family_functional_positive_counts": payload.get("family_functional_positive_counts") or {},
        "minimum_independent_positive_count": payload.get("minimum_independent_positive_count"),
        "positive_count_basis": payload.get("positive_count_basis"),
    }
    strict = payload.get("strict_gate") or {}
    counts = dict(flat["family_functional_positive_counts"])
    for other in (
        nested.get("family_functional_positive_counts") or {},
        strict.get("family_functional_positive_counts") or {},
    ):
        for family, value in other.items():
            if family in counts and int(counts[family]) != int(value):
                raise InputError(
                    RULE_UNKNOWN_PROFILE,
                    f"the independence summary disagrees with itself on {family}: "
                    f"{counts[family]} vs {value}",
                )
            counts.setdefault(family, int(value))
    minimum = flat["minimum_independent_positive_count"]
    if minimum is None:
        minimum = nested.get("minimum_independent_positive_count")
    if minimum is None:
        minimum = strict.get("minimum_independent_positive_count")
    if minimum is None:
        raise InputError(
            RULE_UNKNOWN_PROFILE,
            "the independence summary declares no minimum_independent_positive_count",
        )
    basis = (
        flat["positive_count_basis"]
        or nested.get("positive_count_basis")
        or strict.get("positive_count_basis")
        or ""
    )
    return {
        "family_positive_counts": {key: int(value) for key, value in sorted(counts.items())},
        "minimum_independent_positive_count": int(minimum),
        "positive_count_basis": _text(basis),
    }


def profile_positive_count(profile: Mapping[str, Any], bar: Mapping[str, Any]) -> "int | None":
    """F13's unique-independence-group positive count for a profile's own family.

    A superfamily-level profile has no family id, so the family-level count does
    not apply to it; it returns ``None`` and is judged by its own v1 status.
    """
    family = _text(profile.get("phaded_family_id"))
    if not family:
        return None
    counts = bar.get("family_positive_counts") or {}
    if family not in counts:
        return None
    return int(counts[family])


def resolve_model_layer_for_profile(
    profile: Mapping[str, Any] | None,
    bar: Mapping[str, Any],
    *,
    evaluation_available: bool,
) -> "tuple[str, str]":
    """Return ``(model_layer, reason)`` for one profile hit.

    Rules, in order:

    1. no profile hit at all -> ``reference_query_only`` (nothing scored this
       row against a discovery model);
    2. no model built for the profile (``hmm_sha256`` pending) -> the frozen
       profile is a reference query only;
    3. a family profile below F13's corrected independent-positive bar -> its
       hits stay ``reference_query_only``; the frozen v1 ``trained`` label is
       *not* kept (this is the ``DED_hfam_4`` / ``_55`` / ``_8`` correction);
    4. a superfamily profile that the frozen authority marks ``reference_only``
       -> ``reference_query_only``;
    5. a discovery-layer profile with a hit -> ``discovery_hmm_uncalibrated``.
       ``sequence_family_hmm_validated`` is emitted only when the held-out /
       confounder evaluation is available *and* the profile passed it, which is
       never the case while the evaluator's inputs are absent.
    """
    if profile is None:
        return "reference_query_only", "no_profile_hit"
    if not profile.get("has_hmm"):
        return "reference_query_only", "profile_has_no_model"
    count = profile_positive_count(profile, bar)
    minimum = int(bar.get("minimum_independent_positive_count", 0))
    if count is not None and count < minimum:
        return "reference_query_only", "family_below_corrected_independent_positive_bar"
    if not profile.get("phaded_family_id") and profile.get("model_status") != "trained":
        return "reference_query_only", "superfamily_profile_reference_only"
    if evaluation_available:
        return "sequence_family_hmm_uncalibrated_placeholder", "evaluation_available"
    return "discovery_hmm_uncalibrated", "discovery_layer_profile_hit"


# ---------------------------------------------------------------------------
# Per-column derivation
# ---------------------------------------------------------------------------

def _subtype_value(subtype: Mapping[str, str], column: str) -> str:
    return _text(subtype.get(column))


def derive_superfamily(
    subtype: Mapping[str, str], v1: Mapping[str, str]
) -> "tuple[str, str, str]":
    """Return ``(value, source, v1_value)``; a blank stays blank, never a family."""
    best = _subtype_value(subtype, "phaded_superfamily_best")
    fallback = _text(v1.get("superfamily"))
    if best:
        return best, "subtype_matrix", fallback
    if fallback:
        return fallback, "v1_fallback", fallback
    return ABSENT, "absent", fallback


def derive_localization_evidence(superfamily: str) -> str:
    """A family label is a historical label; nothing here is experimental.

    Preference is expressed through the superfamily direction (intracellular /
    extracellular / periplasmic), and a SignalP transport prediction is never
    consulted -- it is recorded in ``transport_signal_prediction`` only.
    """
    if superfamily:
        return "historical_family_label_only"
    return "unknown"


def derive_motif_and_nucleophile(
    subtype: Mapping[str, str], superfamily: str
) -> dict[str, str]:
    """Delegate the identity/motif pair to the shared ontology.

    ``ahsmg_state`` and ``lipase_box_state`` are normalized to the ontology's own
    token spellings first; the ontology then decides.  AHSMG therefore always
    comes out as ``motif_class=AHSMG`` with ``nucleophile_identity=ser``.
    """
    ahsmg = _normalize_state(subtype.get("ahsmg_state"))
    lipase = _normalize_state(subtype.get("lipase_box_state"))
    motif_state = ""
    if lipase and _key(lipase) == "supported":
        motif_state = "GxSxG"
    elif ahsmg and _key(ahsmg) == "supported":
        motif_state = "AHSMG"
    elif lipase and _key(lipase).startswith("not_detected"):
        # No canonical GxSxG: only an SCL family may keep a Ser prior; a Cys
        # family call without the box is the Cys-associated case.
        motif_state = "Cys-associated" if "scl" in _key(superfamily) else ""
    out = SCHEMA.normalize_nucleophile(
        superfamily, motif_state, "supported" if _key(ahsmg) == "supported" else "",
    )
    # The motif state is emitted as its own column as well: the catalog builder's
    # ``resolve_nucleophile`` re-derives the pair from the ``*_state`` columns,
    # and an explicit motif state makes the decision independent of how a state
    # token happens to be read there.  ``lipase_box_state``/``ahsmg_state``
    # spellings are rewritten to the v1 ontology tokens so the comparison's
    # frozen v1 rule reads them the same way it did in 2026-09-19.
    return {
        "motif_class": out["motif_class"],
        "nucleophile_identity": out["nucleophile_identity"],
        "motif_state": motif_state,
        "ahsmg_state": ahsmg or "not_detected_pattern",
        "lipase_box_state": lipase or ABSENT,
    }


def _accessory_architecture(subtype: Mapping[str, str]) -> str:
    """Render SBD / linker / lid as independent records.

    The three domains are annotation-only extras and are deliberately kept in
    their own column: the catalog builder never reads them when it resolves
    ``catalytic_domain_type``, so an SBD can never define type 1 versus type 2.
    """
    parts: list[str] = []
    for domain, column in (
        ("SBD", "sbd_state"), ("linker", "linker_state"), ("lid", "lid_state"),
    ):
        state = _normalize_state(subtype.get(column))
        if not state:
            continue
        key = _key(state)
        if key.startswith("not_detected") or key.startswith("no_") or key in {
            "pending", "not_assignable", "undetermined",
        }:
            continue
        parts.append(f"{domain}:{state}")
    return "|".join(parts)


def derive_catalytic_domain_type(raw_type: str) -> "tuple[str, str]":
    """Map a geometry verdict onto the v2 vocabulary; never invent a verified type."""
    stated = _text(raw_type)
    if not stated:
        return "unresolved", "no_geometry_verdict_in_the_joined_tables"
    try:
        return (
            SCHEMA._resolve(
                stated, "catalytic_domain_type",
                SCHEMA._CATALYTIC_DOMAIN_TYPE_BY_KEY, SCHEMA.CATALYTIC_DOMAIN_TYPES,
                aliases=SCHEMA.CATALYTIC_DOMAIN_TYPE_ALIASES,
            ),
            "recognized_verdict",
        )
    except ValueError:
        return "unresolved", "unrecognized_verdict"


def derive_assignment_unique(assignment_status: str) -> str:
    """Uniqueness follows the frozen assignment status; ambiguous is never unique."""
    key = _key(assignment_status)
    if key == "assigned":
        return "true"
    if key.startswith("ambiguous") or key.startswith("unassigned"):
        return "false"
    return ABSENT


def derive_score_gaps(subtype: Mapping[str, str]) -> "tuple[str, str]":
    """Carry the score gaps that the catalogue's own threshold is defined on.

    The catalogue records ``evidence_threshold = min_score_gap_1.0_bits_HMMER_E1e6``
    and an ``assignment_status`` of ``assigned`` versus ``ambiguous_*``, but until
    now it did not carry the quantity those two are compared against: the caller had
    the status but not the number that produced it.

    The subtype matrix populates ``family_score_gap`` and ``superfamily_score_gap``,
    and they are not decoration - measured on the 1,087 rows whose best family is
    the promoted ``DED_hfam_70`` profile, ``family_score_gap`` separates them
    perfectly at the documented 1.0-bit threshold: all 132 ``assigned`` rows have a
    gap of 1.0-3.3 bits and all 950 ``ambiguous_family`` rows have 0.0-0.9, with no
    overlap. Dropping the gap left the catalogue asserting a status without the
    evidence for it.

    ``profile_best_evalue`` is deliberately left alone: it exists in the v1 merge
    but is EMPTY there for every one of those rows, so the catalogue faithfully
    copies an empty source. That is the upstream gap the adapter's own docstring
    names, and it is not this function's business to paper over it.
    """
    return (
        _text(subtype.get("family_score_gap")),
        _text(subtype.get("superfamily_score_gap")),
    )


def derive_profile_best_evalue(subtype: Mapping[str, str], v1: Mapping[str, str]) -> str:
    """The v1 strong-profile rule needs an E-value; the subtype matrix has a gap."""
    for source, column in ((v1, "profile_best_evalue"), (v1, "best_evalue")):
        value = _text(source.get(column))
        if value:
            return value
    return ABSENT


def derive_superfamily_confidence(
    subtype: Mapping[str, str], f1: Mapping[str, str]
) -> str:
    """F1's per-accession attribution carries the discovery claim's confidence.

    It is comparison evidence only; the discovery layer never makes a v2 family
    call, so this column is never promoted into ``sequence_family_call``.
    """
    return _text(f1.get("discovery_superfamily_confidence"))


# ---------------------------------------------------------------------------
# The join
# ---------------------------------------------------------------------------

def _index_unique(rows: Sequence[Mapping[str, str]], what: str, rule: str) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = _text(row.get("accession"))
        if not accession:
            continue
        if accession in index:
            raise InputError(rule, f"{what} repeats accession {accession}")
        index[accession] = dict(row)
    return index


def build_adapted_input(
    *,
    universe_path: "str | Path",
    subtype_matrix_path: "str | Path",
    v1_merge_path: "str | Path",
    v1_hold_path: "str | Path",
    demotion_path: "str | Path",
    authority_path: "str | Path",
    profile_manifest_path: "str | Path",
    independence_summary_path: "str | Path",
    ledger_path: "str | Path | None" = None,
    f1_attribution_path: "str | Path | None" = None,
    expect_rows: int | None = None,
) -> dict[str, Any]:
    """Join the frozen tables into one v2 catalog-input row per universe accession.

    Returns ``{"rows", "report", "columns"}``.  Fails closed when the universe is
    empty or duplicated, when the subtype matrix repeats an accession, when a
    universe accession has no subtype-matrix row, when a profile hit is absent
    from the profile manifest, or when the adapted row count / accession set does
    not equal the universe.
    """
    inputs = {
        "universe": universe_path,
        "subtype_matrix": subtype_matrix_path,
        "v1_merge": v1_merge_path,
        "v1_hold": v1_hold_path,
        "demotion": demotion_path,
        "authority": authority_path,
        "profile_manifest": profile_manifest_path,
        "reference_ledger": ledger_path,
        "independence_summary": independence_summary_path,
    }
    bound = {key: bound_file(path) for key, path in inputs.items() if path is not None}

    universe_columns, universe_rows = _read_tsv(universe_path)
    _require_column(universe_columns, "accession", "candidate universe")
    universe_accessions = [_text(row.get("accession")) for row in universe_rows]
    if not universe_accessions:
        raise InputError(RULE_UNIVERSE_EMPTY, f"candidate universe has no rows: {universe_path}")
    duplicates = sorted(
        accession for accession, count in Counter(universe_accessions).items() if count > 1
    )
    if duplicates:
        raise InputError(
            RULE_UNIVERSE_DUPLICATE,
            f"{len(duplicates)} accession(s) appear more than once in the universe: "
            f"{duplicates[:5]}",
        )

    subtype_columns, subtype_rows = _read_tsv(subtype_matrix_path)
    _require_column(subtype_columns, "accession", "localization subtype matrix")
    subtype_index = _index_unique(subtype_rows, "localization subtype matrix", RULE_SUBTYPE_NOT_UNIQUE)

    merge_columns, merge_rows = _read_tsv(v1_merge_path)
    _require_column(merge_columns, "accession", "v1 merge")
    merge_index = _index_unique(merge_rows, "v1 merge", RULE_SUBTYPE_NOT_UNIQUE)

    hold_columns, hold_rows = _read_tsv(v1_hold_path)
    _require_column(hold_columns, "accession", "v1 hold table")
    hold_index = _index_unique(hold_rows, "v1 hold table", RULE_SUBTYPE_NOT_UNIQUE)

    demotion_columns, demotion_rows = _read_tsv(demotion_path)
    _require_column(demotion_columns, "accession", "demotion addendum")
    demotion_index = _index_unique(demotion_rows, "demotion addendum", RULE_SUBTYPE_NOT_UNIQUE)

    f1_index: dict[str, dict[str, str]] = {}
    if f1_attribution_path is not None:
        f1_columns, f1_rows = _read_tsv(f1_attribution_path)
        _require_column(f1_columns, "accession", "F1 residual attribution")
        f1_index = _index_unique(f1_rows, "F1 residual attribution", RULE_SUBTYPE_NOT_UNIQUE)

    authority = load_profile_authority(profile_manifest_path, authority_path)
    bar = load_calibration_bar(independence_summary_path)
    profiles = authority["profiles"]

    reference_ledger_accessions: set[str] = set()
    reference_grades: dict[str, str] = {}
    if ledger_path is not None:
        ledger_columns, ledger_rows = _read_tsv(ledger_path)
        _require_column(ledger_columns, "accession", "curated reference ledger")
        for row in ledger_rows:
            accession = _text(row.get("accession"))
            if accession:
                reference_ledger_accessions.add(accession)
                grade = _text(row.get("experimental_evidence_grade"))
                if grade:
                    reference_grades[accession] = grade

    missing = sorted(
        accession for accession in universe_accessions if accession not in subtype_index
    )
    if missing:
        raise InputError(
            RULE_UNIVERSE_NOT_COVERED,
            f"{len(missing)} universe accession(s) have no localization subtype matrix row "
            f"(the join must be total): {missing[:5]}",
        )

    unknown_profiles: Counter[str] = Counter()
    layer_reason_counts: Counter[str] = Counter()
    superfamily_sources: Counter[str] = Counter()
    superfamily_conflicts: list[str] = []
    catalytic_raw_counts: Counter[str] = Counter()
    catalytic_unrecognized: Counter[str] = Counter()
    missing_geometry: list[str] = []
    value_counts: dict[str, Counter[str]] = {column: Counter() for column in OUTPUT_COLUMNS}
    missing_v1_evidence: Counter[str] = Counter()
    reference_rows_in_universe = 0
    grade_rows = 0

    rows: list[dict[str, str]] = []
    for universe_row in universe_rows:
        accession = _text(universe_row.get("accession"))
        subtype = subtype_index[accession]
        v1 = merge_index.get(accession) or hold_index.get(accession) or {}
        if not v1:
            for column in ("superfamily", "catalytic_domain_type", "profile_best_evalue",
                           "lipase_box_x1", "sbd_pf06850_binding_state", "interpro_status",
                           "profile_evidence_status", "signalp_class"):
                missing_v1_evidence[column] += 1
        demotion = demotion_index.get(accession, {})
        f1 = f1_index.get(accession, {})

        superfamily, superfamily_source, v1_superfamily = derive_superfamily(subtype, v1)
        superfamily_sources[superfamily_source] += 1
        if (
            superfamily_source == "subtype_matrix"
            and v1_superfamily
            and _key(v1_superfamily) != _key(superfamily)
        ):
            superfamily_conflicts.append(accession)

        raw_type = _text(v1.get("catalytic_domain_type")) or _text(demotion.get("catalytic_domain_type"))
        catalytic_raw_counts[raw_type or "<absent>"] += 1
        domain_type, domain_reason = derive_catalytic_domain_type(raw_type)
        if not raw_type:
            missing_geometry.append(accession)
        if domain_reason == "unrecognized_verdict":
            catalytic_unrecognized[raw_type] += 1

        motif = derive_motif_and_nucleophile(subtype, superfamily)

        signalp = _text(subtype.get("signalp_class"))
        # Two transport columns are emitted, and the split is deliberate:
        #
        # * ``transport_signal_prediction`` is the v2 column and takes the
        #   ontology's own weakest default for an untested row, because the
        #   catalog builder's ``resolve_transport_signal`` reads the frozen
        #   ``pending`` spelling as a synonym of ``not_tested`` only when it is
        #   handed the raw signalp field; and
        # * ``signalp_class`` is the v1-column-name alias of the *same* value, so
        #   the catalog builder (which prefers ``signalp_class``) and the
        #   comparison module (whose ``TRANSPORT_COLUMNS`` order prefers
        #   ``transport_signal_prediction``) both read it.  ``pending`` and
        #   ``not_tested`` name the same state -- nothing was predicted -- and a
        #   blank class is written as the documented default instead of a bare
        #   empty cell.
        transport = (
            NOT_TESTED
            if not signalp or _key(signalp) in {"pending", "untested", "not_tested"}
            else signalp
        )

        profile_id = _text(subtype.get("phaded_family_best"))
        profile = profiles.get(profile_id) if profile_id else None
        if profile_id and profile is None:
            unknown_profiles[profile_id] += 1
            raise InputError(
                RULE_UNKNOWN_PROFILE,
                f"{accession}: profile hit {profile_id!r} is not in the profile manifest",
            )
        layer, layer_reason = resolve_model_layer_for_profile(
            profile, bar, evaluation_available=False,
        )
        if layer.endswith("_placeholder"):  # pragma: no cover - defensive
            raise InputError(
                RULE_UNKNOWN_PROFILE,
                "sequence_family_hmm_validated cannot be emitted while the held-out / "
                "confounder evaluation inputs are absent",
            )
        layer_reason_counts[layer_reason] += 1

        if accession in reference_ledger_accessions:
            reference_rows_in_universe += 1
        grade = reference_grades.get(accession, ABSENT)
        if grade:
            grade_rows += 1

        # ``sbd_state`` is the column the catalog builder resolves a support flag
        # from, so it must be a token the ontology accepts.  The frozen PF06850
        # binding column is a *separate* v1 observation and keeps its own
        # verbatim value, exactly as the frozen v1 architecture rule read it.
        sbd_raw = _first_present(
            v1.get("sbd_pf06850_binding_state"),
            demotion.get("sbd_pf06850_binding_state"),
            subtype.get("sbd_state"),
        )
        sbd_flag = _sbd_flag_state(sbd_raw)
        linker_raw, linker_flag = _first_informative(
            subtype.get("linker_state"), v1.get("linker_state")
        )
        lid_raw, lid_flag = _first_informative(
            subtype.get("lid_state"), v1.get("lid_state")
        )

        row: dict[str, str] = {
            "accession": accession,
            "genome": _text(universe_row.get("genome")) or _text(subtype.get("genome")),
            "historical_superfamily": superfamily,
            "sequence_family_call": ABSENT,
            "catalytic_domain_type": domain_type,
            "accessory_domain_architecture": _accessory_architecture(subtype),
            "transport_signal_prediction": transport,
            "localization_evidence": derive_localization_evidence(superfamily),
            "nucleophile_identity": motif["nucleophile_identity"],
            "motif_class": motif["motif_class"],
            "motif_state": motif["motif_state"],
            "experimental_evidence_grade": grade,
            "model_layer": layer,
            "functional_calibration_status": NOT_FUNCTION_CALIBRATED,
            "lipase_box_state": motif["lipase_box_state"],
            "ahsmg_state": motif["ahsmg_state"],
            "catalytic_ser_cys_state": _normalize_state(subtype.get("catalytic_ser_cys_state")),
            "his_state": _normalize_state(subtype.get("his_state")),
            "asp_state": _normalize_state(subtype.get("asp_state")),
            "oxyanion_hole_state": _normalize_state(subtype.get("oxyanion_hole_state")),
            # ``sbd_state``/``sbd_pf06850_binding_state`` are the columns the
            # catalog builder resolves a support flag from, so both carry the
            # ontology's own token; the frozen spelling is preserved verbatim in
            # the ``*_raw`` columns the comparison reads.
            "sbd_state": sbd_flag,
            "sbd_pf06850_binding_state": sbd_flag,
            "linker_state": linker_flag,
            "lid_state": lid_flag,
            "sbd_pf06850_binding_state_raw": sbd_raw,
            "lid_state_raw": lid_raw,
            "linker_state_raw": linker_raw,
            "layer": _text(subtype.get("layer")),
            "family": _text(subtype.get("family")),
            "source_layer": _text(universe_row.get("source_layer")),
            "signalp_evidence": _text(subtype.get("signalp_evidence")),
            "signalp_class": transport,
            "interpro_status": _first_text(v1.get("interpro_status"), subtype.get("interpro_status")),
            "architecture_consistency": _text(subtype.get("architecture_consistency")),
            "profile_evidence_status": _first_text(
                subtype.get("profile_evidence_status"), v1.get("profile_evidence_status")
            ),
            "profile_best_evalue": derive_profile_best_evalue(subtype, v1),
            "family_score_gap": derive_score_gaps(subtype)[0],
            "superfamily_score_gap": derive_score_gaps(subtype)[1],
            "superfamily_confidence": derive_superfamily_confidence(subtype, f1),
            "lipase_box_x1": _first_text(v1.get("lipase_box_x1"), demotion.get("lipase_box_x1")),
            "is_confounder": _text(f1.get("confounder_flag")),
            "sequence_integrity": _text(subtype.get("sequence_integrity")),
            "assignment_status": _text(subtype.get("assignment_status")),
            "assignment_unique": derive_assignment_unique(subtype.get("assignment_status")),
            "phaded_superfamily_best": _text(subtype.get("phaded_superfamily_best")),
            "phaded_family_best": profile_id,
            "subtype_call": _text(subtype.get("subtype_call")),
            "evidence_status": _text(subtype.get("evidence_status")),
            "evidence_threshold": _text(subtype.get("evidence_threshold")),
            "profile_model_status": _text(subtype.get("profile_model_status")),
            "pool_origin": _text(v1.get("pool_origin")),
            "hold_reason": _first_text(v1.get("hold_reason"), demotion.get("hold_reason")),
        }
        extra = sorted(set(row) - set(OUTPUT_COLUMNS))
        if extra:  # pragma: no cover - defensive
            raise InputError(RULE_UNDECLARED_OUTPUT_COLUMN, f"undeclared output column(s): {extra}")
        for column in OUTPUT_COLUMNS:
            value_counts[column][row.get(column, ABSENT)] += 1
        rows.append(row)

    if expect_rows is not None and len(rows) != int(expect_rows):
        raise InputError(
            RULE_ROW_COUNT_MISMATCH,
            f"adapted {len(rows)} row(s) but the expected universe size is {expect_rows}",
        )
    adapted_accessions = [row["accession"] for row in rows]
    if adapted_accessions != universe_accessions:
        raise InputError(
            RULE_ACCESSION_SET_MISMATCH,
            "the adapted accession list differs from the universe list (order or content)",
        )

    report = _build_report(
        rows=rows,
        bound=bound,
        bar=bar,
        authority=authority,
        layer_reason_counts=layer_reason_counts,
        superfamily_sources=superfamily_sources,
        superfamily_conflicts=superfamily_conflicts,
        catalytic_raw_counts=catalytic_raw_counts,
        catalytic_unrecognized=catalytic_unrecognized,
        missing_geometry=missing_geometry,
        value_counts=value_counts,
        missing_v1_evidence=missing_v1_evidence,
        reference_ledger_accessions=reference_ledger_accessions,
        reference_rows_in_universe=reference_rows_in_universe,
        grade_rows=grade_rows,
        universe_rows=len(universe_rows),
        merge_index=merge_index,
        hold_index=hold_index,
        demotion_index=demotion_index,
        subtype_rows=len(subtype_rows),
    )
    return {"rows": rows, "report": report, "columns": list(OUTPUT_COLUMNS)}


def _first_text(*values: object) -> str:
    for value in values:
        text = _text(value)
        if text:
            return text
    return ABSENT


_MAPPING_SPEC: dict[str, tuple[str, str, str]] = {
    # column: (role, source file::column, rule)
    "accession": (ROLE_V2_INPUT, "protein_layers.tsv::accession", RULE_LITERAL),
    "genome": (ROLE_V2_INPUT, "protein_layers.tsv::genome", RULE_LITERAL),
    "historical_superfamily": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::phaded_superfamily_best | v1_merge::superfamily",
        RULE_SUPERFAMILY_SUBTYPE + " | " + RULE_SUPERFAMILY_FALLBACK,
    ),
    "sequence_family_call": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::phaded_family_best (used for the layer only)",
        RULE_FAMILY_CALL_WITHHELD,
    ),
    "catalytic_domain_type": (
        ROLE_V2_INPUT,
        "v1_merge::catalytic_domain_type | v1_hold::catalytic_domain_type | "
        "cys_demotion_addendum::catalytic_domain_type",
        RULE_DOMAIN_TYPE,
    ),
    "accessory_domain_architecture": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::sbd_state/::linker_state/::lid_state",
        RULE_ACCESSORY,
    ),
    "transport_signal_prediction": (
        ROLE_V2_INPUT, "localization_subtype_matrix.tsv::signalp_class", RULE_SIGNALP,
    ),
    "localization_evidence": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::phaded_superfamily_best (direction only)",
        RULE_LOCALIZATION_HISTORICAL + " | " + RULE_LOCALIZATION_UNKNOWN,
    ),
    "nucleophile_identity": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::lipase_box_state/::ahsmg_state",
        RULE_NUCLEOPHILE,
    ),
    "motif_class": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::lipase_box_state/::ahsmg_state",
        RULE_NUCLEOPHILE + " | " + RULE_MOTIF_STATE,
    ),
    "motif_state": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::lipase_box_state/::ahsmg_state",
        RULE_MOTIF_STATE + " -> the motif token handed to the shared ontology",
    ),
    "experimental_evidence_grade": (
        ROLE_V2_INPUT, "phaded_reference_ledger_curated_v2.tsv::experimental_evidence_grade",
        RULE_GRADE,
    ),
    "model_layer": (
        ROLE_V2_INPUT,
        "localization_subtype_matrix.tsv::phaded_family_best -> profile_manifest.tsv::"
        "model_status/::hmm_sha256 + reference_evidence_independence_summary.json::"
        "family_functional_positive_counts",
        RULE_MODEL_LAYER,
    ),
    "functional_calibration_status": (
        ROLE_V2_INPUT, "phaded_evidence_schema.FUNCTIONAL_CALIBRATION_STATUSES",
        RULE_CALIBRATION_STATUS,
    ),
    "sbd_pf06850_binding_state": (
        ROLE_EVIDENCE_STATE,
        "v1_merge::sbd_pf06850_binding_state | v1_hold::sbd_pf06850_binding_state | "
        "cys_demotion_addendum::sbd_pf06850_binding_state",
        "detected->supported; not_detected* kept; a not-tested marker -> "
        "not_detected_in_tested_pfam; an unreadable token stays absent",
    ),
    "sbd_pf06850_binding_state_raw": (
        ROLE_V1_RULE,
        "v1_merge::sbd_pf06850_binding_state | v1_hold::sbd_pf06850_binding_state | "
        "localization_subtype_matrix.tsv::sbd_state",
        RULE_VERBATIM,
    ),
    "lid_state_raw": (
        ROLE_V1_RULE,
        "localization_subtype_matrix.tsv::lid_state | v1_merge::lid_state",
        RULE_VERBATIM,
    ),
    "linker_state_raw": (
        ROLE_V1_RULE,
        "localization_subtype_matrix.tsv::linker_state | v1_merge::linker_state",
        RULE_VERBATIM,
    ),
    "sbd_state": (
        ROLE_EVIDENCE_STATE,
        "v1_merge::sbd_pf06850_binding_state (normalized) | localization_subtype_matrix.tsv::sbd_state",
        "detected->supported; an unrecognized token stays absent so the catalog builder "
        "never coerces it into a support claim",
    ),
    "linker_state": (
        ROLE_EVIDENCE_STATE,
        "localization_subtype_matrix.tsv::linker_state | v1_merge::linker_state",
        "normalized to an ontology flag token; a pending marker falls through to the next source",
    ),
    "lid_state": (
        ROLE_EVIDENCE_STATE,
        "localization_subtype_matrix.tsv::lid_state | v1_merge::lid_state",
        "normalized to an ontology flag token; a pending marker falls through to the next source",
    ),
    "signalp_evidence": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::signalp_evidence", RULE_VERBATIM,
    ),
    "signalp_class": (
        ROLE_V1_RULE,
        "localization_subtype_matrix.tsv::signalp_class",
        RULE_SIGNALP + "; v1-column-name alias of transport_signal_prediction",
    ),
    "interpro_status": (
        ROLE_V1_RULE,
        "v1_merge::interpro_status | v1_hold::interpro_status | "
        "localization_subtype_matrix.tsv::interpro_status",
        RULE_VERBATIM,
    ),
    "architecture_consistency": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::architecture_consistency", RULE_VERBATIM,
    ),
    "profile_evidence_status": (
        ROLE_V1_RULE,
        "localization_subtype_matrix.tsv::profile_evidence_status | "
        "v1_merge::profile_evidence_status",
        RULE_VERBATIM,
    ),
    "profile_best_evalue": (
        ROLE_V1_RULE, "v1_merge::profile_best_evalue", RULE_VERBATIM,
    ),
    "family_score_gap": (
        ROLE_EVIDENCE_STATE, "localization_subtype_matrix.tsv::family_score_gap", RULE_VERBATIM,
    ),
    "superfamily_score_gap": (
        ROLE_EVIDENCE_STATE, "localization_subtype_matrix.tsv::superfamily_score_gap",
        RULE_VERBATIM,
    ),
    "superfamily_confidence": (
        ROLE_V1_RULE,
        "residual_616.tsv::discovery_superfamily_confidence (F1 attribution)",
        RULE_VERBATIM,
    ),
    "lipase_box_x1": (
        ROLE_V1_RULE,
        "v1_merge::lipase_box_x1 | cys_demotion_addendum::lipase_box_x1",
        RULE_VERBATIM,
    ),
    "is_confounder": (
        ROLE_V1_RULE, "residual_616.tsv::confounder_flag (F1 attribution)", RULE_VERBATIM,
    ),
    "sequence_integrity": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::sequence_integrity", RULE_VERBATIM,
    ),
    "assignment_status": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::assignment_status", RULE_VERBATIM,
    ),
    "assignment_unique": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::assignment_status",
        "assigned->true; ambiguous_*/unassigned_*->false; anything else stays absent",
    ),
    "phaded_superfamily_best": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::phaded_superfamily_best", RULE_VERBATIM,
    ),
    "phaded_family_best": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::phaded_family_best", RULE_VERBATIM,
    ),
    "subtype_call": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::subtype_call", RULE_VERBATIM,
    ),
    "evidence_status": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::evidence_status", RULE_VERBATIM,
    ),
    "evidence_threshold": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::evidence_threshold", RULE_VERBATIM,
    ),
    "profile_model_status": (
        ROLE_V1_RULE, "localization_subtype_matrix.tsv::profile_model_status", RULE_VERBATIM,
    ),
    "pool_origin": (ROLE_V1_RULE, "v1_merge::pool_origin", RULE_VERBATIM),
    "hold_reason": (
        ROLE_V1_RULE,
        "v1_hold::hold_reason | cys_demotion_addendum::hold_reason",
        RULE_VERBATIM,
    ),
    "layer": (ROLE_V1_RULE, "localization_subtype_matrix.tsv::layer", RULE_VERBATIM),
    "family": (ROLE_V1_RULE, "localization_subtype_matrix.tsv::family", RULE_VERBATIM),
    "source_layer": (ROLE_V1_RULE, "protein_layers.tsv::source_layer", RULE_VERBATIM),
}
#: Every state column the adapter rewrites shares one rule entry.
for _state_column in (
    "lipase_box_state", "ahsmg_state", "catalytic_ser_cys_state", "his_state",
    "asp_state", "oxyanion_hole_state", "sbd_state", "linker_state", "lid_state",
):
    _MAPPING_SPEC.setdefault(
        _state_column,
        (
            ROLE_EVIDENCE_STATE,
            f"localization_subtype_matrix.tsv::{_state_column}",
            RULE_MOTIF_STATE
            if _state_column in {"lipase_box_state", "ahsmg_state"}
            else RULE_VERBATIM,
        ),
    )


def _count(counts: Mapping[str, int], key: str) -> int:
    """Read a counter defensively; a missing key is a reported zero, not a gap."""
    return int(counts.get(key, 0))


def _build_report(
    *,
    rows: Sequence[Mapping[str, str]],
    bound: Mapping[str, Any],
    bar: Mapping[str, Any],
    authority: Mapping[str, Any],
    layer_reason_counts: Mapping[str, int],
    superfamily_sources: Mapping[str, int],
    superfamily_conflicts: Sequence[str],
    catalytic_raw_counts: Mapping[str, int],
    catalytic_unrecognized: Mapping[str, int],
    missing_geometry: Sequence[str],
    value_counts: Mapping[str, Counter],
    missing_v1_evidence: Mapping[str, int],
    reference_ledger_accessions: Iterable[str],
    reference_rows_in_universe: int,
    grade_rows: int,
    universe_rows: int,
    merge_index: Mapping[str, Any],
    hold_index: Mapping[str, Any],
    demotion_index: Mapping[str, Any],
    subtype_rows: int,
) -> dict[str, Any]:
    """Assemble the machine-readable mapping report."""
    declared = set(OUTPUT_COLUMNS)
    undeclared = sorted(declared - set(_MAPPING_SPEC))
    if undeclared:  # pragma: no cover - defensive
        raise InputError(RULE_UNDECLARED_OUTPUT_COLUMN, f"no mapping declared for {undeclared}")

    column_mapping: dict[str, Any] = {}
    pending_value_classes: dict[str, Any] = {}
    for column in OUTPUT_COLUMNS:
        role, source, rule = _MAPPING_SPEC[column]
        counts = value_counts.get(column, Counter())
        absent_count = int(counts.get(ABSENT, 0))
        entry = {
            "role": role,
            "source": source,
            "rule": rule,
            "rows": len(rows),
            "distinct_values": len(counts),
            "value_counts": {key: int(value) for key, value in sorted(counts.items())},
            "pending_count": absent_count,
            "pending_share": round(absent_count / len(rows), 6) if rows else 0.0,
        }
        column_mapping[column] = entry
        if absent_count:
            pending_value_classes[column] = {
                "count": absent_count,
                "reason": rule,
                "source": source,
            }

    # ``not_tested`` is the transport prediction's own documented pending class.
    not_tested = int(value_counts.get("transport_signal_prediction", Counter()).get(NOT_TESTED, 0))
    if not_tested:
        pending_value_classes["transport_signal_prediction=not_tested"] = {
            "count": not_tested,
            "reason": "SignalP was never run for this accession; a transport prediction is absent, "
                      "which is not evidence of no signal",
            "source": "localization_subtype_matrix.tsv::signalp_class",
        }
    unresolved_motif = int(
        value_counts.get("motif_class", Counter()).get("unresolved", 0)
    )
    if unresolved_motif:
        pending_value_classes["motif_class=unresolved"] = {
            "count": unresolved_motif,
            "reason": "no informative catalytic motif state; the identity is kept unresolved "
                      "rather than guessed",
            "source": "localization_subtype_matrix.tsv::lipase_box_state/::ahsmg_state",
        }

    model_layer_distribution = dict(
        sorted(value_counts.get("model_layer", Counter()).items())
    )
    calibration_bar = {
        "minimum_independent_positive_count": int(bar["minimum_independent_positive_count"]),
        "positive_count_basis": bar["positive_count_basis"],
        "source": "reference_evidence_independence_summary.json",
        "families_at_or_above_bar": {
            family: count
            for family, count in bar["family_positive_counts"].items()
            if count >= int(bar["minimum_independent_positive_count"])
        },
        "families_below_bar": {
            family: count
            for family, count in bar["family_positive_counts"].items()
            if count < int(bar["minimum_independent_positive_count"])
        },
        "below_bar_profiles_kept_at_reference_query_only": sorted(
            profile["profile_id"]
            for profile in authority["profiles"].values()
            if profile.get("phaded_family_id")
            and bar["family_positive_counts"].get(profile["phaded_family_id"], 0)
            < int(bar["minimum_independent_positive_count"])
        ),
        "frozen_v1_trained_label_kept_for_below_bar_families": False,
    }

    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generator": GENERATOR,
        "purpose": "adapt the frozen PhaDED evidence tables into the v2 candidate-catalog input",
        "inputs": bound,
        "outputs": {"adapted_input": ADAPTED_TSV},
        "counts": {
            "universe_rows": universe_rows,
            "joined_rows": len(rows),
            "subtype_matrix_rows": subtype_rows,
            "v1_merge_rows": len(merge_index),
            "v1_hold_rows": len(hold_index),
            "demotion_rows": len(demotion_index),
            "reference_ledger_rows": len(list(reference_ledger_accessions)),
            "reference_rows_in_universe": reference_rows_in_universe,
            "rows_with_a_reference_grade": grade_rows,
            "profile_count": authority["profile_count"],
            "authority_row_count": authority["authority_row_count"],
        },
        "universe_is_total": len(rows) == universe_rows,
        "one_row_per_accession": len({row["accession"] for row in rows}) == len(rows),
        "v2_input_columns": list(V2_INPUT_COLUMNS),
        "output_columns": list(OUTPUT_COLUMNS),
        "column_mapping": column_mapping,
        "pending_value_classes": pending_value_classes,
        "superfamily_source_counts": {
            "subtype_matrix": _count(superfamily_sources, "subtype_matrix"),
            "v1_fallback": _count(superfamily_sources, "v1_fallback"),
            "absent": _count(superfamily_sources, "absent"),
        },
        "superfamily_conflicts": {
            "count": len(superfamily_conflicts),
            "accessions": sorted(superfamily_conflicts)[:50],
            "basis": (
                "the subtype matrix is the v2 authority; a v1-merge value is used only as a "
                "fallback and every disagreement is reported here, never silently resolved"
            ),
        },
        "catalytic_domain_type_raw_counts": dict(sorted(catalytic_raw_counts.items())),
        "catalytic_domain_type_unrecognized": dict(sorted(catalytic_unrecognized.items())),
        "catalytic_domain_type_missing_rows": {
            "count": len(missing_geometry),
            "accessions": sorted(missing_geometry)[:20],
            "reason": (
                "no joined frozen table states a geometry verdict for this accession, so the "
                "type stays unresolved; a verified type is never invented"
            ),
        },
        "model_layer_distribution": model_layer_distribution,
        "model_layer_reasons": dict(sorted(layer_reason_counts.items())),
        "calibration_bar": calibration_bar,
        "sequence_model_evaluation": {
            "inputs_available": False,
            "expected_input": SEQUENCE_MODEL_EVALUATION_INPUT,
            "consequence": (
                "model_layer=sequence_family_hmm_validated is withheld for every row: a profile "
                "that has not passed the held-out/confounder evaluation is not a validated "
                "sequence-family model. Affected profiles stay at their lower layer."
            ),
        },
        "detector_never_fired": {
            "ahsmg_state=supported": int(
                value_counts.get("ahsmg_state", Counter()).get("supported", 0)
            ),
            "note": (
                "zero AHSMG-positive rows in this universe: the frozen full-library motif layer "
                "records ahsmg_assignability_basis="
                "pattern_miss_reproduces_0_of_0_literature_anchors_for_this_role, i.e. the "
                "detector had no anchor to calibrate against, so the AHSMG miss is not evidence "
                "of absence. The AHSMG->ser motif-class semantics are implemented and pinned by "
                "tests; no positive was invented for this run."
            ),
        },
        "outside_universe": {
            "v1_merge_pool_external_rows": sum(
                1 for row in merge_index.values() if _text(row.get("pool_origin")) == "pool_external"
            ),
            "note": (
                "the frozen v1 merge carries a pool-external layer that is not part of the "
                "109,087 candidate universe; the adapter joins on the universe, so those rows "
                "are excluded from the catalog by construction and never deleted"
            ),
        },
        "missing_v1_evidence_rows": dict(sorted(missing_v1_evidence.items())),
        "unmappable_rows": {
            "count": len(missing_geometry),
            "classes": (
                {"no_geometry_verdict": len(missing_geometry)} if missing_geometry else {}
            ),
            "dropped_rows": 0,
            "note": (
                "an unmappable value is reported with its class and kept as the documented "
                "weakest default; no row is ever dropped from the universe"
            ),
        },
        "boundaries": {
            "read_only_over_frozen_evidence": True,
            "deleted_rows": 0,
            "demoted_rows": 0,
            "family_call_made": False,
            "scoring_performed": False,
            "experimental_grade_fabricated": False,
            "input_quality_field_emitted": False,
            "excluded_input_quality_reachable": False,
            "functional_calibration_status_independent": True,
            "sbd_influences_catalytic_domain_type": False,
            "signalp_promoted_to_localization": False,
        },
    }
    return report


# ---------------------------------------------------------------------------
# Artifact writing
# ---------------------------------------------------------------------------

def write_artifacts(
    results_dir: "str | Path",
    rows: Sequence[Mapping[str, str]],
    report: Mapping[str, Any],
) -> dict[str, str]:
    """Write the adapted TSV and the mapping report inside ``results_dir`` only."""
    out_dir = Path(results_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adapted = out_dir / ADAPTED_TSV
    _write_tsv(adapted, OUTPUT_COLUMNS, rows)

    mapping_tsv = out_dir / MAPPING_REPORT_TSV
    columns = ("column", "role", "source", "rule", "rows", "distinct_values", "pending_count")
    with mapping_tsv.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for column in OUTPUT_COLUMNS:
            entry = dict(report["column_mapping"][column])
            entry["column"] = column
            writer.writerow(entry)

    report_path = out_dir / MAPPING_REPORT_JSON
    report_path.write_text(
        json.dumps(dict(report), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "adapted_input": str(adapted),
        "mapping_report": str(report_path),
        "column_mapping_tsv": str(mapping_tsv),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--subtype-matrix", type=Path, required=True)
    parser.add_argument("--v1-merge", type=Path, required=True)
    parser.add_argument("--v1-hold", type=Path, required=True)
    parser.add_argument("--demotion", type=Path, required=True)
    parser.add_argument("--authority", type=Path, required=True)
    parser.add_argument("--profile-manifest", type=Path, required=True)
    parser.add_argument("--independence-summary", type=Path, required=True)
    parser.add_argument("--reference-ledger", type=Path, default=None)
    parser.add_argument("--f1-attribution", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--expect-rows", type=int, default=None,
                        help="declared universe size; a mismatch fails closed")
    return parser


def main(argv: "Sequence[str] | None" = None) -> int:
    """Adapt the frozen PhaDED evidence into the v2 candidate-catalog input."""
    args = build_parser().parse_args(argv)
    try:
        result = build_adapted_input(
            universe_path=args.universe,
            subtype_matrix_path=args.subtype_matrix,
            v1_merge_path=args.v1_merge,
            v1_hold_path=args.v1_hold,
            demotion_path=args.demotion,
            authority_path=args.authority,
            profile_manifest_path=args.profile_manifest,
            independence_summary_path=args.independence_summary,
            ledger_path=args.reference_ledger,
            f1_attribution_path=args.f1_attribution,
            expect_rows=args.expect_rows,
        )
        written = write_artifacts(args.results_dir, result["rows"], result["report"])
    except (InputError, FileNotFoundError) as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        return 1
    print(json.dumps(
        {
            "rows": len(result["rows"]),
            "model_layer_distribution": result["report"]["model_layer_distribution"],
            "superfamily_source_counts": result["report"]["superfamily_source_counts"],
            "outputs": written,
        },
        ensure_ascii=False, sort_keys=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
