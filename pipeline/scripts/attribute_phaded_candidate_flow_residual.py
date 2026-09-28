#!/usr/bin/env python3
"""Attribute the 616-accession residual of the frozen v1 PhaDED candidate flow.

The frozen v1 ledger does not close::

    109,087 (v1 candidate universe)  -  36,611 (pool-in high confidence)
                                     -  71,860 (pool-in hold)          =  616

``reconcile_phaded_candidate_flow.check_hold_residual`` reports that residual as
``status="mismatch"`` with ``residual_unexplained_count=616``.  That loud failure
is **intentional** and this module does **not** silence it: the check is called
with an empty ``other_dispositions`` mapping and its payload is written verbatim
to ``check_hold_residual.json``, so the alarm survives beside the attribution.
Attribution is **not** a disposition: filing the 616 into buckets does not make
them disposed, and only F14 (v2 candidate catalog) may do that.

What this module proves, per accession
--------------------------------------
Step 1 derives the residual by an explicit accession **set difference**
(``universe - high_confidence - hold``) and never by subtracting the three
published aggregates; the aggregate expression is recomputed only to *show*
whether it agrees, and a disagreement is reported, not hidden.

Step 2 attributes every residual accession to the stage that provably dropped it.
The stage is the superfamily-assignment step of the frozen v1 filter
(``filter_phaded_high_confidence.py``, lines 207-228, ``if sf is None: continue``):
a candidate that gets no admissible superfamily call is written to *neither* the
high-confidence table nor the hold table, and therefore carries no ``hold_reason``
at all.  ``check_assignment_completeness`` proves that reading in both directions:
the frozen assignment evidence yields exactly ``|high_confidence ∪ hold|``
admissible accessions, and every admissible accession has a flow row.

Run layout
----------
``--output-dir`` receives every file this module writes; nothing is ever written
outside it, and an existing non-empty output directory is refused (an evidence
directory must never be written into).  A missing input path is always an error --
an absent input is never treated as an empty input.

Attribution vocabulary
----------------------
The vocabulary requested for F1 is kept verbatim and *extended only where the
frozen data shows a mechanism the requested list does not name*.  ``ATTRIBUTION_BASIS``
records for each value whether it was requested or added as an extension, and
``ATTRIBUTION_STATUS`` records whether the value was observed, measured zero, or
is structurally unreachable for a residual accession.  An accession whose frozen
evidence does not prove a stage is reported as ``unexplained`` with its own count:
it is never pushed into the nearest bucket.

Boundaries
----------
Read-only analysis with the Python standard library.  This module reads frozen
evidence, writes only into ``--output-dir``, never touches an existing run, and
never modifies a candidate, a table or a count.  It performs no scoring, no
recalculation and no registry change; the attribution is a statement about the
frozen v1 tables, not a new piece of biological evidence.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

SCHEMA_VERSION = "1.0"

_SCRIPTS_DIR = Path(__file__).resolve().parent


def _load_sibling(module_name: str):
    """Import a sibling script of ``pipeline/scripts`` by path.

    The scripts directory is not a Python package, so a plain ``import`` only
    works when the caller already put it on ``sys.path`` (tests do).  Loading by
    path keeps this module usable both as a script and as a test target.
    """
    if str(_SCRIPTS_DIR) not in sys.path:
        sys.path.insert(0, str(_SCRIPTS_DIR))
    try:
        return importlib.import_module(module_name)
    except ImportError:  # pragma: no cover - exercised only in odd sys.path states
        spec = importlib.util.spec_from_file_location(
            module_name, _SCRIPTS_DIR / ("%s.py" % module_name)
        )
        if spec is None or spec.loader is None:
            raise ImportError("cannot load sibling module %s" % module_name) from None
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return module


#: The frozen v1 ledger/checks module.  Reused instead of re-implemented: the
#: hold-residual check, the accession readers and the hashers are *the same
#: objects* here, which the test suite asserts by identity.
reconcile = _load_sibling("reconcile_phaded_candidate_flow")

#: The frozen v1 high-confidence filter.  Imported only for its superfamily
#: registry range (``SUPERFAMILY_CRITERIA``), so "outside the registry range" is
#: decided by the frozen authority rather than by a local copy of the list.
frozen_filter = _load_sibling("filter_phaded_high_confidence")

check_hold_residual = reconcile.check_hold_residual
read_accession_list = reconcile.read_accession_list
read_reason_table = reconcile.read_reason_table
sha256_file = reconcile.sha256_file
ensure_output_dir = reconcile.ensure_output_dir
write_tsv = reconcile.write_tsv
write_json = reconcile.write_json
bound_file = reconcile.bound_file
STATUS_OK = reconcile.STATUS_OK
STATUS_MISMATCH = reconcile.STATUS_MISMATCH
DISPOSITIONS = reconcile.DISPOSITIONS

#: ``check_hold_residual`` must keep reporting the residual as unexplained until
#: F14 assigns a real ``primary_disposition``: this module deliberately supplies
#: no ``other_dispositions`` bucket.
CHECK_HOLD_RESIDUAL_OTHER_DISPOSITIONS_SOURCE = (
    "not_supplied: F1 attribution is a stage analysis, not a primary_disposition"
)

PROFILE_EVIDENCE_ADMITTED = ("profile_trained_hit", "profile_ambiguous_family")
PROFILE_EVIDENCE_AMBIGUOUS = "profile_ambiguous_superfamily"
ASSIGNMENT_UNASSIGNED = "unassigned_PhaDED_like"
ASSIGNMENT_AMBIGUOUS = "ambiguous_superfamily"
EVIDENCE_STATUS_NO_PROFILE_SCORE = "no_profile_score"
EVIDENCE_STATUS_PROFILE_SCORE = "profile_score"
SCORE_GAP_THRESHOLD_LABEL = "min_score_gap_1.0_bits_HMMER_E1e6"
CONFIDENCE_UNIQUE = "unique"
AMBIGUOUS_CLAIM_PREFIX = "ambiguous:"

#: The frozen score-tier boundaries of the pool-external profile layer
#: (``build_phaded_pool_external_profile_layer.py::score_tier``).  Reused by
#: reference for the ``score_tier`` column so the column cannot drift from the
#: frozen definition.
SCORE_TIER_STRONG_E = 1e-30
SCORE_TIER_MID_E = 1e-10

# ---------------------------------------------------------------------------
# Attribution vocabulary
# ---------------------------------------------------------------------------

ATTRIBUTION_TIER_BELOW_THRESHOLD = "tier_below_threshold"
ATTRIBUTION_OUTSIDE_REGISTRY_FAMILY_RANGE = "outside_registry_family_range"
ATTRIBUTION_POOL_EXTERNAL_BOUNDARY = "pool_external_boundary"
ATTRIBUTION_SCORE_TIER_NOT_UPGRADED = "score_tier_not_upgraded"
ATTRIBUTION_PREFILTER_EXCLUDED = "prefilter_excluded"
ATTRIBUTION_TABLE_CONSTRUCTION_GAP = "table_construction_gap"
ATTRIBUTION_PROFILE_UNASSIGNED_NO_DISCOVERY_CLAIM = (
    "profile_unassigned_no_discovery_claim"
)
ATTRIBUTION_DISCOVERY_CLAIM_NOT_UNIQUE = "discovery_claim_not_unique"
ATTRIBUTION_UNEXPLAINED = "unexplained"

#: Every attribution value, in report order.  The counts of these values sum to
#: the residual size, so no value may ever be dropped from a report.
ATTRIBUTION_ORDER: tuple[str, ...] = (
    ATTRIBUTION_PROFILE_UNASSIGNED_NO_DISCOVERY_CLAIM,
    ATTRIBUTION_DISCOVERY_CLAIM_NOT_UNIQUE,
    ATTRIBUTION_SCORE_TIER_NOT_UPGRADED,
    ATTRIBUTION_POOL_EXTERNAL_BOUNDARY,
    ATTRIBUTION_OUTSIDE_REGISTRY_FAMILY_RANGE,
    ATTRIBUTION_TIER_BELOW_THRESHOLD,
    ATTRIBUTION_PREFILTER_EXCLUDED,
    ATTRIBUTION_TABLE_CONSTRUCTION_GAP,
    ATTRIBUTION_UNEXPLAINED,
)

#: Vocabulary values added here because the frozen data shows a mechanism the
#: requested list does not name.  Each one is justified in ``ATTRIBUTION_BASIS``.
ATTRIBUTION_EXTENSIONS: tuple[str, ...] = (
    ATTRIBUTION_PROFILE_UNASSIGNED_NO_DISCOVERY_CLAIM,
    ATTRIBUTION_DISCOVERY_CLAIM_NOT_UNIQUE,
)

#: ``observed`` -- fired on the frozen data; ``measured_zero`` -- the mechanism
#: exists and was measured to be empty; ``structural_zero`` -- the mechanism
#: cannot fire for a residual accession under the frozen filter, and the reason
#: is recorded; ``pending`` -- the evidence needed was not supplied.
ATTRIBUTION_STATUS: dict[str, str] = {
    ATTRIBUTION_PROFILE_UNASSIGNED_NO_DISCOVERY_CLAIM: "observed",
    ATTRIBUTION_DISCOVERY_CLAIM_NOT_UNIQUE: "observed",
    ATTRIBUTION_SCORE_TIER_NOT_UPGRADED: "observed",
    ATTRIBUTION_POOL_EXTERNAL_BOUNDARY: "measured_zero",
    ATTRIBUTION_OUTSIDE_REGISTRY_FAMILY_RANGE: "measured_zero",
    ATTRIBUTION_TIER_BELOW_THRESHOLD: "structural_zero",
    ATTRIBUTION_PREFILTER_EXCLUDED: "structural_zero",
    ATTRIBUTION_TABLE_CONSTRUCTION_GAP: "measured_zero",
    ATTRIBUTION_UNEXPLAINED: "measured_zero",
}

ATTRIBUTION_BASIS: dict[str, str] = {
    ATTRIBUTION_PROFILE_UNASSIGNED_NO_DISCOVERY_CLAIM: (
        "EXTENSION (not in the requested list): no superfamily call from either "
        "frozen source is admissible.  The trained-profile layer scored the "
        "accession and produced no call (subtype_matrix::assignment_status=unassigned_PhaDED_like, "
        "subtype_matrix::evidence_status=no_profile_score, "
        "subtype_matrix::phaded_superfamily_best is empty, "
        "subtype_matrix::superfamily_score_gap below the frozen "
        "subtype_matrix::evidence_threshold %s), and the discovery layer emitted "
        "no row at all in superfamily_classification.tsv, whose own "
        "claimed_candidates sum equals its row count." % SCORE_GAP_THRESHOLD_LABEL
    ),
    ATTRIBUTION_DISCOVERY_CLAIM_NOT_UNIQUE: (
        "EXTENSION (not in the requested list): the discovery layer did emit a row "
        "but its superfamily_claim is not unique, so the frozen filter's "
        "'superfamily_confidence == \"unique\"' gate refuses it "
        "(superfamily_classification.tsv::superfamily_confidence, "
        "::n_superfamilies_hit).  The trained-profile layer produced no call for "
        "the same accession, so there is no fallback source."
    ),
    ATTRIBUTION_SCORE_TIER_NOT_UPGRADED: (
        "REQUESTED value, observed meaning refined: a trained-profile score exists "
        "(subtype_matrix::evidence_status=profile_score) but it was not upgraded to "
        "a superfamily call (subtype_matrix::profile_evidence_status="
        "profile_ambiguous_superfamily, which the frozen filter's admitted set "
        "{profile_trained_hit, profile_ambiguous_family} excludes), and "
        "subtype_matrix::superfamily_score_gap is below the frozen "
        "subtype_matrix::evidence_threshold %s." % SCORE_GAP_THRESHOLD_LABEL
    ),
    ATTRIBUTION_POOL_EXTERNAL_BOUNDARY: (
        "REQUESTED value: the accession is present in a supplied pool-external "
        "layer (its own accession column), whose layer rule -- score-only, NOT high "
        "confidence -- governs it, so the pool-in tables cannot explain it."
    ),
    ATTRIBUTION_OUTSIDE_REGISTRY_FAMILY_RANGE: (
        "REQUESTED value: a superfamily call exists but names a superfamily outside "
        "the frozen v1 filter registry "
        "(filter_phaded_high_confidence.py::SUPERFAMILY_CRITERIA, cross-checked "
        "against filter_summary.json::per_superfamily)."
    ),
    ATTRIBUTION_TIER_BELOW_THRESHOLD: (
        "REQUESTED value, structurally unreachable for a residual accession: on this "
        "path the frozen gate is a score-gap threshold "
        "(subtype_matrix::evidence_threshold %s), not an E-value tier, and the "
        "frozen discovery branch tests uniqueness BEFORE the score tier, so a "
        "non-unique claim never reaches a tier test.  Every residual accession "
        "carries subtype_matrix::superfamily_score_gap below that threshold; the "
        "numeric fact is recorded in the score_gap_bits column instead of being "
        "reported as a separate bucket." % SCORE_GAP_THRESHOLD_LABEL
    ),
    ATTRIBUTION_PREFILTER_EXCLUDED: (
        "REQUESTED value, structurally unreachable for a residual accession: the "
        "confounder test sits in the frozen filter's *common criteria* block, which "
        "runs after the superfamily-assignment step, so a recorded confounder "
        "without an admissible assignment is dropped before that test and lands in "
        "its stage bucket instead.  Recorded confounders are carried as the "
        "confounder_flag column, never as the primary attribution -- calling them "
        "'prefilter_excluded' would name a stage the frozen code never reached."
    ),
    ATTRIBUTION_TABLE_CONSTRUCTION_GAP: (
        "REQUESTED value: the frozen evidence proves an admissible assignment exists "
        "(subtype_matrix::phaded_superfamily_best inside the registry with "
        "subtype_matrix::profile_evidence_status in {profile_trained_hit, "
        "profile_ambiguous_family}, or a unique in-registry "
        "superfamily_classification.tsv::superfamily_claim) yet the accession is in "
        "neither high_confidence_candidates.tsv nor hold_candidates.tsv.  That is a "
        "defect of the flow tables, not a missing-evidence finding."
    ),
    ATTRIBUTION_UNEXPLAINED: (
        "REQUESTED value: the frozen evidence does not prove a stage for this "
        "accession -- no row in any frozen evidence table, or column values that "
        "match none of the documented mechanisms.  Reported with its own count and "
        "accession list; never folded into a neighbouring bucket."
    ),
}

#: Attribution -> the disposition RULE proposed for that attribution (not a
#: per-accession hand edit).  ``primary_disposition`` is one of the six v2
#: dispositions of the redesign spec.
DISPOSITION_RULES: dict[str, dict[str, str]] = {
    ATTRIBUTION_PROFILE_UNASSIGNED_NO_DISCOVERY_CLAIM: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "subtype_matrix::assignment_status=unassigned_PhaDED_like; "
            "subtype_matrix::evidence_status=no_profile_score; "
            "subtype_matrix::phaded_superfamily_best=''; "
            "subtype_matrix::superfamily_score_gap < subtype_matrix::evidence_threshold; "
            "superfamily_classification.tsv has no row for the accession "
            "(candidate_union.faa carries all 109,087 accessions, so the row is "
            "absent because no family was hit, not because the table was built short)"
        ),
        "justification": (
            "Two independent frozen sources agree that no family/superfamily call is "
            "admissible, so the accession cannot be placed in a homology layer.  It is "
            "also not excludable: AGENTS.md reserves excluded_input_quality for "
            "documented input defects and forbids using it for scientific uncertainty. "
            "The candidates stay in the catalog as unresolved-function candidates with "
            "their existing evidence flags."
        ),
        "action": (
            "Retain as candidates; attach evidence_flags from the frozen columns "
            "(interpro generic domain support, lipase-box/PF06850 motif state) and "
            "leave them recallable by a future calibrated run.  Do not delete, do not "
            "demote."
        ),
    },
    ATTRIBUTION_DISCOVERY_CLAIM_NOT_UNIQUE: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "superfamily_classification.tsv::superfamily_confidence=ambiguous_2; "
            "::n_superfamilies_hit=2; ::superfamily_claim=ambiguous:<sf_a>|<sf_b> "
            "(every named superfamily is inside the v1 registry); "
            "::model_layer=discovery_hmm_uncalibrated; "
            "subtype_matrix::phaded_superfamily_best=''"
        ),
        "justification": (
            "The only family evidence is an uncalibrated discovery-layer claim that "
            "fell into two superfamilies; AGENTS.md forbids using the discovery layer "
            "for screening or demotion, so its ambiguity may neither promote the "
            "accession into a homology layer nor exclude it.  Undetermined family "
            "identity is unresolved function, not a negative result."
        ),
        "action": (
            "Retain as candidates with evidence_flag=superfamily_call_ambiguous_2way "
            "and the two competing superfamilies recorded; resolve by a calibrated "
            "run or by the v2 evidence layer, never by the discovery layer."
        ),
    },
    ATTRIBUTION_SCORE_TIER_NOT_UPGRADED: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "subtype_matrix::evidence_status=profile_score; "
            "subtype_matrix::profile_evidence_status=profile_ambiguous_superfamily; "
            "subtype_matrix::assignment_status=ambiguous_superfamily; "
            "subtype_matrix::superfamily_score_gap < subtype_matrix::evidence_threshold "
            "(min_score_gap_1.0_bits_HMMER_E1e6)"
        ),
        "justification": (
            "A trained-profile score exists but the winning superfamily is ambiguous "
            "below the preregistered 1.0-bit separation, so the published layer must "
            "not carry a family call.  The evidence is real but not decisive: this is "
            "unresolved function with a recorded near-miss, not an exclusion."
        ),
        "action": (
            "Retain as candidates with evidence_flag=profile_superfamily_ambiguous and "
            "the measured score_gap_bits; they are the first re-scoring targets if the "
            "separation rule is ever re-registered (separate authorisation)."
        ),
    },
    ATTRIBUTION_POOL_EXTERNAL_BOUNDARY: {
        "primary_disposition": "remote_homolog_candidate",
        "proving_evidence": (
            "<pool-external layer file>::accession (the layer's own column); "
            "<pool-external layer file>::score_tier and ::evidence_layer=score_only"
        ),
        "justification": (
            "PROVISIONAL rule: an accession that belongs to a pool-external layer is "
            "governed by that layer's published rule (score-only, never called high "
            "confidence).  Within that layer a strong tier denotes candidate remote "
            "homology; mid/weak tiers stay unresolved.  This bucket is empty on the "
            "frozen pool-in residual (measured intersection 0), so the rule is stated "
            "for completeness and must be re-checked when a layer is supplied."
        ),
        "action": (
            "Dispose under the pool-external rule: remote_homolog_candidate only for "
            "score_tier=strong, otherwise function_unresolved; never excluded."
        ),
    },
    ATTRIBUTION_OUTSIDE_REGISTRY_FAMILY_RANGE: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "superfamily_classification.tsv::superfamily_claim (or "
            "subtype_matrix::phaded_superfamily_best) names a superfamily absent from "
            "filter_phaded_high_confidence.py::SUPERFAMILY_CRITERIA; cross-checked "
            "against filter_summary.json::per_superfamily"
        ),
        "justification": (
            "A call outside the frozen registry means the registry range, not the "
            "accession, is the open question.  Silently dropping such an accession is "
            "exactly the v1 defect this task exists to expose, and exclusions are "
            "forbidden for a scope question."
        ),
        "action": (
            "Retain as candidates with evidence_flag=superfamily_outside_v1_registry "
            "and raise a registry-scope decision; do not extend the registry as a "
            "side effect of an attribution."
        ),
    },
    ATTRIBUTION_TIER_BELOW_THRESHOLD: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "subtype_matrix::evidence_threshold=min_score_gap_1.0_bits_HMMER_E1e6 and "
            "subtype_matrix::superfamily_score_gap (below threshold); "
            "superfamily_classification.tsv::best_evalue tier under "
            "build_phaded_pool_external_profile_layer.py::score_tier"
        ),
        "justification": (
            "PROVISIONAL rule for a structurally unreachable bucket: an accession "
            "whose only evidence is a below-threshold score has no admissible call, so "
            "it is unresolved function.  Only a strong tier under the frozen score_tier "
            "definition could ever justify a remote-homology layer, and it does not "
            "occur on this residual."
        ),
        "action": (
            "function_unresolved with evidence_flag=score_below_admission_threshold; "
            "if this bucket ever becomes non-empty, the flow tables disagree with the "
            "assignment evidence and must be re-derived before any disposition."
        ),
    },
    ATTRIBUTION_PREFILTER_EXCLUDED: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "confounder_candidates.tsv::confounder_flag and "
            "::recommended_confidence_cap (recorded per accession, carried as the "
            "confounder_flag column)"
        ),
        "justification": (
            "PROVISIONAL rule for a structurally unreachable bucket: a confounder flag "
            "is a scientific confounder suspicion, not a documented input-quality "
            "defect, so it can never justify excluded_input_quality.  The frozen table "
            "itself only recommends a confidence cap."
        ),
        "action": (
            "Retain as function_unresolved with "
            "evidence_flag=lipase_esterase_confounder_suspected and apply "
            "recommended_confidence_cap as a cap, never as a deletion."
        ),
    },
    ATTRIBUTION_TABLE_CONSTRUCTION_GAP: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "subtype_matrix::phaded_superfamily_best + "
            "subtype_matrix::profile_evidence_status, or "
            "superfamily_classification.tsv::superfamily_claim + "
            "::superfamily_confidence=unique, against the "
            "high_confidence_candidates.tsv / hold_candidates.tsv row sets"
        ),
        "justification": (
            "PROVISIONAL rule for a measured-zero bucket: the assignment evidence says "
            "the accession was admissible while the flow tables do not carry it.  That "
            "is an engineering inconsistency in the frozen tables, so no scientific "
            "disposition may be inferred from it."
        ),
        "action": (
            "function_unresolved and block ledger closure: escalate as a table defect "
            "and re-derive the flow tables before F14 relies on them."
        ),
    },
    ATTRIBUTION_UNEXPLAINED: {
        "primary_disposition": "function_unresolved",
        "proving_evidence": (
            "subtype_matrix::<no row for the accession>; "
            "superfamily_classification::<no row for the accession>; "
            "motif_completion_full::<no row for the accession>; or frozen column "
            "values outside every documented mechanism "
            "(subtype_matrix::assignment_status / ::evidence_status / "
            "superfamily_classification::superfamily_confidence)"
        ),
        "justification": (
            "An unexplained residual record must never be force-fitted: the honest "
            "output is unresolved function plus an explicit pending marker, because "
            "inventing a stage would fabricate evidence."
        ),
        "action": (
            "Retain as function_unresolved with evidence_flag=attribution_pending and "
            "list for human review before F14; excluded_input_quality is explicitly "
            "not available for an unexplained record."
        ),
    },
}

#: Logical evidence-source names used in every citation, so a claim can be
#: resolved to a bound path + SHA-256 in the run summary.
EVIDENCE_SOURCES: dict[str, str] = {
    "protein_layers": "v1 candidate universe (protein_layers.tsv)",
    "high_confidence_candidates": "pool-in high-confidence candidates",
    "hold_candidates": "pool-in hold candidates",
    "filter_summary": "pool-in filter summary (filter_summary.json)",
    "subtype_matrix": "full-library subtype matrix",
    "superfamily_classification": "family discovery superfamily classification",
    "motif_completion_full": "full-library motif completion table",
    "confounder_candidates": "x1 hydrophobicity confounder candidates",
}

#: Evidence columns read out of each frozen table, with the normalized names used
#: internally.  A missing column is an input error (never silently empty).
MATRIX_COLUMNS: dict[str, str] = {
    "accession": "accession",
    "assignment_status": "matrix_assignment_status",
    "evidence_status": "matrix_evidence_status",
    "evidence_threshold": "matrix_evidence_threshold",
    "phaded_superfamily_best": "matrix_phaded_superfamily_best",
    "superfamily_score_gap": "matrix_superfamily_score_gap",
    "profile_evidence_status": "matrix_profile_evidence_status",
    "sequence_integrity": "matrix_sequence_integrity",
    "evidence_grade": "matrix_evidence_grade",
}

CLASSIFICATION_COLUMNS: dict[str, str] = {
    "protein_id": "accession",
    "superfamily_claim": "discovery_superfamily_claim",
    "superfamily_confidence": "discovery_superfamily_confidence",
    "n_superfamilies_hit": "discovery_n_superfamilies_hit",
    "best_evalue": "discovery_best_evalue",
    "model_layer": "discovery_model_layer",
}

MOTIF_COLUMNS: dict[str, str] = {
    "accession": "accession",
    "candidate_prior_superfamily": "motif_candidate_prior_superfamily",
    "candidate_assignment_status": "motif_candidate_assignment_status",
    "prior_anchor_basis": "motif_prior_anchor_basis",
    "lipase_box_state": "motif_lipase_box_state",
    "sbd_pf06850_binding_state": "motif_sbd_pf06850_binding_state",
}

CONFOUNDER_COLUMNS: dict[str, str] = {
    "accession": "accession",
    "confounder_flag": "confounder_flag",
    "recommended_confidence_cap": "confounder_cap",
}

#: Normalized evidence-row shape.  Every key is always present; an empty string
#: means "the frozen source carries no value here", never "not looked at".
EVIDENCE_DEFAULTS: dict[str, str] = {
    "v1_layer": "",
    "v1_family": "",
    "matrix_present": "false",
    "matrix_assignment_status": "",
    "matrix_evidence_status": "",
    "matrix_evidence_threshold": "",
    "matrix_phaded_superfamily_best": "",
    "matrix_superfamily_score_gap": "",
    "matrix_profile_evidence_status": "",
    "matrix_sequence_integrity": "",
    "matrix_evidence_grade": "",
    "discovery_present": "false",
    "discovery_superfamily_claim": "",
    "discovery_superfamily_confidence": "",
    "discovery_n_superfamilies_hit": "",
    "discovery_best_evalue": "",
    "discovery_model_layer": "",
    "motif_present": "false",
    "motif_candidate_prior_superfamily": "",
    "motif_candidate_assignment_status": "",
    "motif_prior_anchor_basis": "",
    "motif_lipase_box_state": "",
    "motif_sbd_pf06850_binding_state": "",
    "confounder_flag": "",
    "confounder_cap": "",
}

RESIDUAL_FIELDS = [
    "accession",
    "layer",
    "family",
    "major_superfamily",
    "major_superfamily_layer",
    "major_superfamily_source",
    "discovery_claim_recall_only",
    "score_tier",
    "best_evalue",
    "score_gap_bits",
    "matrix_assignment_status",
    "matrix_evidence_status",
    "matrix_evidence_threshold",
    "matrix_profile_evidence_status",
    "discovery_row_present",
    "discovery_superfamily_confidence",
    "discovery_n_superfamilies_hit",
    "motif_candidate_assignment_status",
    "confounder_flag",
    "sequence_integrity",
    "attribution",
    "attribution_basis",
    "attribution_evidence",
    "pool_external_layers",
]

ATTRIBUTION_FIELDS = [
    "attribution",
    "count",
    "primary_disposition",
    "status",
    "proving_evidence",
    "accessions",
]

MANIFEST_FIELDS = ["name", "path", "size", "sha256", "status"]


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def load_frozen_filter_registry() -> tuple[str, ...]:
    """Return the superfamily registry range of the frozen v1 filter, sorted."""
    registry = tuple(sorted(getattr(frozen_filter, "SUPERFAMILY_CRITERIA")))
    if not registry:
        raise ValueError("the frozen v1 filter exposes no SUPERFAMILY_CRITERIA registry")
    return registry


frozen_filter_registry_for_testing = load_frozen_filter_registry


def score_tier(evalue: Any) -> str:
    """Bucket an E-value with the frozen pool-external boundaries."""
    try:
        value = float(str(evalue))
    except (TypeError, ValueError):
        return "no_score"
    if value < SCORE_TIER_STRONG_E:
        return "strong"
    if value < SCORE_TIER_MID_E:
        return "mid"
    return "weak"


def claim_members(claim: str) -> list[str]:
    """Split a (possibly ambiguous) superfamily claim into its named members."""
    text = (claim or "").strip()
    if text.startswith(AMBIGUOUS_CLAIM_PREFIX):
        text = text[len(AMBIGUOUS_CLAIM_PREFIX):]
        return [part.strip() for part in text.split("|") if part.strip()]
    return [text] if text else []


def make_evidence_row(accession: str, **overrides: Any) -> dict[str, str]:
    """Build one normalized evidence row; every field defaults to "absent"."""
    row = dict(EVIDENCE_DEFAULTS)
    row["accession"] = str(accession)
    unknown = sorted(set(overrides) - set(EVIDENCE_DEFAULTS) - {"accession"})
    if unknown:
        raise ValueError("unknown evidence field(s): %s" % ", ".join(unknown))
    for key, value in overrides.items():
        if isinstance(value, bool):
            row[key] = "true" if value else "false"
        else:
            row[key] = "" if value is None else str(value)
    return row


def _flag(row: Mapping[str, Any], key: str) -> bool:
    return str(row.get(key, "")).strip().lower() in {"1", "true", "yes"}


def _value(row: Mapping[str, Any], key: str) -> str:
    return str(row.get(key, "") or "").strip()


def check_registry_range(
    observed_names: Iterable[str], registry_names: Iterable[str]
) -> dict[str, Any]:
    """Compare observed superfamily names with the frozen registry range."""
    observed = {str(name).strip() for name in observed_names if str(name).strip()}
    registry = {str(name).strip() for name in registry_names if str(name).strip()}
    if not registry:
        raise ValueError("the frozen registry range is empty")
    outside = sorted(observed - registry)
    return {
        "check": "check_registry_range",
        "status": STATUS_OK if not outside else STATUS_MISMATCH,
        "explanation": (
            "%d observed superfamily name(s) are all inside the frozen registry range"
            % len(observed)
            if not outside
            else "%d observed superfamily name(s) are outside the frozen registry "
            "range: %s" % (len(outside), ", ".join(outside))
        ),
        "registry_count": len(registry),
        "observed_count": len(observed),
        "outside_the_frozen_registry": outside,
    }


def check_assignment_completeness(
    *,
    admissible: Iterable[str],
    high_confidence: Iterable[str],
    hold: Iterable[str],
) -> dict[str, Any]:
    """Prove ``high_confidence ∪ hold`` equals the admissible-assignment set.

    An accession with an admissible superfamily call always reaches one of the two
    frozen flow tables, and an accession carrying a flow row always has an
    admissible call.  Either direction failing is a mismatch: it means the flow
    tables and the assignment evidence disagree, so an attribution built on them
    would be unsound.
    """
    admissible_set = {str(value) for value in admissible if str(value)}
    flow_union = {str(value) for value in high_confidence if str(value)} | {
        str(value) for value in hold if str(value)
    }
    if not admissible_set and not flow_union:
        raise ValueError("assignment completeness check needs non-empty inputs")
    admissible_without_a_flow_row = sorted(admissible_set - flow_union)
    flow_row_without_an_admissible_assignment = sorted(flow_union - admissible_set)
    status = (
        STATUS_OK
        if not admissible_without_a_flow_row and not flow_row_without_an_admissible_assignment
        else STATUS_MISMATCH
    )
    problems: list[str] = []
    if admissible_without_a_flow_row:
        problems.append(
            "%d admissible accession(s) have no flow row: %s"
            % (
                len(admissible_without_a_flow_row),
                ", ".join(admissible_without_a_flow_row[:5]),
            )
        )
    if flow_row_without_an_admissible_assignment:
        problems.append(
            "%d flow row(s) have no admissible assignment: %s"
            % (
                len(flow_row_without_an_admissible_assignment),
                ", ".join(flow_row_without_an_admissible_assignment[:5]),
            )
        )
    return {
        "check": "check_assignment_completeness",
        "status": status,
        "explanation": (
            "the frozen assignment evidence yields %d admissible accession(s) and the "
            "two flow tables carry %d accession(s); the sets are identical, so the "
            "only way to fall out of both tables is to have no admissible superfamily "
            "call" % (len(admissible_set), len(flow_union))
            if status == STATUS_OK
            else "the assignment evidence and the flow tables disagree: "
            + "; ".join(problems)
        ),
        "admissible_count": len(admissible_set),
        "flow_union_count": len(flow_union),
        "admissible_without_a_flow_row": admissible_without_a_flow_row[:50],
        "admissible_without_a_flow_row_count": len(admissible_without_a_flow_row),
        "flow_row_without_an_admissible_assignment": (
            flow_row_without_an_admissible_assignment[:50]
        ),
        "flow_row_without_an_admissible_assignment_count": len(
            flow_row_without_an_admissible_assignment
        ),
        "problems": problems,
    }


def exit_code_for(attribution_closes: bool, hold_residual_status: str, residual_count: int) -> int:
    """Return the CLI exit code.

    ``0`` the attribution closes (and, when a residual exists, the frozen alarm is
    still reported); ``1`` the attribution does not close -- a defect of this
    analysis; ``2`` the residual is non-empty while ``check_hold_residual``
    reported ``ok``, i.e. the intentional alarm was silenced and the run must not
    be trusted.
    """
    if not attribution_closes:
        return 1
    if residual_count and hold_residual_status == STATUS_OK:
        return 2
    return 0


def disposition_rule(attribution: str) -> dict[str, str]:
    """Return a copy of the disposition rule of ``attribution``."""
    if attribution not in DISPOSITION_RULES:
        raise ValueError("no disposition rule for attribution %r" % (attribution,))
    rule = dict(DISPOSITION_RULES[attribution])
    rule["attribution"] = attribution
    rule["status"] = ATTRIBUTION_STATUS[attribution]
    rule["basis"] = ATTRIBUTION_BASIS[attribution]
    return rule


def attribution_basis_tag(attribution: str) -> str:
    """One-line origin tag: requested vocabulary, or an extension of it."""
    if attribution in ATTRIBUTION_EXTENSIONS:
        return "EXTENSION: not in the requested F1 list, justified from the frozen data"
    return "REQUESTED: kept from the F1 vocabulary"


# ---------------------------------------------------------------------------
# The attribution itself (pure)
# ---------------------------------------------------------------------------


def _matrix_branch_admissible(row: Mapping[str, Any], registry: set[str]) -> bool:
    best = _value(row, "matrix_phaded_superfamily_best")
    return bool(best) and best in registry and _value(
        row, "matrix_profile_evidence_status"
    ) in PROFILE_EVIDENCE_ADMITTED


def _discovery_branch_admissible(row: Mapping[str, Any], registry: set[str]) -> bool:
    if not _flag(row, "discovery_present") and not _value(row, "discovery_superfamily_claim"):
        return False
    if _value(row, "discovery_superfamily_confidence") != CONFIDENCE_UNIQUE:
        return False
    members = claim_members(_value(row, "discovery_superfamily_claim"))
    return bool(members) and all(member in registry for member in members)


def _attribute_one(
    row: Mapping[str, Any], registry: set[str], external_layers: Mapping[str, set[str]]
) -> tuple[str, list[str], list[str]]:
    """Return ``(attribution, evidence_citations, pool_external_layer_names)``."""
    accession = row["accession"]
    layers = sorted(
        name for name, members in external_layers.items() if accession in members
    )
    if layers:
        return (
            ATTRIBUTION_POOL_EXTERNAL_BOUNDARY,
            [
                "pool_external_layer::accession contains %s" % accession,
                "pool_external_layer::score_tier (layer %s)" % ", ".join(layers),
            ],
            layers,
        )

    best = _value(row, "matrix_phaded_superfamily_best")
    claim = _value(row, "discovery_superfamily_claim")
    confidence = _value(row, "discovery_superfamily_confidence")
    n_superfamilies = _value(row, "discovery_n_superfamilies_hit")
    evidence_status = _value(row, "matrix_evidence_status")
    profile_status = _value(row, "matrix_profile_evidence_status")
    assignment_status = _value(row, "matrix_assignment_status")
    gap = _value(row, "matrix_superfamily_score_gap")
    threshold = _value(row, "matrix_evidence_threshold")
    members = claim_members(claim)
    outside = [member for member in members if member not in registry]

    # 1. a call exists but names a superfamily outside the frozen registry range
    if (best and best not in registry) or outside:
        citations = []
        if best and best not in registry:
            citations.append(
                "subtype_matrix::phaded_superfamily_best=%s (outside the registry)" % best
            )
        if outside:
            citations.append(
                "superfamily_classification::superfamily_claim=%s (outside the registry)"
                % ";".join(outside)
            )
        return (
            ATTRIBUTION_OUTSIDE_REGISTRY_FAMILY_RANGE,
            citations + ["filter_phaded_high_confidence::SUPERFAMILY_CRITERIA"],
            [],
        )

    # 2. a trained-profile score exists but was not upgraded to a superfamily call
    if evidence_status == EVIDENCE_STATUS_PROFILE_SCORE and not _matrix_branch_admissible(
        row, registry
    ):
        return (
            ATTRIBUTION_SCORE_TIER_NOT_UPGRADED,
            [
                "subtype_matrix::evidence_status=%s" % evidence_status,
                "subtype_matrix::profile_evidence_status=%s" % (profile_status or "<empty>"),
                "subtype_matrix::assignment_status=%s" % (assignment_status or "<empty>"),
                "subtype_matrix::superfamily_score_gap=%s < subtype_matrix::evidence_threshold=%s"
                % (gap or "<empty>", threshold or "<empty>"),
            ],
            [],
        )

    # 3. the evidence proves an admissible assignment, yet there is no flow row
    if _matrix_branch_admissible(row, registry) or _discovery_branch_admissible(
        row, registry
    ):
        return (
            ATTRIBUTION_TABLE_CONSTRUCTION_GAP,
            [
                "subtype_matrix::phaded_superfamily_best=%s with "
                "subtype_matrix::profile_evidence_status=%s"
                % (best or "<empty>", profile_status or "<empty>"),
                "superfamily_classification::superfamily_claim=%s with "
                "::superfamily_confidence=%s" % (claim or "<empty>", confidence or "<empty>"),
            ],
            [],
        )

    # 4. a discovery row exists but its claim is not unique
    if confidence and confidence != CONFIDENCE_UNIQUE:
        return (
            ATTRIBUTION_DISCOVERY_CLAIM_NOT_UNIQUE,
            [
                "superfamily_classification::superfamily_confidence=%s" % confidence,
                "superfamily_classification::n_superfamilies_hit=%s"
                % (n_superfamilies or "<empty>"),
                "superfamily_classification::superfamily_claim=%s" % (claim or "<empty>"),
                "subtype_matrix::phaded_superfamily_best=<empty>",
            ],
            [],
        )

    # 5. no score at all and no admissible discovery claim
    if (
        evidence_status == EVIDENCE_STATUS_NO_PROFILE_SCORE
        and assignment_status == ASSIGNMENT_UNASSIGNED
        and not best
    ):
        return (
            ATTRIBUTION_PROFILE_UNASSIGNED_NO_DISCOVERY_CLAIM,
            [
                "subtype_matrix::assignment_status=%s" % assignment_status,
                "subtype_matrix::evidence_status=%s" % evidence_status,
                "subtype_matrix::phaded_superfamily_best=<empty>",
                "subtype_matrix::superfamily_score_gap=%s < subtype_matrix::evidence_threshold=%s"
                % (gap or "<empty>", threshold or "<empty>"),
                (
                    "superfamily_classification::<no row for the accession>"
                    if not claim
                    else "superfamily_classification::superfamily_claim=%s" % claim
                ),
            ],
            [],
        )

    # 6. the frozen columns prove no stage: say so, with the values that failed
    return (
        ATTRIBUTION_UNEXPLAINED,
        [
            "subtype_matrix::assignment_status=%s" % (assignment_status or "<empty>"),
            "subtype_matrix::evidence_status=%s" % (evidence_status or "<empty>"),
            "superfamily_classification::superfamily_confidence=%s"
            % (confidence or "<empty>"),
            "subfamily evidence matches no documented mechanism",
        ],
        [],
    )


def attribute_residual(
    residual_accessions: Sequence[str],
    layer_rows: Sequence[Mapping[str, Any]],
    filter_summary: Mapping[str, Any],
    external_layer: Mapping[str, Iterable[str]] | Iterable[str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Attribute every residual accession to the stage that provably dropped it.

    Returns ``(rows, summary)``.  ``rows`` carries one entry per residual
    accession (sorted) with exactly one ``attribution``; ``summary`` carries the
    per-attribution counts, the ``unexplained`` count and the accession lists.

    ``layer_rows`` are normalized evidence rows (see :data:`EVIDENCE_DEFAULTS`) --
    typically one per residual accession.  A residual accession without an
    evidence row is reported as ``unexplained``; an evidence row outside the
    residual is counted in ``evidence_rows_outside_residual`` rather than dropped
    silently.

    ``filter_summary`` is the frozen ``filter_summary.json`` payload.  It is used
    for the registry-range cross-check only, so an empty mapping is accepted for a
    caller with no summary at hand; ``filter_summary_used`` records which case
    applied.  The registry range itself always comes from the frozen filter
    (``filter_phaded_high_confidence::SUPERFAMILY_CRITERIA``), never from this
    mapping.

    ``external_layer`` maps a layer name to its accessions, or is a flat iterable
    of accessions (named ``pool_external``).  When nothing is supplied the
    ``pool_external_boundary`` bucket is reported as ``pending`` -- never as a
    proven zero.

    Raises ``ValueError`` for a duplicated residual accession or an evidence row
    without an accession.
    """
    residual_list = [str(value) for value in residual_accessions]
    counts = Counter(residual_list)
    duplicates = sorted(value for value, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(
            "residual accession(s) appear more than once: %s" % ", ".join(duplicates[:5])
        )
    residual_set = set(residual_list)

    external_layers: dict[str, set[str]] = {}
    external_supplied = False
    if external_layer is not None:
        if isinstance(external_layer, Mapping):
            external_layers = {
                str(name): {str(value) for value in values if str(value)}
                for name, values in external_layer.items()
            }
            external_supplied = bool(external_layers)
        else:
            flat = {str(value) for value in external_layer if str(value)}
            external_layers = {"pool_external": flat}
            external_supplied = True

    rows_by_accession: dict[str, Mapping[str, Any]] = {}
    total_evidence_rows = 0
    outside = 0
    for raw in layer_rows:
        if "accession" not in raw or not str(raw.get("accession", "")).strip():
            raise ValueError("evidence row without an accession: %r" % (dict(raw),))
        accession = str(raw["accession"]).strip()
        total_evidence_rows += 1
        if accession not in residual_set:
            outside += 1
            continue
        if accession in rows_by_accession:
            raise ValueError(
                "more than one evidence row for residual accession %s" % accession
            )
        rows_by_accession[accession] = raw

    if residual_list and not filter_summary:
        filter_summary_used = False
    else:
        filter_summary_used = bool(filter_summary)

    registry = set(load_frozen_filter_registry())
    summary_registry = set()
    if isinstance(filter_summary, Mapping):
        per_superfamily = filter_summary.get("per_superfamily") or {}
        if isinstance(per_superfamily, Mapping):
            summary_registry = {str(name) for name in per_superfamily}
    registry_check = (
        check_registry_range(summary_registry, registry) if summary_registry else None
    )

    attributed: list[dict[str, Any]] = []
    attribution_counts: Counter[str] = Counter()
    unexplained: list[str] = []
    confounder_flagged: list[str] = []
    for accession in sorted(residual_set):
        raw = rows_by_accession.get(accession)
        if raw is None:
            attribution = ATTRIBUTION_UNEXPLAINED
            basis = "no evidence row in any frozen evidence table supplied to this analysis"
            citations = [
                "no evidence row in any frozen evidence table supplied to this analysis",
                "attribution_basis: no evidence row",
            ]
            layers: list[str] = []
        else:
            attribution, citations, layers = _attribute_one(raw, registry, external_layers)
            basis = attribution_basis_tag(attribution)
        row = {
            "accession": accession,
            "layer": _value(raw or {}, "v1_layer"),
            "family": _value(raw or {}, "v1_family"),
            "attribution": attribution,
            "attribution_basis": basis,
            "attribution_evidence": "; ".join(citations),
            "pool_external_layers": ";".join(layers),
            "confounder_flag": _value(raw or {}, "confounder_flag"),
            "score_gap_bits": _value(raw or {}, "matrix_superfamily_score_gap"),
            "matrix_assignment_status": _value(raw or {}, "matrix_assignment_status"),
            "matrix_evidence_status": _value(raw or {}, "matrix_evidence_status"),
            "matrix_evidence_threshold": _value(raw or {}, "matrix_evidence_threshold"),
            "matrix_profile_evidence_status": _value(
                raw or {}, "matrix_profile_evidence_status"
            ),
            "sequence_integrity": _value(raw or {}, "matrix_sequence_integrity"),
            "discovery_row_present": (
                "true"
                if _flag(raw or {}, "discovery_present")
                or _value(raw or {}, "discovery_superfamily_claim")
                else "false"
            ),
            "discovery_superfamily_confidence": _value(
                raw or {}, "discovery_superfamily_confidence"
            ),
            "discovery_n_superfamilies_hit": _value(raw or {}, "discovery_n_superfamilies_hit"),
            "motif_candidate_assignment_status": _value(
                raw or {}, "motif_candidate_assignment_status"
            ),
            "best_evalue": _value(raw or {}, "discovery_best_evalue"),
        }
        row.update(major_superfamily_of(raw or {}))
        # Recall-only layer evidence, kept OUT of every family-call field: the
        # discovery layer must not produce a family call, so its raw claim is
        # recorded here under an explicitly recall-only name for audit only.
        row["discovery_claim_recall_only"] = _value(raw or {}, "discovery_superfamily_claim")
        row["score_tier"] = score_tier(row["best_evalue"]) if row["best_evalue"] else "no_score"
        attributed.append(row)
        attribution_counts[attribution] += 1
        if attribution == ATTRIBUTION_UNEXPLAINED:
            unexplained.append(accession)
        if row["confounder_flag"]:
            confounder_flagged.append(accession)

    attribution_counts_all = {
        name: attribution_counts.get(name, 0) for name in ATTRIBUTION_ORDER
    }
    pending_buckets: dict[str, str] = {}
    if not external_supplied:
        pending_buckets[ATTRIBUTION_POOL_EXTERNAL_BOUNDARY] = (
            "no pool-external layer was supplied, so the boundary could not be tested"
        )
    sum_counts = sum(attribution_counts_all.values())
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "residual_count": len(residual_set),
        "attribution_counts": attribution_counts_all,
        "sum_attribution_counts": sum_counts,
        "attribution_counts_sum_to_residual": sum_counts == len(residual_set),
        "unexplained_count": attribution_counts_all[ATTRIBUTION_UNEXPLAINED],
        "unexplained_accessions": unexplained,
        "pending_buckets": pending_buckets,
        "pending_counts": sorted(pending_buckets),
        "pool_external_layers_supplied": sorted(external_layers),
        "pool_external_layer_counts": {
            name: len(members) for name, members in sorted(external_layers.items())
        },
        "pool_external_layer_intersection_with_residual": {
            name: len(members & residual_set)
            for name, members in sorted(external_layers.items())
        },
        "attribution_status": {
            name: ATTRIBUTION_STATUS[name] for name in ATTRIBUTION_ORDER
        },
        "attribution_extensions": list(ATTRIBUTION_EXTENSIONS),
        "vocabulary": list(ATTRIBUTION_ORDER),
        "evidence_rows_total": total_evidence_rows,
        "evidence_rows_outside_residual": outside,
        "evidence_rows_used": len(rows_by_accession),
        "confounder_flagged_count": len(confounder_flagged),
        "confounder_flagged_accessions": confounder_flagged,
        "registry_range_check": registry_check,
        "filter_summary_used": filter_summary_used,
        "note": (
            "attribution names the stage that provably dropped an accession; it is "
            "NOT a primary_disposition -- see the disposition proposal for the rule "
            "that maps each attribution onto one"
        ),
    }
    return attributed, summary


def major_superfamily_of(row: Mapping[str, Any]) -> dict[str, str]:
    """Best available *admissible* superfamily evidence, with its source and layer.

    The v1 candidate universe table carries no superfamily column, so the value is
    read from the frozen assignment evidence.  Only sources that are themselves
    admissible family evidence are eligible: the trained-profile assignment
    (``subtype_matrix::phaded_superfamily_best``) and, failing that, the motif
    table's candidate prior (``motif_completion_full::candidate_prior_superfamily``).
    The source column records which frozen file+column produced the value, so the
    cell can never be mistaken for a value of ``protein_layers.tsv``.

    Discovery-layer rule (AGENTS.md): the discovery layer
    (``model_layer=discovery_hmm_uncalibrated``) is recall-only and **must not
    produce a family call** -- including indirectly, by having its value decide
    which superfamily a row reports.  The discovery layer's
    ``superfamily_claim`` is therefore *not* consulted here at all: it never
    selects, fills or gates ``major_superfamily``.  It is recorded separately by
    the row builder as an explicitly labelled recall-only field
    (``discovery_claim_recall_only``) so the evidence stays visible for audit
    without acting as a family call.  An empty ``major_superfamily`` therefore
    means "no admissible (non-discovery) frozen source carries a family call",
    which is exactly the condition that lets such an accession fall out of both
    the high-confidence and the hold tables.
    """
    best = _value(row, "matrix_phaded_superfamily_best")
    if best:
        return {
            "major_superfamily": best,
            "major_superfamily_layer": "sequence_family",
            "major_superfamily_source": "subtype_matrix::phaded_superfamily_best",
        }
    prior = _value(row, "motif_candidate_prior_superfamily")
    if prior:
        return {
            "major_superfamily": prior,
            "major_superfamily_layer": "candidate_prior",
            "major_superfamily_source": "motif_completion_full::candidate_prior_superfamily",
        }
    return {
        "major_superfamily": "",
        "major_superfamily_layer": "absent",
        "major_superfamily_source": (
            "absent: no admissible (non-discovery) frozen source carries a superfamily "
            "call for this accession"
        ),
    }


# ---------------------------------------------------------------------------
# Frozen-input readers
# ---------------------------------------------------------------------------


def _resolve_input(path: Path, what: str) -> Path:
    candidate = Path(path)
    if not candidate.exists():
        raise ValueError("%s path does not exist: %s" % (what, candidate))
    if not candidate.is_file():
        raise ValueError("%s path is not a regular file: %s" % (what, candidate))
    return candidate


def _open_table(path: Path, what: str, required: Sequence[str]) -> tuple[Any, list[str]]:
    handle = path.open(encoding="utf-8-sig", newline="")
    reader = csv.reader(handle, delimiter="\t")
    try:
        header = [column.strip() for column in next(reader)]
    except StopIteration:
        handle.close()
        raise ValueError("%s is empty: %s" % (what, path)) from None
    missing = [column for column in required if column not in header]
    if missing:
        handle.close()
        raise ValueError(
            "%s is missing columns: %s (found: %s)"
            % (what, ",".join(missing), ",".join(header))
        )
    return handle, header


def read_filter_summary(path: Path) -> dict[str, Any]:
    """Read the frozen ``filter_summary.json`` and require its documented keys."""
    source = _resolve_input(path, "filter summary")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError("filter summary is not valid JSON: %s" % error) from None
    if not isinstance(payload, Mapping):
        raise ValueError("filter summary must be a JSON object: %s" % source)
    missing = [
        key
        for key in ("total_high_confidence", "total_held", "hold_breakdown", "per_superfamily")
        if key not in payload
    ]
    if missing:
        raise ValueError(
            "filter summary is missing keys: %s" % ", ".join(missing)
        )
    return dict(payload)


def read_hold_accessions(path: Path) -> list[str]:
    """Return the sorted unique accessions of the frozen hold table."""
    return read_accession_list(path, "pool-in hold candidates")


def scan_subtype_matrix(
    path: Path,
    detail_accessions: Iterable[str],
    admissible_universe: Iterable[str] | None = None,
    registry: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Stream the subtype matrix once for per-accession detail and admissibility.

    Returns ``{"detail": {accession: normalized row}, "admissible": set,
    "rows_for_universe": set}``.  ``admissible`` implements the frozen filter's
    first superfamily branch: ``phaded_superfamily_best`` inside the registry and
    ``profile_evidence_status`` in the admitted set.
    """
    detail_wanted = {str(value) for value in detail_accessions}
    universe = {str(value) for value in (admissible_universe or ())}
    registry_set = {str(name) for name in (registry or load_frozen_filter_registry())}
    source = _resolve_input(path, "full-library subtype matrix")
    handle, header = _open_table(
        source, "full-library subtype matrix", list(MATRIX_COLUMNS)
    )
    detail: dict[str, dict[str, str]] = {}
    admissible: set[str] = set()
    rows_for_universe: set[str] = set()
    with handle:
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            accession = (raw.get("accession") or "").strip()
            if not accession:
                continue
            if universe:
                if accession not in universe:
                    continue
                rows_for_universe.add(accession)
            if accession in detail_wanted:
                detail[accession] = make_evidence_row(
                    accession,
                    matrix_present=True,
                    **{
                        normalized: (raw.get(column) or "")
                        for column, normalized in MATRIX_COLUMNS.items()
                        if normalized != "accession"
                    },
                )
            if (
                (raw.get("phaded_superfamily_best") or "").strip() in registry_set
                and (raw.get("profile_evidence_status") or "").strip()
                in PROFILE_EVIDENCE_ADMITTED
            ):
                admissible.add(accession)
    return {
        "detail": detail,
        "admissible": admissible,
        "rows_for_universe": rows_for_universe,
    }


def scan_superfamily_classification(
    path: Path,
    detail_accessions: Iterable[str],
    admissible_universe: Iterable[str] | None = None,
    registry: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Stream the discovery classification once for detail and admissibility.

    ``admissible`` implements the frozen filter's second superfamily branch: a
    ``unique`` claim whose members are all inside the registry range.
    """
    detail_wanted = {str(value) for value in detail_accessions}
    universe = {str(value) for value in (admissible_universe or ())}
    registry_set = {str(name) for name in (registry or load_frozen_filter_registry())}
    source = _resolve_input(path, "superfamily classification")
    handle, header = _open_table(
        source, "superfamily classification", list(CLASSIFICATION_COLUMNS)
    )
    detail: dict[str, dict[str, str]] = {}
    admissible: set[str] = set()
    rows_for_universe: set[str] = set()
    with handle:
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            accession = (raw.get("protein_id") or "").strip()
            if not accession:
                continue
            if universe:
                rows_for_universe.add(accession)
            if accession in detail_wanted:
                values = {
                    normalized: (raw.get(column) or "")
                    for column, normalized in CLASSIFICATION_COLUMNS.items()
                    if normalized != "accession"
                }
                detail[accession] = make_evidence_row(
                    accession, discovery_present=True, **values
                )
            claim = (raw.get("superfamily_claim") or "").strip()
            confidence = (raw.get("superfamily_confidence") or "").strip()
            members = claim_members(claim)
            if (
                confidence == CONFIDENCE_UNIQUE
                and members
                and all(member in registry_set for member in members)
            ):
                admissible.add(accession)
    return {
        "detail": detail,
        "admissible": admissible,
        "rows_for_universe": rows_for_universe,
    }


def scan_motif_completion(path: Path, detail_accessions: Iterable[str]) -> dict[str, Any]:
    """Stream the motif completion table and keep the requested accessions."""
    detail_wanted = {str(value) for value in detail_accessions}
    source = _resolve_input(path, "motif completion table")
    handle, header = _open_table(
        source, "motif completion table", list(MOTIF_COLUMNS)
    )
    detail: dict[str, dict[str, str]] = {}
    seen = 0
    with handle:
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            accession = (raw.get("accession") or "").strip()
            if not accession:
                continue
            seen += 1
            if accession in detail_wanted:
                detail[accession] = make_evidence_row(
                    accession,
                    motif_present=True,
                    **{
                        normalized: (raw.get(column) or "")
                        for column, normalized in MOTIF_COLUMNS.items()
                        if normalized != "accession"
                    },
                )
    return {"detail": detail, "rows": seen}


def read_confounder_evidence(path: Path) -> dict[str, dict[str, str]]:
    """Read the confounder table: ``{accession: {confounder_flag, confounder_cap}}``."""
    source = _resolve_input(path, "confounder candidates")
    handle, header = _open_table(
        source, "confounder candidates", list(CONFOUNDER_COLUMNS)
    )
    out: dict[str, dict[str, str]] = {}
    with handle:
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for raw in reader:
            if not raw:
                continue
            accession = (raw.get("accession") or "").strip()
            if not accession:
                continue
            out[accession] = {
                "confounder_flag": (raw.get("confounder_flag") or "").strip(),
                "confounder_cap": (raw.get("recommended_confidence_cap") or "").strip(),
            }
    return out


# ---------------------------------------------------------------------------
# Report writers
# ---------------------------------------------------------------------------


def write_input_manifest(output_dir: Path, sources: Mapping[str, Path]) -> list[dict[str, Any]]:
    """Write ``input_manifest.tsv`` with a real size + SHA-256 per input."""
    rows: list[dict[str, Any]] = []
    for name in sorted(sources):
        path = sources[name]
        bound = bound_file(path)
        rows.append(
            {
                "name": name,
                "path": str(bound["path"]),
                "size": bound["size"] if bound["size"] is not None else "pending",
                "sha256": bound["sha256"] if bound["sha256"] is not None else "pending",
                "status": bound["status"],
            }
        )
    write_tsv(output_dir / "input_manifest.tsv", MANIFEST_FIELDS, rows)
    return rows


def write_attribution_table(
    output_dir: Path, base: str, attribution_summary: Mapping[str, Any]
) -> Path:
    """Write the per-attribution counts plus the accession list per attribution."""
    by_attribution: dict[str, list[str]] = {name: [] for name in ATTRIBUTION_ORDER}
    for accession, attribution in attribution_summary["assignments"]:
        by_attribution.setdefault(attribution, []).append(accession)
    rows = []
    for name in ATTRIBUTION_ORDER:
        rule = disposition_rule(name)
        rows.append(
            {
                "attribution": name,
                "count": attribution_summary["attribution_counts"][name],
                "primary_disposition": rule["primary_disposition"],
                "status": rule["status"],
                "proving_evidence": rule["proving_evidence"],
                "accessions": ";".join(sorted(by_attribution.get(name, []))),
            }
        )
    return write_tsv(output_dir / ("%s_attribution.tsv" % base), ATTRIBUTION_FIELDS, rows)


def write_disposition_proposal(
    output_dir: Path,
    base: str,
    attribution_summary: Mapping[str, Any],
    context: Mapping[str, Any],
) -> Path:
    """Write the disposition RULE per attribution, with its proving evidence."""
    lines: list[str] = []
    lines.append("# %s residual: disposition proposal (rule, not per-accession edits)" % base)
    lines.append("")
    lines.append(
        "Candidate-only analysis of the frozen v1 candidate flow.  Attribution names "
        "the *stage* that provably dropped an accession; it is **not** a "
        "`primary_disposition`.  The rules below map each attribution onto exactly one "
        "of the six v2 dispositions.  No candidate is deleted, demoted or excluded by "
        "this document."
    )
    lines.append("")
    lines.append("## Residual and its check")
    lines.append("")
    lines.append("- residual expression: `%s`" % context["residual_expression"])
    lines.append(
        "- residual (explicit accession set difference): **%d**"
        % attribution_summary["residual_count"]
    )
    lines.append(
        "- attribution counts sum: **%d** (`sum == residual`: %s)"
        % (
            attribution_summary["sum_attribution_counts"],
            attribution_summary["attribution_counts_sum_to_residual"],
        )
    )
    lines.append(
        "- unexplained: **%d**" % attribution_summary["unexplained_count"]
    )
    lines.append(
        "- `check_hold_residual` status: **%s** (silenced: %s)"
        % (context["check_hold_residual_status"], context["check_hold_residual_silenced"])
    )
    lines.append(
        "- assignment completeness: `high_confidence ∪ hold` == admissible assignments: "
        "**%s**" % context["assignment_completeness_status"]
    )
    lines.append("")
    lines.append("## Rule per attribution")
    lines.append("")
    lines.append(
        "| attribution | count | status | target `primary_disposition` | proving evidence (file::column) |"
    )
    lines.append("|---|---:|---|---|---|")
    for name in ATTRIBUTION_ORDER:
        rule = disposition_rule(name)
        lines.append(
            "| `%s` | %d | %s | `%s` | %s |"
            % (
                name,
                attribution_summary["attribution_counts"][name],
                rule["status"],
                rule["primary_disposition"],
                rule["proving_evidence"].replace("|", "\\|"),
            )
        )
    lines.append("")
    lines.append("## Justification and action per attribution")
    lines.append("")
    for name in ATTRIBUTION_ORDER:
        rule = disposition_rule(name)
        lines.append("### `%s` -- %d accession(s)" % (name, attribution_summary["attribution_counts"][name]))
        lines.append("")
        lines.append("- vocabulary basis: %s" % rule["basis"])
        lines.append("- target `primary_disposition`: `%s`" % rule["primary_disposition"])
        lines.append("- justification: %s" % rule["justification"])
        lines.append("- proving evidence: %s" % rule["proving_evidence"])
        lines.append("- action: %s" % rule["action"])
        lines.append("")
    lines.append("## How this closes the flow")
    lines.append("")
    lines.append(
        "The rules above are inputs to F14.  Until F14 emits one "
        "`primary_disposition` per accession and passes them to "
        "`check_hold_residual` as `other_dispositions`, the check must keep "
        "reporting `mismatch`: attribution explains *why* the 616 fell out, it does "
        "not dispose of them.  De-duplicating the residual by renaming a bucket, by "
        "relaxing the expectation, or by dropping the accessions from the universe "
        "are all forbidden."
    )
    lines.append("")
    lines.append("## Evidence sources")
    lines.append("")
    for name in sorted(context["evidence_sources"]):
        bound = context["evidence_sources"][name]
        lines.append(
            "- `%s` -- %s (sha256 `%s`, %s bytes)"
            % (name, bound["path"], bound["sha256"], bound["size"])
        )
    lines.append("")
    lines.append(
        "*Candidate-only: every statement here is about candidate homology or "
        "functional potential, not about a validated PHB/PHA degradation phenotype.*"
    )
    lines.append("")
    path = output_dir / ("%s_disposition_proposal.md" % base)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_namespace(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Attribute the residual of the frozen v1 candidate flow and propose a "
            "disposition rule per attribution."
        )
    )
    parser.add_argument("--universe", type=Path, required=True, help="v1 candidate universe TSV")
    parser.add_argument("--high-confidence", type=Path, required=True)
    parser.add_argument("--hold", type=Path, required=True)
    parser.add_argument("--filter-summary", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--motif", type=Path, required=True)
    parser.add_argument("--confounders", type=Path, required=True)
    parser.add_argument(
        "--pool-external-layer",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="a pool-external layer (repeatable); may be supplied more than once",
    )
    parser.add_argument(
        "--expect-residual-count",
        type=int,
        default=None,
        help="pin the residual size; a mismatch is reported and changes the exit code",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def parse_pool_external_layers(specs: Sequence[str]) -> dict[str, Path]:
    layers: dict[str, Path] = {}
    for spec in specs or ():
        name, separator, raw_path = str(spec).partition("=")
        name = name.strip()
        if not separator or not name or not raw_path.strip():
            raise ValueError(
                "a pool-external layer must be given as NAME=PATH, got: %r" % (spec,)
            )
        if name in layers:
            raise ValueError("pool-external layer %r was supplied twice" % name)
        layers[name] = _resolve_input(Path(raw_path.strip()), "pool-external layer %s" % name)
    return layers


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Execute the attribution and write every report into ``--output-dir``."""
    paths = {
        "universe": _resolve_input(args.universe, "v1 candidate universe"),
        "high_confidence": _resolve_input(args.high_confidence, "pool-in high confidence"),
        "hold": _resolve_input(args.hold, "pool-in hold"),
        "filter_summary": _resolve_input(args.filter_summary, "filter summary"),
        "matrix": _resolve_input(args.matrix, "full-library subtype matrix"),
        "classification": _resolve_input(args.classification, "superfamily classification"),
        "motif": _resolve_input(args.motif, "motif completion table"),
        "confounders": _resolve_input(args.confounders, "confounder candidates"),
    }
    external_paths = parse_pool_external_layers(args.pool_external_layer)

    # the output directory is validated before anything else is written
    output_dir = ensure_output_dir(args.output_dir)

    universe = read_accession_list(paths["universe"], "v1 candidate universe")
    high_confidence = read_accession_list(paths["high_confidence"], "pool-in high confidence")
    hold = read_accession_list(paths["hold"], "pool-in hold candidates")
    hold_reasons = read_reason_table(paths["hold"], "hold_reason", "pool-in hold candidates")
    filter_summary = read_filter_summary(paths["filter_summary"])

    universe_set = set(universe)
    high_confidence_set = set(high_confidence)
    hold_set = set(hold)

    # Step 1: the residual is an explicit accession SET DIFFERENCE.
    residual = sorted(universe_set - high_confidence_set - hold_set)

    # The intentional loud failure is reproduced, never silenced: no
    # other_dispositions bucket is supplied.
    hold_residual = check_hold_residual(
        universe=universe_set,
        high_confidence=high_confidence_set,
        hold=hold_set,
        other_dispositions={},
        expected_residual_count=args.expect_residual_count,
        other_dispositions_source=CHECK_HOLD_RESIDUAL_OTHER_DISPOSITIONS_SOURCE,
    )
    write_json(output_dir / "check_hold_residual.json", hold_residual)

    registry = load_frozen_filter_registry()

    # Step 2: one evidence row per residual accession, plus the admissible set for
    # the completeness check.
    residual_set = set(residual)
    matrix_scan = scan_subtype_matrix(
        paths["matrix"], residual_set, admissible_universe=universe_set, registry=registry
    )
    classification_scan = scan_superfamily_classification(
        paths["classification"],
        residual_set,
        admissible_universe=universe_set,
        registry=registry,
    )
    motif_scan = scan_motif_completion(paths["motif"], residual_set)
    confounders = read_confounder_evidence(paths["confounders"])

    external_layer_accessions: dict[str, list[str]] = {}
    external_layer_sets: dict[str, set[str]] = {}
    for name, path in external_paths.items():
        accessions = read_accession_list(path, "pool-external layer %s" % name)
        external_layer_accessions[name] = accessions
        external_layer_sets[name] = set(accessions)

    universe_layers = {}
    for accession in residual_set:
        universe_layers[accession] = {}
    with paths["universe"].open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            accession = (raw.get("accession") or "").strip()
            if accession in universe_layers:
                universe_layers[accession] = {
                    "v1_layer": (raw.get("layer") or "").strip(),
                    "v1_family": (raw.get("family") or "").strip(),
                }
    if not universe_layers and residual_set:
        raise ValueError("no universe row was found for any residual accession")

    evidence_rows: list[dict[str, str]] = []
    for accession in residual:
        row = dict(matrix_scan["detail"].get(accession) or {})
        if not row:
            row = make_evidence_row(accession)
        discovery = classification_scan["detail"].get(accession)
        if discovery:
            for key in (
                "discovery_present",
                "discovery_superfamily_claim",
                "discovery_superfamily_confidence",
                "discovery_n_superfamilies_hit",
                "discovery_best_evalue",
                "discovery_model_layer",
            ):
                row[key] = discovery[key]
        motif = motif_scan["detail"].get(accession)
        if motif:
            for key in (
                "motif_present",
                "motif_candidate_prior_superfamily",
                "motif_candidate_assignment_status",
                "motif_prior_anchor_basis",
                "motif_lipase_box_state",
                "motif_sbd_pf06850_binding_state",
            ):
                row[key] = motif[key]
        confounder = confounders.get(accession)
        if confounder:
            row["confounder_flag"] = confounder["confounder_flag"]
            row["confounder_cap"] = confounder["confounder_cap"]
        row.update(universe_layers.get(accession) or {})
        evidence_rows.append(row)

    rows, attribution_summary = attribute_residual(
        residual, evidence_rows, filter_summary, external_layer_sets or None
    )

    # Step 3: prove that "no admissible assignment" is the only way out of both
    # flow tables, in both directions.
    admissible = set(matrix_scan["admissible"]) | set(classification_scan["admissible"])
    completeness = check_assignment_completeness(
        admissible=admissible, high_confidence=high_confidence_set, hold=hold_set
    )
    registry_check = check_registry_range(
        set(filter_summary["per_superfamily"]) | set(registry), registry
    )

    # Step 4: the disposition rule per attribution + the run's own summary.
    assignments = [(row["accession"], row["attribution"]) for row in rows]
    attribution_summary_with_assignments = dict(attribution_summary)
    attribution_summary_with_assignments["assignments"] = assignments

    base = "residual_%d" % len(residual)
    write_tsv(output_dir / ("%s.tsv" % base), RESIDUAL_FIELDS, rows)
    write_attribution_table(output_dir, base, attribution_summary_with_assignments)

    hold_reason_counts = {reason: len(ids) for reason, ids in sorted(hold_reasons.items())}
    hold_reason_union = {accession for ids in hold_reasons.values() for accession in ids}
    frozen_breakdown = dict(filter_summary["hold_breakdown"])
    hold_anchor = {
        "observed_hold_reason_counts_nonexclusive": hold_reason_counts,
        "frozen_filter_summary_hold_breakdown": frozen_breakdown,
        "matches_frozen_filter_summary": hold_reason_counts == frozen_breakdown,
        "observed_hold_reason_union_count": len(hold_reason_union),
        "sum_of_observed_reason_counts": sum(hold_reason_counts.values()),
        "hold_table_row_count": len(hold_set),
        "note": (
            "the four reason counts are non-exclusive; the union (not the sum) is the "
            "hold layer size"
        ),
    }
    if not hold_anchor["matches_frozen_filter_summary"]:
        raise ValueError(
            "the hold table's reason counts do not match filter_summary.json: %s vs %s"
            % (hold_reason_counts, frozen_breakdown)
        )

    aggregate_count = len(universe_set) - len(high_confidence_set) - len(hold_set)
    attribution_closes = bool(attribution_summary["attribution_counts_sum_to_residual"])
    if args.expect_residual_count is not None and len(residual) != int(
        args.expect_residual_count
    ):
        attribution_closes = False
    exit_code = exit_code_for(
        attribution_closes, str(hold_residual.get("status")), len(residual)
    )

    hold_residual_problems_source = (
        "outside_universe_count=%d; "
        "overlap_high_confidence_hold_count=%d; "
        "residual_unexplained_count=%d; "
        "residual_bucket_names=%s"
        % (
            hold_residual["outside_universe_count"],
            hold_residual["overlap_high_confidence_hold_count"],
            hold_residual["residual_unexplained_count"],
            hold_residual["residual_bucket_names"],
        )
    )

    evidence_sources = {}
    for name, path in sorted(paths.items()):
        evidence_sources[name] = bound_file(path)
    for name, path in sorted(external_paths.items()):
        evidence_sources["pool_external_layer:%s" % name] = bound_file(path)

    context = {
        "residual_expression": (
            "%d - %d - %d (explicit accession set difference)"
            % (len(universe_set), len(high_confidence_set), len(hold_set))
        ),
        "check_hold_residual_status": str(hold_residual.get("status")),
        "check_hold_residual_silenced": False,
        "assignment_completeness_status": completeness["status"],
        "evidence_sources": {
            name: {
                "path": bound["path"],
                "sha256": bound["sha256"] if bound["sha256"] is not None else "pending",
                "size": bound["size"] if bound["size"] is not None else "pending",
            }
            for name, bound in evidence_sources.items()
        },
    }
    proposal_path = write_disposition_proposal(
        output_dir, base, attribution_summary_with_assignments, context
    )

    manifest_rows = write_input_manifest(
        output_dir, {**paths, **{"pool_external_layer:%s" % k: v for k, v in external_paths.items()}}
    )

    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "tool": "attribute_phaded_candidate_flow_residual",
        "purpose": (
            "attribute the residual of the frozen v1 candidate flow to the stage that "
            "provably dropped each accession and propose a disposition rule per "
            "attribution"
        ),
        "residual_expression": context["residual_expression"],
        "universe_count": len(universe_set),
        "high_confidence_count": len(high_confidence_set),
        "hold_count": len(hold_set),
        "residual_count": len(residual),
        "aggregate_subtraction_count": aggregate_count,
        "aggregate_subtraction_agrees_with_the_set_difference": (
            aggregate_count == len(residual)
        ),
        "aggregate_subtraction_expression": (
            "universe_count - high_confidence_count - hold_count"
        ),
        "expected_residual_count": args.expect_residual_count,
        "expectation_met": (
            None
            if args.expect_residual_count is None
            else len(residual) == int(args.expect_residual_count)
        ),
        "attribution_counts": attribution_summary["attribution_counts"],
        "sum_attribution_counts": attribution_summary["sum_attribution_counts"],
        "attribution_counts_sum_to_residual": attribution_summary[
            "attribution_counts_sum_to_residual"
        ],
        "unexplained_count": attribution_summary["unexplained_count"],
        "unexplained_accessions": attribution_summary["unexplained_accessions"],
        "pending_buckets": attribution_summary["pending_buckets"],
        "attribution_status": attribution_summary["attribution_status"],
        "attribution_extensions": attribution_summary["attribution_extensions"],
        "vocabulary": attribution_summary["vocabulary"],
        "evidence_rows_total": attribution_summary["evidence_rows_total"],
        "evidence_rows_used": attribution_summary["evidence_rows_used"],
        "evidence_rows_outside_residual": attribution_summary[
            "evidence_rows_outside_residual"
        ],
        "confounder_flagged_count": attribution_summary["confounder_flagged_count"],
        "hold_reason_anchor": hold_anchor,
        "sanity_anchor": {
            "hold_reason_counts_nonexclusive": hold_reason_counts,
            "hold_reason_union_count": len(hold_reason_union),
            "residual_count": len(residual),
        },
        "check_hold_residual": hold_residual,
        "check_hold_residual_status": str(hold_residual.get("status")),
        "check_hold_residual_silenced": False,
        "check_hold_residual_other_dispositions_source": (
            CHECK_HOLD_RESIDUAL_OTHER_DISPOSITIONS_SOURCE
        ),
        "check_hold_residual_problems_source": hold_residual_problems_source,
        "assignment_completeness": completeness,
        "registry_range_check": registry_check,
        "registry_range": list(registry),
        "pool_external_layers_supplied": sorted(external_layer_accessions),
        "pool_external_layer_counts": {
            name: len(values) for name, values in sorted(external_layer_accessions.items())
        },
        "pool_external_layer_intersection_with_residual": {
            name: len(values & residual_set)
            for name, values in sorted(external_layer_sets.items())
        },
        "pool_external_layer_intersection_with_universe": {
            name: len(values & universe_set)
            for name, values in sorted(external_layer_sets.items())
        },
        "input_manifest": manifest_rows,
        "evidence_sources": context["evidence_sources"],
        "disposition_proposal": str(proposal_path),
        "exit_code": exit_code,
        "boundary": (
            "read-only attribution of frozen v1 tables: no new scoring, no registry "
            "change, no family call, no candidate deletion, no primary_disposition "
            "assignment, no modification of any existing run"
        ),
    }
    summary["outputs"] = sorted(
        path.name for path in output_dir.iterdir() if path.is_file()
    ) + [Path(base).name + "_summary.json"]
    write_json(output_dir / ("%s_summary.json" % base), summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    """Run the attribution, print the JSON summary and return the exit code.

    Exit codes: ``0`` the attribution closes (and, with a non-empty residual, the
    ``check_hold_residual`` alarm is still reported); ``1`` the attribution does
    not close (an internal or pinned-count inconsistency); ``2`` a non-empty
    residual was reported while ``check_hold_residual`` returned ``ok``, which
    means the intentional alarm was silenced.
    """
    args = build_arg_namespace(argv)
    summary = run(args)
    printable = {
        key: summary[key]
        for key in (
            "residual_expression",
            "residual_count",
            "attribution_counts",
            "sum_attribution_counts",
            "attribution_counts_sum_to_residual",
            "unexplained_count",
            "pending_buckets",
            "check_hold_residual_status",
            "check_hold_residual_silenced",
            "assignment_completeness",
            "registry_range_check",
            "pool_external_layers_supplied",
            "pool_external_layer_intersection_with_residual",
            "hold_reason_anchor",
            "exit_code",
        )
    }
    printable["outputs"] = summary.get("outputs", [])
    print(json.dumps(printable, ensure_ascii=False, sort_keys=True, indent=2))
    return int(summary["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
