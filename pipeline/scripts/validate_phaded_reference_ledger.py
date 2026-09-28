#!/usr/bin/env python3
"""Fail-closed validation for a source-bound PhaDED reference ledger.

Besides the frozen ledger contract (source-bound family assignment, FASTA/SHA-256
binding, strict-training provenance) this module validates the v2 evidence-quality
columns declared in ``pipeline/config/phaded_reference_evidence_columns.tsv`` and
reports duplicate/independence structure.

Evidence-column modes
---------------------
The frozen 723-row ledger was frozen before the evidence-model redesign and has
none of the v2 columns.  Two modes keep that ledger validatable while still
enforcing every value a row *does* declare:

``compat`` (default)
    The new columns are optional.  Columns absent from the ledger are reported as
    ``pending`` in ``report["pending_evidence_columns"]`` and nothing about the
    frozen validation result changes.  Semantic checks still run for every
    declared value, fail-closed:

    * :func:`validate_reference_row` runs on every row, so a row that claims
      grade ``E3``/``E2`` must carry the full direct-experiment provenance even in
      compat mode;
    * unknown vocabulary values (grade, assay directness, sequence integrity,
      the two eligibility booleans) are rejected;
    * a grade ``A`` row may not claim ``functional_calibration_eligible=true``;
    * ``discovery_training_eligible=true`` requires a complete sequence and a
      consistent family call;
    * a contradictory direct experimental negative blocks functional calibration
      without removing the row from the reference audit;
    * contradictions that are *reportable* rather than fatal (for example grade
      ``E1`` declared with ``assay_directness=direct``) are recorded per row in
      ``conflicts`` and surfaced in the independence summary instead of being
      silently dropped.

``strict``
    Requested with ``--evidence-columns-mode strict`` / ``--strict-evidence-columns``
    (or ``validate_ledger(..., evidence_column_mode="strict")``).  Every evidence
    column must be present in the ledger, the curated columns must carry a
    non-``pending`` value of the right vocabulary, direct records must also bind
    ``independence_group``, and declared-field conflicts become fatal.

``pending`` (and empty) is the AGENTS.md sentinel for "not curated yet"; a
``pending`` value is tolerated in compat mode only where the semantics do not
require the field, and is never inferred or fabricated.  ``not_assessed`` stays a
*declared* value (assessed, nothing found) for the free-text evidence columns.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


REQUIRED_COLUMNS = {
    "reference_id", "accession", "sequence_sha256", "sequence_source", "source_database",
    "source_database_version", "retrieval_date", "organism", "taxonomy_id",
    "phaded_superfamily", "phaded_family_id", "reported_localization", "substrate_class",
    "substrate_detail", "experimental_assay", "experimental_result", "evidence_status",
    "positive_negative_control", "catalytic_residues", "lipase_box_or_ahsmg",
    "oxyanion_hole_evidence", "domain_architecture", "primary_doi", "pmid", "pmcid", "notes",
}
ALLOWED_EVIDENCE = {
    "experimental_positive", "experimental_negative", "challenge_control",
    "annotation_only", "pending_review",
}
ALLOWED_SUPERFAMILIES = {
    "intracellular nPHASCL without lipase box",
    "intracellular nPHASCL with lipase box",
    "periplasmic PHA depolymerases",
    "intracellular nPHAMCL",
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
    "extracellular native-SCL/PhaZ7-like",
    "extracellular dPHAMCL",
}
FAMILY_COLUMNS = {"phaded_superfamily", "phaded_family_id"}
STRICT_PROVENANCE_COLUMNS = {
    "sequence_source", "source_database", "source_database_version", "retrieval_date",
    "organism", "taxonomy_id", "experimental_assay", "experimental_result",
}

# --- v2 evidence model (Task 4) -------------------------------------------------

#: Machine-readable contract for the evidence columns; the validator refuses to run
#: when the config disagrees with :data:`EVIDENCE_COLUMNS`.
EVIDENCE_COLUMNS_CONFIG = Path(__file__).resolve().parents[1] / "config" / "phaded_reference_evidence_columns.tsv"

#: Canonical evidence columns, in the order frozen by the redesign plan.
EVIDENCE_COLUMNS = (
    "accession_version", "gene_id", "study_id", "experiment_unit_id", "independence_group",
    "experimental_evidence_grade", "assay_directness", "experimental_substrate_class",
    "experimental_system", "negative_control_description", "catalytic_mutant_evidence",
    "localization_evidence", "sequence_integrity", "discovery_training_eligible",
    "functional_calibration_eligible",
)

#: E3 = purified-protein direct conversion (ideally with a catalytic mutant),
#: E2 = knockout/complement or strong in-vivo causal evidence,
#: E1 = expression, correlative phenotype or other indirect experimental link,
#: A  = database annotation, sequence inference or inherited historical family call.
EXPERIMENTAL_EVIDENCE_GRADES = ("E3", "E2", "E1", "A")
DIRECT_EXPERIMENTAL_GRADES = ("E3", "E2")
#: ``direct`` is the only directness that can carry a direct functional positive.
ASSAY_DIRECTNESS_VOCABULARY = ("direct", "indirect", "annotation")
SEQUENCE_INTEGRITY_VOCABULARY = ("complete", "partial", "truncated", "unresolved")
BOOLEAN_VOCABULARY = ("true", "false")
#: "not curated yet" sentinels; only these are tolerated by the compat semantics.
PENDING_TOKENS = frozenset({"", "pending"})
#: Fields a declared E3/E2 record must bind (plan Task 4: study, experiment unit,
#: substrate, system, a direct result and an evidence source).
DIRECT_REQUIRED_COLUMNS = (
    "study_id", "experiment_unit_id", "experimental_substrate_class", "experimental_system",
)
#: Strict mode additionally requires the experiment to be independence-grouped.
STRICT_DIRECT_REQUIRED_COLUMNS = DIRECT_REQUIRED_COLUMNS + ("independence_group",)
EVIDENCE_SOURCE_COLUMNS = ("primary_doi", "pmid", "pmcid")
#: Columns strict mode requires to be curated (non-pending) on every row.
STRICT_CURATED_COLUMNS = (
    "accession_version", "experimental_evidence_grade", "assay_directness",
    "sequence_integrity", "discovery_training_eligible", "functional_calibration_eligible",
)
#: Retained functional-calibration governance threshold; counted by unique
#: independence group, never by row count.
MINIMUM_INDEPENDENT_POSITIVE_COUNT = 3
EVIDENCE_COLUMN_MODES = ("compat", "strict")
_EVIDENCE_COLUMN_MODE_ALIASES = {
    "compat": "compat", "legacy": "compat", "default": "compat", "lenient": "compat",
    "strict": "strict", "strict-evidence": "strict", "strict-evidence-columns": "strict",
}
#: Conflict tags that are genuine evidence (a direct negative) rather than a
#: declared-field contradiction.
_NON_FATAL_CONFLICTS = frozenset({"contradictory_direct_experimental_negative"})


def _text(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    return str(value).strip()


def _norm(value: Any) -> str:
    return _text(value).lower()


def _is_pending(value: Any) -> bool:
    return _norm(value) in PENDING_TOKENS


def _is_true(value: Any) -> bool:
    return _norm(value) == "true"


def _describe(value: Any) -> str:
    text = _text(value)
    return "empty" if not text else f'"{text}"'


def _has_value(value: Any) -> bool:
    return bool(_text(value)) and not _is_pending(value)


def _regular(path: Path, label: str) -> None:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"{label} is not a regular file: {path}")


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    _regular(path, "input")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return fields, list(reader)


def _fasta(path: Path) -> dict[str, tuple[str, str]]:
    _regular(path, "reference FASTA")
    records: dict[str, tuple[str, str]] = {}
    current = None
    current_accession = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                if current in records:
                    raise ValueError(f"duplicate FASTA reference_id: {current}")
                records[current] = (current_accession or "", "".join(chunks))
            header = line[1:].strip()
            tokens = header.split("|")
            if len(tokens) != 2 or not tokens[0].strip():
                raise ValueError("FASTA header must be >reference_id|accession")
            if not tokens[1].strip():
                raise ValueError("FASTA header has empty accession")
            current = tokens[0].strip()
            current_accession = tokens[1].strip()
            chunks = []
        else:
            if current is None:
                raise ValueError("FASTA sequence appears before header")
            chunks.append(line)
    if current is not None:
        if current in records:
            raise ValueError(f"duplicate FASTA reference_id: {current}")
        records[current] = (current_accession or "", "".join(chunks))
    if not records:
        raise ValueError("reference FASTA contains no records")
    return records


# --- evidence columns config ----------------------------------------------------


def load_evidence_columns(path: str | Path | None = None) -> tuple[str, ...]:
    """Read the evidence-column contract and require it to match :data:`EVIDENCE_COLUMNS`.

    Accepts either one column name per line or a header+rows TSV whose first
    column holds the names.  Blank lines and ``#`` comments are ignored.
    """
    config = Path(path) if path is not None else EVIDENCE_COLUMNS_CONFIG
    _regular(config, "evidence columns config")
    declared: list[str] = []
    for raw in config.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        name = line.split("\t")[0].strip()
        if not name:
            continue
        if not declared and name.lower() in {"column_name", "evidence_column", "column"}:
            continue
        declared.append(name)
    if not declared:
        raise ValueError(f"evidence columns config declares no columns: {config}")
    duplicates = sorted(name for name, count in Counter(declared).items() if count > 1)
    missing = sorted(set(EVIDENCE_COLUMNS) - set(declared))
    unexpected = sorted(set(declared) - set(EVIDENCE_COLUMNS))
    if duplicates or missing or unexpected:
        raise ValueError(
            "evidence columns config mismatch: "
            + f"missing={missing} unexpected={unexpected} duplicated={duplicates}"
        )
    return tuple(declared)


def _normalize_evidence_column_mode(mode: Any, strict_evidence_columns: bool | None = None) -> str:
    text = _norm(mode) or "compat"
    if text not in _EVIDENCE_COLUMN_MODE_ALIASES:
        raise ValueError(
            f"invalid evidence_column_mode: {mode}; allowed: " + ",".join(EVIDENCE_COLUMN_MODES)
        )
    resolved = _EVIDENCE_COLUMN_MODE_ALIASES[text]
    # The boolean alias can only escalate to strict, never silently relax an
    # explicit strict request.
    if strict_evidence_columns and resolved == "compat":
        resolved = "strict"
    return resolved


# --- per-row evidence semantics -------------------------------------------------


def _row_label(row: Mapping[str, Any], index: int | None = None) -> str:
    label = _text(row.get("accession")) or _text(row.get("reference_id")) or "row"
    return f"row {index} ({label})" if index is not None else label


def _family_consistency_problem(row: Mapping[str, Any]) -> str:
    superfamily = _text(row.get("phaded_superfamily"))
    family = _text(row.get("phaded_family_id"))
    if superfamily not in ALLOWED_SUPERFAMILIES:
        return f"phaded_superfamily is {_describe(superfamily)} (not an allowed PhaDED superfamily)"
    if not _has_value(family) or family.lower() in {"unresolved", "unassigned"}:
        return f"phaded_family_id is {_describe(family)}"
    return ""


def row_independence_key(row: Mapping[str, Any]) -> tuple[str, str]:
    """Return the most specific available independence key and the column it came from.

    ``independence_group`` wins; otherwise the key degrades conservatively through
    ``study_id``+``experiment_unit_id``, ``experiment_unit_id``, ``study_id`` and
    finally ``accession_version``/``accession``.  Rows without any usable binding
    return ``("", "unresolved")`` and can never be counted as a functional positive.
    """
    group = _text(row.get("independence_group"))
    if _has_value(group):
        return group, "independence_group"
    study = _text(row.get("study_id"))
    unit = _text(row.get("experiment_unit_id"))
    if _has_value(unit):
        if _has_value(study):
            return f"{study}::{unit}", "study_id+experiment_unit_id"
        return unit, "experiment_unit_id"
    if _has_value(study):
        return study, "study_id"
    for column in ("accession_version", "accession"):
        value = _text(row.get(column))
        if _has_value(value):
            return value, column
    return "", "unresolved"


def _evidence_view(row: Mapping[str, Any]) -> dict:
    """Tolerant, non-raising view of one row's evidence fields."""
    grade_raw = _text(row.get("experimental_evidence_grade"))
    grade = "" if _is_pending(grade_raw) else grade_raw.upper()
    directness = "" if _is_pending(row.get("assay_directness")) else _norm(row.get("assay_directness"))
    integrity = "" if _is_pending(row.get("sequence_integrity")) else _norm(row.get("sequence_integrity"))
    discovery = "" if _is_pending(row.get("discovery_training_eligible")) else _norm(row.get("discovery_training_eligible"))
    calibration = "" if _is_pending(row.get("functional_calibration_eligible")) else _norm(row.get("functional_calibration_eligible"))
    status = _norm(row.get("evidence_status"))
    organism = _text(row.get("organism"))
    genus_token = organism.split()[0] if organism else ""
    genus = genus_token if genus_token else "unresolved"
    independence_group, independence_source = row_independence_key(row)
    is_direct_positive = (
        status == "experimental_positive"
        and grade in DIRECT_EXPERIMENTAL_GRADES
        and directness == "direct"
    )
    is_direct_negative = status == "experimental_negative" and (
        grade in DIRECT_EXPERIMENTAL_GRADES or directness == "direct"
    )
    conflicts: list[str] = []
    if grade in DIRECT_EXPERIMENTAL_GRADES and status not in {
        "experimental_positive", "experimental_negative", "challenge_control",
    }:
        conflicts.append(f"grade_{grade}_with_evidence_status_{status or 'empty'}")
    if grade == "E1" and directness == "direct":
        conflicts.append("grade_E1_declared_as_direct_assay")
    if grade == "A" and directness and directness != "annotation":
        conflicts.append(f"grade_A_declared_as_{directness}_assay")
    if is_direct_negative:
        conflicts.append("contradictory_direct_experimental_negative")
    if calibration == "true" and is_direct_negative:
        blocked_reason: str | None = "contradictory_direct_experimental_negative"
    elif grade == "A":
        blocked_reason = "grade_A_is_annotation_only"
    elif not is_direct_positive:
        blocked_reason = "not_a_direct_functional_positive"
    else:
        blocked_reason = None
    return {
        "reference_id": _text(row.get("reference_id")),
        "accession": _text(row.get("accession")),
        "accession_version": _text(row.get("accession_version")),
        "gene_id": _text(row.get("gene_id")),
        "sequence_sha256": _norm(row.get("sequence_sha256")),
        "study_id": _text(row.get("study_id")),
        "experiment_unit_id": _text(row.get("experiment_unit_id")),
        "family_id": _text(row.get("phaded_family_id")),
        "superfamily": _text(row.get("phaded_superfamily")),
        "organism": organism,
        "genus": genus,
        "evidence_status": status,
        "evidence_grade": grade or "pending",
        "assay_directness": directness or "pending",
        "sequence_integrity": integrity or "pending",
        "discovery_training_eligible": discovery or "pending",
        "functional_calibration_eligible": calibration or "pending",
        "independence_group": independence_group,
        "independence_group_source": independence_source,
        "pending_evidence_columns": sorted(
            column for column in EVIDENCE_COLUMNS if _is_pending(row.get(column))
        ),
        "is_direct_functional_positive": is_direct_positive,
        "counts_as_functional_positive": is_direct_positive and bool(independence_group),
        "is_contradictory_direct_negative": is_direct_negative,
        "functional_calibration_blocked_reason": blocked_reason,
        "uncurated_legacy_positive": status == "experimental_positive" and not grade,
        "conflicts": conflicts,
    }


def validate_reference_row(
    row: Mapping[str, Any], *, strict: bool = False, index: int | None = None
) -> dict:
    """Validate one row's evidence-quality fields and return a normalized view.

    Raises ``ValueError`` naming the offending field.  Fail-closed checks that run
    in both modes: vocabulary of grade/directness/integrity/eligibility, full
    provenance for a declared E3/E2 record, grade ``A`` may not claim functional
    calibration, eligibility preconditions, and a contradictory direct negative
    claiming functional calibration.  In strict mode the curated columns may not
    be ``pending``, direct records must bind ``independence_group``, and declared
    field conflicts become fatal.
    """
    label = _row_label(row, index)
    grade_raw = _text(row.get("experimental_evidence_grade"))
    if not _is_pending(grade_raw) and grade_raw.upper() not in EXPERIMENTAL_EVIDENCE_GRADES:
        raise ValueError(
            f"{label}: invalid experimental_evidence_grade: {grade_raw}; allowed: "
            + ",".join(EXPERIMENTAL_EVIDENCE_GRADES)
        )
    directness_raw = _text(row.get("assay_directness"))
    if not _is_pending(directness_raw) and _norm(directness_raw) not in ASSAY_DIRECTNESS_VOCABULARY:
        raise ValueError(
            f"{label}: invalid assay_directness: {directness_raw}; allowed: "
            + ",".join(ASSAY_DIRECTNESS_VOCABULARY)
        )
    integrity_raw = _text(row.get("sequence_integrity"))
    if not _is_pending(integrity_raw) and _norm(integrity_raw) not in SEQUENCE_INTEGRITY_VOCABULARY:
        raise ValueError(
            f"{label}: invalid sequence_integrity: {integrity_raw}; allowed: "
            + ",".join(SEQUENCE_INTEGRITY_VOCABULARY)
        )
    for column in ("discovery_training_eligible", "functional_calibration_eligible"):
        value = _text(row.get(column))
        if not _is_pending(value) and _norm(value) not in BOOLEAN_VOCABULARY:
            raise ValueError(
                f"{label}: invalid {column}: {value}; allowed: " + ",".join(BOOLEAN_VOCABULARY)
            )

    view = _evidence_view(row)
    grade = "" if view["evidence_grade"] == "pending" else view["evidence_grade"]

    if strict:
        uncurated = [column for column in STRICT_CURATED_COLUMNS if _is_pending(row.get(column))]
        if uncurated:
            raise ValueError(
                f"{label}: strict evidence mode requires curated (non-pending) values for: "
                + ",".join(uncurated)
            )

    if grade in DIRECT_EXPERIMENTAL_GRADES:
        required = STRICT_DIRECT_REQUIRED_COLUMNS if strict else DIRECT_REQUIRED_COLUMNS
        problems = [
            f"{column} is {_describe(row.get(column))}"
            for column in required
            if _is_pending(row.get(column))
        ]
        if view["assay_directness"] != "direct":
            problems.append(
                f"assay_directness is {_describe(row.get('assay_directness'))} instead of direct"
            )
        if _is_pending(row.get("experimental_result")):
            problems.append(f"experimental_result is {_describe(row.get('experimental_result'))}")
        if not any(_has_value(row.get(column)) for column in EVIDENCE_SOURCE_COLUMNS):
            problems.append("no evidence source (primary_doi/pmid/pmcid)")
        if problems:
            raise ValueError(
                f"{label}: experimental evidence grade {grade} requires a direct experimental "
                "record with complete provenance: " + "; ".join(problems)
            )

    if grade == "A" and _is_true(row.get("functional_calibration_eligible")):
        raise ValueError(
            f"{label}: functional_calibration_eligible must be false for experimental evidence "
            "grade A (database annotation cannot calibrate a functional model)"
        )

    if view["is_contradictory_direct_negative"] and _is_true(row.get("functional_calibration_eligible")):
        raise ValueError(
            f"{label}: contradictory direct experimental negative blocks functional calibration; "
            "functional_calibration_eligible must be false (the row stays in the reference audit)"
        )

    if _is_true(row.get("discovery_training_eligible")):
        blockers: list[str] = []
        if view["sequence_integrity"] != "complete":
            blockers.append(
                f"sequence_integrity is {_describe(row.get('sequence_integrity'))} instead of complete"
            )
        family_problem = _family_consistency_problem(row)
        if family_problem:
            blockers.append(family_problem)
        if blockers:
            raise ValueError(
                f"{label}: discovery_training_eligible=true requires a complete sequence with a "
                "consistent family call: " + "; ".join(blockers)
            )

    if strict:
        fatal = [conflict for conflict in view["conflicts"] if conflict not in _NON_FATAL_CONFLICTS]
        if fatal:
            raise ValueError(
                f"{label}: strict evidence mode rejects declared-field conflicts: "
                + "; ".join(fatal)
            )

    view["strict"] = strict
    return view


def summarize_independence(rows: Iterable[Mapping[str, Any]]) -> dict:
    """Report duplicate and independence structure for a set of evidence rows.

    Counts are reported per accession, accession version, sequence SHA-256, gene,
    experiment unit, study, genus, family and independence group.  Functional
    positive counts are always by unique ``independence_group`` (falling back
    conservatively through the other bindings), never by raw row count, so
    repeated accessions, repeated database records of the same gene and repeated
    reports of the same experiment cannot inflate the count.  Every row is still
    individually audited and listed in ``audited_accessions``.
    """
    views = [_evidence_view(row) for row in rows]

    def _counts(selector) -> dict[str, int]:
        counter: Counter[str] = Counter()
        for view in views:
            counter[str(selector(view) or "unresolved")] += 1
        return dict(sorted(counter.items()))

    counts_by_accession = _counts(lambda view: view["accession"])
    counts_by_accession_version = _counts(lambda view: view["accession_version"])
    counts_by_sequence_sha256 = _counts(lambda view: view["sequence_sha256"])
    counts_by_gene = _counts(lambda view: view["gene_id"])
    counts_by_experiment_unit = _counts(lambda view: view["experiment_unit_id"])
    counts_by_study = _counts(lambda view: view["study_id"])
    counts_by_genus = _counts(lambda view: view["genus"])
    counts_by_family = _counts(lambda view: view["family_id"])
    counts_by_independence_group = _counts(lambda view: view["independence_group"])
    named_counts = {
        "accession": counts_by_accession,
        "accession_version": counts_by_accession_version,
        "sequence_sha256": counts_by_sequence_sha256,
        "gene": counts_by_gene,
        "experiment_unit": counts_by_experiment_unit,
        "study": counts_by_study,
        "genus": counts_by_genus,
        "family": counts_by_family,
        "independence_group": counts_by_independence_group,
    }

    positive_views = [view for view in views if view["counts_as_functional_positive"]]
    positive_groups = sorted({view["independence_group"] for view in positive_views if view["independence_group"]})
    positive_studies = sorted({view["study_id"] or "unresolved" for view in positive_views})
    positive_genera = sorted({view["genus"] for view in positive_views})
    positive_families = sorted({view["family_id"] or "unresolved" for view in positive_views})

    def _positive_counts(selector) -> dict[str, int]:
        return dict(sorted(Counter(str(selector(view)) for view in positive_views).items()))

    grade_counts: Counter[str] = Counter(view["evidence_grade"] for view in views)
    summary = {
        "row_count": len(views),
        "audited_row_count": len(views),
        "audited_accessions": sorted({view["accession"] for view in views if view["accession"]}),
        "counts_by_accession": counts_by_accession,
        "counts_by_accession_version": counts_by_accession_version,
        "counts_by_sequence_sha256": counts_by_sequence_sha256,
        "counts_by_gene": counts_by_gene,
        "counts_by_experiment_unit": counts_by_experiment_unit,
        "counts_by_study": counts_by_study,
        "counts_by_genus": counts_by_genus,
        "counts_by_family": counts_by_family,
        "counts_by_independence_group": counts_by_independence_group,
        "unique_accessions": len(counts_by_accession),
        "unique_accession_versions": len(counts_by_accession_version),
        "unique_sequences": len(counts_by_sequence_sha256),
        "unique_genes": len(counts_by_gene),
        "unique_experiment_units": len(counts_by_experiment_unit),
        "unique_studies": len(counts_by_study),
        "unique_genera": len(counts_by_genus),
        "unique_families": len(counts_by_family),
        "unique_independence_groups": len(counts_by_independence_group),
        "duplicate_counts": {
            name: {key: value for key, value in counts.items() if value > 1}
            for name, counts in named_counts.items()
        },
        "grade_counts": dict(sorted(grade_counts.items())),
        "evidence_status_counts": _counts(lambda view: view["evidence_status"]),
        "independence_group_source_counts": _counts(lambda view: view["independence_group_source"]),
        "positive_count_basis": "unique_independence_group",
        "functional_positive_count": len(positive_groups),
        "unique_independence_group_count": len(positive_groups),
        "functional_positive_row_count": len(positive_views),
        "functional_positive_independence_groups": positive_groups,
        "functional_positive_accessions": sorted(
            {view["accession"] for view in positive_views if view["accession"]}
        ),
        "functional_positive_by_group": _positive_counts(lambda view: view["independence_group"]),
        "functional_positive_by_study": _positive_counts(lambda view: view["study_id"] or "unresolved"),
        "functional_positive_by_genus": _positive_counts(lambda view: view["genus"]),
        "functional_positive_by_family": _positive_counts(lambda view: view["family_id"] or "unresolved"),
        "functional_positive_studies": positive_studies,
        "functional_positive_genera": positive_genera,
        "functional_positive_families": positive_families,
        "functional_positive_rows_without_independence_key": sum(
            1 for view in positive_views if not view["independence_group"]
        ),
        "e1_rows": grade_counts.get("E1", 0),
        "e1_positive_rows": sum(
            1 for view in views
            if view["evidence_grade"] == "E1" and view["evidence_status"] == "experimental_positive"
        ),
        "e1_rows_with_direct_assay": sorted(
            {
                view["accession"]
                for view in views
                if view["evidence_grade"] == "E1" and view["assay_directness"] == "direct" and view["accession"]
            }
        ),
        "indirect_experimental_positive_rows": sum(
            1 for view in views
            if view["evidence_status"] == "experimental_positive"
            and not view["counts_as_functional_positive"]
        ),
        "uncurated_experimental_positive_rows": sum(
            1 for view in views if view["uncurated_legacy_positive"]
        ),
        "uncurated_experimental_positive_accessions": sorted(
            {view["accession"] for view in views if view["uncurated_legacy_positive"] and view["accession"]}
        ),
        "contradictory_direct_negative_accessions": sorted(
            {
                view["accession"]
                for view in views
                if view["is_contradictory_direct_negative"] and view["accession"]
            }
        ),
        "functional_calibration_blocked_accessions": sorted(
            {
                view["accession"]
                for view in views
                if view["functional_calibration_blocked_reason"] and view["accession"]
            }
        ),
        "functional_calibration_eligible_count": sum(
            1 for view in views if view["functional_calibration_eligible"] == "true"
        ),
        "discovery_training_eligible_count": sum(
            1 for view in views if view["discovery_training_eligible"] == "true"
        ),
        "evidence_conflicts": {
            view["accession"] or view["reference_id"] or "row": list(view["conflicts"])
            for view in views
            if view["conflicts"]
        },
    }
    return summary


def validate_ledger(
    ledger_path: str | Path,
    family_definitions_path: str | Path,
    fasta_path: str | Path,
    output_path: str | Path | None = None,
    strict_training: bool = False,
    expected_family_count: int = 38,
    evidence_column_mode: str = "compat",
    evidence_columns_path: str | Path | None = None,
    strict_evidence_columns: bool | None = None,
) -> dict:
    mode = _normalize_evidence_column_mode(evidence_column_mode, strict_evidence_columns)
    evidence_config = Path(evidence_columns_path) if evidence_columns_path is not None else EVIDENCE_COLUMNS_CONFIG
    declared_evidence_columns = load_evidence_columns(evidence_config)
    ledger_path, family_definitions_path, fasta_path = map(Path, (ledger_path, family_definitions_path, fasta_path))
    fields, rows = _read_tsv(ledger_path)
    missing = sorted(REQUIRED_COLUMNS - set(fields))
    if missing:
        raise ValueError("ledger missing required columns: " + ",".join(missing))
    pending_evidence_columns = sorted(set(declared_evidence_columns) - set(fields))
    if mode == "strict" and pending_evidence_columns:
        raise ValueError(
            "ledger missing required evidence columns: " + ",".join(pending_evidence_columns)
        )
    def_fields, definitions = _read_tsv(family_definitions_path)
    missing_defs = sorted(FAMILY_COLUMNS - set(def_fields))
    if missing_defs:
        raise ValueError("family definitions missing required columns: " + ",".join(missing_defs))
    if expected_family_count < 1:
        raise ValueError("expected_family_count must be positive")
    allowed_pairs: set[tuple[str, str]] = set()
    family_ids: set[str] = set()
    for definition_index, definition in enumerate(definitions, start=2):
        superfamily = (definition.get("phaded_superfamily") or "").strip()
        family = (definition.get("phaded_family_id") or "").strip()
        if superfamily not in ALLOWED_SUPERFAMILIES:
            raise ValueError(f"family definition row {definition_index}: invalid PhaDED superfamily")
        if not family:
            raise ValueError(f"family definition row {definition_index}: missing phaded_family_id")
        if family in family_ids:
            raise ValueError(f"duplicate phaded_family_id in definitions: {family}")
        family_ids.add(family)
        allowed_pairs.add((superfamily, family))
    if len(family_ids) != expected_family_count:
        raise ValueError(
            f"family definitions contain {len(family_ids)} unique family IDs; "
            f"expected {expected_family_count} unique family IDs"
        )
    fasta = _fasta(fasta_path)
    seen_accessions: set[str] = set()
    seen_hashes: set[str] = set()
    seen_refs: set[str] = set()
    statuses: Counter[str] = Counter()
    families: Counter[str] = Counter()
    failures: list[str] = []
    ledger_refs: set[str] = set()
    evidence_views: list[dict] = []
    for index, row in enumerate(rows, start=2):
        ref = (row.get("reference_id") or "").strip()
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError(f"row {index}: missing accession")
        if not ref:
            raise ValueError(f"row {index}: missing reference_id")
        if not row.get("primary_doi", "").strip() and not row.get("pmid", "").strip() and not row.get("pmcid", "").strip():
            raise ValueError(f"{accession}: missing primary identifier")
        if accession in seen_accessions:
            raise ValueError(f"duplicate accession: {accession}")
        seen_accessions.add(accession)
        if ref in seen_refs:
            raise ValueError(f"duplicate reference_id: {ref}")
        seen_refs.add(ref)
        digest = (row.get("sequence_sha256") or "").strip().lower()
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError(f"{accession}: invalid sequence_sha256")
        if digest in seen_hashes:
            raise ValueError(f"duplicate sequence SHA-256: {digest}")
        seen_hashes.add(digest)
        superfamily = (row.get("phaded_superfamily") or "").strip()
        family = (row.get("phaded_family_id") or "").strip()
        if superfamily not in ALLOWED_SUPERFAMILIES:
            raise ValueError(f"{accession}: invalid PhaDED superfamily")
        if (superfamily, family) not in allowed_pairs:
            raise ValueError(f"{accession}: family assignment is not source-bound")
        status = (row.get("evidence_status") or "").strip()
        if status not in ALLOWED_EVIDENCE:
            raise ValueError(f"{accession}: invalid evidence_status")
        if status in {"experimental_positive", "experimental_negative"}:
            if not row.get("experimental_assay", "").strip() or not row.get("experimental_result", "").strip():
                raise ValueError(f"{accession}: experimental call lacks assay/result")
        if strict_training and status in {"annotation_only", "pending_review"}:
            raise ValueError(f"{accession}: annotation_only/pending_review forbidden in strict-training")
        if strict_training:
            missing_provenance = sorted(
                column for column in STRICT_PROVENANCE_COLUMNS
                if not row.get(column, "").strip()
            )
            if missing_provenance:
                raise ValueError(
                    f"{accession}: strict-training missing source provenance: "
                    + ",".join(missing_provenance)
                )
        evidence_views.append(validate_reference_row(row, strict=mode == "strict", index=index))
        statuses[status] += 1
        families[family] += 1
        ledger_refs.add(ref)
    if ledger_refs != set(fasta):
        missing_fasta = sorted(ledger_refs - set(fasta))
        extra_fasta = sorted(set(fasta) - ledger_refs)
        raise ValueError(f"FASTA-to-ledger mismatch: missing={missing_fasta}, extra={extra_fasta}")
    for row in rows:
        ref = row["reference_id"].strip()
        fasta_accession, sequence = fasta[ref]
        if fasta_accession != row["accession"].strip():
            raise ValueError(
                f"{ref}: FASTA accession mismatch: header={fasta_accession}, "
                f"ledger={row['accession'].strip()}"
            )
        observed = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if observed != row["sequence_sha256"].strip().lower():
            raise ValueError(f"{row['accession']}: FASTA sequence SHA-256 mismatch")
    independence_summary = summarize_independence(rows)
    # Every family in the ledger is reported, so a family with zero curated
    # independent positives is visible instead of silently omitted.
    family_positive_groups: dict[str, set[str]] = {family: set() for family in families}
    family_positive_genera: dict[str, set[str]] = {family: set() for family in families}
    for view in evidence_views:
        if not view["counts_as_functional_positive"]:
            continue
        family_key = view["family_id"] or "unresolved"
        family_positive_groups.setdefault(family_key, set()).add(view["independence_group"])
        family_positive_genera.setdefault(family_key, set()).add(view["genus"])
    family_functional_positive_counts = {
        family_key: len(groups) for family_key, groups in sorted(family_positive_groups.items())
    }
    family_distinct_genus_counts = {
        family_key: len(genera) for family_key, genera in sorted(family_positive_genera.items())
    }
    report = {
        "ledger": str(ledger_path), "family_definitions": str(family_definitions_path), "reference_fasta": str(fasta_path),
        "row_count": len(rows), "fasta_count": len(fasta),
        "evidence_status_counts": dict(sorted(statuses.items())),
        "family_counts": dict(sorted(families.items())), "family_definition_count": len(family_ids),
        "expected_family_count": expected_family_count, "duplicate_accessions": [],
        "validation_failures": failures, "strict_training": strict_training, "status": "valid",
        "source_file_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (ledger_path, family_definitions_path, fasta_path)
        },
        "evidence_column_mode": mode,
        "evidence_columns_config": str(evidence_config),
        "evidence_columns_config_sha256": hashlib.sha256(evidence_config.read_bytes()).hexdigest(),
        "evidence_columns_declared": sorted(set(fields) & set(declared_evidence_columns)),
        "pending_evidence_columns": pending_evidence_columns,
        "evidence_grade_counts": independence_summary["grade_counts"],
        "evidence_conflict_rows": independence_summary["evidence_conflicts"],
        "functional_positive_count": independence_summary["functional_positive_count"],
        "functional_positive_row_count": independence_summary["functional_positive_row_count"],
        "functional_calibration_gate": {
            "minimum_independent_positive_count": MINIMUM_INDEPENDENT_POSITIVE_COUNT,
            "positive_count_basis": "unique_independence_group",
            "minimum_distinct_genus_count": None,
            "genus_breadth_basis": "reported_not_enforced",
            "family_functional_positive_counts": family_functional_positive_counts,
            "family_distinct_genus_counts": family_distinct_genus_counts,
            "families_meeting_minimum": sorted(
                family_key for family_key, count in family_functional_positive_counts.items()
                if count >= MINIMUM_INDEPENDENT_POSITIVE_COUNT
            ),
            "families_below_minimum": sorted(
                family_key for family_key, count in family_functional_positive_counts.items()
                if count < MINIMUM_INDEPENDENT_POSITIVE_COUNT
            ),
        },
        "independence_summary": independence_summary,
    }
    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger")
    parser.add_argument("family_definitions")
    parser.add_argument("reference_fasta")
    parser.add_argument("--output")
    parser.add_argument("--strict-training", action="store_true")
    parser.add_argument("--expected-family-count", type=int, default=38)
    parser.add_argument(
        "--evidence-columns-mode", choices=EVIDENCE_COLUMN_MODES, default="compat",
        help="compat (default) tolerates the pre-redesign ledger; strict requires curated columns",
    )
    parser.add_argument(
        "--strict-evidence-columns", action="store_true",
        help="alias that escalates --evidence-columns-mode to strict",
    )
    parser.add_argument(
        "--evidence-columns-config", default=None,
        help=f"evidence column contract (default: {EVIDENCE_COLUMNS_CONFIG})",
    )
    parser.add_argument(
        "--independence-summary-json", nargs="?", const="-", default=None, metavar="PATH",
        help="print the duplicate/independence summary as JSON ('-' or no PATH prints to stdout)",
    )
    args = parser.parse_args(argv)
    report = validate_ledger(
        args.ledger, args.family_definitions, args.reference_fasta, args.output,
        args.strict_training, args.expected_family_count,
        evidence_column_mode=(
            "strict" if args.strict_evidence_columns else args.evidence_columns_mode
        ),
        evidence_columns_path=args.evidence_columns_config,
    )
    target = args.independence_summary_json
    if target is None:
        print(f"validated {report['row_count']} PhaDED reference rows")
        return 0
    payload = json.dumps(
        report["independence_summary"], indent=2, ensure_ascii=False, sort_keys=True
    ) + "\n"
    if target == "-":
        # JSON-only stdout keeps the flag pipeable into `python -m json.tool`.
        print(payload, end="")
        return 0
    output = Path(target)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(payload, encoding="utf-8")
    print(f"wrote independence summary to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
