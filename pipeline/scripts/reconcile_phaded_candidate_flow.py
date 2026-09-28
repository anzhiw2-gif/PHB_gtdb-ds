#!/usr/bin/env python3
"""Close the PhaDED candidate flow ledger and reconcile every count discrepancy.

This module is the machine-readable answer to acceptance boundary item 5 of the
2026-09-28 evidence-model redesign: *the 616 unexplained flow records, the
4,301/6,214 set relations and the gRodon 66-record difference each have a
machine-readable explanation.*

Two layers, deliberately separated:

``reconcile``
    The **ledger**. One row per accession of the frozen candidate universe with
    exactly one ``primary_disposition`` from the six primary dispositions of the
    redesign spec, plus all applicable hold / demotion / deferred reasons. The
    summary reports ``unaccounted`` and ``multiply_disposed`` separately, so a
    candidate that is missing a disposition can never be silently hidden inside a
    total.

``check_*``
    The **named reconciliation checks**. Each one takes accession-level sets (or
    the raw input rows) and proves a documented count by an explicit set
    difference. No check may be satisfied by arithmetic on aggregates: a check
    that is given only totals is not implemented here, precisely because such a
    check cannot fail.

Failure contract (one rule for every check, no exceptions)
---------------------------------------------------------
Every ``check_*`` function returns a machine-readable ``dict`` that always
carries ``status``: ``"ok"`` or ``"mismatch"``. A check **never** raises for a
scientific inconsistency; it reports ``"mismatch"`` plus the offending
accessions. Two callers turn that into a loud failure:

* :func:`assert_check_ok` raises ``ValueError`` with the check name and the
  explanation -- use it in tests and in any gate that must stop;
* the CLI exits non-zero (``1`` for a mismatching check, ``2`` for a flow that
  does not close) after writing the reports, so the evidence is never lost.

Structural problems (a missing column, a malformed table, a missing path, a
non-positive expectation) raise ``ValueError`` immediately: they are input
errors, not findings.

Boundaries
----------
Candidate accounting only. This module writes nothing outside the caller's
``--output-dir``, creates no analysis run and no dated deploy directory, never
touches ``runs/``, ``results/`` or ``deploy/``, and never modifies a candidate,
a profile or a count. It reads evidence; it does not produce any.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

SCHEMA_VERSION = "1.0"

# ---------------------------------------------------------------------------
# Disposition vocabulary (redesign spec, section "候选目录" / candidate catalog)
# ---------------------------------------------------------------------------

DISPOSITION_CORE_SEQUENCE_HOMOLOG = "core_sequence_homolog"
DISPOSITION_PROBABLE_SEQUENCE_HOMOLOG = "probable_sequence_homolog"
DISPOSITION_REMOTE_HOMOLOG_CANDIDATE = "remote_homolog_candidate"
DISPOSITION_FUNCTION_UNRESOLVED = "function_unresolved"
DISPOSITION_DEFERRED_STRUCTURE_REVIEW = "deferred_structure_review"
DISPOSITION_EXCLUDED_INPUT_QUALITY = "excluded_input_quality"

#: The six primary dispositions. A v2 candidate has exactly one of these.
DISPOSITIONS: tuple[str, ...] = (
    DISPOSITION_CORE_SEQUENCE_HOMOLOG,
    DISPOSITION_PROBABLE_SEQUENCE_HOMOLOG,
    DISPOSITION_REMOTE_HOMOLOG_CANDIDATE,
    DISPOSITION_FUNCTION_UNRESOLVED,
    DISPOSITION_DEFERRED_STRUCTURE_REVIEW,
    DISPOSITION_EXCLUDED_INPUT_QUALITY,
)

#: The legacy v1 vocabulary the projection can read. ``high_confidence`` is a v1
#: selection flag, never a v2 disposition.
LEGACY_DISPOSITION_HIGH_CONFIDENCE = "high_confidence"
LEGACY_DISPOSITION_HOLD = "hold"
HIGH_CONFIDENCE_DISPOSITIONS: tuple[str, ...] = (
    DISPOSITION_CORE_SEQUENCE_HOMOLOG,
    DISPOSITION_PROBABLE_SEQUENCE_HOMOLOG,
)

#: Hold / demotion / deferred reason vocabulary observed in the frozen evidence
#: (``runs/20260919_phaded_nucleophile_unify_01`` and
#: ``runs/20260920_phaded_three_gaps_01/results/cys_demotion``).  The ledger keeps
#: reasons as free strings so a new reason is never silently dropped, but the
#: known reasons are pinned here for review.
HOLD_REASONS: tuple[str, ...] = (
    "localization_conflict",
    "architecture_criteria",
    "common_criteria",
    "nucleophile_conflict",
)
DEMOTE_REASONS: tuple[str, ...] = (
    "localization_conflict_signalp_predicted_export",
)
DEFERRED_REASONS: tuple[str, ...] = (
    "with_lipase_deferred_structure_review",
)

#: The frozen hold split of the 2026-09-21 review, kept as a comment-level
#: expectation only; the checks recompute it, they never hard-code it.
REVIEW_HOLD_BREAKDOWN = {
    "localization_conflict": 48488,
    "architecture_criteria": 19562,
    "common_criteria": 2538,
    "nucleophile_conflict": 1272,
}
#: With-lipase (``DED_hfam_2``) pool-external hits: a deferred structure layer that
#: must stay outside every candidate count and must never be deleted.
REVIEW_WITH_LIPASE_DEFERRED = 1206655

REVIEW_UNIVERSE = 109087
REVIEW_POOL_HIGH_CONFIDENCE = 36611
REVIEW_POOL_HIGH_CONFIDENCE_AFTER_DEMOTION = 36559
REVIEW_HOLD = 71860
REVIEW_MERGED_HIGH_CONFIDENCE = 40690
REVIEW_POOL_EXTERNAL_STRONG = 4131
REVIEW_CYS_DEMOTED = 52
REVIEW_PHAC_FUNNEL_HITS = 1970
REVIEW_IPHAZ_FUNNEL_HITS = 4244
REVIEW_FUNNEL_UNION = 4301
REVIEW_FUNNEL_SUM_OF_REASON_COUNTS = (
    REVIEW_PHAC_FUNNEL_HITS + REVIEW_IPHAZ_FUNNEL_HITS  # 6214
)
REVIEW_GRODON_ELIGIBLE_POSITIVES = 4507
REVIEW_GRODON_MANIFEST_POSITIVES = 4441
REVIEW_GRODON_DELTA = 66
REVIEW_SUPERFAMILIES = 8
REVIEW_FAMILIES = 38
REVIEW_PROFILES = 46
REVIEW_TRAINED_PROFILES = 9
REVIEW_REFERENCE_ONLY_PROFILES = 37

#: ``109087 - 36611 - 71860`` -- the residual set of the 2026-09-21 review.
HOLD_RESIDUAL_EXPRESSION = "%d - %d - %d" % (
    REVIEW_UNIVERSE,
    REVIEW_POOL_HIGH_CONFIDENCE,
    REVIEW_HOLD,
)
RESIDUAL_EXPRESSION = HOLD_RESIDUAL_EXPRESSION

# ---------------------------------------------------------------------------
# Column contracts
# ---------------------------------------------------------------------------

ACCESSION_COLUMN_CANDIDATES = (
    "accession",
    "protein_id",
    "genome_id",
    "target",
    "identifier",
)

RECORDS_REQUIRED_COLUMNS = ("accession", "primary_disposition")
RECORDS_OPTIONAL_COLUMNS = (
    "hold_reasons",
    "demote_reasons",
    "defer_reasons",
    "deferred",
    "disposition_basis",
)

AUTHORITY_REQUIRED_COLUMNS = ("layer_id", "layer_kind")
FAMILY_DEFINITION_REQUIRED_COLUMNS = ("phaded_superfamily", "phaded_family_id")
PROFILE_REQUIRED_COLUMNS = ("profile_id", "model_status")

#: ``profile_kind`` is optional in the profile manifest; when present it must
#: agree with the authority table's ``layer_kind``.
PROFILE_KIND_COLUMN = "profile_kind"
LAYER_ID_COLUMN = "layer_id"
FAMILY_ID_COLUMN = "phaded_family_id"
SUPERFAMILY_COLUMN = "phaded_superfamily"

MODEL_STATUS_TRAINED = "trained"
MODEL_STATUS_REFERENCE_ONLY = "reference_only"
MODEL_STATUSES: tuple[str, ...] = (MODEL_STATUS_TRAINED, MODEL_STATUS_REFERENCE_ONLY)

#: Motif criteria and the columns a completion table must carry for them.  A
#: ``(state, assignability)`` pair means "the run wrote a state for this
#: criterion" -> covered; ``assignability == "assignable"`` -> resolvable.
MOTIF_CRITERIA_COLUMNS: tuple[tuple[str, str, str], ...] = (
    ("lipase_box", "lipase_box_state", "lipase_box_assignability"),
    ("his", "his_state", "his_assignability"),
    ("asp", "asp_state", "asp_assignability"),
    ("oxyanion_hole", "oxyanion_hole_state", "oxyanion_hole_assignability"),
    ("ahsmg", "ahsmg_state", "ahsmg_assignability"),
    ("lid", "lid_state", "lid_assignability"),
    ("linker", "linker_state", "linker_assignability"),
    ("sbd", "sbd_state", "sbd_assignability"),
)
MOTIF_ASSIGNABLE = "assignable"
MOTIF_NOT_ASSIGNABLE = "not_assignable"

STATUS_OK = "ok"
STATUS_MISMATCH = "mismatch"

EXIT_OK = 0
EXIT_CHECK_MISMATCH = 1
EXIT_FLOW_DOES_NOT_CLOSE = 2

FLOW_FIELDS = [
    "accession",
    "primary_disposition",
    "hold_reasons",
    "demote_reasons",
    "defer_reasons",
    "hold_or_demote",
    "demoted_from_v1",
    "deferred",
    "record_count",
    "disposition_conflict",
    "disposition_basis",
    "notes",
]

OVERLAP_FIELDS = [
    "reason_a",
    "reason_b",
    "count_a",
    "count_b",
    "sum_of_counts",
    "overlap_count",
    "union_count",
    "jaccard",
    "overlap_accessions",
    "accessions_only_in_a",
    "accessions_only_in_b",
]

DEFAULT_ACCESSION_SAMPLE = 5


# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------


def _as_set(values: Iterable[str] | None) -> set[str]:
    return {str(value) for value in values or () if str(value) != ""}


def _sorted_sample(values: Iterable[str], limit: int = DEFAULT_ACCESSION_SAMPLE) -> list[str]:
    return sorted(values)[:limit]


def _split_reasons(value: Any) -> list[str]:
    """Split one cell into a de-duplicated, order-preserving reason list."""
    if value is None:
        return []
    if isinstance(value, (list, tuple, set, frozenset)):
        parts = [str(item) for item in value]
    else:
        text = str(value).strip()
        if not text:
            return []
        parts = [part.strip() for part in text.replace(",", ";").split(";")]
    out: list[str] = []
    seen: set[str] = set()
    for part in parts:
        if part and part not in seen:
            seen.add(part)
            out.append(part)
    return out


def _parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    return text in {"1", "true", "yes", "y", "t"}


def _require_columns(rows: Sequence[Mapping[str, Any]], columns: Sequence[str], what: str) -> None:
    if not rows:
        raise ValueError("%s is empty" % what)
    missing = [column for column in columns if column not in rows[0]]
    if missing:
        raise ValueError("%s is missing columns: %s" % (what, ",".join(missing)))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve_input(path: Path, what: str) -> Path:
    """Return ``path`` if it is an existing regular file, else raise.

    A missing path is always an error: an absent input must never be interpreted
    as an empty input.
    """
    candidate = Path(path)
    if not candidate.exists():
        raise ValueError("%s path does not exist: %s" % (what, candidate))
    if not candidate.is_file():
        raise ValueError("%s path is not a regular file: %s" % (what, candidate))
    return candidate


def ensure_output_dir(path: Path) -> Path:
    """Create ``path``, or accept it when it is an existing empty directory.

    An existing non-empty directory is refused so that this tool can never write
    into an evidence directory.
    """
    out_dir = Path(path)
    if out_dir.exists():
        if not out_dir.is_dir():
            raise ValueError("output path is not a directory: %s" % out_dir)
        if any(out_dir.iterdir()):
            raise ValueError("refusing to write into a non-empty directory: %s" % out_dir)
    else:
        out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> Path:
    path = Path(path)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")
    return path


def write_json(path: Path, payload: Any) -> Path:
    path = Path(path)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# Readers (strict: a missing column or an empty table is an error)
# ---------------------------------------------------------------------------


def _looks_like_accession(token: str) -> bool:
    """True when a bare token looks like a protein/genome identifier.

    Real GTDB accessions always carry a digit, a ``|`` or a ``.``
    (``GCA_000012145.1|CP000053.1_249``, ``EAO52570.1``), while column names
    (``accession``, ``protein_id``) carry none.  This is what distinguishes a
    one-column table with a header from a headerless one-identifier-per-line list.
    """
    return any(character.isdigit() or character in "|._" for character in token)


def _read_table(path: Path, what: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV with a header row, or a headerless single-column id list.

    A file whose first line carries no tab is treated as a headerless
    one-identifier-per-line list and is reported with an empty header.  A file
    whose first line is a bare column name (``accession``, ``protein_id``, ...) is
    read as a one-column table with a header; if it then has no data row it is an
    error, so an empty input can never masquerade as an empty-but-valid result.
    """
    source = _resolve_input(path, what)
    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            first = next(reader)
        except StopIteration:
            raise ValueError("%s is empty: %s" % (what, source)) from None
        rows: list[dict[str, str]] = []
        if len(first) > 1:
            header = [column.strip() for column in first]
            for values in reader:
                if not values or all(value == "" for value in values):
                    continue
                rows.append(
                    {
                        header[index]: (values[index] if index < len(values) else "")
                        for index in range(len(header))
                    }
                )
        elif len(first) == 1:
            value = first[0].strip()
            if value.lower() in ACCESSION_COLUMN_CANDIDATES:
                header = [ACCESSION_COLUMN_CANDIDATES[0]]
                for row in reader:
                    if not row or not row[0].strip():
                        continue
                    rows.append({header[0]: row[0].strip()})
            else:
                header = []
                values = ([value] if value else []) + [
                    row[0].strip() for row in reader if row and row[0].strip()
                ]
                rows = [{ACCESSION_COLUMN_CANDIDATES[0]: token} for token in values]
        else:
            raise ValueError("%s has an unreadable first row: %s" % (what, source))
        if not rows:
            raise ValueError("%s has no data rows: %s" % (what, source))
    return header, rows


def _pick_accession_column(header: Sequence[str], what: str) -> str:
    for candidate in ACCESSION_COLUMN_CANDIDATES:
        if candidate in header:
            return candidate
    raise ValueError(
        "%s has no accession column (looked for %s); found: %s"
        % (what, ", ".join(ACCESSION_COLUMN_CANDIDATES), ", ".join(header))
    )


def read_accession_list(path: Path, what: str = "accession list") -> list[str]:
    """Return the sorted unique accessions of ``path``.

    Accepts a TSV with an ``accession`` / ``protein_id`` / ``genome_id`` column or
    a headerless one-identifier-per-line file.  Missing paths, directories and
    empty files are errors.
    """
    header, rows = _read_table(path, what)
    if not header:
        column = ACCESSION_COLUMN_CANDIDATES[0]
        accessions = {row[column] for row in rows if row.get(column)}
    else:
        column = _pick_accession_column(header, what)
        accessions = {row[column] for row in rows if row.get(column)}
    if not accessions:
        raise ValueError("%s has no usable accession value: %s" % (what, path))
    return sorted(accessions)


def read_reason_table(
    path: Path, reason_column: str, what: str | None = None
) -> dict[str, list[str]]:
    """Group accessions by reason: ``{reason: [accession, ...]}``.

    ``reason_column`` defaults to nothing -- pass e.g. ``hold_reason``.  A row
    whose reason cell is empty is an error: an unlabelled reason would silently
    disappear from every nonexcusive count.
    """
    label = what or ("%s table" % reason_column)
    header, rows = _read_table(path, label)
    if not header:
        raise ValueError("%s needs a header row with an accession and a reason column" % label)
    accession_column = _pick_accession_column(header, label)
    if reason_column not in header:
        raise ValueError(
            "%s is missing the reason column %r; found: %s"
            % (label, reason_column, ", ".join(header))
        )
    grouped: dict[str, list[str]] = {}
    for row in rows:
        accession = row.get(accession_column, "")
        if not accession:
            raise ValueError("%s has a row without an accession" % label)
        reason = (row.get(reason_column) or "").strip()
        if not reason:
            raise ValueError("%s has an empty %s for %s" % (label, reason_column, accession))
        grouped.setdefault(reason, []).append(accession)
    if not grouped:
        raise ValueError("%s has no usable reason row: %s" % (label, path))
    return {reason: sorted(set(accessions)) for reason, accessions in sorted(grouped.items())}


def read_records_table(path: Path, what: str = "candidate records") -> list[dict[str, Any]]:
    """Read the candidate record ledger (one row per accession).

    Required columns: ``accession``, ``primary_disposition``.
    Optional columns: ``hold_reasons``, ``demote_reasons``, ``defer_reasons``,
    ``deferred``, ``disposition_basis``.  Reason cells are ``;``-separated.
    """
    header, rows = _read_table(path, what)
    if not header:
        raise ValueError(
            "%s needs a header row with %s"
            % (what, ", ".join(RECORDS_REQUIRED_COLUMNS))
        )
    missing = [column for column in RECORDS_REQUIRED_COLUMNS if column not in header]
    if missing:
        raise ValueError("%s is missing columns: %s" % (what, ",".join(missing)))
    records: list[dict[str, Any]] = []
    for row in rows:
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError("%s has a row without an accession" % what)
        records.append(
            {
                "accession": accession,
                "primary_disposition": (row.get("primary_disposition") or "").strip(),
                "hold_reasons": _split_reasons(row.get("hold_reasons")),
                "demote_reasons": _split_reasons(row.get("demote_reasons")),
                "defer_reasons": _split_reasons(row.get("defer_reasons")),
                "deferred": _parse_bool(row.get("deferred")),
                "disposition_basis": (row.get("disposition_basis") or "").strip(),
            }
        )
    return records


def read_authority_rows(
    path: Path, what: str = "classification authority"
) -> list[dict[str, str]]:
    """Read ``phaded_classification_authority.tsv``; ``layer_id``, ``layer_kind`` required."""
    header, rows = _read_table(path, what)
    if not header:
        raise ValueError("%s needs a header row" % what)
    missing = [column for column in AUTHORITY_REQUIRED_COLUMNS if column not in header]
    if missing:
        raise ValueError("%s is missing columns: %s" % (what, ",".join(missing)))
    return rows


def read_family_definition_rows(
    path: Path, what: str = "family definitions"
) -> list[dict[str, str]]:
    """Read ``phaded_family_definitions.tsv``; superfamily and family id required."""
    header, rows = _read_table(path, what)
    if not header:
        raise ValueError("%s needs a header row" % what)
    missing = [
        column for column in FAMILY_DEFINITION_REQUIRED_COLUMNS if column not in header
    ]
    if missing:
        raise ValueError("%s is missing columns: %s" % (what, ",".join(missing)))
    return rows


def read_profile_rows(path: Path, what: str = "profile manifest") -> list[dict[str, str]]:
    """Read ``profile_manifest.tsv``; ``profile_id`` and ``model_status`` required."""
    header, rows = _read_table(path, what)
    if not header:
        raise ValueError("%s needs a header row" % what)
    missing = [column for column in PROFILE_REQUIRED_COLUMNS if column not in header]
    if missing:
        raise ValueError("%s is missing columns: %s" % (what, ",".join(missing)))
    return rows


def read_motif_layer(path: Path, what: str = "motif completion table") -> dict[str, Any]:
    """Read a motif completion table and derive coverage versus resolvability.

    Returns ``{"covered": set, "resolvable": set, "by_criterion": {...}}`` where
    ``covered`` is every accession the motif run reached for at least one
    criterion and ``resolvable`` is every accession for which at least one
    criterion is ``assignable``.  A state column is required for each criterion
    listed in :data:`MOTIF_CRITERIA_COLUMNS`.
    """
    header, rows = _read_table(path, what)
    if not header:
        raise ValueError("%s needs a header row" % what)
    accession_column = _pick_accession_column(header, what)
    present = [
        (criterion, state_column, assignability_column)
        for criterion, state_column, assignability_column in MOTIF_CRITERIA_COLUMNS
        if state_column in header
    ]
    if not present:
        raise ValueError(
            "%s carries no motif criterion state column (looked for %s)"
            % (what, ", ".join(column for _, column, _ in MOTIF_CRITERIA_COLUMNS))
        )
    covered: set[str] = set()
    resolvable: set[str] = set()
    by_criterion: dict[str, dict[str, set[str]]] = {}
    for criterion, state_column, assignability_column in present:
        criterion_covered: set[str] = set()
        criterion_resolvable: set[str] = set()
        for row in rows:
            accession = (row.get(accession_column) or "").strip()
            if not accession:
                continue
            state = (row.get(state_column) or "").strip()
            assignability = (row.get(assignability_column) or "").strip()
            if state or assignability:
                criterion_covered.add(accession)
                covered.add(accession)
            if assignability == MOTIF_ASSIGNABLE:
                criterion_resolvable.add(accession)
                resolvable.add(accession)
        by_criterion[criterion] = {
            "covered": criterion_covered,
            "resolvable": criterion_resolvable,
        }
    return {"covered": covered, "resolvable": resolvable, "by_criterion": by_criterion}


# ---------------------------------------------------------------------------
# Ledger: the candidate flow
# ---------------------------------------------------------------------------


def to_record(accession: str, disposition: str, **extra: Any) -> dict[str, Any]:
    """Build one ledger record. ``hold_reasons`` etc. accept a list or a ``;``-string."""
    record: dict[str, Any] = {
        "accession": str(accession),
        "primary_disposition": str(disposition or ""),
        "hold_reasons": _split_reasons(extra.get("hold_reasons")),
        "demote_reasons": _split_reasons(extra.get("demote_reasons")),
        "defer_reasons": _split_reasons(extra.get("defer_reasons")),
        "deferred": _parse_bool(extra.get("deferred", False)),
        "disposition_basis": str(extra.get("disposition_basis", "") or ""),
    }
    return record


def _normalize_record(raw: Any) -> dict[str, Any]:
    if isinstance(raw, Mapping):
        accession = str(raw.get("accession", "") or "")
        disposition = str(raw.get("primary_disposition", "") or "")
        if not accession:
            raise ValueError("candidate record without an accession: %r" % (dict(raw),))
        if disposition and disposition not in DISPOSITIONS:
            raise ValueError(
                "unknown primary_disposition %r for %s; allowed: %s"
                % (disposition, accession, ", ".join(DISPOSITIONS))
            )
        return {
            "accession": accession,
            "primary_disposition": disposition,
            "hold_reasons": _split_reasons(raw.get("hold_reasons")),
            "demote_reasons": _split_reasons(raw.get("demote_reasons")),
            "defer_reasons": _split_reasons(raw.get("defer_reasons")),
            "deferred": _parse_bool(raw.get("deferred", False)),
            "disposition_basis": str(raw.get("disposition_basis", "") or ""),
        }
    raise ValueError("candidate record must be a mapping, got %r" % (type(raw).__name__,))


def reconcile(
    universe: Iterable[str],
    records: Iterable[Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build the candidate flow ledger for ``universe``.

    Returns ``(rows, summary)``.  ``rows`` has one entry per universe accession
    (sorted by accession) carrying its single ``primary_disposition`` plus all
    applicable hold / demotion / deferred reasons.  ``summary`` separately
    reports ``unaccounted`` (no valid disposition at all) and
    ``multiply_disposed`` (more than one *distinct* disposition claimed), so that
    a missing disposition can never be hidden inside a total.

    Raises ``ValueError`` when the universe is empty or duplicated, when a record
    falls outside the universe, or when a disposition is not one of the six
    primary dispositions of the redesign spec.
    """
    universe_list = [str(value) for value in universe]
    if not universe_list:
        raise ValueError("candidate universe is empty")
    if len(set(universe_list)) != len(universe_list):
        duplicates = sorted(
            value for value, count in Counter(universe_list).items() if count > 1
        )
        raise ValueError(
            "candidate universe has duplicate accessions: %s"
            % _sorted_sample(duplicates)
        )
    universe_set = set(universe_list)

    normalized = [_normalize_record(raw) for raw in records]
    outside = {record["accession"] for record in normalized} - universe_set
    if outside:
        raise ValueError(
            "candidate record outside the universe: %s" % _sorted_sample(outside)
        )

    by_accession: dict[str, list[dict[str, Any]]] = {}
    for record in normalized:
        by_accession.setdefault(record["accession"], []).append(record)

    rows: list[dict[str, Any]] = []
    disposition_counts: Counter[str] = Counter()
    hold_reason_counts: Counter[str] = Counter()
    demote_reason_counts: Counter[str] = Counter()
    defer_reason_counts: Counter[str] = Counter()
    unaccounted: list[str] = []
    multiply_disposed: list[str] = []
    demoted_accessions: list[str] = []
    hold_accessions: set[str] = set()
    high_confidence_accessions: set[str] = set()
    deferred_accessions: set[str] = set()

    for accession in sorted(universe_set):
        group = by_accession.get(accession, [])
        claimed = [
            record["primary_disposition"]
            for record in group
            if record["primary_disposition"]
        ]
        distinct = sorted(set(claimed))
        disposition = distinct[0] if distinct else ""
        if len(distinct) > 1:
            multiply_disposed.append(accession)
        if not distinct:
            unaccounted.append(accession)

        hold_reasons = sorted(
            {reason for record in group for reason in record["hold_reasons"]}
        )
        demote_reasons = sorted(
            {reason for record in group for reason in record["demote_reasons"]}
        )
        defer_reasons = sorted(
            {reason for record in group for reason in record["defer_reasons"]}
        )
        deferred = any(record["deferred"] for record in group) or bool(defer_reasons)
        basis = next(
            (record["disposition_basis"] for record in group if record["disposition_basis"]),
            "",
        )

        disposition_counts[disposition] += 1
        for reason in hold_reasons:
            hold_reason_counts[reason] += 1
        for reason in demote_reasons:
            demote_reason_counts[reason] += 1
        for reason in defer_reasons:
            defer_reason_counts[reason] += 1
        if hold_reasons:
            hold_accessions.add(accession)
        if demote_reasons:
            demoted_accessions.append(accession)
        if disposition in HIGH_CONFIDENCE_DISPOSITIONS:
            high_confidence_accessions.add(accession)
        if deferred:
            deferred_accessions.add(accession)

        hold_or_demote = (
            "hold+demote" if hold_reasons and demote_reasons
            else "hold" if hold_reasons
            else "demote" if demote_reasons
            else ""
        )
        rows.append(
            {
                "accession": accession,
                "primary_disposition": disposition,
                "hold_reasons": ";".join(hold_reasons),
                "demote_reasons": ";".join(demote_reasons),
                "defer_reasons": ";".join(defer_reasons),
                "hold_or_demote": hold_or_demote,
                "demoted_from_v1": "true" if demote_reasons else "false",
                "deferred": "true" if deferred else "false",
                "record_count": len(group),
                "disposition_conflict": ";".join(distinct) if len(distinct) > 1 else "",
                "disposition_basis": basis,
                "notes": (
                    "no_primary_disposition_recorded" if not distinct else ""
                ),
            }
        )

    disposition_partition = {
        disposition: disposition_counts.get(disposition, 0) for disposition in DISPOSITIONS
    }
    observed_only = {
        disposition: count
        for disposition, count in sorted(disposition_counts.items())
        if disposition and disposition not in DISPOSITIONS
    }
    if observed_only:
        raise ValueError("unknown disposition reached the ledger: %r" % observed_only)

    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "universe_count": len(universe_set),
        "primary_disposition_counts": {
            disposition: count
            for disposition, count in disposition_partition.items()
            if count
        },
        "primary_disposition_counts_all_six": disposition_partition,
        "unaccounted": len(unaccounted),
        "unaccounted_accessions": _sorted_sample(unaccounted),
        "multiply_disposed": len(multiply_disposed),
        "multiply_disposed_accessions": _sorted_sample(multiply_disposed),
        "demoted_from_v1": len(demoted_accessions),
        "demoted_from_v1_accessions": _sorted_sample(demoted_accessions),
        "demote_reason_counts": dict(sorted(demote_reason_counts.items())),
        "demote_reason_union_count": len(demoted_accessions),
        "hold_reason_counts_nonexclusive": dict(sorted(hold_reason_counts.items())),
        "sum_of_hold_reason_counts": sum(hold_reason_counts.values()),
        "hold_reason_union_count": len(hold_accessions),
        "hold_and_demote_overlap_count": len(hold_accessions & set(demoted_accessions)),
        "hold_and_demote_overlap_accessions": _sorted_sample(
            hold_accessions & set(demoted_accessions)
        ),
        "deferred_reason_counts": dict(sorted(defer_reason_counts.items())),
        "deferred_reason_union_count": len(deferred_accessions),
        "deferred_accessions": _sorted_sample(deferred_accessions),
        "high_confidence_dispositions": list(HIGH_CONFIDENCE_DISPOSITIONS),
        "high_confidence_primary_count": len(high_confidence_accessions),
        "records_count": len(normalized),
        "accessions_with_records": len(by_accession),
        "disposition_partition_closes": (
            sum(disposition_partition.values()) + len(unaccounted) == len(universe_set)
        ),
        "disposition_vocabulary": list(DISPOSITIONS),
        "unaccounted_is_a_count_not_a_set": True,
        "note": (
            "primary_disposition is exclusive (exactly one per accession); hold, "
            "demotion and deferred reasons are nonexcusive and are counted as a "
            "union, never as a sum of the per-reason counts"
        ),
    }
    if not summary["disposition_partition_closes"]:
        raise ValueError(
            "disposition partition does not close: %d dispositions + %d unaccounted != %d"
            % (sum(disposition_partition.values()), len(unaccounted), len(universe_set))
        )
    return rows, summary


def summarize_reason_sets(reason_sets: Mapping[str, Iterable[str]]) -> dict[str, Any]:
    """Summarize overlapping reason sets without double-counting proteins.

    The union is the number of *distinct* proteins across all reasons; the sum of
    the per-reason counts is reported beside it so the overlap is never mistaken
    for extra proteins.  Raises ``ValueError`` when no reason set is supplied.
    """
    if not reason_sets:
        raise ValueError("at least one reason set is required")
    sets = {str(reason): _as_set(values) for reason, values in reason_sets.items()}
    if not sets:
        raise ValueError("at least one reason set is required")
    reason_counts = {reason: len(values) for reason, values in sorted(sets.items())}
    union: set[str] = set()
    for values in sets.values():
        union |= values
    membership: Counter[str] = Counter()
    for values in sets.values():
        for accession in values:
            membership[accession] += 1
    overlap = {accession for accession, count in membership.items() if count > 1}
    reasons = sorted(sets)
    pairwise: list[dict[str, Any]] = []
    for index, reason_a in enumerate(reasons):
        for reason_b in reasons[index + 1:]:
            shared = sets[reason_a] & sets[reason_b]
            pairwise.append(
                {
                    "reason_a": reason_a,
                    "reason_b": reason_b,
                    "overlap_count": len(shared),
                    "union_count": len(sets[reason_a] | sets[reason_b]),
                }
            )
    return {
        "reason_names": reasons,
        "reason_counts": reason_counts,
        "sum_of_reason_counts": sum(reason_counts.values()),
        "union_count": len(union),
        "overlap_count": len(overlap),
        "overlap_accessions": sorted(overlap),
        "overlap_is_not_extra_proteins": True,
        "pairwise_overlaps": pairwise,
    }


def overlap_matrix(reason_sets: Mapping[str, Iterable[str]]) -> list[dict[str, Any]]:
    """One row per unordered reason pair (sorted, no self-pairs)."""
    if not reason_sets:
        raise ValueError("at least one reason set is required")
    sets = {str(reason): _as_set(values) for reason, values in reason_sets.items()}
    reasons = sorted(sets)
    rows: list[dict[str, Any]] = []
    for index, reason_a in enumerate(reasons):
        for reason_b in reasons[index + 1:]:
            shared = sets[reason_a] & sets[reason_b]
            union = sets[reason_a] | sets[reason_b]
            rows.append(
                {
                    "reason_a": reason_a,
                    "reason_b": reason_b,
                    "count_a": len(sets[reason_a]),
                    "count_b": len(sets[reason_b]),
                    "sum_of_counts": len(sets[reason_a]) + len(sets[reason_b]),
                    "overlap_count": len(shared),
                    "union_count": len(union),
                    "jaccard": ("%.6f" % (len(shared) / len(union))) if union else "0.000000",
                    "overlap_accessions": ";".join(sorted(shared)),
                    "accessions_only_in_a": ";".join(sorted(sets[reason_a] - sets[reason_b])),
                    "accessions_only_in_b": ";".join(sorted(sets[reason_b] - sets[reason_a])),
                }
            )
    return rows


# ---------------------------------------------------------------------------
# Check machinery
# ---------------------------------------------------------------------------


def assert_check_ok(payload: Mapping[str, Any]) -> Mapping[str, Any]:
    """Raise ``ValueError`` unless ``payload["status"] == "ok"``.

    The message names the check and repeats the explanation, so a failing check
    is loud wherever it is used.
    """
    if payload.get("status") != STATUS_OK:
        raise ValueError(
            "%s: status=%s: %s"
            % (
                payload.get("check", "unnamed check"),
                payload.get("status"),
                payload.get("explanation", "no explanation"),
            )
        )
    return payload


def check_custom_counts(
    observed: Mapping[str, int],
    expected: Mapping[str, int],
) -> dict[str, Any]:
    """Compare pinned expectation keys against observed counts.

    Used for the ``--expect name=value`` pins of the CLI.  A pinned name that has
    no observed count is a mismatch, never a pass.
    """
    mismatches: list[dict[str, Any]] = []
    checked = 0
    for name in sorted(expected):
        want = int(expected[name])
        if name not in observed:
            mismatches.append(
                {"name": name, "expected": want, "observed": None, "difference": None}
            )
            continue
        got = int(observed[name])
        checked += 1
        if got != want:
            mismatches.append(
                {"name": name, "expected": want, "observed": got, "difference": got - want}
            )
    unpinned = sorted(name for name in observed if name not in expected)
    for name in unpinned:
        mismatches.append(
            {
                "name": name,
                "expected": None,
                "observed": int(observed[name]),
                "difference": None,
            }
        )
    return {
        "check": "check_custom_counts",
        "status": STATUS_OK if not mismatches else STATUS_MISMATCH,
        "explanation": (
            "%d pinned count(s) match and every reported count is pinned" % checked
            if not mismatches
            else "%d count(s) do not reconcile: %s"
            % (
                len(mismatches),
                ", ".join(
                    (
                        "%s has no observed value (expected %s)" % (item["name"], item["expected"])
                        if item["observed"] is None
                        else "%s observed %s but expected %s"
                        % (item["name"], item["observed"], item["expected"])
                    )
                    for item in mismatches
                ),
            )
        ),
        "checked_count": checked,
        "mismatches": mismatches,
        "unpinned_names": unpinned,
        "counts": {name: int(value) for name, value in sorted(observed.items())},
        "expectations": {name: int(value) for name, value in sorted(expected.items())},
    }


# ---------------------------------------------------------------------------
# Check 1 -- the hold residual (the "616 unexplained flow records")
# ---------------------------------------------------------------------------


def check_hold_residual(
    *,
    universe: Iterable[str],
    high_confidence: Iterable[str],
    hold: Iterable[str],
    other_dispositions: Mapping[str, Iterable[str]],
    expected_residual_count: int | None = None,
    expression: str = HOLD_RESIDUAL_EXPRESSION,
    other_dispositions_source: str = "caller_supplied",
) -> dict[str, Any]:
    """Explain ``universe - high_confidence - hold`` as an explicit set difference.

    The residual is **never** computed as a subtraction of three aggregates: it is
    the accession set difference, and every residual accession must be named by an
    ``other_dispositions`` bucket.  A non-empty ``residual_unexplained`` set, any
    accession outside the universe, or any overlap between the high-confidence
    and hold layers makes the check a mismatch.

    Expected columns
    ----------------
    ``universe``/``high_confidence``/``hold``: accession sets (``accession``,
    ``protein_id`` or ``genome_id`` column, or a headerless one-per-line list).
    ``other_dispositions``: ``{bucket_name: accession set}``, e.g.
    ``{"excluded_input_quality": {...}, "deferred_structure_review": {...}}``.
    ``other_dispositions_source`` names where those buckets came from
    (``caller_supplied`` by default; the CLI passes
    ``ledger_primary_disposition``, because then the residual is closed by
    construction and the check's value is proving that no accession was left
    unnamed -- a residual accession that the ledger also does not name still
    fails).
    """
    universe_set = _as_set(universe)
    if not universe_set:
        raise ValueError("hold residual check needs a non-empty universe")
    high_confidence_set = _as_set(high_confidence)
    hold_set = _as_set(hold)
    buckets = {str(name): _as_set(values) for name, values in (other_dispositions or {}).items()}

    outside = (high_confidence_set | hold_set | set().union(*buckets.values()) if buckets
               else high_confidence_set | hold_set) - universe_set
    overlap = high_confidence_set & hold_set
    residual = universe_set - high_confidence_set - hold_set
    explained = set().union(*buckets.values()) if buckets else set()
    residual_explained = residual & explained
    residual_unexplained = residual - explained
    bucket_counts = {
        name: len(values & residual) for name, values in sorted(buckets.items())
    }
    bucket_counts["unexplained"] = len(residual_unexplained)

    status = STATUS_OK
    problems: list[str] = []
    if outside:
        status = STATUS_MISMATCH
        problems.append(
            "%d accession(s) are not in the universe: %s"
            % (len(outside), ", ".join(_sorted_sample(outside)))
        )
    if overlap:
        status = STATUS_MISMATCH
        problems.append(
            "%d accession(s) are both high-confidence and held: %s"
            % (len(overlap), ", ".join(_sorted_sample(overlap)))
        )
    if residual_unexplained:
        status = STATUS_MISMATCH
        problems.append(
            "%d residual accession(s) have no disposition: %s"
            % (len(residual_unexplained), ", ".join(_sorted_sample(residual_unexplained)))
        )
    if expected_residual_count is not None and len(residual) != int(expected_residual_count):
        status = STATUS_MISMATCH
        problems.append(
            "residual %d != expected %d" % (len(residual), int(expected_residual_count))
        )
    if sum(bucket_counts.values()) != len(residual):
        status = STATUS_MISMATCH
        problems.append(
            "residual bucket counts sum to %d but the residual is %d"
            % (sum(bucket_counts.values()), len(residual))
        )

    explanation = (
        "residual set difference %s gives %d accession(s) (|universe| %d - "
        "|high_confidence| %d - |hold| %d = %d); %d of them are named by %d disposition "
        "bucket(s) and %d stayed unexplained"
        % (
            expression,
            len(residual),
            len(universe_set),
            len(high_confidence_set),
            len(hold_set),
            len(residual),
            len(residual_explained),
            len(buckets),
            len(residual_unexplained),
        )
    )
    if problems:
        explanation = explanation + "; MISMATCH: " + "; ".join(problems)
    return {
        "check": "check_hold_residual",
        "status": status,
        "explanation": explanation,
        "expression": expression,
        "universe_count": len(universe_set),
        "high_confidence_count": len(high_confidence_set),
        "hold_count": len(hold_set),
        "residual_count": len(residual),
        "residual_accessions": _sorted_sample(residual, 50),
        "residual_explained_count": len(residual_explained),
        "residual_unexplained_count": len(residual_unexplained),
        "residual_unexplained_accessions": _sorted_sample(residual_unexplained, 50),
        "residual_bucket_counts": bucket_counts,
        "residual_bucket_names": sorted(buckets),
        "overlap_high_confidence_hold_count": len(overlap),
        "overlap_high_confidence_hold_accessions": _sorted_sample(overlap),
        "outside_universe_count": len(outside),
        "outside_universe_accessions": _sorted_sample(outside),
        "expected_residual_count": (
            None if expected_residual_count is None else int(expected_residual_count)
        ),
        "other_dispositions_source": other_dispositions_source,
        "problems": problems,
    }


# ---------------------------------------------------------------------------
# Check 2 -- the 52 Cys demotions
# ---------------------------------------------------------------------------


def check_cys_demotion(
    *,
    demoted: Iterable[str],
    hold: Iterable[str],
    high_confidence_after: Iterable[str] | None = None,
    high_confidence_before: Iterable[str] | None = None,
    expected_count: int | None = None,
    hold_reason: str = "localization_conflict",
    required: bool = False,
) -> dict[str, Any]:
    """State whether the demoted records join the hold set or form a separate disposition.

    Every demoted accession must appear exactly once and must be absent from the
    post-demotion high-confidence layer; it is then either in the hold layer
    (``join_hold``) or recorded outside it (``separate_disposition``).  Both sets
    are reported by accession, so the reader never has to reconstruct them, and the
    explanation states in words which of the two the evidence shows.

    A demoted accession that is still high confidence is a mismatch, and a demoted
    accession listed twice is a mismatch.  Demotions held outside the hold layer
    are *not* a defect by themselves: the question the check answers is which layer
    the demotion joined, and a demotion can legitimately be recorded as its own
    disposition rather than inside the hold set.

    Expected columns
    ----------------
    ``demoted``: one accession per row (``accession``/``protein_id``); a repeated
    accession is a mismatch.  ``hold``/``high_confidence_after``/
    ``high_confidence_before``: accession sets.  ``hold_reason`` names the reason
    the demotion is expected to carry in the hold layer (default
    ``localization_conflict``, the reason published by
    ``runs/20260920_phaded_three_gaps_01/results/cys_demotion``).
    """
    demoted_list = [str(value) for value in demoted if str(value)]
    demoted_counts = Counter(demoted_list)
    duplicates = sorted(value for value, count in demoted_counts.items() if count > 1)
    demoted_set = set(demoted_list)
    hold_set = _as_set(hold)
    after_set = None if high_confidence_after is None else _as_set(high_confidence_after)
    before_set = None if high_confidence_before is None else _as_set(high_confidence_before)

    join_hold = demoted_set & hold_set
    still_high_confidence = (
        set() if after_set is None else demoted_set & after_set
    )
    separate = demoted_set - hold_set
    unaccounted = separate

    status = STATUS_OK
    problems: list[str] = []
    if duplicates:
        status = STATUS_MISMATCH
        problems.append(
            "demoted accession(s) listed more than once: %s"
            % ", ".join(_sorted_sample(duplicates))
        )
    if still_high_confidence:
        status = STATUS_MISMATCH
        problems.append(
            "%d demoted accession(s) are still high confidence: %s"
            % (len(still_high_confidence), ", ".join(_sorted_sample(still_high_confidence)))
        )
    if unaccounted is None:  # pragma: no cover - defensive, unaccounted is always a set
        unaccounted = separate
    if expected_count is not None and len(demoted_set) != int(expected_count):
        status = STATUS_MISMATCH
        problems.append(
            "demoted count %d != expected %d" % (len(demoted_set), int(expected_count))
        )
    if required and not demoted_set:
        status = STATUS_MISMATCH
        problems.append("the demotion layer is empty but a demotion is required")

    joins = bool(join_hold) and not separate
    separate_layer = bool(separate)
    if problems:
        explanation = "the %d demoted accession(s) do not reconcile: %s" % (
            len(demoted_set),
            "; ".join(problems),
        )
    elif joins:
        explanation = (
            "the %d demoted accession(s) join the hold set (reason %s); they form no "
            "separate disposition and are absent from the post-demotion high-confidence "
            "layer" % (len(demoted_set), hold_reason)
        )
    else:
        explanation = (
            "the %d demoted accession(s) form a separate disposition: none of them is in "
            "the hold layer; the demotion is recorded outside the hold set and the "
            "accessions are absent from the post-demotion high-confidence layer"
            % len(demoted_set)
        )

    added_to_hold = (
        None
        if before_set is None
        else sorted(((before_set - after_set) if after_set is not None else set()) & hold_set)
    )
    return {
        "check": "check_cys_demotion",
        "status": status,
        "explanation": explanation,
        "demoted_count": len(demoted_set),
        "demoted_accessions": sorted(demoted_set),
        "duplicate_demoted_accessions": duplicates,
        "expected_count": None if expected_count is None else int(expected_count),
        "join_hold": joins,
        "join_hold_count": len(sorted(join_hold)),
        "join_hold_accessions": sorted(join_hold),
        "separate_disposition": separate_layer,
        "separate_disposition_count": len(sorted(separate)),
        "separate_disposition_accessions": sorted(separate),
        "hold_reason": hold_reason,
        "hold_count": len(hold_set),
        "high_confidence_after_count": None if after_set is None else len(after_set),
        "high_confidence_before_count": None if before_set is None else len(before_set),
        "added_to_hold_accessions": added_to_hold,
        "still_high_confidence_count": len(still_high_confidence),
        "still_high_confidence_accessions": _sorted_sample(still_high_confidence),
        "unaccounted_count": len(unaccounted),
        "unaccounted_accessions": sorted(unaccounted),
        "problems": problems,
    }


# ---------------------------------------------------------------------------
# Check 3 -- the Cys funnel union
# ---------------------------------------------------------------------------


def check_cys_funnel_union(
    *,
    phac: Iterable[str],
    iphaz: Iterable[str],
    universe: Iterable[str] | None = None,
    merged_high_confidence: Iterable[str] | None = None,
    expected_union_count: int | None = None,
    expected_sum_of_reason_counts: int | None = None,
) -> dict[str, Any]:
    """Explain why ``|PhaC| + |iPhaZ|`` is larger than their distinct union.

    The union is the number of distinct proteins; ``sum_of_reason_counts`` counts
    the same protein once per reason it matches.  The overlap (``PhaC ∩ iPhaZ``)
    is exactly the difference, and it is reported by accession.  When a candidate
    universe is supplied, the union is additionally split into in-universe and
    out-of-universe parts, because the documented 4,301-protein union is the
    out-of-universe part of the 4,214-hit funnel union.

    Expected columns
    ----------------
    ``phac``/``iphaz``: accession sets (tblout ``target``-derived lists are fine:
    an ``accession``/``protein_id`` column or one identifier per line).
    ``universe``: the frozen candidate universe.
    ``merged_high_confidence``: the merged high-confidence set (40,690 when the
    frozen merge is supplied).  Their intersection with the universe must be
    exactly the in-universe part of the funnel union: a funnel hit that is inside
    the universe must be a merged high-confidence candidate, so any leftover is a
    mismatch.  This is the cross-table reconciliation the 2026-09-21 review
    requested for the 4,301/6,214 relation.
    """
    phac_set = _as_set(phac)
    iphaz_set = _as_set(iphaz)
    if not phac_set and not iphaz_set:
        raise ValueError("cys funnel check needs at least one non-empty hit set")
    summary = summarize_reason_sets({"PhaC": phac_set, "iPhaZ": iphaz_set})
    overlap = phac_set & iphaz_set
    union = phac_set | iphaz_set

    funnel_in_universe: set[str] | None = None
    funnel_outside_universe: set[str] | None = None
    universe_set: set[str] | None = None
    if universe is not None:
        universe_set = _as_set(universe)
        if not universe_set:
            raise ValueError("cys funnel check needs a non-empty universe when one is given")
        funnel_in_universe = union & universe_set
        funnel_outside_universe = union - universe_set

    merged_set: set[str] | None = None
    in_universe_not_merged: set[str] = set()
    if merged_high_confidence is not None:
        merged_set = _as_set(merged_high_confidence)
        if not merged_set:
            raise ValueError(
                "cys funnel check needs a non-empty merged high-confidence set when one is given"
            )
        if funnel_in_universe is None:
            raise ValueError(
                "the merged high-confidence cross-check needs the candidate universe"
            )
        in_universe_not_merged = funnel_in_universe - merged_set

    status = STATUS_OK
    problems: list[str] = []
    if summary["sum_of_reason_counts"] - summary["overlap_count"] != summary["union_count"]:
        status = STATUS_MISMATCH
        problems.append("sum of reason counts minus overlap does not equal the union")
    if expected_union_count is not None and summary["union_count"] != int(expected_union_count):
        status = STATUS_MISMATCH
        problems.append(
            "union %d != expected %d" % (summary["union_count"], int(expected_union_count))
        )
    if (
        expected_sum_of_reason_counts is not None
        and summary["sum_of_reason_counts"] != int(expected_sum_of_reason_counts)
    ):
        status = STATUS_MISMATCH
        problems.append(
            "sum of reason counts %d != expected %d"
            % (summary["sum_of_reason_counts"], int(expected_sum_of_reason_counts))
        )
    if in_universe_not_merged:
        status = STATUS_MISMATCH
        problems.append(
            "%d in-universe funnel hit(s) are not merged high-confidence: %s"
            % (
                len(in_universe_not_merged),
                ", ".join(_sorted_sample(in_universe_not_merged)),
            )
        )

    explanation = (
        "PhaC %d + iPhaZ %d = %d reported rows, but %d protein(s) are in both reasons, so "
        "the distinct union is %d (= %d - %d), not the sum"
        % (
            summary["reason_counts"]["PhaC"],
            summary["reason_counts"]["iPhaZ"],
            summary["sum_of_reason_counts"],
            summary["overlap_count"],
            summary["union_count"],
            summary["sum_of_reason_counts"],
            summary["overlap_count"],
        )
    )
    if funnel_outside_universe is not None:
        explanation += (
            "; %d of the union is inside the candidate universe and %d is outside it"
            % (len(funnel_in_universe or ()), len(funnel_outside_universe))
        )
    if merged_set is not None:
        explanation += (
            "; all %d in-universe funnel hit(s) are in the %d merged high-confidence set"
            % (len(funnel_in_universe or ()), len(merged_set))
        )
    if problems:
        explanation = explanation + "; MISMATCH: " + "; ".join(problems)

    return {
        "check": "check_cys_funnel_union",
        "status": status,
        "explanation": explanation,
        "reason_counts": summary["reason_counts"],
        "sum_of_reason_counts": summary["sum_of_reason_counts"],
        "union_count": summary["union_count"],
        "overlap_count": summary["overlap_count"],
        "overlap_accessions": _sorted_sample(overlap, 50),
        "pairwise_overlaps": summary["pairwise_overlaps"],
        "phac_count": len(phac_set),
        "iphaz_count": len(iphaz_set),
        "universe_count": None if universe_set is None else len(universe_set),
        "funnel_in_universe_count": (
            None if funnel_in_universe is None else len(funnel_in_universe)
        ),
        "funnel_in_universe_accessions": (
            None if funnel_in_universe is None else _sorted_sample(funnel_in_universe, 50)
        ),
        "funnel_outside_universe_count": (
            None if funnel_outside_universe is None else len(funnel_outside_universe)
        ),
        "funnel_outside_universe_accessions": (
            None
            if funnel_outside_universe is None
            else _sorted_sample(funnel_outside_universe, 50)
        ),
        "merged_high_confidence_count": (
            None if merged_set is None else len(merged_set)
        ),
        "in_universe_not_merged_count": (
            None if merged_set is None else len(in_universe_not_merged)
        ),
        "in_universe_not_merged_accessions": (
            None
            if merged_set is None
            else _sorted_sample(in_universe_not_merged, 50)
        ),
        "expected_union_count": (
            None if expected_union_count is None else int(expected_union_count)
        ),
        "expected_sum_of_reason_counts": (
            None
            if expected_sum_of_reason_counts is None
            else int(expected_sum_of_reason_counts)
        ),
        "problems": problems,
    }


# ---------------------------------------------------------------------------
# Check 4 -- the gRodon manifest delta
# ---------------------------------------------------------------------------


def check_grodon_manifest_delta(
    *,
    eligible_positives: Iterable[str],
    manifest_positives: Iterable[str],
    exclusion_reasons: Mapping[str, str],
    prediction_failures: Mapping[str, str],
    expected_eligible_count: int | None = None,
    expected_manifest_count: int | None = None,
    expected_difference: int | None = None,
) -> dict[str, Any]:
    """Explain the eligible-positive to manifest-positive difference by reason bucket.

    Every eligible positive is either included in the manifest or carries exactly
    one reason it is not: a selection exclusion (no same-genus control, missing
    genome FASTA, ...) or a prediction failure (``too_few_ribosomal_hits``,
    ``hmmsearch_failed``, ...).  The per-reason bucket counts must sum to the
    difference, no accession may be explained twice, and a manifest positive that
    is not an eligible positive is a mismatch.

    Expected columns
    ----------------
    ``exclusion_reasons``: ``{genome_id: reason}`` for eligible positives that the
    manifest selection dropped.  ``prediction_failures``:
    ``{genome_id: error}`` for genomes that were selected but whose prediction
    failed.  Both keys are genome ids; both value columns are free strings that
    become the bucket names.
    """
    eligible = _as_set(eligible_positives)
    manifest = _as_set(manifest_positives)
    if not eligible:
        raise ValueError("gRodon delta check needs a non-empty eligible-positive set")
    exclusions = {str(key): str(value) for key, value in (exclusion_reasons or {}).items()}
    failures = {str(key): str(value) for key, value in (prediction_failures or {}).items()}

    difference = len(eligible) - len(manifest)
    excluded = eligible - manifest
    explained_by = Counter[str]()
    multiply_explained: list[str] = []
    for accession in excluded:
        reasons = [
            reason
            for reason in (exclusions.get(accession), failures.get(accession))
            if reason
        ]
        if len(reasons) > 1:
            multiply_explained.append(accession)
        if reasons:
            explained_by[reasons[0]] += 1
    explained = {
        accession
        for accession in excluded
        if exclusions.get(accession) or failures.get(accession)
    }
    unexplained = excluded - explained
    manifest_not_eligible = manifest - eligible
    unexpected_keys = sorted(
        (set(exclusions) | set(failures)) - eligible
    )

    status = STATUS_OK
    problems: list[str] = []
    if manifest_not_eligible:
        status = STATUS_MISMATCH
        problems.append(
            "%d manifest positive(s) are not eligible positives: %s"
            % (len(manifest_not_eligible), ", ".join(_sorted_sample(manifest_not_eligible)))
        )
    if unexplained:
        status = STATUS_MISMATCH
        problems.append(
            "%d eligible positive(s) have no disposition: %s"
            % (len(unexplained), ", ".join(_sorted_sample(unexplained)))
        )
    if multiply_explained:
        status = STATUS_MISMATCH
        problems.append(
            "%d accession(s) carry more than one disposition: %s"
            % (len(multiply_explained), ", ".join(_sorted_sample(multiply_explained)))
        )
    if sum(explained_by.values()) != difference:
        status = STATUS_MISMATCH
        problems.append(
            "reason buckets sum to %d but the difference is %d"
            % (sum(explained_by.values()), difference)
        )
    for name, want in (
        ("eligible_positive_count", expected_eligible_count),
        ("manifest_positive_count", expected_manifest_count),
        ("difference", expected_difference),
    ):
        if want is None:
            continue
        got = {
            "eligible_positive_count": len(eligible),
            "manifest_positive_count": len(manifest),
            "difference": difference,
        }[name]
        if got != int(want):
            status = STATUS_MISMATCH
            problems.append("%s %d != expected %d" % (name, got, int(want)))

    explanation = (
        "%d eligible positive(s) - %d manifest positive(s) = %d record(s), explained by "
        "%s"
        % (
            len(eligible),
            len(manifest),
            difference,
            ", ".join("%s %d" % (reason, count) for reason, count in sorted(explained_by.items()))
            or "no bucket",
        )
    )
    if problems:
        explanation = explanation + "; MISMATCH: " + "; ".join(problems)

    count_expectation_failures = [
        name
        for name, want, got in (
            ("eligible_positive_count", expected_eligible_count, len(eligible)),
            ("manifest_positive_count", expected_manifest_count, len(manifest)),
            ("difference", expected_difference, difference),
        )
        if want is not None and got != int(want)
    ]
    return {
        "check": "check_grodon_manifest_delta",
        "status": status,
        "explanation": explanation,
        "eligible_positive_count": len(eligible),
        "manifest_positive_count": len(manifest),
        "difference": difference,
        "reason_bucket_counts": dict(sorted(explained_by.items())),
        "delta_explained_count": len(explained),
        "delta_unexplained_count": len(unexplained),
        "delta_unexplained_accessions": _sorted_sample(unexplained, 50),
        "manifest_not_eligible_count": len(manifest_not_eligible),
        "manifest_not_eligible_accessions": _sorted_sample(manifest_not_eligible, 50),
        "multiply_explained_count": len(multiply_explained),
        "multiply_explained_accessions": _sorted_sample(multiply_explained),
        "no_genus_control_count": explained_by.get("no_genus_control", 0),
        "missing_fasta_count": explained_by.get("missing_fasta", 0),
        "failed_prediction_count": sum(
            count
            for reason, count in explained_by.items()
            if reason in set(failures.values())
        ),
        "exclusion_reason_keys_outside_eligible": unexpected_keys,
        "expected_eligible_count": (
            None if expected_eligible_count is None else int(expected_eligible_count)
        ),
        "expected_manifest_count": (
            None if expected_manifest_count is None else int(expected_manifest_count)
        ),
        "expected_difference": (
            None if expected_difference is None else int(expected_difference)
        ),
        "count_expectation_failures": count_expectation_failures,
        "problems": problems,
    }


# ---------------------------------------------------------------------------
# Check 5 -- motif coverage versus motif resolvability
# ---------------------------------------------------------------------------


def check_motif_coverage_vs_resolvability(
    *,
    covered: Iterable[str],
    resolvable: Iterable[str],
    criteria_resolved: Mapping[str, Iterable[str]] | None = None,
    covered_but_unresolvable_is_expected: bool = True,
) -> dict[str, Any]:
    """Separate "the motif run reached this protein" from "the motif is resolvable here".

    These are two different numbers and the difference is a finding, not a
    rounding: a protein the run reached but for which no criterion is assignable
    (``covered_but_unresolvable``) carries no motif information at all, and a
    protein with an assignable criterion that the run never reached
    (``resolvable_but_uncovered``) is a coverage defect.  Both sets are listed
    separately.

    ``resolvable_but_uncovered`` is **always** a mismatch: a resolvable motif row
    outside the run coverage means the run missed a protein it should have
    processed.  ``covered_but_unresolvable`` is a mismatch unless the caller states
    that the state is expected for this layer -- the default is ``True`` because
    the frozen motif layer genuinely carries criteria (``linker``, ``sbd``) that
    are unassignable for every candidate, so covering a protein is not the same as
    resolving it.  Pass ``False`` to make even that state a hard finding.

    Expected columns
    ----------------
    ``covered``/``resolvable``: accession sets.  ``criteria_resolved``:
    ``{criterion: accession set with that criterion assignable}``.
    """
    covered_set = _as_set(covered)
    resolvable_set = _as_set(resolvable)
    if not covered_set:
        raise ValueError("motif check needs a non-empty coverage layer")
    criteria = {
        str(name): _as_set(values) for name, values in (criteria_resolved or {}).items()
    }

    covered_but_unresolvable = covered_set - resolvable_set
    resolvable_but_uncovered = resolvable_set - covered_set
    criteria_payload: dict[str, Any] = {}
    criteria_uncovered: set[str] = set()
    for name, values in sorted(criteria.items()):
        uncovered = values - covered_set
        criteria_uncovered |= uncovered
        criteria_payload[name] = {
            "covered_count": len(covered_set),
            "resolvable_count": len(values),
            "unresolvable_count": len(covered_set - values),
            "uncovered_count": len(uncovered),
            "resolvable_but_uncovered_accessions": _sorted_sample(uncovered),
        }

    status = STATUS_OK
    problems: list[str] = []
    if resolvable_but_uncovered:
        status = STATUS_MISMATCH
        problems.append(
            "%d resolvable accession(s) were never covered: %s"
            % (len(resolvable_but_uncovered), ", ".join(_sorted_sample(resolvable_but_uncovered)))
        )
    if criteria_uncovered:
        status = STATUS_MISMATCH
        problems.append(
            "%d criterion row(s) are resolvable outside the coverage layer: %s"
            % (len(criteria_uncovered), ", ".join(_sorted_sample(criteria_uncovered)))
        )
    if covered_but_unresolvable and not covered_but_unresolvable_is_expected:
        status = STATUS_MISMATCH
        problems.append(
            "%d covered accession(s) have no assignable criterion: %s"
            % (
                len(covered_but_unresolvable),
                ", ".join(_sorted_sample(covered_but_unresolvable)),
            )
        )

    explanation = (
        "the motif run covered %d protein(s) and resolved %d of them, so the two numbers "
        "differ by %d covered-but-unresolvable protein(s); %d protein(s) are resolvable but "
        "uncovered"
        % (
            len(covered_set),
            len(resolvable_set),
            len(covered_but_unresolvable),
            len(resolvable_but_uncovered),
        )
    )
    if problems:
        explanation = explanation + "; MISMATCH: " + "; ".join(problems)

    return {
        "check": "check_motif_coverage_vs_resolvability",
        "status": status,
        "explanation": explanation,
        "covered_count": len(covered_set),
        "resolvable_count": len(resolvable_set),
        "covered_but_unresolvable_count": len(covered_but_unresolvable),
        "covered_but_unresolvable_accessions": _sorted_sample(covered_but_unresolvable, 50),
        "resolvable_but_uncovered_count": len(resolvable_but_uncovered),
        "resolvable_but_uncovered_accessions": _sorted_sample(resolvable_but_uncovered, 50),
        "covered_but_unresolvable_is_expected": bool(covered_but_unresolvable_is_expected),
        "criteria": criteria_payload,
        "problem_count": len(problems),
        "problems": problems,
    }


# ---------------------------------------------------------------------------
# Check 6 -- the profile inventory
# ---------------------------------------------------------------------------


def check_profile_inventory(
    *,
    authority_rows: Sequence[Mapping[str, Any]],
    family_rows: Sequence[Mapping[str, Any]],
    profile_rows: Sequence[Mapping[str, Any]],
    expected_superfamilies: int | None = None,
    expected_families: int | None = None,
    expected_trained: int | None = None,
    expected_reference_only: int | None = None,
) -> dict[str, Any]:
    """Count superfamilies, families, trained HMMs and reference-only profiles from the tables.

    Nothing is hard-coded: every count is recomputed from the supplied rows, and
    the expectations default to ``None`` (recompute only).  The CLI pins the
    documented values (8 / 38 / 9 / 37) so the real inventory is checked against
    the review, while a caller with other tables is not forced to fake them.

    The totals must reconcile -- superfamily definitions + family definitions must
    equal the profile count, and trained + reference-only must equal it too -- and
    every family definition must name a superfamily that the ledger carries.

    Expected columns
    ----------------
    ``authority_rows``: ``layer_id``, ``layer_kind`` in
    ``{"superfamily", "family"}`` (optional ``gate_profile_id``, which must be
    covered by the profile manifest when present).
    ``family_rows``: ``phaded_superfamily``, ``phaded_family_id``.
    ``profile_rows``: ``profile_id``, ``model_status`` in
    ``{"trained", "reference_only"}``, optional matching ``profile_kind``.
    """
    _require_columns(authority_rows, AUTHORITY_REQUIRED_COLUMNS, "classification authority")
    if not family_rows:
        raise ValueError("family definition table is empty")
    _require_columns(
        family_rows, FAMILY_DEFINITION_REQUIRED_COLUMNS, "family definition table"
    )
    if not profile_rows:
        raise ValueError("profile manifest is empty")
    _require_columns(profile_rows, PROFILE_REQUIRED_COLUMNS, "profile manifest")

    layer_kind_counts = Counter(
        str(row.get("layer_kind", "")).strip() for row in authority_rows
    )
    layer_kind_counts.pop("", None)
    superfamily_count = layer_kind_counts.get("superfamily", 0)
    family_count_from_authority = layer_kind_counts.get("family", 0)
    definition_counts = dict(sorted(layer_kind_counts.items()))
    definition_count = sum(definition_counts.values())

    duplicate_layer_ids = sorted(
        value
        for value, count in Counter(
            str(row.get("layer_id", "")).strip() for row in authority_rows
        ).items()
        if value and count > 1
    )
    empty_layer_ids = [
        index
        for index, row in enumerate(authority_rows)
        if not str(row.get("layer_id", "")).strip()
    ]
    unknown_layer_kinds = sorted(
        kind for kind in layer_kind_counts if kind not in {"superfamily", "family"}
    )

    superfamily_ids = {
        str(row.get("layer_id", "")).strip()
        for row in authority_rows
        if str(row.get("layer_kind", "")).strip() == "superfamily"
    }
    family_ids = {
        str(row.get("phaded_family_id", "")).strip() for row in family_rows
    }
    family_rows_bound = [
        row
        for row in family_rows
        if str(row.get(SUPERFAMILY_COLUMN, "")).strip() in superfamily_ids
    ]
    families_with_unknown_superfamily = sorted(
        str(row.get(FAMILY_ID_COLUMN, "")).strip()
        for row in family_rows
        if str(row.get(SUPERFAMILY_COLUMN, "")).strip() not in superfamily_ids
    )
    duplicate_family_ids = sorted(
        value
        for value, count in Counter(
            str(row.get(FAMILY_ID_COLUMN, "")).strip() for row in family_rows
        ).items()
        if value and count > 1
    )
    covered_superfamilies = sorted(
        {
            str(row.get(SUPERFAMILY_COLUMN, "")).strip()
            for row in family_rows
            if str(row.get(SUPERFAMILY_COLUMN, "")).strip() in superfamily_ids
        }
    )
    uncovered_superfamilies = sorted(superfamily_ids - set(covered_superfamilies))

    model_status_counts = Counter(
        str(row.get("model_status", "")).strip() for row in profile_rows
    )
    model_status_counts.pop("", None)
    trained_count = model_status_counts.get(MODEL_STATUS_TRAINED, 0)
    reference_only_count = model_status_counts.get(MODEL_STATUS_REFERENCE_ONLY, 0)
    profile_count = len(profile_rows)
    unknown_model_statuses = sorted(
        status for status in model_status_counts if status not in MODEL_STATUSES
    )
    duplicate_profile_ids = sorted(
        value
        for value, count in Counter(
            str(row.get("profile_id", "")).strip() for row in profile_rows
        ).items()
        if value and count > 1
    )
    profile_ids = {str(row.get("profile_id", "")).strip() for row in profile_rows}
    gate_profile_ids = {
        str(row.get("gate_profile_id", "")).strip()
        for row in authority_rows
        if str(row.get("gate_profile_id", "")).strip()
    }
    gate_profiles_missing = sorted(gate_profile_ids - profile_ids)
    profile_kind_mismatches = sorted(
        str(row.get("profile_id", "")).strip()
        for row in profile_rows
        if PROFILE_KIND_COLUMN in row
        and str(row.get(PROFILE_KIND_COLUMN, "")).strip()
        and str(row.get(PROFILE_KIND_COLUMN, "")).strip()
        not in {"superfamily", "family"}
    )
    profile_definition_mismatch = abs(profile_count - definition_count)

    status = STATUS_OK
    problems: list[str] = []
    for name, want, got in (
        ("superfamily_count", expected_superfamilies, superfamily_count),
        ("family_count", expected_families, family_count_from_authority),
        ("trained_count", expected_trained, trained_count),
        ("reference_only_count", expected_reference_only, reference_only_count),
    ):
        if want is not None and got != int(want):
            status = STATUS_MISMATCH
            problems.append("%s %d != expected %d" % (name, got, int(want)))
    if family_count_from_authority != len(family_rows):
        status = STATUS_MISMATCH
        problems.append(
            "authority family rows %d != family definition rows %d"
            % (family_count_from_authority, len(family_rows))
        )
    if profile_definition_mismatch:
        status = STATUS_MISMATCH
        problems.append(
            "profile count %d != definition count %d" % (profile_count, definition_count)
        )
    if trained_count + reference_only_count != profile_count:
        status = STATUS_MISMATCH
        problems.append(
            "trained %d + reference_only %d != profile count %d"
            % (trained_count, reference_only_count, profile_count)
        )
    if unknown_model_statuses:
        status = STATUS_MISMATCH
        problems.append(
            "unknown model_status value(s): %s" % ", ".join(unknown_model_statuses)
        )
    if families_with_unknown_superfamily:
        status = STATUS_MISMATCH
        problems.append(
            "%d family definition(s) name an unknown superfamily: %s"
            % (
                len(families_with_unknown_superfamily),
                ", ".join(_sorted_sample(families_with_unknown_superfamily)),
            )
        )
    if duplicate_layer_ids or empty_layer_ids:
        status = STATUS_MISMATCH
        problems.append(
            "layer_id is not a unique non-empty key (%d duplicate, %d empty)"
            % (len(duplicate_layer_ids), len(empty_layer_ids))
        )
    if duplicate_family_ids:
        status = STATUS_MISMATCH
        problems.append(
            "duplicate phaded_family_id: %s"
            % ", ".join(_sorted_sample(duplicate_family_ids))
        )
    if duplicate_profile_ids:
        status = STATUS_MISMATCH
        problems.append(
            "duplicate profile_id: %s"
            % ", ".join(_sorted_sample(duplicate_profile_ids))
        )
    if unknown_layer_kinds:
        status = STATUS_MISMATCH
        problems.append("unknown layer_kind value(s): %s" % ", ".join(unknown_layer_kinds))
    if gate_profiles_missing:
        status = STATUS_MISMATCH
        problems.append(
            "%d authority gate_profile_id(s) have no profile row: %s"
            % (len(gate_profiles_missing), ", ".join(_sorted_sample(gate_profiles_missing)))
        )
    if profile_kind_mismatches:
        status = STATUS_MISMATCH
        problems.append(
            "profile_kind outside {superfamily, family}: %s"
            % ", ".join(_sorted_sample(profile_kind_mismatches))
        )
    if uncovered_superfamilies:
        status = STATUS_MISMATCH
        problems.append(
            "%d superfamily(ies) have no family definition: %s"
            % (len(uncovered_superfamilies), ", ".join(_sorted_sample(uncovered_superfamilies)))
        )

    explanation = (
        "%d superfamily definition(s) + %d family definition(s) = %d; the profile manifest "
        "carries %d profile(s) = %d trained + %d reference-only; %d family definition(s) are "
        "bound to a known superfamily"
        % (
            superfamily_count,
            family_count_from_authority,
            definition_count,
            profile_count,
            trained_count,
            reference_only_count,
            len(family_rows_bound),
        )
    )
    if problems:
        explanation = explanation + "; MISMATCH: " + "; ".join(problems)

    return {
        "check": "check_profile_inventory",
        "status": status,
        "explanation": explanation,
        "superfamily_count": superfamily_count,
        "family_count": family_count_from_authority,
        "family_definition_row_count": len(family_rows),
        "definition_count": definition_count,
        "definition_counts": definition_counts,
        "profile_count": profile_count,
        "trained_count": trained_count,
        "reference_only_count": reference_only_count,
        "trained_plus_reference_only": trained_count + reference_only_count,
        "profile_definition_mismatch": profile_definition_mismatch,
        "model_status_counts": dict(sorted(model_status_counts.items())),
        "unknown_model_statuses": unknown_model_statuses,
        "families_bound_to_a_superfamily_count": len(family_rows_bound),
        "families_with_unknown_superfamily": _sorted_sample(families_with_unknown_superfamily),
        "superfamily_ids": sorted(superfamily_ids),
        "superfamilies_with_a_family_definition": covered_superfamilies,
        "superfamilies_without_a_family_definition": _sorted_sample(uncovered_superfamilies),
        "duplicate_layer_ids": _sorted_sample(duplicate_layer_ids),
        "duplicate_family_ids": _sorted_sample(duplicate_family_ids),
        "duplicate_profile_ids": _sorted_sample(duplicate_profile_ids),
        "empty_layer_id_rows": len(empty_layer_ids),
        "unknown_layer_kinds": unknown_layer_kinds,
        "gate_profiles_missing": _sorted_sample(gate_profiles_missing),
        "profile_kind_mismatches": _sorted_sample(profile_kind_mismatches),
        "expected": {
            "superfamilies": expected_superfamilies,
            "families": expected_families,
            "trained": expected_trained,
            "reference_only": expected_reference_only,
        },
        "problems": problems,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


CHECK_EXPECTATION_KEYS: dict[str, tuple[str, ...]] = {
    "check_hold_residual": ("hold_residual_count", "hold_residual_unexplained_count"),
    "check_cys_demotion": ("cys_demoted_count",),
    "check_cys_funnel_union": ("cys_funnel_union_count", "cys_funnel_sum_of_reason_counts"),
    "check_grodon_manifest_delta": (
        "grodon_eligible_positive_count",
        "grodon_manifest_positive_count",
        "grodon_delta",
    ),
    "check_motif_coverage_vs_resolvability": (
        "motif_covered_count",
        "motif_resolvable_count",
    ),
    "check_profile_inventory": (
        "profile_superfamily_count",
        "profile_family_count",
        "profile_trained_count",
        "profile_reference_only_count",
    ),
}


def build_arg_namespace(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Close the PhaDED candidate flow ledger and emit the machine-readable "
            "reconciliation reports. Writes only inside --output-dir."
        ),
        epilog=(
            "checks and their inputs: check_hold_residual (--universe --records "
            "--hold-reasons --high-confidence) | check_cys_demotion (--cys-demoted "
            "--hold-reasons --high-confidence --v1-high-confidence) | "
            "check_cys_funnel_union (--phac-hits --iphaz-hits --high-confidence) | "
            "check_grodon_manifest_delta (--grodon-eligible --grodon-manifest "
            "--grodon-dispositions) | check_motif_coverage_vs_resolvability "
            "(--motif-table) | check_profile_inventory (--authority "
            "--family-definitions --profile-manifest). A missing input skips its "
            "check and is reported in checks_skipped; a supplied input that cannot "
            "be read is an error."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--universe", type=Path, required=True,
                       help="candidate universe accessions (TSV with an accession column or one per line)")
    parser.add_argument("--records", type=Path, required=True,
                       help="candidate record ledger: accession, primary_disposition, "
                            "and optional hold_reasons/demote_reasons/defer_reasons/deferred columns")
    parser.add_argument("--output-dir", type=Path, required=True,
                       help="must not exist or must be empty; nothing is written outside it")
    parser.add_argument("--hold-reasons", type=Path, default=None,
                       help="hold layer table: accession + hold_reason (checks 1 and 2)")
    parser.add_argument("--v1-high-confidence", type=Path, default=None,
                       help="v1 high-confidence accessions, used only for the v1/v2 impact summary")
    parser.add_argument("--high-confidence", type=Path, default=None,
                       help="post-demotion high-confidence accessions (checks 1 and 2)")
    parser.add_argument("--cys-demoted", type=Path, default=None,
                       help="demoted accession table: accession (+ optional demote_reason) (check 2)")
    parser.add_argument("--phac-hits", type=Path, default=None,
                       help="PhaC funnel hit accessions (check 3)")
    parser.add_argument("--iphaz-hits", type=Path, default=None,
                       help="iPhaZ funnel hit accessions (check 3)")
    parser.add_argument("--grodon-eligible", type=Path, default=None,
                       help="eligible-positive genome ids (check 4)")
    parser.add_argument("--grodon-manifest", type=Path, default=None,
                       help="manifest table: genome_id + phaZ_status (check 4)")
    parser.add_argument("--grodon-dispositions", type=Path, default=None,
                       help="genome_id + reason for records that are not manifest positives (check 4)")
    parser.add_argument("--motif-table", type=Path, default=None,
                       help="motif completion table with <criterion>_state and "
                            "<criterion>_assignability columns (check 5)")
    parser.add_argument("--authority", type=Path, default=None,
                       help="phaded_classification_authority.tsv: layer_id + layer_kind (check 6)")
    parser.add_argument("--family-definitions", type=Path, default=None,
                       help="phaded_family_definitions.tsv: phaded_superfamily + phaded_family_id (check 6)")
    parser.add_argument("--profile-manifest", type=Path, default=None,
                       help="profile_manifest.tsv: profile_id + model_status (check 6)")
    parser.add_argument("--expect", action="append", default=[], metavar="NAME=VALUE",
                       help="pin a reported count; a mismatch makes the CLI exit non-zero")
    return parser.parse_args(argv)


def parse_expectations(values: Iterable[str]) -> dict[str, int]:
    """Parse repeated ``NAME=VALUE`` pins into ``{name: int}``."""
    pins: dict[str, int] = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--expect needs NAME=VALUE, got %r" % value)
        name, raw = value.split("=", 1)
        name = name.strip()
        if not name:
            raise ValueError("--expect needs a non-empty name, got %r" % value)
        try:
            pins[name] = int(raw.strip())
        except ValueError:
            raise ValueError("--expect %s needs an integer, got %r" % (name, raw)) from None
    return pins


def expectation_key_map(pins: Mapping[str, int]) -> dict[str, dict[str, int]]:
    """Route pinned expectation names to the check they belong to."""
    routed: dict[str, dict[str, int]] = {}
    for name, value in pins.items():
        for check_name, keys in CHECK_EXPECTATION_KEYS.items():
            if name in keys:
                routed.setdefault(check_name, {})[name] = value
                break
        else:
            known = sorted(
                key for keys in CHECK_EXPECTATION_KEYS.values() for key in keys
            )
            raise ValueError(
                "unknown expectation %r; known names: %s" % (name, ", ".join(known))
            )
    return routed


def observed_check_counts(payload: Mapping[str, Any]) -> dict[str, int]:
    """The counts a check exposes to ``--expect``."""
    mapping = {
        "check_hold_residual": (
            ("hold_residual_count", "residual_count"),
            ("hold_residual_unexplained_count", "residual_unexplained_count"),
        ),
        "check_cys_demotion": (("cys_demoted_count", "demoted_count"),),
        "check_cys_funnel_union": (
            ("cys_funnel_union_count", "union_count"),
            ("cys_funnel_sum_of_reason_counts", "sum_of_reason_counts"),
        ),
        "check_grodon_manifest_delta": (
            ("grodon_eligible_positive_count", "eligible_positive_count"),
            ("grodon_manifest_positive_count", "manifest_positive_count"),
            ("grodon_delta", "difference"),
        ),
        "check_motif_coverage_vs_resolvability": (
            ("motif_covered_count", "covered_count"),
            ("motif_resolvable_count", "resolvable_count"),
        ),
        "check_profile_inventory": (
            ("profile_superfamily_count", "superfamily_count"),
            ("profile_family_count", "family_count"),
            ("profile_trained_count", "trained_count"),
            ("profile_reference_only_count", "reference_only_count"),
        ),
    }
    out: dict[str, int] = {}
    for name, key in mapping.get(str(payload.get("check", "")), ()):
        if key in payload and payload[key] is not None:
            out[name] = int(payload[key])
    return out


def read_genome_status_table(
    path: Path, what: str = "gRodon manifest"
) -> tuple[list[str], dict[str, list[str]], dict[str, str]]:
    """Read a genome table with a status column.

    Returns ``(genome_ids, {status: [genome_id, ...]}, {genome_id: status})``.  The
    accession column is ``genome_id``/``accession``/``protein_id``; the status
    column is ``phaZ_status`` when present, otherwise the caller's ``what`` label
    must say which status column is used -- this reader requires ``phaZ_status``.
    """
    header, rows = _read_table(path, what)
    if not header:
        raise ValueError("%s needs a header row" % what)
    accession_column = _pick_accession_column(header, what)
    if "phaZ_status" not in header:
        raise ValueError("%s is missing the phaZ_status column" % what)
    by_status: dict[str, list[str]] = {}
    status_by_id: dict[str, str] = {}
    for row in rows:
        genome_id = (row.get(accession_column) or "").strip()
        if not genome_id:
            raise ValueError("%s has a row without a genome id" % what)
        status = (row.get("phaZ_status") or "").strip()
        by_status.setdefault(status, []).append(genome_id)
        status_by_id[genome_id] = status
    return sorted(status_by_id), {
        status: sorted(set(ids)) for status, ids in sorted(by_status.items())
    }, status_by_id


def read_disposition_reasons(path: Path, what: str = "disposition reasons") -> dict[str, str]:
    """Read ``genome_id`` + reason (``reason``/``status``/``error``/``disposition``)."""
    header, rows = _read_table(path, what)
    if not header:
        raise ValueError("%s needs a header row" % what)
    accession_column = _pick_accession_column(header, what)
    reason_column = None
    for candidate in ("reason", "disposition", "status", "error"):
        if candidate in header:
            reason_column = candidate
            break
    if reason_column is None:
        raise ValueError(
            "%s needs a reason column (reason/disposition/status/error); found: %s"
            % (what, ", ".join(header))
        )
    out: dict[str, str] = {}
    for row in rows:
        genome_id = (row.get(accession_column) or "").strip()
        reason = (row.get(reason_column) or "").strip()
        if not genome_id:
            raise ValueError("%s has a row without a genome id" % what)
        if not reason:
            raise ValueError("%s has an empty %s for %s" % (what, reason_column, genome_id))
        out[genome_id] = reason
    return out


def run(args: argparse.Namespace) -> dict[str, Any]:
    """Execute the ledger and every check whose inputs were supplied."""
    universe_path = _resolve_input(args.universe, "candidate universe")
    records_path = _resolve_input(args.records, "candidate records")
    output_dir = ensure_output_dir(args.output_dir)

    universe = read_accession_list(universe_path, "candidate universe")
    records = read_records_table(records_path)
    universe_set = set(universe)

    rows, summary = reconcile(universe=universe, records=records)
    write_tsv(output_dir / "candidate_flow.tsv", FLOW_FIELDS, rows)
    write_json(output_dir / "candidate_flow_summary.json", summary)

    pins = parse_expectations(args.expect)
    routed = expectation_key_map(pins)

    checks: dict[str, Any] = {}
    skipped: dict[str, str] = {}

    # check 1: hold residual
    if args.hold_reasons is None:
        skipped["check_hold_residual"] = "no --hold-reasons table supplied"
    else:
        hold_reasons_by_reason = read_reason_table(args.hold_reasons, "hold_reason")
        hold_set = {accession for ids in hold_reasons_by_reason.values() for accession in ids}
        if args.high_confidence is None:
            high_confidence_set = {
                row["accession"]
                for row in rows
                if row["primary_disposition"] in HIGH_CONFIDENCE_DISPOSITIONS
            }
        else:
            high_confidence_set = set(
                read_accession_list(args.high_confidence, "post-demotion high confidence")
            )
        other: dict[str, list[str]] = {}
        for disposition in DISPOSITIONS:
            members = [
                row["accession"]
                for row in rows
                if row["primary_disposition"] == disposition
                and row["accession"] not in high_confidence_set
                and row["accession"] not in hold_set
            ]
            if members:
                other[disposition] = members
        hold_residual_pins = routed.get("check_hold_residual", {})
        checks["check_hold_residual"] = check_hold_residual(
            universe=universe_set,
            high_confidence=high_confidence_set,
            hold=hold_set,
            other_dispositions=other,
            expected_residual_count=hold_residual_pins.get("hold_residual_count"),
            other_dispositions_source="ledger_primary_disposition",
        )

    # check 2: the Cys demotions
    if args.cys_demoted is None:
        skipped["check_cys_demotion"] = "no --cys-demoted table supplied"
    else:
        demoted = read_accession_list(args.cys_demoted, "demoted accessions")
        hold_set = (
            {
                accession
                for ids in read_reason_table(args.hold_reasons, "hold_reason").values()
                for accession in ids
            }
            if args.hold_reasons is not None
            else set()
        )
        after = (
            set(read_accession_list(args.high_confidence, "post-demotion high confidence"))
            if args.high_confidence is not None
            else None
        )
        before = (
            set(read_accession_list(args.v1_high_confidence, "v1 high confidence"))
            if args.v1_high_confidence is not None
            else None
        )
        checks["check_cys_demotion"] = check_cys_demotion(
            demoted=demoted,
            hold=hold_set,
            high_confidence_after=after,
            high_confidence_before=before,
            expected_count=routed.get("check_cys_demotion", {}).get("cys_demoted_count"),
            required=True,
        )

    # check 3: the Cys funnel union
    if args.phac_hits is None or args.iphaz_hits is None:
        skipped["check_cys_funnel_union"] = (
            "needs both --phac-hits and --iphaz-hits accession lists"
        )
    else:
        funnel_pins = routed.get("check_cys_funnel_union", {})
        merged = (
            set(read_accession_list(args.high_confidence, "post-demotion high confidence"))
            if args.high_confidence is not None
            else None
        )
        checks["check_cys_funnel_union"] = check_cys_funnel_union(
            phac=read_accession_list(args.phac_hits, "PhaC funnel hits"),
            iphaz=read_accession_list(args.iphaz_hits, "iPhaZ funnel hits"),
            universe=universe_set,
            merged_high_confidence=merged,
            expected_union_count=funnel_pins.get("cys_funnel_union_count"),
            expected_sum_of_reason_counts=funnel_pins.get(
                "cys_funnel_sum_of_reason_counts"
            ),
        )

    # check 4: the gRodon manifest delta
    if args.grodon_eligible is None or args.grodon_manifest is None:
        skipped["check_grodon_manifest_delta"] = (
            "needs --grodon-eligible and --grodon-manifest (and optionally "
            "--grodon-dispositions)"
        )
    else:
        eligible = set(
            read_accession_list(args.grodon_eligible, "gRodon eligible positives")
        )
        _, by_status, _ = read_genome_status_table(args.grodon_manifest)
        manifest_positives = set(by_status.get("phaZ_positive", []))
        reasons = (
            read_disposition_reasons(args.grodon_dispositions)
            if args.grodon_dispositions is not None
            else {}
        )
        excluded = eligible - manifest_positives
        exclusion_reasons = {
            genome_id: reason
            for genome_id, reason in reasons.items()
            if genome_id in excluded
        }
        failure_names = {"too_few_ribosomal_hits", "hmmsearch_failed", "prediction_failed"}
        prediction_failures = {
            genome_id: reason
            for genome_id, reason in reasons.items()
            if genome_id in excluded and reason.split(":")[0] in failure_names
        }
        exclusion_reasons = {
            genome_id: reason
            for genome_id, reason in exclusion_reasons.items()
            if genome_id not in prediction_failures
        }
        grodon_pins = routed.get("check_grodon_manifest_delta", {})
        checks["check_grodon_manifest_delta"] = check_grodon_manifest_delta(
            eligible_positives=eligible,
            manifest_positives=manifest_positives,
            exclusion_reasons=exclusion_reasons,
            prediction_failures=prediction_failures,
            expected_eligible_count=grodon_pins.get("grodon_eligible_positive_count"),
            expected_manifest_count=grodon_pins.get("grodon_manifest_positive_count"),
            expected_difference=grodon_pins.get("grodon_delta"),
        )

    # check 5: motif coverage versus resolvability
    if args.motif_table is None:
        skipped["check_motif_coverage_vs_resolvability"] = "no --motif-table supplied"
    else:
        layer = read_motif_layer(args.motif_table)
        motif_pins = routed.get("check_motif_coverage_vs_resolvability", {})
        checks["check_motif_coverage_vs_resolvability"] = (
            check_motif_coverage_vs_resolvability(
                covered=layer["covered"],
                resolvable=layer["resolvable"],
                criteria_resolved={
                    criterion: payload["resolvable"]
                    for criterion, payload in layer["by_criterion"].items()
                },
                covered_but_unresolvable_is_expected=True,
            )
        )

    # check 6: the profile inventory
    if (
        args.authority is None
        or args.family_definitions is None
        or args.profile_manifest is None
    ):
        skipped["check_profile_inventory"] = (
            "needs --authority, --family-definitions and --profile-manifest"
        )
    else:
        profile_pins = routed.get("check_profile_inventory", {})
        checks["check_profile_inventory"] = check_profile_inventory(
            authority_rows=read_authority_rows(args.authority),
            family_rows=read_family_definition_rows(args.family_definitions),
            profile_rows=read_profile_rows(args.profile_manifest),
            expected_superfamilies=profile_pins.get(
                "profile_superfamily_count", REVIEW_SUPERFAMILIES
            ),
            expected_families=profile_pins.get("profile_family_count", REVIEW_FAMILIES),
            expected_trained=profile_pins.get(
                "profile_trained_count", REVIEW_TRAINED_PROFILES
            ),
            expected_reference_only=profile_pins.get(
                "profile_reference_only_count", REVIEW_REFERENCE_ONLY_PROFILES
            ),
        )

    # expectation pins that belong to a check that could not run are an error:
    # a pin that was never evaluated must not look satisfied.
    unresolved_pins = {
        name: value
        for name, value in pins.items()
        if any(
            check_name not in checks
            for check_name, keys in CHECK_EXPECTATION_KEYS.items()
            if name in keys
        )
    }

    # the reason sets of the flow itself, always reported as an overlap matrix
    hold_sets: dict[str, set[str]] = {}
    for row in rows:
        if row["hold_reasons"]:
            for reason in row["hold_reasons"].split(";"):
                hold_sets.setdefault(reason, set()).add(row["accession"])
    demote_sets: dict[str, set[str]] = {}
    for row in rows:
        if row["demote_reasons"]:
            for reason in row["demote_reasons"].split(";"):
                demote_sets.setdefault("demote:" + reason, set()).add(row["accession"])
    reason_sets: dict[str, set[str]] = {}
    reason_sets.update(hold_sets)
    reason_sets.update(demote_sets)
    if args.phac_hits is not None and args.iphaz_hits is not None:
        reason_sets["funnel_PhaC"] = set(
            read_accession_list(args.phac_hits, "PhaC funnel hits")
        )
        reason_sets["funnel_iPhaZ"] = set(
            read_accession_list(args.iphaz_hits, "iPhaZ funnel hits")
        )
    if reason_sets:
        matrix_rows = overlap_matrix(reason_sets)
        write_tsv(output_dir / "overlap_matrix.tsv", OVERLAP_FIELDS, matrix_rows)
        reason_summary = summarize_reason_sets(reason_sets)
    else:
        matrix_rows = []
        write_tsv(output_dir / "overlap_matrix.tsv", OVERLAP_FIELDS, [])
        reason_summary = None

    for name, payload in checks.items():
        write_json(output_dir / ("%s.json" % name), payload)

    statuses = {name: payload.get("status") for name, payload in checks.items()}
    check_mismatches = sorted(name for name, status in statuses.items() if status != STATUS_OK)

    if unresolved_pins:
        unresolved = check_custom_counts({}, unresolved_pins)
        checks["check_custom_counts"] = unresolved
        write_json(output_dir / "check_custom_counts.json", unresolved)
        statuses["check_custom_counts"] = unresolved["status"]
        check_mismatches.append("check_custom_counts")

    flow_closes = (
        int(summary["unaccounted"]) == 0 and int(summary["multiply_disposed"]) == 0
    )
    if flow_closes and not check_mismatches:
        exit_code = EXIT_OK
    elif not flow_closes:
        exit_code = EXIT_FLOW_DOES_NOT_CLOSE
    else:
        exit_code = EXIT_CHECK_MISMATCH

    checks_summary = {
        "schema_version": SCHEMA_VERSION,
        "universe_count": summary["universe_count"],
        "unaccounted": summary["unaccounted"],
        "multiply_disposed": summary["multiply_disposed"],
        "demoted_from_v1": summary["demoted_from_v1"],
        "hold_reason_counts_nonexclusive": summary["hold_reason_counts_nonexclusive"],
        "hold_reason_union_count": summary["hold_reason_union_count"],
        "reason_set_summary": reason_summary,
        "check_statuses": statuses,
        "check_mismatches": sorted(set(check_mismatches)),
        "checks_skipped": skipped,
        "exit_code": exit_code,
        "failure_contract": (
            "every check returns status ok|mismatch and never raises for a scientific "
            "inconsistency; assert_check_ok() raises ValueError and the CLI exits "
            "non-zero (2 when the flow does not close, 1 when a check mismatches)"
        ),
    }
    write_json(output_dir / "checks_summary.json", checks_summary)

    v1_impact = None
    if args.v1_high_confidence is not None:
        v1 = set(read_accession_list(args.v1_high_confidence, "v1 high confidence"))
        v2_high = {row["accession"] for row in rows
                   if row["primary_disposition"] in HIGH_CONFIDENCE_DISPOSITIONS}
        v1_impact = {
            "v1_high_confidence_count": len(v1),
            "v2_high_confidence_count": len(v2_high),
            "still_high_confidence_count": len(v1 & v2_high),
            "demoted_from_v1_count": len(v1 - v2_high),
            "promoted_from_v1_count": len(v2_high - v1),
            "v1_outside_universe_count": len(v1 - universe_set),
        }
        write_json(output_dir / "v1_impact.json", v1_impact)

    inputs = {
        "universe": bound_file(universe_path),
        "records": bound_file(records_path),
    }
    for name, path in (
        ("hold_reasons", args.hold_reasons),
        ("v1_high_confidence", args.v1_high_confidence),
        ("high_confidence", args.high_confidence),
        ("cys_demoted", args.cys_demoted),
        ("phac_hits", args.phac_hits),
        ("iphaz_hits", args.iphaz_hits),
        ("grodon_eligible", args.grodon_eligible),
        ("grodon_manifest", args.grodon_manifest),
        ("grodon_dispositions", args.grodon_dispositions),
        ("motif_table", args.motif_table),
        ("authority", args.authority),
        ("family_definitions", args.family_definitions),
        ("profile_manifest", args.profile_manifest),
    ):
        if path is not None:
            inputs[name] = bound_file(path)

    run_summary = {
        "schema_version": SCHEMA_VERSION,
        "tool": "reconcile_phaded_candidate_flow",
        "purpose": (
            "close the candidate flow ledger and explain every documented count "
            "discrepancy as an explicit accession-level set difference"
        ),
        "universe_count": summary["universe_count"],
        "primary_disposition_counts": summary["primary_disposition_counts"],
        "unaccounted": summary["unaccounted"],
        "unaccounted_accessions": summary["unaccounted_accessions"],
        "multiply_disposed": summary["multiply_disposed"],
        "multiply_disposed_accessions": summary["multiply_disposed_accessions"],
        "demoted_from_v1": summary["demoted_from_v1"],
        "hold_reason_counts_nonexclusive": summary["hold_reason_counts_nonexclusive"],
        "sum_of_hold_reason_counts": summary["sum_of_hold_reason_counts"],
        "hold_reason_union_count": summary["hold_reason_union_count"],
        "deferred_reason_counts": summary["deferred_reason_counts"],
        "reason_set_summary": reason_summary,
        "checks": checks,
        "checks_skipped": skipped,
        "v1_impact": v1_impact,
        "exit_code": exit_code,
        "expectations": pins,
        "inputs": inputs,
        "outputs": sorted(
            path.name for path in output_dir.iterdir() if path.is_file()
        ),
        "boundary": (
            "candidate accounting only: no new scoring, no registry change, no family "
            "call, no deletion, no analysis run, no dated deploy directory"
        ),
    }
    write_json(output_dir / "reconcile_summary.json", run_summary)
    return run_summary


def bound_file(path: Path) -> dict[str, Any]:
    """Bind one input file with its size and SHA-256 (never fabricate a hash)."""
    source = Path(path)
    if not source.is_file():
        return {"path": str(source), "status": "pending", "sha256": None, "size": None}
    return {
        "path": str(source),
        "status": "verified",
        "sha256": sha256_file(source),
        "size": source.stat().st_size,
    }


def main(argv: list[str] | None = None) -> int:
    """Run the ledger and checks, print the JSON summary and return the exit code.

    Exit codes: ``0`` all supplied checks reconcile and the flow closes, ``1`` a
    supplied check reports ``status="mismatch"``, ``2`` the flow itself does not
    close (``unaccounted`` or ``multiply_disposed`` is non-zero).  The reports are
    written before the exit code is returned, so evidence is never lost.
    """
    args = build_arg_namespace(argv)
    payload = run(args)
    printable = {
        key: payload[key]
        for key in (
            "universe_count",
            "primary_disposition_counts",
            "unaccounted",
            "multiply_disposed",
            "demoted_from_v1",
            "hold_reason_counts_nonexclusive",
            "hold_reason_union_count",
            "checks_skipped",
            "exit_code",
        )
    }
    printable["checks"] = payload["checks"]
    print(json.dumps(printable, ensure_ascii=False, sort_keys=True))
    return int(payload["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
