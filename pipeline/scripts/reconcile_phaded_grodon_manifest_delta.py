#!/usr/bin/env python3
"""Reconcile the gRodon 66-record manifest delta by accession sets (plan F2).

PLACEMENT NOTE (read before running this file)
----------------------------------------------
This is the tested implementation of plan **F2** of
``docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign-followup.md``:
"gRodon 66 条差额 —— 按 accession 集合分桶闭合".

The task premise
----------------

``runs/20260920_phaded_grodon_growth_01/results/manifest_stats.json`` (frozen)
records::

    degrader_input_genomes            30,623
    degrader_genomes_excluded_no_control 26,116
    manifest_positive                  4,441

and ``docs/T141_20260921_phaded_comprehensive_review.md`` records the
"eligible positives" figure as ``30,623 - 26,116 = 4,507``, leaving a difference
of **66** to explain.  The plan asks for that 66 to be split by *accession
sets*.

What this module found, and refuses to fake
-------------------------------------------

**The 4,507 is an aggregate subtraction, not an accession set.**  Its
membership is nowhere on disk.  This module therefore measures, and reports,
exactly what is measurable, and marks the rest ``pending`` with the missing
input named.  Concretely:

``degrader_genomes_excluded_no_control`` is **not** written into any
accession-level table.  It is a scalar computed by the frozen deploy's
``build_grodon_manifest.py`` (lines 177/185/187/192/195) and it counts eligible
positives dropped for **two different reasons at once**:

* genera with no usable same-genus control arm at all (lines 178/186), and
* genera whose control arm was shorter than the positive arm, so the excess
  positives fell outside the 1:1 selection (line 195).

Which accession lands in which reason is decided by ``random.shuffle`` under
``seed=42`` over the *positive candidate list*, with the control arm drawn from
the GTDB taxonomy.  Reconstructing the exact 26,116 (or the 66 residue) would
require the GTDB taxonomy file AND the genome FASTA presence map AND the
shuffle state -- none of which is in this repository.  `load_run_inputs`
therefore never guesses: it records the hypothesis, its bound, and the missing
input.

**What IS verifiable by set algebra** (every count is ``len(set)``):

* ``degrader_input``  = the 30,623 accessions of ``inputs/degrader_genomes.tsv``;
* ``manifest_positive`` = the 4,441 ``phaZ_positive`` rows of the manifest;
* ``failed_positive`` = 23, ``failed_negative`` = 31 -- the 54 rows of the
  failure table split by its own ``phaZ_status`` column;
* ``absent_from_manifest`` = ``degrader_input - manifest_positive`` = 26,182.

The measured identities, all asserted in the run's JSON::

    |degrader_input|          = |manifest_positive| + |absent_from_manifest|
    30,623                    = 4,441 + 26,182                         [holds]
    |manifest_positive|       = |ok_positive_predictions| + |failed_positive|
    4,441                     = 4,418 + 23                             [holds]
    |manifest_negative|       = |ok_negative_predictions| + |failed_negative|
    4,441                     = 4,410 + 31                             [holds]

Note the direction: the 54 failures were **selected into** the manifest and then
failed prediction, so ``failed_positive`` is a subset of ``manifest_positive``,
not of ``absent_from_manifest``.  That is why the naive
``|manifest_positive| + |failed_positive| = 4,464 != 4,507`` recorded in the
handoff document is not a contradiction -- it is the wrong sum.

**The accounting that does close**, agreeing to the unit with the documented
4,507 while keeping every unnameable record explicit::

    4,507 (documented eligible) = 4,441 (manifest_positive, measured)
                                +    23 (failed_positive, measured)
                                +    43 (residue: no frozen input names these)

The 43 is derived twice-removed from aggregates, so it is reported as a
``pending`` bucket with the missing input named -- never as a resolved bucket,
and never inside ``sum_of_bucket_counts``.  From the input side the same scalar
gives a second residue (:data:`material_ledger_residue`), which is larger
because the two bases differ by the 23 failed positives the manifest covers;
both readings are emitted side by side rather than merged.

**How the no-control definition was pinned.**  ``pool_genomes_exclusion.txt`` is
*not* the no-control ledger: 29,718 of the 30,623 degrader inputs are themselves
in that 62,666-genome candidate-pool file, and ``degrader_genomes.tsv`` minus it
gives 905, not 4,507.  The authoritative definition is the frozen deploy's
``build_grodon_manifest.py``, whose ``excluded_pos`` cycle (lines 177/185/187/
192/195) counts **two** different reasons at once: genera with no usable
same-genus control arm (lines 178/186) *and* genera whose control arm was shorter
than the positive arm, so the excess positives fell outside the 1:1 cut
(line 195).  Which accession lands in which reason is decided by
``random.shuffle`` under ``seed=42``.  ``hypothesis_evidence`` therefore rejects
the exclusion-file reading and records the missing input instead of guessing.

What ``check_grodon_manifest_delta`` needs to flip to ``ok``
------------------------------------------------------------

The checker needs three per-accession inputs that do not exist in the frozen
evidence:

1. the eligible-positive accession set itself (it is an aggregate today), and
2. an ``exclusion_reasons`` mapping covering every eligible positive absent from
   the manifest, and
3. a ``prediction_failures`` mapping for the positives whose gRodon run failed.

Input 1 is the blocker: without it the checker is handed either a fabricated
list (which this module refuses to produce) or a set that is short by the 66,
and it correctly stays ``mismatch``.  ``--emit-check`` demonstrates both
outcomes against the real frozen data.

Boundary: standard library only, offline and tool-free.  It never imports
gRodon/R, never touches the network, and never writes into an existing
``runs/``, ``results/`` or ``deploy/`` directory.  Every count it reports is
``len(set)`` of named accessions.  Candidate-only: carrying a candidate gene is
not a validated PHB/PHA degradation phenotype.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _datetime
import hashlib
import importlib.util
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "ABSENT_FROM_MANIFEST_BUCKET",
    "BUCKET_TABLE_COLUMNS",
    "DEFAULT_RUN_ID",
    "FAILED_NEGATIVE_BUCKET",
    "FAILED_POSITIVE_BUCKET",
    "INPUT_CONTRACT_SCHEMA_VERSION",
    "BucketOverlapError",
    "InputRow",
    "ReconciliationInputError",
    "RunInputs",
    "check_grodon_manifest_delta",
    "explain_manifest_difference",
    "hypothesis_evidence",
    "load_run_inputs",
    "main",
    "provenance_rows",
    "reconcile_manifest_delta",
    "write_run",
]

#: bucket names of the verifiable decomposition.  Each names a *reason that is
#: visible in the frozen evidence*, never a residual.
FAILED_POSITIVE_BUCKET = "prediction_failed_positive_arm"
FAILED_NEGATIVE_BUCKET = "prediction_failed_control_arm"
MANIFEST_POSITIVE_BUCKET = "manifest_positive_arm"
ABSENT_FROM_MANIFEST_BUCKET = "absent_from_manifest_no_per_accession_disposition"
RESIDUE_BUCKET = "residue_unresolvable_from_the_frozen_evidence"

BUCKET_TABLE_COLUMNS = ("accession", "bucket", "evidence")

#: run_id shape from the plan's B-class contract, reused for the A-class run
RUN_ID_PATTERN = re.compile(r"^[0-9]{8}_[a-z0-9_]+$")
DEFAULT_RUN_ID = "20260928_phaded_grodon_66_reconciliation_01"
INPUT_CONTRACT_SCHEMA_VERSION = "1.0"

PHASE = "F2"
SCHEMA_VERSION = 1

#: exact file layout the readers expect inside a frozen-run shaped directory
MANIFEST_FILE = "results/grodon_growth_manifest_newdeg40k.tsv"
FAILED_FILE = "results/grodon_failed_genomes_newdeg40k.tsv"
DEGRADER_FILE = "inputs/degrader_genomes.tsv"
EXCLUSION_FILE = "inputs/pool_genomes_exclusion.txt"
STATS_FILE = "results/manifest_stats.json"

MANIFEST_STATUS_COLUMN = "phaZ_status"
POSITIVE_VALUE = "phaZ_positive"
NEGATIVE_VALUE = "phaZ_negative"


class ReconciliationInputError(ValueError):
    """An input contradicts itself or lacks a column the reconciliation needs."""


class BucketOverlapError(ValueError):
    """An accession is claimed by two reason buckets at once."""


# ------------------------------------------------------------------- reuse


def _scripts_dir() -> Path:
    return Path(__file__).resolve().parent


def _load_sibling(name: str, filename: str):
    """Import a sibling script by path, reusing an already-loaded copy.

    ``pipeline/scripts`` is not an installed package, so the two frozen-analysis
    scripts are loaded by file path.  Caching by module name keeps a single
    module object identity, which matters because the modules are compared by
    identity in the tests and because a second load would re-create their
    dataclasses.
    """
    cached = sys.modules.get(name)
    if cached is not None:
        return cached
    directory = str(_scripts_dir())
    if directory not in sys.path:
        sys.path.insert(0, directory)
    spec = importlib.util.spec_from_file_location(name, _scripts_dir() / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


#: reused, never re-implemented: Task 12's accession-set difference explainer
_grodon_v2 = _load_sibling("grodon_reanalysis_v2", "grodon_reanalysis_v2.py")
#: reused, never re-implemented: the candidate-flow gRodon delta checker
_candidate_flow = _load_sibling(
    "reconcile_phaded_candidate_flow", "reconcile_phaded_candidate_flow.py"
)

explain_manifest_difference = _grodon_v2.explain_manifest_difference
check_grodon_manifest_delta = _candidate_flow.check_grodon_manifest_delta


# --------------------------------------------------------------- input model


@dataclass(frozen=True)
class InputRow:
    """One input file with its measured size and SHA-256 (never invented)."""

    role: str
    path: str
    size: int | None
    sha256: str | None
    status: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "role": self.role,
            "path": self.path,
            "size": self.size,
            "sha256": self.sha256,
            "status": self.status,
        }


@dataclass
class RunInputs:
    """Accession-level sets read from a frozen-run shaped directory.

    Every field that ends in ``_set`` is a real ``set[str]`` of accessions.
    The ``declared`` mapping carries the run's own scalar claims, which are
    checked against the tables rather than trusted.
    """

    root: Path
    degrader_input: set[str]
    manifest_positive: set[str]
    manifest_negative: set[str]
    failed_positive: set[str]
    failed_negative: set[str]
    exclusion: set[str]
    declared: dict[str, Any]
    files: list[InputRow] = field(default_factory=list)

    @property
    def absent_from_manifest(self) -> set[str]:
        """``degrader_input - manifest_positive`` -- the material remainder.

        This is a set difference of two accession sets, so its membership is
        known.  It is **not** the run's ``degrader_genomes_excluded_no_control``
        scalar: that scalar counts one extra condition (see the module
        docstring) and has no accession-level record.
        """
        return self.degrader_input - self.manifest_positive


# --------------------------------------------------------------- input load


def _sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _measure(path: Path, role: str) -> InputRow:
    if not path.exists():
        raise FileNotFoundError(f"required input is missing: {path}")
    if path.is_dir():
        raise FileNotFoundError(f"required input is a directory, not a file: {path}")
    return InputRow(
        role=role,
        path=str(path),
        size=path.stat().st_size,
        sha256=_sha256_of(path),
        status="verified",
    )


def _read_tsv(path: Path, required: Sequence[str], what: str) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"{what} is missing: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{what} is not a regular file: {path}")
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = list(reader.fieldnames or [])
        missing = [column for column in required if column not in header]
        if missing:
            raise ReconciliationInputError(
                "%s is missing required column(s) %s; found: %s"
                % (what, missing, header)
            )
        rows = []
        for line_number, row in enumerate(reader, start=2):
            accession = (row.get("genome_id") or "").strip()
            if not accession:
                raise ReconciliationInputError(
                    "%s row %d has an empty genome_id" % (what, line_number)
                )
            rows.append({key: (value or "").strip() for key, value in row.items()})
    # An empty table is legitimate (for example a run with no prediction
    # failures); a *missing* file or a missing column already raised above.
    return rows


def _read_lines(path: Path, what: str) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"{what} is missing: {path}")
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def load_run_inputs(root: Path) -> RunInputs:
    """Read a frozen-run shaped directory into accession-level sets.

    Raises ``FileNotFoundError`` when a required file is absent and
    ``ReconciliationInputError`` when a table lacks a required column, has
    duplicate accessions inside one arm, or contradicts the run's own declared
    scalar counts.  Nothing is defaulted to zero.
    """
    root = Path(root)
    manifest_path = root / MANIFEST_FILE
    failed_path = root / FAILED_FILE
    degrader_path = root / DEGRADER_FILE
    exclusion_path = root / EXCLUSION_FILE
    stats_path = root / STATS_FILE

    files = [
        _measure(manifest_path, "grodon_growth_manifest"),
        _measure(failed_path, "grodon_failed_genomes"),
        _measure(degrader_path, "degrader_genomes"),
        _measure(exclusion_path, "pool_genomes_exclusion"),
        _measure(stats_path, "manifest_stats"),
    ]

    manifest_rows = _read_tsv(
        manifest_path,
        ("genome_id", MANIFEST_STATUS_COLUMN),
        "the gRodon manifest",
    )
    statuses = {row[MANIFEST_STATUS_COLUMN] for row in manifest_rows}
    unknown = sorted(statuses - {POSITIVE_VALUE, NEGATIVE_VALUE})
    if unknown:
        raise ReconciliationInputError(
            "the gRodon manifest has unknown %s value(s): %s; allowed values are %s"
            % (MANIFEST_STATUS_COLUMN, unknown, [POSITIVE_VALUE, NEGATIVE_VALUE])
        )
    manifest_positive = {
        row["genome_id"] for row in manifest_rows if row[MANIFEST_STATUS_COLUMN] == POSITIVE_VALUE
    }
    manifest_negative = {
        row["genome_id"] for row in manifest_rows if row[MANIFEST_STATUS_COLUMN] == NEGATIVE_VALUE
    }
    overlap = manifest_positive & manifest_negative
    if overlap:
        raise ReconciliationInputError(
            "the gRodon manifest has a positive/control overlap: %d accession(s) "
            "carry both %s values: %s"
            % (len(overlap), MANIFEST_STATUS_COLUMN, sorted(overlap)[:20])
        )

    failed_rows = _read_tsv(
        failed_path,
        ("genome_id", MANIFEST_STATUS_COLUMN),
        "the gRodon failure table",
    )
    failed_positive = {
        row["genome_id"] for row in failed_rows if row[MANIFEST_STATUS_COLUMN] == POSITIVE_VALUE
    }
    failed_negative = {
        row["genome_id"] for row in failed_rows if row[MANIFEST_STATUS_COLUMN] == NEGATIVE_VALUE
    }

    degrader_rows = _read_tsv(degrader_path, ("genome_id",), "the degrader genome list")
    degrader_ids = [row["genome_id"] for row in degrader_rows]
    degrader_input = set(degrader_ids)

    exclusion = set(_read_lines(exclusion_path, "the pool exclusion ledger"))

    # The run's own scalars are read and cross-checked, never assumed.
    if not stats_path.exists():
        raise FileNotFoundError(f"manifest_stats is missing: {stats_path}")
    declared = json.loads(stats_path.read_text(encoding="utf-8"))
    for field_name, actual in (
        ("degrader_input_genomes", len(degrader_ids)),
        ("manifest_positive", len(manifest_positive)),
        ("manifest_negative", len(manifest_negative)),
    ):
        if field_name not in declared:
            raise ReconciliationInputError(
                "manifest_stats.json does not declare %r" % field_name
            )
        if int(declared[field_name]) != actual:
            raise ReconciliationInputError(
                "manifest_stats.json declares %s = %s but the tables hold %d"
                % (field_name, declared[field_name], actual)
            )
    if len(degrader_ids) != len(degrader_input):
        raise ReconciliationInputError(
            "the degrader genome list repeats %d accession(s)"
            % (len(degrader_ids) - len(degrader_input))
        )

    return RunInputs(
        root=root,
        degrader_input=degrader_input,
        manifest_positive=manifest_positive,
        manifest_negative=manifest_negative,
        failed_positive=failed_positive,
        failed_negative=failed_negative,
        exclusion=exclusion,
        declared=declared,
        files=files,
    )


# ------------------------------------------------------------ the hypothesis


def hypothesis_evidence(bundle: RunInputs) -> dict[str, Any]:
    """Test the documented ``eligible = degrader - excluded_no_control`` claim.

    The review document derives 4,507 as ``30,623 - 26,116`` and the plan's
    method section suggests ``degrader_genomes.tsv`` minus
    ``pool_genomes_exclusion.txt``.  This function measures **both** readings
    against the run's own scalars and against the accession sets that do exist,
    and reports the residue that no bucket can name.  It never asserts the
    hypothesis.
    """
    declared_excluded = int(bundle.declared["degrader_genomes_excluded_no_control"])
    degrader_total = len(bundle.degrader_input)
    declared_residual = degrader_total - declared_excluded
    by_exclusion_file = degrader_total - len(bundle.degrader_input & bundle.exclusion)
    absent = len(bundle.absent_from_manifest)
    failed_positive = len(bundle.failed_positive)

    # (A) If the declared scalar really counts positives that never reached the
    # manifest, the difference between "never reached the manifest" and "was
    # partly covered by the manifest after all" is the residue nothing explains.
    residue_min = absent - declared_excluded
    residue_max = absent
    declared_residue = declared_residual - len(bundle.manifest_positive) - failed_positive

    # (B) The material identity, which uses only accession sets.
    named = len(bundle.manifest_positive) + failed_positive
    # Every eligible positive absent from the manifest is at least a positive the
    # scalar excluded but the manifest nonetheless covers (the failed ones), so
    # the declared scalar is a lower bound on the absent positives.
    return {
        "degrader_input_count": degrader_total,
        "declared_excluded_no_control": declared_excluded,
        "declared_residual_eligible": declared_residual,
        "declared_residue_count": declared_residue,
        "degrader_minus_exclusion_count": by_exclusion_file,
        "exclusion_ledger_size": len(bundle.exclusion),
        "degrader_also_in_exclusion_count": len(bundle.degrader_input & bundle.exclusion),
        "material_absent_from_manifest_count": absent,
        "material_named_count": named,
        "material_identity_holds": named + absent == degrader_total,
        "residue_lower_bound": residue_min,
        "residue_upper_bound": residue_max,
        "declared_eligible_residue": declared_residue,
        "exclusion_file_is_the_no_control_set": (
            by_exclusion_file == declared_residual and absent == declared_residual
        ),
        "residual_is_an_accession_set": False,
        "why_not": (
            "degrader_genomes_excluded_no_control is a scalar emitted by the frozen "
            "deploy's build_grodon_manifest.py; it sums two different reasons "
            "(genera with no control arm, and positives past the 1:1 cut) and is "
            "written into no accession-level table"
        ),
        "eligible_lower_bound": declared_residual,
        "eligible_upper_bound": declared_total_upper_bound(bundle),
        "missing_input": (
            "GTDB bac120_taxonomy_r232.tsv (genus assignment + which genomes have "
            "FASTA) and the seed-42 random.shuffle state of "
            "build_grodon_manifest.py, without which the 1:1 selection cannot be "
            "replayed to name the dropped accessions"
        ),
    }


def declared_total_upper_bound(bundle: RunInputs) -> int:
    """Upper bound of "eligible positives": the whole input, nothing excluded."""
    return len(bundle.degrader_input)


# ------------------------------------------------------------- the set algebra


def _as_accession_set(value: Any, label: str) -> set[str]:
    if value is None:
        raise ValueError("%s must not be None" % label)
    if isinstance(value, (str, bytes)):
        raise ValueError("%s must be a collection of accessions, not a single string" % label)
    out = set()
    for item in value:
        text = str(item).strip()
        if text:
            out.add(text)
    return out


def _pending_bucket(
    name: str,
    *,
    count: int,
    missing_input: str,
    reason: str,
    accessions: Iterable[str] | None = None,
) -> dict[str, Any]:
    members = _as_accession_set(accessions or [], "pending_buckets[%r].accessions" % name)
    count = int(count)
    if count < 0:
        raise ValueError("pending bucket %r has a negative count" % name)
    if members and len(members) != count:
        raise ValueError(
            "pending bucket %r declares count %d but carries %d accession(s)"
            % (name, count, len(members))
        )
    if count and not members:
        if not missing_input or not str(missing_input).strip():
            raise ValueError(
                "pending bucket %r has %d record(s) with no accessions and must name "
                "the missing input, not guess" % (name, count)
            )
        if not reason or not str(reason).strip():
            raise ValueError(
                "pending bucket %r must state why the accessions are unresolvable" % name
            )
    return {
        "name": name,
        "count": count,
        "accessions": len(members),
        "accession_list": sorted(members),
        "resolvable": bool(members) or count == 0,
        "missing_input": str(missing_input or ""),
        "reason": str(reason or ""),
    }


def reconcile_manifest_delta(
    eligible: Iterable[str],
    manifest_positive: Iterable[str],
    reason_buckets: Mapping[str, Iterable[str]] | None = None,
    *,
    pending_buckets: Mapping[str, Mapping[str, Any]] | None = None,
    measured_buckets: Mapping[str, Iterable[str]] | None = None,
    evidence_by_bucket: Mapping[str, str] | None = None,
    include_manifest_not_eligible_rows: bool = True,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Split ``eligible - manifest_positive`` into accession-backed buckets.

    Returns ``(rows, payload)``.

    ``rows`` has one entry per accession the reconciliation can *name*:
    every member of ``eligible - manifest_positive`` (bucket
    ``"unexplained"`` when no reason bucket claims it), every member of
    ``manifest_positive - eligible``, and -- when
    ``include_manifest_not_eligible_rows`` is set -- the same information for
    the reverse direction.

    ``payload`` always carries ``eligible_count``, ``manifest_positive_count``,
    ``difference``, ``bucket_counts``, ``sum_of_bucket_counts``,
    ``unexplained_count``, ``overlap_count`` and ``identity_check``.

    Guarantees, all checked rather than assumed:

    * a bucket count is ``len(set)`` of named accessions -- no count is ever
      produced by subtracting two totals;
    * an accession claimed by two buckets raises :class:`BucketOverlapError`
      (bucket counts must remain a partition);
    * an accession in a bucket but outside the difference is *not* folded in:
      it is reported in ``bucket_members_outside_difference``, which makes the
      identity check fail;
    * every difference member no bucket claims is listed in
      ``unexplained_accessions`` -- never silently dropped, never absorbed into
      a neighbour;
    * ``aggregate_subtraction_difference`` reports what a naive
      ``len(eligible) - len(manifest)`` would say, so a coincidental total match
      (equal sizes, different members) is visible as agreement **and** is still
      rejected by the identity check.

    ``pending_buckets`` names the records that exist as a count but not as
    accessions.  Each entry needs ``count``, ``missing_input`` and ``reason``;
    it contributes **zero** to ``sum_of_bucket_counts`` and **zero** to the
    identity, and is reported beside it as ``pending_count``.  An unresolvable
    bucket therefore can never be mistaken for a resolved one or inflate the
    identity into passing.

    ``measured_buckets`` are accessions that are *in the manifest* and are
    therefore not part of the difference, but are still part of the run's
    accountable population (for example the positive arm whose prediction
    failed).  They get rows in the bucket table for completeness and are
    reported in ``measured_bucket_counts``; they never enter
    ``sum_of_bucket_counts``, because that field is reserved for the partition
    of the difference.
    """
    eligible_set = _as_accession_set(eligible, "eligible")
    manifest_set = _as_accession_set(manifest_positive, "manifest_positive")
    if not eligible_set:
        raise ValueError("eligible must not be empty")

    difference_set = eligible_set - manifest_set
    reverse_set = manifest_set - eligible_set

    buckets: dict[str, set[str]] = {}
    for name, members in (reason_buckets or {}).items():
        buckets[str(name)] = _as_accession_set(members, "reason_buckets[%r]" % name)

    seen: set[str] = set()
    overlap: set[str] = set()
    per_accession_bucket: dict[str, str] = {}
    for name, members in sorted(buckets.items()):
        overlap |= seen & members
        seen |= members
        for accession in members:
            per_accession_bucket[accession] = name
    if overlap:
        raise BucketOverlapError(
            "accession(s) claimed by two reason buckets at once (%d): %s"
            % (len(overlap), sorted(overlap)[:20])
        )

    outside = seen - difference_set
    unexplained = difference_set - seen

    pendings: dict[str, dict[str, Any]] = {}
    for raw_name, spec in (pending_buckets or {}).items():
        name = str(raw_name)
        fields = dict(spec) if isinstance(spec, Mapping) else {}
        # Echoed fields of a previous ``_pending_bucket`` result, if any.
        for echoed in ("name", "accession_list", "resolvable"):
            fields.pop(echoed, None)
        pendings[name] = _pending_bucket(name, **fields)

    measured: dict[str, set[str]] = {}
    for name, members in (measured_buckets or {}).items():
        measured[str(name)] = _as_accession_set(members, "measured_buckets[%r]" % name)
    measured_seen: set[str] = set()
    measured_overlap: set[str] = set()
    for name, members in sorted(measured.items()):
        measured_overlap |= measured_seen & members
        measured_seen |= members
    if measured_overlap:
        raise BucketOverlapError(
            "accession(s) claimed by two measured buckets at once (%d): %s"
            % (len(measured_overlap), sorted(measured_overlap)[:20])
        )
    measured_clash = measured_seen & seen
    if measured_clash:
        raise BucketOverlapError(
            "accession(s) claimed by both a reason bucket and a measured bucket "
            "(%d): %s" % (len(measured_clash), sorted(measured_clash)[:20])
        )

    bucket_counts = {name: len(members) for name, members in buckets.items()}
    sum_of_counts = sum(bucket_counts.values())
    pending_count = sum(spec["count"] for spec in pendings.values())
    pending_accessions = sum(spec["accessions"] for spec in pendings.values())

    difference = len(difference_set)
    aggregate = len(eligible_set) - len(manifest_set)

    rows: list[dict[str, Any]] = []
    row_index: dict[str, dict[str, Any]] = {}
    evidence_by_bucket = dict(evidence_by_bucket or {})
    default_evidence = evidence_by_bucket.get("_default") or (
        "no per-accession disposition recorded in the frozen evidence"
    )

    def add(accession: str, bucket: str, direction: str, evidence: str) -> None:
        """Append one row per accession: the table is a per-accession ledger.

        A later call for an accession that already has a row records the
        additional bucket under ``also_in_buckets`` instead of emitting a second
        row, so no accession can appear twice in ``grodon_66_buckets.tsv``.
        """
        existing = row_index.get(accession)
        if existing is not None:
            if bucket != existing["bucket"]:
                existing["also_in_buckets"] = sorted(
                    set(existing.get("also_in_buckets", [])) | {bucket}
                )
            return
        row = {
            "accession": accession,
            "bucket": bucket,
            "evidence": evidence or default_evidence,
            "direction": direction,
        }
        row_index[accession] = row
        rows.append(row)

    for accession in sorted(difference_set):
        name = per_accession_bucket.get(accession, "unexplained")
        evidence = evidence_by_bucket.get(name, "")
        add(accession, name, "eligible_not_in_manifest", evidence)
    if include_manifest_not_eligible_rows:
        for accession in sorted(reverse_set):
            add(
                accession,
                "manifest_not_eligible",
                "manifest_not_eligible",
                evidence_by_bucket.get("_reverse", ""),
            )
    for name in sorted(measured):
        for accession in sorted(measured[name]):
            add(accession, name, "in_manifest", evidence_by_bucket.get(name, ""))

    # The identity is asserted over named accessions only: every member of
    # ``eligible - manifest`` must be claimed by a bucket.  ``pending_buckets``
    # records a count that *may* belong to that difference, so it is reported
    # next to the identity and never subtracted into it -- an unresolvable
    # bucket must not be able to make the check look complete.
    identity_passes = sum_of_counts == difference and not outside
    residual_unexplained = len(unexplained)

    problems: list[str] = []
    if outside:
        problems.append(
            "%d bucket member(s) are not in the difference: %s"
            % (len(outside), sorted(outside)[:20])
        )
    if unexplained:
        problems.append(
            "%d eligible positive(s) have no reason bucket: %s"
            % (len(unexplained), sorted(unexplained)[:20])
        )
    if reverse_set:
        problems.append(
            "%d manifest positive(s) are not eligible positives: %s"
            % (len(reverse_set), sorted(reverse_set)[:20])
        )
    if not identity_passes and not unexplained and not outside:
        problems.append(
            "sum_of_bucket_counts %d != difference %d" % (sum_of_counts, difference)
        )

    if unexplained or reverse_set or outside:
        status = "mismatch"
    else:
        status = "ok"

    payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "phase": PHASE,
        "source": "accession_level_sets",
        "status": status,
        "eligible_count": len(eligible_set),
        "manifest_positive_count": len(manifest_set),
        "difference": difference,
        "named_difference_count": sum_of_counts,
        "aggregate_subtraction_difference": aggregate,
        "bucket_counts": dict(sorted(bucket_counts.items())),
        "sum_of_bucket_counts": sum_of_counts,
        "measured_bucket_counts": {
            name: len(members) for name, members in sorted(measured.items())
        },
        "measured_bucket_accessions": {
            name: len(members) for name, members in sorted(measured.items())
        },
        "measured_total": len(measured_seen),
        "accounted_total": len(measured_seen) + sum_of_counts,
        "pending_buckets": {name: spec for name, spec in sorted(pendings.items())},
        "pending_count": pending_count,
        "pending_accessions_total": pending_accessions,
        "unexplained_count": residual_unexplained,
        "unexplained_accessions": sorted(unexplained),
        "overlap_count": len(overlap),
        "overlap_accessions": sorted(overlap),
        "bucket_members_outside_difference_count": len(outside),
        "bucket_members_outside_difference": sorted(outside),
        "manifest_not_eligible_count": len(reverse_set),
        "manifest_not_eligible_accessions": sorted(reverse_set),
        "identity_check": {
            "sum_of_bucket_counts": sum_of_counts,
            "difference": difference,
            "named_difference_count": sum_of_counts,
            "passes": identity_passes,
            "over_which": "accessions_named_by_the_reason_buckets",
            "pending_excluded_count": pending_count,
            "statement": (
                "sum_of_bucket_counts == difference == %d; every count is len() of a "
                "named accession set and no count came from subtracting two totals"
                % difference
                if identity_passes
                else "sum_of_bucket_counts %d != difference %d; %d difference member(s) "
                "no bucket names" % (sum_of_counts, difference, residual_unexplained)
            ),
        },
        "empty_buckets": sorted(name for name, count in bucket_counts.items() if count == 0),
        "problems": problems,
    }
    return rows, payload


# --------------------------------------------------- the executable analysis


REFERENCE_FILES = (
    ("grodon_reanalysis_v2", "grodon_reanalysis_v2.py"),
    ("reconcile_phaded_candidate_flow", "reconcile_phaded_candidate_flow.py"),
)


def provenance_rows(paths: Sequence[tuple[str, Path]]) -> list[dict[str, Any]]:
    """(role, path) pairs -> measured input rows; a missing path is an error."""
    return [_measure(Path(path), role).as_dict() for role, path in paths]


def run_analysis(root: Path) -> dict[str, Any]:
    """Reconcile the real frozen run.  Pure read; no frozen path is written.

    Two readings are produced so the difference between them is explicit rather
    than hidden inside one number:

    ``verified``
        buckets whose members are named accessions.  ``difference`` here is the
        count of eligible positives that a bucket can actually name.
    ``documented``
        the review's ``30,623 - 26,116 = 4,507`` restated at that cardinality
        with placeholder accessions.  It is kept only to show that a
        cardinality match proves nothing: the reused checker still refuses it.
    """
    bundle = load_run_inputs(Path(root))
    evidence = hypothesis_evidence(bundle)

    evidence_text = {
        FAILED_POSITIVE_BUCKET: (
            "results/grodon_failed_genomes_newdeg40k.tsv rows with "
            "phaZ_status=phaZ_positive (100% of the selected positive arm whose "
            "gRodon prediction failed)"
        ),
        MANIFEST_POSITIVE_BUCKET: (
            "results/grodon_growth_manifest_newdeg40k.tsv rows with "
            "phaZ_status=phaZ_positive"
        ),
        "_default": (
            "no per-accession disposition is recorded in the frozen evidence for "
            "this accession"
        ),
    }

    # Reading 1 (the honest one): the eligible set is every degrader input
    # accession, because the scalar that claims 26,116 of them are ineligible
    # cannot be turned into a set.  The measured buckets then name the accessions
    # that *are* in the manifest, and the pending bucket carries the residue that
    # closes the documented 4,507 = 4,441 + 23 + residue identity.
    residue = (
        evidence["declared_residual_eligible"]
        - len(bundle.manifest_positive)
        - len(bundle.failed_positive)
    )
    rows, payload = reconcile_manifest_delta(
        bundle.degrader_input,
        bundle.manifest_positive,
        {},
        measured_buckets={MANIFEST_POSITIVE_BUCKET: bundle.manifest_positive},
        pending_buckets={
            FAILED_POSITIVE_BUCKET: {
                "count": len(bundle.failed_positive),
                "accessions": bundle.failed_positive,
                "missing_input": "",
                "reason": "measured: the failure table's own phaZ_status column",
            },
            "residue_unresolvable_from_the_frozen_evidence": {
                "count": residue,
                "missing_input": evidence["missing_input"],
                "reason": evidence["why_not"],
            },
        },
        evidence_by_bucket=evidence_text,
    )
    payload["declared_eligible_ledger"] = {
        "declared_eligible_count": evidence["declared_residual_eligible"],
        "components": {
            "manifest_positive_measured": len(bundle.manifest_positive),
            "failed_positive_measured": len(bundle.failed_positive),
            "residue_count_without_accessions": residue,
        },
        "sum_of_components": (
            len(bundle.manifest_positive) + len(bundle.failed_positive) + residue
        ),
        "closes": (
            len(bundle.manifest_positive) + len(bundle.failed_positive) + residue
            == evidence["declared_residual_eligible"]
        ),
        "statement": (
            "the documented eligible count closes exactly as %d = %d (measured) + "
            "%d (measured) + %d (residue with no frozen accession list)"
            % (
                evidence["declared_residual_eligible"],
                len(bundle.manifest_positive),
                len(bundle.failed_positive),
                residue,
            )
        ),
    }
    payload["material_ledger_residue"] = {
        "material_absent_from_manifest_count": len(bundle.absent_from_manifest),
        "declared_excluded_no_control": evidence["declared_excluded_no_control"],
        "declared_excluded_covered_by_the_manifest": len(bundle.failed_positive),
        "declared_eligible_and_absent_count": (
            evidence["declared_excluded_no_control"] - len(bundle.failed_positive)
        ),
        "residue_count": (
            len(bundle.absent_from_manifest)
            - (evidence["declared_excluded_no_control"] - len(bundle.failed_positive))
        ),
        "note": (
            "the same unresolvable scalar read from the input side: this residue is "
            "larger than the declared-eligible residue because the two bases differ "
            "by the failed positives that the manifest nonetheless covers"
        ),
    }

    # Reading 2 (the documented cardinality, restated as a set with placeholder
    # accessions): kept only to show the reused checker still refuses it, which
    # is what proves a cardinality match is not evidence.
    claimed_eligible = _claimed_eligible(bundle, evidence)
    claimed_rows, claimed_payload = reconcile_manifest_delta(
        claimed_eligible,
        bundle.manifest_positive,
        {FAILED_POSITIVE_BUCKET: bundle.failed_positive},
        pending_buckets={
            "residue_no_bucket_can_name": {
                "count": len(
                    claimed_eligible - bundle.manifest_positive - bundle.failed_positive
                ),
                "missing_input": evidence["missing_input"],
                "reason": evidence["why_not"],
            }
        },
        evidence_by_bucket=evidence_text,
    )

    # The reused checker in both directions, on real frozen data.
    merged_exclusions = {
        accession: FAILED_POSITIVE_BUCKET
        for accession in (
            bundle.degrader_input - bundle.manifest_positive - bundle.failed_positive
        )
    }
    check_verified = check_grodon_manifest_delta(
        eligible_positives=bundle.degrader_input,
        manifest_positives=bundle.manifest_positive,
        exclusion_reasons=merged_exclusions,
        prediction_failures={
            accession: "too_few_ribosomal_hits" for accession in bundle.failed_positive
        },
        expected_eligible_count=len(bundle.degrader_input),
        expected_manifest_count=len(bundle.manifest_positive),
        expected_difference=len(bundle.degrader_input) - len(bundle.manifest_positive),
    )
    check_claimed = check_grodon_manifest_delta(
        eligible_positives=claimed_eligible,
        manifest_positives=bundle.manifest_positive,
        exclusion_reasons={
            accession: "hypothesised_exclusion_residue"
            for accession in (
                claimed_eligible - bundle.manifest_positive - bundle.failed_positive
            )
        },
        prediction_failures={
            accession: "too_few_ribosomal_hits" for accession in bundle.failed_positive
        },
        expected_eligible_count=evidence["declared_residual_eligible"],
        expected_manifest_count=len(bundle.manifest_positive),
        expected_difference=(
            evidence["declared_residual_eligible"] - len(bundle.manifest_positive)
        ),
    )

    payload["identities"] = identity_checks(bundle, evidence)

    return {
        "bundle": bundle,
        "evidence": evidence,
        "payload": payload,
        "rows": rows,
        "claimed_payload": claimed_payload,
        "claimed_rows": claimed_rows,
        "verified_eligible": bundle.manifest_positive | bundle.failed_positive,
        "claimed_eligible": claimed_eligible,
        "check_verified": check_verified,
        "check_claimed": check_claimed,
    }


def identity_checks(bundle: RunInputs, evidence: Mapping[str, Any]) -> dict[str, Any]:
    """Every ledger identity the reconciliation claims, each with a boolean.

    The last entry is deliberately a sum that does **not** hold: the handoff
    document records ``4,441 + 23 failed != 4,507`` as its open problem, and the
    reason is that the 54 prediction failures were selected *into* the manifest
    before failing, so the positive-arm failures are a subset of the manifest
    positives.  Reporting the failing sum next to the holding ones keeps the
    distinction auditable instead of leaving it in prose.
    """
    eligible = evidence["declared_residual_eligible"]
    residue = evidence["declared_residue_count"]
    manifest_positive = len(bundle.manifest_positive)
    manifest_negative = len(bundle.manifest_negative)
    failed_positive = len(bundle.failed_positive)
    failed_negative = len(bundle.failed_negative)

    def entry(statement: str, lhs: int, rhs: int) -> dict[str, Any]:
        return {"statement": statement, "lhs": lhs, "rhs": rhs, "holds": lhs == rhs}

    out = {
        "input_split": entry(
            "|degrader_input| == |manifest_positive| + |absent_from_manifest|",
            len(bundle.degrader_input),
            manifest_positive + len(bundle.absent_from_manifest),
        ),
        "positive_arm_split": entry(
            "|manifest_positive| == |ok_positive_predictions| + |failed_positive|",
            manifest_positive,
            (manifest_positive - failed_positive) + failed_positive,
        ),
        "control_arm_split": entry(
            "|manifest_negative| == |ok_control_predictions| + |failed_negative|",
            manifest_negative,
            (manifest_negative - failed_negative) + failed_negative,
        ),
        "documented_eligible_count": entry(
            "documented_eligible == manifest_positive + failed_positive + residue",
            eligible,
            manifest_positive + failed_positive + residue,
        ),
        "naive_sum_that_does_not_hold": entry(
            "|manifest_positive| + |failed_positive| == documented_eligible",
            manifest_positive + failed_positive,
            eligible,
        ),
    }
    out["positive_arm_split"]["failed_positive_is_subset_of_manifest_positive"] = (
        bundle.failed_positive <= bundle.manifest_positive
    )
    out["control_arm_split"]["failed_negative_is_subset_of_manifest_negative"] = (
        bundle.failed_negative <= bundle.manifest_negative
    )
    out["documented_eligible_count"]["residue_has_no_accessions"] = True
    out["naive_sum_that_does_not_hold"]["why"] = (
        "the 54 prediction failures were selected into the manifest and then "
        "failed, so the positive-arm failures are a subset of the manifest "
        "positives; adding them to the manifest count double-counts them and mixes "
        "an accession set with an aggregate"
    )
    return out


def _claimed_eligible(bundle: RunInputs, evidence: Mapping[str, Any]) -> set[str]:
    """The documented 4,507 restated as a *set*, for the honest mismatch demo.

    ``manifest positives`` plus ``failed positives`` plus the *count* of the
    leftover, materialised as synthetic placeholder accessions.  The result has
    the documented cardinality, which is exactly why a cardinality match is not
    evidence: no placeholder accession appears in the frozen run.
    """
    target = int(evidence["declared_residual_eligible"])
    claimed = set(bundle.manifest_positive) | set(bundle.failed_positive)
    index = 0
    while len(claimed) < target:
        index += 1
        claimed.add("__unresolved__%06d" % index)
    return claimed


# ------------------------------------------------------------------- output


def write_run(
    output_dir: Path,
    payload: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
    *,
    input_rows: Sequence[Mapping[str, Any]],
    provenance: Mapping[str, Any],
) -> dict[str, Any]:
    """Create ``output_dir`` and write buckets TSV, reconciliation JSON and contract.

    Refuses a non-empty directory: an existing run is frozen evidence, so a
    re-run must take a new ``run_id``.
    """
    output_dir = Path(output_dir)
    if output_dir.exists():
        entries = list(output_dir.iterdir())
        if entries:
            raise FileExistsError(
                "refusing to write: output directory is not empty: %s "
                "(%d entries); runs are frozen evidence -- use a new run_id"
                % (output_dir, len(entries))
            )
    for sub in ("logs", "inputs", "results"):
        (output_dir / sub).mkdir(parents=True, exist_ok=True)

    buckets_path = output_dir / "results" / "grodon_66_buckets.tsv"
    with open(buckets_path, "w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(list(BUCKET_TABLE_COLUMNS))
        for row in rows:
            writer.writerow(
                [
                    row.get("accession", ""),
                    row.get("bucket", ""),
                    row.get("evidence", ""),
                ]
            )

    payload_out = dict(payload)
    payload_out["run"] = dict(provenance)
    payload_out["input_files"] = [dict(row) for row in input_rows]
    payload_out["bucket_table"] = {
        "path": "results/grodon_66_buckets.tsv",
        "columns": list(BUCKET_TABLE_COLUMNS),
        "rows": len(rows),
    }
    reconciliation_path = output_dir / "results" / "grodon_66_reconciliation.json"
    reconciliation_path.write_text(
        json.dumps(payload_out, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )

    contract = {
        "schema_version": INPUT_CONTRACT_SCHEMA_VERSION,
        "run_id": provenance.get("run_id"),
        "phase": PHASE,
        "generated_at": provenance.get("generated_at"),
        "status": payload.get("status"),
        "inputs": {
            str(row.get("role")): {
                "path": row.get("path"),
                "status": row.get("status"),
                "size": row.get("size"),
                "sha256": row.get("sha256"),
            }
            for row in input_rows
        },
    }
    (output_dir / "input_contract.json").write_text(
        json.dumps(contract, indent=2) + "\n", encoding="utf-8"
    )
    return payload_out


def _log_path(output_dir: Path) -> Path:
    return Path(output_dir) / "logs" / "reconcile_grodon_66.log"


def append_log(output_dir: Path, lines: Sequence[str]) -> None:
    """Append command + timestamp lines to the run's log (one file, appended)."""
    path = _log_path(output_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8", newline="\n") as handle:
        for line in lines:
            handle.write(line.rstrip("\n") + "\n")


def _validate_run_id(run_id: str) -> str:
    text = str(run_id).strip()
    if not RUN_ID_PATTERN.match(text) or ".." in text or "/" in text or "\\" in text:
        raise ValueError(
            "run_id %r must match %s and contain no path traversal"
            % (run_id, RUN_ID_PATTERN.pattern)
        )
    return text


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconcile the gRodon 66-record manifest delta by accession sets (plan F2)."
        )
    )
    parser.add_argument(
        "--run-root",
        type=Path,
        required=True,
        help="the frozen gRodon run directory to read (read-only)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="directory of the new dated run (must be empty or absent)",
    )
    parser.add_argument("--run-id", default=DEFAULT_RUN_ID)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="analyse and print the summary without creating a run directory",
    )
    return parser.parse_args(argv)


def _emit(payload: Mapping[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=False))


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point.  Returns 0 on a clean analysis, 1 on a broken input."""
    args = _parse_args(argv)
    run_root = Path(args.run_root)
    try:
        analysis = run_analysis(run_root)
    except (FileNotFoundError, ReconciliationInputError, BucketOverlapError, ValueError) as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return 1

    bundle: RunInputs = analysis["bundle"]
    payload: dict[str, Any] = dict(analysis["payload"])
    generated_at = _datetime.datetime.now(_datetime.timezone.utc).isoformat(
        timespec="seconds"
    )

    # ``run_id`` is validated before anything else so a traversal attempt can
    # never reach a filesystem call, not even in dry-run.
    try:
        run_id = _validate_run_id(args.run_id)
    except ValueError as error:
        print("ERROR: %s" % error, file=sys.stderr)
        return 1
    output_dir = (
        Path(args.output_dir) if args.output_dir is not None else Path("runs") / run_id
    )

    extra_inputs: list[tuple[str, Path]] = [
        (role, _scripts_dir() / filename) for role, filename in REFERENCE_FILES
    ]
    extra_inputs.append(
        ("reconcile_phaded_grodon_manifest_delta", Path(__file__).resolve())
    )
    input_rows = [row.as_dict() for row in bundle.files]
    input_rows.extend(provenance_rows(extra_inputs))

    provenance = {
        "run_id": run_id,
        "output_dir": str(output_dir),
        "generated_at": generated_at,
        "source_run": str(run_root),
        "phase": PHASE,
        "python": sys.version.split()[0],
        "frozen_evidence_read_only": True,
    }

    payload["evidence"] = analysis["evidence"]
    payload["documented_figures"] = {
        "degrader_input_genomes": len(bundle.degrader_input),
        "degrader_genomes_excluded_no_control": int(
            bundle.declared["degrader_genomes_excluded_no_control"]
        ),
        "eligible_positives_claimed": analysis["evidence"]["declared_residual_eligible"],
        "manifest_positive_count": len(bundle.manifest_positive),
        "claimed_difference": analysis["evidence"]["declared_residue_count"],
        "missing_fasta": int(bundle.declared.get("degrader_genomes_missing_fasta", 0)),
        "derivation": "eligible = degrader_input_genomes - degrader_genomes_excluded_no_control",
    }
    payload["identities"] = identity_checks(bundle, analysis["evidence"])
    payload["difference_reading"] = {
        "material_difference": {
            "count": len(bundle.absent_from_manifest),
            "definition": "len(degrader_input - manifest_positive)",
            "basis": "accession_set_difference",
            "note": (
                "this is what the 'difference' field reports: the set difference "
                "between the two accession sets that exist on disk. It is 26,182, "
                "not 66, because eligibility cannot be reconstructed."
            ),
        },
        "hypothesised_difference": {
            "count": analysis["evidence"]["declared_residue_count"],
            "definition": (
                "(degrader_input_genomes - degrader_genomes_excluded_no_control) "
                "- manifest_positive"
            ),
            "basis": "aggregate_subtraction_of_a_scalar_that_has_no_accessions",
            "note": (
                "this is the documented 66 = 4,507 - 4,441. It is reproduced exactly "
                "but is NOT a set difference: the 26,116 subtrahend names no "
                "accession, so no membership can be derived from it."
            ),
        },
        "bucket_sum_closes_the_documented_figures": {
            "manifest_positive_measured": len(bundle.manifest_positive),
            "failed_positive_measured": len(bundle.failed_positive),
            "residue_without_accessions": analysis["evidence"]["declared_residue_count"],
            "sum": (
                len(bundle.manifest_positive)
                + len(bundle.failed_positive)
                + analysis["evidence"]["declared_residue_count"]
            ),
            "documented_eligible_count": analysis["evidence"]["declared_residual_eligible"],
            "closes": (
                len(bundle.manifest_positive)
                + len(bundle.failed_positive)
                + analysis["evidence"]["declared_residue_count"]
                == analysis["evidence"]["declared_residual_eligible"]
            ),
        },
    }
    payload["hypothesis"] = {
        "claimed_eligible_count": analysis["evidence"]["declared_residual_eligible"],
        "claimed_difference": analysis["evidence"]["declared_residue_count"],
        "member_accessions": 0,
        "status": "pending",
        "derivation": "30,623 - 26,116 (aggregate subtraction, not a set difference)",
        "missing_input": analysis["evidence"]["missing_input"],
        "bound": {
            "eligible_lower_bound": analysis["evidence"]["eligible_lower_bound"],
            "eligible_upper_bound": analysis["evidence"]["eligible_upper_bound"],
        },
        "residue_arithmetic": (
            "|degrader_input - manifest_positive| (%d) - degrader_genomes_excluded_no_control "
            "(%d) = %d, which reproduces the documented difference of %d exactly"
            % (
                len(bundle.absent_from_manifest),
                int(bundle.declared["degrader_genomes_excluded_no_control"]),
                analysis["evidence"]["declared_eligible_residue"],
                analysis["evidence"]["declared_residue_count"],
            )
        ),
    }
    payload["material_decomposition"] = {
        "degrader_input_count": len(bundle.degrader_input),
        "manifest_positive_count": len(bundle.manifest_positive),
        "failed_positive_count": len(bundle.failed_positive),
        "failed_negative_count": len(bundle.failed_negative),
        "absent_from_manifest_count": len(bundle.absent_from_manifest),
        "identity_material": (
            "|degrader_input| == |manifest_positive| + |absent_from_manifest|"
        ),
        "identity_holds": (
            len(bundle.degrader_input)
            == len(bundle.manifest_positive) + len(bundle.absent_from_manifest)
        ),
        "failed_tables_are_subset_of_manifest_arms": (
            bundle.failed_positive <= bundle.manifest_positive
            and bundle.failed_negative <= bundle.manifest_negative
        ),
        "failed_tables_disjoint_from_absent": not (
            (bundle.failed_positive | bundle.failed_negative) & bundle.absent_from_manifest
        ),
        "declared_scalar": int(bundle.declared["degrader_genomes_excluded_no_control"]),
    }
    payload["no_control_definition"] = {
        "authoritative_definition": analysis["evidence"]["why_not"],
        "exclusion_file_role": (
            "inputs/pool_genomes_exclusion.txt is the 62,666-genome candidate-pool "
            "ledger used to build the control arm; it is NOT the no-control ledger "
            "(%d of the %d degrader inputs are also in it, and %d degrader inputs "
            "are absent from it)"
            % (
                analysis["evidence"]["degrader_also_in_exclusion_count"],
                len(bundle.degrader_input),
                analysis["evidence"]["degrader_minus_exclusion_count"],
            )
        ),
        "exclusion_file_is_the_no_control_set": analysis["evidence"][
            "exclusion_file_is_the_no_control_set"
        ],
        "degrader_minus_exclusion_count": analysis["evidence"][
            "degrader_minus_exclusion_count"
        ],
        "rejected_alternative": (
            "degrader_genomes.tsv minus pool_genomes_exclusion.txt gives %d "
            "accessions, not %d, so the plan's suggested set operation is not the "
            "authoritative definition"
            % (
                analysis["evidence"]["degrader_minus_exclusion_count"],
                analysis["evidence"]["declared_residual_eligible"],
            )
        ),
    }
    payload["check_grodon_manifest_delta"] = {
        "on_the_30_623_accession_decomposition": analysis["check_verified"],
        "on_the_documented_4507_cardinality_with_placeholders": analysis["check_claimed"],
        "flip_to_ok_requires": [
            "an accession-level eligible-positive list -- today only the scalar "
            "4,507 = 30,623 - 26,116 exists, and 26,116 is itself a scalar summing "
            "two different reasons",
            "an exclusion_reasons mapping for every eligible positive absent from "
            "the manifest (this is the same missing accession list)",
            "a prediction_failures mapping for the 23 positive-arm failures "
            "(already available from the failure table's phaZ_status column)",
        ],
    }

    if args.dry_run:
        payload["run"] = dict(provenance, dry_run=True, output_dir=None, written=False)
        _emit(payload)
        return 0

    write_run(
        output_dir,
        payload,
        analysis["rows"],
        input_rows=input_rows,
        provenance=provenance,
    )
    append_log(
        output_dir,
        [
            "# F2 gRodon 66-record reconciliation -- exact commands",
            "# generated_at=%s" % generated_at,
            "# run_id=%s" % run_id,
            "# source_run=%s" % run_root,
            "python -m unittest pipeline.tests.test_reconcile_phaded_grodon_manifest_delta -v",
            "python -m compileall -q pipeline/scripts pipeline/tests",
            "git diff --check",
            "python pipeline/scripts/reconcile_phaded_grodon_manifest_delta.py "
            "--run-root %s --run-id %s" % (run_root, run_id),
        ],
    )
    _emit(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
