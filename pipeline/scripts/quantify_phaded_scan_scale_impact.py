#!/usr/bin/env python3
"""Quantify how the corrected scan-13 E-value scale hits the frozen v1 universe.

The question this module answers
--------------------------------
Task F7 rescaled the frozen scan-13 ``hits_all.tsv`` off the per-shard HMMER
``Z`` and onto one measured full-library ``Z`` (``615,969,589``), producing
``hits_all_rescaled.tsv`` — the frozen 6,743,197 rows plus ``E-value_shard_z``,
``E-value_rescaled``, ``z_shard``, ``z_total`` and
``passes_rescaled_threshold``.  The headline was:

    rows_passing_original_threshold = 6,742,621
    rows_passing_rescaled_threshold = 6,305,428
    delta_rows                      = -437,193   (newly rejected)

The frozen v1 candidate universe — the **109,087** proteins of
``runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv``
— was built from those very hits, at the *per-shard* scale.  So for every
universe ``accession`` this module asks:

1. does it keep **at least one** row with ``passes_rescaled_threshold`` true?
2. how many accessions lose **all** support (``rescaled_support_rows == 0``)
   and were therefore supported only by newly rejected rows?
3. who are they, decomposed by ``family``, ``source_layer`` and ``layer``, and
   cross-tabulated with the universe's own ``layer`` column?
4. of the 437,193 newly rejected rows, how many were the **only** support for
   some universe accession (consequential) and how many were redundant?
5. how many universe accessions appear in **no** hit row at all (a
   join-completeness check that must be 0 or explicitly explained — a prior
   frozen measurement put the pool inside the old hit table at
   109,034 / 109,087, i.e. 53 accessions unaccounted for).

The join key is measured, never assumed
---------------------------------------
The universe carries ``accession``; the rescaled table carries both ``protein``
and ``tacc``.  ``accession == protein`` is *expected*, so this module counts, in
one streaming pass, how many universe accessions are found in each column and in
neither, and reports the counts per column in ``join.per_key``.  A row whose
``protein`` cell and ``tacc`` cell name **different** universe accessions is not
resolved by preference: the row is recorded and the run fails closed with
``JoinKeyAmbiguityError`` unless one column is pinned with ``--join-key``.

Bounded memory on a 985 MB table
--------------------------------
``hits_all_rescaled.tsv`` is **985,026,355 bytes** (6,743,197 data rows) on the
server and must never be loaded.  The table is read once, one **size-bounded**
line at a time (``readline(MAX_LINE_LENGTH)``; an over-long line is refused, not
buffered), and the only retained state is O(U) in the universe size U: a hash
map with one small slot record per universe accession plus one aggregate
counter per cross-tab bucket.  Rows that claim no universe accession are counted
and dropped immediately, which is what makes the state independent of the 6.7 M
row count — the table can grow without the process growing.
``memory_bound_bytes`` states the bound and every run reports it as
``bounded_state``.

Why a hash join instead of an external sort/merge: a merge join would have to
spill and sort 985 MB of hit rows purely to shed state that a hash map already
bounds at O(U) — the universe is 109,087 entries, three orders of magnitude
smaller than the table — and it would add a large temp-file failure mode on a
server whose ``runs/`` tree is already 1 TB.  The state that remains is
proportional to the universe, not to the hits table.

Threshold semantics
-------------------
``passes_rescaled_threshold`` is read verbatim; pass/fail is **never** inferred
from ``E-value_rescaled`` while the column exists (the pre-flight fails closed
and names the observed values otherwise, and the only fallback is the explicit
``--infer-pass-from-rescaled``, recorded as
``scale_column.inferred_from_rescaled = true``).  The frozen per-shard state of
a row is recomputed from the preserved ``E-value`` cell with the project rule
``E < threshold`` (strict) at ``--threshold`` (default ``1e-5``, the frozen
``hmmsearch -E``), because F7's table carries no frozen-pass column.  A row is
**newly rejected** when it passes the frozen scale and fails the corrected one,
and the accounting is exact:

    rows_newly_rejected = rows_newly_rejected_consequential
                          + rows_newly_rejected_redundant

A newly rejected row is *consequential* when its accession keeps no passing row
at all, and *redundant* when that accession keeps at least one.

What is written (only ever inside ``--output-dir``)
---------------------------------------------------
* ``accession_scale_impact.tsv`` — one row per universe accession:
  ``accession, genome, family, source_layer, layer, total_hit_rows,
  rescaled_support_rows, frozen_support_rows, lost_all_support,
  rejection_consequential``;
* ``scale_impact_summary.json`` — every count, the cross-tabs, the join-key
  evidence, the frozen-delta cross-check and the bounded-state bound;
* ``accession_scale_impact_affected.tsv`` — the affected accession list
  (``lost_all_support == true``) with the same columns.

Failure modes (all fail closed)
-------------------------------
* missing/unreadable input, a missing required column, an unparsable frozen
  E-value → usage error, nothing written;
* a ``passes_rescaled_threshold`` cell outside the accepted true/false
  spellings, or a missing pass column without ``--infer-pass-from-rescaled``
  → inconsistent, nothing written, observed values named;
* differing universe accessions claimed by ``protein`` vs ``tacc`` → join-key
  ambiguity, nothing written unless one column is pinned;
* a universe size other than ``--expected-universe-size`` (default 109,087) or
  a duplicate accession → inconsistent, nothing written;
* a supplied ``--delta-report`` whose row counts disagree with the recomputed
  ones → inconsistent, nothing written;
* accessions in **no** hit row: the evidence tables are written (they are the
  evidence for the claim) and the exit code is non-zero unless the count is
  explained with ``--expected-no-hit N`` or downgraded with
  ``--allow-no-hit-accessions``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Sequence, Tuple

#: The frozen v1 universe columns that this module needs.
UNIVERSE_REQUIRED_COLUMNS = [
    "accession",
    "genome",
    "family",
    "source_layer",
    "layer",
]

#: The join-key candidate columns of the F7 rescaled table.
KEY_COLUMNS = ("protein", "tacc")

#: Column names produced by F7's ``rescale_phaded_scan_hits.py``.
DEFAULT_PASS_COLUMN = "passes_rescaled_threshold"
RESCALED_EVALUE_COLUMN = "E-value_rescaled"
FROZEN_EVALUE_COLUMN = "E-value_shard_z"

#: Spellings accepted for the frozen pass column.  ``pending`` (a row whose
#: shard had no scale) is deliberately **not** accepted: a pending row is not a
#: decision, and counting it either way would silently change the answer.
TRUE_VALUES = frozenset({"true"})
FALSE_VALUES = frozenset({"false"})

#: The project threshold rule (see the docstring of rescale_phaded_scan_hits.py).
THRESHOLD_RULE = "E < threshold (strict)"

#: The frozen v1 universe size, measured (109,087 unique accessions).
DEFAULT_EXPECTED_UNIVERSE_SIZE = 109_087

#: Rows scanned before any output exists, to catch an unusable scale column
#: early instead of after a multi-minute pass over 985 MB.
DEFAULT_PREFLIGHT_ROWS = 100_000

#: Hard cap on the number of join-key slots retained.  The universe is 109,087
#: accessions, so anything far above that means the "universe" is not the
#: universe and the run must stop instead of growing without bound.
MAX_JOIN_KEYS = 2_000_000

#: Hard cap on one line of either TSV.  The widest frozen hit row is ~200 bytes;
#: 1 MiB cannot appear in a well-formed table and refusing it keeps the read
#: bound honest (no silently materialised multi-megabyte "line").
MAX_LINE_LENGTH = 1024 * 1024

#: Per-accession state model, in bytes: an interned accession string plus the
#: slot record (see ``AccessionInfo``/``AccessionTally``) and the dict entry.
#: Used only to *state* the bound; a test measures the real peak.
BYTES_PER_UNIVERSE_ACCESSION = 272
FIXED_STATE_BYTES = 3 * 1024 * 1024

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_INCONSISTENT = 3
EXIT_REFUSED = 4

AUTO_JOIN_KEY = "auto"

SCALE_IMPACT_COLUMNS = [
    "accession",
    "genome",
    "family",
    "source_layer",
    "layer",
    "total_hit_rows",
    "rescaled_support_rows",
    "frozen_support_rows",
    "lost_all_support",
    "rejection_consequential",
]

IMPACT_FILE = "accession_scale_impact.tsv"
SUMMARY_FILE = "scale_impact_summary.json"
AFFECTED_FILE = "accession_scale_impact_affected.tsv"
KNOWN_OUTPUT_FILES = (IMPACT_FILE, SUMMARY_FILE, AFFECTED_FILE)

#: Buckets every cross-tab entry carries.
CROSS_TAB_FIELDS = (
    "universe_accessions",
    "accessions_with_hits",
    "accessions_in_no_hit_row",
    "accessions_retaining_rescaled_support",
    "accessions_lost_all_support",
    "accessions_lost_all_support_consequential",
    "rows_joining_universe",
    "rows_passing_frozen",
    "rows_passing_rescaled",
    "rows_newly_rejected",
    "rows_newly_rejected_consequential",
    "rows_newly_rejected_redundant",
)

#: Frozen evidence trees and deployment trees this module must never write into.
PROTECTED_ROOTS = ("runs", "results", "deploy")


class JoinKeyAmbiguityError(ValueError):
    """Raised when two key columns claim different universe accessions."""

    def __init__(self, conflicts: Mapping[str, Sequence[str]], limit: int = 10) -> None:
        self.conflicts = {key: tuple(values) for key, values in conflicts.items()}
        shown = sorted(self.conflicts)[:limit]
        detail = "; ".join(
            f"row {key}: {', '.join(self.conflicts[key][:5])}" for key in shown
        )
        if not self.conflicts:
            super().__init__(
                "ambiguous join key: no universe accession is claimed by any key "
                f"column ({', '.join(KEY_COLUMNS)}); nothing can be joined"
            )
            return
        super().__init__(
            f"ambiguous join key: {len(self.conflicts)} row(s) whose "
            f"{' and '.join(KEY_COLUMNS)} cells name different universe "
            f"accessions ({detail}). The join key must be measured, not "
            "preferred: re-check the table or pin one column with "
            "--join-key protein (or tacc)."
        )


class ScaleColumnError(ValueError):
    """Raised when the frozen pass column is missing or holds foreign values."""


class ScaleImpactInconsistentError(ValueError):
    """Raised when two measurements of the same quantity disagree."""


def _open_text(path: Path, mode: str):
    """Open a TSV for reading/writing (``.gz`` supported on the read side)."""
    if mode.startswith("r") and path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def memory_bound_bytes(universe_size: int) -> int:
    """The stated memory bound: O(U) state plus a fixed interpreter allowance."""
    if universe_size < 0:
        raise ValueError(f"universe_size must be non-negative, got {universe_size!r}")
    return FIXED_STATE_BYTES + BYTES_PER_UNIVERSE_ACCESSION * int(universe_size)


# --------------------------------------------------------------------------
# streaming TSV reading
# --------------------------------------------------------------------------
def iter_tsv_rows(
    path: Path | str,
    *,
    required_columns: Sequence[str],
    limit: int | None = None,
    max_line_length: int = MAX_LINE_LENGTH,
) -> Iterator[Dict[str, str]]:
    """Yield one dict per data row, streaming and never materialising the file.

    Reads are explicitly bounded (``readline(max_line_length)``) rather than
    ``for line in handle``: the frozen table is 985 MB and the bound has to be
    observable at the handle — a test can wrap the stream and refuse an unsized
    read — not merely intended.  A line that fills the bound is refused instead
    of being concatenated into an unbounded buffer.

    The header is validated against ``required_columns`` before any row is
    yielded, so a wrong table fails before any output is produced.  Extra
    columns are tolerated (the frozen table gained five after F7).  ``limit``
    caps the number of data rows yielded.
    """
    path = Path(path)
    if max_line_length < 1:
        raise ValueError(f"max_line_length must be positive, got {max_line_length!r}")
    handle = _open_text(path, "rt")
    try:
        raw_header = handle.readline(max_line_length)
        if not raw_header:
            raise ValueError(f"{path}: the table is empty")
        if not raw_header.endswith("\n"):
            raise ValueError(
                f"{path} line 1: more than {max_line_length} bytes; refusing to "
                "buffer it"
            )
        header = [cell.strip() for cell in raw_header.rstrip("\r\n").split("\t")]
        missing = [column for column in required_columns if column not in header]
        if missing:
            raise ValueError(
                f"{path}: missing required column(s): {', '.join(missing)} "
                f"(found: {', '.join(header)})"
            )
        index = {column: header.index(column) for column in required_columns}
        width = len(header)
        seen = 0
        number = 1
        while True:
            line = handle.readline(max_line_length)
            if not line:
                return
            number += 1
            if not line.endswith("\n"):
                raise ValueError(
                    f"{path} line {number}: more than {max_line_length} bytes or a "
                    "final line without a terminator; refusing to buffer it"
                )
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) < width:
                raise ValueError(
                    f"{path} line {number}: expected {width} columns, got "
                    f"{len(cells)}"
                )
            yield {column: cells[position] for column, position in index.items()}
            seen += 1
            if limit is not None and seen >= limit:
                return
    finally:
        handle.close()


def sha256_file(path: Path | str, chunk_bytes: int = 1024 * 1024) -> str:
    """Return the SHA-256 of ``path``, streamed in fixed-size blocks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# value parsing
# --------------------------------------------------------------------------
def pass_value(cell: Any) -> bool:
    """Return the boolean of a ``passes_rescaled_threshold`` cell.

    Only the exact spellings produced by F7 are accepted (``true``/``false``,
    case-insensitive).  ``pending``, ``1``/``0`` and an empty cell all raise:
    guessing one of them would silently move rows between "kept" and "lost".
    """
    text = ("" if cell is None else str(cell)).strip()
    lowered = text.lower()
    if lowered in TRUE_VALUES:
        return True
    if lowered in FALSE_VALUES:
        return False
    raise ValueError(
        f"{DEFAULT_PASS_COLUMN} has the unexpected value {cell!r}; accepted "
        "values are 'true' and 'false' (a 'pending' row means its shard had no "
        "scale and cannot be counted either way)"
    )


def parse_evalue(cell: Any, field: str) -> float:
    """Parse a frozen E-value cell into a finite, non-negative float."""
    text = ("" if cell is None else str(cell)).strip()
    if not text:
        raise ValueError(f"{field} is empty; no E-value was recorded")
    try:
        number = float(text)
    except ValueError as error:
        raise ValueError(f"{field} is not a number: {cell!r}") from error
    if math.isnan(number) or math.isinf(number):
        raise ValueError(f"{field} must be finite, got {cell!r}")
    if number < 0:
        raise ValueError(f"{field} must be non-negative, got {cell!r}")
    return number


def parse_positive_float(cell: Any, field: str) -> float:
    """Parse a strictly positive float."""
    number = parse_evalue(cell, field)
    if number <= 0:
        raise ValueError(f"{field} must be positive, got {cell!r}")
    return number


def format_threshold(value: float) -> str:
    """Render a threshold the way the frozen reports do (``%.6g``)."""
    return "%.6g" % float(value)


def non_negative_int(value: Any, field: str) -> int:
    """Parse a non-negative integer, rejecting bools and junk."""
    if isinstance(value, bool):
        raise ValueError(f"{field} must be an integer, got {value!r}")
    if isinstance(value, int):
        number = value
    else:
        text = ("" if value is None else str(value)).strip()
        if not text or not text.lstrip("+").isdigit():
            raise ValueError(f"{field} must be an integer, got {value!r}")
        number = int(text, 10)
    if number < 0:
        raise ValueError(f"{field} must be non-negative, got {value!r}")
    return number


# --------------------------------------------------------------------------
# per-accession state
# --------------------------------------------------------------------------
class AccessionInfo:
    """The immutable per-accession coordinates copied from the universe."""

    __slots__ = ("genome", "family", "source_layer", "layer")

    def __init__(self, genome: str, family: str, source_layer: str, layer: str) -> None:
        self.genome = genome
        self.family = family
        self.source_layer = source_layer
        self.layer = layer


class AccessionTally:
    """The mutable per-accession counters; O(1) per accession, O(U) overall."""

    __slots__ = (
        "total_hit_rows",
        "rescaled_support_rows",
        "frozen_support_rows",
        "newly_rejected_rows",
    )

    def __init__(self) -> None:
        self.total_hit_rows = 0
        self.rescaled_support_rows = 0
        self.frozen_support_rows = 0
        self.newly_rejected_rows = 0


class KeyStats:
    """Per-key-column join evidence: how many accessions and rows it resolves."""

    __slots__ = ("found", "rows_matched", "users")

    def __init__(self) -> None:
        self.found = 0
        self.rows_matched = 0
        self.users: set = set()


def read_universe(path: Path | str) -> Dict[str, AccessionInfo]:
    """Read ``protein_layers.tsv`` into an accession-keyed map of coordinates.

    A repeated accession is refused: the universe is documented as 109,087
    unique accessions, so a duplicate means the wrong (or corrupted) table, and
    a silent last-wins would quietly change the join.
    """
    path = Path(path)
    universe: Dict[str, AccessionInfo] = {}
    for number, row in enumerate(
        iter_tsv_rows(path, required_columns=UNIVERSE_REQUIRED_COLUMNS), start=2
    ):
        accession = (row.get("accession") or "").strip()
        if not accession:
            raise ValueError(f"{path} line {number}: empty accession")
        if accession in universe:
            raise ValueError(
                f"{path} line {number}: duplicate accession {accession!r}; the "
                "universe must have exactly one row per accession"
            )
        universe[accession] = AccessionInfo(
            genome=(row.get("genome") or "").strip(),
            family=(row.get("family") or "").strip(),
            source_layer=(row.get("source_layer") or "").strip(),
            layer=(row.get("layer") or "").strip(),
        )
    if not universe:
        raise ValueError(f"{path}: the universe has no data rows")
    return universe


class _CrossTab:
    """A per-dimension-value rollup of the accession outcomes."""

    def __init__(self, dimension: str) -> None:
        self.dimension = dimension
        self.buckets: Dict[str, Dict[str, int]] = {}

    def bucket(self, value: str) -> Dict[str, int]:
        key = value if value else "<empty>"
        entry = self.buckets.get(key)
        if entry is None:
            entry = {field: 0 for field in CROSS_TAB_FIELDS}
            self.buckets[key] = entry
        return entry

    def as_dict(self) -> Dict[str, Dict[str, int]]:
        return {key: dict(self.buckets[key]) for key in sorted(self.buckets)}


class ScaleImpactAccumulator:
    """The single streaming pass: one bounded update per hit row."""

    def __init__(
        self,
        universe: Mapping[str, AccessionInfo],
        *,
        join_key: str,
        pass_column: str = DEFAULT_PASS_COLUMN,
        infer_pass: bool = False,
        threshold: float | None = None,
        max_join_keys: int = MAX_JOIN_KEYS,
    ) -> None:
        self.universe = universe
        self.join_key = join_key
        self.pass_column = pass_column
        self.infer_pass = bool(infer_pass)
        self.threshold = threshold
        self.max_join_keys = int(max_join_keys)
        self.tallies: Dict[str, AccessionTally] = {}
        self.key_stats: Dict[str, KeyStats] = {
            key: KeyStats() for key in KEY_COLUMNS
        }
        self.conflicts: Dict[str, Tuple[str, ...]] = {}
        self.pass_values: Dict[str, int] = {}
        # row-level accounting
        self.rows_total = 0
        self.rows_joining_universe = 0
        self.rows_not_in_universe = 0
        self.rows_passing_rescaled = 0
        self.rows_failing_rescaled = 0
        self.rows_passing_frozen = 0
        self.rows_failing_frozen = 0
        self.rows_newly_rejected = 0
        # cross-tabs
        self.tabs = {
            "family": _CrossTab("family"),
            "source_layer": _CrossTab("source_layer"),
            "layer": _CrossTab("layer"),
            "layer_x_family": _CrossTab("layer_x_family"),
        }

    # -- streaming -------------------------------------------------------
    def observe(
        self,
        row: Mapping[str, str],
        *,
        row_number: int,
        path: Path | str,
        on_conflict: str = "raise",
    ) -> None:
        """Fold one hit row into the per-accession state.

        ``on_conflict`` is ``"raise"`` for a pre-flight (stop at the first
        ambiguity) or ``"collect"`` for the full pass (record every conflicted
        row so the failure names all of them).
        """
        self.rows_total += 1
        try:
            passes = self._passes(row)
        except ValueError as error:
            raise ScaleColumnError(f"{path} line {row_number}: {error}") from error

        if self.infer_pass:
            frozen = passes
        else:
            frozen_cell = row.get(FROZEN_EVALUE_COLUMN)
            if frozen_cell is None:
                frozen_cell = row.get("E-value")
            frozen = (
                parse_evalue(
                    frozen_cell,
                    f"{path} line {row_number}: frozen E-value "
                    f"({FROZEN_EVALUE_COLUMN})",
                )
                < self.threshold
            )

        claims: Dict[str, str | None] = {}
        for key in KEY_COLUMNS:
            cell = (row.get(key) or "").strip()
            claims[key] = cell if cell in self.universe else None

        if self.join_key == AUTO_JOIN_KEY:
            effective = next((key for key in KEY_COLUMNS if claims[key]), None)
        else:
            effective = self.join_key
        primary: List[str] = []
        if effective is not None and claims.get(effective):
            primary = [str(claims[effective])]

        if self._claims_conflict(claims, primary):
            self.conflicts[str(row_number)] = tuple(
                sorted({value for value in claims.values() if value})
            )
            if on_conflict == "raise":
                raise JoinKeyAmbiguityError(self.conflicts)

        if effective is not None and primary:
            stats = self.key_stats[effective]
            if primary[0] not in stats.users:
                stats.users.add(primary[0])
                stats.found += 1
            stats.rows_matched += 1

        if not primary:
            self.rows_not_in_universe += 1
            return
        self.rows_joining_universe += 1

        # Row-level counters are only advanced for rows that actually claim a
        # universe accession: the hit table is far wider than the universe
        # (6,743,197 rows vs 109,087 accessions), and its counts must reconcile
        # with the per-accession sums.
        if passes:
            self.rows_passing_rescaled += 1
        else:
            self.rows_failing_rescaled += 1
        if frozen:
            self.rows_passing_frozen += 1
        else:
            self.rows_failing_frozen += 1
        newly_rejected = bool(frozen and not passes)
        if newly_rejected:
            self.rows_newly_rejected += 1

        for accession in primary:
            tally = self.tally(accession)
            tally.total_hit_rows += 1
            if passes:
                tally.rescaled_support_rows += 1
            if frozen:
                tally.frozen_support_rows += 1
            if newly_rejected:
                tally.newly_rejected_rows += 1

    def _passes(self, row: Mapping[str, str]) -> bool:
        """Return the row's rescaled pass state, from the column or by rule."""
        cell = row.get(self.pass_column)
        if cell is None:
            if self.infer_pass:
                return (
                    parse_evalue(
                        row.get(RESCALED_EVALUE_COLUMN), RESCALED_EVALUE_COLUMN
                    )
                    < self.threshold
                )
            raise ScaleColumnError(
                f"{self.pass_column} is absent from the table and pass/fail is "
                f"never inferred from {RESCALED_EVALUE_COLUMN} without "
                "--infer-pass-from-rescaled"
            )
        raw = str(cell).strip()
        self.pass_values[raw] = self.pass_values.get(raw, 0) + 1
        return pass_value(raw)

    def _claims_conflict(
        self, claims: Mapping[str, str | None], primary: Sequence[str]
    ) -> bool:
        """True when two key columns name different universe accessions.

        A column that names nothing (``tacc`` is ``-`` for most frozen rows)
        never conflicts; only a cell that resolves to a *different* accession
        does.  A pinned ``--join-key`` is an explicit decision, so the other
        column is then only reported, not treated as a conflict.
        """
        if not primary or self.join_key != AUTO_JOIN_KEY:
            return False
        for key in KEY_COLUMNS:
            value = claims.get(key)
            if value is not None and value not in primary:
                return True
        return False

    # -- state -----------------------------------------------------------
    def tally(self, accession: str) -> AccessionTally:
        tally = self.tallies.get(accession)
        if tally is None:
            if len(self.tallies) >= self.max_join_keys:
                raise ScaleColumnError(
                    f"more than {self.max_join_keys} universe accessions have "
                    "been claimed by hit rows; the universe is far smaller than "
                    "that, so the join is wrong and the state is not bounded"
                )
            tally = AccessionTally()
            self.tallies[accession] = tally
        return tally

    def observed_pass_values(self) -> List[str]:
        return sorted(self.pass_values)

    def columns_agree(self) -> str:
        """How the key columns relate: proven, one-sided, disjoint or overlapping."""
        protein = self.key_stats["protein"].users
        tacc = self.key_stats["tacc"].users
        if not protein and not tacc:
            return "neither"
        if not tacc:
            return "protein_only"
        if not protein:
            return "tacc_only"
        if protein == tacc:
            return "identical"
        return "overlap" if protein & tacc else "disjoint"


def classify_accessions(
    universe: Mapping[str, AccessionInfo], accumulator: ScaleImpactAccumulator
) -> Dict[str, Dict[str, Any]]:
    """Turn the tallies into one outcome record per universe accession."""
    outcomes: Dict[str, Dict[str, Any]] = {}
    for accession, info in universe.items():
        tally = accumulator.tallies.get(accession)
        total = tally.total_hit_rows if tally else 0
        rescaled = tally.rescaled_support_rows if tally else 0
        frozen = tally.frozen_support_rows if tally else 0
        rejected = tally.newly_rejected_rows if tally else 0
        # ``lost_all_support`` needs a row to lose: an accession that never
        # appears on this table is reported separately as a join-completeness
        # failure, not as a scale casualty.
        lost = bool(total > 0 and rescaled == 0)
        outcomes[accession] = {
            "accession": accession,
            "genome": info.genome,
            "family": info.family,
            "source_layer": info.source_layer,
            "layer": info.layer,
            "total_hit_rows": total,
            "rescaled_support_rows": rescaled,
            "frozen_support_rows": frozen,
            "newly_rejected_rows": rejected,
            "lost_all_support": lost,
            "rejection_consequential": bool(lost and rejected > 0),
        }
    return outcomes


# --------------------------------------------------------------------------
# argument surface
# --------------------------------------------------------------------------
def build_arg_namespace(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Quantify how the corrected (full-library Z) scan-13 E-value scale "
            "affects the frozen v1 candidate universe in protein_layers.tsv. "
            "Streams the 985 MB rescaled hits table; writes only new files."
        )
    )
    parser.add_argument(
        "--universe",
        type=Path,
        required=True,
        help="frozen v1 universe (protein_layers.tsv)",
    )
    parser.add_argument(
        "--hits",
        type=Path,
        required=True,
        help="frozen F7 hits_all_rescaled.tsv (985 MB on the server; streamed)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="NEW directory for this run's outputs (never runs/, results/ or deploy/)",
    )
    parser.add_argument(
        "--affected-out",
        type=Path,
        default=None,
        help=f"affected-accession TSV name inside --output-dir (default {AFFECTED_FILE})",
    )
    parser.add_argument(
        "--delta-report",
        type=Path,
        default=None,
        help="optional F7 scale_delta_report.json to cross-check the row counts",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1e-5,
        help="frozen reporting threshold used to recompute the per-shard pass "
        "state of the preserved E-value column (default 1e-05)",
    )
    parser.add_argument(
        "--pass-column",
        default=DEFAULT_PASS_COLUMN,
        help=f"the rescaled pass/fail column (default {DEFAULT_PASS_COLUMN})",
    )
    parser.add_argument(
        "--infer-pass-from-rescaled",
        action="store_true",
        help="explicitly allow deriving pass/fail from E-value_rescaled when the "
        "pass column is absent; recorded in the output",
    )
    parser.add_argument(
        "--join-key",
        choices=(AUTO_JOIN_KEY, *KEY_COLUMNS),
        default=AUTO_JOIN_KEY,
        help="pin the join column; 'auto' resolves it from the measured matches "
        "and fails closed when the two columns disagree",
    )
    parser.add_argument(
        "--expected-universe-size",
        type=int,
        default=DEFAULT_EXPECTED_UNIVERSE_SIZE,
        help="fail closed unless the universe has exactly this many accessions "
        f"(default {DEFAULT_EXPECTED_UNIVERSE_SIZE}; 0 disables the check)",
    )
    parser.add_argument(
        "--expected-no-hit",
        type=int,
        default=None,
        help="accept exactly this many universe accessions being absent from "
        "every hit row; without it that count is an unexplained inconsistency",
    )
    parser.add_argument(
        "--allow-no-hit-accessions",
        action="store_true",
        help="downgrade the no-hit join-completeness failure to a warning "
        "(the count is still reported)",
    )
    parser.add_argument(
        "--preflight-rows",
        type=int,
        default=DEFAULT_PREFLIGHT_ROWS,
        help="rows scanned before writing anything, to catch an unusable scale "
        f"column early (default {DEFAULT_PREFLIGHT_ROWS})",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="read at most N hit rows (smoke run); reported as limited=true",
    )
    parser.add_argument(
        "--sha256",
        action="store_true",
        help="record the SHA-256 of both inputs (a second full read of the hits)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace this run's own output files if they already exist",
    )
    parser.add_argument("--quiet", action="store_true", help="suppress the stderr line")
    return parser.parse_args(argv)


def parse_args_never_exits(argv: Sequence[str] | None = None):
    """Parse ``argv``; return ``(status, namespace)`` instead of exiting."""
    try:
        return EXIT_OK, build_arg_namespace(argv)
    except SystemExit as exit_error:
        code = exit_error.code
        return (code if isinstance(code, int) else EXIT_USAGE), None


# --------------------------------------------------------------------------
# safety checks
# --------------------------------------------------------------------------
def _validate_input(path: Path, field: str) -> Path:
    if not path.exists():
        raise ValueError(f"{field} does not exist: {path}")
    if path.is_dir():
        raise ValueError(f"{field} is a directory, not a file: {path}")
    if not path.is_file():
        raise ValueError(f"{field} is not a regular file: {path}")
    return path


def _same_file(left: Path, right: Path) -> bool:
    try:
        if left.resolve() == right.resolve():
            return True
        return left.exists() and right.exists() and left.samefile(right)
    except OSError:  # pragma: no cover - defensive
        return False


def _check_output_dir(output_dir: Path, inputs: Sequence[Path], force: bool) -> None:
    """Refuse an unsafe output directory before anything is read or written."""
    resolved = output_dir.resolve()
    for name in PROTECTED_ROOTS:
        if resolved.name == name or name in resolved.parts:
            raise FileExistsError(
                f"--output-dir {output_dir} is inside the frozen evidence tree "
                f"{name!r}; runs/, results/ and deploy/ are read-only"
            )
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"--output-dir is not a directory: {output_dir}")
    for path in inputs:
        if _same_file(output_dir, path):
            raise FileExistsError(
                f"--output-dir resolves to the input {path}; refusing to write "
                "over an input"
            )
        if path.exists() and path.resolve().parent == resolved:
            raise FileExistsError(
                f"--output-dir {output_dir} is the directory holding the input "
                f"{path}; refusing to create outputs beside frozen inputs"
            )
    if output_dir.exists():
        entries = sorted(entry.name for entry in output_dir.iterdir())
        unknown = [name for name in entries if name not in KNOWN_OUTPUT_FILES]
        if unknown:
            raise FileExistsError(
                "--output-dir is not empty and holds unrelated entries: "
                f"{', '.join(unknown)}. Pass a NEW empty directory (--force only "
                "replaces this run's own output files)."
            )
        if entries and not force:
            raise FileExistsError(
                f"--output-dir already holds {', '.join(entries)}; pass --force "
                "to replace this run's own output files."
            )


def _write_text_atomically(path: Path, text: str, force: bool) -> None:
    """Publish ``text`` at ``path`` via a sibling temp file, never partially."""
    if path.exists() and not force:
        raise FileExistsError(f"{path} already exists; pass --force to replace it")
    temporary = path.with_name(path.name + ".partial")
    if temporary.exists():
        raise FileExistsError(
            f"a partial output is already present: {temporary}; remove it first"
        )
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


# --------------------------------------------------------------------------
# the run
# --------------------------------------------------------------------------
def read_tsv_header(path: Path | str, *, max_line_length: int = MAX_LINE_LENGTH) -> List[str]:
    """Return the header cells of a TSV, reading one bounded line."""
    path = Path(path)
    with _open_text(path, "rt") as handle:
        raw = handle.readline(max_line_length)
    if not raw:
        raise ValueError(f"{path}: the table is empty")
    if not raw.endswith("\n"):
        raise ValueError(
            f"{path} line 1: more than {max_line_length} bytes; refusing to buffer it"
        )
    return [cell.strip() for cell in raw.rstrip("\r\n").split("\t")]


def _preflight_scale_column(
    hits: Path, *, pass_column: str, rows: int, infer_pass: bool = False
) -> Dict[str, Any]:
    """Read the first ``rows`` rows and refuse an unusable scale column.

    Runs before any output exists: a table whose pass column is missing, or
    holds a value outside ``true``/``false``, must never reach the accumulator,
    where a guess could be mistaken for a measurement.
    """
    header = read_tsv_header(hits)
    if pass_column not in header and not infer_pass:
        raise ScaleColumnError(
            f"{hits}: the pass column {pass_column!r} is absent (found: "
            f"{', '.join(header)}). Pass/fail is never inferred from "
            f"{RESCALED_EVALUE_COLUMN} without --infer-pass-from-rescaled, so "
            "this run refuses rather than guessing."
        )
    for required in (
        (FROZEN_EVALUE_COLUMN, *KEY_COLUMNS)
        if infer_pass
        else (pass_column, FROZEN_EVALUE_COLUMN, *KEY_COLUMNS)
    ):
        if required not in header:
            raise ValueError(
                f"{hits}: missing required column(s): {required} (found: "
                f"{', '.join(header)})"
            )
    values: Dict[str, int] = {}
    checked = 0
    column = RESCALED_EVALUE_COLUMN if infer_pass else pass_column
    required = (
        [RESCALED_EVALUE_COLUMN, FROZEN_EVALUE_COLUMN]
        if infer_pass
        else [pass_column, FROZEN_EVALUE_COLUMN]
    )
    for row in iter_tsv_rows(hits, required_columns=required, limit=rows):
        checked += 1
        raw = str(row.get(column) or "").strip()
        values[raw] = values.get(raw, 0) + 1
        try:
            if infer_pass:
                parse_evalue(raw, RESCALED_EVALUE_COLUMN)
            else:
                pass_value(raw)
        except ValueError as error:
            if infer_pass:
                raise ValueError(
                    f"{hits}: {error} (observed so far: "
                    f"{', '.join(repr(value) for value in sorted(values))})"
                ) from error
            raise ScaleColumnError(
                f"{hits}: {error} (observed so far: "
                f"{', '.join(repr(value) for value in sorted(values))})"
            ) from error
    if checked == 0:
        raise ValueError(f"{hits}: the rescaled hits table has no data rows")
    return {"rows_checked": checked, "observed_values": sorted(values)}


def _cross_check_delta_report(
    path: Path, summary: Mapping[str, Any], *, limited: bool
) -> Dict[str, Any]:
    """Compare the recomputed row counts against F7's frozen delta report."""
    if limited:
        return {
            "status": "skipped_limited",
            "path": str(path),
            "note": "a --limit run cannot reproduce the whole-table counts",
            "mismatches": None,
        }
    report = json.loads(Path(path).read_text(encoding="utf-8"))
    recomputed_delta = summary["rows_passing_rescaled"] - summary["rows_passing_frozen"]
    pairs = [
        ("total_rows", summary["rows_total"]),
        ("rows_passing_original_threshold", summary["rows_passing_frozen"]),
        ("rows_passing_rescaled_threshold", summary["rows_passing_rescaled"]),
        ("delta_rows", recomputed_delta),
    ]
    mismatches: Dict[str, Dict[str, Any]] = {}
    for key, recomputed in pairs:
        recorded = report.get(key)
        if recorded != recomputed:
            mismatches[key] = {"report": recorded, "recomputed": recomputed}
    return {
        "status": "mismatch" if mismatches else "match",
        "path": str(path),
        "mismatches": mismatches or None,
        "report_total_rows": report.get("total_rows"),
        "report_rows_passing_original_threshold": report.get(
            "rows_passing_original_threshold"
        ),
        "report_rows_passing_rescaled_threshold": report.get(
            "rows_passing_rescaled_threshold"
        ),
        "report_delta_rows": report.get("delta_rows"),
        "recomputed_total_rows": summary["rows_total"],
        "recomputed_rows_passing_original_threshold": summary["rows_passing_frozen"],
        "recomputed_rows_passing_rescaled_threshold": summary["rows_passing_rescaled"],
        "recomputed_delta_rows": recomputed_delta,
        "threshold_rule": THRESHOLD_RULE,
    }


def _impact_tsv(
    outcomes: Mapping[str, Mapping[str, Any]], accessions: Iterable[str]
) -> str:
    lines = ["\t".join(SCALE_IMPACT_COLUMNS)]
    for accession in accessions:
        record = outcomes[accession]
        lines.append(
            "\t".join(
                [
                    record["accession"],
                    record["genome"],
                    record["family"],
                    record["source_layer"],
                    record["layer"],
                    str(record["total_hit_rows"]),
                    str(record["rescaled_support_rows"]),
                    str(record["frozen_support_rows"]),
                    "true" if record["lost_all_support"] else "false",
                    "true" if record["rejection_consequential"] else "false",
                ]
            )
        )
    return "\n".join(lines) + "\n"


def _resolve_affected_target(
    args: argparse.Namespace, output_dir: Path
) -> Path:
    """Resolve ``--affected-out`` and refuse anything outside ``--output-dir``."""
    if args.affected_out is None:
        return output_dir / AFFECTED_FILE
    candidate = Path(args.affected_out)
    target = candidate if candidate.is_absolute() else output_dir / candidate
    if target.resolve().parent != output_dir.resolve():
        raise FileExistsError(
            f"--affected-out {target} is not directly inside --output-dir "
            f"{output_dir}; outputs are never written elsewhere"
        )
    if target.name not in (IMPACT_FILE, AFFECTED_FILE):
        raise FileExistsError(
            f"--affected-out {target.name} is not one of this run's accession "
            f"table names ({IMPACT_FILE}, {AFFECTED_FILE})"
        )
    return target


def run(args: argparse.Namespace) -> Dict[str, Any]:
    """Measure the scale impact; write the outputs; return the summary.

    Raises ``ValueError`` for an unusable input, ``FileExistsError`` for an
    unsafe output target, ``ScaleColumnError`` for an unusable pass column,
    ``JoinKeyAmbiguityError`` for disagreeing key columns and
    ``ScaleImpactInconsistentError`` for an unexplained count mismatch.
    """
    universe_path = Path(args.universe)
    hits_path = Path(args.hits)
    output_dir = Path(args.output_dir)
    delta_path = Path(args.delta_report) if args.delta_report else None

    _validate_input(universe_path, "--universe")
    _validate_input(hits_path, "--hits")
    if delta_path is not None:
        _validate_input(delta_path, "--delta-report")

    threshold = parse_positive_float(args.threshold, "--threshold")
    pass_column = str(args.pass_column or "").strip()
    if not pass_column:
        raise ValueError("--pass-column must not be empty")
    if args.infer_pass_from_rescaled and pass_column != DEFAULT_PASS_COLUMN:
        raise ValueError(
            "--infer-pass-from-rescaled derives pass/fail from "
            f"{RESCALED_EVALUE_COLUMN}; it cannot be combined with --pass-column "
            f"{pass_column!r}"
        )
    expected_size = non_negative_int(
        args.expected_universe_size, "--expected-universe-size"
    )
    preflight_rows = non_negative_int(args.preflight_rows, "--preflight-rows")
    if preflight_rows < 1:
        raise ValueError(
            f"--preflight-rows must be positive, got {args.preflight_rows!r}"
        )
    limit = None if args.limit is None else non_negative_int(args.limit, "--limit")
    if limit is not None and limit < 1:
        raise ValueError(f"--limit must be a positive integer, got {args.limit!r}")
    expected_no_hit = (
        None
        if args.expected_no_hit is None
        else non_negative_int(args.expected_no_hit, "--expected-no-hit")
    )
    join_key = str(args.join_key or AUTO_JOIN_KEY)
    if join_key not in (AUTO_JOIN_KEY, *KEY_COLUMNS):
        raise ValueError(f"unknown --join-key {join_key!r}")

    _check_output_dir(output_dir, (universe_path, hits_path), bool(args.force))
    affected_target = _resolve_affected_target(args, output_dir)

    # 1. the universe.  Read whole: 109,087 rows of coordinates *is* the O(U)
    #    state the memory bound is stated for.
    universe = read_universe(universe_path)
    if expected_size and len(universe) != expected_size:
        raise ScaleImpactInconsistentError(
            f"the universe has {len(universe)} accessions but "
            f"--expected-universe-size is {expected_size}; refusing to measure a "
            "scale impact on an unexpected candidate universe"
        )

    # 2. pre-flight the scale column, before anything is read in full
    preflight = _preflight_scale_column(
        hits_path,
        pass_column=pass_column,
        rows=min(preflight_rows, limit) if limit is not None else preflight_rows,
        infer_pass=bool(args.infer_pass_from_rescaled),
    )
    if args.infer_pass_from_rescaled:
        scale_meta: Dict[str, Any] = {
            "column": RESCALED_EVALUE_COLUMN,
            "true_value": sorted(TRUE_VALUES)[0],
            "false_value": sorted(FALSE_VALUES)[0],
            "observed_values": [],
            "inferred_from_rescaled": True,
            "inference_rule": (
                f"{RESCALED_EVALUE_COLUMN} < {format_threshold(threshold)}"
            ),
            "rows_checked": preflight["rows_checked"],
        }
    else:
        scale_meta = {
            "column": pass_column,
            "true_value": sorted(TRUE_VALUES)[0],
            "false_value": sorted(FALSE_VALUES)[0],
            "observed_values": preflight["observed_values"],
            "inferred_from_rescaled": False,
            "inference_rule": None,
            "rows_checked": preflight["rows_checked"],
        }

    accumulator = ScaleImpactAccumulator(
        universe,
        join_key=join_key,
        pass_column=pass_column,
        infer_pass=bool(args.infer_pass_from_rescaled),
        threshold=threshold,
    )

    # 3. the single streaming pass over the frozen rescaled table
    scale_columns = (
        [RESCALED_EVALUE_COLUMN] if args.infer_pass_from_rescaled else [pass_column]
    )
    for number, row in enumerate(
        iter_tsv_rows(
            hits_path,
            required_columns=[
                *KEY_COLUMNS,
                FROZEN_EVALUE_COLUMN,
                *scale_columns,
            ],
            limit=limit,
        ),
        start=2,
    ):
        accumulator.observe(
            row, row_number=number, path=hits_path, on_conflict="collect"
        )

    if accumulator.rows_total == 0:
        raise ValueError(f"{hits_path}: no hit rows were read")

    # 4. per-accession outcomes and the accounting identities
    outcomes = classify_accessions(universe, accumulator)
    if not scale_meta["inferred_from_rescaled"]:
        scale_meta["observed_values"] = accumulator.observed_pass_values()
        scale_meta["rows_checked"] = accumulator.rows_total

    total_hit_rows = sum(record["total_hit_rows"] for record in outcomes.values())
    rescaled_support_rows = sum(
        record["rescaled_support_rows"] for record in outcomes.values()
    )
    frozen_support_rows = sum(
        record["frozen_support_rows"] for record in outcomes.values()
    )
    for label, per_accession, row_level in (
        ("rows", total_hit_rows, accumulator.rows_joining_universe),
        (
            "rescaled support",
            rescaled_support_rows,
            accumulator.rows_passing_rescaled,
        ),
        ("frozen support", frozen_support_rows, accumulator.rows_passing_frozen),
    ):
        if per_accession != row_level:
            raise ScaleImpactInconsistentError(
                f"per-accession {label} ({per_accession}) does not sum to the "
                f"row-level {label} ({row_level})"
            )

    rows_newly_rejected = accumulator.rows_newly_rejected
    rows_newly_rejected_consequential = sum(
        record["newly_rejected_rows"]
        for record in outcomes.values()
        if record["rejection_consequential"]
    )
    rows_newly_rejected_redundant = sum(
        record["newly_rejected_rows"]
        for record in outcomes.values()
        if record["total_hit_rows"] > 0 and not record["rejection_consequential"]
    )
    if rows_newly_rejected != (
        rows_newly_rejected_consequential + rows_newly_rejected_redundant
    ):
        raise ScaleImpactInconsistentError(
            "a newly rejected row is neither consequential nor redundant: "
            f"{rows_newly_rejected} != {rows_newly_rejected_consequential} + "
            f"{rows_newly_rejected_redundant}"
        )

    accessions_lost = sorted(
        accession
        for accession, record in outcomes.items()
        if record["lost_all_support"]
    )
    no_hit_accessions = sorted(
        accession
        for accession, record in outcomes.items()
        if record["total_hit_rows"] == 0
    )
    lost_consequential = sum(
        1 for record in outcomes.values() if record["rejection_consequential"]
    )
    retained = len(outcomes) - len(accessions_lost) - len(no_hit_accessions)

    overall = {
        "universe_accessions": len(outcomes),
        "accessions_with_hits": len(outcomes) - len(no_hit_accessions),
        "accessions_in_no_hit_row": len(no_hit_accessions),
        "accessions_retaining_rescaled_support": retained,
        "accessions_lost_all_support": len(accessions_lost),
        "accessions_lost_all_support_consequential": lost_consequential,
        "rows_joining_universe": accumulator.rows_joining_universe,
        "rows_passing_frozen": accumulator.rows_passing_frozen,
        "rows_passing_rescaled": accumulator.rows_passing_rescaled,
        "rows_newly_rejected": rows_newly_rejected,
        "rows_newly_rejected_consequential": rows_newly_rejected_consequential,
        "rows_newly_rejected_redundant": rows_newly_rejected_redundant,
    }

    # 5. cross-tabs: every accession lands in exactly one bucket per dimension
    for accession, record in outcomes.items():
        family = record["family"]
        for tab, value in (
            (accumulator.tabs["family"], family),
            (accumulator.tabs["source_layer"], record["source_layer"]),
            (accumulator.tabs["layer"], record["layer"]),
            (
                accumulator.tabs["layer_x_family"],
                f"{record['layer']}|{family}",
            ),
        ):
            entry = tab.bucket(value)
            entry["universe_accessions"] += 1
            entry["rows_joining_universe"] += record["total_hit_rows"]
            entry["rows_passing_frozen"] += record["frozen_support_rows"]
            entry["rows_passing_rescaled"] += record["rescaled_support_rows"]
            entry["rows_newly_rejected"] += record["newly_rejected_rows"]
            if record["total_hit_rows"] > 0:
                entry["accessions_with_hits"] += 1
            else:
                entry["accessions_in_no_hit_row"] += 1
            if record["rescaled_support_rows"] > 0:
                entry["accessions_retaining_rescaled_support"] += 1
            if record["lost_all_support"]:
                entry["accessions_lost_all_support"] += 1
            if record["rejection_consequential"]:
                entry["accessions_lost_all_support_consequential"] += 1
                entry["rows_newly_rejected_consequential"] += record[
                    "newly_rejected_rows"
                ]
            else:
                entry["rows_newly_rejected_redundant"] += record[
                    "newly_rejected_rows"
                ]

    for tab in accumulator.tabs.values():
        for field, expected in overall.items():
            total = sum(entry[field] for entry in tab.buckets.values())
            if total != expected:
                raise ScaleImpactInconsistentError(
                    f"cross-tab {tab.dimension}.{field} sums to {total}, expected "
                    f"{expected}"
                )

    # 6. the join key: measured per column, ambiguity fatal, nothing assumed
    per_key: Dict[str, Dict[str, Any]] = {}
    for key in KEY_COLUMNS:
        stats = accumulator.key_stats[key]
        others = set().union(
            *(accumulator.key_stats[other].users for other in KEY_COLUMNS if other != key)
        )
        per_key[key] = {
            "found": stats.found,
            "rows_matched": stats.rows_matched,
            "found_not_in_other_keys": sorted(stats.users - others)[:20],
        }
    found_in_any = set().union(
        *(accumulator.key_stats[key].users for key in KEY_COLUMNS)
    )
    if accumulator.conflicts:
        raise JoinKeyAmbiguityError(accumulator.conflicts)
    if join_key == AUTO_JOIN_KEY:
        resolved = next(
            (key for key in KEY_COLUMNS if accumulator.key_stats[key].users == found_in_any),
            None,
        )
        if resolved is None or not found_in_any:
            raise JoinKeyAmbiguityError({})
        join_basis = "auto"
    else:
        resolved = join_key
        join_basis = "pinned"
    join_report = {
        "join_key": resolved,
        "join_key_basis": join_basis,
        "key_columns": list(KEY_COLUMNS),
        "universe_accessions": len(universe),
        "found_in_protein": per_key["protein"]["found"],
        "found_in_tacc": per_key["tacc"]["found"],
        "found_in_neither": len(universe) - len(found_in_any),
        "columns_agree": accumulator.columns_agree(),
        "conflicting_rows": len(accumulator.conflicts),
        "per_key": per_key,
        "proven": (
            "accession == protein"
            if per_key["protein"]["found"] == len(universe) and resolved == "protein"
            else "not_proven"
        ),
        "join_key_ambiguity": bool(accumulator.conflicts),
    }

    # 7. the join-completeness scope of this run
    if limit is None:
        no_hit_in_scope = True
        out_of_limit_scope: List[str] = []
    else:
        # Only rows up to the limit were seen, so an accession with no row may
        # simply not have been reached yet: report it, but do not count it as a
        # missing join.
        seen_accessions = set(accumulator.tallies)
        out_of_limit_scope = [
            accession
            for accession in no_hit_accessions
            if accession not in seen_accessions
        ]
        no_hit_in_scope = False
    no_hit_measured = [
        accession
        for accession in no_hit_accessions
        if accession not in set(out_of_limit_scope)
    ]

    summary: Dict[str, Any] = {
        "schema_version": "1.0",
        "task": "F7 follow-up: corrected E-value scale vs the frozen v1 universe",
        "inputs": {
            "universe": str(universe_path),
            "hits": str(hits_path),
            "delta_report": str(delta_path) if delta_path else None,
            "universe_sha256": sha256_file(universe_path) if args.sha256 else None,
            "hits_sha256": sha256_file(hits_path) if args.sha256 else None,
        },
        "join": join_report,
        "scale_column": scale_meta,
        "threshold": format_threshold(threshold),
        "threshold_rule": THRESHOLD_RULE,
        "rows_passing_original_basis": (
            "recomputed from the preserved per-shard E-value cell with "
            f"{THRESHOLD_RULE} at {format_threshold(threshold)}"
        ),
        "universe_accessions": len(universe),
        "accessions_processed": len(outcomes),
        "expected_universe_size": expected_size,
        "accessions_with_hits": overall["accessions_with_hits"],
        "accessions_retaining_rescaled_support": retained,
        "accessions_with_zero_rescaled_support": len(accessions_lost),
        "accessions_lost_all_support": len(accessions_lost),
        "accessions_lost_all_support_consequential": lost_consequential,
        "accessions_in_no_hit_row": len(no_hit_accessions),
        "no_hit_accessions": no_hit_accessions,
        "no_hit_accessions_in_scope": len(no_hit_measured) if no_hit_in_scope else None,
        "no_hit_check_applicable": no_hit_in_scope,
        "expected_no_hit": expected_no_hit,
        "no_hit_accessions_within_expected": bool(
            expected_no_hit is not None
            and no_hit_in_scope
            and len(no_hit_accessions) == expected_no_hit
        ),
        "no_hit_check_downgraded": bool(args.allow_no_hit_accessions),
        "accessions_out_of_limit_scope": len(out_of_limit_scope),
        "join_completeness_provable": no_hit_in_scope,
        "rows_total": accumulator.rows_total,
        "rows_joining_universe": accumulator.rows_joining_universe,
        "rows_not_in_universe": accumulator.rows_not_in_universe,
        "rows_limit_scope": accumulator.rows_total,
        "rows_passing_frozen": accumulator.rows_passing_frozen,
        "rows_failing_frozen": accumulator.rows_failing_frozen,
        "rows_passing_rescaled": accumulator.rows_passing_rescaled,
        "rows_failing_rescaled": accumulator.rows_failing_rescaled,
        "rows_passing_original_threshold": accumulator.rows_passing_frozen,
        "rows_passing_rescaled_threshold": accumulator.rows_passing_rescaled,
        "delta_rows": (
            accumulator.rows_passing_rescaled - accumulator.rows_passing_frozen
        ),
        "rows_newly_rejected": rows_newly_rejected,
        "rows_newly_rejected_consequential": rows_newly_rejected_consequential,
        "rows_newly_rejected_redundant": rows_newly_rejected_redundant,
        "affected_accessions": len(accessions_lost),
        "affected_accessions_consequential": lost_consequential,
        "question_answers": {
            "1_keeps_at_least_one_passing_row": retained,
            "1_universe_accessions": len(outcomes),
            "2_loses_all_support": len(accessions_lost),
            "2_loses_all_support_where_the_loss_is_consequential": lost_consequential,
            "2_kept_by_at_least_one_other_row": len(accessions_lost)
            - lost_consequential,
            "3_decomposition": (
                "cross_tabs.{family,source_layer,layer,layer_x_family} carry "
                "universe/passing/lost/no-hit/rejected counters per value"
            ),
            "4_newly_rejected_rows_consequential": rows_newly_rejected_consequential,
            "4_newly_rejected_rows_redundant": rows_newly_rejected_redundant,
            "4_newly_rejected_rows_total": rows_newly_rejected,
            "5_accessions_in_no_hit_row": len(no_hit_accessions),
            "5_explained": bool(
                expected_no_hit is not None
                and no_hit_in_scope
                and len(no_hit_accessions) == expected_no_hit
            ),
        },
        "interpretation_boundaries": [
            "A lost accession is one whose every frozen hit row fails E < "
            f"{format_threshold(threshold)} on the corrected full-library scale; "
            "it is not a claim that the protein is absent from GTDB.",
            "The loss means the sequence search evidence for that candidate does "
            "not survive the longer-library scale; any other evidence layer is "
            "untouched by this measurement.",
            "Newly rejected rows can only appear where the corrected scale is "
            "stricter, which holds whenever every per-shard Z is below Z_total.",
        ],
        "cross_tabs": {name: tab.as_dict() for name, tab in accumulator.tabs.items()},
        "bounded_state": {
            "universe_size": len(universe),
            "estimated_bytes": memory_bound_bytes(len(universe)),
            "model": (
                "O(universe) state: one slot record per universe accession "
                "(counters plus a has-support flag) plus one aggregate counter "
                "per cross-tab bucket; rows claiming no universe accession are "
                "counted and dropped, so the state does not grow with the "
                "6,743,197-row / 985 MB hit table"
            ),
            "bytes_per_universe_accession": BYTES_PER_UNIVERSE_ACCESSION,
            "fixed_state_bytes": FIXED_STATE_BYTES,
            "max_join_keys": MAX_JOIN_KEYS,
            "max_line_length": MAX_LINE_LENGTH,
        },
        "limited": limit is not None,
        "limit": limit,
        "limited_note": (
            f"limited=true: only the first {limit} hit row(s) were read, so no "
            "whole-table count and no join-completeness claim is made"
            if limit is not None
            else None
        ),
        "outputs": {
            "impact_tsv": str(output_dir / IMPACT_FILE),
            "summary_json": str(output_dir / SUMMARY_FILE),
            "affected_tsv": str(affected_target),
        },
        "frozen_evidence_modified": False,
    }

    if delta_path is not None:
        cross_check = _cross_check_delta_report(
            delta_path, summary, limited=limit is not None
        )
        summary["frozen_delta_cross_check"] = cross_check
        if cross_check["status"] == "mismatch":
            raise ScaleImpactInconsistentError(
                "the recomputed row counts disagree with the supplied F7 delta "
                f"report {delta_path}: "
                f"{json.dumps(cross_check['mismatches'], sort_keys=True)}"
            )
    else:
        summary["frozen_delta_cross_check"] = {
            "status": "not_supplied",
            "path": None,
            "mismatches": None,
        }

    # 8. write.  Only after every check that does not need the tables on disk.
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_text_atomically(
        output_dir / IMPACT_FILE,
        _impact_tsv(outcomes, sorted(outcomes)),
        bool(args.force),
    )
    _write_text_atomically(
        affected_target,
        _impact_tsv(outcomes, accessions_lost),
        bool(args.force),
    )

    # 9. join completeness.  The evidence tables are already published because
    #    they *are* the evidence for the claim; the exit code is still non-zero.
    unexplained = no_hit_in_scope and len(no_hit_accessions) > 0 and not (
        expected_no_hit is not None and len(no_hit_accessions) == expected_no_hit
    )
    _write_text_atomically(
        output_dir / SUMMARY_FILE,
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        bool(args.force),
    )
    if unexplained and not args.allow_no_hit_accessions:
        raise ScaleImpactInconsistentError(
            f"{len(no_hit_accessions)} universe accession(s) appear in no hit row "
            "at all (first: "
            + ", ".join(no_hit_accessions[:5])
            + "); the universe was derived from these hits, so this join is "
            "incomplete. Explain it with --expected-no-hit N (the measured count) "
            "or downgrade it with --allow-no-hit-accessions."
        )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    status, args = parse_args_never_exits(argv)
    if args is None:
        return status
    try:
        payload = run(args)
    except (
        JoinKeyAmbiguityError,
        ScaleColumnError,
        ScaleImpactInconsistentError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_INCONSISTENT
    except FileExistsError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_REFUSED
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_USAGE
    if not args.quiet:
        print(
            "universe %d: %d keep rescaled support, %d lose all support "
            "(%d consequential), %d in no hit row; newly rejected rows %d "
            "(%d consequential / %d redundant)"
            % (
                payload["universe_accessions"],
                payload["accessions_retaining_rescaled_support"],
                payload["accessions_lost_all_support"],
                payload["accessions_lost_all_support_consequential"],
                payload["accessions_in_no_hit_row"],
                payload["rows_newly_rejected"],
                payload["rows_newly_rejected_consequential"],
                payload["rows_newly_rejected_redundant"],
            ),
            file=sys.stderr,
        )
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
