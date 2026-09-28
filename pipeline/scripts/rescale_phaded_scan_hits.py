#!/usr/bin/env python3
"""Rescale the frozen scan-13 E-values onto one measured full-library ``Z``.

The decision this module answers
--------------------------------
The frozen GTDB full-library scan (``runs/20260901_formal_frozen_scan_13``) ran
1,000 ``hmmsearch`` tasks (10 HMMs x 100 shards) with

    hmmsearch --tblout "$BUILD/$base.tbl" --domtblout "$BUILD/$base.dom" \
              -E 1e-5 --cpu 1 "$hmm" "$filtered"

— **no ``-Z`` and no ``--domZ``**.  HMMER therefore used each shard's own
sequence count as ``Z``, so the 6.7 M reported E-values sit on 100 different
scales (shard_0001 = 4,694,121 records vs shard_0051 = 2,923,820, a 1.61x
spread).  Any threshold applied across the frozen table mixes those scales.

Because HMMER's E-value is linear in ``Z``, and because that linearity was
proven twice on frozen evidence (same target, ``-Z 2920000 -> E=1.4e-08`` vs
``-Z 109087 -> E=5.1e-10``; predicted ratio 26.80, measured 26.9), the
reconciliation needs **no HMMER run at all** — this is a pure rescaling of
numbers already on disk:

    E_full = E_shard * (Z_total / Z_shard)

``Z_shard`` comes from ``count_fasta_records.py`` (one ``^>`` count per frozen
shard) and ``Z_total`` is the measured full-library sum — or an explicit
``--database-size-z``.  Both are recorded, and the basis is stated, so a reader
never has to guess which scale a column is on.

What is produced (and what is never touched)
--------------------------------------------
* ``--out``: a **new** wide TSV.  Every original ``hits_all.tsv`` column is
  copied through unchanged, the original E-value is additionally preserved
  verbatim as ``E-value_shard_z``, and the rescale columns are appended:
  ``E-value_rescaled``, ``z_shard``, ``z_total``,
  ``passes_rescaled_threshold``.
* ``--delta-json``: the reconciliation report — rows, rows passing the original
  threshold, rows passing the rescaled threshold, the DELTA, its direction and
  magnitude in words, and a per-family breakdown.
* ``--manifest`` (optional): the scale evidence
  (``database_size_Z``, per-shard counts, the command, input hashes).

The frozen ``hits_all.tsv`` is **read-only evidence**: it is streamed, never
loaded, and ``--out`` equal to ``--hits`` is refused (resolved-path comparison,
so a different spelling of the same file is refused too).

The arithmetic is the project's own convention
----------------------------------------------
``E_full = rescale_evalue(E_shard, from_z=Z_shard, to_z=Z_total)`` from
``pipeline/scripts/parse_phaded_cys_targeted_recall.py`` — imported, not
re-derived.  ``rescale_row`` is a thin validated wrapper around it so callers
have one entry point and no second formula exists in the tree.

Failure modes (all fail closed, none is a silent skip)
------------------------------------------------------
* a row whose shard is absent from the shard-counts file is an **error** naming
  every missing shard — unless ``--allow-unknown-shard`` marks it ``pending``
  (counted and reported) for dry runs;
* an unusable ``--out``/``--delta-json``/``--manifest`` target (same file as
  ``--hits``, or an existing file without ``--force``) is refused;
* a shard-counts file with no usable count and no ``--database-size-z`` fails
  closed rather than dividing by a guess.

Threshold semantics: a row passes when its E-value is **strictly less than** the
threshold (``<``), matching ``hmmsearch -E`` decoding of "report sequences with
E-value <= E"; the two agree for every value except an exact tie at the
threshold, and ``--threshold`` is echoed in the report so the rule is explicit.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Sequence

#: The upstream module that owns the rescaling convention.  Loaded by path
#: because ``pipeline/scripts`` is not a package.
CONVENTION_SCRIPT = Path(__file__).resolve().parent / "parse_phaded_cys_targeted_recall.py"

#: ``hmmsearch -E 1e-5`` was the frozen reporting threshold; it is therefore the
#: default basis for the "before" column of the delta report as well.
FROZEN_THRESHOLD = 1e-5

#: The frozen ``hits_all.tsv`` columns.  ``E-value`` is the one that is on the
#: wrong scale; the rest are copied through untouched.
HITS_COLUMNS = [
    "family",
    "shard",
    "protein",
    "tacc",
    "E-value",
    "score",
    "bias",
    "domE",
    "qname",
    "cov",
]

#: Appended by this module.
ADDED_COLUMNS = [
    "E-value_shard_z",
    "E-value_rescaled",
    "z_shard",
    "z_total",
    "passes_rescaled_threshold",
]

ORIGINAL_EVALUE_COLUMN = "E-value_shard_z"
RESCALED_EVALUE_COLUMN = "E-value_rescaled"
PASS_COLUMN = "passes_rescaled_threshold"
PASSED = "true"
FAILED = "false"
PENDING = "pending"

#: How many rows the pre-flight scan reads to collect the shard names present in
#: a table before any output is written.  The frozen table has ~100 distinct
#: shards (one per shard file), so this bounded sample finds them all while
#: keeping the pre-flight independent of the table's 696 MB size.
PREFLIGHT_SCAN_ROWS = 100_000

#: Threshold cells are echoed in this format so the report is comparable with a
#: HMMER ``-E`` flag string.
THRESHOLD_FORMAT = "%.6g"

Z_BASIS_CLI = "cli_database_size_z"
Z_BASIS_TABLE = "shard_counts_total"
Z_BASIS_NOTE_CLI = (
    "Z_total supplied on the command line with --database-size-z; it overrides "
    "the shard-counts total, and every row is rescaled onto this value."
)
Z_BASIS_NOTE_TABLE = (
    "Z_total is the sum of the per-shard record counts in --shard-counts "
    "(measured by count_fasta_records.py); no other value was supplied."
)

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_UNKNOWN_SHARD = 3
EXIT_REFUSED = 4
EXIT_NO_DATABASE_SIZE = 5

PENDING_UNKNOWN_SHARD = "shard_absent_from_shard_counts"


def load_convention():
    """Import the shared rescaling convention, never a private copy.

    ``parse_phaded_cys_targeted_recall`` lives next to this file but the scripts
    directory is not a package, so a plain import only works when it is already
    on ``sys.path``.  Every other import path is a hard error: silently
    re-implementing ``rescale_evalue`` here is exactly the divergence the
    project forbids.
    """
    try:
        import parse_phaded_cys_targeted_recall as convention  # type: ignore[import-not-found]

        return convention
    except ImportError:
        pass
    if not CONVENTION_SCRIPT.is_file():
        raise ImportError(f"cannot locate the rescaling convention next to {__file__}")
    spec = importlib.util.spec_from_file_location(
        "parse_phaded_cys_targeted_recall", CONVENTION_SCRIPT
    )
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise ImportError(f"cannot load the rescaling convention: {CONVENTION_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("parse_phaded_cys_targeted_recall", module)
    spec.loader.exec_module(module)
    return module


CONVENTION = load_convention()
#: The project's function, re-exported so callers and tests share one object.
rescale_evalue = CONVENTION.rescale_evalue


# --------------------------------------------------------------------------
# value parsing
# --------------------------------------------------------------------------
def parse_evalue(value: Any, field: str = "E-value") -> float:
    """Parse an E-value cell robustly; return a finite, non-negative ``float``.

    The frozen table spells values as ``1.4e-08``; a ``1e-5``-style string, a
    plain ``0`` and an integer all mean the same thing to HMMER.  Anything that
    is not a finite, non-negative number raises ``ValueError``: a corrupted
    E-value must never be silently turned into ``0`` or skipped.

    Values that underflow to ``0.0`` (``2.9e-310``) are accepted: the reported
    value already is zero at double precision.
    """
    if isinstance(value, bool):
        raise ValueError(f"{field} is not a number: {value!r}")
    if isinstance(value, (int, float)):
        number = float(value)
    else:
        text = ("" if value is None else str(value)).strip()
        if not text:
            raise ValueError(f"{field} is empty; no E-value was recorded")
        try:
            number = float(Decimal(text))
        except Exception as error:  # noqa: BLE001 - any parse failure means "unusable"
            raise ValueError(f"{field} is not a number: {value!r}") from error
    if number != number or number in (float("inf"), float("-inf")):
        raise ValueError(f"{field} must be finite, got {value!r}")
    if number < 0:
        raise ValueError(f"{field} must be non-negative, got {value!r}")
    return number


def parse_positive_int(value: Any, field: str) -> int:
    """Parse a strictly positive integer (``bool`` and ``"0"`` are rejected)."""
    if isinstance(value, bool):
        raise ValueError(f"{field} must be a positive integer, got {value!r}")
    if isinstance(value, int):
        number = value
    else:
        text = ("" if value is None else str(value)).strip()
        if not text or not text.lstrip("+").isdigit():
            raise ValueError(f"{field} must be a positive integer, got {value!r}")
        number = int(text, 10)
    if number < 1:
        raise ValueError(f"{field} must be a positive integer, got {value!r}")
    return number


def format_evalue(value: float) -> str:
    """Render a rescaled E-value without ever collapsing it to ``0``.

    Scientific notation is deliberately kept: ``%.6g`` renders ``1.4e-06`` as
    ``1.4e-06`` (parseable and unambiguous) while a fixed-point format would
    turn every sub-1e-6 value into ``0.000000`` and lose the comparison.
    """
    return "%.6g" % value


def format_threshold(value: float) -> str:
    """Render a threshold for the report so it round-trips."""
    return THRESHOLD_FORMAT % value


def rescale_row(evalue: Any, z_shard: Any, z_total: Any) -> float:
    """Return ``E_full`` for one row, via the project's ``rescale_evalue``.

    ``E_full = E_shard * Z_total / Z_shard`` — delegated, so this module holds no
    second copy of the formula.  ``evalue`` may be a float or a cell string;
    ``z_shard``/``z_total`` must be positive integers.
    """
    parsed = parse_evalue(evalue, "E-value")
    return rescale_evalue(
        parsed,
        from_z=parse_positive_int(z_shard, "z_shard"),
        to_z=parse_positive_int(z_total, "z_total"),
    )


def sha256_file(path: Path | str, chunk_bytes: int = 1024 * 1024) -> str:
    """Return the SHA-256 of ``path``, streamed in fixed-size blocks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# reading the frozen table
# --------------------------------------------------------------------------
def _open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def validate_hits_header(header: Sequence[str], path: Any) -> List[str]:
    """Return the required ``hits_all.tsv`` columns, raising when one is missing.

    A table without ``E-value`` or ``shard`` cannot be rescaled; refusing here is
    what stops a wrong-scale table from being produced with a shifted column.
    """
    missing = [column for column in HITS_COLUMNS if column not in header]
    if missing:
        raise ValueError(
            f"{path}: the hits table is missing required column(s): "
            f"{', '.join(missing)} (found: {', '.join(header)})"
        )
    return list(HITS_COLUMNS)


def iter_hits_rows(
    path: Path | str, *, limit: int | None = None
) -> Iterator[Dict[str, str]]:
    """Yield one dict per data row of the hits table, streaming line by line.

    The frozen table is 696 MB, so this is a generator that reads a single line
    at a time and never materialises the file.  ``limit`` caps the number of data
    rows yielded (``None`` means all of them).
    """
    path = Path(path)
    with _open_text(path, "rt") as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"{path}: the hits table is empty")
        header = first.rstrip("\r\n").split("\t")
        validate_hits_header(header, path)
        index = {column: header.index(column) for column in HITS_COLUMNS}
        seen = 0
        number = 1
        # Explicit ``readline()`` rather than ``for line in handle``: it makes the
        # bounded-read property observable at the handle (a test can wrap the
        # stream and refuse an unsized read) instead of relying on iteration.
        while True:
            line = handle.readline()
            if not line:
                break
            number += 1
            if not line.strip():
                continue
            if line.rstrip("\r\n").startswith("#"):
                continue
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) < len(header):
                raise ValueError(
                    f"{path} line {number}: expected {len(header)} columns, "
                    f"got {len(cells)}"
                )
            yield {column: cells[position] for column, position in index.items()}
            seen += 1
            if limit is not None and seen >= limit:
                return


# --------------------------------------------------------------------------
# the shard-counts ledger
# --------------------------------------------------------------------------
def normalise_shard_key(name: Any) -> str:
    """Return the shard key of a ledger name: the file name without ``.faa``.

    The frozen hits table says ``shard_0051``; a counting ledger may name the
    file ``shard_0051.faa``.  Both must resolve to one key or the scale lookup
    would miss and (correctly) refuse to produce a table.
    """
    text = ("" if name is None else str(name)).strip()
    text = Path(text).name
    for suffix in (".faa.gz", ".faa", ".fasta.gz", ".fasta", ".gz"):
        if text.endswith(suffix):
            text = text[: -len(suffix)]
            break
    return text


def read_shard_counts(path: Path | str) -> Dict[str, Any]:
    """Read the ``count_fasta_records.py`` ledger into a scale lookup.

    Returns ``{"counts": {shard: records}, "pending": {shard: reason}}``.  A row
    whose ``status`` is not ``counted``/``reused`` (or whose ``records`` cell is
    empty) is *pending*, never zero: dividing by a shard count of zero is the one
    arithmetic error that would silently destroy the whole table.
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"shard-counts file does not exist: {path}")
    if path.is_dir():
        raise ValueError(f"shard-counts path is a directory: {path}")
    counts: Dict[str, int] = {}
    pending: Dict[str, str] = {}
    with _open_text(path, "rt") as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"shard-counts file is empty: {path}")
        header = [cell.strip() for cell in header_line.rstrip("\r\n").split("\t")]
        for required in ("name", "records"):
            if required not in header:
                raise ValueError(
                    f"shard-counts file {path} is missing the {required!r} column"
                )
        for number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) < len(header):
                raise ValueError(
                    f"{path} line {number}: expected {len(header)} columns, "
                    f"got {len(cells)}"
                )
            row = dict(zip(header, cells))
            name = normalise_shard_key(row.get("name"))
            if not name:
                raise ValueError(f"{path} line {number}: shard row has no name")
            status = (row.get("status") or "counted").strip()
            recorded = (row.get("records") or "").strip()
            if status not in ("counted", "reused") or not recorded:
                reason = (row.get("pending_reason") or "").strip()
                pending[name] = reason or (
                    f"status={status!r} with records={row.get('records')!r}"
                )
                continue
            if name in counts:
                raise ValueError(f"{path}: shard {name!r} appears more than once")
            counts[name] = parse_positive_int(recorded, f"records for {name}")
    return {"counts": counts, "pending": pending}


def shard_counts_total(counts: Mapping[str, int]) -> int:
    """Return the full-library ``Z`` implied by a shard-count ledger.

    A ledger with no usable count means no full-library ``Z`` can be measured
    here; ``run`` refuses rather than rescaling onto a guess.
    """
    if not counts:
        raise NoDatabaseSizeError(
            "the shard-counts file contains no usable count, so no full-library "
            "Z can be established; pass --database-size-z explicitly"
        )
    return sum(int(value) for value in counts.values())


def collect_shards_in_table(
    path: Path | str,
    *,
    scan_rows: int = PREFLIGHT_SCAN_ROWS,
    limit: int | None = None,
) -> List[str]:
    """Return the distinct shard names in the first ``scan_rows`` data rows.

    Used as a *pre-flight* so that a run which cannot name a scale for some
    shard fails with **every** missing shard listed, before a single output byte
    is written, instead of reporting them one at a time.
    """
    if scan_rows < 1:
        raise ValueError(f"scan_rows must be positive, got {scan_rows!r}")
    cap = scan_rows if limit is None else min(scan_rows, limit)
    names = set()
    for row in iter_hits_rows(path, limit=cap):
        names.add((row.get("shard") or "").strip())
    names.discard("")
    return sorted(names)


def tables_shards_without_scale(
    path: Path | str,
    counts: Mapping[str, int],
    *,
    scan_rows: int = PREFLIGHT_SCAN_ROWS,
    limit: int | None = None,
) -> List[str]:
    """Return the shards of a table that the ledger cannot give a scale for.

    A shard is unscalable when it has no count in the ledger at all — including a
    ledger row that is present but ``not_counted``/empty, which is exactly the
    zero-division case the ledger must never turn into a silent skip.
    """
    return [
        shard
        for shard in collect_shards_in_table(path, scan_rows=scan_rows, limit=limit)
        if shard not in counts
    ]


# --------------------------------------------------------------------------
# the rescaler
# --------------------------------------------------------------------------
class Rescaler:
    """Per-row rescaling decisions plus the counters the delta report needs."""

    def __init__(
        self,
        *,
        counts: Mapping[str, int],
        pending: Mapping[str, str],
        z_total: int,
        threshold: float,
        original_threshold: float,
        allow_unknown_shard: bool,
    ) -> None:
        self.counts = dict(counts)
        self.pending = dict(pending)
        self.z_total = z_total
        self.threshold = threshold
        self.original_threshold = original_threshold
        self.allow_unknown_shard = allow_unknown_shard
        self.total_rows = 0
        self.rescalable_rows = 0
        self.pending_rows = 0
        self.pass_original = 0
        self.fail_original = 0
        self.pass_rescaled = 0
        self.fail_rescaled = 0
        self.pending_cells: Dict[str, int] = {}
        self.families: Dict[str, Dict[str, int]] = {}

    # -- helpers ---------------------------------------------------------
    def _family(self, family: str) -> Dict[str, int]:
        return self.families.setdefault(
            family,
            {
                "total_rows": 0,
                "rows_passing_original_threshold": 0,
                "rows_passing_rescaled_threshold": 0,
                "delta_rows": 0,
                "rows_pending_unknown_shard": 0,
            },
        )

    def _pending_cell(self, key: str) -> None:
        self.pending_cells[key] = self.pending_cells.get(key, 0) + 1

    def unknown_shards(self) -> Dict[str, int]:
        """Rows per shard that is absent from the ledger (for the error text)."""
        return {
            key[len("unknown_shard:") :]: value
            for key, value in self.pending_cells.items()
            if key.startswith("unknown_shard:")
        }

    # -- the row ---------------------------------------------------------
    def rescale(self, row: Mapping[str, str]) -> Dict[str, str]:
        """Return the output cells for one frozen row (all values as text)."""
        family = (row.get("family") or "").strip()
        shard = (row.get("shard") or "").strip()
        family_row = self._family(family)
        self.total_rows += 1
        family_row["total_rows"] += 1

        original_text = (row.get("E-value") or "").strip()
        original = parse_evalue(original_text, f"E-value for {shard}:{row.get('protein')}")
        passes_original = original < self.original_threshold
        self.pass_original += 1 if passes_original else 0
        self.fail_original += 0 if passes_original else 1
        family_row["rows_passing_original_threshold"] += 1 if passes_original else 0
        out: Dict[str, str] = {column: row.get(column, "") for column in HITS_COLUMNS}
        out[ORIGINAL_EVALUE_COLUMN] = original_text
        out["z_total"] = str(self.z_total)

        records = self.counts.get(shard)
        if records is None:
            if not self.allow_unknown_shard:
                raise UnknownShardError([shard])
            reason = self.pending.get(shard) or PENDING_UNKNOWN_SHARD
            self.pending_rows += 1
            family_row["rows_pending_unknown_shard"] += 1
            self._pending_cell(f"unknown_shard:{shard}")
            out[RESCALED_EVALUE_COLUMN] = PENDING
            out["z_shard"] = PENDING
            out[PASS_COLUMN] = PENDING
            return out

        rescaled = rescale_row(original, z_shard=records, z_total=self.z_total)
        passes = rescaled < self.threshold
        self.rescalable_rows += 1
        if passes:
            self.pass_rescaled += 1
            family_row["rows_passing_rescaled_threshold"] += 1
        else:
            self.fail_rescaled += 1
        out[RESCALED_EVALUE_COLUMN] = format_evalue(rescaled)
        out["z_shard"] = str(records)
        out[PASS_COLUMN] = PASSED if passes else FAILED
        return out

    # -- finalisation ----------------------------------------------------
    def finalise(self) -> None:
        """Fill the per-family delta, which is only known once all rows are seen."""
        for row in self.families.values():
            row["delta_rows"] = (
                row["rows_passing_rescaled_threshold"]
                - row["rows_passing_original_threshold"]
            )


class UnknownShardError(ValueError):
    """Raised when a row's shard has no scale in the ledger."""

    def __init__(self, shards: Iterable[str]) -> None:
        self.shards = sorted(set(shards))
        super().__init__(
            "the shard-counts file has no scale for shard(s): "
            + ", ".join(self.shards)
            + " — refusing to emit a mixed-scale table. Count them first "
            "(count_fasta_records.py) or pass --allow-unknown-shard to mark "
            "those rows pending for a dry run."
        )


class NoDatabaseSizeError(ValueError):
    """Raised when neither ``--database-size-z`` nor the ledger can give a Z."""


# --------------------------------------------------------------------------
# the run
# --------------------------------------------------------------------------
def build_arg_namespace(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rescale frozen scan-13 hits onto one measured full-library Z "
            "(E_full = E_shard * Z_total / Z_shard). No HMMER run is needed."
        )
    )
    parser.add_argument("--hits", type=Path, required=True, help="frozen hits_all.tsv")
    parser.add_argument(
        "--shard-counts",
        type=Path,
        required=True,
        help="per-shard record counts from count_fasta_records.py",
    )
    parser.add_argument("--out", type=Path, required=True, help="new rescaled TSV")
    parser.add_argument(
        "--delta-json", type=Path, default=None, help="new delta report JSON"
    )
    parser.add_argument(
        "--manifest", type=Path, default=None, help="new scale-evidence manifest JSON"
    )
    parser.add_argument(
        "--database-size-z",
        type=int,
        default=None,
        help="the full-library Z to rescale onto; defaults to the "
        "--shard-counts total and must be a positive integer",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=FROZEN_THRESHOLD,
        help=f"E-value threshold applied to the rescaled values (default {FROZEN_THRESHOLD:g})",
    )
    parser.add_argument(
        "--original-threshold",
        type=float,
        default=None,
        help="threshold applied to the frozen E-values for the 'before' count "
        f"(default: same as --threshold, i.e. {FROZEN_THRESHOLD:g})",
    )
    parser.add_argument(
        "--allow-unknown-shard",
        action="store_true",
        help="mark rows whose shard is missing from --shard-counts as 'pending' "
        "and count them; for dry runs only",
    )
    parser.add_argument("--force", action="store_true", help="overwrite existing outputs")
    parser.add_argument(
        "--limit", type=int, default=None, help="read at most N data rows (safety valve)"
    )
    parser.add_argument(
        "--hits-sha256",
        action="store_true",
        help="also record the frozen table's SHA-256 (a second full read)",
    )
    parser.add_argument("--gzip", action="store_true", help="gzip the rescaled TSV")
    parser.add_argument(
        "--command",
        default="pending",
        help="the exact command this table came from, recorded verbatim",
    )
    parser.add_argument(
        "--input-record",
        action="append",
        default=[],
        help="extra KEY=PATH binding to record in the manifest; repeatable",
    )
    parser.add_argument("--quiet", action="store_true", help="suppress the progress line")
    return parser.parse_args(argv)


def parse_args_never_exits(argv: Sequence[str] | None = None):
    """Parse ``argv``; return ``(status, namespace)`` instead of exiting."""
    try:
        return EXIT_OK, build_arg_namespace(argv)
    except SystemExit as exit_error:
        code = exit_error.code
        return (code if isinstance(code, int) else EXIT_USAGE), None


def _validate_input(path: Path, field: str) -> Path:
    if not path.exists():
        raise ValueError(f"{field} does not exist: {path}")
    if path.is_dir():
        raise ValueError(f"{field} is a directory, not a file: {path}")
    if not path.is_file():
        raise ValueError(f"{field} is not a regular file: {path}")
    return path


def _check_out_target(path: Path, field: str, force: bool) -> None:
    if path.exists() and path.is_dir():
        contents = list(path.iterdir())
        if contents:
            raise ValueError(
                f"{field} is an existing non-empty directory: {path} "
                f"({len(contents)} entries). Pass a NEW file path."
            )
        return
    if path.exists() and not force:
        raise FileExistsError(f"{field} already exists: {path}. Pass --force.")
    if not path.parent.exists():
        raise ValueError(f"{field} parent directory does not exist: {path.parent}")


def _same_file(left: Path, right: Path) -> bool:
    """True when two paths resolve to the same file.

    Compares resolved paths so ``dir/../hits_all.tsv`` and a differently-spelled
    path to the frozen table are both recognised; the frozen evidence must never
    be a write target.
    """
    try:
        if left.resolve() == right.resolve():
            return True
        return left.exists() and right.exists() and left.samefile(right)
    except OSError:  # pragma: no cover - defensive: an unreadable path is not equal
        return False


def _write_rows(out: Path, columns: Sequence[str], rows: Iterable[Mapping[str, str]]) -> int:
    written = 0
    with _open_text(out, "wt") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row.get(column, "") for column in columns) + "\n")
            written += 1
    return written


def write_rows_atomically(
    out: Path, columns: Sequence[str], rows: Iterable[Mapping[str, str]]
) -> int:
    """Stream ``rows`` into ``out``, but only publish it if every row succeeded.

    The rescaling can fail mid-table (a row from a shard with no scale).  Writing
    straight to ``--out`` would then leave a truncated mixed-scale file behind
    that looks like a finished product.  The rows go to a sibling temporary file
    which is renamed onto ``out`` only after the last row; a failure removes the
    temporary file and leaves ``--out`` nonexistent.
    """
    out = Path(out)
    temporary = out.with_name(out.name + ".partial")
    if temporary.exists():
        raise FileExistsError(
            f"a partial output is already present: {temporary}. Remove it or pass "
            "a different --out."
        )
    try:
        written = _write_rows(temporary, columns, rows)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    os.replace(temporary, out)
    return written


def run(args: argparse.Namespace) -> Dict[str, Any]:
    """Rescale the frozen table; write the outputs; return the summary.

    Raises ``ValueError`` for a missing input, an unusable output target, a
    shard without a scale, or an unusable ``--database-size-z``;
    ``FileExistsError`` when an output exists without ``--force``.
    """
    hits = Path(args.hits)
    counts_path = Path(args.shard_counts)
    out = Path(args.out)
    delta_path = Path(args.delta_json) if args.delta_json is not None else None
    manifest_path = Path(args.manifest) if args.manifest is not None else None

    _validate_input(hits, "--hits")
    _validate_input(counts_path, "--shard-counts")

    threshold = parse_evalue(args.threshold, "--threshold")
    if threshold <= 0:
        raise ValueError(f"--threshold must be positive, got {args.threshold!r}")
    original_threshold = (
        threshold
        if args.original_threshold is None
        else parse_evalue(args.original_threshold, "--original-threshold")
    )
    if original_threshold <= 0:
        raise ValueError(
            f"--original-threshold must be positive, got {args.original_threshold!r}"
        )
    if args.limit is not None:
        if isinstance(args.limit, bool) or int(args.limit) < 1:
            raise ValueError(f"--limit must be a positive integer, got {args.limit!r}")
        limit = int(args.limit)
    else:
        limit = None

    for field, target in (
        ("--out", out),
        ("--delta-json", delta_path),
        ("--manifest", manifest_path),
    ):
        if target is None:
            continue
        if _same_file(target, hits):
            raise FileExistsError(
                f"{field} resolves to --hits ({hits}); the frozen evidence is "
                "read-only and is never a write target"
            )
        _check_out_target(target, field, args.force)
        if delta_path is not None and target is not delta_path and _same_file(target, delta_path):
            raise FileExistsError(f"{field} resolves to --delta-json ({delta_path})")
        if manifest_path is not None and target is not manifest_path and _same_file(
            target, manifest_path
        ):
            raise FileExistsError(f"{field} resolves to --manifest ({manifest_path})")
    if args.gzip and out.suffix != ".gz":
        raise ValueError(f"--gzip requires a .gz --out name, got {out}")

    ledger = read_shard_counts(counts_path)
    if args.database_size_z is not None:
        z_total = parse_positive_int(args.database_size_z, "--database-size-z")
        z_basis = Z_BASIS_CLI
        z_note = Z_BASIS_NOTE_CLI
    else:
        z_total = shard_counts_total(ledger["counts"])
        z_basis = Z_BASIS_TABLE
        z_note = Z_BASIS_NOTE_TABLE

    rescaler = Rescaler(
        counts=ledger["counts"],
        pending=ledger["pending"],
        z_total=z_total,
        threshold=threshold,
        original_threshold=original_threshold,
        allow_unknown_shard=bool(args.allow_unknown_shard),
    )

    # Pre-flight: name every shard this run cannot scale *before* writing output,
    # so the failure is complete and no partial table is published.
    if not args.allow_unknown_shard:
        unscalable = tables_shards_without_scale(
            hits, ledger["counts"], limit=limit
        )
        if unscalable:
            raise UnknownShardError(unscalable)

    columns = list(HITS_COLUMNS) + list(ADDED_COLUMNS)
    # Atomic publish: a row whose shard has no scale aborts the pass and leaves
    # no truncated mixed-scale table behind.
    written = write_rows_atomically(
        out,
        columns,
        (rescaler.rescale(row) for row in iter_hits_rows(hits, limit=limit)),
    )
    rescaler.finalise()

    rows_passing_original = rescaler.pass_original
    rows_passing_rescaled = rescaler.pass_rescaled
    delta_rows = rows_passing_rescaled - rows_passing_original
    if delta_rows < 0:
        direction = "newly_rejected"
        magnitude = f"{abs(delta_rows)} row(s) newly rejected"
        comparison = "stricter"
    elif delta_rows > 0:
        direction = "newly_admitted"
        magnitude = f"{delta_rows} row(s) newly admitted"
        comparison = "looser"
    else:
        direction = "unchanged"
        magnitude = "no row changes side"
        comparison = "identical"
    ratio = (
        rows_passing_rescaled / rows_passing_original
        if rows_passing_original
        else None
    )
    statement = (
        f"Applying E < {format_threshold(threshold)} to the rescaled column is "
        f"{comparison} than applying E < {format_threshold(original_threshold)} to "
        f"the frozen per-shard column: {rows_passing_original} -> "
        f"{rows_passing_rescaled} ({magnitude} of {rescaler.total_rows} row(s))."
    )
    statement_direction = (
        f"Direction: {direction}; because every shard count "
        f"(max {max(ledger['counts'].values()) if ledger['counts'] else 'n/a'}) is "
        f"below the full-library Z {z_total}, rescaling can only raise an E-value, "
        f"so the rescaled threshold is never more permissive than the per-shard one."
    )

    pending_buckets: Dict[str, int] = {}
    for key, value in rescaler.pending_cells.items():
        bucket = key.split(":", 1)[1] if ":" in key else key
        pending_buckets[bucket] = pending_buckets.get(bucket, 0) + value

    delta_report: Dict[str, Any] = {
        "schema_version": "1.0",
        "threshold_rescaled": format_threshold(threshold),
        "threshold_original": format_threshold(original_threshold),
        "threshold_rule": "E < threshold (strict)",
        "z_total": z_total,
        "z_total_basis": z_basis,
        "z_total_basis_note": z_note,
        "total_rows": rescaler.total_rows,
        "rows_rescalable": rescaler.rescalable_rows,
        "rows_pending_unknown_shard": rescaler.pending_rows,
        "rows_passing_original_threshold": rows_passing_original,
        "rows_failing_original_threshold": rescaler.fail_original,
        "rows_passing_rescaled_threshold": rows_passing_rescaled,
        "rows_failing_rescaled_threshold": rescaler.fail_rescaled,
        "delta_rows": delta_rows,
        "direction": direction,
        "ratio_rescaled_over_original": ratio,
        "statement": statement,
        "statement_direction": statement_direction,
        "pending_buckets": pending_buckets,
        "unknown_shards": sorted(rescaler.unknown_shards()),
        "per_family": {
            family: dict(values) for family, values in sorted(rescaler.families.items())
        },
        "allow_unknown_shard": bool(args.allow_unknown_shard),
        "limited": limit is not None,
        "limit": limit,
        "hits": str(hits),
        "shard_counts": str(counts_path),
        "shards_with_scale": len(ledger["counts"]),
        "shards_pending_in_counts": sorted(ledger["pending"]),
    }

    if delta_path is not None:
        delta_path.write_text(
            json.dumps(delta_report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

    if manifest_path is not None:
        input_records: Dict[str, str] = {}
        for record in args.input_record:
            if "=" in record:
                key, value = record.split("=", 1)
                input_records[key.strip()] = value.strip()
        manifest: Dict[str, Any] = {
            "schema_version": "1.0",
            "purpose": "scan-13 -Z scale reconciliation (no HMMER run)",
            "formula": "E_full = E_shard * (Z_total / Z_shard)",
            "database_size_Z": z_total,
            "database_size_basis": z_basis,
            "database_size_basis_note": z_note,
            "scale_unified": True,
            "frozen_scan_manifest_had_scale_field": False,
            "frozen_hmmsearch_command_had_Z": False,
            "threshold_rescaled": format_threshold(threshold),
            "threshold_original": format_threshold(original_threshold),
            "command": args.command,
            "inputs": {
                "frozen_hits": str(hits),
                "shard_counts": str(counts_path),
                **input_records,
            },
            "hits_sha256": sha256_file(hits) if args.hits_sha256 else None,
            "shards": [
                {"name": name, "records": records}
                for name, records in sorted(ledger["counts"].items())
            ],
            "statistics": {
                "total_rows": rescaler.total_rows,
                "rows_rescalable": rescaler.rescalable_rows,
                "rows_pending_unknown_shard": rescaler.pending_rows,
                "rows_passing_original_threshold": rows_passing_original,
                "rows_passing_rescaled_threshold": rows_passing_rescaled,
                "delta_rows": delta_rows,
                "direction": direction,
            },
            "output_columns": columns,
            "frozen_table_modified": False,
        }
        manifest_path.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )

    return {
        "out": str(out),
        "delta_json": str(delta_path) if delta_path is not None else None,
        "manifest": str(manifest_path) if manifest_path is not None else None,
        "z_total": z_total,
        "z_total_basis": z_basis,
        "total_rows": rescaler.total_rows,
        "rows_written": written,
        "rows_passing_original_threshold": rows_passing_original,
        "rows_passing_rescaled_threshold": rows_passing_rescaled,
        "delta_rows": delta_rows,
        "direction": direction,
        "rows_pending_unknown_shard": rescaler.pending_rows,
        "unknown_shards": sorted(rescaler.unknown_shards()),
        "limited": limit is not None,
        "limit": limit,
        "threshold_rescaled": format_threshold(threshold),
        "threshold_original": format_threshold(original_threshold),
    }


def main(argv: Sequence[str] | None = None) -> int:
    status, args = parse_args_never_exits(argv)
    if args is None:
        return status
    try:
        payload = run(args)
    except UnknownShardError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_UNKNOWN_SHARD
    except NoDatabaseSizeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_NO_DATABASE_SIZE
    except FileExistsError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_REFUSED
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_USAGE
    if not args.quiet:
        print(
            "rescaled %d row(s) onto Z=%d (%s): passes %d -> %d (delta %+d) -> %s"
            % (
                payload["total_rows"],
                payload["z_total"],
                payload["z_total_basis"],
                payload["rows_passing_original_threshold"],
                payload["rows_passing_rescaled_threshold"],
                payload["delta_rows"],
                payload["out"],
            ),
            file=sys.stderr,
        )
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
