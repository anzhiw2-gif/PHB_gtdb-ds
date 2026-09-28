#!/usr/bin/env python3
"""Count FASTA records per frozen scan shard and measure the full-library ``Z``.

Why this module exists (task F7, scan-13 ``-Z`` reconciliation)
--------------------------------------------------------------
The frozen GTDB full-library scan (``runs/20260901_formal_frozen_scan_13``) ran
``hmmsearch`` **without** ``-Z`` and without ``--domZ``:

    hmmsearch --tblout "$BUILD/$base.tbl" --domtblout "$BUILD/$base.dom" \
              -E 1e-5 --cpu 1 "$hmm" "$filtered"

HMMER therefore used each shard's *own* sequence count as ``Z``, so every shard
reported E-values on a different scale (``shard_0001`` = 4,694,121 records vs
``shard_0051`` = 2,923,820, a 1.61x spread).  The frozen manifest records no
per-shard count at all, and the long-standing full-library total was only an
approximation ("~2.92e8").  This module closes both gaps by *measuring* them:

* one record count per shard, and
* the sum of those counts = the measured full-library ``Z``.

Counting semantics: ``grep -c '^>'`` equivalence
-----------------------------------------------
The project's established convention for a shard record count is
``grep -c '^>' <shard>.faa`` (used for the two counts already on record,
shard_0001 = 4,694,121 and shard_0051 = 2,923,820).  To keep a count produced
here **comparable with every historical count**, the default mode is defined as
exactly that operation and nothing more:

* a record is a **line whose first byte is ``>``** — i.e. a line matching ``^>``;
* ``grep`` splits lines on ``\\n`` only, so a ``\\r\\n`` file counts identically
  (the trailing ``\\r`` stays part of the line);
* a file whose last line has no trailing newline still contributes that line —
  just like ``grep``;
* a ``>`` that is *not* the first byte of a line (for example inside a sequence
  line, ``MK>AV``) is **not** a header and is not counted;
* an empty file counts ``0``.

Nothing else is validated in the default mode, so it can never be stricter than
``grep``.  ``--strict`` adds validation *on top* (a sequence line before any
header, a control character or non-sequence character in a sequence line, an
empty header, an underivable header identifier, an over-long line).  Use
``--strict`` on the real 267 GB shards only when the extra pass is affordable:
strict mode raises on the first malformed line instead of reporting a possibly
wrong count.

Memory: shards are up to ~2.3 GB each and the whole library is ~267 GB, so
counting streams each file in fixed 1 MiB binary chunks and never materialises
the file.  The only growth allowed is a single line; a line longer than
``DEFAULT_MAX_LINE_BYTES`` (64 MiB) is treated as corruption and raises instead
of growing without bound.

Output contract
---------------
A TSV with columns ``name``, ``path``, ``bytes``, ``records``, ``status``,
``pending_reason`` (plus ``sha256`` when ``--sha256`` is passed) and a JSON
summary on stdout:

    {"shards": N, "total_records": Z, "database_size_Z": Z, "out": ...}

``database_size_Z`` is the sum of the per-shard counts and
``database_size_basis`` says so in words; the full record is written with
``--summary-json``.  ``--sha256`` is **off** by default because hashing the
whole library costs a second full read of 267 GB.

Resumability
------------
``--resume-from <tsv>`` reuses counts already present in an earlier table and
counts only the shards that are missing from it.  Reused rows keep their
``status``/``path``/``bytes``/``sha256`` and are labelled ``reused`` (or
``reused_unhashed`` when the requested ``--sha256`` column is not available in
the resume file), so which numbers were measured *now* is never ambiguous.
Rows in the resume file that are not part of this input set are reported as
``reused_not_in_input`` rather than silently dropped.

Write safety
------------
The output is written **after** every shard has been counted, so a failure
leaves no partial ledger.  An existing ``--out`` file is refused unless
``--force`` is given.  A non-empty existing *directory* passed as ``--out`` is
refused outright (the coordinator may only pass an empty directory or a file
path).

This module never writes to ``runs/``, ``results/`` or ``deploy/`` by itself;
the caller supplies ``--out`` inside a new dated run.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Sequence

#: Stream chunk.  1 MiB keeps a 267 GB library bounded while amortising syscalls.
DEFAULT_CHUNK_BYTES = 1024 * 1024

#: A single FASTA line (a wrapped sequence line, realistically <= 80 bytes) may
#: never exceed this.  Beyond it the file is corrupt and growing the buffer
#: would defeat the streaming bound, so it raises instead.
DEFAULT_MAX_LINE_BYTES = 64 * 1024 * 1024

#: The record marker, as in ``grep '^>'``.
HEADER_BYTE = b">"

#: Characters allowed on a sequence line in ``--strict`` mode.  The IUPAC amino
#: acid alphabet (including the ambiguity codes), nucleotide codes, the gap
#: characters ``-`` and ``.``, and ``*`` (stop).  ASCII-only, case-insensitive.
SEQUENCE_ALPHABET = frozenset(
    "ACDEFGHIKLMNPQRSTVWYBXZJUO*-.acdefghiklmnpqrstvwybxzjuo"
)

BASE_COLUMNS = ["name", "path", "bytes", "records", "status", "pending_reason"]
SHA_COLUMN = "sha256"
COUNTED = "counted"
REUSED = "reused"
REUSED_UNHASHED = "reused_unhashed"
#: A row that carries no count: excluded from Z and reported as pending.
NOT_COUNTED = "not_counted"
PENDING_NOT_COUNTED = "not_counted_by_choice"
#: A shard with zero records is a defect, not a measurement of zero.
PENDING_EMPTY = "zero_records_possible_truncation_or_wrong_input"

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_REFUSED = 4

#: Stated verbatim in the summary so the basis of ``database_size_Z`` is never
#: guessed from the number alone.
DATABASE_SIZE_BASIS = "sum of per-shard record counts (one count per frozen scan shard)"


# --------------------------------------------------------------------------
# pure helpers
# --------------------------------------------------------------------------
def count_records(text: Any, strict: bool = False) -> int:
    """Count FASTA records in ``text``, matching ``grep -c '^>'``.

    ``text`` is either a string (split on ``\\n``) or an iterable of lines
    (line terminators optional).  A record is a line whose **first character**
    is ``>``; the trailing ``\\r`` of a CRLF file, and a final line without a
    newline, change nothing.  A ``>`` inside a sequence line is not a record.

    With ``strict=True`` the same count is produced, but a malformed file
    raises ``ValueError`` instead of being counted silently.  The default
    (``strict=False``) is deliberately exactly ``grep``'s semantics.
    """
    if isinstance(text, str):
        lines: Iterable[str] = text.split("\n")
    elif isinstance(text, (bytes, bytearray)):
        raise ValueError("count_records takes text, not bytes; use count_fasta_file")
    else:
        lines = text

    records = 0
    for index, line in enumerate(lines, start=1):
        if line.startswith(">"):
            if strict:
                _validate_header(line[1:], index)
            records += 1
            continue
        if strict and line.strip():
            if records == 0:
                raise ValueError(
                    f"line {index}: sequence data before any FASTA header"
                )
            _validate_sequence_line(line, index)
    return records


def _validate_header(header: str, index: int) -> None:
    """Reject a header whose identifier cannot be derived (strict mode only)."""
    body = header.rstrip("\r\n")
    if not body.strip():
        raise ValueError(f"line {index}: FASTA header is empty")
    if "\x00" in body:
        raise ValueError(f"line {index}: FASTA header contains a NUL byte")
    identifier = body.split()[0]
    if not identifier:
        raise ValueError(f"line {index}: FASTA header has no identifier")


def _validate_sequence_line(line: str, index: int) -> None:
    """Reject a sequence line carrying control or non-sequence characters.

    A ``>`` on a sequence line (``MK>AV``) is one of those characters: a ``>``
    may only ever start a header line.  The default mode still counts such a
    line the way ``grep -c '^>'`` does — it is simply not a record.
    """
    body = line.rstrip("\r\n")
    if "\x00" in body:
        raise ValueError(f"line {index}: sequence line contains a NUL byte")
    if not body.isascii():
        raise ValueError(f"line {index}: sequence line is not ASCII")
    bad = sorted({char for char in body if char not in SEQUENCE_ALPHABET})
    if bad:
        preview = "".join(repr(char) for char in bad[:5])
        raise ValueError(
            f"line {index}: sequence line has an invalid character {preview}"
        )


def validate_record_count(value: Any, *, require_positive: bool = False) -> int:
    """Return ``value`` as a non-negative ``int`` or raise ``ValueError``.

    ``bool`` is rejected explicitly: ``True`` is an ``int`` in Python but is
    never a sequence count.  With ``require_positive=True`` a count of ``0`` is
    rejected as well — an empty shard is a defect, never a measurement, and must
    never be merged into a scale ledger as if it were one.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"record count must be an integer, got {value!r}")
    if value < 0:
        raise ValueError(f"record count must be non-negative, got {value!r}")
    if require_positive and value < 1:
        raise ValueError(f"record count must be at least 1, got {value!r}")
    return value


def merge_counts(
    existing: Mapping[str, Any] | None,
    new: Mapping[str, Any] | None,
) -> Dict[str, Any]:
    """Merge two ``{name: records}`` mappings without silently overwriting.

    Returns ``{"counts": {name: records}, "new": [names], "reused": [names]}``.
    An entry already present in ``existing`` keeps its **existing** value and is
    listed under ``reused`` — that is what makes ``--resume-from`` idempotent
    and leaves the already-published number untouched.  ``new`` lists the names
    whose value came from this call, in insertion order.
    """
    counts: Dict[str, int] = {}
    reused: List[str] = []
    added: List[str] = []
    for name, value in (existing or {}).items():
        counts[str(name)] = validate_record_count(value, require_positive=True)
        reused.append(str(name))
    for name, value in (new or {}).items():
        key = str(name)
        merged = validate_record_count(value, require_positive=True)
        if key in counts:
            continue  # never overwrite an already-published measurement
        counts[key] = merged
        added.append(key)
    return {"counts": counts, "new": added, "reused": reused}


def split_complete_lines(buffer: bytes) -> tuple[list[bytes], bytes]:
    """Split ``buffer`` into complete ``\\n``-terminated lines plus the remainder.

    Lines are split on ``b"\\n"`` only, exactly like ``grep``.  A CRLF file
    yields lines ending in ``\\r`` and therefore counts identically.
    """
    parts = buffer.split(b"\n")
    return parts[:-1], parts[-1]


def iter_lines_bounded(
    handle,
    *,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    max_line_bytes: int = DEFAULT_MAX_LINE_BYTES,
) -> Iterator[bytes]:
    """Yield ``\\n``-terminated lines (terminator stripped) from a binary handle.

    The handle is read in ``chunk_bytes`` chunks and at most one line is held in
    memory, so a 2.3 GB shard is never materialised.  A line longer than
    ``max_line_bytes`` raises ``ValueError`` rather than growing the buffer
    without bound.
    """
    if chunk_bytes < 1:
        raise ValueError(f"chunk_bytes must be positive, got {chunk_bytes!r}")
    if max_line_bytes < 1:
        raise ValueError(f"max_line_bytes must be positive, got {max_line_bytes!r}")
    buffer = b""
    while True:
        chunk = handle.read(chunk_bytes)
        if chunk:
            buffer += chunk
        lines, buffer = split_complete_lines(buffer)
        for line in lines:
            yield line
        if len(buffer) > max_line_bytes:
            raise ValueError(
                f"a line exceeds max_line_bytes={max_line_bytes}; the file looks "
                "unwrapped or corrupt (refusing to buffer it without bound)"
            )
        if not chunk:
            if buffer:
                yield buffer
            return


def count_headered_lines(
    handle,
    *,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    max_line_bytes: int = DEFAULT_MAX_LINE_BYTES,
) -> int:
    """Count ``^>`` lines in a binary stream without building every line.

    The default counting mode needs only the count, so it avoids materialising
    each line: a line starts with ``>`` exactly when the ``>`` follows a ``\\n``
    (or begins the file), so per chunk this is ``chunk.count(b"\\n>")`` plus a
    check of the chunk's first byte when the previous chunk ended a line.

    Every chunk is scanned exactly once and no byte is scanned twice, so a header
    that lands on a chunk boundary is counted exactly once.  The line-length
    guard of :func:`iter_lines_bounded` is applied here too: the running
    single-line length may not exceed ``max_line_bytes``.
    """
    if chunk_bytes < 1:
        raise ValueError(f"chunk_bytes must be positive, got {chunk_bytes!r}")
    if max_line_bytes < 1:
        raise ValueError(f"max_line_bytes must be positive, got {max_line_bytes!r}")
    records = 0
    at_line_start = True
    line_length = 0
    while True:
        chunk = handle.read(chunk_bytes)
        if not chunk:
            return records
        if at_line_start and chunk[0:1] == HEADER_BYTE:
            records += 1
        records += chunk.count(b"\n>")
        newline = chunk.rfind(b"\n")
        if newline == -1:
            line_length += len(chunk)
        else:
            # The tail after the last newline is the incomplete line so far.
            line_length = len(chunk) - newline - 1
        if line_length > max_line_bytes:
            raise ValueError(
                f"a line exceeds max_line_bytes={max_line_bytes}; the file looks "
                "unwrapped or corrupt (refusing to buffer it without bound)"
            )
        at_line_start = chunk.endswith(b"\n")


def sha256_file(path: Path) -> str:
    """Return the SHA-256 of ``path``, streamed in 1 MiB blocks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(DEFAULT_CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# file-level counting
# --------------------------------------------------------------------------
def _strict_scan_binary(path: Path, chunk_bytes: int, max_line_bytes: int) -> int:
    """Count and validate ``path`` from its bytes; reject non-UTF-8 (strict)."""
    with path.open("rb") as handle:
        return _count_binary_stream(
            handle, strict=True, chunk_bytes=chunk_bytes, max_line_bytes=max_line_bytes
        )


def _count_binary_stream(
    handle, *, strict: bool, chunk_bytes: int, max_line_bytes: int
) -> int:
    """Count ``^>`` lines in a binary stream, optionally validating each line.

    Default mode takes the vectorised path (no per-line objects); ``strict``
    walks the lines so every malformed line can be named.
    """
    if not strict:
        return count_headered_lines(
            handle, chunk_bytes=chunk_bytes, max_line_bytes=max_line_bytes
        )
    records = 0
    line_number = 0
    for raw in iter_lines_bounded(
        handle, chunk_bytes=chunk_bytes, max_line_bytes=max_line_bytes
    ):
        line_number += 1
        if not raw:
            continue
        if raw.startswith(HEADER_BYTE):
            try:
                header = raw.decode("utf-8")
            except UnicodeDecodeError as error:
                raise ValueError(
                    f"{handle_name(handle)} line {line_number}: FASTA header "
                    f"is not valid UTF-8 ({error})"
                ) from error
            _validate_header(header, line_number)
            records += 1
            continue
        if not raw.strip():
            continue
        if records == 0:
            raise ValueError(
                f"{handle_name(handle)} line {line_number}: sequence data "
                "before any FASTA header"
            )
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError(
                f"{handle_name(handle)} line {line_number}: sequence line is "
                f"not valid UTF-8 ({error})"
            ) from error
        _validate_sequence_line(text, line_number)
    return records


def handle_name(handle) -> str:
    """Best-effort name of an open handle, for error messages."""
    name = getattr(handle, "name", None)
    return str(name) if name is not None else "<stream>"


def _validate_input_path(path: Path) -> Path:
    """Return ``path`` resolved, or raise ``ValueError`` when it is unusable."""
    if not path.exists():
        raise ValueError(f"input FASTA does not exist: {path}")
    if path.is_dir():
        raise ValueError(f"input FASTA is a directory, not a file: {path}")
    if not path.is_file():
        raise ValueError(f"input FASTA is not a regular file: {path}")
    try:
        with path.open("rb"):
            pass
    except OSError as error:
        raise ValueError(f"input FASTA is not readable: {path} ({error})") from error
    return path.resolve()


def count_fasta_file(
    path: Path | str,
    *,
    strict: bool = False,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    max_line_bytes: int = DEFAULT_MAX_LINE_BYTES,
    with_sha256: bool = False,
) -> Dict[str, Any]:
    """Count records in one FASTA file; return ``{records, bytes[, sha256]}``.

    Streams the file in fixed ``chunk_bytes`` chunks (never a whole-file read).
    ``strict`` adds the validations documented in the module docstring and
    rejects non-UTF-8 bytes.  ``bytes`` is the file size actually read.
    """
    path = Path(path)
    _validate_input_path(path)
    seen = 0
    with path.open("rb") as handle:
        records = _count_binary_stream(
            handle,
            strict=strict,
            chunk_bytes=chunk_bytes,
            max_line_bytes=max_line_bytes,
        )
        seen = handle.tell()
    result: Dict[str, Any] = {"records": records, "bytes": seen}
    if with_sha256:
        result["sha256"] = sha256_file(path)
    return result


# --------------------------------------------------------------------------
# counting-table IO
# --------------------------------------------------------------------------
def counts_columns(with_sha256: bool) -> List[str]:
    """Return the TSV column order for a record-count table."""
    columns = list(BASE_COLUMNS)
    if with_sha256:
        columns.insert(4, SHA_COLUMN)
    return columns


def _open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def parse_record_count(value: Any, field: str = "records") -> int:
    """Parse a record count from a TSV cell (``"123"`` / ``123``)."""
    if isinstance(value, bool):
        raise ValueError(f"{field} must be an integer, got {value!r}")
    if isinstance(value, int):
        return validate_record_count(value)
    text = ("" if value is None else str(value)).strip()
    if not text:
        raise ValueError(f"{field} is empty; a blank count is never 0")
    try:
        number = int(text, 10)
    except ValueError as error:
        raise ValueError(f"{field} is not an integer: {value!r}") from error
    return validate_record_count(number)


def read_counts_table(path: Path | str) -> List[Dict[str, str]]:
    """Read a record-count TSV (plain or ``.gz``) into a list of row dicts.

    Raises ``ValueError`` for a missing file or a missing required column: a
    scale ledger without counts must fail closed rather than default to zero.
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"record-count table does not exist: {path}")
    if path.is_dir():
        raise ValueError(f"record-count table is a directory: {path}")
    rows: List[Dict[str, str]] = []
    with _open_text(path, "rt") as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"record-count table is empty: {path}")
        header = [cell.strip() for cell in header_line.rstrip("\r\n").split("\t")]
        for required in ("name", "records"):
            if required not in header:
                raise ValueError(
                    f"record-count table {path} is missing the {required!r} column"
                )
        for number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) != len(header):
                raise ValueError(
                    f"{path} line {number}: expected {len(header)} columns, "
                    f"got {len(cells)}"
                )
            rows.append(dict(zip(header, cells)))
    return rows


def counts_by_name(rows: Iterable[Mapping[str, str]]) -> Dict[str, Any]:
    """Index count rows by ``name`` with duplicates rejected.

    Returns ``{"counts": {name: records}, "rows": {name: row}, "unreadable":
    {name: row}}``.  A duplicate ``name`` is an error: two different counts for
    one shard would make the rescaling scale ambiguous.  Rows whose status is
    ``not_counted`` are collected separately (as ``unreadable``) and are never
    treated as a count of zero.
    """
    counts: Dict[str, int] = {}
    rows_by_name: Dict[str, Dict[str, str]] = {}
    unreadable: Dict[str, Dict[str, str]] = {}
    for row in rows:
        name = (row.get("name") or "").strip()
        if not name:
            raise ValueError("record-count table has a row without a name")
        if name in rows_by_name:
            raise ValueError(f"record-count table has a duplicate shard {name!r}")
        rows_by_name[name] = dict(row)
        status = (row.get("status") or COUNTED).strip()
        if status == NOT_COUNTED:
            unreadable[name] = dict(row)
            continue
        counts[name] = parse_record_count(row.get("records"), f"records for {name}")
    return {"counts": counts, "rows": rows_by_name, "unreadable": unreadable}


# --------------------------------------------------------------------------
# shard discovery
# --------------------------------------------------------------------------
def discover_shard_files(shard_dir: Path | str) -> List[Path]:
    """Return the ``*.faa`` files of ``shard_dir``, sorted deterministically.

    ``shard_0001.faa`` sorts before ``shard_0002.faa``; the sort is by name so
    the emitted table order does not depend on the filesystem or on ``--jobs``.
    """
    directory = Path(shard_dir)
    if not directory.exists():
        raise ValueError(f"shard directory does not exist: {directory}")
    if not directory.is_dir():
        raise ValueError(f"shard directory is not a directory: {directory}")
    files = [path for path in directory.iterdir() if path.is_file() and path.suffix == ".faa"]
    if not files:
        raise ValueError(
            f"shard directory contains no FASTA (*.faa) files: {directory}"
        )
    return sorted(files, key=lambda path: (path.name, str(path)))


def shard_key(path: Path | str) -> str:
    """Return the shard key of a FASTA path: its file name without ``.faa``.

    The frozen ``hits_all.tsv`` identifies a shard as ``shard_0051`` (no
    extension); the counting table names the file ``shard_0051.faa``.  This
    helper is the single place that bridges the two vocabularies.
    """
    name = Path(path).name
    if name.endswith(".faa"):
        name = name[: -len(".faa")]
    elif name.endswith(".faa.gz"):
        name = name[: -len(".faa.gz")]
    return name


def _display_path(path: Path, *, absolute: bool = False) -> str:
    """Render a path for the ledger: relative when possible, absolute otherwise.

    With ``absolute=True`` the ledger records the server-side path verbatim, so
    it stays self-describing when read from another working directory; the
    default records a path relative to the working directory, which keeps the
    table portable between the server and a local checkout.  Only a path on
    another drive (``relpath`` is impossible) falls back to absolute.
    """
    if absolute:
        return str(path)
    try:
        return os.path.relpath(path, Path.cwd()).replace(os.sep, "/")
    except ValueError:  # different drive on Windows
        return str(path)


# --------------------------------------------------------------------------
# the run
# --------------------------------------------------------------------------
def build_arg_namespace(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count FASTA records per scan shard (grep -c '^>' semantics) and "
            "report the measured full-library Z."
        )
    )
    parser.add_argument(
        "--shard-dir",
        type=Path,
        default=None,
        help="directory of *.faa shards (all of them are counted)",
    )
    parser.add_argument(
        "--file",
        dest="files",
        action="append",
        default=[],
        help="explicit FASTA file; repeatable. Overrides --shard-dir for the same path",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="record-count TSV to write (a new path; an existing non-empty "
        "directory is refused)",
    )
    parser.add_argument(
        "--resume-from",
        type=Path,
        default=None,
        help="earlier record-count TSV; its counts are reused and those shards "
        "are not re-counted",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="shards counted concurrently (caller-supplied; no default fan-out)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="also validate structure (sequence before a header, control "
        "characters, over-long lines) instead of matching grep exactly",
    )
    parser.add_argument(
        "--sha256",
        action="store_true",
        help="also record each shard's SHA-256 (a second full read; default off)",
    )
    parser.add_argument(
        "--absolute-paths",
        action="store_true",
        help="record absolute shard paths in the ledger (default: paths relative "
        "to the working directory)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing --out file",
    )
    parser.add_argument(
        "--gzip",
        action="store_true",
        help="gzip the output table (name it *.tsv.gz); readers accept both",
    )
    parser.add_argument(
        "--chunk-bytes",
        type=int,
        default=DEFAULT_CHUNK_BYTES,
        help=f"stream chunk size (default {DEFAULT_CHUNK_BYTES})",
    )
    parser.add_argument(
        "--max-line-bytes",
        type=int,
        default=DEFAULT_MAX_LINE_BYTES,
        help=f"reject a longer single line (default {DEFAULT_MAX_LINE_BYTES})",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="also write the full JSON summary to this new file",
    )
    parser.add_argument(
        "--command",
        default="pending",
        help="the exact command this table came from, recorded verbatim",
    )
    parser.add_argument(
        "--input-record",
        action="append",
        default=[],
        help="extra KEY=PATH binding to record in the summary; repeatable",
    )
    parser.add_argument("--quiet", action="store_true", help="suppress the progress line")
    return parser.parse_args(argv)


def validate_jobs(jobs: Any) -> int:
    """Return ``jobs`` as a positive ``int`` (the caller decides the value)."""
    if isinstance(jobs, bool) or not isinstance(jobs, int) or jobs < 1:
        raise ValueError(f"--jobs must be a positive integer, got {jobs!r}")
    return jobs


def _resolve_inputs(args: argparse.Namespace) -> List[Path]:
    """Return the de-duplicated input FASTA list in deterministic order."""
    explicit: List[Path] = []
    for raw in args.files:
        explicit.append(Path(raw))
    if args.shard_dir is None and not explicit:
        raise ValueError("no input: pass --shard-dir and/or --file")
    discovered: List[Path] = (
        discover_shard_files(args.shard_dir) if args.shard_dir is not None else []
    )

    ordered: List[Path] = []
    seen: set[str] = set()
    for path in list(discovered) + explicit:
        _validate_input_path(path)
        key = str(path.resolve()).lower()
        if key in seen:
            continue
        seen.add(key)
        ordered.append(path)
    return ordered


def _check_output_target(out: Path, force: bool) -> None:
    """Refuse an unusable ``--out`` before any work is done."""
    if out.exists() and out.is_dir():
        contents = list(out.iterdir())
        if contents:
            raise ValueError(
                f"--out is an existing non-empty directory: {out} "
                f"({len(contents)} entries). Pass a NEW file path."
            )
        return
    if out.exists() and not force:
        raise FileExistsError(
            f"--out already exists: {out}. Pass --force to overwrite it."
        )
    if not out.parent.exists():
        raise ValueError(f"--out parent directory does not exist: {out.parent}")


def _write_rows(
    out: Path,
    columns: Sequence[str],
    rows: Sequence[Mapping[str, Any]],
    *,
    gzip_output: bool,
) -> None:
    with _open_text(out, "wt") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write(
                "\t".join("" if row.get(c) is None else str(row.get(c)) for c in columns)
                + "\n"
            )


def _count_one(
    path: Path,
    *,
    strict: bool,
    with_sha256: bool,
    chunk_bytes: int,
    max_line_bytes: int,
    absolute_paths: bool = False,
) -> Dict[str, Any]:
    result = count_fasta_file(
        path,
        strict=strict,
        chunk_bytes=chunk_bytes,
        max_line_bytes=max_line_bytes,
        with_sha256=with_sha256,
    )
    records = int(result["records"])
    row: Dict[str, Any] = {
        "name": shard_key(path),
        "path": _display_path(path, absolute=absolute_paths),
        "bytes": result["bytes"],
        "records": records,
        "status": COUNTED,
        "pending_reason": "",
    }
    if with_sha256:
        row[SHA_COLUMN] = result["sha256"]
    if records < 1:
        # A zero-record shard is a defect (truncation, wrong filter, wrong path),
        # never a measurement.  It is still reported, but never summed into Z.
        row["status"] = NOT_COUNTED
        row["pending_reason"] = PENDING_EMPTY
    return row


def _reuse_row(
    name: str,
    previous: Mapping[str, str],
    *,
    with_sha256: bool,
    default_path: str | None,
) -> Dict[str, Any]:
    """Convert a resume-file row into an output row labelled ``reused``."""
    row: Dict[str, Any] = {
        "name": name,
        "path": (previous.get("path") or default_path or "").strip(),
        "bytes": (previous.get("bytes") or "").strip(),
        "records": previous.get("records", ""),
        "pending_reason": "",
    }
    if with_sha256 and (previous.get(SHA_COLUMN) or "").strip():
        row[SHA_COLUMN] = previous[SHA_COLUMN].strip()
        row["status"] = REUSED
    elif with_sha256:
        row[SHA_COLUMN] = ""
        row["status"] = REUSED_UNHASHED
    else:
        row["status"] = REUSED
    return row


def run(args: argparse.Namespace) -> Dict[str, Any]:
    """Count the inputs and write the table; return the summary payload.

    Raises ``ValueError`` for an unusable path or argument (missing shard dir,
    missing explicit file, bad ``--jobs``, unusable ``--out``) and
    ``FileExistsError`` when ``--out`` exists without ``--force``.
    """
    jobs = validate_jobs(args.jobs)
    if args.chunk_bytes < 1:
        raise ValueError(f"--chunk-bytes must be positive, got {args.chunk_bytes!r}")
    if args.max_line_bytes < 1:
        raise ValueError(
            f"--max-line-bytes must be positive, got {args.max_line_bytes!r}"
        )

    inputs = _resolve_inputs(args)
    if not inputs:
        raise ValueError("no FASTA inputs to count")

    out = Path(args.out)
    _check_output_target(out, args.force)
    if args.summary_json is not None and Path(args.summary_json).exists() and not args.force:
        raise FileExistsError(
            f"--summary-json already exists: {args.summary_json}. Pass --force."
        )
    if args.gzip and out.suffix != ".gz":
        raise ValueError(f"--gzip requires a .gz --out name, got {out}")

    previous = {"counts": {}, "rows": {}, "unreadable": {}}
    if args.resume_from is not None:
        previous = counts_by_name(read_counts_table(args.resume_from))

    rows: List[Dict[str, Any]] = []
    pending: List[Dict[str, Any]] = []
    reused_count = 0
    recount_note: List[Dict[str, str]] = []
    for path in inputs:
        key = shard_key(path)
        entry = previous["rows"].get(key)
        if entry is not None and str(entry.get("records", "")).strip():
            rows.append(
                _reuse_row(
                    key,
                    entry,
                    with_sha256=args.sha256,
                    default_path=_display_path(
                        path, absolute=bool(args.absolute_paths)
                    ),
                )
            )
            reused_count += 1
        elif entry is not None:
            # Present but without a usable count (an earlier unreadable row):
            # count it now and say so, instead of reusing an empty cell.
            recount_note.append(
                {
                    "name": key,
                    "reason": (
                        "earlier row had status=%s and no count"
                        % (entry.get("status") or "")
                    ),
                }
            )
        elif key in previous["unreadable"]:
            reason = (previous["unreadable"][key].get("pending_reason") or "").strip()
            rows.append(
                {
                    "name": key,
                    "path": _display_path(path, absolute=bool(args.absolute_paths)),
                    "bytes": (previous["unreadable"][key].get("bytes") or "").strip(),
                    "records": "",
                    "status": NOT_COUNTED,
                    "pending_reason": reason or PENDING_NOT_COUNTED,
                    **({SHA_COLUMN: ""} if args.sha256 else {}),
                }
            )
            pending.append(
                {
                    "name": key,
                    "reason": reason or PENDING_NOT_COUNTED,
                    "records": None,
                }
            )

    to_count = [
        path
        for path in inputs
        if not (
            shard_key(path) in previous["rows"]
            and str(previous["rows"][shard_key(path)].get("records", "")).strip()
        )
        and shard_key(path) not in previous["unreadable"]
    ]
    counted_rows: List[Dict[str, Any]] = []
    if to_count:
        if jobs == 1:
            results = [
                _count_one(
                    path,
                    strict=args.strict,
                    with_sha256=args.sha256,
                    chunk_bytes=args.chunk_bytes,
                    max_line_bytes=args.max_line_bytes,
                    absolute_paths=bool(args.absolute_paths),
                )
                for path in to_count
            ]
        else:
            # Threads, not processes: this is I/O and byte counting, and the
            # interpreter stays in one process so no pickling of paths.
            with ThreadPoolExecutor(max_workers=jobs) as pool:
                results = list(
                    pool.map(
                        lambda path: _count_one(
                            path,
                            strict=args.strict,
                            with_sha256=args.sha256,
                            chunk_bytes=args.chunk_bytes,
                            max_line_bytes=args.max_line_bytes,
                            absolute_paths=bool(args.absolute_paths),
                        ),
                        to_count,
                    )
                )
        by_name = {row["name"]: row for row in results}
        for path in to_count:
            counted_rows.append(by_name[shard_key(path)])

    rows.extend(counted_rows)
    rows.sort(key=lambda row: (row["name"], str(row.get("path", ""))))
    unreadable = [row for row in rows if row.get("status") == NOT_COUNTED]
    for row in unreadable:
        if not any(entry["name"] == row["name"] for entry in pending):
            pending.append(
                {
                    "name": row["name"],
                    "reason": row.get("pending_reason") or PENDING_NOT_COUNTED,
                    "records": None,
                }
            )

    total = sum(int(row["records"]) for row in rows if row.get("status") != NOT_COUNTED)
    shards = len(rows)
    complete = not pending
    pending_reason = (
        ""
        if complete
        else "not_counted/unreadable shards present (%d): %s"
        % (len(pending), ", ".join(entry["name"] for entry in pending))
    )

    columns = counts_columns(args.sha256)
    _write_rows(out, columns, rows, gzip_output=args.gzip)

    reused_not_in_input = sorted(
        name
        for name in previous["counts"]
        if name not in {row["name"] for row in rows}
    )

    input_records: Dict[str, str] = {}
    for record in args.input_record:
        if "=" in record:
            key, value = record.split("=", 1)
            input_records[key.strip()] = value.strip()

    summary: Dict[str, Any] = {
        "shards": shards,
        "total_records": total,
        "database_size_Z": total,
        "database_size_basis": DATABASE_SIZE_BASIS,
        "out": str(out),
        "counted": len(counted_rows),
        "reused": reused_count,
        "pending": [entry["name"] for entry in pending],
        "pending_reason": pending_reason,
        "complete": complete,
        "status_counts": _status_counts(rows),
        "sha256_recorded": bool(args.sha256),
        "strict": bool(args.strict),
        "jobs": jobs,
        "reused_not_in_input": reused_not_in_input,
        "recounted_over_empty_resume_rows": recount_note,
        "resume_from": str(args.resume_from) if args.resume_from else None,
        "command": args.command,
        "inputs": input_records,
        "shards_detail": [
            {
                "name": row["name"],
                "records": (
                    None if row.get("status") == NOT_COUNTED else int(row["records"])
                ),
                "status": row.get("status", ""),
                "bytes": row.get("bytes", ""),
                "path": row.get("path", ""),
            }
            for row in rows
        ],
    }

    if args.summary_json is not None:
        path = Path(args.summary_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
            encoding="utf-8",
        )
    return summary


def _status_counts(rows: Iterable[Mapping[str, Any]]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for row in rows:
        status = str(row.get("status") or "")
        counts[status] = counts.get(status, 0) + 1
    return counts


def parse_args_never_exits(argv: Sequence[str] | None = None) -> tuple[int, Any]:
    """Parse ``argv``; return ``(status, namespace)`` instead of exiting.

    ``argparse`` raises ``SystemExit`` for an unusable argv.  ``main`` is
    documented to *return* its status, so it catches that here and reduces it to
    an exit code; the usage message argparse already wrote to stderr is kept.
    """
    try:
        return EXIT_OK, build_arg_namespace(argv)
    except SystemExit as exit_error:
        code = exit_error.code
        return (code if isinstance(code, int) else EXIT_USAGE), None


def main(argv: Sequence[str] | None = None) -> int:
    status, args = parse_args_never_exits(argv)
    if args is None:
        return status
    try:
        summary = run(args)
    except FileExistsError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_REFUSED
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return EXIT_USAGE
    if not args.quiet:
        print(
            "counted %d shard(s): counted=%d reused=%d total_records=%d -> %s"
            % (
                summary["shards"],
                summary["counted"],
                summary["reused"],
                summary["total_records"],
                summary["out"],
            ),
            file=sys.stderr,
        )
    payload = {
        key: summary[key]
        for key in (
            "shards",
            "total_records",
            "database_size_Z",
            "database_size_basis",
            "out",
            "counted",
            "reused",
            "pending",
            "pending_reason",
            "complete",
            "reused_not_in_input",
        )
    }
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
