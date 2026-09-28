"""Task F7 — the authoritative per-shard FASTA record counts and full-library ``Z``.

The frozen GTDB scan 13 ran ``hmmsearch`` **without** ``-Z``, so HMMER used each
shard's own sequence count as ``Z`` and every shard reported E-values on its own
scale (``shard_0001`` = 4,694,121 vs ``shard_0051`` = 2,923,820 records).  The
frozen manifest records no per-shard count at all, so the reconciliation run
begins here: count ``^>`` lines in every frozen shard, sum them, and that sum is
the measured full-library ``database_size_Z``.

These tests are the failing test for ``pipeline/scripts/count_fasta_records.py``.
They pin:

* the counting semantics — byte-for-byte the same convention as the project's
  established ``grep -c '^>'``, so a count produced here is comparable with every
  historical shard count that was taken with ``grep``;
* the edge cases that decide whether that equivalence actually holds (empty
  file, CRLF, header-only file, a file whose last line has no newline, and a
  ``>`` that appears *inside* a sequence line and is therefore not a header);
* ``--strict`` validation (a sequence line before any header, binary/garbage
  lines, an over-long line) — off by default, so the default is never stricter
  than ``grep``;
* resumability/idempotence (``--resume-from`` reuses counts and skips work) and
  the refusal to overwrite an existing output file without ``--force``;
* ``main(argv=None) -> int`` exit codes.  A missing input path raises
  ``ValueError`` from ``run()``; ``main()`` turns it into a non-zero exit.

Nothing here touches ``runs/``, ``results/`` or ``deploy/``: every fixture is a
synthetic file inside a ``tempfile`` directory.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import threading
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "pipeline" / "scripts"
MODULE_PATH = SCRIPTS / "count_fasta_records.py"

#: Must match the default in the script; only a small synthetic fixture is
#: written at this size, never a real shard.
DEFAULT_MAX_LINE_BYTES = 64 * 1024 * 1024


def load_module():
    """Import ``pipeline/scripts/count_fasta_records.py`` by path.

    ``pipeline/scripts`` is not a package (it is a flat directory of standalone
    entrypoints, run as ``python pipeline/scripts/<name>.py``), so the module is
    loaded by file path exactly as the shell wrappers do it.
    """
    if not MODULE_PATH.is_file():
        raise AssertionError(f"missing record-counting script: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location("count_fasta_records_f7", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("count_fasta_records_f7", module)
    spec.loader.exec_module(module)
    return module


class _RecordingFile:
    """Wrap a binary handle and record the size of every ``read`` request.

    ``read()`` called with no argument would slurp the whole file, so a bounded
    stream must never do it.  The wrapper turns "did this stream the file?"
    into an assertion on observed read sizes instead of a guess.
    """

    def __init__(self, handle, sizes):
        self._handle = handle
        self._sizes = sizes

    def read(self, *args, **kwargs):
        self._sizes.append(args[0] if args else None)
        return self._handle.read(*args, **kwargs)

    def __iter__(self):
        return iter(self._handle)

    def __getattr__(self, name):
        return getattr(self._handle, name)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._handle.close()
        return False


def _record_block(name: str, sequence: str) -> bytes:
    """Return one FASTA record as bytes, in the project's ``grep -c '^>'`` shape."""
    return f">{name}\n{sequence}\n".encode("ascii")


class CountRecordsPureTests(unittest.TestCase):
    """``count_records`` is the pure core: text in, count out."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_empty_file_counts_zero(self):
        self.assertEqual(0, self.module.count_records(""))

    def test_counts_header_lines_only(self):
        text = ">a\nMKAV\n>b\nMK\n"
        self.assertEqual(2, self.module.count_records(text))

    def test_last_line_without_newline_still_counts(self):
        self.assertEqual(1, self.module.count_records(">only_header"))

    def test_crlf_headers_count(self):
        text = ">a\r\nMKAV\r\n>b\r\nMK\r\n"
        self.assertEqual(2, self.module.count_records(text))

    def test_headers_only_file_counts_every_header(self):
        self.assertEqual(3, self.module.count_records(">a\n>b\n>c\n"))

    def test_gt_inside_a_sequence_line_is_not_a_header(self):
        """``grep -c '^>'`` anchors the ``>``; a mid-line ``>`` is not a record."""
        text = ">a\nMK>AV\n>b\nMK\n"
        self.assertEqual(2, self.module.count_records(text))

    def test_gt_after_leading_whitespace_is_not_a_header(self):
        self.assertEqual(1, self.module.count_records(">a\nMK\n  >b\nMK\n"))

    def test_blank_and_comment_lines_are_not_records(self):
        self.assertEqual(1, self.module.count_records("\n>a\nMK\n\n#note\n"))

    def test_accepts_an_iterable_of_lines_without_line_endings(self):
        self.assertEqual(2, self.module.count_records([">a", "MK", ">b", "MK"]))

    def test_strict_rejects_a_sequence_line_before_any_header(self):
        with self.assertRaises(ValueError) as caught:
            self.module.count_records("MKAV\n>a\nMK\n", strict=True)
        self.assertIn("before any FASTA header", str(caught.exception))

    def test_strict_accepts_a_well_formed_file(self):
        text = ">a desc here\nMKAV\n>b\nMK\n"
        self.assertEqual(2, self.module.count_records(text, strict=True))

    def test_strict_rejects_a_mid_sequence_gt_line(self):
        with self.assertRaises(ValueError) as caught:
            self.module.count_records(">a\nMK>AV\n", strict=True)
        self.assertIn(">", str(caught.exception))

    def test_strict_rejects_garbage_characters_in_a_sequence_line(self):
        with self.assertRaises(ValueError) as caught:
            self.module.count_records(">a\nMK\x01AV\n", strict=True)
        self.assertIn("invalid character", str(caught.exception))

    def test_strict_rejects_a_non_ascii_sequence_line(self):
        with self.assertRaises(ValueError) as caught:
            self.module.count_records(">a\nMKAVé\n", strict=True)
        self.assertIn("not ASCII", str(caught.exception))

    def test_default_mode_tolerates_garbage_like_grep_does(self):
        """Default must never be stricter than ``grep -c '^>'``."""
        self.assertEqual(1, self.module.count_records(">a\nMK\x01AV\n"))

    def test_strict_tolerates_the_iupac_and_gap_alphabet(self):
        text = ">a\nMK-X*BJZOU\n"
        self.assertEqual(1, self.module.count_records(text, strict=True))


class CountFastaFileTests(unittest.TestCase):
    """``count_fasta_file`` reads a real path and reports bytes as well."""

    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def _write(self, name: str, data: bytes) -> Path:
        path = self.dir / name
        path.write_bytes(data)
        return path

    def test_empty_file(self):
        path = self._write("empty.faa", b"")
        self.assertEqual({"records": 0, "bytes": 0}, self.module.count_fasta_file(path))

    def test_counts_and_byte_size_match_the_file(self):
        data = _record_block("a", "MK") + _record_block("b", "MKAV")
        path = self._write("shard.faa", data)
        result = self.module.count_fasta_file(path)
        self.assertEqual(2, result["records"])
        self.assertEqual(len(data), result["bytes"])
        self.assertEqual(path.stat().st_size, result["bytes"])

    def test_crlf_file(self):
        path = self._write("crlf.faa", b">a\r\nMK\r\n>b\r\nMK\r\n")
        self.assertEqual(2, self.module.count_fasta_file(path)["records"])

    def test_last_line_without_a_newline(self):
        path = self._write("tail.faa", b">a\nMK\n>b\nMK")
        self.assertEqual(2, self.module.count_fasta_file(path)["records"])

    def test_header_only_file(self):
        path = self._write("headers.faa", b">a\n>b\n")
        self.assertEqual(2, self.module.count_fasta_file(path)["records"])

    def test_mid_sequence_gt_is_not_counted(self):
        path = self._write("mid.faa", b">a\nMK>AV\n>b\nMK\n")
        self.assertEqual(2, self.module.count_fasta_file(path)["records"])

    def test_missing_path_raises_value_error(self):
        with self.assertRaises(ValueError) as caught:
            self.module.count_fasta_file(self.dir / "nope.faa")
        self.assertIn("nope.faa", str(caught.exception))

    def test_directory_path_raises_value_error(self):
        with self.assertRaises(ValueError):
            self.module.count_fasta_file(self.dir)

    def test_strict_rejects_binary_content(self):
        path = self._write("binary.faa", b">a\nMK\x00\x01\x02AV\n")
        with self.assertRaises(ValueError):
            self.module.count_fasta_file(path, strict=True)

    def test_strict_accepts_an_ordinary_file(self):
        path = self._write("ok.faa", _record_block("a", "MKAV"))
        self.assertEqual(1, self.module.count_fasta_file(path, strict=True)["records"])


class FastPathEquivalenceTests(unittest.TestCase):
    """The vectorised default path must agree with the line-by-line path.

    ``count_records`` (per line, and the only path strict mode uses) is the
    reference; ``count_headered_lines`` is the chunk-based fast path used by
    default on real shards.  They must never disagree — including when a header
    lands exactly on a chunk boundary, which is the only place the fast path's
    overlap window matters.
    """

    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def _both(self, data: bytes, chunk: int) -> tuple[int, int]:
        """Return (reference count, fast-path count) for the same bytes.

        The reference is ``count_records`` — the per-line implementation that
        strict mode uses — so the comparison is between two independent code
        paths, not between the fast path and itself.
        """
        path = self.dir / "case.faa"
        path.write_bytes(data)
        text = data.decode("utf-8")
        reference = self.module.count_records(text)
        with path.open("rb") as handle:
            fast = self.module.count_headered_lines(
                handle, chunk_bytes=chunk, max_line_bytes=1024 * 1024
            )
        # the file-level entry point must agree with the fast path too
        self.assertEqual(
            reference, self.module.count_fasta_file(path, chunk_bytes=chunk)["records"]
        )
        return reference, fast

    def test_agrees_on_the_edge_cases_at_every_chunk_size(self):
        cases = {
            "empty": b"",
            "header_only": b">a\n",
            "header_without_newline": b">a",
            "headers_only": b">a\n>b\n>c\n",
            "crlf": b">a\r\nMK\r\n>b\r\nMK\r\n",
            "mid_sequence_gt": b">a\nMK>AV\n>b\nMK\n",
            "gt_after_indent": b">a\nMK\n >b\nMK\n",
            "blank_lines": b"\n>a\nMK\n\n#note\n",
            "trailing_gt": b">a\nMK\n>",
            "long_records": b"".join(
                b">rec_%04d\n%s\n" % (i, b"M" * 79) for i in range(200)
            ),
        }
        for name, data in cases.items():
            for chunk in (1, 2, 3, 7, 16, 64, 1000):
                with self.subTest(case=name, chunk=chunk):
                    reference, fast = self._both(data, chunk)
                    self.assertEqual(reference, fast)

    def test_counts_a_header_exactly_on_a_chunk_boundary_once(self):
        # chunk = 4 puts the boundary so that ">b" starts a chunk, and the
        # preceding newline sits in the previous chunk.
        data = b">aa\n" + b"M" * 8 + b"\n>b\nMK\n"
        self.assertEqual(2, self._both(data, 4)[0])
        self.assertEqual(2, self._both(data, 4)[1])
        # also when the whole window is exactly one chunk
        for chunk in range(1, 20):
            with self.subTest(chunk=chunk):
                reference, fast = self._both(data, chunk)
                self.assertEqual(2, reference)
                self.assertEqual(reference, fast)

    def test_fast_path_still_enforces_the_line_budget(self):
        path = self.dir / "runaway.faa"
        path.write_bytes(b">a\n" + b"M" * (2 * 1024 * 1024) + b"\n")
        with self.assertRaises(ValueError) as caught:
            self.module.count_fasta_file(path, max_line_bytes=1024 * 1024)
        self.assertIn("max_line_bytes", str(caught.exception))

    def test_fast_path_reads_bounded_chunks(self):
        path = self.dir / "big.faa"
        with path.open("wb") as handle:
            for i in range(20_000):
                handle.write(b">rec_%05d\nMKAV\n" % i)
        sizes: list = []
        real_open = Path.open
        chunk = 8 * 1024

        def patched_open(self, *args, **kwargs):
            handle = real_open(self, *args, **kwargs)
            if Path(self) == path:
                return _RecordingFile(handle, sizes)
            return handle

        with mock.patch.object(Path, "open", patched_open):
            result = self.module.count_fasta_file(
                path, chunk_bytes=chunk, max_line_bytes=1024 * 1024
            )
        self.assertEqual(20_000, result["records"])
        for size in sizes:
            self.assertIsNotNone(size, "the fast path must never read() unsized")
            self.assertLessEqual(size, chunk)
        self.assertGreater(len(sizes), 1)


class LineStreamingTests(unittest.TestCase):
    """A 267 GB shard must stream: the reader never asks for the whole file."""
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.max_line_bytes = 64 * 1024

    def test_multi_chunk_file_reads_bounded_chunks(self):
        self.size = 2_000
        block = _record_block("contig_1|gene_1", "M" * 999)
        path = self.dir / "big.faa"
        with path.open("wb") as handle:
            for _ in range(self.size):
                handle.write(block)
        sizes: list = []
        real_open = Path.open
        chunk = 64 * 1024

        def patched_open(self, *args, **kwargs):
            handle = real_open(self, *args, **kwargs)
            if Path(self) == path:
                return _RecordingFile(handle, sizes)
            return handle

        # Patch ``pathlib.Path.open``: ``Path.open`` resolves ``io.open`` at call
        # time, so patching ``builtins.open`` would not intercept it.
        with mock.patch.object(Path, "open", patched_open):
            result = self.module.count_fasta_file(
                path, chunk_bytes=chunk, max_line_bytes=self.max_line_bytes
            )
        self.assertEqual(self.size, result["records"])
        self.assertTrue(sizes, "the reader never called read()")
        for size in sizes:
            self.assertIsNotNone(
                size, "read() without a size would load the whole shard into memory"
            )
            self.assertLessEqual(
                size,
                chunk,
                f"read({size}) exceeds the bounded chunk budget ({chunk})",
            )
        self.assertGreater(
            len(sizes), 1, "a multi-chunk file must take more than one read"
        )

    def test_over_long_line_raises_instead_of_growing_without_bound(self):
        path = self.dir / "runaway.faa"
        path.write_bytes(b">a\n" + b"M" * (3 * 1024 * 1024) + b"\n")
        with self.assertRaises(ValueError) as caught:
            self.module.count_fasta_file(path, max_line_bytes=1024 * 1024)
        self.assertIn("max_line_bytes", str(caught.exception))


class MergeCountsTests(unittest.TestCase):
    """``merge_counts`` is the resumable/idempotent-friendly core."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_new_entries_are_added_and_existing_ones_reused(self):
        merged = self.module.merge_counts(
            {"shard_0001": 10}, {"shard_0001": 999, "shard_0002": 20}
        )
        self.assertEqual(
            {"shard_0001": 10, "shard_0002": 20}, dict(merged["counts"])
        )
        self.assertEqual(["shard_0002"], list(merged["new"]))
        self.assertEqual(["shard_0001"], list(merged["reused"]))

    def test_merging_is_idempotent(self):
        first = self.module.merge_counts({}, {"a": 1, "b": 2})
        second = self.module.merge_counts(first["counts"], {"a": 1, "b": 2})
        self.assertEqual(dict(first["counts"]), dict(second["counts"]))
        self.assertEqual(["a", "b"], list(second["reused"]))
        self.assertEqual([], list(second["new"]))

    def test_non_positive_counts_are_rejected(self):
        for bad in (0, -1, "5", 1.5, True, None):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    self.module.merge_counts({}, {"a": bad})


class _CountFixture(unittest.TestCase):
    """Shared synthetic shard directory for the CLI-level tests."""

    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.shards = self.dir / "scan_shards"
        self.shards.mkdir(parents=True, exist_ok=True)
        # Keyed by the module's shard key (the file name without ".faa"), which
        # is what a row's "name" column carries.
        self.expected_records = {}
        self._add("shard_0001.faa", 3)
        self._add("shard_0002.faa", 1)
        self._add("shard_0003.faa", 2)
        self._add("shard_0004.faa", 4)
        self.total = sum(self.expected_records.values())

    def _add(self, name: str, records: int, sequence: str = "MKAV") -> None:
        path = self.shards / name
        payload = b"".join(
            _record_block(f"{name}|contig_{index}", sequence)
            for index in range(records)
        )
        path.write_bytes(payload)
        self.expected_records[path.stem] = records

    def _main(self, argv):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = self.module.main(argv)
            except SystemExit as exit_error:  # argparse rejects an unusable argv
                code = exit_error.code if isinstance(exit_error.code, int) else 1
        return code, stdout.getvalue(), stderr.getvalue()

    def _read_tsv(self, path: Path):
        lines = path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        rows = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]
        return header, rows


class CountCliTests(_CountFixture):
    """End-to-end behaviour of ``main(argv=None) -> int``."""

    def test_absolute_paths_flag_records_server_side_paths(self):
        relative_out = self.dir / "relative.tsv"
        absolute_out = self.dir / "absolute.tsv"
        self.assertEqual(
            0,
            self._main(["--shard-dir", str(self.shards), "--out", str(relative_out)])[0],
        )
        code, _, err = self._main(
            [
                "--shard-dir",
                str(self.shards),
                "--out",
                str(absolute_out),
                "--absolute-paths",
            ]
        )
        self.assertEqual(0, code, err)
        _, relative_rows = self._read_tsv(relative_out)
        _, absolute_rows = self._read_tsv(absolute_out)
        self.assertEqual(
            {str(self.shards / f"{r['name']}.faa") for r in absolute_rows},
            {r["path"] for r in absolute_rows},
        )
        # Both spellings must resolve to the same real file.  The default cell is
        # relative to the working directory, except when the shards live on
        # another drive, where ``relpath`` is impossible and absolute is the only
        # correct answer (a Windows temp dir on C: with the repo on D: is exactly
        # that case).
        for relative, absolute in zip(relative_rows, absolute_rows):
            self.assertEqual(relative["name"], absolute["name"])
            cell = Path(relative["path"])
            self.assertTrue(
                (Path.cwd() / cell).resolve() == Path(absolute["path"]).resolve()
            )
            if cell.is_absolute():
                self.assertEqual(Path(absolute["path"]), cell)
            else:
                self.assertEqual(
                    "..", cell.parts[0], "a relative cell must not be drive-bound"
                )

    def test_counts_a_shard_directory_and_prints_a_json_summary(self):
        out = self.dir / "shard_counts.tsv"
        code, stdout, stderr = self._main(
            ["--shard-dir", str(self.shards), "--out", str(out)]
        )
        self.assertEqual(0, code, stderr)
        payload = json.loads(stdout)
        self.assertEqual(4, payload["shards"])
        self.assertEqual(self.total, payload["total_records"])
        self.assertEqual(self.total, payload["database_size_Z"])
        self.assertEqual(str(out), payload["out"])
        self.assertEqual("", payload["pending_reason"])

        header, rows = self._read_tsv(out)
        self.assertEqual(
            ["name", "path", "bytes", "records", "status", "pending_reason"], header
        )
        self.assertEqual(
            ["shard_0001", "shard_0002", "shard_0003", "shard_0004"],
            [r["name"] for r in rows],
        )
        for row in rows:
            self.assertEqual(
                str(self.expected_records[row["name"]]), row["records"]
            )
            self.assertEqual("counted", row["status"])
            self.assertEqual("", row["pending_reason"])
        self.assertEqual(self.total, sum(int(r["records"]) for r in rows))

    def test_tsv_carries_the_counted_path_and_byte_size(self):
        out = self.dir / "shard_counts.tsv"
        code, _, err = self._main(
            ["--shard-dir", str(self.shards), "--out", str(out)]
        )
        self.assertEqual(0, code, err)
        _, rows = self._read_tsv(out)
        for row in rows:
            path = Path(row["path"])
            self.assertTrue(path.is_file())
            self.assertEqual(path.stat().st_size, int(row["bytes"]))
            self.assertEqual(
                self.module.count_records(path.read_text(encoding="utf-8")),
                int(row["records"]),
            )

    def test_sha256_column_is_absent_by_default_and_correct_when_requested(self):
        default_out = self.dir / "default.tsv"
        self.assertEqual(
            0, self._main(["--shard-dir", str(self.shards), "--out", str(default_out)])[0]
        )
        header, _ = self._read_tsv(default_out)
        self.assertNotIn("sha256", header)

        hashed_out = self.dir / "hashed.tsv"
        code, _, err = self._main(
            ["--shard-dir", str(self.shards), "--out", str(hashed_out), "--sha256"]
        )
        self.assertEqual(0, code, err)
        header, rows = self._read_tsv(hashed_out)
        self.assertIn("sha256", header)
        for row in rows:
            expected = hashlib.sha256(Path(row["path"]).read_bytes()).hexdigest()
            self.assertEqual(expected, row["sha256"])

    def test_explicit_files_win_over_the_shard_directory(self):
        out = self.dir / "mixed.tsv"
        explicit = self.shards / "shard_0001.faa"
        code, stdout, err = self._main(
            [
                "--shard-dir",
                str(self.shards),
                "--file",
                str(explicit),
                "--out",
                str(out),
            ]
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(4, payload["shards"])

    def test_duplicate_inputs_are_counted_once(self):
        out = self.dir / "dedup.tsv"
        explicit = str(self.shards / "shard_0001.faa")
        code, stdout, err = self._main(
            ["--file", explicit, "--file", explicit, "--out", str(out)]
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(1, payload["shards"])
        self.assertEqual(self.expected_records["shard_0001"], payload["total_records"])

    def test_output_is_written_only_after_a_successful_count(self):
        out = self.dir / "should_not_exist.tsv"
        bad = self.dir / "missing_shards"
        code, _, _ = self._main(["--shard-dir", str(bad), "--out", str(out)])
        self.assertNotEqual(0, code)
        self.assertFalse(out.exists())

    def test_resume_from_reuses_counts_and_skips_the_work(self):
        """Only shards absent from the resume table may be counted again.

        The prior table deliberately covers three of the four shards, so a pass
        that ignores ``--resume-from`` cannot produce the same result: it would
        recount the now-empty shard_0002 as 0 and shard_0003 as 1, changing the
        total.  The prior table is itself produced by the real counter, so it is
        also a round-trip of the writer and the reader.
        """
        prior = self.dir / "prior.tsv"
        code, _, err = self._main(
            [
                "--file",
                str(self.shards / "shard_0001.faa"),
                "--file",
                str(self.shards / "shard_0002.faa"),
                "--file",
                str(self.shards / "shard_0004.faa"),
                "--out",
                str(prior),
            ]
        )
        self.assertEqual(0, code, err)
        prior_bytes = prior.read_bytes()

        # Corrupt the two shapes a resumed run must not re-measure.
        (self.shards / "shard_0002.faa").write_bytes(b"")
        (self.shards / "shard_0003.faa").write_bytes(_record_block("x", "MK"))
        second = self.dir / "second.tsv"
        code, stdout, err = self._main(
            [
                "--shard-dir",
                str(self.shards),
                "--out",
                str(second),
                "--resume-from",
                str(prior),
            ]
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        # 3 + 1 + 4 reused, plus shard_0003 recounted as 1 -> 9, NOT the 10 a
        # blind recount would give (which would also make shard_0002 empty).
        self.assertEqual(9, payload["total_records"])
        self.assertEqual(3, payload["reused"])
        self.assertEqual(1, payload["counted"])
        self.assertEqual(4, payload["shards"])
        self.assertEqual([], payload["reused_not_in_input"])

        second_rows = {
            Path(row["path"]).name: row for row in self.module.read_counts_table(second)
        }
        # shard_0002 is now an EMPTY file: its unchanged count can only have come
        # from the resume table, which is exactly what reuse must guarantee.
        for name in ("shard_0001.faa", "shard_0002.faa", "shard_0004.faa"):
            self.assertEqual(
                self.expected_records[Path(name).stem],
                int(second_rows[name]["records"]),
                f"{name} must be reused from the resume file",
            )
            self.assertEqual("reused", second_rows[name]["status"])
        self.assertEqual(1, int(second_rows["shard_0003.faa"]["records"]))
        self.assertEqual("counted", second_rows["shard_0003.faa"]["status"])
        # the resumed table itself must not have been silently rewritten
        self.assertEqual(prior_bytes, prior.read_bytes())

    def test_resume_rows_outside_the_input_set_are_reported(self):
        """A resume row for a shard that is not an input is reported, not hidden."""
        carried = self.dir / "carried.tsv"
        carried.write_text(
            "name\tpath\tbytes\trecords\tstatus\tpending_reason\n"
            "shard_0001\tshards/shard_0001.faa\t10\t3\tcounted\t\n"
            "shard_0099\tshards/shard_0099.faa\t99\t42\tcounted\t\n",
            encoding="utf-8",
        )
        out = self.dir / "out.tsv"
        code, stdout, err = self._main(
            [
                "--file",
                str(self.shards / "shard_0001.faa"),
                "--out",
                str(out),
                "--resume-from",
                str(carried),
            ]
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(3, payload["total_records"])
        self.assertEqual(1, payload["reused"])
        self.assertEqual(["shard_0099"], payload["reused_not_in_input"])

    def test_zero_record_shard_is_counted_but_reported_as_not_measured(self):
        """An empty shard is a defect, never a measurement of 0 for the scale."""
        empty = self.dir / "empty.faa"
        empty.write_bytes(b"")
        out = self.dir / "out.tsv"
        code, stdout, err = self._main(
            ["--file", str(self.shards / "shard_0001.faa"), "--file", str(empty),
             "--out", str(out)]
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(self.expected_records["shard_0001"], payload["total_records"])
        self.assertEqual(["empty"], payload["pending"])
        self.assertFalse(payload["complete"])
        self.assertIn("empty", payload["pending_reason"])
        rows = {row["name"]: row for row in self.module.read_counts_table(out)}
        self.assertEqual("0", rows["empty"]["records"])
        self.assertNotEqual("counted", rows["empty"]["status"])

    def test_resume_from_missing_file_is_an_error(self):
        out = self.dir / "out.tsv"
        code, _, err = self._main(
            [
                "--shard-dir",
                str(self.shards),
                "--out",
                str(out),
                "--resume-from",
                str(self.dir / "absent.tsv"),
            ]
        )
        self.assertNotEqual(0, code)
        self.assertIn("absent.tsv", err)
        self.assertFalse(out.exists())

    def test_refuses_to_overwrite_an_existing_output(self):
        out = self.dir / "keep.tsv"
        out.write_text("PRECIOUS\n", encoding="utf-8")
        code, _, err = self._main(
            ["--shard-dir", str(self.shards), "--out", str(out)]
        )
        self.assertNotEqual(0, code)
        self.assertIn("--force", err)
        self.assertEqual("PRECIOUS\n", out.read_text(encoding="utf-8"))

    def test_force_overwrites(self):
        out = self.dir / "keep.tsv"
        out.write_text("PRECIOUS\n", encoding="utf-8")
        code, stdout, err = self._main(
            ["--shard-dir", str(self.shards), "--out", str(out), "--force"]
        )
        self.assertEqual(0, code, err)
        self.assertIn("shard_0001", out.read_text(encoding="utf-8"))

    def test_resume_in_place_is_recounted_without_data_loss(self):
        """``--resume-from X --out X`` is legitimate; ``--force`` acknowledges it."""
        ledger = self.dir / "ledger.tsv"
        self.assertEqual(
            0, self._main(["--shard-dir", str(self.shards), "--out", str(ledger)])[0]
        )
        code, stdout, err = self._main(
            [
                "--shard-dir",
                str(self.shards),
                "--out",
                str(ledger),
                "--resume-from",
                str(ledger),
                "--force",
            ]
        )
        self.assertEqual(0, code, err)
        self.assertEqual(self.total, json.loads(stdout)["total_records"])

    def test_missing_shard_dir_raises_then_exits_non_zero(self):
        out = self.dir / "out.tsv"
        with self.assertRaises(ValueError):
            self.module.run(
                self.module.build_arg_namespace(
                    ["--shard-dir", str(self.dir / "nope"), "--out", str(out)]
                )
            )
        code, _, err = self._main(["--shard-dir", str(self.dir / "nope"), "--out", str(out)])
        self.assertEqual(2, code)
        self.assertIn("nope", err)

    def test_missing_explicit_file_raises_then_exits_non_zero(self):
        out = self.dir / "out.tsv"
        missing = str(self.dir / "nope.faa")
        with self.assertRaises(ValueError):
            self.module.run(
                self.module.build_arg_namespace(["--file", missing, "--out", str(out)])
            )
        code, _, err = self._main(["--file", missing, "--out", str(out)])
        self.assertEqual(2, code)
        self.assertIn("nope.faa", err)

    def test_no_inputs_at_all_is_an_error(self):
        code, _, err = self._main(["--out", str(self.dir / "out.tsv")])
        self.assertEqual(2, code)
        self.assertTrue(err.strip())

    def test_empty_shard_directory_is_an_error(self):
        empty = self.dir / "empty_dir"
        empty.mkdir()
        out = self.dir / "out.tsv"
        code, _, err = self._main(["--shard-dir", str(empty), "--out", str(out)])
        self.assertEqual(2, code)
        self.assertIn("no FASTA", err)
        self.assertFalse(out.exists())

    def test_jobs_must_be_a_positive_integer(self):
        for bad in ("0", "-1", "two"):
            with self.subTest(jobs=bad):
                # A distinct --out per case: the first iteration must not make
                # the next one fail for the wrong reason.
                out = self.dir / f"out_jobs_{bad}.tsv"
                code, _, err = self._main(
                    [
                        "--shard-dir",
                        str(self.shards),
                        "--out",
                        str(out),
                        "--jobs",
                        bad,
                    ]
                )
                self.assertNotEqual(0, code)
                self.assertTrue(err.strip())
                self.assertFalse(out.exists())

    def test_jobs_are_used_in_parallel(self):
        """--jobs N really runs N workers; counts stay deterministic."""
        gate = threading.Barrier(4, timeout=30)
        real_count = self.module.count_fasta_file

        def slow_count(path, **kwargs):
            gate.wait()
            return real_count(path, **kwargs)

        out = self.dir / "parallel.tsv"
        with mock.patch.object(self.module, "count_fasta_file", slow_count):
            code, stdout, err = self._main(
                ["--shard-dir", str(self.shards), "--out", str(out), "--jobs", "4"]
            )
        self.assertEqual(0, code, err)
        self.assertEqual(self.total, json.loads(stdout)["total_records"])

    def test_output_order_is_deterministic_regardless_of_jobs(self):
        first = self.dir / "j1.tsv"
        second = self.dir / "j8.tsv"
        self.assertEqual(
            0,
            self._main(["--shard-dir", str(self.shards), "--out", str(first), "--jobs", "1"])[0],
        )
        self.assertEqual(
            0,
            self._main(["--shard-dir", str(self.shards), "--out", str(second), "--jobs", "8"])[0],
        )
        self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_gzip_output_is_readable_by_the_loader(self):
        import gzip

        out = self.dir / "counts.tsv.gz"
        code, stdout, err = self._main(
            ["--shard-dir", str(self.shards), "--out", str(out), "--gzip"]
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(self.total, payload["total_records"])
        with gzip.open(out, "rt", encoding="utf-8", newline="") as handle:
            self.assertIn("shard_0001", handle.read())
        rows = self.module.read_counts_table(out)
        self.assertEqual(self.total, sum(int(r["records"]) for r in rows))

    def test_summary_and_ledger_bind_the_same_total(self):
        out = self.dir / "counts.tsv"
        ledger = self.dir / "manifest.json"
        code, stdout, err = self._main(
            [
                "--shard-dir",
                str(self.shards),
                "--out",
                str(out),
                "--summary-json",
                str(ledger),
                "--command",
                "python pipeline/scripts/count_fasta_records.py --shard-dir shards "
                "--out counts.tsv",
            ]
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        recorded = json.loads(ledger.read_text(encoding="utf-8"))
        self.assertEqual(payload["total_records"], recorded["database_size_Z"])
        self.assertEqual(
            "sum of per-shard record counts (one count per frozen scan shard)",
            recorded["database_size_basis"],
        )
        self.assertIn("count_fasta_records.py", recorded["command"])
        self.assertEqual(4, len(recorded["shards_detail"]))
        self.assertEqual(self.total, sum(s["records"] for s in recorded["shards_detail"]))


class ContractAndCliSurfaceTests(_CountFixture):
    """Argument surface and module-level contract the coordinator relies on."""

    def test_module_exposes_the_documented_helpers(self):
        for name in (
            "count_records",
            "count_fasta_file",
            "merge_counts",
            "read_counts_table",
            "build_arg_namespace",
            "run",
            "main",
        ):
            with self.subTest(name=name):
                self.assertTrue(
                    callable(getattr(self.module, name)), f"missing helper {name}"
                )

    def test_documented_flags_exist(self):
        args = self.module.build_arg_namespace(
            [
                "--shard-dir",
                "d",
                "--file",
                "f.faa",
                "--out",
                "o.tsv",
                "--jobs",
                "8",
                "--strict",
                "--sha256",
                "--resume-from",
                "prev.tsv",
                "--force",
                "--gzip",
            ]
        )
        self.assertEqual(8, args.jobs)
        self.assertTrue(args.strict)
        self.assertTrue(args.sha256)
        self.assertTrue(args.force)
        self.assertTrue(args.gzip)
        self.assertEqual("prev.tsv", str(args.resume_from))

    def test_sha256_defaults_off_and_jobs_defaults_to_one(self):
        args = self.module.build_arg_namespace(["--shard-dir", "d", "--out", "o.tsv"])
        self.assertFalse(args.sha256)
        self.assertFalse(args.strict)
        self.assertFalse(args.force)
        self.assertEqual(1, args.jobs)

    def test_strict_flows_through_the_cli(self):
        bad = self.dir / "bad"
        bad.mkdir()
        (bad / "shard_0001.faa").write_bytes(b"MKAV\n>a\nMK\n")
        out = self.dir / "strict.tsv"
        code, _, err = self._main(["--shard-dir", str(bad), "--out", str(out)])
        self.assertEqual(0, code, err)
        with self.assertRaises(ValueError):
            self.module.run(
                self.module.build_arg_namespace(
                    ["--shard-dir", str(bad), "--out", str(self.dir / "s2.tsv"), "--strict"]
                )
            )
        code, _, err = self._main(
            ["--shard-dir", str(bad), "--out", str(self.dir / "s3.tsv"), "--strict"]
        )
        self.assertEqual(2, code)
        self.assertIn("before any FASTA header", err)

    def test_max_line_bytes_default_matches_the_documented_budget(self):
        self.assertEqual(DEFAULT_MAX_LINE_BYTES, self.module.DEFAULT_MAX_LINE_BYTES)

    def test_docstring_documents_the_grep_equivalence(self):
        doc = self.module.__doc__ or ""
        self.assertIn("grep -c", doc)
        self.assertIn(">", doc)


if __name__ == "__main__":
    unittest.main()
