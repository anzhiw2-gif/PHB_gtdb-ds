"""Task F7 — rescale the frozen scan-13 E-values onto one full-library ``Z``.

The frozen scan 13 ran ``hmmsearch`` without ``-Z``, so each of the 1,000 shard
searches used its own sequence count as ``Z`` and the 6.7 M reported E-values sit
on 100 different scales.  Linearity of ``E`` in ``Z`` has already been proven on
frozen evidence (``-Z 2920000 -> E=1.4e-08`` vs ``-Z 109087 -> E=5.1e-10``;
predicted ratio 26.80, measured 26.9), so the reconciliation needs **no HMMER
run** at all: ``E_full = E_shard * (Z_total / Z_shard)``.

These tests are the failing test for
``pipeline/scripts/rescale_phaded_scan_hits.py``.  They pin:

* the arithmetic, delegated to the project's existing
  ``parse_phaded_cys_targeted_recall.rescale_evalue`` (never re-derived here),
  and its agreement with the frozen measured ratio;
* monotonicity (a larger ``Z_shard`` means a smaller rescaled E-value) and the
  identity case ``Z_shard == Z_total``;
* exact reproduction of the two frozen measurements (26.80 predicted vs 26.9
  measured) and the working example ``1.4e-08 with 2920000/292000000 -> 1.4e-06``;
* that the frozen input is never the output, that a row from a shard absent from
  the counts file is an ERROR rather than a silent skip, and that
  ``--allow-unknown-shard`` downgrades it to a counted ``pending`` bucket;
* that the 696 MB hits table is **streamed**: whole lines, bounded reads,
  never loaded into memory;
* fail-closed behaviour when neither ``--database-size-z`` nor a usable
  shard-count total exists, and ``main(argv=None) -> int`` exit codes.

Nothing here touches ``runs/``, ``results/`` or ``deploy/``: every fixture is a
synthetic file inside a ``tempfile`` directory.
"""

from __future__ import annotations

import importlib.util
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "pipeline" / "scripts"
MODULE_PATH = SCRIPTS / "rescale_phaded_scan_hits.py"
RESCALE_SOURCE = SCRIPTS / "parse_phaded_cys_targeted_recall.py"

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

#: Frozen measurements (docs/T141_20260921_phaded_comprehensive_review.md and the
#: F7 plan): shard_0001 = 4,694,121 records, shard_0051 = 2,923,820 records.
Z_SHARD_0001 = 4_694_121
Z_SHARD_0051 = 2_923_820
#: The long-standing approximate full-library count, used only as a synthetic
#: stand-in here.  The real ``Z_total`` is measured by count_fasta_records.py.
Z_TOTAL_APPROX = 292_000_000
#: From the frozen check: 2923820 / 109087 = 26.803.
FROZEN_POOL_Z = 109_087
FROZEN_EXPECTED_RATIO = 26.80


def load_module():
    """Import the rescaling script by path.

    ``pipeline/scripts`` is not a package, so the module is loaded by file path
    the way the shell wrappers run it.  The reused convention module
    (``parse_phaded_cys_targeted_recall``) is loaded and registered first, so the
    script's own import resolves to the *project* function rather than to a copy.
    """
    if not MODULE_PATH.is_file():
        raise AssertionError(f"missing rescaling script: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "parse_phaded_cys_targeted_recall", RESCALE_SOURCE
    )
    convention = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("parse_phaded_cys_targeted_recall", convention)
    spec.loader.exec_module(convention)
    sys.modules.setdefault("rescale_convention_f7", convention)

    spec = importlib.util.spec_from_file_location("rescale_hits_f7", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("rescale_hits_f7", module)
    spec.loader.exec_module(module)
    return module


def hits_line(
    *,
    family: str = "DED_hfam_2",
    shard: str = "shard_0051",
    protein: str = "GCA_1|contig_2",
    tacc: str = "GCA_1_2",
    evalue: str = "1.4e-08",
    score: str = "45.6",
    bias: str = "0.1",
    domE: str = "1.4e-08",
    qname: str = "hfam_2",
    cov: str = "0.87",
) -> str:
    """Return one synthetic ``hits_all.tsv`` data row (tab separated)."""
    return "\t".join(
        [family, shard, protein, tacc, evalue, score, bias, domE, qname, cov]
    )


class RescaleRowTests(unittest.TestCase):
    """``rescale_row`` delegates to the project's ``rescale_evalue``."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_reuses_the_project_convention_and_not_a_copy(self):
        """The arithmetic must come from the existing project function."""
        self.assertIs(self.module.rescale_evalue, sys.modules[
            "parse_phaded_cys_targeted_recall"
        ].rescale_evalue)
        with mock.patch.object(
            self.module, "rescale_evalue", side_effect=AssertionError("called")
        ):
            with self.assertRaises(AssertionError):
                self.module.rescale_row(1e-5, 10, 100)

    def test_shared_upstream_function_still_has_the_expected_shape(self):
        upstream = sys.modules["parse_phaded_cys_targeted_recall"]
        self.assertEqual(2e-10, upstream.rescale_evalue(1e-10, 109087, 218174))
        with self.assertRaises(ValueError):
            upstream.rescale_evalue(1e-10, 0, 10)

    def test_working_example_from_the_plan(self):
        self.assertAlmostEqual(
            1.4e-06,
            self.module.rescale_row(1.4e-08, z_shard=2_920_000, z_total=292_000_000),
            places=18,
        )

    def test_frozen_measured_ratio_is_reproduced(self):
        """26.80 predicted from Z_full/Z_recall vs 26.9 measured on frozen data.

        Frozen evidence (same target): ``-Z 2920000 -> E=1.4e-08`` and
        ``-Z 109087 -> E=5.1e-10``.  E is proportional to Z, so the ratio between
        the two reported E-values is the ratio of the two Z values:
        2923820 / 109087 = 26.803 predicted, 26.9 measured.
        """
        predicted = self.module.rescale_row(5.1e-10, FROZEN_POOL_Z, Z_SHARD_0051) / 5.1e-10
        self.assertAlmostEqual(FROZEN_EXPECTED_RATIO, predicted, places=2)
        self.assertAlmostEqual(26.9, predicted, delta=0.15)

    def test_recall_pool_ratio_is_reproduced(self):
        """Rescaling the 109,087-protein pool E onto the shard Z gives ~1.4e-08.

        The frozen check recorded both numbers to two significant figures, so the
        reconstructed value must match within that rounding, not exactly.
        """
        ratio = Z_SHARD_0051 / FROZEN_POOL_Z
        self.assertAlmostEqual(FROZEN_EXPECTED_RATIO, ratio, places=2)
        reconstructed = self.module.rescale_row(5.1e-10, FROZEN_POOL_Z, Z_SHARD_0051)
        self.assertAlmostEqual(1.4e-08, reconstructed, delta=0.05 * 1.4e-08)
        self.assertLessEqual(abs(reconstructed - 1.4e-08) / 1.4e-08, 0.03)

    def test_plan_ratio_is_also_what_the_broken_scale_produced(self):
        """The inverse direction: rescaling shard -> full multiplies by 1/26.8."""
        self.assertAlmostEqual(
            1.0 / FROZEN_EXPECTED_RATIO,
            self.module.rescale_row(1.0, Z_SHARD_0051, FROZEN_POOL_Z),
            places=4,
        )

    def test_monotonic_in_z_shard(self):
        """Larger Z_shard means the shard under-reported, so E grows."""
        values = [
            self.module.rescale_row(1e-5, z, Z_TOTAL_APPROX)
            for z in (1_000, 100_000, 2_000_000, 4_694_121, 100_000_000)
        ]
        for smaller, larger in zip(values, values[1:]):
            self.assertLess(larger, smaller)

    def test_identity_when_scales_match(self):
        self.assertEqual(1e-5, self.module.rescale_row(1e-5, 7, 7))

    def test_larger_z_total_increases_e(self):
        small = self.module.rescale_row(1e-5, 1_000, 1_000_000)
        large = self.module.rescale_row(1e-5, 1_000, 100_000_000)
        self.assertGreater(large, small)

    def test_non_positive_or_unusable_scales_raise(self):
        for z_shard, z_total in ((0, 10), (10, 0), (-1, 10), (10, -1), (True, 10)):
            with self.subTest(z_shard=z_shard, z_total=z_total):
                with self.assertRaises(ValueError):
                    self.module.rescale_row(1e-5, z_shard, z_total)

    def test_evalue_must_be_finite_and_non_negative(self):
        for bad in (float("nan"), float("inf"), -1e-5, "abc", None):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    self.module.rescale_row(bad, 10, 100)


class ParseEvalueTests(unittest.TestCase):
    """The frozen ``E-value`` column must parse in every spelling it can take."""

    def setUp(self):
        self.module = load_module()

    def test_accepts_the_spellings_the_frozen_table_can_carry(self):
        for text, expected in (
            ("1.4e-08", 1.4e-08),
            ("1e-5", 1e-5),
            ("1.4E-08", 1.4e-08),
            ("0", 0.0),
            ("0.0", 0.0),
            ("3", 3.0),
            # below the smallest positive double (~4.94e-324): underflows to 0.0
            # rather than crashing.  HMMER never reports a value this small, but
            # a corrupted cell must not take the whole run down either.
            ("1e-400", 0.0),
        ):
            with self.subTest(text=text):
                self.assertEqual(expected, self.module.parse_evalue(text))

    def test_rejects_a_missing_or_garbage_value(self):
        for bad in ("", "  ", "abc", "1e", "nan", "inf", "-1e-5", None):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    self.module.parse_evalue(bad)

    def test_accepts_a_float(self):
        self.assertEqual(1.4e-08, self.module.parse_evalue(1.4e-08))


class _RescaleFixture(unittest.TestCase):
    """A synthetic hits table plus a matching shard-count ledger."""

    #: shard -> (records, rows written).  Z_shard < Z_total for every shard, so
    #: rescaling can only push E-values UP and therefore only reject rows: the
    #: fixture is built to make that direction observable.
    SHARDS = {
        "shard_0001": (4_694_121, 2),
        "shard_0051": (2_923_820, 3),
        "shard_0099": (1_000_000, 1),
    }
    Z_TOTAL = 8_617_941

    #: Frozen E-value per shard.  Every value passes the frozen ``-E 1e-5`` on
    #: that shard's own scale, and exactly one shard (``shard_0099``, the
    #: smallest count and therefore the largest inflation) is pushed back over
    #: the threshold by rescaling onto Z_total — so the fixture's delta is
    #: non-zero and its direction is observable.
    EVALUES = {
        "shard_0001": "4.6e-06",
        "shard_0051": "2.6e-09",
        "shard_0099": "9.2e-06",
    }

    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.counts = self.dir / "shard_counts.tsv"
        self.hits = self.dir / "hits_all.tsv"
        self.out = self.dir / "hits_rescaled.tsv"
        self.delta = self.dir / "delta.json"
        self._write_counts()
        self._write_hits()

    def _write_counts(
        self,
        rows=None,
        *,
        include_header: bool = True,
        columns=("name", "path", "bytes", "records", "status", "pending_reason"),
    ):
        rows = self.SHARDS if rows is None else rows
        lines = []
        if include_header:
            lines.append("\t".join(columns))
        for name, (records, _) in rows.items():
            cells = {
                "name": name,
                "path": f"SCAN_SHARDS_DIR/{name}.faa",
                "bytes": "1000",
                "records": str(records),
                "status": "counted",
                "pending_reason": "",
            }
            lines.append("\t".join(cells.get(c, "") for c in columns))
        self.counts.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return self.counts

    def _write_hits(self, rows=None):
        """Write a synthetic hits table; ``rows`` is a list of column dicts."""
        if rows is None:
            rows = []
            for shard, (_, count) in self.SHARDS.items():
                for index in range(count):
                    rows.append(
                        {
                            "family": "DED_hfam_2" if shard != "shard_0099" else "DED_hfam_9",
                            "shard": shard,
                            "protein": f"{shard}|p{index}",
                            "tacc": f"{shard}_t{index}",
                            "E-value": self.EVALUES[shard],
                            "score": "45.6",
                            "bias": "0.1",
                            "domE": "1.4e-08",
                            "qname": "hfam_2",
                            "cov": "0.87",
                        }
                    )
        lines = ["\t".join(HITS_COLUMNS)]
        for row in rows:
            lines.append("\t".join(str(row.get(c, "")) for c in HITS_COLUMNS))
        self.hits.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return self.hits

    def _main(self, argv):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = self.module.main(argv)
            except SystemExit as exit_error:
                code = exit_error.code if isinstance(exit_error.code, int) else 1
        return code, stdout.getvalue(), stderr.getvalue()

    def _read_out(self):
        lines = self.out.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        return header, [
            dict(zip(header, line.split("\t"))) for line in lines[1:] if line
        ]

    def _argv(self, *extra):
        return [
            "--hits",
            str(self.hits),
            "--shard-counts",
            str(self.counts),
            "--out",
            str(self.out),
            "--delta-json",
            str(self.delta),
            *extra,
        ]


class RescaleRunTests(_RescaleFixture):
    """End-to-end rescaling of a synthetic frozen table."""

    def test_rescales_every_row_and_preserves_the_original_value(self):
        code, stdout, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        header, rows = self._read_out()
        # every original column survives, plus the new scale columns
        for column in HITS_COLUMNS:
            self.assertIn(column, header)
        self.assertIn("E-value_shard_z", header)
        self.assertIn("E-value_rescaled", header)
        self.assertIn("z_shard", header)
        self.assertIn("z_total", header)
        self.assertIn("passes_rescaled_threshold", header)
        self.assertEqual(sum(count for _, count in self.SHARDS.values()), len(rows))

        by_shard = {}
        for row in rows:
            by_shard.setdefault(row["shard"], []).append(row)
        for shard, (records, count) in self.SHARDS.items():
            self.assertEqual(count, len(by_shard[shard]))
            for row in by_shard[shard]:
                self.assertEqual(str(records), row["z_shard"])
                self.assertEqual(str(self.Z_TOTAL), row["z_total"])
                # the original cell is preserved byte for byte
                self.assertEqual(row["E-value"], row["E-value_shard_z"])
                expected = self.module.rescale_row(
                    self.module.parse_evalue(row["E-value_shard_z"]), records,
                    self.Z_TOTAL,
                )
                self.assertAlmostEqual(
                    expected, float(row["E-value_rescaled"]), places=9
                )

    def test_threshold_is_applied_to_the_rescaled_value(self):
        code, _, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        _, rows = self._read_out()
        for row in rows:
            rescaled = float(row["E-value_rescaled"])
            expected = "true" if rescaled < 1e-5 else "false"
            self.assertEqual(expected, row["passes_rescaled_threshold"])
            if row["shard"] == "shard_0051":
                # 2.6e-09 * (8617941/2923820) = 7.7e-09: still far under 1e-5
                self.assertEqual("true", row["passes_rescaled_threshold"])
            if row["shard"] == "shard_0001":
                # 4.6e-06 * (8617941/4694121) = 8.4e-06: still under 1e-5
                self.assertEqual("true", row["passes_rescaled_threshold"])
            if row["shard"] == "shard_0099":
                # 9.2e-06 * (8617941/1000000) = 7.9e-05: the scale finally bites
                self.assertEqual("false", row["passes_rescaled_threshold"])

    def test_custom_threshold_changes_the_pass_column_only(self):
        code, _, err = self._main(self._argv("--threshold", "1e-3"))
        self.assertEqual(0, code, err)
        _, rows = self._read_out()
        for row in rows:
            rescaled = float(row["E-value_rescaled"])
            self.assertEqual(
                "true" if rescaled < 1e-3 else "false",
                row["passes_rescaled_threshold"],
            )

    def test_z_total_is_derived_from_the_counts_file_and_stated(self):
        code, stdout, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(self.Z_TOTAL, payload["z_total"])
        self.assertEqual("shard_counts_total", payload["z_total_basis"])
        delta = json.loads(self.delta.read_text(encoding="utf-8"))
        self.assertEqual(self.Z_TOTAL, delta["z_total"])
        self.assertEqual("shard_counts_total", delta["z_total_basis"])
        self.assertIn("sum", delta["z_total_basis_note"])

    def test_explicit_database_size_z_overrides_and_is_reported(self):
        code, stdout, err = self._main(self._argv("--database-size-z", "100000000"))
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(100_000_000, payload["z_total"])
        self.assertEqual("cli_database_size_z", payload["z_total_basis"])
        _, rows = self._read_out()
        self.assertEqual({"100000000"}, {row["z_total"] for row in rows})

    def test_database_size_z_must_be_a_positive_integer(self):
        for bad in ("0", "-5", "abc", "1.5"):
            with self.subTest(bad=bad):
                out = self.dir / f"out_{bad.replace('-', 'm').replace('.', '_')}.tsv"
                code, _, err = self._main(
                    [
                        "--hits",
                        str(self.hits),
                        "--shard-counts",
                        str(self.counts),
                        "--out",
                        str(out),
                        "--database-size-z",
                        bad,
                    ]
                )
                self.assertNotEqual(0, code)
                self.assertFalse(out.exists())
                self.assertTrue(err.strip())

    def test_fails_closed_when_no_z_can_be_established(self):
        empty = self.dir / "empty_counts.tsv"
        empty.write_text(
            "name\tpath\tbytes\trecords\tstatus\tpending_reason\n"
            "shard_0051\tp\t10\t\tnot_counted\tzero_records\n",
            encoding="utf-8",
        )
        out = self.dir / "never.tsv"
        code, _, err = self._main(
            ["--hits", str(self.hits), "--shard-counts", str(empty), "--out", str(out)]
        )
        self.assertEqual(5, code)
        self.assertIn("database-size-z", err)
        self.assertFalse(out.exists())
        # it is a typed failure, not a string match in main()
        with self.assertRaises(self.module.NoDatabaseSizeError):
            self.module.shard_counts_total({})

    def test_delta_report_arithmetic_is_consistent(self):
        code, stdout, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        delta = json.loads(self.delta.read_text(encoding="utf-8"))
        _, rows = self._read_out()

        self.assertEqual(len(rows), delta["total_rows"])
        self.assertEqual(len(rows), payload["total_rows"])
        passes_original = sum(
            1 for row in rows if self.module.parse_evalue(row["E-value_shard_z"]) < 1e-5
        )
        passes_rescaled = sum(
            1 for row in rows if row["passes_rescaled_threshold"] == "true"
        )
        self.assertEqual(passes_original, delta["rows_passing_original_threshold"])
        self.assertEqual(passes_rescaled, delta["rows_passing_rescaled_threshold"])
        self.assertEqual(
            passes_rescaled - passes_original, delta["delta_rows"]
        )
        self.assertEqual(
            delta["rows_passing_rescaled_threshold"]
            - delta["rows_passing_original_threshold"],
            delta["delta_rows"],
        )
        # direction and magnitude are stated, not implied
        self.assertIn(delta["direction"], {"newly_rejected", "newly_admitted", "unchanged"})
        if delta["delta_rows"] < 0:
            self.assertEqual("newly_rejected", delta["direction"])
            self.assertLess(delta["ratio_rescaled_over_original"], 1.0)
        elif delta["delta_rows"] > 0:
            self.assertEqual("newly_admitted", delta["direction"])
        else:
            self.assertEqual("unchanged", delta["direction"])
        self.assertIn("stricter", delta["statement"] + delta["statement_direction"])
        self.assertEqual("newly_rejected", delta["direction"])
        self.assertEqual(-1, delta["delta_rows"])
        self.assertLess(delta["ratio_rescaled_over_original"], 1.0)
        self.assertEqual("1e-05", delta["threshold_rescaled"])
        self.assertEqual("1e-05", delta["threshold_original"])
        self.assertEqual(
            delta["rows_passing_original_threshold"]
            + delta["rows_failing_original_threshold"],
            delta["rows_rescalable"],
        )

    def test_per_family_breakdown_sums_to_the_global_delta(self):
        code, _, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        delta = json.loads(self.delta.read_text(encoding="utf-8"))
        families = delta["per_family"]
        self.assertTrue(families)
        self.assertEqual(
            delta["rows_passing_original_threshold"],
            sum(entry["rows_passing_original_threshold"] for entry in families.values()),
        )
        self.assertEqual(
            delta["rows_passing_rescaled_threshold"],
            sum(entry["rows_passing_rescaled_threshold"] for entry in families.values()),
        )
        self.assertEqual(
            delta["total_rows"], sum(entry["total_rows"] for entry in families.values())
        )
        self.assertEqual(
            delta["delta_rows"], sum(entry["delta_rows"] for entry in families.values())
        )

    def test_shard_columns_and_provenance_are_preserved(self):
        code, _, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        _, rows = self._read_out()
        input_rows = [
            line.split("\t")
            for line in self.hits.read_text(encoding="utf-8").splitlines()[1:]
            if line
        ]
        self.assertEqual(len(input_rows), len(rows))
        for source, row in zip(input_rows, rows):
            for index, column in enumerate(HITS_COLUMNS):
                self.assertEqual(source[index], row[column])

    def test_shard_extension_spellings_resolve_to_the_same_key(self):
        """A counts row may name the file ``shard_0051.faa``; hits say ``shard_0051``."""
        self._write_counts(
            rows={
                "shard_0001.faa": (4_694_121, 2),
                "shard_0051.faa": (2_923_820, 3),
                "shard_0099.faa": (1_000_000, 1),
            }
        )
        code, _, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        _, rows = self._read_out()
        self.assertEqual(6, len(rows))
        self.assertEqual(
            {"4694121", "2923820", "1000000"}, {row["z_shard"] for row in rows}
        )

    def test_manifest_sidecar_records_the_scale_contract(self):
        manifest = self.dir / "rescale_manifest.json"
        code, _, err = self._main(
            self._argv(
                "--manifest",
                str(manifest),
                "--command",
                "python pipeline/scripts/rescale_phaded_scan_hits.py ...",
                "--input-record",
                "frozen_scan_run=SOME_RUN_DIR",
            )
        )
        self.assertEqual(0, code, err)
        recorded = json.loads(manifest.read_text(encoding="utf-8"))
        self.assertEqual(self.Z_TOTAL, recorded["database_size_Z"])
        self.assertEqual("shard_counts_total", recorded["database_size_basis"])
        self.assertEqual(len(self.SHARDS), len(recorded["shards"]))
        self.assertEqual(
            self.Z_TOTAL, sum(entry["records"] for entry in recorded["shards"])
        )
        self.assertIn("rescale_phaded_scan_hits.py", recorded["command"])
        self.assertIn("SOME_RUN_DIR", recorded["inputs"]["frozen_scan_run"])
        self.assertIn("hits_sha256", recorded)
        self.assertIn("statistics", recorded)
        # six rows; every one passes the frozen -E 1e-5 on its own shard's scale,
        # and shard_0099 alone is pushed over the threshold by rescaling
        self.assertEqual(6, recorded["statistics"]["rows_passing_original_threshold"])
        self.assertEqual(5, recorded["statistics"]["rows_passing_rescaled_threshold"])
        self.assertEqual(-1, recorded["statistics"]["delta_rows"])

    def test_hits_sha256_flag_binds_the_frozen_input(self):
        manifest = self.dir / "manifest.json"
        code, _, err = self._main(
            self._argv("--manifest", str(manifest), "--hits-sha256")
        )
        self.assertEqual(0, code, err)
        recorded = json.loads(manifest.read_text(encoding="utf-8"))
        self.assertIsNotNone(recorded["hits_sha256"])
        self.assertEqual(64, len(recorded["hits_sha256"]))
        # the frozen input must be byte-identical after the run
        self.assertEqual(
            recorded["hits_sha256"], self.module.sha256_file(self.hits)
        )


class RescaleUnknownShardTests(_RescaleFixture):
    """A shard missing from the counts file is an ERROR, never a silent skip."""

    def setUp(self):
        super().setUp()
        self._write_hits(
            rows=[
                {
                    "family": "DED_hfam_2",
                    "shard": "shard_0051",
                    "protein": "p1",
                    "tacc": "t1",
                    "E-value": "1.4e-08",
                    "score": "45.6",
                    "bias": "0.1",
                    "domE": "1.4e-08",
                    "qname": "hfam_2",
                    "cov": "0.87",
                },
                {
                    "family": "DED_hfam_2",
                    "shard": "shard_0777",
                    "protein": "p2",
                    "tacc": "t2",
                    "E-value": "1.4e-06",
                    "score": "40.0",
                    "bias": "0.2",
                    "domE": "1.4e-06",
                    "qname": "hfam_2",
                    "cov": "0.80",
                },
                {
                    "family": "DED_hfam_3",
                    "shard": "shard_0777",
                    "protein": "p3",
                    "tacc": "t3",
                    "E-value": "1.4e-07",
                    "score": "41.0",
                    "bias": "0.2",
                    "domE": "1.4e-07",
                    "qname": "hfam_3",
                    "cov": "0.81",
                },
            ]
        )

    def test_unknown_shard_raises_and_names_the_shard(self):
        out = self.dir / "out.tsv"
        argv = [
            "--hits",
            str(self.hits),
            "--shard-counts",
            str(self.counts),
            "--out",
            str(out),
        ]
        with self.assertRaises(ValueError) as caught:
            self.module.run(self.module.build_arg_namespace(argv))
        self.assertIn("shard_0777", str(caught.exception))
        code, _, err = self._main(argv)
        self.assertNotEqual(0, code)
        self.assertIn("shard_0777", err)
        self.assertFalse(out.exists())

    def test_allow_unknown_shard_counts_them_as_pending(self):
        code, stdout, err = self._main(self._argv("--allow-unknown-shard"))
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(2, payload["rows_pending_unknown_shard"])
        self.assertEqual(["shard_0777"], payload["unknown_shards"])
        header, rows = self._read_out()
        self.assertEqual(3, len(rows))
        pending = [row for row in rows if row["shard"] == "shard_0777"]
        self.assertEqual(2, len(pending))
        for row in pending:
            self.assertEqual("pending", row["E-value_rescaled"])
            self.assertEqual("pending", row["z_shard"])
            self.assertEqual(str(self.Z_TOTAL), row["z_total"])
            self.assertEqual("pending", row["passes_rescaled_threshold"])
            # the original value is still preserved for the dry run
            self.assertEqual(row["E-value"], row["E-value_shard_z"])

        delta = json.loads(self.delta.read_text(encoding="utf-8"))
        self.assertEqual(2, delta["rows_pending_unknown_shard"])
        self.assertEqual(["shard_0777"], delta["unknown_shards"])
        self.assertEqual(3, delta["total_rows"])
        self.assertEqual(
            2,
            delta["per_family"]["DED_hfam_2"]["rows_pending_unknown_shard"]
            + delta["per_family"]["DED_hfam_3"]["rows_pending_unknown_shard"],
        )
        # the pending rows must not be counted as passing or failing
        self.assertEqual(
            3,
            delta["rows_passing_original_threshold"]
            + delta["rows_failing_original_threshold"],
        )
        self.assertEqual(
            delta["rows_passing_rescaled_threshold"]
            + delta["rows_failing_rescaled_threshold"],
            delta["rows_rescalable"],
        )

    def test_error_names_every_missing_shard(self):
        self._write_hits(
            rows=[
                {
                    "family": "DED_hfam_2",
                    "shard": shard,
                    "protein": f"p{n}",
                    "tacc": f"t{n}",
                    "E-value": "1.4e-08",
                    "score": "45.6",
                    "bias": "0.1",
                    "domE": "1.4e-08",
                    "qname": "hfam_2",
                    "cov": "0.87",
                }
                for n, shard in enumerate(["shard_0777", "shard_0888", "shard_0777"])
            ]
        )
        code, _, err = self._main(
            [
                "--hits",
                str(self.hits),
                "--shard-counts",
                str(self.counts),
                "--out",
                str(self.dir / "out.tsv"),
            ]
        )
        self.assertNotEqual(0, code)
        self.assertIn("shard_0777", err)
        self.assertIn("shard_0888", err)


class RescaleWriteSafetyTests(_RescaleFixture):
    """The frozen table is read-only evidence and must never be the output."""

    def test_out_equal_to_hits_is_refused(self):
        argv = [
            "--hits",
            str(self.hits),
            "--shard-counts",
            str(self.counts),
            "--out",
            str(self.hits),
        ]
        with self.assertRaises(FileExistsError) as caught:
            self.module.run(self.module.build_arg_namespace(argv))
        self.assertIn("--hits", str(caught.exception))
        before = self.hits.read_bytes()
        code, _, err = self._main(argv)
        self.assertEqual(4, code)
        self.assertIn("--hits", err)
        self.assertEqual(before, self.hits.read_bytes())

    def test_out_equal_to_hits_is_refused_for_a_different_spelling(self):
        spelling = self.dir / "sub" / ".." / self.hits.name
        code, _, err = self._main(
            [
                "--hits",
                str(self.hits),
                "--shard-counts",
                str(self.counts),
                "--out",
                str(spelling),
            ]
        )
        self.assertEqual(4, code)
        self.assertIn("--hits", err)

    def test_delta_json_equal_to_hits_is_refused(self):
        code, _, err = self._main(
            [
                "--hits",
                str(self.hits),
                "--shard-counts",
                str(self.counts),
                "--out",
                str(self.out),
                "--delta-json",
                str(self.hits),
            ]
        )
        self.assertEqual(4, code)
        self.assertIn("--hits", err)

    def test_existing_output_is_refused_without_force(self):
        self.out.write_text("PRECIOUS\n", encoding="utf-8")
        code, _, err = self._main(self._argv())
        self.assertEqual(4, code)
        self.assertIn("--force", err)
        self.assertEqual("PRECIOUS\n", self.out.read_text(encoding="utf-8"))

    def test_force_overwrites(self):
        self.out.write_text("PRECIOUS\n", encoding="utf-8")
        code, _, err = self._main(self._argv("--force"))
        self.assertEqual(0, code, err)
        self.assertIn("E-value_rescaled", self.out.read_text(encoding="utf-8"))

    def test_missing_inputs_raise_value_error_then_exit_non_zero(self):
        out = self.dir / "out.tsv"
        for flag, name in (
            ("--hits", "absent_hits.tsv"),
            ("--shard-counts", "absent_counts.tsv"),
        ):
            with self.subTest(flag=flag):
                argv = [
                    "--hits",
                    str(self.hits),
                    "--shard-counts",
                    str(self.counts),
                    "--out",
                    str(out),
                ]
                argv[argv.index(flag) + 1] = str(self.dir / name)
                with self.assertRaises(ValueError) as caught:
                    self.module.run(self.module.build_arg_namespace(argv))
                self.assertIn(name, str(caught.exception))
                code, _, err = self._main(argv)
                self.assertNotEqual(0, code)
                self.assertIn(name, err)
                self.assertFalse(out.exists())

    def test_counts_file_missing_the_records_column_fails_closed(self):
        self._write_counts(columns=("name", "path", "bytes", "status"))
        code, _, err = self._main(self._argv())
        self.assertNotEqual(0, code)
        self.assertIn("records", err)
        self.assertFalse(self.out.exists())


class _RecordingFile:
    """Wrap a text handle and record every ``read``/``readline`` size."""

    def __init__(self, handle, sizes):
        self._handle = handle
        self._sizes = sizes

    def read(self, *args, **kwargs):
        self._sizes.append(args[0] if args else None)
        if not args:
            raise AssertionError("read() without a size would load the whole table")
        return self._handle.read(*args, **kwargs)

    def readline(self, *args, **kwargs):
        size = args[0] if args else -1
        self._sizes.append(size)
        return self._handle.readline(*args, **kwargs)

    def __iter__(self):
        raise AssertionError("iterating the handle streams an unbounded buffer")

    def __getattr__(self, name):
        return getattr(self._handle, name)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._handle.close()
        return False


class RescaleStreamingTests(_RescaleFixture):
    """The 696 MB hits table must stream; it must never be read whole."""

    def test_large_hits_table_streams_with_bounded_reads(self):
        """Bounded-memory construction, verified by observed read sizes.

        A multi-megabyte synthetic table is built row by row, then the real
        reader runs while every ``read``/``readline`` size is recorded.  Any
        attempt to slurp the file (``read()`` with no size, or ``for line in
        handle``) raises inside the wrapper.
        """
        block = hits_line(protein="G" * 400, shard="shard_0051", evalue="2.6e-09")
        with self.hits.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(HITS_COLUMNS) + "\n")
            for _ in range(6000):
                handle.write(block + "\n")
        self.assertGreater(self.hits.stat().st_size, 2 * 1024 * 1024)

        sizes: list = []
        real_open = Path.open
        hits_path = Path(self.hits)

        def patched_open(self, *args, **kwargs):
            handle = real_open(self, *args, **kwargs)
            if Path(self) == hits_path:
                return _RecordingFile(handle, sizes)
            return handle

        # The pre-flight and the main pass each open the frozen table; both must
        # use bounded reads.  A ``_RecordingFile`` raises if either tries to
        # slurp the file or to iterate it line by line.
        with mock.patch.object(Path, "open", patched_open):
            code, _, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        self.assertTrue(sizes, "the hits table was never read")
        rows = self._read_out()[1]
        self.assertEqual(6000, len(rows))
        self.assertEqual(
            6000, json.loads(self.delta.read_text(encoding="utf-8"))["total_rows"]
        )

    def test_reader_yields_one_row_at_a_time_and_rejects_a_bad_header(self):
        self.hits.write_text("not\tthe\tright\theader\n", encoding="utf-8")
        with self.assertRaises(ValueError) as caught:
            list(self.module.iter_hits_rows(self.hits))
        message = str(caught.exception)
        self.assertIn("E-value", message)

    def test_iter_hits_rows_is_a_generator(self):
        import inspect

        self.assertTrue(inspect.isgeneratorfunction(self.module.iter_hits_rows))


class RescaleCliSurfaceTests(_RescaleFixture):
    """Argument surface and the documented contract."""

    def test_documented_flags_exist(self):
        args = self.module.build_arg_namespace(
            [
                "--hits",
                "h.tsv",
                "--shard-counts",
                "c.tsv",
                "--out",
                "o.tsv",
                "--delta-json",
                "d.json",
                "--manifest",
                "m.json",
                "--database-size-z",
                "292000000",
                "--threshold",
                "1e-5",
                "--allow-unknown-shard",
                "--force",
                "--limit",
                "10",
            ]
        )
        self.assertEqual(292_000_000, args.database_size_z)
        self.assertEqual(1e-5, args.threshold)
        self.assertTrue(args.allow_unknown_shard)
        self.assertTrue(args.force)
        self.assertEqual(10, args.limit)

    def test_threshold_defaults_to_the_frozen_reporting_threshold(self):
        args = self.module.build_arg_namespace(
            ["--hits", "h.tsv", "--shard-counts", "c.tsv", "--out", "o.tsv"]
        )
        self.assertEqual(1e-5, args.threshold)
        self.assertIsNone(args.database_size_z)
        self.assertFalse(args.allow_unknown_shard)

    def test_limit_stops_early_and_is_reported(self):
        code, stdout, err = self._main(self._argv("--limit", "2"))
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(2, payload["total_rows"])
        self.assertTrue(payload["limited"])

    def test_module_exposes_the_documented_helpers(self):
        for name in (
            "rescale_row",
            "rescale_evalue",
            "parse_evalue",
            "iter_hits_rows",
            "build_arg_namespace",
            "run",
            "main",
        ):
            with self.subTest(name=name):
                self.assertTrue(
                    callable(getattr(self.module, name)), f"missing helper {name}"
                )

    def test_docstring_states_the_formula_and_the_no_hmmer_claim(self):
        doc = (self.module.__doc__ or "").lower()
        self.assertIn("e_full", doc)
        self.assertIn("z_total", doc)
        self.assertIn("no hmmer run at all", doc)

    def test_frozen_table_is_not_touched_by_a_successful_run(self):
        before = self.hits.read_bytes()
        code, _, err = self._main(self._argv())
        self.assertEqual(0, code, err)
        self.assertEqual(before, self.hits.read_bytes())


if __name__ == "__main__":
    unittest.main()
