"""Task F7 follow-up — does the corrected E-value scale keep the v1 universe?

Task F7 rescaled the frozen scan-13 E-values onto one measured full-library
``Z`` and found that **437,193 of the 6,742,621** previously passing hit rows no
longer pass ``E < 1e-5`` on the corrected scale.  The frozen v1 candidate
universe (109,087 proteins, ``protein_layers.tsv``) was constructed from those
hits at the *per-shard* scale, so the open question is whether any universe
accession is left with **zero** surviving rows — i.e. was supported only by rows
that the corrected scale rejects.

These are the failing tests for
``pipeline/scripts/quantify_phaded_scan_scale_impact.py``.  They pin:

* the **join key**.  The universe carries ``accession``; the rescaled hits table
  carries ``protein`` **and** ``tacc``.  ``accession == protein`` is *expected*
  but is never assumed: the run reports accessions found in ``protein``, in
  ``tacc`` and in neither, and a row whose ``protein`` and ``tacc`` cells name
  *different* universe accessions fails closed (the tested ambiguity path);
* the four per-accession outcomes: one passing row (kept), only rejected rows
  (lost, and the rejection is **consequential**), rejected + passing (kept, and
  the rejected row is **redundant**, not consequential), and no hit row at all
  (counted separately, never silently dropped);
* the row-level accounting identity
  ``rejected == consequential + redundant`` and
  ``frozen_support == rescaled_support + rejected``;
* **pre-flight fail-closed** on the scale column: a missing
  ``passes_rescaled_threshold`` or an unexpected cell value names the observed
  values and refuses to infer pass/fail from ``E-value_rescaled``;
* **bounded memory / streaming**: the 985 MB hits table is read one bounded
  line at a time, never slurped, and the retained state is O(universe), not
  O(rows) — proven by an observed peak-allocation measurement;
* ``main(argv=None) -> int``, ``--limit`` reporting, output-dir refusal, and an
  optional cross-check against the frozen F7 ``scale_delta_report.json``.

Nothing here touches ``runs/``, ``results/`` or ``deploy/``: every fixture is a
synthetic file inside a ``tempfile`` directory.
"""

from __future__ import annotations

import importlib.util
import io
import json
import sys
import tempfile
import tracemalloc
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "pipeline" / "scripts"
MODULE_PATH = SCRIPTS / "quantify_phaded_scan_scale_impact.py"

#: The frozen v1 universe columns (``protein_layers.tsv``).
UNIVERSE_COLUMNS = [
    "accession",
    "genome",
    "family",
    "source_layer",
    "signal_type",
    "assignment",
    "length",
    "layer",
]

#: The frozen rescaled-hits columns: the original ten plus F7's five.
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
    "E-value_shard_z",
    "E-value_rescaled",
    "z_shard",
    "z_total",
    "passes_rescaled_threshold",
]

#: Default frozen reporting threshold; the per-shard ``E-value`` column is on
#: the frozen scale, so its own pass state is recomputed with this threshold.
FROZEN_THRESHOLD = 1e-5

#: Frozen F7 facts, reused here only as *synthetic* fixture magnitudes.
F7_TOTAL_ROWS = 6_743_197
F7_PASSING_ORIGINAL = 6_742_621
F7_PASSING_RESCALED = 6_305_428
F7_DELTA_ROWS = -437_193


def load_module():
    """Import the quantification script by file path (scripts is not a package)."""
    if not MODULE_PATH.is_file():
        raise AssertionError(f"missing quantification script: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location("quantify_scale_impact_f7", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("quantify_scale_impact_f7", module)
    spec.loader.exec_module(module)
    return module


def universe_row(
    *,
    accession: str,
    genome: str = "GCA_000000001.1",
    family: str = "ePhaZ",
    source_layer: str = "ePhaZ_curated_core",
    signal_type: str = "OTHER",
    assignment: str = "ePhaZ_like",
    length: str = "382",
    layer: str = "ePhaZ_curated_core_nonsecreted",
) -> dict:
    """Return one synthetic ``protein_layers.tsv`` row."""
    return {
        "accession": accession,
        "genome": genome,
        "family": family,
        "source_layer": source_layer,
        "signal_type": signal_type,
        "assignment": assignment,
        "length": length,
        "layer": layer,
    }


def hits_row(
    *,
    protein: str,
    tacc: str = "-",
    family: str = "ePhaZ_curated_core",
    shard: str = "shard_0051",
    evalue: str = "1.4e-08",
    rescaled: str = "3.0e-07",
    passes: str = "true",
    score: str = "45.6",
    bias: str = "0.1",
    domE: str = "1.4e-08",
    qname: str = "ePhaZ_curated_core_aln",
    cov: str = "0.87",
) -> dict:
    """Return one synthetic ``hits_all_rescaled.tsv`` row.

    ``E-value`` (frozen per-shard scale) and ``E-value_rescaled`` are spelled
    independently on purpose: the frozen pass state is recovered from ``E-value``
    while the surviving state comes from ``passes_rescaled_threshold``.
    """
    return {
        "family": family,
        "shard": shard,
        "protein": protein,
        "tacc": tacc,
        "E-value": evalue,
        "score": score,
        "bias": bias,
        "domE": domE,
        "qname": qname,
        "cov": cov,
        "E-value_shard_z": evalue,
        "E-value_rescaled": rescaled,
        "z_shard": "2923820",
        "z_total": "615969589",
        "passes_rescaled_threshold": passes,
    }


class _RecordingFile:
    """Wrap a text handle and refuse every unbounded read.

    The project idiom (``test_rescale_phaded_scan_hits.py``): patching
    ``builtins.open`` does **not** intercept ``Path.open``, so the recorder wraps
    the handle returned by a patched ``Path.open``.  An unsized ``read()`` or a
    line iteration would buffer the whole 985 MB table, so both raise here.
    """

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
        if size is None or (isinstance(size, int) and size < 0):
            raise AssertionError("readline() without a bound streams an unbounded buffer")
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


class _ScaleImpactFixture(unittest.TestCase):
    """A synthetic universe plus a synthetic rescaled-hits table."""

    #: accession -> (genome, family, source_layer, layer)
    UNIVERSE = {
        "ACC_KEPT": ("G_1", "ePhaZ", "ePhaZ_curated_core", "ePhaZ_curated_core_secreted"),
        "ACC_LOST": ("G_2", "iPhaZ", "iPhaZ_curated_core", "iPhaZ_curated_core_secreted"),
        "ACC_MIXED": ("G_3", "PhaC", "PhaC_registry", "PhaC_nonsecreted"),
        "ACC_ABSENT": ("G_4", "BdhA", "BdhA_registry", "BdhA_nonsecreted"),
        "ACC_MULTIROW": ("G_5", "ePhaZ", "ePhaZ_broad_discovery", "ePhaZ_broad_nonsecreted"),
    }

    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.universe = self.dir / "protein_layers.tsv"
        self.hits = self.dir / "hits_all_rescaled.tsv"
        self.out_dir = self.dir / "out"
        self._write_universe()
        self._write_hits()

    # -- fixture writers -------------------------------------------------
    def _write_universe(self, rows=None):
        rows = self.UNIVERSE if rows is None else rows
        if isinstance(rows, dict):
            rows = [
                universe_row(
                    accession=accession,
                    genome=payload[0],
                    family=payload[1],
                    source_layer=payload[2],
                    layer=payload[3],
                )
                for accession, payload in rows.items()
            ]
        lines = ["\t".join(UNIVERSE_COLUMNS)]
        for row in rows:
            lines.append("\t".join(str(row.get(c, "")) for c in UNIVERSE_COLUMNS))
        self.universe.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return self.universe

    def _write_hits(self, rows=None):
        if rows is None:
            rows = [
                # one passing row -> keeps support
                hits_row(protein="ACC_KEPT"),
                # only rejected rows -> loses everything, consequential
                hits_row(
                    protein="ACC_LOST",
                    family="iPhaZ",
                    evalue="4.6e-06",
                    rescaled="9.7e-05",
                    passes="false",
                ),
                # rejected + passing -> keeps support, rejected row redundant
                hits_row(
                    protein="ACC_MIXED",
                    family="PhaC",
                    evalue="4.6e-06",
                    rescaled="9.7e-05",
                    passes="false",
                ),
                hits_row(
                    protein="ACC_MIXED",
                    family="PhaC",
                    evalue="2.6e-09",
                    rescaled="5.5e-08",
                    passes="true",
                ),
                # two passing rows -> keeps support twice over
                hits_row(protein="ACC_MULTIROW", family="ePhaZ_broad_discovery"),
                hits_row(protein="ACC_MULTIROW", family="ePhaZ_broad_discovery"),
                # a row that no universe accession claims (the table is wider
                # than the universe: 6.7 M rows vs 109,087 accessions)
                hits_row(protein="ACC_NOT_IN_UNIVERSE"),
            ]
        lines = ["\t".join(HITS_COLUMNS)]
        for row in rows:
            lines.append("\t".join(str(row.get(c, "")) for c in HITS_COLUMNS))
        self.hits.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return self.hits

    # -- invocation helpers ---------------------------------------------
    def _argv(self, *extra, out_dir=None):
        return [
            "--universe",
            str(self.universe),
            "--hits",
            str(self.hits),
            "--output-dir",
            str(out_dir if out_dir is not None else self.out_dir),
            "--expected-universe-size",
            str(len(self.UNIVERSE)),
            *extra,
        ]

    def _main(self, argv):
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = self.module.main(argv)
            except SystemExit as exit_error:
                code = exit_error.code if isinstance(exit_error.code, int) else 1
        return code, stdout.getvalue(), stderr.getvalue()

    def _summary(self):
        return json.loads(
            (self.out_dir / "scale_impact_summary.json").read_text(encoding="utf-8")
        )

    def _impact(self):
        path = self.out_dir / "accession_scale_impact.tsv"
        lines = path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        return header, [
            dict(zip(header, line.split("\t"))) for line in lines[1:] if line
        ]

    def _by_accession(self):
        _, rows = self._impact()
        return {row["accession"]: row for row in rows}


class PerAccessionOutcomeTests(_ScaleImpactFixture):
    """The four outcomes the task names, one test each."""

    def setUp(self):
        super().setUp()
        # NAME_3 has no row at all in this fixture; explain that one absence so
        # the join check does not mask the per-accession outcomes under test.
        code, stdout, err = self._main(self._argv("--expected-no-hit", "1"))
        self.assertEqual(0, code, err)

    def test_one_passing_row_keeps_support(self):
        row = self._by_accession()["ACC_KEPT"]
        self.assertEqual("1", row["total_hit_rows"])
        self.assertEqual("1", row["rescaled_support_rows"])
        self.assertEqual("1", row["frozen_support_rows"])
        self.assertEqual("false", row["lost_all_support"])
        self.assertEqual("false", row["rejection_consequential"])

    def test_only_rejected_rows_loses_all_support_consequentially(self):
        row = self._by_accession()["ACC_LOST"]
        self.assertEqual("1", row["total_hit_rows"])
        self.assertEqual("0", row["rescaled_support_rows"])
        self.assertEqual("1", row["frozen_support_rows"])
        self.assertEqual("true", row["lost_all_support"])
        self.assertEqual("true", row["rejection_consequential"])

    def test_rejected_plus_passing_keeps_support_and_the_rejection_is_redundant(self):
        row = self._by_accession()["ACC_MIXED"]
        self.assertEqual("2", row["total_hit_rows"])
        self.assertEqual("1", row["rescaled_support_rows"])
        self.assertEqual("2", row["frozen_support_rows"])
        self.assertEqual("false", row["lost_all_support"])
        self.assertEqual("false", row["rejection_consequential"])
        summary = self._summary()
        # the fixture has exactly two rejected rows: ACC_LOST's (consequential,
        # it keeps nothing) and ACC_MIXED's (redundant, ACC_MIXED keeps a pass)
        self.assertEqual(2, summary["rows_newly_rejected"])
        self.assertEqual(1, summary["rows_newly_rejected_consequential"])
        self.assertEqual(1, summary["rows_newly_rejected_redundant"])

    def test_multi_row_support_keeps_support(self):
        row = self._by_accession()["ACC_MULTIROW"]
        self.assertEqual("2", row["total_hit_rows"])
        self.assertEqual("2", row["rescaled_support_rows"])
        self.assertEqual("false", row["lost_all_support"])

    def test_accession_with_no_hit_row_is_counted_separately_not_dropped(self):
        row = self._by_accession()["ACC_ABSENT"]
        self.assertEqual("0", row["total_hit_rows"])
        self.assertEqual("0", row["rescaled_support_rows"])
        self.assertEqual("0", row["frozen_support_rows"])
        # nothing was lost: there was never any row on this table
        self.assertEqual("false", row["lost_all_support"])
        self.assertEqual("false", row["rejection_consequential"])
        summary = self._summary()
        self.assertEqual(1, summary["accessions_in_no_hit_row"])
        self.assertEqual(["ACC_ABSENT"], summary["no_hit_accessions"])
        self.assertTrue(summary["no_hit_accessions_within_expected"])

    def test_every_universe_accession_appears_in_the_per_accession_table(self):
        self.assertEqual(
            set(self.UNIVERSE), set(self._by_accession())
        )
        summary = self._summary()
        self.assertEqual(len(self.UNIVERSE), summary["universe_accessions"])
        self.assertEqual(len(self.UNIVERSE), summary["accessions_processed"])

    def test_headline_counts(self):
        summary = self._summary()
        self.assertEqual(1, summary["accessions_lost_all_support"])
        self.assertEqual(1, summary["accessions_lost_all_support_consequential"])
        self.assertEqual(3, summary["accessions_retaining_rescaled_support"])
        self.assertEqual(7, summary["rows_total"])
        self.assertEqual(6, summary["rows_joining_universe"])

class NoHitRowExitTests(_ScaleImpactFixture):
    """A join-completeness failure is loud, explainable and never silent."""

    def test_no_hit_accessions_exit_non_zero_by_default(self):
        code, stdout, err = self._main(self._argv())
        # the fixture has one accession in no row at all
        self.assertNotEqual(0, code)
        self.assertIn("ACC_ABSENT", err)
        self.assertIn("no hit row", err)
        # the evidence tables are still published: they *are* the evidence
        self.assertTrue(
            (self.out_dir / "accession_scale_impact.tsv").is_file()
        )
        summary = self._summary()
        self.assertEqual(1, summary["accessions_in_no_hit_row"])
        self.assertEqual(["ACC_ABSENT"], summary["no_hit_accessions"])
        self.assertFalse(summary["no_hit_accessions_within_expected"])

    def test_explained_no_hit_count_is_accepted(self):
        code, stdout, err = self._main(self._argv("--expected-no-hit", "1"))
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertTrue(payload["no_hit_accessions_within_expected"])
        self.assertEqual(1, payload["accessions_in_no_hit_row"])

    def test_wrong_expected_no_hit_still_fails_closed(self):
        code, stdout, err = self._main(self._argv("--expected-no-hit", "0"))
        self.assertNotEqual(0, code)
        self.assertIn("ACC_ABSENT", err)

    def test_allow_flag_downgrades_but_still_reports(self):
        code, stdout, err = self._main(self._argv("--allow-no-hit-accessions"))
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(1, payload["accessions_in_no_hit_row"])
        self.assertTrue(payload["no_hit_check_downgraded"])
        self.assertFalse(payload["no_hit_accessions_within_expected"])


class JoinKeyTests(_ScaleImpactFixture):
    """``accession == protein`` is proven, never assumed."""

    def test_join_key_is_reported_per_column_and_auto_resolves_to_protein(self):
        code, stdout, err = self._main(self._argv("--expected-no-hit", "1"))
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        join = payload["join"]
        self.assertEqual("protein", join["join_key"])
        self.assertEqual("auto", join["join_key_basis"])
        self.assertEqual(5, join["universe_accessions"])
        # ACC_ABSENT has no row at all, so protein resolves the other four
        self.assertEqual(4, join["found_in_protein"])
        self.assertEqual(0, join["found_in_tacc"])
        self.assertEqual(1, join["found_in_neither"])
        # the per-column detail is the evidence for the claim
        self.assertEqual(4, join["per_key"]["protein"]["found"])
        self.assertEqual(0, join["per_key"]["tacc"]["found"])
        self.assertEqual(6, join["per_key"]["protein"]["rows_matched"])
        self.assertEqual(0, join["per_key"]["tacc"]["rows_matched"])
        self.assertEqual("protein_only", join["columns_agree"])
        self.assertEqual(0, join["conflicting_rows"])

    def test_tacc_join_is_reported_when_it_matches(self):
        """A table that spells the key in ``tacc`` is *also* counted."""
        self._write_hits(
            rows=[hits_row(protein="some_other_protein", tacc="ACC_KEPT")]
        )
        code, stdout, err = self._main(
            self._argv("--allow-no-hit-accessions")
        )
        self.assertEqual(0, code, err)
        join = json.loads(stdout)["join"]
        self.assertEqual(0, join["found_in_protein"])
        self.assertEqual(1, join["found_in_tacc"])
        self.assertEqual(4, join["found_in_neither"])
        self.assertEqual("tacc_only", join["columns_agree"])
        self.assertEqual(1, join["per_key"]["tacc"]["rows_matched"])

    def test_a_row_whose_protein_and_tacc_name_different_accessions_fails_closed(self):
        self._write_hits(
            rows=[
                hits_row(protein="ACC_KEPT", tacc="ACC_LOST"),
                hits_row(protein="ACC_MIXED"),
                hits_row(protein="ACC_LOST"),
            ]
        )
        out = self.dir / "ambiguous_out"
        code, stdout, err = self._main(
            self._argv("--expected-no-hit", "1", out_dir=out)
        )
        self.assertEqual(3, code)
        self.assertIn("ambiguous", err.lower())
        self.assertIn("ACC_LOST", err + stdout)
        self.assertIn("ACC_KEPT", err + stdout)
        self.assertFalse(out.exists())

    def test_ambiguous_join_is_typed_not_a_string_match(self):
        self._write_hits(
            rows=[hits_row(protein="ACC_KEPT", tacc="ACC_LOST")]
        )
        args = self.module.build_arg_namespace(self._argv())
        with self.assertRaises(self.module.JoinKeyAmbiguityError) as caught:
            self.module.run(args)
        self.assertIn("ACC_KEPT", str(caught.exception))
        self.assertIn("ACC_LOST", str(caught.exception))

    def test_protein_column_must_exist(self):
        self.hits.write_text(
            "\t".join(c for c in HITS_COLUMNS if c != "protein") + "\n",
            encoding="utf-8",
        )
        code, _, err = self._main(self._argv())
        self.assertNotEqual(0, code)
        self.assertIn("protein", err)

    def test_auto_join_key_requires_both_key_columns(self):
        """Without ``tacc`` the run cannot *prove* the key; it says so."""
        self.hits.write_text(
            "\t".join(c for c in HITS_COLUMNS if c != "tacc")
            + "\n"
            + "\t".join(
                [
                    "ePhaZ",
                    "shard_0051",
                    "ACC_KEPT",
                    "1e-08",
                    "1",
                    "0",
                    "1e-08",
                    "q",
                    "0.9",
                    "1e-08",
                    "3e-07",
                    "2923820",
                    "615969589",
                    "true",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        out = self.dir / "no_tacc"
        code, _, err = self._main(self._argv(out_dir=out))
        self.assertNotEqual(0, code)
        self.assertIn("tacc", err)
        self.assertFalse(out.exists())

    def test_universe_must_carry_its_columns(self):
        self.universe.write_text("accession\tfamily\nX\ty\n", encoding="utf-8")
        code, _, err = self._main(self._argv())
        self.assertNotEqual(0, code)
        self.assertIn("genome", err)
        self.assertFalse((self.out_dir / "accession_scale_impact.tsv").exists())


class ScaleColumnPreflightTests(_ScaleImpactFixture):
    """Never infer pass/fail when the frozen column exists."""

    def test_missing_pass_column_fails_closed_and_names_the_fallback(self):
        self.hits.write_text(
            "\t".join(c for c in HITS_COLUMNS if c != "passes_rescaled_threshold")
            + "\n"
            + "ePhaZ\tshard_0051\tACC_KEPT\t-\t1e-9\t1\t0\t1e-9\tq\t1\n",
            encoding="utf-8",
        )
        out = self.dir / "nofallback"
        code, stdout, err = self._main(self._argv(out_dir=out))
        self.assertNotEqual(0, code)
        self.assertIn("passes_rescaled_threshold", err)
        self.assertFalse(out.exists())
        # the message must say the run refuses rather than inventing the column
        self.assertIn("without", err.lower())

    def test_explicit_fallback_is_loudly_recorded(self):
        self.hits.write_text(
            "\t".join(c for c in HITS_COLUMNS if c != "passes_rescaled_threshold")
            + "\n"
            + "\t".join(
                [
                    "ePhaZ",
                    "shard_0051",
                    "ACC_KEPT",
                    "-",
                    "1e-09",
                    "1",
                    "0",
                    "1e-09",
                    "q",
                    "0.9",
                    "1e-09",
                    "5e-07",
                    "2923820",
                    "615969589",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        out = self.dir / "fallback"
        code, stdout, err = self._main(
            self._argv(
                "--infer-pass-from-rescaled",
                "--threshold",
                "1e-5",
                "--allow-no-hit-accessions",
                out_dir=out,
            )
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertTrue(payload["scale_column"]["inferred_from_rescaled"])
        self.assertEqual(
            "E-value_rescaled < 1e-05", payload["scale_column"]["inference_rule"]
        )

    def test_unexpected_pass_values_fail_closed_naming_them(self):
        self._write_hits(
            rows=[hits_row(protein="ACC_KEPT", passes="pending")]
        )
        code, stdout, err = self._main(self._argv())
        self.assertNotEqual(0, code)
        self.assertIn("pending", err)
        self.assertFalse((self.out_dir / "accession_scale_impact.tsv").exists())

    def test_other_false_spellings_are_named_too(self):
        self._write_hits(rows=[hits_row(protein="ACC_KEPT", passes="0")])
        code, _, err = self._main(self._argv())
        self.assertNotEqual(0, code)
        self.assertIn("0", err)
        self.assertIn("passes_rescaled_threshold", err)

    def test_accepted_spellings_are_reported_verbatim(self):
        rows = [
            hits_row(protein="ACC_KEPT", passes="true"),
            hits_row(protein="ACC_MULTIROW", passes="TRUE"),
        ]
        self._write_hits(rows=rows)
        code, stdout, err = self._main(self._argv("--allow-no-hit-accessions"))
        self.assertEqual(0, code, err)
        scale = json.loads(stdout)["scale_column"]
        self.assertEqual(["TRUE", "true"], scale["observed_values"])
        self.assertEqual("true", scale["true_value"])
        self.assertEqual("false", scale["false_value"])
        self.assertFalse(scale["inferred_from_rescaled"])

    def test_preflight_samples_before_any_output_is_written(self):
        """A bad value deep in the table is caught before output exists."""
        rows = [hits_row(protein="ACC_KEPT") for _ in range(50)]
        rows.append(hits_row(protein="ACC_KEPT", passes="maybe"))
        self._write_hits(rows=rows)
        out = self.dir / "late_bad"
        code, _, err = self._main(
            self._argv("--preflight-rows", "100", out_dir=out)
        )
        self.assertNotEqual(0, code)
        self.assertIn("maybe", err)
        self.assertFalse(out.exists())

    def test_inference_mode_needs_no_pass_column_at_all(self):
        """``--infer-pass-from-rescaled`` must not demand the absent column.

        Regression: the pre-flight used to read ``--pass-column`` even when the
        caller had explicitly asked for the rule-based fallback, so a table
        without the column failed with a misleading "missing column" error
        instead of being rescaled from ``E-value_rescaled``.
        """
        self.hits.write_text(
            "\t".join(c for c in HITS_COLUMNS if c != "passes_rescaled_threshold")
            + "\n"
            + "\t".join(
                [
                    "ePhaZ",
                    "shard_0051",
                    "ACC_KEPT",
                    "-",
                    "1e-09",
                    "1",
                    "0",
                    "1e-09",
                    "q",
                    "0.9",
                    "1e-09",
                    "5e-07",
                    "2923820",
                    "615969589",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        code, stdout, err = self._main(
            self._argv(
                "--infer-pass-from-rescaled",
                "--allow-no-hit-accessions",
            )
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertTrue(payload["scale_column"]["inferred_from_rescaled"])
        # 5e-07 < 1e-05, so the rule keeps this row
        self.assertEqual(1, payload["rows_passing_rescaled"])
        self.assertEqual(1, payload["rows_joining_universe"])
        rows = self._by_accession()
        self.assertEqual("1", rows["ACC_KEPT"]["rescaled_support_rows"])

    def test_inference_mode_still_refuses_a_bad_rescaled_cell(self):
        self.hits.write_text(
            "\t".join(c for c in HITS_COLUMNS if c != "passes_rescaled_threshold")
            + "\n"
            + "\t".join(
                [
                    "ePhaZ",
                    "shard_0051",
                    "ACC_KEPT",
                    "-",
                    "1e-09",
                    "1",
                    "0",
                    "1e-09",
                    "q",
                    "0.9",
                    "1e-09",
                    "not_a_number",
                    "2923820",
                    "615969589",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        out = self.dir / "bad_rescaled"
        code, _, err = self._main(
            self._argv("--infer-pass-from-rescaled", out_dir=out)
        )
        self.assertNotEqual(0, code)
        self.assertIn("not_a_number", err)
        self.assertFalse(out.exists())


class StreamingTests(_ScaleImpactFixture):
    """The 985 MB hits table must stream and the retained state is O(universe)."""

    #: Rows in the large synthetic table.  Each row is >400 bytes so the file is
    #: several MB, like the frozen table is hundreds of MB.
    BIG_ROWS = 6000

    #: Every big row names this universe accession (so the join is exercised);
    #: the padding lives in the *family* column, which is not a join key.
    BIG_ACCESSION = "ACC_KEPT"

    def _big_hits_row(self, index: int, filler: str) -> dict:
        return hits_row(
            protein=self.BIG_ACCESSION,
            family=f"ePhaZ_curated_core|padding_{filler}_{index}",
            tacc="-",
        )

    def _write_big_hits(self, rows: int | None = None) -> Path:
        count = self.BIG_ROWS if rows is None else rows
        filler = "G" * 400
        with self.hits.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(HITS_COLUMNS) + "\n")
            for index in range(count):
                row = self._big_hits_row(index, filler)
                handle.write("\t".join(str(row[c]) for c in HITS_COLUMNS) + "\n")
        return self.hits

    def test_large_hits_table_streams_with_bounded_reads(self):
        self._write_big_hits()
        self.assertGreater(self.hits.stat().st_size, 2 * 1024 * 1024)
        sizes: list = []
        real_open = Path.open
        hits_path = Path(self.hits)

        def patched_open(self, *args, **kwargs):
            handle = real_open(self, *args, **kwargs)
            if Path(self) == hits_path:
                return _RecordingFile(handle, sizes)
            return handle

        with mock.patch.object(Path, "open", patched_open):
            code, stdout, err = self._main(
                self._argv("--allow-no-hit-accessions")
            )
        self.assertEqual(0, code, err)
        self.assertTrue(sizes, "the hits table was never read")
        payload = json.loads(stdout)
        self.assertEqual(self.BIG_ROWS, payload["rows_total"])
        self.assertEqual(self.BIG_ROWS, payload["rows_joining_universe"])
        self.assertEqual(0, payload["rows_newly_rejected"])
        self.assertLessEqual(self.BIG_ROWS + 100, len(sizes))
        # the five-accession universe is the only retained state
        self.assertEqual(5, payload["bounded_state"]["universe_size"])

    def test_retained_state_does_not_grow_with_the_row_count(self):
        """Same universe, 4x the rows, same-order memory: O(universe) not O(rows)."""

        def peak(rows: int, out_name: str):
            self._write_big_hits(rows=rows)
            out_dir = self.dir / out_name
            tracemalloc.start()
            try:
                tracemalloc.reset_peak()
                code, stdout, err = self._main(
                    self._argv(
                        "--allow-no-hit-accessions",
                        out_dir=out_dir,
                    )
                )
                self.assertEqual(0, code, err)
                _, high = tracemalloc.get_traced_memory()
            finally:
                tracemalloc.stop()
            return high, json.loads(stdout)

        small_peak, small_payload = peak(1500, "mem_small")
        big_peak, big_payload = peak(6000, "mem_big")
        self.assertEqual(1500, small_payload["rows_total"])
        self.assertEqual(6000, big_payload["rows_total"])
        # 4x the rows for the same universe: the delta stays far below the
        # 4x-sized table itself (which is >8 MB of line text at 400 B/row)
        self.assertLess(abs(big_peak - small_peak), 8 * 1024 * 1024)
        # and the module states its own bound rather than leaving it implicit:
        # a fixed interpreter allowance plus a per-accession slot estimate
        self.assertEqual(
            self.module.FIXED_STATE_BYTES
            + self.module.BYTES_PER_UNIVERSE_ACCESSION * len(self.UNIVERSE),
            small_payload["bounded_state"]["estimated_bytes"],
        )
        self.assertEqual(
            small_payload["bounded_state"]["estimated_bytes"],
            self.module.memory_bound_bytes(len(self.UNIVERSE)),
        )
        # the model's per-accession term is a real byte-scale figure ...
        self.assertLess(self.module.BYTES_PER_UNIVERSE_ACCESSION, 4096)
        # ... and the measured peak for this 5-accession universe stays inside
        # the stated ceiling (the ceiling is a bound, not a prediction)
        self.assertLess(small_peak, self.module.memory_bound_bytes(len(self.UNIVERSE)))
        self.assertLess(big_peak, self.module.memory_bound_bytes(len(self.UNIVERSE)))
        # the ceiling never decreases as the universe grows, and the reported
        # figures are documented as a ceiling rather than an estimate
        self.assertLess(
            self.module.memory_bound_bytes(1000),
            self.module.memory_bound_bytes(1001),
        )
        model = small_payload["bounded_state"]["model"].lower()
        self.assertIn("ceiling", model)
        self.assertIn("measured_note", small_payload["bounded_state"])

    def test_iter_tsv_rows_is_a_generator_with_bounded_line_reads(self):
        self.assertTrue(
            __import__("inspect").isgeneratorfunction(self.module.iter_tsv_rows)
        )
        sizes: list = []
        real_open = Path.open
        hits_path = Path(self.hits)

        def patched_open(self, *args, **kwargs):
            handle = real_open(self, *args, **kwargs)
            if Path(self) == hits_path:
                return _RecordingFile(handle, sizes)
            return handle

        with mock.patch.object(Path, "open", patched_open):
            rows = list(
                self.module.iter_tsv_rows(self.hits, required_columns=HITS_COLUMNS)
            )
        # six rows name universe accessions, one names nobody (no header)
        self.assertEqual(6, len([r for r in rows if r["protein"] in self.UNIVERSE]))
        self.assertEqual(7, len(rows))
        self.assertEqual(
            {"ACC_NOT_IN_UNIVERSE"},
            {row["protein"] for row in rows if row["protein"] not in self.UNIVERSE},
        )
        self.assertTrue(all(size > 0 for size in sizes))

    def test_over_long_line_is_refused_rather_than_buffered(self):
        self.hits.write_text(
            "\t".join(HITS_COLUMNS) + "\n" + ("x" * 10) + "\n", encoding="utf-8"
        )
        with self.assertRaises(ValueError) as caught:
            list(
                self.module.iter_tsv_rows(
                    self.hits, required_columns=[], max_line_length=5
                )
            )
        self.assertIn("more than 5 bytes", str(caught.exception))


class CliAndSafetyTests(_ScaleImpactFixture):
    """Argument surface, exit codes and write safety."""

    def test_documented_flags_exist(self):
        args = self.module.build_arg_namespace(
            self._argv(
                "--limit",
                "10",
                "--threshold",
                "1e-5",
                "--expected-universe-size",
                "5",
                "--expected-no-hit",
                "0",
                "--join-key",
                "protein",
                "--pass-column",
                "passes_rescaled_threshold",
                "--preflight-rows",
                "1000",
                "--sha256",
                "--allow-no-hit-accessions",
            )
        )
        self.assertEqual(10, args.limit)
        self.assertEqual(1e-5, args.threshold)
        self.assertEqual(5, args.expected_universe_size)
        self.assertEqual(0, args.expected_no_hit)
        self.assertEqual("protein", args.join_key)
        self.assertEqual("passes_rescaled_threshold", args.pass_column)
        self.assertEqual(1000, args.preflight_rows)
        self.assertTrue(args.sha256)
        self.assertTrue(args.allow_no_hit_accessions)

    def test_defaults(self):
        args = self.module.build_arg_namespace(self._argv())
        # the fixture pins its own size; the module default is the frozen v1 one
        self.assertEqual(len(self.UNIVERSE), args.expected_universe_size)
        self.assertEqual(109_087, self.module.DEFAULT_EXPECTED_UNIVERSE_SIZE)
        default_args = self.module.build_arg_namespace(
            ["--universe", "u.tsv", "--hits", "h.tsv", "--output-dir", "o"]
        )
        self.assertEqual(109_087, default_args.expected_universe_size)
        self.assertEqual(1e-5, default_args.threshold)
        self.assertEqual("passes_rescaled_threshold", default_args.pass_column)
        self.assertEqual("auto", default_args.join_key)
        self.assertIsNone(default_args.limit)
        self.assertIsNone(default_args.expected_no_hit)
        self.assertFalse(default_args.allow_no_hit_accessions)
        self.assertFalse(default_args.sha256)
        self.assertFalse(default_args.infer_pass_from_rescaled)
        self.assertEqual(100_000, default_args.preflight_rows)

    def test_module_exposes_the_documented_helpers(self):
        for name in (
            "iter_tsv_rows",
            "read_universe",
            "pass_value",
            "memory_bound_bytes",
            "build_arg_namespace",
            "run",
            "main",
        ):
            with self.subTest(name=name):
                self.assertTrue(
                    callable(getattr(self.module, name)), f"missing helper {name}"
                )

    def test_limit_is_reported_and_said_to_be_incomplete(self):
        code, stdout, err = self._main(
            self._argv("--limit", "3", "--allow-no-hit-accessions")
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertTrue(payload["limited"])
        self.assertEqual(3, payload["limit"])
        self.assertEqual(3, payload["rows_total"])
        self.assertEqual(3, payload["rows_limit_scope"])
        # a bounded read cannot prove join completeness
        self.assertFalse(payload["join_completeness_provable"])
        self.assertIn("limit", (payload["limited_note"] or "").lower())

    def test_limit_zero_or_negative_is_refused(self):
        for bad in ("0", "-3"):
            with self.subTest(bad=bad):
                code, _, err = self._main(
                    self._argv("--limit", bad, "--allow-no-hit-accessions")
                )
                self.assertNotEqual(0, code)
                self.assertTrue(err.strip())

    def test_missing_inputs_raise_then_exit_non_zero(self):
        for flag, name in (
            ("--universe", "absent_universe.tsv"),
            ("--hits", "absent_hits.tsv"),
        ):
            with self.subTest(flag=flag):
                argv = self._argv("--allow-no-hit-accessions")
                argv[argv.index(flag) + 1] = str(self.dir / name)
                with self.assertRaises(ValueError) as caught:
                    self.module.run(self.module.build_arg_namespace(argv))
                self.assertIn(name, str(caught.exception))
                code, _, err = self._main(argv)
                self.assertNotEqual(0, code)
                self.assertIn(name, err)

    def test_non_empty_output_dir_is_refused(self):
        self.out_dir.mkdir()
        (self.out_dir / "PRECIOUS.txt").write_text("keep\n", encoding="utf-8")
        code, _, err = self._main(self._argv())
        self.assertEqual(4, code)
        self.assertIn("--force", err)
        self.assertEqual(
            "keep\n", (self.out_dir / "PRECIOUS.txt").read_text(encoding="utf-8")
        )
        self.assertFalse((self.out_dir / "accession_scale_impact.tsv").exists())

    def test_force_refuses_to_delete_unrelated_files(self):
        self.out_dir.mkdir()
        (self.out_dir / "PRECIOUS.txt").write_text("keep\n", encoding="utf-8")
        code, _, err = self._main(
            self._argv("--force", "--expected-no-hit", "1")
        )
        self.assertEqual(4, code)
        self.assertIn("PRECIOUS.txt", err)
        self.assertEqual(
            "keep\n", (self.out_dir / "PRECIOUS.txt").read_text(encoding="utf-8")
        )

    def test_force_overwrites_only_this_run_outputs(self):
        self.out_dir.mkdir()
        for name in (
            "accession_scale_impact.tsv",
            "scale_impact_summary.json",
            "accession_scale_impact_affected.tsv",
        ):
            (self.out_dir / name).write_text("stale\n", encoding="utf-8")
        code, stdout, err = self._main(
            self._argv("--force", "--expected-no-hit", "1")
        )
        self.assertEqual(0, code, err)
        self.assertEqual(1, json.loads(stdout)["accessions_in_no_hit_row"])
        self.assertIn(
            "accession",
            (self.out_dir / "accession_scale_impact.tsv")
            .read_text(encoding="utf-8")
            .splitlines()[0],
        )

    def test_universe_path_cannot_be_an_output_target(self):
        argv = self._argv("--expected-no-hit", "1")
        argv[argv.index("--output-dir") + 1] = str(self.dir)
        code, _, err = self._main(argv)
        self.assertEqual(4, code)
        self.assertIn("output-dir", err)

    def test_output_dir_must_not_be_a_file(self):
        target = self.dir / "not_a_dir"
        target.write_text("x\n", encoding="utf-8")
        code, _, err = self._main(self._argv(out_dir=target))
        self.assertNotEqual(0, code)
        self.assertIn("not_a_dir", err)

    def test_affected_list_contains_exactly_the_lost_accessions(self):
        code, _, err = self._main(self._argv("--expected-no-hit", "1"))
        self.assertEqual(0, code, err)
        path = self.out_dir / "accession_scale_impact_affected.tsv"
        lines = path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        rows = [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]
        self.assertEqual(["ACC_LOST"], [row["accession"] for row in rows])
        self.assertEqual("iPhaZ", rows[0]["family"])
        self.assertEqual("1", rows[0]["frozen_support_rows"])

    def test_universe_size_check_fails_closed(self):
        code, stdout, err = self._main(
            self._argv("--expected-universe-size", "7")
        )
        self.assertNotEqual(0, code)
        self.assertIn("7", err)
        self.assertIn(str(len(self.UNIVERSE)), err)
        with self.assertRaises(ValueError):
            self.module.run(
                self.module.build_arg_namespace(
                    self._argv("--expected-universe-size", "7")
                )
            )

    def test_universe_size_check_can_be_disabled_with_zero(self):
        code, stdout, err = self._main(
            self._argv(
                "--expected-universe-size",
                "0",
                "--expected-no-hit",
                "1",
            )
        )
        self.assertEqual(0, code, err)
        self.assertEqual(0, json.loads(stdout)["expected_universe_size"])

    def test_duplicate_accession_rows_are_refused(self):
        rows = [
            universe_row(accession="ACC_KEPT"),
            universe_row(accession="ACC_KEPT"),
        ]
        self._write_universe(rows=rows)
        code, _, err = self._main(self._argv())
        self.assertNotEqual(0, code)
        self.assertIn("ACC_KEPT", err)
        self.assertIn("duplicate", err.lower())


class CrossTabTests(_ScaleImpactFixture):
    """The decomposition the task asks for, by family and by layer."""

    def setUp(self):
        super().setUp()
        code, _, err = self._main(self._argv("--expected-no-hit", "1"))
        self.assertEqual(0, code, err)

    def test_per_family_counts_cover_every_universe_accession(self):
        families = self._summary()["cross_tabs"]["family"]
        self.assertEqual(
            len(self.UNIVERSE),
            sum(entry["universe_accessions"] for entry in families.values()),
        )
        self.assertEqual(
            1, families["iPhaZ"]["accessions_lost_all_support"]
        )
        self.assertEqual(
            1, families["iPhaZ"]["accessions_lost_all_support_consequential"]
        )
        self.assertEqual(
            0, families["ePhaZ"]["accessions_lost_all_support"]
        )
        self.assertEqual(
            1, families["BdhA"]["accessions_in_no_hit_row"]
        )

    def test_per_layer_counts_cover_every_universe_accession(self):
        layers = self._summary()["cross_tabs"]["layer"]
        self.assertEqual(
            len(self.UNIVERSE),
            sum(entry["universe_accessions"] for entry in layers.values()),
        )
        self.assertEqual(
            1,
            layers["iPhaZ_curated_core_secreted"]["accessions_lost_all_support"],
        )

    def test_source_layer_and_layer_are_both_reported(self):
        tabs = self._summary()["cross_tabs"]
        self.assertIn("source_layer", tabs)
        self.assertEqual(
            len(self.UNIVERSE),
            sum(
                entry["universe_accessions"]
                for entry in tabs["source_layer"].values()
            ),
        )
        combined = tabs["layer_x_family"]
        self.assertEqual(
            len(self.UNIVERSE),
            sum(entry["universe_accessions"] for entry in combined.values()),
        )
        self.assertIn(
            "ePhaZ_curated_core_secreted",
            {key.split("|")[0] for key in combined},
        )

    def test_row_accounting_identity_holds(self):
        summary = self._summary()
        self.assertEqual(
            summary["rows_newly_rejected"],
            summary["rows_newly_rejected_consequential"]
            + summary["rows_newly_rejected_redundant"],
        )
        self.assertEqual(
            summary["rows_passing_frozen"] - summary["rows_passing_rescaled"],
            summary["rows_newly_rejected"],
        )
        self.assertEqual(
            summary["rows_passing_frozen"],
            summary["rows_newly_rejected"] + summary["rows_passing_rescaled"],
        )

    def test_rejected_rows_are_counted_once_each(self):
        summary = self._summary()
        # ACC_LOST's single rejected row is consequential; ACC_MIXED's is not
        self.assertEqual(2, summary["rows_newly_rejected"])
        self.assertEqual(1, summary["rows_newly_rejected_consequential"])
        self.assertEqual(1, summary["rows_newly_rejected_redundant"])
        by_family = summary["cross_tabs"]["family"]
        self.assertEqual(1, by_family["iPhaZ"]["rows_newly_rejected_consequential"])
        self.assertEqual(1, by_family["PhaC"]["rows_newly_rejected_redundant"])
        self.assertEqual(0, by_family["PhaC"]["rows_newly_rejected_consequential"])
        self.assertEqual(
            2,
            by_family["iPhaZ"]["rows_newly_rejected"]
            + by_family["PhaC"]["rows_newly_rejected"],
        )

    def test_rows_claiming_no_accession_are_excluded_from_every_row_counter(self):
        """Regression: row counters must not advance for unclaimed rows.

        The frozen table is far wider than the universe (6,743,197 rows vs
        109,087 accessions).  If an unclaimed row moved ``rows_passing_frozen``
        or ``rows_newly_rejected``, the summary would no longer be the
        per-accession picture of the universe and the row-level/per-accession
        identities would silently disagree.
        """
        universe_rows = [
            universe_row(accession="ACC_ONE"),
            universe_row(accession="ACC_TWO"),
        ]
        self._write_universe(rows=universe_rows)
        self._write_hits(
            rows=[
                hits_row(protein="ACC_ONE"),
                hits_row(
                    protein="NOT_IN_UNIVERSE_A",
                    evalue="1e-09",
                    rescaled="1e-08",
                    passes="false",
                ),
                hits_row(
                    protein="NOT_IN_UNIVERSE_B",
                    evalue="1e-09",
                    rescaled="1e-08",
                    passes="true",
                ),
            ]
        )
        code, stdout, err = self._main(
            self._argv(
                "--force",
                "--allow-no-hit-accessions",
                "--expected-universe-size",
                "2",
            )
        )
        self.assertEqual(0, code, err)
        payload = json.loads(stdout)
        self.assertEqual(3, payload["rows_total"])
        self.assertEqual(1, payload["rows_joining_universe"])
        self.assertEqual(2, payload["rows_not_in_universe"])
        # nothing about the two unclaimed rows reaches any row-level counter
        self.assertEqual(1, payload["rows_passing_frozen"])
        self.assertEqual(0, payload["rows_failing_frozen"])
        self.assertEqual(1, payload["rows_passing_rescaled"])
        self.assertEqual(0, payload["rows_failing_rescaled"])
        self.assertEqual(0, payload["rows_newly_rejected"])
        self.assertEqual(0, payload["accessions_lost_all_support"])
        rows = self._by_accession()
        self.assertEqual("1", rows["ACC_ONE"]["frozen_support_rows"])
        self.assertEqual("0", rows["ACC_TWO"]["total_hit_rows"])

    def test_rows_outside_the_universe_are_counted_not_ignored(self):
        summary = self._summary()
        self.assertEqual(7, summary["rows_total"])
        self.assertEqual(6, summary["rows_joining_universe"])
        self.assertEqual(1, summary["rows_not_in_universe"])
        # the row claiming no universe accession must not move any row-level
        # counter: six joined rows, six frozen passes, four rescaled passes
        self.assertEqual(
            6, summary["rows_passing_frozen"] + summary["rows_failing_frozen"]
        )
        self.assertEqual(
            6, summary["rows_passing_rescaled"] + summary["rows_failing_rescaled"]
        )
        self.assertEqual(6, summary["rows_passing_frozen"])
        self.assertEqual(4, summary["rows_passing_rescaled"])
        # and the per-accession sums equal those row-level counters
        self.assertEqual(
            6,
            sum(int(row["frozen_support_rows"]) for row in self._by_accession().values()),
        )
        self.assertEqual(
            4,
            sum(
                int(row["rescaled_support_rows"])
                for row in self._by_accession().values()
            ),
        )


class DeltaReportCrossCheckTests(_ScaleImpactFixture):
    """When the frozen F7 delta report is supplied, the numbers must match."""

    def _write_report(self, **overrides):
        """A frozen-delta-report stand-in whose numbers match the hits fixture.

        The fixture is a seven-row table; six rows claim a universe accession
        and the seventh names nobody.  Those six carry six frozen per-shard
        passes and four corrected-scale passes (two rows are rejected by the
        longer-library scale), so the report says exactly that and
        ``overrides`` breaks one field at a time.
        """
        payload = {
            "z_total": 615_969_589,
            "total_rows": 7,
            "rows_passing_original_threshold": 6,
            "rows_passing_rescaled_threshold": 4,
            "delta_rows": -2,
            "threshold_rescaled": "1e-05",
            "threshold_original": "1e-05",
        }
        payload.update(overrides)
        path = self.dir / "scale_delta_report.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_matching_report_is_recorded(self):
        report = self._write_report()
        code, stdout, err = self._main(
            self._argv(
                "--delta-report",
                str(report),
                "--allow-no-hit-accessions",
            )
        )
        self.assertEqual(0, code, err)
        recorded = json.loads(stdout)["frozen_delta_cross_check"]
        self.assertEqual("match", recorded["status"])
        self.assertEqual(7, recorded["report_total_rows"])
        self.assertEqual(6, recorded["report_rows_passing_original_threshold"])
        self.assertIsNone(recorded["mismatches"])
        self.assertEqual(7, recorded["recomputed_total_rows"])
        self.assertEqual(6, recorded["recomputed_rows_passing_original_threshold"])
        self.assertEqual(4, recorded["recomputed_rows_passing_rescaled_threshold"])
        self.assertEqual(-2, recorded["recomputed_delta_rows"])

    def test_mismatching_report_fails_closed(self):
        report = self._write_report(
            rows_passing_rescaled_threshold=5,
            delta_rows=-1,
        )
        out = self.dir / "mismatch_out"
        code, stdout, err = self._main(
            self._argv(
                "--delta-report",
                str(report),
                "--allow-no-hit-accessions",
                out_dir=out,
            )
        )
        self.assertEqual(3, code)
        self.assertIn("rows_passing_rescaled_threshold", err)
        self.assertIn("delta_rows", err)
        self.assertFalse(out.exists())
    def test_limited_run_skips_the_cross_check_explicitly(self):
        report = self._write_report()
        code, stdout, err = self._main(
            self._argv(
                "--delta-report",
                str(report),
                "--limit",
                "3",
                "--allow-no-hit-accessions",
            )
        )
        self.assertEqual(0, code, err)
        recorded = json.loads(stdout)["frozen_delta_cross_check"]
        self.assertEqual("skipped_limited", recorded["status"])


class ModuleContractTests(_ScaleImpactFixture):
    """The module's documented contract."""

    def test_docstring_states_the_question_and_the_join_decision(self):
        doc = (self.module.__doc__ or "").lower()
        self.assertIn("passes_rescaled_threshold", doc)
        self.assertIn("protein_layers.tsv", doc)
        self.assertIn("985", doc)
        self.assertIn("join", doc)

    def test_threshold_rule_is_the_project_rule(self):
        self.assertEqual("E < threshold (strict)", self.module.THRESHOLD_RULE)

    def test_pass_value_accepts_only_the_documented_spellings(self):
        self.assertTrue(self.module.pass_value("true"))
        self.assertTrue(self.module.pass_value("TRUE"))
        self.assertFalse(self.module.pass_value("false"))
        self.assertFalse(self.module.pass_value("FALSE"))
        for bad in ("pending", "", "0", "1", "yes", None):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    self.module.pass_value(bad)

    def test_exit_codes_are_named(self):
        self.assertEqual(0, self.module.EXIT_OK)
        self.assertEqual(2, self.module.EXIT_USAGE)
        self.assertEqual(3, self.module.EXIT_INCONSISTENT)
        self.assertEqual(4, self.module.EXIT_REFUSED)


if __name__ == "__main__":
    unittest.main()
