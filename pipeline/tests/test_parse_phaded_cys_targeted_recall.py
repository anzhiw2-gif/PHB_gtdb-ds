"""Tests for the bounded Cys targeted-recall parser (2026-09-17, Step 5).

Written before the implementation (test-first).  The parser scores a *bounded*
existing GTDB protein pool with the uncalibrated Cys discovery-layer HMM,
rebuilds the recall tables, and quantifies how many discovery hits fall
**outside** the frozen 109,087-protein candidate universe.

Two invariants are load-bearing and must never be relaxed:

* the discovery threshold is frozen at ``E < 1e-5`` (identical to the
  pre-registered run) and ``classify_hits`` refuses any other value;
* an output row never carries a family call and never changes a
  ``subtype_call``; every row is labelled ``discovery_hmm_uncalibrated``.

Because the discovery HMM has no discriminative power (it hits 443/563 of the
x1 lipase/esterase suspicion set), every recalled row is a *review priority
queue* entry, never specificity evidence.  The parser therefore must not
fabricate Pfam/lipase-box annotations for proteins that no existing annotation
table covers: those fields are ``not_available``.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "parse_phaded_cys_targeted_recall.py"
)

TABLOUT_HEADER = (
    "# hmmsearch :: search profile(s) against a sequence database\n"
    "# target name        accession  query name  accession  E-value  score  bias"
    "  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to  from  to  acc description\n"
    "#---\n"
)

MOTIF_HEADER = [
    "accession", "genome", "candidate_prior_superfamily", "candidate_assignment_status",
    "sbd_pf06850_binding_state", "lipase_box_state",
]


def load_module():
    spec = importlib.util.spec_from_file_location(
        "parse_phaded_cys_targeted_recall", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def write_fasta(path, records):
    """records: list of (id, sequence) or (id, sequence, description)."""
    chunks = []
    for record in records:
        if len(record) == 3:
            identifier, sequence, description = record
            chunks.append(">%s %s\n%s\n" % (identifier, description, sequence))
        else:
            identifier, sequence = record
            chunks.append(">%s\n%s\n" % (identifier, sequence))
    write_text(path, "".join(chunks))


def tblout_lines(records):
    """records: list of (target, evalue, score)."""
    out = [TABLOUT_HEADER]
    for target, evalue, score in records:
        out.append(
            "%-40s %-10s %-12s %-10s %9s %6s %5s %3d %3d %9s %9s %6s %5s %4d %4d %4d %4d %4d %4d %4s %s\n"
            % (target, "-", "cys_discovery", "-", evalue, score, "0.1", 1, 1,
               evalue, evalue, score, "0.1", 1, 10, 1, 10, 1, 10, "0.95", "desc")
        )
    out.append("# Program:         hmmsearch\n")
    out.append("# [ok]\n")
    return "".join(out)


class TbloutParsingTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_parse_tblout_keeps_best_hit_per_target(self):
        path = self.root / "a.tbl"
        write_text(path, tblout_lines([
            ("GCA_1|CTX_1", "1e-30", "120.0"),
            ("GCA_1|CTX_1", "3e-40", "150.0"),
            ("GCA_2|CTX_9", "0.5", "12.0"),
        ]))
        parsed = self.module.parse_tblout(path)
        self.assertEqual(parsed["GCA_1|CTX_1"]["best_e"], 3e-40)
        self.assertEqual(parsed["GCA_1|CTX_1"]["best_score"], 150.0)
        self.assertEqual(parsed["GCA_1|CTX_1"]["n_reported_rows"], 2)
        self.assertEqual(parsed["GCA_2|CTX_9"]["best_e"], 0.5)

    def test_parse_tblout_does_not_expand_pipe_aliases(self):
        """GTDB protein ids contain '|' *as part of the identifier*.

        Unlike the frozen reference panel (``>reference_id|accession``), the GTDB
        headers are ``>GCA_xxxx|contig_N``; splitting on ``|`` would create
        spurious keys and could silently merge unrelated records.
        """
        path = self.root / "b.tbl"
        write_text(path, tblout_lines([("GCA_000008085.1|AE017199.1_1", "1e-9", "60.0")]))
        parsed = self.module.parse_tblout(path)
        self.assertIn("GCA_000008085.1|AE017199.1_1", parsed)
        self.assertNotIn("GCA_000008085.1", parsed)
        self.assertNotIn("AE017199.1_1", parsed)

    def test_parse_tblout_skips_comments_and_short_rows(self):
        path = self.root / "c.tbl"
        write_text(path, TABLOUT_HEADER + "\n" + "x y\n" + "# trailing\n")
        self.assertEqual(self.module.parse_tblout(path), {})

    def test_unreported_target_has_no_entry(self):
        path = self.root / "d.tbl"
        write_text(path, tblout_lines([("GCA_1|CTX_1", "1e-9", "60.0")]))
        parsed = self.module.parse_tblout(path)
        self.assertNotIn("GCA_2|CTX_2", parsed)


class ThresholdTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_constants_are_pinned(self):
        self.assertEqual(self.module.DISCOVERY_EVALUE_THRESHOLD, 1e-5)
        self.assertEqual(self.module.MODEL_LAYER, "discovery_hmm_uncalibrated")

    def test_classify_hits_refuses_any_other_threshold(self):
        with self.assertRaises(ValueError):
            self.module.classify_hits({"x": {"best_e": 1e-9}}, threshold=1e-3)
        with self.assertRaises(ValueError):
            self.module.classify_hits({"x": {"best_e": 1e-9}}, threshold=1e-6)

    def test_hit_requires_strictly_less_than_threshold(self):
        classified = self.module.classify_hits(
            {
                "below": {"best_e": 9.9e-6},
                "exactly": {"best_e": 1e-5},
                "above": {"best_e": 1.1e-5},
            },
            threshold=self.module.DISCOVERY_EVALUE_THRESHOLD,
        )
        self.assertEqual(set(classified["hits"]), {"below"})
        self.assertEqual(set(classified["non_hits"]), {"exactly", "above"})


class EvalueScaleTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_recall_pool_z_is_pinned_to_the_preregistered_scale(self):
        self.assertEqual(self.module.RECALL_POOL_Z, 109087)

    def test_rescale_evalue_is_linear_in_database_size(self):
        self.assertAlmostEqual(
            self.module.rescale_evalue(1e-10, from_z=109087, to_z=218174), 2e-10
        )
        self.assertAlmostEqual(
            self.module.rescale_evalue(1e-10, from_z=109087, to_z=109087), 1e-10
        )

    def test_rescale_rejects_non_positive_z(self):
        with self.assertRaises(ValueError):
            self.module.rescale_evalue(1e-10, from_z=0, to_z=10)
        with self.assertRaises(ValueError):
            self.module.rescale_evalue(1e-10, from_z=10, to_z=0)


class RowBuildingTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _hits(self):
        path = self.root / "hits.tbl"
        write_text(path, tblout_lines([
            ("GCA_1|CTX_1", "1e-30", "120.0"),
            ("GCA_2|CTX_2", "1e-7", "80.0"),
            ("GCA_3|CTX_3", "1e-4", "40.0"),
        ]))
        return self.module.parse_tblout(path)

    def test_hit_rows_mark_pool_membership(self):
        rows = self.module.build_hit_rows(
            pool_ids=["GCA_1|CTX_1", "GCA_2|CTX_2", "GCA_3|CTX_3"],
            tblout_hits=self._hits(),
            candidate_universe={"GCA_1|CTX_1"},
            motif_rows={},
            tag_layer=set(),
        )
        by_id = {row["protein_id"]: row for row in rows}
        self.assertEqual(len(rows), 2)
        self.assertEqual(by_id["GCA_1|CTX_1"]["in_existing_candidate_pool"], "true")
        self.assertEqual(by_id["GCA_2|CTX_2"]["in_existing_candidate_pool"], "false")

    def test_hit_rows_carry_only_boundary_declarations(self):
        rows = self.module.build_hit_rows(
            pool_ids=["GCA_1|CTX_1"],
            tblout_hits=self._hits(),
            candidate_universe=set(),
            motif_rows={},
            tag_layer=set(),
        )
        self.assertEqual(rows[0]["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(rows[0]["family_call_made"], "false")
        self.assertEqual(rows[0]["subtype_call_impact"], "none")

    def test_hit_rows_are_sorted_by_protein_id(self):
        rows = self.module.build_hit_rows(
            pool_ids=["GCA_1|CTX_1", "GCA_2|CTX_2"],
            tblout_hits=self._hits(),
            candidate_universe=set(),
            motif_rows={},
            tag_layer=set(),
        )
        self.assertEqual(
            [row["protein_id"] for row in rows], ["GCA_1|CTX_1", "GCA_2|CTX_2"]
        )

    def test_in_pool_hits_take_annotations_from_the_motif_table(self):
        motif_rows = {
            "GCA_1|CTX_1": {
                "accession": "GCA_1|CTX_1",
                "candidate_assignment_status": "unassigned_PhaDED_like",
                "candidate_prior_superfamily": "",
                "sbd_pf06850_binding_state": "detected",
                "lipase_box_state": "not_detected_pattern",
            }
        }
        rows = self.module.build_hit_rows(
            pool_ids=["GCA_1|CTX_1"],
            tblout_hits=self._hits(),
            candidate_universe={"GCA_1|CTX_1"},
            motif_rows=motif_rows,
            tag_layer={"GCA_1|CTX_1"},
        )
        self.assertEqual(rows[0]["lipase_box_state"], "not_detected_pattern")
        self.assertEqual(rows[0]["sbd_pf06850_binding_state"], "detected")
        self.assertEqual(rows[0]["in_task6_tag_layer"], "true")

    def test_out_of_pool_hits_must_not_fabricate_annotations(self):
        rows = self.module.build_hit_rows(
            pool_ids=["GCA_2|CTX_2"],
            tblout_hits=self._hits(),
            candidate_universe=set(),
            motif_rows={},
            tag_layer=set(),
        )
        self.assertEqual(rows[0]["lipase_box_state"], self.module.NOT_AVAILABLE)
        self.assertEqual(rows[0]["sbd_pf06850_binding_state"], self.module.NOT_AVAILABLE)
        self.assertEqual(rows[0]["candidate_assignment_status"], self.module.NOT_AVAILABLE)

    def test_new_candidate_rows_are_exactly_hits_outside_the_pool(self):
        hit_rows = self.module.build_hit_rows(
            pool_ids=["GCA_1|CTX_1", "GCA_2|CTX_2"],
            tblout_hits=self._hits(),
            candidate_universe={"GCA_1|CTX_1"},
            motif_rows={},
            tag_layer=set(),
        )
        new_rows = self.module.build_new_candidate_rows(hit_rows)
        self.assertEqual([row["protein_id"] for row in new_rows], ["GCA_2|CTX_2"])
        for row in new_rows:
            self.assertEqual(row["in_existing_candidate_pool"], "false")
            self.assertEqual(row["task6_tag_layer_claimed"], "false")

    def test_new_candidate_rows_declare_a_definitional_zero_tag_overlap(self):
        hit_rows = self.module.build_hit_rows(
            pool_ids=["GCA_2|CTX_2"],
            tblout_hits=self._hits(),
            candidate_universe=set(),
            motif_rows={},
            tag_layer={"GCA_2|CTX_2"},
        )
        new_rows = self.module.build_new_candidate_rows(hit_rows)
        self.assertEqual(new_rows[0]["task6_tag_layer_claimed"], "false")
        self.assertEqual(new_rows[0]["task6_overlap_basis"], self.module.DEFINITIONAL_ZERO)

    def test_new_candidate_rows_are_review_priority_only(self):
        hit_rows = self.module.build_hit_rows(
            pool_ids=["GCA_2|CTX_2"],
            tblout_hits=self._hits(),
            candidate_universe=set(),
            motif_rows={},
            tag_layer=set(),
        )
        new_rows = self.module.build_new_candidate_rows(hit_rows)
        self.assertEqual(new_rows[0]["review_priority_only"], "true")
        self.assertEqual(new_rows[0]["discriminative_power_caveat"], self.module.NO_DISCRIMINATION_NOTE)
        self.assertEqual(new_rows[0]["family_call_made"], "false")
        self.assertEqual(new_rows[0]["new_family_call_made"], "false")
        self.assertEqual(new_rows[0]["subtype_call_impact"], "none")


class SummaryTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def _summary(self, **overrides):
        kwargs = dict(
            pool_size=1000,
            pool_scope="full_pool",
            sampling_statement="",
            hits=200,
            new_candidates=150,
            in_pool_hits=50,
            tag_layer_size=40,
            new_candidate_tag_overlap=0,
            strict_hits_at_pool_z=180,
            pool_z=1000,
            scored_from_annotation_table=50,
        )
        kwargs.update(overrides)
        return self.module.build_summary(**kwargs)

    def test_total_pool_z_is_reported(self):
        self.assertEqual(self._summary()["hmmsearch_database_size_Z"], 109087)

    def test_scope_and_sampling_statement_are_recorded(self):
        summary = self._summary()
        self.assertEqual(summary["pool_scope"], "full_pool")
        self.assertEqual(summary["sampling_declared"], False)
        self.assertEqual(summary["is_full_library_rescan"], False)

    def test_sampled_scope_requires_a_sampling_statement(self):
        with self.assertRaises(ValueError):
            self._summary(pool_scope="declared_shard_sample", sampling_statement="")
        summary = self._summary(
            pool_scope="declared_shard_sample",
            sampling_statement="uniform 20 of 100 shards",
        )
        self.assertEqual(summary["sampling_declared"], True)

    def test_unknown_scope_is_rejected(self):
        with self.assertRaises(ValueError):
            self._summary(pool_scope="whatever")

    def test_counts_and_rates_are_consistent(self):
        summary = self._summary()
        self.assertEqual(summary["counts"]["recall_pool_proteins_scored"], 1000)
        self.assertEqual(summary["counts"]["recall_pool_hits"], 200)
        self.assertEqual(summary["counts"]["new_cys_candidates"], 150)
        self.assertAlmostEqual(summary["hit_rate"], 0.2)
        self.assertAlmostEqual(summary["new_candidate_share_of_hits"], 0.75)

    def test_strict_variant_never_exceeds_the_reported_hit_count(self):
        summary = self._summary()
        self.assertEqual(summary["counts"]["strict_hits_at_actual_pool_z"], 180)
        self.assertLessEqual(
            summary["counts"]["strict_hits_at_actual_pool_z"],
            summary["counts"]["recall_pool_hits"],
        )
        self.assertEqual(summary["hmmsearch_database_size_Z"], 109087)
        self.assertEqual(summary["actual_pool_size_z"], 1000)

    def test_boundary_flags_are_always_false(self):
        summary = self._summary()
        self.assertFalse(summary["family_call_made"])
        self.assertFalse(summary["new_family_call_made"])
        self.assertEqual(summary["subtype_call_rows_changed"], 0)
        self.assertFalse(summary["calibration_gate_met"])
        self.assertTrue(summary["calibration_gate_unchanged"])
        self.assertEqual(summary["model_layer"], "discovery_hmm_uncalibrated")


class EndToEndTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.run_dir = self.root / "runs" / "20260917_phaded_cys_targeted_recall_01"
        (self.run_dir / "results").mkdir(parents=True)
        (self.run_dir / "inputs").mkdir(parents=True)
        (self.run_dir / "logs").mkdir(parents=True)

        self.pool = self.run_dir / "inputs" / "recall_pool.faa"
        write_fasta(self.pool, [
            ("GCA_1|CTX_1", "MKV", "first"),
            ("GCA_2|CTX_2", "MKV"),
            ("GCA_3|CTX_3", "MKV"),
        ])
        self.candidates = self.run_dir / "inputs" / "candidate_union.faa"
        write_fasta(self.candidates, [("GCA_1|CTX_1", "MKV")])
        self.tblout = self.run_dir / "results" / "recall_pool_threshold.tblout"
        write_text(self.tblout, tblout_lines([
            ("GCA_1|CTX_1", "1e-30", "120.0"),
            ("GCA_2|CTX_2", "1e-8", "70.0"),
            ("GCA_3|CTX_3", "1e-4", "30.0"),
        ]))
        self.motif = self.run_dir / "inputs" / "motif_completion_full.tsv"
        write_text(
            self.motif,
            "\t".join(MOTIF_HEADER) + "\n"
            + "\t".join([
                "GCA_1|CTX_1", "GCA_1", "", "unassigned_PhaDED_like",
                "detected", "not_detected_pattern",
            ]) + "\n",
        )
        self.snapshot = self.run_dir / "logs" / "server_load_snapshot.json"
        write_text(self.snapshot, json.dumps({"load_average": "0.5 0.5 0.5"}))

    def tearDown(self):
        self.tmp.cleanup()

    def _args(self, **overrides):
        args = self.module.build_arg_namespace([
            "--run-dir", str(self.run_dir),
            "--pool-faa", str(self.pool),
            "--tblout", str(self.tblout),
            "--candidate-faa", str(self.candidates),
            "--motif-tsv", str(self.motif),
            "--pool-scope", "full_pool",
            "--actual-pool-size-z", "3",
            "--server-snapshot", str(self.snapshot),
            "--tool", "hmmsearch=HMMER 3.4 (Aug 2023)",
            "--command", "hmmsearch --cpu 40 -E 1e-5 -Z 109087 ...",
            "--input-record", "cys_discovery_hmm=%s" % (self.run_dir / "inputs" / "hmm.hmm"),
        ])
        for key, value in overrides.items():
            setattr(args, key, value)
        return args

    def test_run_writes_the_four_required_outputs(self):
        write_text(self.run_dir / "inputs" / "hmm.hmm", "HMMER3/f\n")
        payload = self.module.run(self._args())
        for name in (
            "recall_hits.tsv",
            "new_cys_candidates.tsv",
            "recall_summary.json",
            "recall_manifest.json",
        ):
            self.assertTrue((self.run_dir / "results" / name).is_file(), name)
        self.assertEqual(payload["counts"]["recall_pool_hits"], 2)
        self.assertEqual(payload["counts"]["new_cys_candidates"], 1)

    def test_manifest_records_hashes_and_pending_for_missing_inputs(self):
        manifest = self.module.run(self._args())["manifest"]
        self.assertEqual(
            manifest["inputs"]["cys_discovery_hmm"]["sha256"], self.module.sha256_file(
                self.run_dir / "inputs" / "hmm.hmm"
            ) if (self.run_dir / "inputs" / "hmm.hmm").is_file() else None
        )
        self.assertEqual(manifest["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(manifest["threads"]["max_single_task_threads_allowed"], 40)
        self.assertFalse(manifest["registry_modified"])

    def test_run_is_fail_closed_on_an_empty_pool(self):
        write_fasta(self.pool, [])
        with self.assertRaises(ValueError):
            self.module.run(self._args())

    def test_run_is_fail_closed_when_no_pool_record_is_reported(self):
        write_text(self.tblout, tblout_lines([]))
        with self.assertRaises(ValueError):
            self.module.run(self._args())

    def test_outputs_declare_no_call_columns(self):
        write_text(self.run_dir / "inputs" / "hmm.hmm", "HMMER3/f\n")
        self.module.run(self._args())
        scan = self.module.scan_output_call_columns(
            scores=self.module.RECALL_HIT_FIELDS,
            claims=self.module.NEW_CANDIDATE_FIELDS,
        )
        self.assertEqual(scan["call_columns"], [])
        self.assertTrue(scan["boundary_declaration_columns"])


if __name__ == "__main__":
    unittest.main()
