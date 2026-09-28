"""Tests for the Cys discovery-HMM score parser and gate evaluator
(2026-09-17, Step 7).

Written before the implementation (test-first).  The parser consumes the
archived HMMER tblout files, rebuilds the 109,087-row score table, derives the
unassigned-claim list, and evaluates the four pre-registered gates G1-G4
without ever lowering the pre-registered threshold.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "parse_phaded_cys_discovery_hmm_scores.py"

CYS = "intracellular nPHASCL without lipase box"

TABLOUT_HEADER = (
    "# hmmsearch :: search profile(s) against a sequence database\n"
    "# target name        accession  query name  accession  E-value  score  bias"
    "  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to  from  to  acc description\n"
    "#---\n"
)

LEDGER_HEADER = [
    "reference_id", "accession", "phaded_superfamily", "phaded_family_id",
    "reported_localization", "evidence_status",
]

MOTIF_HEADER = [
    "accession", "genome", "candidate_prior_superfamily", "candidate_assignment_status",
    "sbd_pf06850_binding_state", "lipase_box_state",
]


def load_module():
    spec = importlib.util.spec_from_file_location(
        "parse_phaded_cys_discovery_hmm_scores", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def tblout_lines(records):
    """records: list of (target, accession, evalue, score, desc)."""
    out = [TABLOUT_HEADER]
    for target, accession, evalue, score, desc in records:
        out.append(
            "%-24s %-10s %-12s %-10s %9s %6s %5s %3d %3d %9s %9s %6s %5s %4d %4d %4d %4d %4d %4d %4s %s\n"
            % (target, accession, "cys_discovery", "-", evalue, score, "0.1", 1, 1,
               evalue, evalue, score, "0.1", 1, 10, 1, 10, 1, 10, "0.95", desc)
        )
    out.append("# Program:         hmmsearch\n")
    return "".join(out)


class TbloutParserTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_parse_tblout_keeps_best_hit_per_target(self):
        path = self.root / "a.tbl"
        write_text(path, tblout_lines([
            ("candA", "-", "1e-30", "120.0", "desc"),
            ("candA", "-", "3e-40", "150.0", "desc"),
            ("candB", "-", "0.5", "12.0", "desc"),
        ]))
        parsed = self.module.parse_tblout(path)
        self.assertEqual(parsed["candA"]["best_e"], 3e-40)
        self.assertEqual(parsed["candA"]["best_score"], 150.0)
        self.assertEqual(parsed["candB"]["best_e"], 0.5)
        # hmmsearch --tblout emits one row per target sequence, so this counter
        # is "reported rows", not a domain count
        self.assertEqual(parsed["candA"]["n_reported_rows"], 2)
        self.assertNotIn("n_domains", parsed["candA"])

    def test_confounder_burden_counts_real_records_not_alias_keys(self):
        records = ["GCA_1|CTX1", "GCA_2|CTX2", "GCA_3|CTX3"]
        path = self.root / "conf.tbl"
        write_text(path, tblout_lines([
            ("GCA_1|CTX1", "-", "1e-9", "60.0", "desc"),
            ("GCA_2|CTX2", "-", "0.5", "9.0", "desc"),
        ]))
        burden = self.module.evaluate_confounder_burden(
            records=records, hit_sources=[self.module.parse_tblout(path)])
        self.assertEqual(burden["records"], 3)
        self.assertEqual(burden["hits"], 1)
        self.assertAlmostEqual(burden["hit_rate"], 1 / 3)
        self.assertEqual([entry["accession"] for entry in burden["hit_list"]], ["GCA_1|CTX1"])

    def test_discovery_threshold_is_pinned_and_not_adjustable(self):
        self.assertEqual(self.module.DISCOVERY_EVALUE_THRESHOLD, 1e-5)
        self.assertEqual(self.module.MODEL_LAYER, "discovery_hmm_uncalibrated")
        with self.assertRaises(ValueError):
            self.module.classify_hits({"x": {"best_e": 1e-4}}, threshold=1e-3)

    def test_parse_tblout_indexes_pipe_delimited_name_aliases(self):
        path = self.root / "aliases.tbl"
        write_text(path, tblout_lines([
            ("DED_hfam_57_0038|CAJ96855.1", "-", "1e-9", "60.0", "desc"),
            ("Q84C08|mcc_control", "-", "0.4", "8.0", "desc"),
        ]))
        parsed = self.module.parse_tblout(path)
        for alias in ("DED_hfam_57_0038|CAJ96855.1", "CAJ96855.1", "DED_hfam_57_0038",
                      "Q84C08|mcc_control", "Q84C08"):
            self.assertIn(alias, parsed)
        self.assertEqual(parsed["CAJ96855.1"]["best_e"], 1e-9)
        self.assertEqual(parsed["CAJ96855.1"]["n_reported_rows"], 1)
        self.assertEqual(parsed["Q84C08"]["best_e"], 0.4)

    def test_parse_tblout_does_not_duplicate_counts_through_aliases(self):
        path = self.root / "dup.tbl"
        write_text(path, tblout_lines([
            ("DED_hfam_61_0001|BAF86293.1", "-", "1e-30", "120.0", "desc"),
            ("DED_hfam_61_0001|BAF86293.1", "-", "1e-40", "150.0", "desc"),
        ]))
        parsed = self.module.parse_tblout(path)
        self.assertEqual(parsed["BAF86293.1"]["n_reported_rows"], 2)
        self.assertEqual(parsed["BAF86293.1"]["best_e"], 1e-40)

    def test_scores_table_covers_every_candidate_and_marks_hits(self):
        candidates = ["candA", "candB", "candC"]
        primary = {"candA": {"best_e": 1e-30, "best_score": 120.0, "n_domains": 1}}
        supplementary = {
            "candA": {"best_e": 1e-30, "best_score": 120.0, "n_domains": 1},
            "candB": {"best_e": 0.4, "best_score": 11.0, "n_domains": 1},
            "candC": {"best_e": 9.9, "best_score": -3.0, "n_domains": 0},
        }
        rows = self.module.build_score_rows(candidates, primary, supplementary)
        self.assertEqual(len(rows), 3)
        by_accession = {r["accession"]: r for r in rows}
        self.assertEqual(by_accession["candA"]["discovery_hit"], "true")
        self.assertEqual(by_accession["candA"]["score_source"], "pre_registered_threshold_run")
        self.assertEqual(by_accession["candB"]["discovery_hit"], "false")
        self.assertEqual(by_accession["candB"]["score_source"], "supplementary_unfiltered_run")
        self.assertEqual(by_accession["candC"]["score_source"], "supplementary_unfiltered_run")
        for row in rows:
            self.assertEqual(row["model_layer"], "discovery_hmm_uncalibrated")
            self.assertEqual(row["family_call_made"], "false")
            self.assertEqual(row["subtype_call_impact"], "none")

    def test_build_score_rows_rejects_missing_candidate(self):
        with self.assertRaises(ValueError):
            self.module.build_score_rows(["candA"], {}, {})


class GateEvaluationTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _reference_fixture(self):
        ledger = self.root / "ledger.tsv"
        rows = [
            {"reference_id": "DED_hfam_65_0001", "accession": "CAJ92291.1", "phaded_superfamily": CYS,
             "phaded_family_id": "DED_hfam_65", "reported_localization": "intracellular",
             "evidence_status": "experimental_positive"},
            {"reference_id": "DED_hfam_61_0001", "accession": "BAF86293.1", "phaded_superfamily": CYS,
             "phaded_family_id": "DED_hfam_61", "reported_localization": "intracellular",
             "evidence_status": "annotation_only"},
            {"reference_id": "DED_hfam_2_0001", "accession": "AAT62319.1",
             "phaded_superfamily": "intracellular nPHASCL with lipase box",
             "phaded_family_id": "DED_hfam_2", "reported_localization": "intracellular",
             "evidence_status": "annotation_only"},
            {"reference_id": "DED_hfam_57_0038", "accession": "CAJ96855.1",
             "phaded_superfamily": "extracellular dPHASCL type 1",
             "phaded_family_id": "DED_hfam_57", "reported_localization": "extracellular",
             "evidence_status": "annotation_only"},
        ]
        lines = ["\t".join(LEDGER_HEADER)]
        for row in rows:
            lines.append("\t".join(row[c] for c in LEDGER_HEADER))
        write_text(ledger, "\n".join(lines) + "\n")
        return rows

    def test_g1_passes_when_no_extracellular_hit_and_counts_lipase_box_cross_talk(self):
        ledger_rows = self._reference_fixture()
        reference_hits = {
            "CAJ92291.1": {"best_e": 1e-200, "best_score": 600.0},
            "BAF86293.1": {"best_e": 1e-150, "best_score": 500.0},
            "AAT62319.1": {"best_e": 1e-9, "best_score": 60.0},
        }
        gate = self.module.evaluate_g1(
            reference_hits=reference_hits,
            ledger_rows=ledger_rows,
            detected_controls={"Q84C08": {"best_e": 0.7, "best_score": 5.0}},
        )
        self.assertEqual(gate["gate_id"], "G1")
        self.assertEqual(gate["unexplained_hits"], 0)
        self.assertEqual(gate["outcome"], "passed")
        self.assertEqual(gate["explained_cross_talk_count"], 1)
        self.assertEqual(gate["explained_cross_talk"][0]["accession"], "AAT62319.1")
        self.assertEqual(gate["population_size"], 2, "1 extracellular reference + Q84C08")

    def test_g1_fails_when_an_extracellular_reference_is_hit(self):
        ledger_rows = self._reference_fixture()
        reference_hits = {"CAJ96855.1": {"best_e": 1e-12, "best_score": 80.0}}
        gate = self.module.evaluate_g1(
            reference_hits=reference_hits, ledger_rows=ledger_rows, detected_controls={})
        self.assertEqual(gate["unexplained_hits"], 1)
        self.assertEqual(gate["outcome"], "failed")
        self.assertEqual(gate["threshold"], "unexplained_hits == 0")

    def test_g1_fails_when_q84c08_is_hit(self):
        ledger_rows = self._reference_fixture()
        gate = self.module.evaluate_g1(
            reference_hits={}, ledger_rows=ledger_rows,
            detected_controls={"Q84C08": {"best_e": 1e-11, "best_score": 70.0}})
        self.assertEqual(gate["unexplained_hits"], 1)
        self.assertEqual(gate["outcome"], "failed")
        self.assertEqual(gate["unexplained_hit_list"][0]["accession"], "Q84C08")

    def test_g2_uses_integrity_passing_denominator_and_labels_in_sample(self):
        gate = self.module.evaluate_g2(
            reference_hits={"CAJ92291.1": {"best_e": 1e-200}, "BAF86293.1": {"best_e": 0.5}},
            training_accessions=["CAJ92291.1", "BAF86293.1"],
            threshold=0.90,
        )
        self.assertEqual(gate["recall"], 0.5)
        self.assertEqual(gate["outcome"], "failed")
        self.assertEqual(gate["measurement_scope"], "in_sample_not_held_out")

    def test_g3_requires_byte_identical_formal_scan_models(self):
        gate = self.module.evaluate_g3(before_sha256="a" * 64, after_sha256="a" * 64,
                                       labelled_outputs=["x", "y"], unlabelled_outputs=[])
        self.assertEqual(gate["outcome"], "passed")
        gate2 = self.module.evaluate_g3(before_sha256="a" * 64, after_sha256="b" * 64,
                                        labelled_outputs=["x"], unlabelled_outputs=[])
        self.assertEqual(gate2["outcome"], "failed")

    def test_g4_requires_no_family_or_subtype_column(self):
        gate = self.module.evaluate_g4(
            score_columns=["accession", "best_E", "discovery_hit"],
            claim_columns=["accession", "candidate_assignment_status"],
            subtype_call_rows_changed=0,
        )
        self.assertEqual(gate["outcome"], "passed")
        bad = self.module.evaluate_g4(
            score_columns=["accession", "family_call"],
            claim_columns=["accession"], subtype_call_rows_changed=0)
        self.assertEqual(bad["outcome"], "failed")

    def test_scan_labelled_outputs_flags_untagged_files(self):
        tagged = self.root / "scores.tsv"
        write_text(tagged, "accession\tmodel_layer\nx\tdiscovery_hmm_uncalibrated\n")
        hmm = self.root / "cys_discovery.hmm"
        write_text(hmm, "HMMER3/f\nNAME  cys_discovery_uncalibrated\n")
        untagged = self.root / "legacy.tsv"
        write_text(untagged, "accession\nx\n")
        labelled, unlabelled = self.module.scan_labelled_outputs([tagged, hmm, untagged])
        self.assertEqual(labelled, [str(tagged), str(hmm)])
        self.assertEqual(unlabelled, [str(untagged)])
        gate = self.module.evaluate_g3(before_sha256="a" * 64, after_sha256="a" * 64,
                                       labelled_outputs=labelled, unlabelled_outputs=unlabelled)
        self.assertEqual(gate["outcome"], "failed")

    def test_g1_detects_hit_reached_through_reference_header_alias(self):
        ledger_rows = self._reference_fixture()
        path = self.root / "ref_alias.tbl"
        write_text(path, tblout_lines([("DED_hfam_57_0038|CAJ96855.1", "-", "1e-12", "80.0", "d")]))
        reference_hits = self.module.parse_tblout(path)
        gate = self.module.evaluate_g1(reference_hits=reference_hits, ledger_rows=ledger_rows,
                                       detected_controls={})
        self.assertEqual(gate["unexplained_hits"], 1)
        self.assertEqual(gate["outcome"], "failed")
        self.assertEqual(gate["unexplained_hit_list"][0]["accession"], "CAJ96855.1")

    def test_g2_recall_is_measured_through_reference_header_aliases(self):
        path = self.root / "g2_alias.tbl"
        write_text(path, tblout_lines([
            ("DED_hfam_65_0001|CAJ92291.1", "-", "1e-200", "700.0", "d"),
            ("DED_hfam_61_0001|BAF86293.1", "-", "0.5", "9.0", "d"),
        ]))
        gate = self.module.evaluate_g2(
            reference_hits=self.module.parse_tblout(path),
            training_accessions=["CAJ92291.1", "BAF86293.1"], threshold=0.90)
        self.assertEqual(gate["recall"], 0.5)
        self.assertEqual(gate["outcome"], "failed")

    def test_g4_does_not_flag_boundary_declaration_columns(self):
        gate = self.module.evaluate_g4(
            score_columns=["accession", "family_call_made", "subtype_call_impact"],
            claim_columns=["accession", "new_family_call_made", "candidate_assignment_status"],
            subtype_call_rows_changed=0,
        )
        self.assertEqual(gate["outcome"], "passed")
        self.assertEqual(gate["call_columns"], [])
        self.assertEqual(gate["boundary_declaration_columns"],
                         ["family_call_made", "new_family_call_made", "subtype_call_impact"])
        # the literal token scan is still reported, so the observation is not hidden
        self.assertEqual(gate["literal_token_columns"],
                         ["family_call_made", "new_family_call_made", "subtype_call_impact"])

    def test_claim_rows_join_assignment_status_and_tag_layer_overlap(self):
        claims = self.module.build_claim_rows(
            hit_rows=[{"accession": "candA", "best_E": "1e-30", "best_score": "120.0"}],
            motif_rows={"candA": {"candidate_assignment_status": "unassigned_PhaDED_like",
                                  "candidate_prior_superfamily": ""}},
            tag_layer_accessions={"candA", "other"},
        )
        self.assertEqual(len(claims), 1)
        row = claims[0]
        self.assertEqual(row["candidate_assignment_status"], "unassigned_PhaDED_like")
        self.assertEqual(row["task6_tag_layer_claimed"], "true")
        self.assertEqual(row["subtype_call_impact"], "none")
        self.assertEqual(row["new_family_call_made"], "false")


if __name__ == "__main__":
    unittest.main()
