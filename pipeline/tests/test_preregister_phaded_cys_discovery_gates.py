"""Tests for the pre-registered falsification gates of the Cys discovery HMM
(2026-09-17, Step 4).

Written before the implementation (test-first).  The pre-registration must be
written *before* training and scoring, must carry a UTC/local timestamp, must
pin the discovery threshold and the HMMER ``-Z`` comparison count, and must
fail closed if training artefacts already exist in the run directory.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "preregister_phaded_cys_discovery_gates.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "preregister_phaded_cys_discovery_gates", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PreregistrationTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.run_dir = self.root / "runs" / "20260917_phaded_cys_discovery_hmm_01"
        (self.run_dir / "inputs").mkdir(parents=True)
        (self.run_dir / "logs").mkdir()
        (self.run_dir / "results").mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def test_gates_cover_g1_to_g4_with_thresholds_and_fallbacks(self):
        gates = self.module.gate_definitions()
        ids = [g["gate_id"] for g in gates]
        self.assertEqual(ids, ["G1", "G2", "G3", "G4"])
        for gate in gates:
            self.assertTrue(gate["statement"])
            self.assertTrue(gate["metric"])
            self.assertTrue(gate["expected_observation"])
            self.assertTrue(gate["fallback_action"])
            self.assertIn(gate["gate_id"], {"G1", "G2", "G3", "G4"})
        by_id = {g["gate_id"]: g for g in gates}
        self.assertEqual(by_id["G1"]["threshold"], "unexplained_hits == 0")
        self.assertEqual(by_id["G2"]["threshold"], "recall >= 0.90")
        self.assertEqual(by_id["G4"]["threshold"], "subtype_call_rows_changed == 0")

    def test_preregistration_pins_threshold_and_database_size(self):
        payload = self.module.write_preregistration(self.run_dir, formal_scan_models_sha256="abc123")
        self.assertEqual(payload["discovery_evalue_threshold"], "1e-05")
        self.assertEqual(payload["hmmsearch_database_size_Z"], 109087)
        self.assertEqual(payload["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(payload["threshold_lowering_permitted"], False)
        self.assertEqual(payload["calibration_gate_unchanged"], True)
        self.assertEqual(payload["formal_scan_models_tsv_sha256_before"], "abc123")
        self.assertTrue(payload["registered_at_utc"])
        self.assertTrue(payload["registered_at_local"])
        self.assertEqual(payload["registered_before_training"], "true")

    def test_g1_declares_extracellular_panel_and_explained_cross_talk(self):
        payload = self.module.write_preregistration(self.run_dir, formal_scan_models_sha256="abc123")
        g1 = payload["gates"][0]
        self.assertIn("Q84C08", g1["population"])
        self.assertIn("365", g1["population"])
        self.assertIn("29", g1["explained_cross_talk_population"])
        self.assertEqual(g1["explained_cross_talk_counted_as_gate_failure"], "false")

    def test_writes_tsv_markdown_and_json_with_same_timestamp(self):
        payload = self.module.write_preregistration(self.run_dir, formal_scan_models_sha256="abc123")
        tsv = self.run_dir / "inputs" / "falsification_preregistration.tsv"
        md = self.run_dir / "inputs" / "falsification_preregistration.md"
        js = self.run_dir / "inputs" / "falsification_preregistration.json"
        for path in (tsv, md, js):
            self.assertTrue(path.is_file(), path)
        lines = tsv.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(len(lines), 5, "header + G1..G4")
        self.assertIn("G4", tsv.read_text(encoding="utf-8"))
        self.assertEqual(json.loads(js.read_text(encoding="utf-8"))["registered_at_utc"],
                         payload["registered_at_utc"])
        self.assertIn("E-value < 1e-5", md.read_text(encoding="utf-8"))

    def test_refuses_to_preregister_after_training_started(self):
        (self.run_dir / "results" / "cys_discovery.hmm").write_text("HMMER3/f\n", encoding="utf-8")
        with self.assertRaises(RuntimeError):
            self.module.write_preregistration(self.run_dir, formal_scan_models_sha256="abc123")

    def test_refuses_when_scoring_output_already_exists(self):
        (self.run_dir / "results" / "cys_discovery_scores.tsv").write_text("accession\n", encoding="utf-8")
        with self.assertRaises(RuntimeError):
            self.module.write_preregistration(self.run_dir, formal_scan_models_sha256="abc123")


if __name__ == "__main__":
    unittest.main()
