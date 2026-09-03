import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class FormalScanContractTests(unittest.TestCase):
    def test_registry_has_unique_models_and_frozen_thresholds(self):
        rows = (ROOT / "config" / "formal_scan_models.tsv").read_text(encoding="utf-8").splitlines()
        header = rows[0].split("\t")
        self.assertEqual(header, ["model", "hmm_source", "threshold", "min_cov", "report_group"])
        models = [line.split("\t") for line in rows[1:] if line.strip()]
        self.assertEqual(len({row[0] for row in models}), 10)
        self.assertIn("ePhaZ_broad_discovery", {row[0] for row in models})
        self.assertEqual(next(row[1] for row in models if row[0] == "iPhaZ"), "frozen_data_root")
        self.assertTrue(all(row[2] == "e-5" for row in models))
        self.assertEqual(next(row[3] for row in models if row[0] == "OH"), "0.6")

    def test_script_requires_dated_run_and_frozen_registry(self):
        script = (ROOT / "scripts" / "formal_frozen_screen.sh").read_text(encoding="utf-8")
        self.assertIn("formal_scan_models.tsv", script)
        self.assertIn("RUN_ID", script)
        self.assertIn("mkdir \"$RUN_ROOT\"", script)
        self.assertIn("hmmsearch", script)
        self.assertIn("expected_count", script)
        self.assertIn("frozen_data_root", script)
        self.assertIn("FROZEN_MODEL_ROOT", script)
        self.assertIn("--preflight-only", script)
        self.assertIn("planned_not_run", script)
        self.assertIn('cp "$SCRIPT_DIR/formal_frozen_screen.sh"', script)
        self.assertIn('"formal_frozen_screen": pathlib.Path', script)
        self.assertIn("filter_hmmsearch_shard.py", script)
        self.assertIn("overlength_exclusions.tsv", script)

    def test_overlength_filter_has_a_targeted_contract(self):
        script = (ROOT / "scripts" / "filter_hmmsearch_shard.py")
        self.assertTrue(script.exists())
        source = script.read_text(encoding="utf-8")
        self.assertIn("100000", source)
        self.assertIn("overlength_exclusions.tsv", source)

    def test_parallel_runner_supports_resume_and_caps_concurrency(self):
        script = (ROOT / "scripts" / "formal_frozen_screen_parallel.sh")
        self.assertTrue(script.exists())
        source = script.read_text(encoding="utf-8")
        self.assertIn("--resume-from", source)
        self.assertIn("parallel -j", source)
        self.assertIn("THREADS -le 60", source)
        self.assertIn("task_status", source)
        self.assertIn("parent_run", source)
        self.assertIn("scan_manifest.json", source)

    def test_monitor_script_is_read_only_and_reports_progress_signals(self):
        script = (ROOT / "scripts" / "monitor_formal_scan.sh")
        self.assertTrue(script.exists())
        source = script.read_text(encoding="utf-8")
        self.assertIn("--run-dir", source)
        self.assertIn("--interval", source)
        self.assertIn("failed_tasks.tsv", source)
        self.assertIn("task_status", source)
        self.assertIn("hmmsearch", source)
        self.assertNotIn("rm -", source)
        self.assertNotIn("kill ", source)


if __name__ == "__main__":
    unittest.main()
