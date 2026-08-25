"""Contract tests for isolated dated runs and GTDB input provenance."""
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_context.py"
SPEC = importlib.util.spec_from_file_location("run_context", SCRIPT)
run_context = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run_context)


class RunContextTests(unittest.TestCase):
    def test_validate_run_id_accepts_safe_dated_identifier(self):
        self.assertEqual(
            run_context.validate_run_id("20260824T120000Z_9a7d02d"),
            "20260824T120000Z_9a7d02d",
        )

    def test_validate_run_id_rejects_path_traversal_and_empty_values(self):
        for run_id in ("", "..", "../escape", "a/b", "a\\b", "run id"):
            with self.subTest(run_id=run_id):
                with self.assertRaises(ValueError):
                    run_context.validate_run_id(run_id)

    def test_create_run_layout_makes_isolated_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = run_context.create_run_layout(tmp, "20260824T120000Z_test")
            self.assertEqual(run_dir, Path(tmp) / "runs" / "20260824T120000Z_test")
            for name in ("logs", "inputs", "results"):
                self.assertTrue((run_dir / name).is_dir(), name)
            with self.assertRaises(FileExistsError):
                run_context.create_run_layout(tmp, "20260824T120000Z_test")

    def test_input_contract_can_bind_explicit_legacy_run_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "project"
            run_dir.mkdir()
            contract = run_context.build_input_contract(run_dir, run_id="legacy")
            self.assertEqual(contract["run_id"], "legacy")

    def test_sha256_file_matches_standard_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.tsv"
            path.write_bytes(b"taxonomy\n")
            expected = hashlib.sha256(b"taxonomy\n").hexdigest()
            self.assertEqual(run_context.sha256_file(path), expected)

    def test_input_contract_hashes_verified_inputs_and_marks_missing_gtdb_pending(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = run_context.create_run_layout(tmp, "20260824T120000Z_test")
            taxonomy = Path(tmp) / "taxonomy.tsv"
            taxonomy.write_text("GCA_1\td__Bacteria\n", encoding="utf-8")
            contract = run_context.build_input_contract(
                run_dir,
                gtdb_inputs={
                    "taxonomy": taxonomy,
                    "metadata": Path(tmp) / "missing_metadata.tsv",
                    "tree": None,
                },
                inputs={"parameters": taxonomy},
            )
            self.assertEqual(contract["status"], "pending")
            self.assertEqual(contract["gtdb"]["taxonomy"]["status"], "verified")
            self.assertEqual(
                contract["gtdb"]["taxonomy"]["sha256"],
                run_context.sha256_file(taxonomy),
            )
            self.assertEqual(contract["gtdb"]["metadata"]["status"], "pending")
            self.assertIsNone(contract["gtdb"]["metadata"]["sha256"])
            self.assertEqual(contract["gtdb"]["tree"]["status"], "pending")
            self.assertEqual(contract["inputs"]["parameters"]["status"], "verified")

    def test_non_gtdb_missing_input_is_missing_and_contract_writes_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = run_context.create_run_layout(tmp, "20260824T120000Z_test")
            contract = run_context.write_input_contract(
                run_dir,
                inputs={"hmm": Path(tmp) / "missing.hmm"},
            )
            output = run_dir / "input_contract.json"
            self.assertTrue(output.is_file())
            loaded = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(loaded, contract)
            self.assertEqual(contract["status"], "pending")
            self.assertEqual(contract["inputs"]["hmm"]["status"], "missing")


if __name__ == "__main__":
    unittest.main()
