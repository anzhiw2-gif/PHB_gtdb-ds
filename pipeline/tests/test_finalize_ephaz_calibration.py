import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "finalize_ephaz_calibration.py"


def load_module():
    spec = importlib.util.spec_from_file_location("finalize_ephaz_calibration", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FinalizeEphazCalibrationTests(unittest.TestCase):
    def test_writes_audit_report_and_hash_manifest_without_model_promotion(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            inputs = root / "inputs"
            results = root / "results" / "calibration"
            inputs.mkdir(parents=True)
            results.mkdir(parents=True)
            (root / "input_contract.json").write_text(
                json.dumps({"status": "planned_candidate_only", "authorization": {}}),
                encoding="utf-8",
            )
            panel = inputs / "subtype_training_panels.tsv"
            fields = [
                "subtype", "accession", "role", "substrate_class", "localization",
                "phaded_superfamily", "genus", "sequence_sha256", "evidence_status",
                "threshold", "coverage_threshold", "heldout_hit", "negative_hit",
            ]
            row = {
                "subtype": "e_dPHASCL_type1", "accession": "A1", "role": "train_positive",
                "substrate_class": "SCL", "localization": "extracellular",
                "phaded_superfamily": "extracellular dPHASCL type 1", "genus": "GenusA",
                "sequence_sha256": hashlib.sha256(b"A1").hexdigest(),
                "evidence_status": "experimental_positive", "threshold": "1e-5",
                "coverage_threshold": "0.0", "heldout_hit": "false", "negative_hit": "false",
            }
            with panel.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(row)
            (results / "leave_one_out_metrics.tsv").write_text(
                "subtype\ttrain_positive\theldout_positive\tgenera\theldout_recovered\tnegative_hits\tthreshold\tcoverage_threshold\tdecision\n"
                "e_dPHASCL_type1\t1\t0\t1\tfalse\t0\t1e-5\t0.0\treference_only_insufficient_panel\n",
                encoding="utf-8",
            )
            decision = {
                "status": "candidate-only",
                "decisions": {"e_dPHASCL_type1": "reference_only_insufficient_panel"},
                "formal_registry_modified": False,
                "historical_ephaz_model_overwritten": False,
            }
            (results / "subtype_model_decision.json").write_text(
                json.dumps(decision), encoding="utf-8"
            )
            output = module.finalize(root)
            self.assertEqual(output["status"], "completed_candidate_only")
            self.assertFalse(output["model_promotion_authorized"])
            self.assertTrue((results / "ephaz_seed_audit.tsv").is_file())
            self.assertTrue((results / "calibration_status.md").is_file())
            manifest = json.loads((results / "calibration_sha256_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("inputs/subtype_training_panels.tsv", manifest["files"])
            self.assertIn("results/calibration/ephaz_seed_audit.tsv", manifest["files"])


if __name__ == "__main__":
    unittest.main()
