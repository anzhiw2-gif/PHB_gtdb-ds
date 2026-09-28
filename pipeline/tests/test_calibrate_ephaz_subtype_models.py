import csv
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "calibrate_ephaz_subtype_models.py"


def load_module():
    spec = importlib.util.spec_from_file_location("calibrate_ephaz_subtype_models", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FIELDS = [
    "subtype", "accession", "role", "substrate_class", "localization", "phaded_superfamily",
    "genus", "sequence_sha256", "evidence_status", "threshold", "coverage_threshold",
    "heldout_hit", "negative_hit",
]


def write_manifest(root, rows):
    path = root / "panels.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = {field: "" for field in FIELDS}
            output.update(row)
            writer.writerow(output)
    return path


def positive(accession, subtype="e_dPHASCL_type1", genus="GenusA", role="train_positive", substrate="SCL"):
    return {
        "subtype": subtype, "accession": accession, "role": role, "substrate_class": substrate,
        "localization": "extracellular", "phaded_superfamily": "extracellular dPHASCL type 1",
        "genus": genus, "sequence_sha256": hashlib.sha256(accession.encode("ascii")).hexdigest(), "evidence_status": "experimental_positive",
        "threshold": "1e-5", "coverage_threshold": "0.6", "heldout_hit": "true", "negative_hit": "false",
    }


class EphaZSubtypeCalibrationTests(unittest.TestCase):
    def test_rejects_q51718_from_phb_scl_training_panel(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = [positive("B2NHN2"), positive("Q51718", genus="GenusQ", substrate="MCL")]
            with self.assertRaisesRegex(ValueError, "substrate-class mismatch"):
                module.calibrate_subtype_models(write_manifest(root, rows), root / "out")

    def test_rejects_heldout_accession_present_in_training(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = [positive("B2NHN2"), positive("B2NHN2", role="heldout_positive")]
            with self.assertRaisesRegex(ValueError, "train-heldout leakage"):
                module.calibrate_subtype_models(write_manifest(root, rows), root / "out")

    def test_insufficient_panel_is_explicitly_reference_only(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            result = module.calibrate_subtype_models(
                write_manifest(root, [positive("B2NHN2")]), root / "out"
            )
            self.assertEqual(result["decisions"]["e_dPHASCL_type1"], "reference_only_insufficient_panel")

    def test_rejects_threshold_changes_within_subtype(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = positive("B2NHN2")
            second = positive("O05527", genus="GenusB")
            second["threshold"] = "1e-10"
            with self.assertRaisesRegex(ValueError, "threshold change"):
                module.calibrate_subtype_models(write_manifest(root, [first, second]), root / "out")

    def test_rejects_nonexperimental_training_evidence(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            rows = [positive(f"P{i}", genus=f"Genus{i}") for i in range(1, 4)]
            rows[0]["evidence_status"] = "annotation_only"
            with self.assertRaisesRegex(ValueError, "evidence_status"):
                module.calibrate_subtype_models(write_manifest(root, rows), root / "out")


if __name__ == "__main__":
    unittest.main()
