import csv
import importlib.util
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "calibrate_ephaz_control_classes.py"


class CalibrationClassTests(unittest.TestCase):
    def _load(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT}")
        spec = importlib.util.spec_from_file_location("calibrate_ephaz_control_classes", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module

    def test_summarize_metrics_keeps_control_classes_separate(self):
        module = self._load()
        classes = {
            "CORE": "ePhaZ_curated_core",
            "REMOTE": "ePhaZ_architecture_remote",
            "CHALLENGE": "iPhaZ_like_challenge",
            "NEG": "negative",
        }
        hits = {
            "CORE": {"evalue": 1e-40, "bitscore": 100.0, "coverage": 0.9},
            "REMOTE": {"evalue": 1e-8, "bitscore": 70.0, "coverage": 0.8},
            "NEG": {"evalue": 1e-30, "bitscore": 120.0, "coverage": 0.9},
        }
        rows = module.summarize_by_class(classes, hits, threshold=1e-10, min_cov=0.5)
        by_class = {row["control_class"]: row for row in rows}
        self.assertEqual((by_class["ePhaZ_curated_core"]["TP"], by_class["ePhaZ_curated_core"]["FN"]), (1, 0))
        self.assertEqual((by_class["ePhaZ_architecture_remote"]["TP"], by_class["ePhaZ_architecture_remote"]["FN"]), (0, 1))
        self.assertEqual(by_class["iPhaZ_like_challenge"]["tested"], 1)
        self.assertEqual(by_class["iPhaZ_like_challenge"]["detected"], 0)
        self.assertEqual(by_class["negative"]["FP"], 1)

    def test_unknown_control_class_is_rejected(self):
        module = self._load()
        with self.assertRaises(ValueError):
            module.validate_control_classes({"X": "not_a_supported_class"})

    def test_challenge_fasta_is_separate_from_ephaz_positive_fasta(self):
        module = self._load()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            positive = root / "positive.faa"
            challenge = root / "challenge.faa"
            negative = root / "negative.faa"
            classes = root / "classes.tsv"
            for path, records in (
                (positive, [("CORE", "M" * 100)]),
                (challenge, [("CHALLENGE", "M" * 90)]),
                (negative, [("NEG", "M" * 80)]),
            ):
                with path.open("w", encoding="utf-8") as handle:
                    for accession, sequence in records:
                        handle.write(f">{accession}|test\n{sequence}\n")
            for name in ("ephaz.hmm", "iphaz.hmm"):
                (root / name).write_text("HMMER3/f [test]\n", encoding="utf-8")
            classes.write_text(
                "accession\tcontrol_class\n"
                "CORE\tePhaZ_curated_core\n"
                "CHALLENGE\tiPhaZ_like_challenge\n",
                encoding="utf-8",
            )
            def fake_hmmsearch(command, check, capture_output, text):
                (Path(command[command.index("--tblout") + 1])).write_text("", encoding="utf-8")
                (Path(command[command.index("--domtblout") + 1])).write_text("", encoding="utf-8")

            with patch.object(module.subprocess, "run", side_effect=fake_hmmsearch):
                module.calibrate_control_classes(
                    root / "ephaz.hmm", root / "iphaz.hmm", positive, negative,
                    classes, root / "out", challenge_faa=challenge,
                    hmmsearch_bin="hmmsearch",
                )
            with (root / "out" / "calibration_by_class.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_class = {(row["model"], row["control_class"]): row for row in rows}
            self.assertEqual(by_class[("ePhaZ", "ePhaZ_curated_core")]["tested"], "1")
            self.assertEqual(by_class[("ePhaZ", "iPhaZ_like_challenge")]["tested"], "1")
            self.assertEqual(by_class[("ePhaZ", "iPhaZ_like_challenge")]["challenge_detected"], "0")


if __name__ == "__main__":
    unittest.main()
