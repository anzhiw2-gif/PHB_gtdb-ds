import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "reclassify_ephaz_controls.py"


class ReclassifyControlTests(unittest.TestCase):
    def _load(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT}")
        spec = importlib.util.spec_from_file_location("reclassify_ephaz_controls", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module

    def test_six_controls_are_challenges_and_excluded_from_ephaz_panel(self):
        module = self._load()
        six = {
            "A0ABY7N197", "A0A3Q9BTL6", "A0ABY9RN41",
            "A0A3N5XWX8", "A0ABT5L4W6", "A0ABT7T1G3",
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            controls = root / "controls.tsv"
            positive = root / "positive.faa"
            classes = root / "classes.tsv"
            outdir = root / "out"
            rows = [
                ("B2NHN2", "positive", "e-PhaZ_EC", "ePhaZ_curated_core"),
                *[(acc, "positive", "e-PhaZ_remote", "iPhaZ_like_challenge") for acc in sorted(six)],
            ]
            with controls.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["accession", "label", "query_group"])
                writer.writerows(row[:3] for row in rows)
            with classes.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["accession", "control_class", "rationale"])
                writer.writerows((row[0], row[3], "test") for row in rows)
            with positive.open("w", encoding="utf-8") as handle:
                for accession, _, _, _ in rows:
                    handle.write(f">{accession}|test\n{'M' * 100}\n")

            result = module.reclassify_controls(controls, positive, classes, outdir)
            self.assertEqual(result["counts"]["iPhaZ_like_challenge"], 6)
            self.assertEqual(result["counts"]["ePhaZ_curated_core"], 1)
            challenge = {
                line[1:].split("|", 1)[0]
                for line in (outdir / "iPhaZ_like_challenge.faa").read_text(encoding="utf-8").splitlines()
                if line.startswith(">")
            }
            ephaz = {
                line[1:].split("|", 1)[0]
                for line in (outdir / "ephaz_positive_controls.faa").read_text(encoding="utf-8").splitlines()
                if line.startswith(">")
            }
            self.assertEqual(challenge, six)
            self.assertEqual(ephaz, {"B2NHN2"})
            self.assertTrue((outdir / "control_governance.json").is_file())

    def test_missing_class_accession_fails_closed(self):
        module = self._load()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            controls = root / "controls.tsv"
            positive = root / "positive.faa"
            classes = root / "classes.tsv"
            controls.write_text("accession\tlabel\tquery_group\nX\tpositive\te-PhaZ\n", encoding="utf-8")
            positive.write_text(">X|test\nMMMM\n", encoding="utf-8")
            classes.write_text("accession\tcontrol_class\trationale\nY\tePhaZ_curated_core\ttest\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                module.reclassify_controls(controls, positive, classes, root / "out")


if __name__ == "__main__":
    unittest.main()
