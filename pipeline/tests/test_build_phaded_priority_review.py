import importlib.util
import unittest
import csv
import tempfile

from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_phaded_priority_review.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_phaded_priority_review", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PriorityReviewTests(unittest.TestCase):
    def test_main_writes_manifest_for_interpro_output(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "input.tsv"
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["accession", "layer", "architecture_consistency"], delimiter="\t")
                writer.writeheader()
                writer.writerow({"accession": "A", "layer": "ePhaZ_competition_review", "architecture_consistency": "partial"})
            self.assertEqual(module.main(["--input", str(source), "--output-dir", str(root / "out")]), 0)
            self.assertTrue((root / "out" / "priority_review_manifest.json").is_file())

    def test_queue_order_deduplicates_with_conflict_precedence(self):
        module = load_module()
        rows = [
            {"accession": "C", "layer": "iPhaZ_tier1", "architecture_consistency": "partial"},
            {"accession": "A", "layer": "ePhaZ_curated_core_secreted", "architecture_consistency": "partial"},
            {"accession": "B", "layer": "ePhaZ_competition_review", "architecture_consistency": "partial"},
            {"accession": "A", "layer": "ePhaZ_curated_core_secreted", "architecture_consistency": "conflicting"},
        ]
        queued = module.build_review_queue(rows)
        self.assertEqual([r["accession"] for r in queued], ["A", "B", "C"])
        self.assertEqual(queued[0]["review_stage"], "conflicting")
        self.assertEqual(queued[1]["review_stage"], "ePhaZ_competition_review")

    def test_interpro_forces_conflicts_and_competition_then_caps_strata(self):
        module = load_module()
        rows = [
            {"accession": "X", "layer": "ePhaZ_curated_core_secreted", "architecture_consistency": "partial", "phaded_superfamily_best": "sf", "assignment_status": "assigned", "priority_score": "5"},
            {"accession": "Y", "layer": "ePhaZ_competition_review", "architecture_consistency": "partial", "phaded_superfamily_best": "sf", "assignment_status": "assigned", "priority_score": "1"},
            {"accession": "Z", "layer": "iPhaZ_tier1", "architecture_consistency": "conflicting", "phaded_superfamily_best": "sf2", "assignment_status": "assigned", "priority_score": "1"},
        ]
        selected = module.select_interpro_candidates(rows, cap_per_stratum=1)
        self.assertEqual([r["accession"] for r in selected], ["Z", "Y", "X"])
        self.assertTrue(all(r["interpro_selection_reason"] for r in selected))


if __name__ == "__main__":
    unittest.main()
