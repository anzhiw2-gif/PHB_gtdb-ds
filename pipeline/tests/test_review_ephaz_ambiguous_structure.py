import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "review_ephaz_ambiguous_structure.py"


class StructuralReviewTests(unittest.TestCase):
    def _load(self):
        spec = importlib.util.spec_from_file_location("review_ephaz_ambiguous_structure", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module

    def test_merge_intervals_and_coverage(self):
        module = self._load()
        self.assertEqual(module.merge_intervals([(5, 10), (8, 20), (30, 35)]), [(5, 20), (30, 35)])
        self.assertAlmostEqual(module.coverage_from_intervals([(5, 10), (8, 20), (30, 35)], 100), 0.22)

    def test_classification_is_conservative_and_integrity_gated(self):
        module = self._load()
        ephaz = [(10, 130), (185, 280)]
        iphaz = [(175, 330)]
        self.assertEqual(module.classify_architecture(ephaz, iphaz, 0.45, 0.35), "mixed_cross_family")
        self.assertEqual(module.classify_architecture([], [(175, 330)], 0.0, 0.34), "iPhaZ_consistent")
        self.assertEqual(module.classify_architecture([(10, 130)], [], 0.38, 0.0), "partial_ePhaZ_signal")
        self.assertEqual(module.review_decision("complete", "iPhaZ_consistent"), "provisional_iPhaZ_challenge")
        self.assertEqual(module.review_decision("possible_N_truncation", "iPhaZ_consistent"), "pending_manual")


if __name__ == "__main__":
    unittest.main()
