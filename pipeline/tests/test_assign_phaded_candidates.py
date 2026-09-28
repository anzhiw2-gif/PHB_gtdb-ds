import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "assign_phaded_candidates.py"


def load_module():
    spec = importlib.util.spec_from_file_location("assign_phaded_candidates", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PhaDEDAssignmentTests(unittest.TestCase):
    def test_tied_superfamily_scores_remain_ambiguous(self):
        module = load_module()
        result = module.assign_candidate(
            accession="GCF_000001|gene_1",
            scores={"e_dPHASCL_type1": 102.0, "e_dPHASCL_type2": 102.0},
            profile_status={"e_dPHASCL_type1": "calibrated", "e_dPHASCL_type2": "calibrated"},
        )
        self.assertEqual(result["assignment_status"], "ambiguous_superfamily")
        self.assertEqual(result["phaded_superfamily_best"], "")

    def test_reference_only_family_is_not_presented_as_calibrated(self):
        module = load_module()
        result = module.assign_candidate(
            "GCF_000001|gene_2", {"family_17": 88.0}, {"family_17": "reference_only"}
        )
        self.assertEqual(result["assignment_status"], "reference_only_nearest")
        self.assertEqual(result["phaded_family_best"], "family_17")

    def test_small_gap_is_ambiguous_and_retains_second_profile(self):
        module = load_module()
        result = module.assign_candidate(
            "candidate", {"sf_a": 100.0, "sf_b": 99.5},
            {"sf_a": "calibrated", "sf_b": "calibrated"}, min_score_gap=1.0,
        )
        self.assertEqual(result["assignment_status"], "ambiguous_superfamily")
        self.assertEqual(result["phaded_superfamily_second"], "sf_b")
        self.assertAlmostEqual(result["superfamily_score_gap"], 0.5)

    def test_missing_scores_are_not_forced(self):
        module = load_module()
        result = module.assign_candidate("candidate", {}, {})
        self.assertEqual(result["assignment_status"], "unassigned_PhaDED_like")
        self.assertEqual(result["evidence_status"], "no_profile_score")

    def test_family_ambiguity_is_retained_inside_assigned_superfamily(self):
        module = load_module()
        metadata = {
            "family_a": {"superfamily": "sf", "family": "family_a"},
            "family_b": {"superfamily": "sf", "family": "family_b"},
        }
        result = module.assign_candidate(
            "candidate", {"family_a": 100.0, "family_b": 100.0},
            {"family_a": "calibrated", "family_b": "calibrated"},
            profile_metadata=metadata,
        )
        self.assertEqual(result["assignment_status"], "ambiguous_family")
        self.assertEqual(result["phaded_superfamily_best"], "sf")

    def test_reference_only_loser_does_not_downgrade_calibrated_winner(self):
        module = load_module()
        metadata = {
            "calibrated": {"superfamily": "sf_calibrated", "family": "family_calibrated"},
            "nearest": {"superfamily": "sf_reference", "family": "family_reference"},
        }
        result = module.assign_candidate(
            "candidate", {"calibrated": 100.0, "nearest": 80.0},
            {"calibrated": "calibrated", "nearest": "reference_only"},
            profile_metadata=metadata,
        )
        self.assertEqual(result["assignment_status"], "assigned")
        self.assertEqual(result["phaded_family_best"], "calibrated")

    def test_rejects_duplicate_profile_score_rows(self):
        module = load_module()
        candidates = [{"accession": "candidate", "genome": "G1"}]
        scores = [
            {"accession": "candidate", "profile_id": "profile", "score": "10"},
            {"accession": "candidate", "profile_id": "profile", "score": "9"},
        ]
        with self.assertRaisesRegex(ValueError, "duplicate profile score"):
            module.map_candidates(candidates, scores, {"profile": "calibrated"}, {"profile": {"superfamily": "sf"}})


if __name__ == "__main__":
    unittest.main()
