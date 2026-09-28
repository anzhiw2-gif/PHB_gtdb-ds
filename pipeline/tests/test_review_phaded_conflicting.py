import importlib.util
import unittest

from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "review_phaded_conflicting.py"


def load_module():
    spec = importlib.util.spec_from_file_location("review_phaded_conflicting", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ConflictingReviewTests(unittest.TestCase):
    def test_terminal_stop_is_not_invalid_biology(self):
        module = load_module()
        self.assertEqual(module.sequence_integrity("MAAA*"), "terminal_stop_only")
        self.assertEqual(module.sequence_integrity("MA*AA"), "invalid_internal_character")

    def test_classifies_type2_explicit_and_nphamcl_synthase_separately(self):
        module = load_module()
        explicit = {"phaded_superfamily_best": "extracellular dPHASCL type 2", "interpro_secondary_state": "explicit_pha_related_domain", "pfam_accessions": "PF00756", "assignment_status": "assigned", "family_score_gap": "3.5"}
        synthase = {"phaded_superfamily_best": "intracellular nPHAMCL", "interpro_secondary_state": "pha_synthase_like_alternative", "pfam_accessions": "PF06850", "assignment_status": "ambiguous_family", "family_score_gap": "0"}
        self.assertEqual(module.provisional_label(explicit)[0], "A_retain_high_priority_architecture")
        self.assertEqual(module.provisional_label(synthase)[0], "E_hold_synthase_like_structure")


if __name__ == "__main__":
    unittest.main()
