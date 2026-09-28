import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_phaded_evidence_tiers.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_phaded_evidence_tiers", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class EvidenceTierTests(unittest.TestCase):
    def test_terminal_stop_overrides_stale_invalid_character_label(self):
        module = load_module()
        self.assertEqual(
            module.normalized_sequence_integrity("invalid_character", "MSTAG*"),
            "terminal_stop_only",
        )

    def test_internal_stop_remains_a_sequence_integrity_hold(self):
        module = load_module()
        self.assertEqual(
            module.normalized_sequence_integrity("valid", "MST*AG"),
            "invalid_internal_character",
        )

    def test_conflict_has_hold_grade(self):
        module = load_module()
        row = {"assignment_status": "assigned", "architecture_consistency": "conflicting", "sequence_integrity": "valid"}
        result = module.grade_row(row)
        self.assertEqual(result["evidence_grade"], "HOLD_architecture_conflict")

    def test_invalid_sequence_is_not_called_negative(self):
        module = load_module()
        row = {"assignment_status": "assigned", "architecture_consistency": "partial", "sequence_integrity": "invalid_internal_character"}
        result = module.grade_row(row)
        self.assertEqual(result["evidence_grade"], "HOLD_sequence_integrity")
        self.assertNotIn("negative", result["evidence_grade"].lower())

    def test_partial_architecture_is_retained_with_explicit_grade(self):
        module = load_module()
        row = {"assignment_status": "assigned", "architecture_consistency": "partial", "sequence_integrity": "valid"}
        result = module.grade_row(row)
        self.assertEqual(result["evidence_grade"], "L1_profile_plus_partial_architecture")

    def test_unassigned_remains_pending(self):
        module = load_module()
        row = {"assignment_status": "unassigned_PhaDED_like", "architecture_consistency": "pending", "sequence_integrity": "valid"}
        result = module.grade_row(row)
        self.assertEqual(result["evidence_grade"], "PENDING_profile_unassigned")

    def test_unrun_evidence_layers_are_explicitly_not_tested(self):
        module = load_module()
        result = module.grade_row({"assignment_status": "assigned", "architecture_consistency": "pending", "sequence_integrity": "valid"})
        self.assertEqual(result["structure_evidence_status"], "not_tested_full_library")
        self.assertEqual(result["phylogeny_evidence_status"], "not_tested_full_library")
