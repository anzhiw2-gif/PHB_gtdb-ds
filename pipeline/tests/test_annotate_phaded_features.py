import importlib.util
import unittest
from pathlib import Path
import hashlib


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_phaded_features.py"


def load_module():
    spec = importlib.util.spec_from_file_location("annotate_phaded_features", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PhaDEDFeatureTests(unittest.TestCase):
    def test_phaz7_like_ahsmg_is_not_failed_for_lacking_gx_sx_g(self):
        module = load_module()
        features = module.annotate_sequence(
            "candidate_1", "M" + "A" * 80 + "AHSMG" + "A" * 180,
            "extracellular native-SCL/PhaZ7-like",
        )
        self.assertEqual(features["ahsmg_state"], "present")
        self.assertEqual(features["lipase_box_state"], "not_expected_for_superfamily")

    def test_unavailable_domain_database_stays_pending(self):
        module = load_module()
        features = module.annotate_sequence(
            "candidate_2", "M" + "A" * 240,
            "intracellular nPHAMCL", domain_records=None,
        )
        self.assertEqual(features["pfam_interpro_state"], "pending")
        self.assertEqual(features["lid_state"], "pending")

    def test_coordinates_are_one_based_and_signalp_is_imported(self):
        module = load_module()
        sequence = "M" + "A" * 9 + "GASAG" + "A" * 20
        features = module.annotate_sequence(
            "candidate_3", sequence, "extracellular dPHASCL type 1",
            signalp_record={"type": "SP", "status": "accepted"},
        )
        self.assertEqual(features["lipase_box_state"], "present")
        self.assertEqual(features["lipase_box_coordinates"], "11-15")
        self.assertEqual(features["signalp_class"], "SP")
        self.assertEqual(features["signalp_evidence"], "accepted")

    def test_terminal_stop_is_recorded_without_calling_the_sequence_invalid(self):
        module = load_module()
        features = module.annotate_sequence(
            "candidate_4", "M" + "A" * 30 + "*", "extracellular dPHASCL type 1",
        )
        self.assertEqual(features["sequence_integrity"], "terminal_stop_only")
        self.assertEqual(features["signalp_evidence"], "no_record")
        self.assertNotEqual(features["signalp_class"], "negative")

    def test_terminal_stop_with_missing_initiator_remains_possible_n_truncation(self):
        module = load_module()
        features = module.annotate_sequence(
            "candidate_n_truncated", "AAA*", "extracellular dPHASCL type 1",
        )
        self.assertEqual(features["sequence_integrity"], "possible_N_truncation")

    def test_internal_stop_or_other_non_iupac_residue_is_invalid_internal_character(self):
        module = load_module()
        internal_stop = module.annotate_sequence(
            "candidate_internal_stop", "MAA*AA", "extracellular dPHASCL type 1",
        )
        non_iupac = module.annotate_sequence(
            "candidate_non_iupac", "MAAXAA", "extracellular dPHASCL type 1",
        )
        self.assertEqual(internal_stop["sequence_integrity"], "invalid_internal_character")
        self.assertEqual(non_iupac["sequence_integrity"], "invalid_internal_character")

    def test_feature_merge_preserves_assignment_and_flags_conflict_for_review(self):
        module = load_module()
        assignments = [{"accession": "candidate_5", "assignment_status": "assigned"}]
        features = [{"accession": "candidate_5", "architecture_consistency": "conflicting", "feature_evidence_status": "pending"}]
        merged = module.merge_assignment_evidence(assignments, features)
        self.assertEqual(merged[0]["assignment_status"], "assigned")
        self.assertEqual(merged[0]["assignment_review"], "review_required")

    def test_sequence_fasta_binding_rejects_extra_record_and_hash_mismatch(self):
        module = load_module()
        sequence = "M" + "A" * 20
        rows = [{"accession": "A", "phaded_superfamily": "extracellular dPHASCL type 1", "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest()}]
        with self.assertRaisesRegex(ValueError, "FASTA accession set mismatch"):
            module.annotate_rows(rows, {"A": sequence, "EXTRA": sequence})
        rows[0]["sequence_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "sequence SHA-256 mismatch"):
            module.annotate_rows(rows, {"A": sequence})


if __name__ == "__main__":
    unittest.main()
