import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "finalize_phaded_evidence_completion.py"


def load_module():
    spec = importlib.util.spec_from_file_location("finalize_phaded_evidence_completion", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class EvidenceCompletionTests(unittest.TestCase):
    def test_summary_preserves_candidate_boundary_and_layer_counts(self):
        module = load_module()
        rows = [
            {"accession": "A", "subtype_call": "dPHASCL1_like_candidate", "subtype_confidence": "moderate_candidate_only", "evidence_grade": "L1", "motif_evidence_status": "partial_catalytic_pattern_panel", "domain_evidence_status": "domain_supported_partial", "profile_evidence_status": "profile_trained_hit", "interpro_status": "interpro_supported"},
            {"accession": "B", "subtype_call": "hold_architecture_conflict", "subtype_confidence": "hold", "evidence_grade": "HOLD", "motif_evidence_status": "conflict_relative_position", "domain_evidence_status": "domain_conflict", "profile_evidence_status": "profile_ambiguous_superfamily", "interpro_status": "interpro_supported"},
        ]
        summary = module.summarize_rows(rows)
        self.assertEqual(summary["candidate_count"], 2)
        self.assertEqual(summary["subtype_call_counts"]["hold_architecture_conflict"], 1)
        self.assertEqual(summary["evidence_layers"]["motif_evidence_status"]["conflict_relative_position"], 1)
        self.assertIn("candidate-only", summary["phenotype_boundary"])


if __name__ == "__main__":
    unittest.main()
