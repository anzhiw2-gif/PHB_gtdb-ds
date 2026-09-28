import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "reconcile_phaded_conflict_focus.py"


def load_module():
    spec = importlib.util.spec_from_file_location("reconcile_phaded_conflict_focus", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ConflictFocusTests(unittest.TestCase):
    def test_architecture_conflict_has_priority_over_superfamily_ambiguity(self):
        module = load_module()
        rows = module.select_focus_rows([
            {"accession": "A", "domain_evidence_status": "domain_conflict", "profile_evidence_status": "profile_ambiguous_superfamily", "motif_evidence_status": "conflict_relative_position"},
            {"accession": "B", "domain_evidence_status": "domain_supported_partial", "profile_evidence_status": "profile_ambiguous_superfamily", "motif_evidence_status": "partial_catalytic_pattern_panel"},
            {"accession": "C", "domain_evidence_status": "domain_supported_partial", "profile_evidence_status": "profile_trained_hit", "motif_evidence_status": "partial_catalytic_pattern_panel"},
        ])
        by = {row["accession"]: row for row in rows}
        self.assertEqual(by["A"]["focus_decision"], "hold_architecture_conflict")
        self.assertEqual(by["B"]["focus_decision"], "ambiguous_superfamily_unresolved")
        self.assertNotIn("C", by)


if __name__ == "__main__":
    unittest.main()
