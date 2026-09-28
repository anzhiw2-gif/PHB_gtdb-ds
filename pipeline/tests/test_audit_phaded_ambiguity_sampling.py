import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_phaded_ambiguity_sampling.py"


def load_module():
    spec = importlib.util.spec_from_file_location("audit_phaded_ambiguity_sampling", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AmbiguitySamplingTests(unittest.TestCase):
    def test_gap_bins_and_stratum_are_deterministic(self):
        module = load_module()
        row = {
            "accession": "A", "assignment_status": "ambiguous_family",
            "phaded_superfamily_best": "extracellular dPHASCL type 1",
            "family_score_gap": "0.7", "domain_evidence_status": "domain_conflict",
            "motif_evidence_status": "partial_lipase_box_only",
        }
        self.assertEqual(module.gap_bin(row["family_score_gap"]), "0-1")
        self.assertEqual(module.stratum_key(row), "extracellular dPHASCL type 1|0-1|domain_conflict|partial_lipase_box_only")

    def test_sampling_never_changes_assignment_or_claims_negative(self):
        module = load_module()
        rows = [
            {"accession": f"A{i}", "assignment_status": "ambiguous_family", "phaded_superfamily_best": "sf", "family_score_gap": "0.5", "domain_evidence_status": "domain_pending", "motif_evidence_status": "not_tested_full_library"}
            for i in range(20)
        ]
        sampled = module.sample_rows(rows, fraction=0.1, max_per_stratum=25)
        self.assertEqual(len(sampled), 2)
        self.assertTrue(all(row["assignment_status"] == "ambiguous_family" for row in sampled))
        self.assertTrue(all(row["review_interpretation"] == "candidate_only_unresolved" for row in sampled))


if __name__ == "__main__":
    unittest.main()
