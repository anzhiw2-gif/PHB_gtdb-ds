import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "derive_phaded_pfam_architecture.py"


def load_module():
    spec = importlib.util.spec_from_file_location("derive_phaded_pfam_architecture", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PhaDEDPfamArchitectureTests(unittest.TestCase):
    def test_parses_hmmscan_domtblout_with_query_accession_and_coordinates(self):
        module = load_module()
        row = "Lipase_3 PF01764.21 278 candidate_1 - 350 1e-20 80.0 0.0 1 1 1e-20 1e-20 80.0 0.0 1 278 15 300 12 304 0.98 lipase"
        result = module.parse_domtblout([row])
        self.assertEqual(result[0]["accession"], "candidate_1")
        self.assertEqual(result[0]["pfam_accession"], "PF01764")
        self.assertEqual(result[0]["coordinates"], "15-300")

    def test_only_matching_reference_fingerprint_gives_partial_support(self):
        module = load_module()
        expectation = {"sf_a": {"PF00001"}, "sf_b": {"PF00002"}}
        result = module.classify_candidate_domains("sf_a", {"PF00001"}, expectation)
        self.assertEqual(result["architecture_consistency"], "partial")
        self.assertEqual(result["pfam_interpro_state"], "tested_targeted_pfam_hit")

    def test_reference_fingerprint_normalizes_reference_id_prefix(self):
        module = load_module()
        hits = [{"accession": "DED_hfam_52_0001|P12625.1", "pfam_accession": "PF06850"}]
        result = module.reference_fingerprint(hits, {"DED_hfam_52_0001": "sf_a"}, min_reference_count=1)
        self.assertEqual(result, {"sf_a": {"PF06850"}})

    def test_selected_name_rows_preserve_pfam_accession_name_binding(self):
        module = load_module()
        hits = [{"pfam_accession": "PF06850", "pfam_name": "PHB_depo_C"}]
        self.assertEqual(module.selected_name_rows(hits, {"PF06850"}), [("PF06850", "PHB_depo_C")])

    def test_exclusive_other_superfamily_fingerprint_is_review_conflict(self):
        module = load_module()
        expectation = {"sf_a": {"PF00001"}, "sf_b": {"PF00002"}}
        result = module.classify_candidate_domains("sf_a", {"PF00002"}, expectation)
        self.assertEqual(result["architecture_consistency"], "conflicting")
        self.assertEqual(result["assignment_review"], "review_required")


if __name__ == "__main__":
    unittest.main()
