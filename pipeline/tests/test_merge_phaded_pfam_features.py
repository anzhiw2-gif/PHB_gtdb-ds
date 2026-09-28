import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "merge_phaded_pfam_features.py"


def load_module():
    spec = importlib.util.spec_from_file_location("merge_phaded_pfam_features", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MergePhaDEDPfamFeatureTests(unittest.TestCase):
    def test_pfam_architecture_supersedes_pending_baseline_without_losing_motif_fields(self):
        module = load_module()
        feature = [{"accession": "candidate_1", "lipase_box_state": "present", "architecture_consistency": "pending"}]
        pfam = [{"accession": "candidate_1", "pfam_accessions": "PF00561", "architecture_consistency": "partial", "pfam_interpro_state": "tested_targeted_pfam_hit", "assignment_review": "no_conflict_recorded"}]
        result = module.merge(feature, pfam)
        self.assertEqual(result[0]["lipase_box_state"], "present")
        self.assertEqual(result[0]["pfam_accessions"], "PF00561")
        self.assertEqual(result[0]["architecture_consistency"], "partial")

    def test_merged_header_does_not_duplicate_existing_evidence_fields(self):
        module = load_module()
        fields = module.merged_fields(["accession", "architecture_consistency", "pfam_interpro_state"])
        self.assertEqual(len(fields), len(set(fields)))
