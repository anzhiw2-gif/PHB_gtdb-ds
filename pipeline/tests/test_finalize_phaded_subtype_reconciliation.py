import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "finalize_phaded_subtype_reconciliation.py"


def load_module():
    spec = importlib.util.spec_from_file_location("finalize_phaded_subtype_reconciliation", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FinalizeSubtypeReconciliationTests(unittest.TestCase):
    def test_summarize_rows_reports_required_conflicts_and_evidence_layers(self):
        module = load_module()
        rows = [
            {
                "accession": "A",
                "assignment_status": "assigned",
                "profile_evidence_status": "profile_trained_hit",
                "domain_evidence_status": "domain_conflict",
                "motif_evidence_status": "partial_lipase_box_only",
                "localization_evidence_status": "export_signal_supported",
                "structure_evidence_status_v2": "structure_conflict_candidate_only",
                "phylogeny_evidence_status_v2": "phylogeny_support_candidate_only",
                "interpro_status": "interpro_supported",
                "subtype_call": "hold_architecture_conflict",
                "subtype_confidence": "hold",
                "manual_review_status": "hold_architecture_conflict",
            },
            {
                "accession": "B",
                "assignment_status": "ambiguous_superfamily",
                "profile_evidence_status": "profile_ambiguous_superfamily",
                "domain_evidence_status": "domain_pending",
                "motif_evidence_status": "not_tested_full_library",
                "localization_evidence_status": "not_tested_full_library",
                "structure_evidence_status_v2": "not_tested_full_library",
                "phylogeny_evidence_status_v2": "not_tested_full_library",
                "interpro_status": "interpro_not_reported",
                "subtype_call": "ambiguous_superfamily",
                "subtype_confidence": "unresolved",
                "manual_review_status": "not_in_19_candidate_review",
            },
            {
                "accession": "C",
                "assignment_status": "ambiguous_family",
                "profile_evidence_status": "profile_ambiguous_family",
                "domain_evidence_status": "domain_supported_generic_interpro",
                "motif_evidence_status": "partial_lipase_box_only",
                "localization_evidence_status": "localization_unassigned",
                "structure_evidence_status_v2": "not_tested_full_library",
                "phylogeny_evidence_status_v2": "not_tested_full_library",
                "interpro_status": "interpro_tool_input_excluded",
                "subtype_call": "ambiguous_family_within_superfamily",
                "subtype_confidence": "unresolved",
                "manual_review_status": "not_in_19_candidate_review",
            },
        ]
        summary = module.summarize_rows(rows)
        self.assertEqual(summary["candidate_count"], 3)
        self.assertEqual(summary["architecture_conflict"], 1)
        self.assertEqual(summary["superfamily_ambiguity"], 1)
        self.assertEqual(summary["family_ambiguity"], 1)
        self.assertEqual(summary["profile_evidence_status_counts"]["profile_ambiguous_family"], 1)
        self.assertEqual(summary["interpro_status_counts"]["interpro_tool_input_excluded"], 1)
        self.assertEqual(summary["evidence_layer_counts"]["domain_conflict"], 1)

    def test_write_reports_creates_machine_readable_conflict_and_focus_outputs(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            matrix = output_dir / "matrix.tsv"
            profile_gap = output_dir / "profile_gap.tsv"
            manual = output_dir / "manual.tsv"
            matrix.write_text(
                "\t".join([
                    "accession", "assignment_status", "profile_evidence_status", "domain_evidence_status",
                    "motif_evidence_status", "localization_evidence_status", "structure_evidence_status_v2",
                    "phylogeny_evidence_status_v2", "interpro_status", "subtype_call", "subtype_confidence",
                    "manual_review_status", "evidence_grade", "literature_evidence_tier",
                ]) + "\n" + "\t".join([
                    "A", "assigned", "profile_trained_hit", "domain_supported_partial", "partial_lipase_box_only",
                    "export_signal_supported", "not_tested_full_library", "not_tested_full_library",
                    "interpro_supported", "dPHASCL1_like_candidate", "moderate_candidate_only",
                    "retain_candidate_only", "L1_profile_plus_partial_architecture", "none_exact_candidate_match",
                ]) + "\n",
                encoding="utf-8",
            )
            profile_gap.write_text("profile_id\tcoverage_status\nP\treference_only_family\n", encoding="utf-8")
            manual.write_text(
                "candidate_id\tmanual_review_decision\tnext_action\tsuperfamily\nA\tretain_candidate_only\tinspect\tmanual_sf\n",
                encoding="utf-8",
            )
            result = module.write_reports(matrix, profile_gap, manual, output_dir, "test_run")
            self.assertEqual(result["candidate_count"], 1)
            self.assertTrue((output_dir / "conflict_ambiguity_summary.tsv").exists())
            self.assertTrue((output_dir / "priority_19_candidate_summary.tsv").exists())
            self.assertTrue((output_dir / "reconciliation_report.md").exists())
            focus = module.read_tsv(output_dir / "priority_19_candidate_summary.tsv")
            self.assertEqual(focus[0]["superfamily"], "manual_sf")
            gap = module.read_tsv(output_dir / "profile_gap_plan.tsv")
            self.assertTrue(gap[0]["planned_action"].startswith("planned_not_run"))
            self.assertTrue((output_dir.parent / "input_contract.json").exists())
            self.assertTrue((output_dir.parent / "run_manifest.json").exists())
            self.assertTrue((output_dir.parent / "inputs" / "source_manifest.tsv").exists())
            contract = json.loads((output_dir.parent / "input_contract.json").read_text(encoding="utf-8"))
            self.assertIn("finalizer", contract["source_snapshot"])
            self.assertIn("agents_rules", contract["source_snapshot"])


if __name__ == "__main__":
    unittest.main()
