import importlib.util
import unittest
from pathlib import Path
import tempfile


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_phaded_motifs.py"


def load_module():
    spec = importlib.util.spec_from_file_location("audit_phaded_motifs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MotifAuditTests(unittest.TestCase):
    def test_type1_requires_upstream_oxyanion_pattern_and_separate_catalytic_residues(self):
        module = load_module()
        seq = "M" + "A" * 20 + "HAGCQ" + "A" * 60 + "GLSAG" + "A" * 14 + "GXXDYTV" + "A" * 45 + "GMXHXXPXXG"
        result = module.audit_sequence("A", seq, "extracellular dPHASCL type 1", pfam_accessions="PF10503;PF06850")
        self.assertEqual(result["lipase_box_state"], "supported")
        self.assertEqual(result["catalytic_ser_cys_state"], "supported")
        self.assertEqual(result["asp_state"], "supported")
        self.assertEqual(result["his_state"], "supported")
        self.assertEqual(result["oxyanion_hole_state"], "supported")
        self.assertEqual(result["sbd_state"], "supported")
        self.assertEqual(result["motif_evidence_level"], "motif_pattern_only")

    def test_nphascl_without_lipase_box_does_not_turn_missing_box_into_negative(self):
        module = load_module()
        result = module.audit_sequence("B", "M" + "A" * 400, "intracellular nPHASCL without lipase box", pfam_accessions="PF00561")
        self.assertEqual(result["lipase_box_state"], "not_expected_for_subtype")
        self.assertEqual(result["catalytic_ser_cys_state"], "pending_reference_annotation")
        self.assertIn("cys_his_asp_reference_mapping_pending", result["motif_pending_reason"])

    def test_dphamcl_marks_sbd_and_linker_not_expected_but_lid_is_pending(self):
        module = load_module()
        result = module.audit_sequence("C", "M" + "A" * 300, "extracellular dPHAMCL", pfam_accessions="PF10503")
        self.assertEqual(result["sbd_state"], "not_expected_for_subtype")
        self.assertEqual(result["linker_state"], "not_expected_for_subtype")
        self.assertEqual(result["lid_state"], "pending_reference_annotation")

    def test_type2_flags_oxyanion_on_wrong_side_as_conflict(self):
        module = load_module()
        seq = "M" + "A" * 10 + "HAGCQ" + "A" * 20 + "GLSAG" + "A" * 250
        result = module.audit_sequence("D", seq, "extracellular dPHASCL type 2", pfam_accessions="PF10503")
        self.assertEqual(result["oxyanion_hole_state"], "conflict_relative_position")

    def test_alignment_quality_distinguishes_header_only_files(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            header_only = tmp_path / "aln44.aln"
            header_only.write_text("CLUSTAL W (1.8) multiple sequence alignment\n\n", encoding="ascii")
            usable = tmp_path / "aln55.aln"
            usable.write_text("CLUSTAL W (1.8) multiple sequence alignment\n\nA  ABC\nB  ABD\n", encoding="ascii")
            self.assertEqual(module.inspect_alignment(header_only)["alignment_status"], "empty_or_header_only")
            self.assertEqual(module.inspect_alignment(usable)["alignment_status"], "usable_reference_alignment")


if __name__ == "__main__":
    unittest.main()
