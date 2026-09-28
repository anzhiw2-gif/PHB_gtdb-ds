import importlib.util
import unittest

from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_phaded_subtype_matrix.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_phaded_subtype_matrix", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def feature(accession="A", signalp="accepted_layer_signalp", signal_class="SP", lipase="present"):
    return {
        "accession": accession, "sequence_integrity": "terminal_stop_only",
        "signalp_evidence": signalp, "signalp_class": signal_class,
        "lipase_box_state": lipase, "catalytic_ser_cys_state": "supported" if lipase == "present" else "pending",
        "his_state": "pending", "asp_state": "pending", "oxyanion_hole_state": "pending",
        "sbd_state": "pending", "linker_state": "pending", "lid_state": "pending",
    }


class SubtypeMatrixTests(unittest.TestCase):
    def test_profile_domain_and_localization_are_explicit_and_motif_is_partial(self):
        module = load_module()
        rows = module.build_matrix(
            [{"accession": "A", "genome": "G", "layer": "ePhaZ_curated_core_secreted", "phaded_superfamily_best": "extracellular dPHASCL type 1", "phaded_family_best": "family_DED_hfam_55_hash", "assignment_status": "assigned", "superfamily_score_gap": "10", "family_score_gap": "5"}],
            [feature()],
            [{"accession": "A", "architecture_consistency": "partial", "pfam_interpro_state": "tested_targeted_pfam_hit", "pfam_accessions": "PF10503"}],
            [{"accession": "A", "interpro_status": "interpro_supported", "interpro_secondary_state": "explicit_pha_related_domain", "interpro_hit_count": "2", "interpro_ids": "IPR010126"}],
            [{"profile_kind": "family", "phaded_family_id": "DED_hfam_55", "phaded_superfamily": "extracellular dPHASCL type 1", "model_status": "trained"}],
        )
        row = rows[0]
        self.assertEqual(row["profile_evidence_status"], "profile_trained_hit")
        self.assertEqual(row["domain_evidence_status"], "domain_supported_partial")
        self.assertEqual(row["interpro_status"], "interpro_supported")
        self.assertEqual(row["interpro_evidence_status"], "interpro_supported")
        self.assertEqual(row["motif_evidence_status"], "partial_lipase_box_only")
        self.assertEqual(row["localization_evidence_status"], "export_signal_supported")
        self.assertEqual(row["subtype_call"], "dPHASCL1_like_candidate")
        self.assertIn(row["subtype_confidence"], {"moderate_candidate_only", "low_candidate_only"})

    def test_full_library_signalp_supported_status_counts_as_localization_evidence(self):
        module = load_module()
        rows = module.build_matrix(
            [{"accession": "A", "genome": "G", "layer": "ePhaZ", "phaded_superfamily_best": "extracellular dPHASCL type 1", "assignment_status": "assigned"}],
            [feature(signalp="signalp_supported", signal_class="SP")],
            [{"accession": "A", "architecture_consistency": "partial"}],
            [{"accession": "A", "interpro_status": "interpro_supported", "interpro_secondary_state": "explicit_pha_related_domain"}],
            [],
        )
        self.assertEqual(rows[0]["localization_evidence_status"], "export_signal_supported")

    def test_conflict_and_unassigned_are_not_promoted(self):
        module = load_module()
        rows = module.build_matrix(
            [
                {"accession": "A", "genome": "G", "layer": "ePhaZ_tier2_review", "phaded_superfamily_best": "extracellular dPHASCL type 1", "phaded_family_best": "", "assignment_status": "assigned", "superfamily_score_gap": "2", "family_score_gap": "0"},
                {"accession": "B", "genome": "G2", "layer": "iPhaZ_tier1", "phaded_superfamily_best": "", "phaded_family_best": "", "assignment_status": "unassigned_PhaDED_like", "superfamily_score_gap": "0", "family_score_gap": "0"},
            ],
            [feature("A", signalp="not_bound_in_layer_contract", signal_class="pending", lipase="present"), feature("B", signalp="not_bound_in_layer_contract", signal_class="pending", lipase="not_detected")],
            [{"accession": "A", "architecture_consistency": "conflicting", "pfam_interpro_state": "tested_other_superfamily_hit", "pfam_accessions": "PF00756"}],
            [{"accession": "A", "interpro_status": "interpro_supported", "interpro_secondary_state": "explicit_pha_related_domain", "interpro_hit_count": "1", "interpro_ids": "IPR010126"}, {"accession": "B", "interpro_status": "interpro_not_reported", "interpro_secondary_state": "no_interpro_hit", "interpro_hit_count": "0", "interpro_ids": ""}],
            [],
        )
        by = {row["accession"]: row for row in rows}
        self.assertEqual(by["A"]["subtype_call"], "hold_architecture_conflict")
        self.assertEqual(by["A"]["subtype_confidence"], "hold")
        self.assertEqual(by["B"]["subtype_call"], "unassigned_PhaDED_like")
        self.assertEqual(by["B"]["profile_evidence_status"], "profile_unassigned")
        self.assertEqual(by["B"]["motif_evidence_status"], "not_tested_full_library")

    def test_base_ledger_fields_are_preserved_while_interpro_status_is_replaced(self):
        module = load_module()
        rows = module.build_matrix(
            [{"accession": "A", "genome": "G", "layer": "ePhaZ_tier2_review", "phaded_superfamily_best": "extracellular dPHASCL type 1", "phaded_family_best": "", "assignment_status": "assigned"}],
            [feature()],
            [{"accession": "A", "architecture_consistency": "partial", "pfam_interpro_state": "tested_targeted_pfam_hit", "pfam_accessions": "PF10503"}],
            [{"accession": "A", "interpro_status": "interpro_supported", "interpro_secondary_state": "explicit_pha_related_domain", "interpro_hit_count": "1", "interpro_ids": "IPR010126"}],
            [],
            base_ledger=[{"accession": "A", "evidence_grade": "L1_profile_plus_partial_architecture", "literature_evidence_tier": "none_exact_candidate_match"}],
        )
        self.assertEqual(rows[0]["evidence_grade"], "L1_profile_plus_partial_architecture")
        self.assertEqual(rows[0]["literature_evidence_tier"], "none_exact_candidate_match")
        self.assertEqual(rows[0]["interpro_evidence_status"], "interpro_supported")

    def test_motif_pattern_panel_states_are_preserved_without_calling_them_complete_catalysis(self):
        module = load_module()
        row = module.build_matrix(
            [{"accession": "A", "genome": "G", "layer": "ePhaZ", "phaded_superfamily_best": "extracellular dPHASCL type 1", "assignment_status": "assigned"}],
            [{
                **feature(), "motif_reference_subtype": "extracellular dPHASCL type 1",
                "motif_panel_status": "partial_catalytic_pattern_panel",
                "motif_completeness": "pattern_partial",
                "lipase_box_state": "supported", "catalytic_ser_cys_state": "supported",
                "his_state": "supported", "asp_state": "not_detected_pattern",
                "oxyanion_hole_state": "supported",
            }],
            [{"accession": "A", "architecture_consistency": "partial"}],
            [{"accession": "A", "interpro_status": "interpro_supported", "interpro_secondary_state": "explicit_pha_related_domain"}],
            [],
        )[0]
        self.assertEqual(row["motif_evidence_status"], "partial_catalytic_pattern_panel")
        self.assertEqual(row["motif_completeness"], "pattern_partial")


if __name__ == "__main__":
    unittest.main()
