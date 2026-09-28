import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "integrate_phaded_structure_evidence.py"
spec = importlib.util.spec_from_file_location("integrate_phaded_structure_evidence", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_review_tier_keeps_structure_quality_separate_from_architecture_conflict():
    row = module.build_review_row(
        {
            "source_accession": "A",
            "panel": "extracellular_explicit",
            "phylo_decision": "retain_extracellular_dPHASCL2_like_candidate",
            "superfamily": "extracellular dPHASCL type 2",
            "best_model_mean_plddt": "90.0",
            "best_model_ptm": "0.9",
            "structure_confidence_band": "high_structure_confidence",
        },
        {
            "accession": "A",
            "sequence_integrity": "terminal_stop_only",
            "signalp_class": "SP",
            "lipase_box_state": "present",
            "catalytic_ser_cys_state": "supported",
        },
        {
            "accession": "A",
            "architecture_consistency": "conflicting",
            "pfam_accessions": "PF00756",
            "assignment_review": "review_required",
        },
        {
            "source_accession": "A",
            "panel": "extracellular_explicit",
            "phylo_decision": "retain_extracellular_dPHASCL2_like_candidate",
            "nearest_reference_family": "DED_hfam_70",
            "nearest_reference_distance": "0.4",
            "candidate_reference_mrca_support": "98.0",
        },
    )
    assert row["structure_quality"] == "high"
    assert row["architecture_evidence"] == "conflicting"
    assert row["review_tier"] == "architecture_conflict_manual_review"
    assert row["phenotype_claim"] == "candidate_only"


def test_review_tier_holds_low_confidence_and_n_truncation():
    row = module.build_review_row(
        {
            "source_accession": "B",
            "panel": "ephaz_competition",
            "phylo_decision": "retain_dPHASCL1_like_candidate",
            "superfamily": "extracellular dPHASCL type 1",
            "best_model_mean_plddt": "64.0",
            "best_model_ptm": "0.67",
            "structure_confidence_band": "low_structure_confidence",
        },
        {"accession": "B", "sequence_integrity": "possible_N_truncation", "signalp_class": "OTHER"},
        {"accession": "B", "architecture_consistency": "partial", "pfam_accessions": "PF00326"},
        {"source_accession": "B", "nearest_reference_family": "DED_hfam_45", "nearest_reference_distance": "0.8", "candidate_reference_mrca_support": "92"},
    )
    assert row["review_tier"] == "low_confidence_or_truncation_hold"
    assert row["phenotype_claim"] == "candidate_only"


def test_summarize_pdb_reports_residue_confidence_and_low_run(tmp_path):
    pdb = tmp_path / "x.pdb"
    lines = []
    for i, b in enumerate((90.0, 65.0, 60.0, 85.0), 1):
        lines.append(f"ATOM  {i:5d}  CA  ALA A{i:4d}    {i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{b:6.2f}           C")
    pdb.write_text("\n".join(lines) + "\n", encoding="ascii")
    result = module.summarize_pdb(pdb)
    assert result["residue_count"] == "4"
    assert result["mean_plddt"] == "75.00"
    assert result["low_fraction"] == "0.500"
    assert result["longest_low_run"] == "2"


def test_terminal_stop_and_unbound_signalp_are_informational_only():
    row = module.build_review_row(
        {"source_accession": "C", "panel": "ephaz_competition", "superfamily": "x", "best_model_mean_plddt": "85", "structure_confidence_band": "high_structure_confidence"},
        {"accession": "C", "sequence_integrity": "terminal_stop_only", "signalp_class": "OTHER", "lipase_box_state": "present", "catalytic_ser_cys_state": "supported"},
        {"accession": "C", "architecture_consistency": "partial", "pfam_accessions": "PF10503", "assignment_status": "assigned", "family_score_gap": "20"},
        {"source_accession": "C", "nearest_reference_family": "DED_hfam_45", "candidate_reference_mrca_support": "95", "nearest_reference_distance": "0.5"},
        {"accession": "C", "interpro_secondary_state": "explicit_pha_related_domain"},
    )
    assert row["merged_priority"] == "P1_high_architecture_phylo_candidate"
    assert "TERMINAL_STOP_CLEANED" in row["conflict_codes"]


def test_pfam_span_is_derived_from_coordinates():
    row = module.build_review_row(
        {"source_accession": "D", "panel": "ephaz_competition", "superfamily": "x", "sequence_length_cleaned": "200", "best_model_mean_plddt": "85", "structure_confidence_band": "high_structure_confidence"},
        {"accession": "D", "sequence_integrity": "terminal_stop_only"},
        {"accession": "D", "architecture_consistency": "partial", "pfam_accessions": "PF10503", "pfam_coordinates": "PF10503:10-109"},
        {"source_accession": "D", "candidate_reference_mrca_support": "95"},
    )
    assert row["pfam_span_start"] == "10"
    assert row["pfam_span_end"] == "109"
    assert row["pfam_span_coverage"] == "0.500"
