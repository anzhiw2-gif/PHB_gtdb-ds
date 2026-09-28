import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "review_phaded_candidate_architecture.py"
spec = importlib.util.spec_from_file_location("review_phaded_candidate_architecture", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_classify_dphascl2_pfam_conflict_as_hold():
    result = module.classify_candidate({
        "source_accession": "A",
        "superfamily": "extracellular dPHASCL type 2",
        "pfam_accessions": "PF00756",
        "architecture_evidence": "conflicting",
        "interpro_secondary_state": "explicit_pha_related_domain",
        "structure_quality": "high",
        "pdb_n_terminal_20_mean": "26.46",
        "pdb_low_fraction": "0.192",
        "sequence_integrity": "terminal_stop_only",
    })
    assert result["review_decision"] == "hold_architecture_conflict"
    assert "PFAM_OTHER_SUPERFAMILY" in result["decision_reason"]


def test_classify_short_low_confidence_record_as_gene_model_hold():
    result = module.classify_candidate({
        "source_accession": "B",
        "superfamily": "extracellular dPHASCL type 1",
        "pfam_accessions": "PF10503",
        "architecture_evidence": "partial",
        "interpro_secondary_state": "explicit_pha_related_domain",
        "structure_quality": "low",
        "pdb_low_fraction": "0.490",
        "pdb_longest_low_run": "85",
        "sequence_integrity": "terminal_stop_only",
        "sequence_length_cleaned": "202",
    })
    assert result["review_decision"] == "hold_gene_model_or_structure"


def test_classify_complete_high_candidate_for_comparison():
    result = module.classify_candidate({
        "source_accession": "C",
        "superfamily": "extracellular dPHASCL type 1",
        "pfam_accessions": "PF00326;PF10503",
        "architecture_evidence": "partial",
        "interpro_secondary_state": "explicit_pha_related_domain",
        "structure_quality": "high",
        "pdb_low_fraction": "0.196",
        "sequence_integrity": "terminal_stop_only",
        "sequence_length_cleaned": "363",
    })
    assert result["review_decision"] == "retain_candidate_comparison"
