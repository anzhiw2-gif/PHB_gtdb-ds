import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "pipeline" / "scripts" / "overall_pha_phb_archaea_audit.py"


def load():
    spec = importlib.util.spec_from_file_location("overall_pha_phb_archaea_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_overall_layer_is_preserved_and_phb_focus_is_conservative():
    module = load()
    row = {
        "family": "ePhaZ",
        "layer": "ePhaZ_curated_core_secreted",
        "assignment": "ePhaZ_like",
        "signal_type": "SP",
        "length": "430",
    }
    result = module.classify_pha_row(row)
    assert result["overall_pha_layer"] == "ePhaZ_curated_core_secreted"
    assert result["pha_subtype_evidence"] == "ePhaZ_SCL_secreted_evidence"
    assert result["phb_focus"] == "PHB_candidate_review"
    assert result["phb_exclusion_reason"] == "no_subtype_evidence"


def test_nonsecreted_layer_is_not_classified_as_secreted():
    module = load()
    result = module.classify_pha_row({
        "family": "ePhaZ",
        "layer": "ePhaZ_curated_core_nonsecreted",
        "assignment": "ePhaZ_like",
    })
    assert result["pha_subtype_evidence"] == "ePhaZ_core_nonsecreted_evidence"


def test_mcl_assignment_is_excluded_from_phb_main_result():
    module = load()
    result = module.classify_pha_row({
        "family": "iPhaZ",
        "layer": "iPhaZ_tier1",
        "assignment": "MCL_lipase_associated",
        "signal_type": "",
        "length": "410",
    })
    assert result["phb_focus"] == "MCL_or_nonPHB_review"
    assert "MCL" in result["phb_exclusion_reason"]
    assert result["pha_subtype_evidence"] == "iPhaZ_MCL_competition_evidence"


def test_dual_profile_mcl_like_label_is_loaded_as_competition_evidence(tmp_path):
    module = load()
    dual = tmp_path / "dual_profile.tsv"
    dual.write_text(
        "accession\tclassification\nMCL1\tMCL_like\nPHB1\tPHB_like\n",
        encoding="utf-8",
    )
    assert module._accessions(dual, "MCL_lipase_associated") == {"MCL1"}


def test_review_layers_keep_their_evidence_boundary_in_subtype_label():
    module = load()
    result = module.classify_pha_row({
        "family": "ePhaZ",
        "layer": "ePhaZ_tier2_review",
        "assignment": "",
    })
    assert result["pha_subtype_evidence"] == "ePhaZ_tier2_review_evidence"


def test_archaeal_high_hit_requires_full_audit_fields():
    module = load()
    result = module.audit_archaeal_row({
        "genome": "G1", "locus": "L1", "domain": "Archaea",
        "coverage": "0.91", "evalue": "1e-30", "complete": "true",
        "nearby_markers": "BdhA", "patatin_conflict": "false",
        "phylum": "Halobacteriota", "class": "Halobacteria",
        "sequence_present": "true", "domain_annotation": "Patatin",
    })
    assert result["layer"] == "PhaZh1_like_high"
    assert result["audit_status"] == "review_required"
    assert result["phenotype_claim"] == "candidate_only"


def test_archaeal_nonarchaea_and_conflict_are_not_high():
    module = load()
    result = module.audit_archaeal_row({
        "genome": "G2", "locus": "L2", "domain": "Bacteria",
        "coverage": "0.99", "evalue": "1e-50", "complete": "true",
        "nearby_markers": "PhaC", "patatin_conflict": "true",
        "phylum": "Pseudomonadota", "class": "Gammaproteobacteria",
        "sequence_present": "false", "domain_annotation": "pending",
    })
    assert result["layer"] == "non_PhaZh1_patatin"
    assert result["audit_status"] == "excluded"


def test_fasta_accessions_are_parsed_for_full_archaeal_audit(tmp_path):
    module = load()
    fasta = tmp_path / "arch.faa"
    fasta.write_text(">G1|L1\nMPEPTIDE*\n>G2|L2\nMSEQ\n", encoding="utf-8")
    assert module.fasta_accessions(fasta) == {"G1|L1", "G2|L2"}


def test_n_terminal_patatin_motif_is_audit_field_not_a_phenotype_call(tmp_path):
    module = load()
    fasta = tmp_path / "arch.faa"
    fasta.write_text(
        ">G1|L1\nMNNGASAGAAAA\n>G2|L2\nMNNAAAAAA\n",
        encoding="utf-8",
    )
    assert module.n_terminal_gxsxg(fasta) == {"G1|L1": "present", "G2|L2": "absent"}


def test_run_accepts_an_empty_precreated_results_directory(tmp_path):
    module = load()
    layers = tmp_path / "layers.tsv"
    layers.write_text(
        "accession\tfamily\tlayer\tassignment\nA1\tePhaZ\tePhaZ_curated_core_secreted\tePhaZ_like\n",
        encoding="utf-8",
    )
    arch = tmp_path / "arch.tsv"
    arch.write_text(
        "genome\tlocus\tdomain\tcoverage\tevalue\tcomplete\tnearby_markers\tpatatin_conflict\tlayer\n",
        encoding="utf-8",
    )
    taxonomy = tmp_path / "taxonomy.tsv"
    taxonomy.write_text("genome\tphylum\tclass\tgtdb_taxonomy\n", encoding="utf-8")
    context = tmp_path / "context.tsv"
    context.write_text("hit_locus\tnearby_markers\n", encoding="utf-8")
    output = tmp_path / "results"
    output.mkdir()
    manifest = module.run(layers, arch, taxonomy, context, output)
    assert manifest["counts"]["pha_records"] == 1
    assert (output / "overall_pha_phb.tsv").exists()
