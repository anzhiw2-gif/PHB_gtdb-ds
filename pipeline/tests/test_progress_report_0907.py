import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "pipeline" / "scripts" / "plot_progress_report_0907.py"


def load_module():
    spec = importlib.util.spec_from_file_location("plot_progress_report_0907", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_report_scope_is_limited_to_the_two_main_families():
    module = load_module()
    assert module.MAIN_FAMILIES == ("ePhaZ", "iPhaZ")
    assert module.validate_main_family_scope(pd.DataFrame({"family": ["ePhaZ", "iPhaZ"]})) == []
    assert module.validate_main_family_scope(pd.DataFrame({"family": ["ArchPhaZ_patatin"]})) == [
        "ArchPhaZ_patatin"
    ]


def test_report_has_five_main_figures_and_no_supplementary_branch():
    module = load_module()
    stems = module.report_figure_stems()
    assert len(stems) == 5
    assert all("supplementary" not in stem for stem in stems)
    assert all("archae" not in stem.lower() for stem in stems)


def test_workflow_excludes_archaeal_side_branch():
    module = load_module()
    workflow = module.workflow_table()
    assert not workflow["step"].str.contains("archae", case=False).any()
    assert set(workflow["scope"]) == {"main_line"}


def test_signalp_scope_is_ephaz_curated_core_only():
    module = load_module()
    frame = pd.DataFrame(
        {
            "family": ["ePhaZ", "ePhaZ", "iPhaZ"],
            "source_layer": ["ePhaZ_curated_core", "ePhaZ_tier2", "iPhaZ_tier1"],
        }
    )
    result = module.signalp_scope_counts(frame)
    assert result == {"input_count": 1, "excluded_ephaz_review": 1, "excluded_iphaz": 1}


def test_initial_vs_tier1_counts_keep_family_and_evidence_boundaries():
    module = load_module()
    frame = pd.DataFrame(
        {
            "family": ["ePhaZ", "ePhaZ", "iPhaZ", "iPhaZ"],
            "source_layer": ["ePhaZ_curated_core", "ePhaZ_tier2", "iPhaZ_tier1", "iPhaZ_tier2"],
        }
    )
    result = module.initial_vs_tier1_counts(frame)
    assert result.to_dict("records") == [
        {"family": "ePhaZ", "initial_candidates": 2, "tier1_candidates": 1},
        {"family": "iPhaZ", "initial_candidates": 2, "tier1_candidates": 1},
    ]


def test_genus_distribution_preserves_family_breakdown():
    module = load_module()
    frame = pd.DataFrame(
        {
            "genome": ["g1", "g2", "g2", "g3"],
            "family": ["ePhaZ", "ePhaZ", "iPhaZ", "iPhaZ"],
            "gtdb_taxonomy": [
                "d__Bacteria;p__P;c__C;o__O;f__F;g__Alpha;s__a",
                "d__Bacteria;p__P;c__C;o__O;f__F;g__Alpha;s__b",
                "d__Bacteria;p__P;c__C;o__O;f__F;g__Alpha;s__b",
                "d__Bacteria;p__P;c__C;o__O;f__F;g__Beta;s__c",
            ],
        }
    )
    result = module.genus_distribution_table(frame)
    assert result.to_dict("records") == [
        {"genus": "Alpha", "family": "ePhaZ", "genomes": 2},
        {"genus": "Alpha", "family": "iPhaZ", "genomes": 1},
        {"genus": "Beta", "family": "iPhaZ", "genomes": 1},
    ]


def test_figure4_uses_genome_family_assignment_unit():
    module = load_module()
    assert module.figure4_axis_label() == "Tier1 genome-family assignments"


def test_figure4_marks_gtdb_placeholder_and_candidate_only_boundary():
    module = load_module()
    assert module.format_figure4_genus_label("Palsa-465") == "Palsa-465 (GTDB placeholder)"
    assert module.format_figure4_genus_label("Alpha") == "Alpha"
    caption = module.figure4_caption()
    assert "candidate-only" in caption
    assert "taxonomy" in caption.lower()
    assert "HMM" in caption
