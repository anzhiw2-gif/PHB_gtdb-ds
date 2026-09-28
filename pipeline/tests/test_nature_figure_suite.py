import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "pipeline" / "scripts" / "plot_reframed_nature_figure_suite.py"


def load_module():
    spec = importlib.util.spec_from_file_location("plot_nature_figure_suite", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_subtype_counts_are_deterministic_and_preserve_candidate_layers():
    module = load_module()
    frame = pd.DataFrame(
        {
            "pha_subtype_evidence": [
                "ePhaZ_core_nonsecreted_evidence",
                "ePhaZ_core_nonsecreted_evidence",
                "iPhaZ_intracellular_evidence",
            ]
        }
    )
    result = module.subtype_counts(frame)
    assert list(result["category"]) == [
        "ePhaZ core nonsecreted",
        "iPhaZ intracellular",
    ]
    assert list(result["count"]) == [2, 1]


def test_taxonomy_species_parser_uses_species_rank_only():
    module = load_module()
    assert module.species_from_taxonomy(
        "d__Archaea;p__Halobacteriota;c__Halobacteria;g__Haloferax;s__Haloferax mediterranei"
    ) == "Haloferax mediterranei"
    assert module.species_from_taxonomy("d__Archaea;p__Thermoproteota;c__Nitrososphaeria") == "unclassified"


def test_validate_export_bundle_requires_four_formats_per_figure(tmp_path):
    module = load_module()
    for stem in ("figure_1", "figure_2"):
        for suffix in (".svg", ".pdf", ".tiff", ".png"):
            (tmp_path / f"{stem}{suffix}").write_bytes(b"x")
    assert module.validate_export_bundle(tmp_path, ["figure_1", "figure_2"]) == []
    (tmp_path / "figure_2.tiff").unlink()
    assert module.validate_export_bundle(tmp_path, ["figure_1", "figure_2"]) == ["figure_2.tiff"]


def test_resolve_input_tables_binds_renderer_to_this_run_input_snapshot(tmp_path):
    module = load_module()
    for name in (
        "overall_pha_phb.tsv",
        "archaeal_phazh1_full_audit.tsv",
        "domain_annotation.tsv",
        "positive_negative_calibration.tsv",
        "cross_phylum_manual_review.tsv",
        "tier1_genome_family.tsv",
        "no_signal_competition.tsv",
    ):
        (tmp_path / "inputs").mkdir(exist_ok=True)
        (tmp_path / "inputs" / name).write_text("header\nvalue\n", encoding="utf-8")
    paths = module.resolve_input_tables(tmp_path)
    assert paths["overall"].parent == tmp_path / "inputs"
    assert paths["archaea"].name == "archaeal_phazh1_full_audit.tsv"
    assert paths["calibration"].name == "positive_negative_calibration.tsv"


def test_bar_label_uses_inside_white_text_for_large_bar():
    module = load_module()
    assert module.bar_label(32226, 5000) == ("32,226", "inside", "white")


def test_bar_label_uses_outside_dark_text_for_small_bar():
    module = load_module()
    assert module.bar_label(18, 5000) == ("18", "outside", module.PALETTE["dark"])


def test_main_family_scope_is_only_ephaz_and_iphaz():
    module = load_module()
    assert module.MAIN_FAMILIES == ("ePhaZ", "iPhaZ")
    frame = pd.DataFrame({"family": ["ePhaZ", "iPhaZ"]})
    assert module.validate_main_family_scope(frame) == []
    contaminated = pd.DataFrame({"family": ["ePhaZ", "ArchPhaZ_patatin"]})
    assert module.validate_main_family_scope(contaminated)


def test_signalp_scope_is_curated_ephaz_only():
    module = load_module()
    frame = pd.DataFrame(
        {
            "family": ["ePhaZ", "ePhaZ", "iPhaZ"],
            "source_layer": ["ePhaZ_curated_core", "tier2", "tier1"],
            "signal_type": ["SP", "OTHER", ""],
        }
    )
    result = module.signalp_scope_counts(frame)
    assert result["input_count"] == 1
    assert result["excluded_iPhaZ"] == 1
    assert result["excluded_ePhaZ_tier2"] == 1


def test_taxonomy_distribution_uses_tier1_genome_table():
    module = load_module()
    frame = pd.DataFrame(
        {
            "genome": ["G1", "G1", "G2"],
            "family": ["ePhaZ", "iPhaZ", "ePhaZ"],
            "phylum": ["P1", "P1", "P2"],
            "gtdb_taxonomy": [
                "d__Bacteria;p__P1;g__A;s__A one",
                "d__Bacteria;p__P1;g__A;s__A one",
                "d__Bacteria;p__P2;g__B;s__B one",
            ],
        }
    )
    result = module.genome_distribution(frame)
    assert set(result["family"]) == {"ePhaZ", "iPhaZ"}
    assert result["genomes"].sum() == 3


def test_workflow_main_line_excludes_archaeal_side_branch():
    module = load_module()
    workflow = module.workflow_table()
    assert not workflow["step"].str.contains("Archaeal", case=False).any()
    assert "side_branch" in workflow.columns


def test_competition_display_keeps_minority_assignments_readable_without_rescaling():
    module = load_module()
    frame = pd.DataFrame(
        {
            "assignment": ["ePhaZ_like", "iPhaZ_like", "ambiguous"],
            "proteins": [16809, 9, 18],
        }
    )
    majority, minority = module.competition_display_groups(frame)
    assert majority.to_dict("records") == [{"assignment": "ePhaZ_like", "proteins": 16809}]
    assert minority.to_dict("records") == [
        {"assignment": "iPhaZ_like", "proteins": 9},
        {"assignment": "ambiguous", "proteins": 18},
    ]
