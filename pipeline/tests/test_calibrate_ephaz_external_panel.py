import importlib.util
from pathlib import Path


def load_module():
    path = Path(__file__).parents[1] / "scripts" / "calibrate_ephaz_external_panel.py"
    spec = importlib.util.spec_from_file_location("calibrate_ephaz_external_panel", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_panel_summary_keeps_four_denominators_and_challenge_separate():
    module = load_module()
    panels = {
        "PHB_experimental": ["P1", "P2"],
        "MCL-PHA_experimental": ["M1"],
        "intracellular_PHB": ["I1", "I2"],
        "MCL-PHA_non_PHB": ["N1"],
        "annotation_only": ["A1", "A2"],
        "fragment_challenge": ["F1"],
    }
    hits = {acc: {"evalue": value} for acc, value in {
        "P1": 1e-40, "M1": 1e-20, "I1": 1e-30, "A1": 1e-10, "F1": 1e-50
    }.items()}
    rows = module.summarize_panel_hits(panels, hits, threshold=1e-5)
    by_panel = {row["panel"]: row for row in rows}
    assert by_panel["PHB_experimental"]["tested"] == 2
    assert by_panel["PHB_experimental"]["detected"] == 1
    assert by_panel["MCL-PHA_experimental"]["detected"] == 1
    assert by_panel["intracellular_PHB"]["FP"] == 1
    assert by_panel["MCL-PHA_non_PHB"]["FP"] == 0
    assert by_panel["annotation_only"]["FP"] == 1
    assert by_panel["fragment_challenge"]["challenge_detected"] == 1
    assert by_panel["fragment_challenge"]["formal_denominator"] is False


def test_parse_domtblout_reports_hmm_and_target_coverage():
    module = load_module()
    line = "target  -  100 query - 50 1e-20 100.0 0.0 1 1 1e-20 1e-20 100.0 0.0 1 50 1 50 1 100 0.99 desc"
    result = module.parse_domtblout([line], {"target": 100})
    assert result["target"]["hmm_coverage"] == 1.0
    assert result["target"]["target_coverage"] == 0.5


def test_sensitivity_grid_emits_each_threshold_and_coverage_for_each_panel():
    module = load_module()
    panels = {"PHB_experimental": ["P1"], "intracellular_PHB": ["I1"]}
    hits = {
        "P1": {"evalue": 1e-8, "hmm_coverage": 0.7},
        "I1": {"evalue": 1e-8, "hmm_coverage": 0.3},
    }
    rows = module.summarize_sensitivity_grid(panels, hits, [1e-5, 1e-10], [0.0, 0.6])
    assert len(rows) == 8
    strict = [r for r in rows if r["threshold"] == 1e-10 and r["min_hmm_coverage"] == 0.6]
    assert strict[0]["detected"] == 0
    assert strict[1]["detected"] == 0
