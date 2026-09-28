import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "pipeline" / "scripts" / "build_archaea_manual_review.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_archaea_manual_review", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_build_selects_all_review_and_cross_phylum_records_with_context(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text(
        "genome\tlocus\tlayer\tphylum\tclass\tcoverage\tevalue\tnterm_gxsxg_motif\tregistered_model_hits\tregistered_model_best_bitscore\tnearby_markers\n"
        "g1\tp1\tPhaZh1_like_high\tHalobacteriota\tHalobacteria\t0.9\t1e-20\tpresent\tArchPhaZ_patatin\tArchPhaZ_patatin:100\tBdhA\n"
        "g2\tp2\tPhaZh1_like_review\tHalobacteriota\tHalobacteria\t0.7\t1e-8\tabsent\tArchPhaZ_patatin\tArchPhaZ_patatin:50\tBdhA\n"
        "g3\tp3\tPhaZh1_like_high\tThermoproteota\tNitrososphaeria\t0.8\t1e-10\tpresent\tArchPhaZ_patatin\tArchPhaZ_patatin:60\tBdhA\n",
        encoding="utf-8",
    )
    context = tmp_path / "context.tsv"
    context.write_text(
        "genome\tcontig\thit_locus\thit_family\tstart\tend\tstrand\tnearby_markers\n"
        "g2\tc2\tp2\tArchPhaZ_patatin\t10\t20\t-\tBdhA\n"
        "g3\tc3\tp3\tArchPhaZ_patatin\t30\t40\t+\tBdhA\n",
        encoding="utf-8",
    )
    output = tmp_path / "review.tsv"
    assert load_module().build(annotation, context, output) == 2
    rows = output.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 3
    assert "g2" in rows[1] and "manual_review_required" in rows[1]
    assert "g3" in rows[2] and "c3" in rows[2]
