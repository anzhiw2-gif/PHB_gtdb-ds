from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "runs" / "20260906_nature_figure_suite_03" / "results" / "nature_figure_suite_explanation.docx"


def test_figure_suite_explanation_contains_all_eight_figures_and_boundary():
    assert DOC.is_file() and DOC.stat().st_size > 0
    document = Document(DOC)
    text = "\n".join(paragraph.text for paragraph in document.paragraphs)
    assert "6 张主图 + 2 张补充图" in text
    assert "candidate-only" in text
    assert all(f"图 {index}" in text for index in range(1, 9))
    assert len(document.inline_shapes) == 8
