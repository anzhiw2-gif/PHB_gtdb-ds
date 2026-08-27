"""Tests for Scheme A figure quality-assurance checks."""
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def load_module(name: str):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FigureQATests(unittest.TestCase):
    def test_qa_rejects_svg_without_editable_text(self):
        figure_qa = load_module("check_scheme_a_figures")
        with tempfile.TemporaryDirectory() as tmp:
            svg = Path(tmp) / "broken.svg"
            svg.write_text("<svg><path /></svg>", encoding="utf-8")
            with self.assertRaises(figure_qa.FigureQAError):
                figure_qa.check_svg_text(svg)


if __name__ == "__main__":
    unittest.main()
