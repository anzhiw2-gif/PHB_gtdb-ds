"""Export tests for Scheme A matplotlib figures."""
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


class FigureExportTests(unittest.TestCase):
    def test_figure1_separates_genome_and_protein_levels(self):
        figures = load_module("plot_scheme_a_figures")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data_dir = root / "source_data"
            data_dir.mkdir()
            (data_dir / "figure1_genome_coverage.tsv").write_text(
                "stage\tvalue\tunit\tnote\n"
                "GTDB representative genomes\t199923\tgenomes\tInput census\n"
                "Genomes with >=1 tier1 core candidate\t44814\tgenomes\tCore-family union\n",
                encoding="utf-8",
            )
            (data_dir / "figure1_funnel.tsv").write_text(
                "stage\tvalue\tunit\tnote\n"
                "HMM hit rows\t6769772\thit rows\tAll family-model hit records\n"
                "Strict tier1 core candidates\t74339\tsequences\tFour core candidate families\n",
                encoding="utf-8",
            )
            figures.render_figure_1(data_dir, root / "out")
            svg = (root / "out" / "figure_1_workflow_funnel.svg").read_text(encoding="utf-8")
            self.assertIn("Genome-level coverage", svg)
            self.assertIn("Protein-level screening", svg)
            self.assertIn("199,923", svg)

    def test_renderer_exports_editable_svg_and_all_requested_formats(self):
        figures = load_module("plot_scheme_a_figures")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data_dir = root / "source_data"
            data_dir.mkdir()
            (data_dir / "figure2_core_scale.tsv").write_text(
                "family\tsequences\tgenomes\n"
                "ePhaZ\t38692\t27839\n"
                "iPhaZ\t32926\t25920\n"
                "OH\t1429\t1410\n"
                "ArchPhaZ_hydrolase\t1292\t1236\n",
                encoding="utf-8",
            )
            (data_dir / "figure2_union.tsv").write_text(
                "core_family_union_genomes\tgtdb_representatives\tcoverage_fraction\n"
                "44814\t199923\t0.22416\n",
                encoding="utf-8",
            )
            outputs = figures.render_figure_2(data_dir, root / "out")
            self.assertEqual({path.suffix for path in outputs}, {".svg", ".pdf", ".tiff", ".png"})
            svg = (root / "out" / "figure_2_core_scale.svg").read_text(encoding="utf-8")
            self.assertIn("<text", svg)
            self.assertIn("44,814", svg)


if __name__ == "__main__":
    unittest.main()
