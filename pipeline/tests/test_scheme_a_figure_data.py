"""Unit tests for deterministic Scheme A figure source-data preparation."""
import importlib.util
from pathlib import Path
import unittest

import pandas as pd


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"


def load_module(name: str):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FigureDataTests(unittest.TestCase):
    def test_core_union_excludes_patatin_and_counts_unique_genomes(self):
        figure_data = load_module("build_scheme_a_figure_data")
        tier = pd.DataFrame([
            {"genome": "G1", "family": "ePhaZ"},
            {"genome": "G1", "family": "iPhaZ"},
            {"genome": "G2", "family": "OH"},
            {"genome": "G3", "family": "ArchPhaZ_patatin"},
        ])
        self.assertEqual(figure_data.core_union(tier), 2)

    def test_neighborhood_rate_uses_unique_hit_loci_as_denominator(self):
        figure_data = load_module("build_scheme_a_figure_data")
        audit = pd.DataFrame([
            {"genome": "G1", "locus": "L1", "family": "ePhaZ", "status": "analyzed"},
            {"genome": "G1", "locus": "L2", "family": "ePhaZ", "status": "analyzed"},
        ])
        context = pd.DataFrame([
            {"genome": "G1", "hit_locus": "L1", "hit_family": "ePhaZ", "marker_family": "BdhA"},
            {"genome": "G1", "hit_locus": "L1", "hit_family": "ePhaZ", "marker_family": "BdhA"},
        ])
        result = figure_data.neighborhood_rates(context, audit)
        row = result.iloc[0]
        self.assertEqual(row["candidate_loci"], 2)
        self.assertEqual(row["supported_loci"], 1)
        self.assertEqual(row["support_rate"], 0.5)

    def test_top_signal_phyla_sorts_genome_counts_numerically(self):
        figure_data = load_module("build_scheme_a_figure_data")
        phyla = pd.DataFrame([
            {"phylum": "Low", "genomes": "9"},
            {"phylum": "High", "genomes": "100"},
            {"phylum": "Middle", "genomes": "20"},
        ])
        self.assertEqual(
            list(figure_data.top_signal_phyla(phyla, limit=2)["phylum"]),
            ["High", "Middle"],
        )

    def test_core_phylum_totals_count_each_genome_once(self):
        figure_data = load_module("build_scheme_a_figure_data")
        tier = pd.DataFrame([
            {"genome": "G1", "family": "ePhaZ", "phylum": "P1"},
            {"genome": "G1", "family": "iPhaZ", "phylum": "P1"},
            {"genome": "G2", "family": "OH", "phylum": "P1"},
            {"genome": "G3", "family": "ArchPhaZ_patatin", "phylum": "P1"},
            {"genome": "G4", "family": "ArchPhaZ_hydrolase", "phylum": "P2"},
        ])
        totals = figure_data.core_phylum_totals(tier)
        self.assertEqual(dict(zip(totals["phylum"], totals["genomes"])), {"P1": 2, "P2": 1})

    def test_figure1_genome_coverage_keeps_genomes_separate_from_protein_counts(self):
        figure_data = load_module("build_scheme_a_figure_data")
        tier = pd.DataFrame([
            {"genome": "G1", "family": "ePhaZ"},
            {"genome": "G1", "family": "iPhaZ"},
            {"genome": "G2", "family": "OH"},
            {"genome": "G3", "family": "ArchPhaZ_patatin"},
        ])
        coverage = figure_data.figure1_genome_coverage(tier, total_genomes=10)
        self.assertEqual(list(coverage["value"]), [10, 2])
        self.assertTrue((coverage["unit"] == "genomes").all())


if __name__ == "__main__":
    unittest.main()
