import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_ephaz_iphaZ_competition.py"


class CompetitionAuditTests(unittest.TestCase):
    def _load(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT}")
        spec = importlib.util.spec_from_file_location("audit_ephaz_iphaZ_competition", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module

    def test_select_other_predictions_and_assign_competition_winner(self):
        module = self._load()
        predictions = {"A": "OTHER", "B": "SP", "C": "OTHER"}
        self.assertEqual(module.select_other_predictions(predictions), {"A", "C"})
        scores = {
            "A": {"ephaz": {"bitscore": 20.0, "evalue": 1e-4}, "iphaz": {"bitscore": 100.0, "evalue": 1e-30}},
            "C": {"ephaz": {"bitscore": 80.0, "evalue": 1e-20}, "iphaz": {"bitscore": 82.0, "evalue": 1e-21}},
            "D": {"ephaz": {}, "iphaz": {}},
        }
        result = module.classify_competition(scores, margin_bits=10.0)
        self.assertEqual(result["A"], "iPhaZ_like")
        self.assertEqual(result["C"], "ambiguous")
        self.assertEqual(result["D"], "no_reportable_hit")

    def test_parse_tblout_keeps_best_hit(self):
        module = self._load()
        lines = [
            "A - model - 1e-5 40.0 0.0 1e-5 40.0 0.0",
            "A - model - 1e-20 100.0 0.0 1e-20 100.0 0.0",
        ]
        hits = module.parse_tblout(lines)
        self.assertEqual(hits["A"]["bitscore"], 100.0)
        self.assertEqual(hits["A"]["evalue"], 1e-20)

    def test_preserves_pipe_containing_candidate_ids(self):
        module = self._load()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "tier1.faa"
            predictions = root / "signalp.txt"
            fasta.write_text(
                ">GCA_1|contig_1\nMMMM\n>GCA_1|contig_2\nMMMMM\n",
                encoding="utf-8",
            )
            predictions.write_text(
                "# ID\tPrediction\nGCA_1|contig_1\tOTHER\nGCA_1|contig_2\tOTHER\n",
                encoding="utf-8",
            )
            self.assertEqual(set(module._read_fasta(fasta)), {"GCA_1|contig_1", "GCA_1|contig_2"})
            self.assertEqual(module.select_other_predictions(module._read_predictions(predictions)), {"GCA_1|contig_1", "GCA_1|contig_2"})

    def test_parse_domtblout_keeps_best_domain_coordinates_and_coverages(self):
        module = self._load()
        line = "A - 100 hmm - 50 1e-20 50.0 0.0 1 1 1e-20 1e-20 50.0 0.0 5 45 10 90 10 90 0.99 desc"
        result = module.parse_domtblout([line], {"A": 100})
        self.assertEqual(result["A"]["hmm_from"], 5)
        self.assertEqual(result["A"]["hmm_to"], 45)
        self.assertAlmostEqual(result["A"]["hmm_coverage"], 0.82)
        self.assertAlmostEqual(result["A"]["target_coverage"], 0.81)


if __name__ == "__main__":
    unittest.main()
