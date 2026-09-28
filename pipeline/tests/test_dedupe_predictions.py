import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "deploy" / "20260920_phaded_grodon_growth_01" / "scripts" / "dedupe_predictions.py"


def load_module():
    spec = importlib.util.spec_from_file_location("dedupe_predictions", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class DedupePredictionsTests(unittest.TestCase):
    def test_collapses_duplicates_and_prefers_ok(self):
        module = load_module()
        rows = [
            {"genome_id": "g1", "status": "ok", "growth_rate_per_h": "1.0"},
            {"genome_id": "g2", "status": "failed", "growth_rate_per_h": ""},
            {"genome_id": "g2", "status": "failed", "growth_rate_per_h": ""},
            {"genome_id": "g3", "status": "failed", "growth_rate_per_h": ""},
            {"genome_id": "g3", "status": "ok", "growth_rate_per_h": "0.5"},
            {"genome_id": "g1", "status": "ok", "growth_rate_per_h": "1.0"},
        ]
        out, stats = module.dedupe(rows)
        self.assertEqual([r["genome_id"] for r in out], ["g1", "g2", "g3"])
        self.assertEqual(stats["input_rows"], 6)
        self.assertEqual(stats["unique_genomes"], 3)
        self.assertEqual(stats["duplicates_removed"], 3)
        self.assertEqual(stats["ok"], 2)
        self.assertEqual(stats["failed"], 1)
        # order preserved, ok preferred for g3
        self.assertEqual([r["status"] for r in out], ["ok", "failed", "ok"])
        self.assertEqual(out[2]["growth_rate_per_h"], "0.5")

    def test_no_duplicates_is_identity(self):
        module = load_module()
        rows = [
            {"genome_id": "g1", "status": "ok", "growth_rate_per_h": "1.0"},
            {"genome_id": "g2", "status": "ok", "growth_rate_per_h": "2.0"},
        ]
        out, stats = module.dedupe(rows)
        self.assertEqual(len(out), 2)
        self.assertEqual(stats["duplicates_removed"], 0)


if __name__ == "__main__":
    unittest.main()
