"""Regression tests for domtblout aggregation."""
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "06b_aggregate_hits.py"
SPEC = importlib.util.spec_from_file_location("aggregate", SCRIPT)
aggregate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(aggregate)


def dom_line(target, hmm_from, hmm_to, model_length=100):
    fields = [
        target, "-", "300", "OH", "-", str(model_length), "1e-20",
        "100", "0", "1", "2", "1e-20", "1e-20", "90", "0",
        str(hmm_from), str(hmm_to), "1", "50", "1", "50", "0.99", "desc",
    ]
    return " ".join(fields) + "\n"


class AggregateCoverageTests(unittest.TestCase):
    def test_read_dom_cov_uses_union_of_non_overlapping_hmm_intervals(self):
        with tempfile.TemporaryDirectory() as tmp:
            dom = Path(tmp) / "OH__shard_0001.dom"
            dom.write_text(
                dom_line("protein_1", 1, 30)
                + dom_line("protein_1", 51, 80),
                encoding="utf-8",
            )
            coverage = aggregate.read_dom_cov(str(dom))
        self.assertAlmostEqual(coverage["protein_1"], 0.60)


if __name__ == "__main__":
    unittest.main()
