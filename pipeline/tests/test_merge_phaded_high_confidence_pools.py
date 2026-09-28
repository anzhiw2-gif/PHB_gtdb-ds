import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "merge_phaded_high_confidence_pools.py"


def load_module():
    spec = importlib.util.spec_from_file_location("merge_phaded_high_confidence_pools", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COLS = ["accession", "genome", "superfamily", "nucleophile_type", "high_confidence"]


class MergeHighConfidencePoolsTests(unittest.TestCase):
    def test_merge_adds_pool_origin_after_accession_and_preserves_fields(self):
        module = load_module()
        pool_in = [
            {"accession": "a1", "genome": "g1", "superfamily": "Cys", "nucleophile_type": "cys", "high_confidence": "true"},
        ]
        pool_ext = [
            {"accession": "b1", "genome": "g2", "superfamily": "type2", "nucleophile_type": "ser", "high_confidence": "true"},
        ]
        out_cols, merged = module.merge(pool_in, pool_ext, COLS)
        self.assertEqual(out_cols, ["accession", "pool_origin", "genome", "superfamily", "nucleophile_type", "high_confidence"])
        self.assertEqual(len(merged), 2)
        by_acc = {r["accession"]: r for r in merged}
        self.assertEqual(by_acc["a1"]["pool_origin"], "pool_in")
        self.assertEqual(by_acc["b1"]["pool_origin"], "pool_external")
        self.assertEqual(by_acc["a1"]["superfamily"], "Cys")
        self.assertEqual(by_acc["a1"]["nucleophile_type"], "cys")

    def test_merge_rejects_overlapping_accessions(self):
        module = load_module()
        row = {"accession": "x", "genome": "g", "superfamily": "", "nucleophile_type": "", "high_confidence": ""}
        with self.assertRaises(ValueError):
            module.merge([dict(row)], [dict(row)], COLS)

    def test_merge_sorts_by_pool_origin_then_accession(self):
        module = load_module()
        pool_in = [
            {"accession": "a2", "genome": "g", "superfamily": "", "nucleophile_type": "", "high_confidence": ""},
            {"accession": "a1", "genome": "g", "superfamily": "", "nucleophile_type": "", "high_confidence": ""},
        ]
        pool_ext = [
            {"accession": "b1", "genome": "g", "superfamily": "", "nucleophile_type": "", "high_confidence": ""},
        ]
        _, merged = module.merge(pool_in, pool_ext, COLS)
        self.assertEqual([r["accession"] for r in merged], ["a1", "a2", "b1"])

    def test_merge_requires_accession_as_first_column(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.merge([], [], ["genome", "accession"])


if __name__ == "__main__":
    unittest.main()
