"""Tests for the P5 with-lipase deferred-pool stratified sampler.

Synthetic fixtures only; the real 1,206,655-row pool is streamed by the run.
"""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import sample_phaded_deferred_pool as module  # noqa: E402


FIELDS = list(module.COLUMNS)


def write_pool(path: Path, rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})
    return path


def pool_row(protein_id: str, *, model: str = "DED_hfam_2", evalue: str = "1e-60",
             override: str = "false", genome: str = "G1") -> dict[str, str]:
    return {
        "protein_id": protein_id, "genome": genome, "superfamily": "intracellular nPHASCL with lipase box",
        "discovery_best_family": "family_DED_hfam_2_x", "discovery_families_hit": model,
        "discovery_best_evalue": evalue, "trained_best_model": model,
        "trained_best_evalue": evalue, "override_by_trained": override,
        "model_layer": "discovery_hmm_uncalibrated",
    }


class EvalueBucketTests(unittest.TestCase):
    def test_bucket_edges_follow_the_frozen_ranges(self):
        self.assertEqual(module.evalue_bucket("1e-60"), "lt_1e-50")
        self.assertEqual(module.evalue_bucket("1e-40"), "1e-50_to_1e-30")
        self.assertEqual(module.evalue_bucket("1e-20"), "1e-30_to_1e-10")
        self.assertEqual(module.evalue_bucket("1e-7"), "1e-10_to_1e-5")
        self.assertEqual(module.evalue_bucket("1e-3"), "gte_1e-5")

    def test_unparsable_and_negative_are_reported_not_guessed(self):
        self.assertEqual(module.evalue_bucket(""), "unparsable")
        self.assertEqual(module.evalue_bucket("abc"), "unparsable")
        self.assertEqual(module.evalue_bucket("-1"), "unparsable")


class StrataTests(unittest.TestCase):
    def test_stratum_uses_only_columns_the_pool_actually_has(self):
        key = module.strata_of(pool_row("p1", model="DED_hfam_70", evalue="1e-60"))
        self.assertEqual(key, ("DED_hfam_70", "lt_1e-50", "false"))

    def test_missing_trained_claim_is_not_blank(self):
        key = module.strata_of(pool_row("p1", model="", evalue="1e-3"))
        self.assertEqual(key[0], "no_trained_claim")

    def test_missing_columns_raise(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.tsv"
            path.write_text("protein_id\tgenome\np1\tG1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing columns"):
                list(module.read_pool(path))


class AllocationTests(unittest.TestCase):
    def test_proportional_allocation_and_every_stratum_represented(self):
        counts = {("a", "lt_1e-50", "false"): 900, ("b", "gte_1e-5", "false"): 100}
        allocation = module.allocate(counts, 1000, 100)
        self.assertEqual(sum(allocation.values()), 100)
        self.assertEqual(allocation[("a", "lt_1e-50", "false")], 90)
        self.assertEqual(allocation[("b", "gte_1e-5", "false")], 10)

    def test_tiny_stratum_still_gets_one_draw(self):
        counts = {("a", "lt_1e-50", "false"): 9999, ("b", "gte_1e-5", "false"): 1}
        allocation = module.allocate(counts, 10000, 100)
        self.assertEqual(sum(allocation.values()), 100)
        self.assertGreaterEqual(allocation[("b", "gte_1e-5", "false")], 1)

    def test_sample_below_stratum_count_is_refused(self):
        counts = {("a", "lt_1e-50", "false"): 5, ("b", "gte_1e-5", "false"): 5}
        with self.assertRaisesRegex(ValueError, "below the stratum count"):
            module.allocate(counts, 10, 1)


class SamplingTests(unittest.TestCase):
    def _fixture(self, tmp: Path):
        rows = (
            [pool_row(f"a{i}", model="DED_hfam_2", evalue="1e-60") for i in range(50)]
            + [pool_row(f"b{i}", model="DED_hfam_70", evalue="1e-3", override="true") for i in range(50)]
        )
        return write_pool(tmp / "pool.tsv", rows)

    def test_same_seed_reproduces_the_same_sample(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            pool = self._fixture(tmp)
            counts, total = module.collect_strata(pool)
            allocation = module.allocate(counts, total, 20)
            first = [row["protein_id"] for row in module.draw_sample(pool, allocation, 7)]
            second = [row["protein_id"] for row in module.draw_sample(pool, allocation, 7)]
            self.assertEqual(first, second)
            self.assertEqual(len(first), len(set(first)))

    def test_every_stratum_is_present_in_the_sample(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            pool = self._fixture(tmp)
            counts, total = module.collect_strata(pool)
            allocation = module.allocate(counts, total, 20)
            sample = module.draw_sample(pool, allocation, 42)
            seen = {module.strata_of(row) for row in sample}
            self.assertEqual(seen, set(counts))

    def test_cli_writes_artifacts_and_refuses_a_non_empty_dir(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            pool = self._fixture(tmp)
            out = tmp / "out"
            self.assertEqual(
                module.main(["--pool", str(pool), "--out-dir", str(out),
                             "--sample-size", "10", "--expect-pool-size", "100"]), 0,
            )
            summary = json.loads((out / "p5_stratified_sample_summary.json").read_text("utf-8"))
            self.assertEqual(summary["pool_total"], 100)
            self.assertEqual(summary["sample_size"], 10)
            self.assertEqual(summary["seed"], 42)
            self.assertIn("inherited from the frozen pilot", summary["sample_size_basis"])
            self.assertIn("Deletes nothing", summary["boundary"])
            self.assertEqual(len(summary["keys_requested_but_unavailable"]), 3)
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--pool", str(pool), "--out-dir", str(out), "--sample-size", "10"])

    def test_pool_size_mismatch_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            pool = self._fixture(tmp)
            with self.assertRaisesRegex(ValueError, "does not match the declared"):
                module.main(["--pool", str(pool), "--out-dir", str(tmp / "out"),
                             "--sample-size", "10", "--expect-pool-size", "999999"])

    def test_pilot_composition_reader_buckets_the_blastp_table(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = tmp / "pilot.tsv"
            path.write_text("q1\t30.0\t200\t1e-60\t100\nq2\t25.0\t150\t1e-3\t40\n", encoding="utf-8")
            buckets = module.pilot_strata_from_accessions(path)
            self.assertEqual(buckets["lt_1e-50"], 1)
            self.assertEqual(buckets["gte_1e-5"], 1)


if __name__ == "__main__":
    unittest.main()
