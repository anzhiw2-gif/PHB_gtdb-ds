import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[2] / "deploy" / "20260920_phaded_grodon_growth_01" / "scripts" / "grodon_group_stats.py"
V2_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "grodon_reanalysis_v2.py"


def load_module():
    spec = importlib.util.spec_from_file_location("grodon_group_stats", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_v2_module():
    """Load the v2 module by path, like the other focused test modules do.

    ``importlib`` (rather than ``from grodon_reanalysis_v2 import ...``) keeps
    the loader identical to the dated deploy, which executes the script directly
    instead of importing it as a package module.
    """
    spec = importlib.util.spec_from_file_location("grodon_reanalysis_v2_stats", V2_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V2 = load_v2_module()
CANDIDATE_CARRIER = V2.CANDIDATE_CARRIER
CANDIDATE_NOT_DETECTED = V2.CANDIDATE_NOT_DETECTED
GROWTH_RATE_UNIT = V2.GROWTH_RATE_UNIT


class GrodonGroupStatsTests(unittest.TestCase):
    def test_bh_adjust_matches_manual_bh(self):
        module = load_module()
        rows = [{"p": "0.01"}, {"p": "0.04"}, {"p": "0.03"}, {"p": "0.9"}]
        module.bh_adjust(rows, "p", "q")
        # sorted p: 0.01(rank1), 0.03(rank2), 0.04(rank3), 0.9(rank4)
        # raw q:     0.04,       0.06,       0.05333,   0.9
        # monotone:  0.04,       0.05333,    0.05333,   0.9
        qs = {i: r["q"] for i, r in enumerate(rows)}
        self.assertAlmostEqual(qs[0], 0.04)
        self.assertAlmostEqual(qs[2], 0.053333, places=5)
        self.assertAlmostEqual(qs[1], 0.053333, places=5)
        self.assertAlmostEqual(qs[3], 0.9)

    def test_paired_stats_basic(self):
        module = load_module()
        stats_row = module.paired_stats([0.1, -0.2, 0.05])
        self.assertEqual(stats_row["n_genera"], 3)
        self.assertEqual(stats_row["positive_genera"], 2)
        self.assertEqual(stats_row["negative_genera"], 1)
        self.assertAlmostEqual(stats_row["mean_delta"], -0.05 / 3)
        self.assertEqual(module.paired_stats([])["n_genera"], 0)

    def test_genus_deltas_balance_within_genus(self):
        module = load_module()
        df = pd.DataFrame(
            {
                "genome_id": ["p1", "p2", "p3", "n1", "n2"],
                "phaZ_status": ["phaZ_positive", "phaZ_positive", "phaZ_positive", "phaZ_negative", "phaZ_negative"],
                "major_subtype": ["intracellular nPHASCL without lipase box"] * 3 + ["none"] * 2,
                "group": ["intracellular"] * 3 + ["control"] * 2,
                "genus": ["G1", "G1", "G2", "G1", "G1"],
                "growth_rate_per_h": [0.5, 0.3, 0.9, 0.2, 0.4],
                "status": ["ok"] * 5,
            }
        )
        mask = df["phaZ_status"].eq("phaZ_positive") & df["group"].eq("intracellular")
        neg = df[df["phaZ_status"].eq("phaZ_negative")]
        deltas = module.genus_deltas(df, mask, neg)
        # genus G1: 2 pos (0.5,0.3), 2 neg (0.2,0.4) -> delta = 0.4-0.3 = 0.1
        self.assertEqual(len(deltas), 1)
        self.assertEqual(deltas[0]["genus"], "G1")
        self.assertEqual(deltas[0]["n_phaZ_positive"], 2)
        self.assertAlmostEqual(deltas[0]["delta_positive_minus_negative"], 0.1)

    def test_between_group_test_outputs(self):
        module = load_module()
        result = module.between_group_test([0.2, 0.1, 0.3], [-0.1, 0.0, -0.2], seed=7, n_perm=200, n_boot=200)
        self.assertEqual(result["n_intra_genera"], 3)
        self.assertEqual(result["n_extra_genera"], 3)
        self.assertAlmostEqual(result["mean_difference_intra_minus_extra"], 0.2 - (-0.1))
        self.assertGreaterEqual(result["permutation_p"], 0.0)
        self.assertLessEqual(result["permutation_p"], 1.0)
        self.assertLess(result["bootstrap_ci_low"], result["bootstrap_ci_high"])

    def test_genus_deltas_accepts_mask_from_filtered_subframe(self):
        """Regression: masks built from a filtered sub-frame must not raise IndexingError."""
        module = load_module()
        df = pd.DataFrame(
            {
                "genome_id": ["p1", "p2", "n1", "n2", "x1"],
                "phaZ_status": ["phaZ_positive", "phaZ_positive", "phaZ_negative", "phaZ_negative", "phaZ_positive"],
                "major_subtype": ["A", "B", "none", "none", "A"],
                "group": ["intracellular", "extracellular", "control", "control", "intracellular"],
                "genus": ["G1", "G1", "G1", "G1", "G2"],
                "growth_rate_per_h": [0.5, 0.3, 0.2, 0.4, 0.9],
                "status": ["ok"] * 5,
            }
        )
        df = df[df["status"].eq("ok")].copy()
        pos_rows = df[df["phaZ_status"].eq("phaZ_positive")]
        neg = df[df["phaZ_status"].eq("phaZ_negative")]
        # mask derived from the filtered sub-frame (shorter index) - must be reindexed internally
        mask = pos_rows["major_subtype"].eq("A")
        deltas = module.genus_deltas(df, mask, neg)
        self.assertEqual(len(deltas), 1)
        self.assertEqual(deltas[0]["genus"], "G1")

    def test_grodon_group_stats_empty(self):
        module = load_module()
        result = module.between_group_test([], [0.1], seed=1, n_perm=10, n_boot=10)
        self.assertEqual(result["n_intra_genera"], 0)
        self.assertEqual(result["mean_difference_intra_minus_extra"], "")


class GrodonGroupStatsV2Tests(unittest.TestCase):
    """Task 12 additions: genus dependence, units and the primary estimand.

    These assertions target the canonical
    ``pipeline/scripts/grodon_reanalysis_v2.py`` module.  Every assertion above
    still targets the frozen 2026-09-20 deploy script, which stays untouched.
    """

    def frame(self):
        """4 genera; G3 has 3x the genomes of the others (cluster sizes differ)."""
        rows = []
        for genus, sizes in (("G1", 1), ("G2", 2), ("G3", 3), ("G4", 4)):
            for index in range(sizes):
                rows.append(
                    {
                        "genome_id": f"{genus}_C{index}",
                        "candidate_detection_status": CANDIDATE_CARRIER,
                        "group": "intracellular",
                        "genus": genus,
                        "growth_rate_per_h": 0.60,
                        "doubling_time_h": math.log(2.0) / 0.60,
                    }
                )
                rows.append(
                    {
                        "genome_id": f"{genus}_N{index}",
                        "candidate_detection_status": CANDIDATE_NOT_DETECTED,
                        "group": "control",
                        "genus": genus,
                        "growth_rate_per_h": 0.50,
                        "doubling_time_h": math.log(2.0) / 0.50,
                    }
                )
        return pd.DataFrame(rows)

    def test_v2_genus_resampling_unit_is_the_genus_not_the_genome(self):
        module = V2
        rng = np.random.default_rng(20260928)
        bootstrap = module.cluster_bootstrap_mean_difference(self.frame(), rng, n_boot=25)
        self.assertEqual(bootstrap["resampling_unit"], "genus")
        self.assertEqual(bootstrap["n_genera"], 4)
        self.assertEqual(bootstrap["resample_size"], 4)
        self.assertEqual(len(bootstrap["resampled_units"]), 25 * 4)
        self.assertTrue(set(bootstrap["resampled_units"]) <= {"G1", "G2", "G3", "G4"})
        # 4 genera, not the 20 genomes of the table.
        self.assertNotEqual(bootstrap["resample_size"], len(self.frame()))

    def test_v2_permutation_resamples_whole_genera(self):
        module = V2
        rng = np.random.default_rng(7)
        result = module.genus_level_permutation_test(self.frame(), rng, n_permutations=50)
        self.assertEqual(result["resampling_unit"], "genus")
        self.assertEqual(result["n_genera"], 4)
        self.assertEqual(len(result["resampled_units"]), 50 * 4)
        self.assertTrue(set(result["resampled_units"]) <= {"G1", "G2", "G3", "G4"})

    def test_v2_genus_mean_difference_is_the_primary_estimand(self):
        module = V2
        estimates = module.genus_level_estimates(
            self.frame(), np.random.default_rng(20260928), n_boot=200, n_permutations=200
        )
        self.assertEqual(estimates["primary_estimand"], "genus_level_mean_difference")
        self.assertEqual(estimates["resampling_unit"], "genus")
        self.assertAlmostEqual(estimates["mean_delta"], 0.10, places=9)
        self.assertEqual(estimates["effect_unit"], GROWTH_RATE_UNIT)
        self.assertIn("genus_median_difference", estimates["secondary_analyses"])
        self.assertIn("genus_sign_test", estimates["secondary_analyses"])
        self.assertIn("genome_level_mann_whitney_u", estimates["secondary_analyses"])
        self.assertIn("intracellular_vs_extracellular", estimates["secondary_analyses"])

    def test_v2_null_wording_is_bounded_without_a_preregistered_margin(self):
        module = V2
        outcome = module.equivalence_statement(
            (-0.01, 0.02), mean_delta=0.005, equivalence_margin=None
        )
        self.assertEqual(outcome["basis"], "not_detected")
        self.assertEqual(
            outcome["statement"], "no difference detected under the current design"
        )

    def test_v2_growth_rate_is_ln2_over_doubling_time_in_per_hour(self):
        module = V2
        self.assertEqual(GROWTH_RATE_UNIT, "1/h")
        self.assertAlmostEqual(module.growth_rate_per_h(1.0), math.log(2.0))
        self.assertAlmostEqual(
            module.resolve_growth_rate({"doubling_time_h": 2.0}), math.log(2.0) / 2.0
        )
        # A table that already carries a rate is not converted a second time.
        self.assertAlmostEqual(
            module.resolve_growth_rate({"growth_rate_per_h": 0.25}), 0.25
        )


if __name__ == "__main__":
    unittest.main()
