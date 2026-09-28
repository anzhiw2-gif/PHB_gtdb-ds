"""Tests for the nature-style gRodon growth-comparison figure (plan Task 12).

The figure is the public face of the genus-level analysis, so its *vocabulary*
is part of the contract: group labels must be "candidate gene carrier" and
"candidate not detected under the defined search", the intracellular /
extracellular comparison must be labelled secondary, and the historical
phenotype labels must never appear in any rendered axis, legend, title or
caption.

The pre-existing helper tests stay as they are; the new tests build a synthetic
six-panel figure in a temporary directory so they run offline and never touch
``runs/`` or ``results/``.
"""

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "plot_phaded_grodon_growth_nature.py"

FORBIDDEN = ("degrader", "non-degrader", "non degrader", "non_degrader")

CARRIER = "candidate_gene_carrier"
CONTROL = "candidate_not_detected_under_defined_search"


def load_module():
    spec = importlib.util.spec_from_file_location("plot_phaded_grodon_growth_nature", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_inputs(module, n_genera=6, per_genus=3):
    """Synthetic tables in the v2 column shape (no real run data)."""
    rng = np.random.default_rng(20260928)
    genomes = []
    deltas = []
    subtypes = [
        "extracellular dPHASCL type 1",
        "intracellular nPHASCL without lipase box",
    ]
    for index in range(n_genera):
        genus = f"Genus{index:02d}"
        delta = 0.01 * (index - n_genera / 2)
        deltas.append(delta)
        for arm, status, group in (
            ("carrier", CARRIER, "intracellular" if index % 2 else "extracellular"),
            ("control", CONTROL, "control"),
        ):
            for replicate in range(per_genus):
                growth = 0.5 + delta / 2 + 0.02 * replicate
                growth += 0.001 * rng.random()
                if arm == "control":
                    growth -= delta / 2 + delta / 2
                genomes.append(
                    {
                        "genome_id": f"{genus}_{arm}_{replicate}",
                        "candidate_detection_status": status,
                        "group": group,
                        "major_subtype": subtypes[index % 2],
                        "phylum": "Pseudomonadota",
                        "genus": genus,
                        "species": "sp",
                        "growth_rate_per_h": growth,
                        "doubling_time_h": float(np.log(2.0) / growth),
                    }
                )
    balanced = pd.DataFrame(genomes)

    summary_rows = []
    effect_rows = []
    for index in range(n_genera):
        genus = f"Genus{index:02d}"
        delta = 0.01 * (index - n_genera / 2)
        summary_rows.append(
            {
                "genus": genus,
                "n_candidate_carrier": per_genus,
                "n_candidate_not_detected": per_genus,
                "mean_growth_candidate_carrier": 0.5 + delta / 2,
                "mean_growth_candidate_not_detected": 0.5 - delta / 2,
                "delta_carrier_minus_control": delta,
                "median_growth_candidate_carrier": 0.5 + delta / 2,
                "median_growth_candidate_not_detected": 0.5 - delta / 2,
            }
        )
        effect_rows.append(
            {
                "subtype": subtypes[index % 2],
                "genus": genus,
                "delta_carrier_minus_control": delta,
            }
        )
    summary = pd.DataFrame(summary_rows)
    sf_effects = pd.DataFrame(effect_rows)
    sf_tests = (
        sf_effects.groupby("subtype", sort=True)["delta_carrier_minus_control"]
        .mean()
        .reset_index(name="mean_delta")
    )
    sf_tests["wilcoxon_q_bh"] = 0.42

    tests = pd.DataFrame(
        [
            {"metric": "n_balanced_candidate_carrier", "value": n_genera * per_genus},
            {"metric": "n_balanced_candidate_not_detected", "value": n_genera * per_genus},
            {"metric": "n_balanced_genera", "value": n_genera},
            {"metric": "genus_mean_delta_growth_rate_per_h", "value": float(np.mean(deltas))},
            {"metric": "genus_median_delta_growth_rate_per_h", "value": float(np.median(deltas))},
            {"metric": "wilcoxon_signed_rank_two_sided_p", "value": 0.031},
            {"metric": "stratified_permutation_unweighted_p", "value": 0.004},
            {"metric": "exact_sign_test_two_sided_p", "value": 0.115},
            {"metric": "wilcoxon_effect_r", "value": 0.084},
        ]
    )

    group_effect_rows = []
    for index in range(n_genera):
        group_effect_rows.append(
            {
                "group": "intracellular" if index % 2 == 0 else "extracellular",
                "genus": f"Genus{index:02d}",
                "delta_carrier_minus_control": 0.01 * (index - n_genera / 2),
            }
        )
    grp_effects = pd.DataFrame(group_effect_rows)
    grp_tests = pd.DataFrame(
        [
            {
                "group": "intracellular_vs_extracellular",
                "n_intra_genera": int((grp_effects["group"] == "intracellular").sum()),
                "n_extra_genera": int((grp_effects["group"] == "extracellular").sum()),
                "mean_intra_delta": 0.0,
                "mean_extra_delta": 0.0,
                "mean_difference_intra_minus_extra": -0.0135,
                "mannwhitney_u_p": 0.92,
                "permutation_p": 0.48,
                "p_value": 0.48,
                "bootstrap_ci_low": -0.05,
                "bootstrap_ci_high": 0.03,
            }
        ]
    )

    manifest_stats = {
        "manifest_rows": n_genera * per_genus * 2,
        "manifest_carrier": n_genera * per_genus,
        "manifest_control": n_genera * per_genus,
        "manifest_genera": n_genera,
        "carrier_genomes_input": n_genera * per_genus * 4,
        "carrier_genomes_without_a_same_genus_control": n_genera,
    }
    dedup_stats = {"ok": n_genera * per_genus * 2, "failed": 1, "unique_genomes": n_genera * per_genus * 2 + 1}
    return module.FigureInputs(
        balanced=balanced,
        summary=summary,
        tests=tests,
        sf_tests=sf_tests,
        sf_effects=sf_effects,
        grp_tests=grp_tests,
        grp_effects=grp_effects,
        manifest_stats=manifest_stats,
        dedup_stats=dedup_stats,
    )


class NatureGrowthFigureTests(unittest.TestCase):
    def test_cliff_delta_bounds_and_sign(self):
        module = load_module()
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([0.0, 0.0, 0.0])
        # every a > every b -> delta = +1
        self.assertAlmostEqual(module.cliff_delta(a, b), 1.0)
        self.assertAlmostEqual(module.cliff_delta(b, a), -1.0)
        self.assertAlmostEqual(module.cliff_delta(a, a), 0.0)

    def test_cliff_delta_empty_is_nan(self):
        module = load_module()
        self.assertTrue(np.isnan(module.cliff_delta(np.array([]), np.array([1.0]))))

    def test_bootstrap_mean_ci_brackets_mean(self):
        module = load_module()
        rng = np.random.default_rng(7)
        values = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        lo, hi = module.bootstrap_mean_ci(values, rng, n_boot=400)
        self.assertLessEqual(lo, values.mean())
        self.assertLessEqual(values.mean(), hi)
        self.assertLess(lo, hi)

    def test_genus_bootstrap_ci_brackets_group_difference(self):
        module = load_module()
        rng = np.random.default_rng(11)
        a = np.full(8, 0.20)
        b = np.full(8, 0.30)
        lo, hi = module.genus_bootstrap_ci(a, b, rng, n_boot=200)
        self.assertAlmostEqual(lo, 0.10, places=6)
        self.assertAlmostEqual(hi, 0.10, places=6)

    def test_metric_lookup_raises_on_missing(self):
        module = load_module()
        import pandas as pd

        tests = pd.DataFrame({"metric": ["n_balanced_genera"], "value": [1216]})
        self.assertEqual(module.metric(tests, "n_balanced_genera"), 1216.0)
        with self.assertRaises(KeyError):
            module.metric(tests, "not_a_metric")


class FigureVocabularyTests(unittest.TestCase):
    """Requirement: candidate-carrier wording, secondary labels, no phenotype."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.module = load_module()
        self.module.PANEL_LETTERS.clear()
        self.inputs = synthetic_inputs(self.module)
        self.fig = self.module.build_figure(
            self.inputs,
            rng=np.random.default_rng(20260928),
            out_dir=self.tmp / "figures",
            source_dir=self.tmp / "source_data",
            write_audits=True,
        )
        self.texts = self.module.collect_figure_text(self.fig)
        self.blob = "\n".join(self.texts).lower()
        self.addCleanup(self.module.plt.close, self.fig)

    def test_forbidden_legacy_labels_do_not_appear_in_rendered_text(self):
        call = lambda: [
            token for token in FORBIDDEN if token in self.blob
        ]
        self.assertEqual(call(), [])

    def test_label_audit_helper_reports_pass(self):
        report = self.module.audit_figure_labels(self.fig)
        self.assertEqual(report["verdict"], "PASS")
        self.assertEqual(report["forbidden_labels_found"], [])
        self.assertEqual(report["required_vocabulary_missing"], [])
        self.assertGreater(report["texts_measured"], 20)
        self.module.assert_figure_labels(self.fig)

    def test_label_audit_helper_detects_a_forbidden_label(self):
        module = load_module()
        fig = module.plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(f"{FORBIDDEN[1]} mean")
        with self.assertRaisesRegex(ValueError, "forbidden legacy phenotype label"):
            module.assert_figure_labels(fig)
        module.plt.close(fig)

    def test_required_candidate_carrier_vocabulary_is_rendered(self):
        self.assertIn("candidate carrier", self.blob)
        self.assertIn("candidate not detected", self.blob)
        self.assertIn("no candidate detected under the defined search", self.blob)

    def test_compartment_and_per_genome_analyses_are_labelled_secondary(self):
        self.assertIn("secondary", self.blob)
        self.assertIn("intracellular vs extracellular (secondary)", self.blob)
        self.assertIn("genome-level distributions (secondary)", self.blob)

    def test_primary_estimand_is_rendered_as_the_genus_level_mean_difference(self):
        self.assertIn("genus-level mean difference (primary estimand)", self.blob)
        self.assertIn("resampling unit: genus", self.blob)

    def test_null_wording_is_the_bounded_form(self):
        self.assertIn("no difference detected under the current design", self.blob)

    def test_growth_rate_axis_units_are_per_hour(self):
        self.assertIn("(1/h)", self.blob)

    def test_group_label_helper_maps_statuses_and_rejects_unknown(self):
        module = load_module()
        self.assertEqual(module.group_label(CARRIER), "Candidate carrier")
        self.assertEqual(module.group_label(CONTROL), "Candidate not detected")
        with self.assertRaises(ValueError):
            module.group_label("degrader")

    def test_secondary_label_helper_marks_secondary_text(self):
        module = load_module()
        self.assertEqual(module.secondary_label("Foo"), "Foo (secondary)")

    def test_audits_and_source_data_are_written_to_the_given_directories(self):
        collisions = json_load(self.tmp / "figures" / "growth_comparison_nature.collision-audit.json")
        self.assertIn("verdict", collisions)
        for name in (
            "source_per_genome.tsv",
            "source_genus_deltas.tsv",
            "source_main_tests.tsv",
            "source_superfamily_tests.tsv",
            "source_group_tests.tsv",
            "source_ledger.tsv",
        ):
            self.assertTrue((self.tmp / "source_data" / name).is_file(), name)
        per_genome = pd.read_csv(self.tmp / "source_data" / "source_per_genome.tsv", sep="\t")
        self.assertIn("candidate_detection_status", per_genome.columns)
        self.assertNotIn("phaZ_status", per_genome.columns)


class GrowthRateUnitVerificationTests(unittest.TestCase):
    def test_consistent_table_passes(self):
        module = load_module()
        balanced = pd.DataFrame(
            {"doubling_time_h": [1.0, 2.0], "growth_rate_per_h": [np.log(2.0), np.log(2.0) / 2.0]}
        )
        report = module.verify_growth_rate_unit(balanced)
        self.assertEqual(report["checked_rows"], 2)
        self.assertEqual(report["unit"], "1/h")
        self.assertLess(report["max_relative_error"], 1e-9)

    def test_double_converted_rate_is_rejected(self):
        module = load_module()
        rate = np.log(2.0) / 2.0
        balanced = pd.DataFrame(
            {"doubling_time_h": [2.0], "growth_rate_per_h": [np.log(2.0) / rate]}
        )
        with self.assertRaisesRegex(ValueError, "ln2/doubling_time_h"):
            module.verify_growth_rate_unit(balanced)

    def test_wrong_unit_is_rejected(self):
        module = load_module()
        balanced = pd.DataFrame(
            {"doubling_time_h": [24.0], "growth_rate_per_h": [np.log(2.0) / 24.0 * 24.0]}
        )
        with self.assertRaises(ValueError):
            module.verify_growth_rate_unit(balanced)

    def test_missing_rate_column_is_rejected(self):
        module = load_module()
        with self.assertRaisesRegex(ValueError, "growth_rate_per_h"):
            module.verify_growth_rate_unit(pd.DataFrame({"doubling_time_h": [1.0]}))


def json_load(path):
    import json

    return json.loads(Path(path).read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
