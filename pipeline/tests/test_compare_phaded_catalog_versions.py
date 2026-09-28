"""Tests for the PhaDED v1/v2 catalog comparison and preregistered rule scenarios.

Written before the implementation (test-first).  The module under test answers two
questions without touching frozen evidence:

* *impact*: which accessions the v2 candidate catalog adds to, removes from, or
  relabels relative to the frozen v1 high-confidence table;
* *rule sensitivity*: how the candidate composition changes under the three
  preregistered scenarios (``current_v1_rules``,
  ``localization_uncertain_not_excluded``, ``sequence_validated_models_only``).

Load-bearing invariants pinned here:

* ``compare`` returns **sets** for ``removed`` / ``added`` / ``relabelled`` so the
  answer cannot depend on input insertion order;
* the scenario set is preregistered data inside the module: an unknown scenario
  name raises and the message lists the valid names, it is never ignored;
* the frozen v1 rule constants re-encoded by the comparison module are asserted
  equal to the frozen v1 source module (``filter_phaded_high_confidence.py``) so
  the two cannot drift apart;
* a v1 rule that cannot be reconstructed from the supplied columns yields
  ``pending`` for that row, never a guess, while an explicitly blank column stays
  fail-closed (blank evidence never passes a v1 criterion);
* every composition block carries the scenario name, the rule text and the fixed
  caveat string: the 76:24 intracellular:extracellular split is reported only as
  *candidate composition under the named rule* and never as a biological ratio.

All fixtures are synthetic and in-memory / tempfile: no project data is read and
no catalog is recomputed.
"""

from __future__ import annotations

import csv
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import compare_phaded_catalog_versions as module  # noqa: E402
import filter_phaded_high_confidence as frozen_v1  # noqa: E402


# --- synthetic fixture vocabulary (mirrors, never reads, project values) ------

TYPE1 = "extracellular dPHASCL type 1"
TYPE2 = "extracellular dPHASCL type 2"
CYS = "intracellular nPHASCL without lipase box"
NPHAMCL = "intracellular nPHAMCL"
PHAZ7 = "extracellular native-SCL/PhaZ7-like"

MISSING = object()
"""Sentinel: the column is absent from the row mapping (unknown, not blank)."""


def row(accession: str = "A1", **overrides) -> dict[str, str]:
    """Build one synthetic evidence row; each test overrides only its own field.

    Passing ``MISSING`` removes the column entirely, which is what distinguishes
    "cannot be reconstructed" (pending) from "explicitly blank" (fail-closed).
    """
    base: dict[str, str] = {
        "accession": accession,
        "genome": "GCA_000000001.1",
        "superfamily": TYPE1,
        "interpro_status": "interpro_supported",
        "architecture_consistency": "consistent",
        "profile_evidence_status": "profile_trained_hit",
        "profile_best_evalue": "1e-100",
        "is_confounder": "false",
        "signalp_class": "SP",
        "transport_signal_prediction": "SP",
        "lipase_box_state": "supported",
        "lipase_box_x1": "A",
        "catalytic_domain_type": "type1_verified",
        "lid_state": "not_detected_pattern",
        "sbd_pf06850_binding_state": "not_detected",
        "ahsmg_state": "not_detected_pattern",
        "nucleophile_family": "GxSxG",
        "sequence_family_call": "hfam_test",
        "model_layer": "sequence_family_hmm_validated",
    }
    for key, value in overrides.items():
        if value is MISSING:
            base.pop(key, None)
        else:
            base[key] = value
    return base


def cys_row(accession: str, **overrides) -> dict[str, str]:
    """A Cys-superfamily row that satisfies its own v1 architecture criteria."""
    defaults = {
        "superfamily": CYS,
        "signalp_class": "OTHER",
        "transport_signal_prediction": "OTHER",
        "lipase_box_state": "not_detected_pattern",
        "lipase_box_x1": "",
        "catalytic_domain_type": "undetermined_no_oxyanion",
        "sbd_pf06850_binding_state": "detected",
        "nucleophile_family": "hfam_cys",
        "sequence_family_call": "hfam_cys",
    }
    defaults.update(overrides)
    return row(accession, **defaults)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> Path:
    fields: list[str] = []
    for item in rows:
        for key in item:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for item in rows:
            writer.writerow({field: item.get(field, "") for field in fields})
    return path


# --- the required comparison interface ---------------------------------------


class CompareApiTests(unittest.TestCase):
    """``compare`` reports accession-level changes between two label maps."""

    def test_comparison_reports_added_removed_and_relabelled(self):
        result = module.compare(
            v1={"A": "high_confidence", "B": "high_confidence"},
            v2={"B": "function_unresolved", "C": "core_sequence_homolog"},
        )
        self.assertEqual(result["removed"], {"A"})
        self.assertEqual(result["added"], {"C"})
        self.assertEqual(result["relabelled"], {"B"})

    def test_identical_inputs_report_all_three_keys_empty(self):
        labels = {"A": "high_confidence", "B": "function_unresolved"}
        result = module.compare(v1=labels, v2=dict(labels))
        for key in ("removed", "added", "relabelled"):
            self.assertEqual(result[key], set(), key)
        self.assertEqual(result["retained_unchanged"], {"A", "B"})

    def test_v1_only_subset_reports_removals_only(self):
        v1 = {"A": "high_confidence", "B": "high_confidence", "C": "high_confidence"}
        result = module.compare(v1=v1, v2={"A": "high_confidence"})
        self.assertEqual(result["removed"], {"B", "C"})
        self.assertEqual(result["added"], set())
        self.assertEqual(result["relabelled"], set())
        self.assertEqual(result["retained_unchanged"], {"A"})

    def test_v2_only_superset_reports_additions_only(self):
        v2 = {
            "A": "high_confidence",
            "B": "remote_homolog_candidate",
            "C": "function_unresolved",
        }
        result = module.compare(v1={"A": "high_confidence"}, v2=v2)
        self.assertEqual(result["removed"], set())
        self.assertEqual(result["added"], {"B", "C"})
        self.assertEqual(result["relabelled"], set())
        self.assertEqual(result["retained_unchanged"], {"A"})

    def test_result_is_stable_under_input_insertion_order(self):
        v1 = {"A": "high_confidence", "B": "high_confidence", "C": "high_confidence"}
        v2 = {
            "C": "core_sequence_homolog",
            "B": "probable_sequence_homolog",
            "D": "function_unresolved",
        }
        forward = module.compare(v1=v1, v2=v2)
        backward = module.compare(
            v1=dict(reversed(list(v1.items()))),
            v2=dict(reversed(list(v2.items()))),
        )
        self.assertEqual(forward, backward)
        self.assertEqual(
            json.dumps(module.impact_summary(forward), sort_keys=True),
            json.dumps(module.impact_summary(backward), sort_keys=True),
        )

    def test_impact_summary_reports_counts_transitions_and_deltas(self):
        result = module.compare(
            v1={"A": "high_confidence", "B": "high_confidence", "C": "high_confidence"},
            v2={"A": "high_confidence", "B": "function_unresolved", "D": "core_sequence_homolog"},
        )
        summary = module.impact_summary(result)
        self.assertEqual(summary["v1_count"], 3)
        self.assertEqual(summary["v2_count"], 3)
        self.assertEqual(summary["removed_count"], 1)
        self.assertEqual(summary["added_count"], 1)
        self.assertEqual(summary["relabelled_count"], 1)
        self.assertEqual(
            summary["label_transitions"],
            {"high_confidence->function_unresolved": 1},
        )
        self.assertEqual(summary["label_counts_v1"]["high_confidence"], 3)
        self.assertEqual(summary["label_counts_v2"]["high_confidence"], 1)
        self.assertEqual(summary["label_deltas"]["high_confidence"], -2)


# --- preregistered scenario registry -----------------------------------------


class ScenarioRegistryTests(unittest.TestCase):
    """The scenario set is fixed in the module; unknown names raise."""

    def test_preregistered_scenario_names_are_exactly_the_three(self):
        self.assertEqual(
            set(module.SCENARIOS),
            {
                "current_v1_rules",
                "localization_uncertain_not_excluded",
                "sequence_validated_models_only",
            },
        )

    def test_scenarios_are_data_with_rule_text_criteria_and_source(self):
        for name, scenario in module.SCENARIOS.items():
            self.assertEqual(scenario.name, name)
            self.assertTrue(scenario.description.strip(), name)
            self.assertTrue(scenario.rule.strip(), name)
            self.assertTrue(scenario.source.strip(), name)
            self.assertTrue(scenario.criteria, name)
            for criterion in scenario.criteria:
                self.assertIn(criterion, module.CRITERIA, name)
            for criterion in scenario.removed_criteria:
                self.assertIn(criterion, module.CRITERIA, name)
                self.assertNotIn(criterion, scenario.criteria, name)

    def test_unknown_scenario_raises_and_lists_every_valid_name(self):
        with self.assertRaises(ValueError) as caught:
            module.resolve_scenario_names(["not_a_scenario"])
        message = str(caught.exception)
        for name in module.SCENARIOS:
            self.assertIn(name, message)

    def test_unknown_scenario_also_raises_from_the_cli(self):
        with self.assertRaises(ValueError):
            module.resolve_scenario_names(["current_v1_rules", "typo_scenario"])

    def test_localization_scenario_removes_exactly_the_export_requirement(self):
        strict = set(module.SCENARIOS["current_v1_rules"].criteria)
        relaxed = set(module.SCENARIOS["localization_uncertain_not_excluded"].criteria)
        self.assertEqual(strict - relaxed, {"v1_localization_export_requirement"})
        self.assertEqual(relaxed - strict, set())

    def test_only_the_no_export_rule_survives_in_the_relaxed_scenario(self):
        relaxed = set(module.SCENARIOS["localization_uncertain_not_excluded"].criteria)
        self.assertIn("v1_localization_no_export_requirement", relaxed)
        self.assertNotIn("v1_localization_export_requirement", relaxed)

    def test_sequence_validated_scenario_is_the_pure_model_layer_restriction(self):
        scenario = module.SCENARIOS["sequence_validated_models_only"]
        self.assertEqual(scenario.criteria, ("sequence_family_hmm_validated",))
        self.assertIn("sequence_family_hmm_validated", scenario.rule)

    def test_v1_scenarios_cite_the_frozen_source_file(self):
        for name in ("current_v1_rules", "localization_uncertain_not_excluded"):
            self.assertIn("filter_phaded_high_confidence.py", module.SCENARIOS[name].source)

    def test_module_docstring_documents_the_encoded_v1_rules(self):
        doc = module.__doc__ or ""
        self.assertIn("filter_phaded_high_confidence.py", doc)
        self.assertIn("localization", doc.lower())
        self.assertIn("48,038", doc)
        self.assertIn("pending", doc.lower())

    def test_public_api_is_declared_and_present(self):
        self.assertTrue(module.__all__)
        for name in (
            "compare",
            "impact_summary",
            "SCENARIOS",
            "CRITERIA",
            "resolve_scenario_names",
            "evaluate_row",
            "evaluate_scenario",
            "analyze",
            "build_report",
            "composition_blocks",
            "validate_report_guard",
            "main",
        ):
            self.assertIn(name, module.__all__, name)
            self.assertTrue(hasattr(module, name), name)
        for name in module.__all__:
            self.assertTrue(hasattr(module, name), name)


# --- frozen v1 rule parity (drift guard) -------------------------------------


class FrozenV1RuleParityTests(unittest.TestCase):
    """The re-encoded v1 constants must equal the frozen v1 source constants."""

    def test_superfamily_criteria_match_the_frozen_source(self):
        self.assertEqual(dict(module.V1_SUPERFAMILY_CRITERIA), dict(frozen_v1.SUPERFAMILY_CRITERIA))

    def test_transport_and_hydrophobic_sets_match_the_frozen_source(self):
        self.assertEqual(set(module.V1_EXPORT_CLASSES), set(frozen_v1.EXPORT))
        self.assertEqual(set(module.V1_NO_EXPORT_CLASSES), set(frozen_v1.NO_EXPORT))
        self.assertEqual(set(module.V1_HYDROPHOBIC), set(frozen_v1.HYDROPHOBIC))

    def test_export_requirement_partition_matches_the_frozen_source(self):
        for superfamily, (requirement, _criteria) in frozen_v1.SUPERFAMILY_CRITERIA.items():
            if requirement == "export":
                self.assertIn(superfamily, module.V1_EXPORT_REQUIRING)
                self.assertNotIn(superfamily, module.V1_NO_EXPORT_REQUIRED)
            elif requirement == "no_export":
                self.assertIn(superfamily, module.V1_NO_EXPORT_REQUIRED)
                self.assertNotIn(superfamily, module.V1_EXPORT_REQUIRING)
            else:
                self.assertNotIn(superfamily, module.V1_EXPORT_REQUIRING)
                self.assertNotIn(superfamily, module.V1_NO_EXPORT_REQUIRED)

    def test_nucleophile_typing_matches_the_frozen_source(self):
        for superfamily in frozen_v1.SUPERFAMILY_CRITERIA:
            for lipase_box in ("supported", "not_detected_pattern", "pending"):
                for ahsmg in ("supported", "not_detected_pattern", ""):
                    self.assertEqual(
                        module.v1_nucleophile_type(superfamily, lipase_box, ahsmg),
                        frozen_v1.nucleophile_type(superfamily, lipase_box, ahsmg),
                        (superfamily, lipase_box, ahsmg),
                    )

    def test_v1_hold_reason_labels_match_the_frozen_source_labels(self):
        self.assertEqual(
            set(module.V1_HOLD_REASONS),
            {
                "not_in_v1_scope",
                "common_criteria",
                "localization_conflict",
                "nucleophile_conflict",
                "architecture_criteria",
            },
        )
        self.assertLessEqual(
            set(module.V1_HOLD_REASON_BY_CRITERION.values()),
            set(module.V1_HOLD_REASONS),
        )


# --- scenario behaviour -------------------------------------------------------


class ScenarioBehaviourTests(unittest.TestCase):
    """Each scenario selects exactly the rows its preregistered rule describes."""

    def test_type1_without_a_secretion_signal_is_held_only_by_current_v1_rules(self):
        rows = [row("T1_NO_SIGNAL", signalp_class="OTHER", transport_signal_prediction="OTHER")]
        strict = module.evaluate_scenario(rows, "current_v1_rules")
        relaxed = module.evaluate_scenario(rows, "localization_uncertain_not_excluded")
        self.assertEqual(strict.included, ())
        self.assertEqual(strict.excluded, ("T1_NO_SIGNAL",))
        self.assertEqual(strict.excluded_reasons, {"localization_conflict": 1})
        self.assertEqual(relaxed.included, ("T1_NO_SIGNAL",))
        self.assertEqual(relaxed.pending, ())

    def test_the_no_export_side_of_the_localization_gate_is_retained(self):
        rows = [cys_row("CYS_WITH_SIGNAL", signalp_class="SP", transport_signal_prediction="SP")]
        for name in ("current_v1_rules", "localization_uncertain_not_excluded"):
            result = module.evaluate_scenario(rows, name)
            self.assertEqual(result.included, (), name)
            self.assertEqual(result.excluded_reasons, {"localization_conflict": 1}, name)

    def test_sequence_validated_scenario_selects_only_validated_model_layers(self):
        rows = [
            row("VALIDATED", model_layer="sequence_family_hmm_validated"),
            row("DISCOVERY", model_layer="discovery_hmm_uncalibrated"),
            row("REFERENCE_ONLY", model_layer="reference_query_only"),
        ]
        result = module.evaluate_scenario(rows, "sequence_validated_models_only")
        self.assertEqual(result.included, ("VALIDATED",))
        self.assertEqual(len(result.excluded), 2)

    def test_absent_model_layer_column_is_pending_not_excluded(self):
        rows = [row("NO_LAYER", model_layer=MISSING)]
        result = module.evaluate_scenario(rows, "sequence_validated_models_only")
        self.assertEqual(result.included, ())
        self.assertEqual(result.excluded, ())
        self.assertEqual(result.pending, ("NO_LAYER",))
        self.assertEqual(result.pending_criteria, {"sequence_family_hmm_validated": 1})

    def test_absent_evidence_column_is_pending_not_guessed(self):
        rows = [row("NO_PROFILE", profile_evidence_status=MISSING)]
        result = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(result.pending, ("NO_PROFILE",))
        self.assertEqual(result.pending_criteria, {"v1_strong_profile": 1})
        self.assertEqual(result.included, ())

    def test_blank_evidence_column_stays_fail_closed_under_v1(self):
        rows = [row("BLANK_PROFILE", profile_evidence_status="")]
        result = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(result.excluded, ("BLANK_PROFILE",))
        self.assertEqual(result.excluded_reasons, {"common_criteria": 1})

    def test_confounder_membership_is_pending_without_a_confounder_input(self):
        rows = [row("C1", is_confounder=MISSING)]
        unresolved = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(unresolved.pending_criteria, {"v1_not_confounder": 1})
        flagged = module.evaluate_scenario(
            rows,
            "current_v1_rules",
            context=module.ScenarioContext(confounders=frozenset({"C1"})),
        )
        self.assertEqual(flagged.excluded, ("C1",))
        self.assertEqual(flagged.excluded_reasons, {"common_criteria": 1})
        clean = module.evaluate_scenario(
            rows,
            "current_v1_rules",
            context=module.ScenarioContext(confounders=frozenset()),
        )
        self.assertEqual(clean.included, ("C1",))

    def test_confounder_file_overrides_a_row_level_negative_flag(self):
        rows = [row("C1", is_confounder="false")]
        result = module.evaluate_scenario(
            rows,
            "current_v1_rules",
            context=module.ScenarioContext(confounders=frozenset({"C1"})),
        )
        self.assertEqual(result.excluded, ("C1",))
        self.assertEqual(result.first_failing_criteria, {"C1": "v1_not_confounder"})

    def test_common_criteria_failures_use_the_frozen_hold_reason(self):
        rows = [
            row("NO_INTERPRO", interpro_status="interpro_unsupported"),
            row("ARCH_CONFLICT", architecture_consistency="conflicting"),
            row("CONFOUNDER", is_confounder="true"),
            row("WEAK_PROFILE", profile_evidence_status="discovery_hit", profile_best_evalue="1e-3"),
        ]
        result = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(len(result.excluded), 4)
        self.assertEqual(result.excluded_reasons, {"common_criteria": 4})
        self.assertEqual(
            result.first_failing_criteria,
            {
                "ARCH_CONFLICT": "v1_no_architecture_conflict",
                "CONFOUNDER": "v1_not_confounder",
                "NO_INTERPRO": "v1_interpro_supported",
                "WEAK_PROFILE": "v1_strong_profile",
            },
        )

    def test_discovery_hit_needs_the_frozen_evalue_threshold(self):
        rows = [
            row("DISCOVERY_STRONG", profile_evidence_status="discovery_hit", profile_best_evalue="1e-30"),
            row("DISCOVERY_BORDERLINE", profile_evidence_status="discovery_hit", profile_best_evalue="1e-10"),
        ]
        result = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(result.included, ("DISCOVERY_STRONG",))
        self.assertEqual(result.excluded, ("DISCOVERY_BORDERLINE",))

    def test_architecture_criteria_are_applied_per_superfamily(self):
        rows = [
            row("TYPE1_WRONG_GEOMETRY", catalytic_domain_type="type2_verified"),
            row(
                "NPHAMCL_NO_LID",
                superfamily=NPHAMCL,
                signalp_class="OTHER",
                transport_signal_prediction="OTHER",
                lid_state="not_detected",
            ),
            row("PHAZ7_WITHOUT_AHSMG", superfamily=PHAZ7, ahsmg_state="not_detected_pattern"),
            cys_row("CYS_WITH_LIPASE_BOX", lipase_box_state="supported"),
        ]
        result = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(len(result.excluded), 4)
        # a PhaZ7 row without AHSMG and a Cys row carrying a GxSxG are nucleophile
        # conflicts in the frozen v1 typing; the others are architecture failures
        self.assertEqual(
            result.excluded_reasons,
            {"architecture_criteria": 2, "nucleophile_conflict": 2},
        )
        self.assertEqual(
            result.first_failing_criteria["TYPE1_WRONG_GEOMETRY"], "v1_architecture_criteria"
        )

    def test_unrecognized_superfamily_is_out_of_v1_scope_not_pending(self):
        rows = [row("UNASSIGNED", superfamily="not a v1 superfamily")]
        result = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(result.excluded, ("UNASSIGNED",))
        self.assertEqual(result.excluded_reasons, {"not_in_v1_scope": 1})
        self.assertEqual(result.pending, ())

    def test_absent_superfamily_column_is_pending(self):
        rows = [row("NO_SUPERFAMILY", superfamily=MISSING)]
        result = module.evaluate_scenario(rows, "current_v1_rules")
        self.assertEqual(result.pending, ("NO_SUPERFAMILY",))
        self.assertEqual(
            set(result.pending_criteria),
            {
                "v1_superfamily_recognized",
                "v1_localization_export_requirement",
                "v1_localization_no_export_requirement",
                "v1_architecture_criteria",
            },
        )
        self.assertTrue(all(count == 1 for count in result.pending_criteria.values()))

    def test_a_row_without_an_accession_is_rejected(self):
        with self.assertRaises(ValueError):
            module.evaluate_scenario([{"superfamily": TYPE1}], "current_v1_rules")


# --- report shape, determinism, guard ----------------------------------------


class ReportContentTests(unittest.TestCase):
    """Per-scenario totals, composition and accession-level differences."""

    def fixture_rows(self) -> list[dict[str, str]]:
        return [
            row("T1_A"),
            cys_row("CYS_A"),
            row(
                "T1_NO_SIGNAL",
                signalp_class="OTHER",
                transport_signal_prediction="OTHER",
            ),
            row("DISCOVERY_ONLY", model_layer="discovery_hmm_uncalibrated"),
        ]

    def build(self):
        return module.build_report(
            rows=self.fixture_rows(),
            v1_labels={"T1_A": "high_confidence", "CYS_A": "high_confidence"},
            v2_labels={
                "T1_A": "core_sequence_homolog",
                "CYS_A": "core_sequence_homolog",
                "T1_NO_SIGNAL": "probable_sequence_homolog",
            },
        )

    def test_report_sections_are_present(self):
        report = self.build()
        self.assertIn("v1_v2_impact", report)
        self.assertIn("scenario_definitions", report)
        self.assertEqual(set(report["scenarios"]), set(module.SCENARIOS))
        for name, block in report["scenarios"].items():
            self.assertEqual(block["scenario"], name)
            for key in (
                "rule",
                "total_proteins",
                "genomes",
                "pending_proteins",
                "pending_criteria",
                "excluded_reason_counts",
                "compositions",
                "differences_vs_other_scenarios",
            ):
                self.assertIn(key, block, name)
            for key in ("superfamily", "family", "transport_signal", "candidate_localization"):
                self.assertIn(key, block["compositions"], name)

    def test_totals_and_compositions_count_the_selected_rows(self):
        report = self.build()
        strict = report["scenarios"]["current_v1_rules"]
        self.assertEqual(strict["total_proteins"], 3)
        relaxed = report["scenarios"]["localization_uncertain_not_excluded"]
        self.assertEqual(relaxed["total_proteins"], 4)
        validated = report["scenarios"]["sequence_validated_models_only"]
        self.assertEqual(validated["total_proteins"], 3)
        self.assertEqual(
            relaxed["compositions"]["superfamily"]["counts"],
            {TYPE1: 3, CYS: 1},
        )
        self.assertEqual(
            relaxed["compositions"]["transport_signal"]["counts"],
            {"OTHER": 2, "SP": 2},
        )
        self.assertEqual(
            relaxed["compositions"]["family"]["counts"],
            {"hfam_cys": 1, "hfam_test": 3},
        )
        self.assertEqual(
            relaxed["compositions"]["candidate_localization"]["counts"],
            {"extracellular": 3, "intracellular": 1},
        )

    def test_genomes_are_counted_distinctly(self):
        rows = [
            row("A", genome="G1"),
            row("B", genome="G1"),
            row("C", genome="G2", superfamily=TYPE2, catalytic_domain_type="type2_verified"),
        ]
        report = module.build_report(
            rows=rows, v1_labels={}, v2_labels={}
        )
        block = report["scenarios"]["current_v1_rules"]
        self.assertEqual(block["total_proteins"], 3)
        self.assertEqual(block["genomes"], 2)

    def test_pending_rows_are_counted_separately_from_the_totals(self):
        rows = [row("OK"), row("UNKNOWN", profile_evidence_status=MISSING)]
        report = module.build_report(rows=rows, v1_labels={}, v2_labels={})
        block = report["scenarios"]["current_v1_rules"]
        self.assertEqual(block["total_proteins"], 1)
        self.assertEqual(block["pending_proteins"], 1)
        self.assertEqual(block["pending_criteria"], {"v1_strong_profile": 1})
        # a count that is not decidable may never be reported as final
        self.assertEqual(block["total_proteins_status"], "lower_bound_1_rows_pending_not_decidable")
        self.assertEqual(block["genomes_status"], "lower_bound_1_rows_pending_not_decidable")

    def test_decidable_counts_are_marked_final(self):
        report = module.build_report(rows=[row("OK")], v1_labels={}, v2_labels={})
        block = report["scenarios"]["current_v1_rules"]
        self.assertEqual(block["total_proteins_status"], "final")
        self.assertEqual(block["genomes_status"], "final")

    def test_scenario_difference_counts_are_reported_in_the_json(self):
        report = self.build()
        strict = report["scenarios"]["current_v1_rules"]["differences_vs_other_scenarios"]
        relaxed = report["scenarios"]["localization_uncertain_not_excluded"][
            "differences_vs_other_scenarios"
        ]
        validated = report["scenarios"]["sequence_validated_models_only"][
            "differences_vs_other_scenarios"
        ]
        # the relaxed scenario is a superset here, so nothing is unique to it
        self.assertEqual(relaxed["only_in_scenario_count"], 0)
        self.assertEqual(relaxed["not_in_scenario_count"], 0)
        # the frozen v1 rule holds the signal-less type-1 row that others keep
        self.assertEqual(strict["not_in_scenario_count"], 1)
        self.assertEqual(strict["only_in_scenario_count"], 0)
        # the validated-only scenario drops the discovery-layer-only row
        self.assertEqual(validated["not_in_scenario_count"], 1)
        self.assertEqual(validated["only_in_scenario_count"], 0)
        self.assertEqual(
            relaxed["compared_against"],
            ["current_v1_rules", "sequence_validated_models_only"],
        )
        self.assertIn("accession_differences.tsv", report["accession_list_location"])

    def test_accession_level_differences_are_available_for_the_tsv(self):
        analysis = module.analyze(
            rows=self.fixture_rows(),
            v1_labels={"T1_A": "high_confidence"},
            v2_labels={"T1_A": "core_sequence_homolog"},
            scenario_names=["current_v1_rules", "localization_uncertain_not_excluded"],
        )
        scenario_rows = [r for r in analysis.accession_rows if r["section"] == "scenario"]
        only = {r["accession"] for r in scenario_rows if r["change_kind"] == "only_in_scenario"}
        missing = {r["accession"] for r in scenario_rows if r["change_kind"] == "not_in_scenario"}
        self.assertIn("T1_NO_SIGNAL", only)
        self.assertIn("T1_NO_SIGNAL", missing)

    def test_v1_v2_changes_are_available_for_the_tsv(self):
        analysis = module.analyze(
            rows=self.fixture_rows(),
            v1_labels={"GONE": "high_confidence", "T1_A": "high_confidence"},
            v2_labels={"T1_A": "core_sequence_homolog", "NEW": "core_sequence_homolog"},
        )
        impact_rows = {r["accession"]: r for r in analysis.accession_rows if r["section"] == "v1_v2"}
        self.assertEqual(impact_rows["GONE"]["change_kind"], "removed")
        self.assertEqual(impact_rows["GONE"]["v1_label"], "high_confidence")
        self.assertEqual(impact_rows["NEW"]["change_kind"], "added")
        self.assertEqual(impact_rows["NEW"]["v2_label"], "core_sequence_homolog")
        self.assertEqual(impact_rows["T1_A"]["change_kind"], "relabelled")
        self.assertEqual(
            (impact_rows["T1_A"]["v1_label"], impact_rows["T1_A"]["v2_label"]),
            ("high_confidence", "core_sequence_homolog"),
        )

    def test_report_is_deterministic_across_row_order(self):
        rows = self.fixture_rows()
        forward = module.build_report(rows=rows, v1_labels={"A": "high_confidence"}, v2_labels={"A": "core_sequence_homolog"})
        backward = module.build_report(
            rows=list(reversed(rows)),
            v1_labels={"A": "high_confidence"},
            v2_labels={"A": "core_sequence_homolog"},
        )
        self.assertEqual(forward, backward)
        self.assertEqual(
            json.dumps(forward, sort_keys=True),
            json.dumps(backward, sort_keys=True),
        )

    def test_scenario_subset_selection_is_honoured(self):
        report = module.build_report(
            rows=self.fixture_rows(),
            v1_labels={},
            v2_labels={},
            scenario_names=["current_v1_rules"],
        )
        self.assertEqual(list(report["scenarios"]), ["current_v1_rules"])

    def test_unknown_scenario_in_the_report_builder_raises(self):
        with self.assertRaises(ValueError):
            module.build_report(
                rows=self.fixture_rows(),
                v1_labels={},
                v2_labels={},
                scenario_names=["does_not_exist"],
            )


class CompositionGuardTests(unittest.TestCase):
    """The 76:24 candidate composition can never be read as a biological ratio."""

    def build(self) -> dict:
        rows = [
            row("T1_A"),
            row("T1_B", superfamily=TYPE2, catalytic_domain_type="type2_verified"),
            cys_row("CYS_A"),
        ]
        return module.build_report(rows=rows, v1_labels={}, v2_labels={})

    def test_every_composition_block_carries_caveat_rule_and_scenario(self):
        blocks = module.composition_blocks(self.build())
        self.assertTrue(blocks)
        self.assertGreaterEqual(len(blocks), 4 * len(module.SCENARIOS))
        for block in blocks:
            self.assertEqual(block["caveat"], module.COMPOSITION_CAVEAT)
            self.assertEqual(block["basis"], module.COMPOSITION_BASIS)
            self.assertIn(block["scenario"], module.SCENARIOS)
            self.assertTrue(block["rule"].strip())

    def test_ratio_label_names_the_scenario_and_stays_a_candidate_composition(self):
        report = self.build()
        block = report["scenarios"]["current_v1_rules"]["compositions"]["candidate_localization"]
        self.assertIn("current_v1_rules", block["ratio_label"])
        self.assertIn("candidate composition", block["ratio_label"])
        self.assertIn("intracellular:extracellular", block["ratio_label"])
        self.assertIn("3", block["ratio_label"])

    def test_the_caveat_itself_is_free_of_forbidden_fragments(self):
        for fragment in module.FORBIDDEN_LABEL_FRAGMENTS:
            self.assertNotIn(fragment, module.COMPOSITION_CAVEAT.lower(), fragment)

    def test_the_whole_report_json_text_has_no_forbidden_label(self):
        text = json.dumps(self.build(), ensure_ascii=False).lower()
        for fragment in module.FORBIDDEN_LABEL_FRAGMENTS:
            self.assertNotIn(fragment, text, fragment)

    def test_a_biological_ratio_label_is_rejected(self):
        poisoned = json.loads(json.dumps(self.build()))
        poisoned["scenarios"]["current_v1_rules"]["compositions"]["candidate_localization"][
            "ratio_label"
        ] = "true_ratio 76:24"
        with self.assertRaisesRegex(ValueError, "forbidden"):
            module.validate_report_guard(poisoned)

    def test_a_biological_ratio_key_is_rejected(self):
        poisoned = json.loads(json.dumps(self.build()))
        poisoned["scenarios"]["current_v1_rules"]["compositions"]["candidate_localization"][
            "counts"
        ]["in_vivo_ratio"] = 1
        with self.assertRaisesRegex(ValueError, "forbidden"):
            module.validate_report_guard(poisoned)

    def test_a_nested_biological_ratio_key_is_rejected(self):
        poisoned = json.loads(json.dumps(self.build()))
        poisoned["scenarios"]["current_v1_rules"]["compositions"]["candidate_localization"][
            "nested"
        ] = {"biological_ratio": {"intracellular": 1}}
        with self.assertRaisesRegex(ValueError, "forbidden"):
            module.validate_report_guard(poisoned)

    def test_a_composition_block_without_caveat_or_scenario_is_rejected(self):
        for field in ("caveat", "scenario", "rule"):
            poisoned = json.loads(json.dumps(self.build()))
            del poisoned["scenarios"]["current_v1_rules"]["compositions"]["family"][field]
            with self.assertRaises(ValueError, msg=field):
                module.validate_report_guard(poisoned)

    def test_a_composition_block_with_a_foreign_scenario_name_is_rejected(self):
        poisoned = json.loads(json.dumps(self.build()))
        poisoned["scenarios"]["current_v1_rules"]["compositions"]["family"]["scenario"] = "unknown_rule"
        with self.assertRaises(ValueError):
            module.validate_report_guard(poisoned)

    def test_a_composition_moved_outside_the_compositions_key_is_rejected(self):
        poisoned = json.loads(json.dumps(self.build()))
        block = poisoned["scenarios"]["current_v1_rules"]["compositions"].pop("superfamily")
        poisoned["scenarios"]["current_v1_rules"]["loose_composition"] = block
        with self.assertRaises(ValueError):
            module.validate_report_guard(poisoned)


# --- CLI ----------------------------------------------------------------------


class CliTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="dsh_compare_catalog_")
        self.tmp = Path(self._tmp.name)
        self.v1_path = write_tsv(
            self.tmp / "v1_high_confidence.tsv",
            [
                {"accession": "T1_A", "high_confidence": "1"},
                {"accession": "GONE", "high_confidence": "1"},
            ],
        )
        self.v2_path = write_tsv(
            self.tmp / "phaded_candidate_catalog_v2.tsv",
            [
                dict(row("T1_A"), primary_disposition="core_sequence_homolog"),
                dict(
                    row("T1_NO_SIGNAL", signalp_class="OTHER", transport_signal_prediction="OTHER"),
                    primary_disposition="probable_sequence_homolog",
                ),
                dict(
                    row("DISCOVERY_ONLY", model_layer="discovery_hmm_uncalibrated"),
                    primary_disposition="remote_homolog_candidate",
                ),
            ],
        )
        self.out_dir = self.tmp / "out"

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def run_cli(self, *extra: str):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = module.main(
                [
                    "--v1", str(self.v1_path),
                    "--v2", str(self.v2_path),
                    "--output-dir", str(self.out_dir),
                    *extra,
                ]
            )
        return code, stdout.getvalue(), stderr.getvalue()

    def test_cli_writes_json_and_tsv_and_prints_a_json_summary(self):
        code, stdout, stderr = self.run_cli()
        self.assertEqual(code, 0, stderr)
        report_path = self.out_dir / "catalog_version_comparison.json"
        tsv_path = self.out_dir / "accession_differences.tsv"
        self.assertTrue(report_path.is_file())
        self.assertTrue(tsv_path.is_file())
        written = json.loads(report_path.read_text(encoding="utf-8"))
        printed = json.loads(stdout)
        self.assertEqual(written, printed)
        self.assertEqual(set(written["scenarios"]), set(module.SCENARIOS))
        with tsv_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(rows)
        self.assertTrue(all(r["section"] in {"v1_v2", "scenario"} for r in rows))

    def test_cli_writes_nothing_outside_the_output_dir(self):
        before = {p for p in self.tmp.rglob("*")}
        code, _stdout, stderr = self.run_cli()
        self.assertEqual(code, 0, stderr)
        after = {p for p in self.tmp.rglob("*")}
        created = after - before
        self.assertTrue(created)
        for path in created:
            self.assertTrue(
                str(path).startswith(str(self.out_dir)),
                f"wrote outside --output-dir: {path}",
            )

    def test_cli_refuses_an_existing_non_empty_output_dir(self):
        self.out_dir.mkdir(parents=True)
        keep = self.out_dir / "keep.txt"
        keep.write_text("frozen evidence\n", encoding="utf-8")
        code, _stdout, stderr = self.run_cli()
        self.assertEqual(code, 2)
        self.assertIn("not empty", stderr)
        self.assertEqual(keep.read_text(encoding="utf-8"), "frozen evidence\n")
        self.assertEqual(sorted(p.name for p in self.out_dir.iterdir()), ["keep.txt"])

    def test_cli_refuses_a_missing_input_path(self):
        code, _stdout, stderr = self.run_cli_and_report_missing()
        self.assertEqual(code, 2)
        self.assertIn("does not exist", stderr)

    def run_cli_and_report_missing(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = module.main(
                [
                    "--v1", str(self.tmp / "missing.tsv"),
                    "--v2", str(self.v2_path),
                    "--output-dir", str(self.out_dir),
                ]
            )
        return code, stdout.getvalue(), stderr.getvalue()

    def test_cli_rejects_an_unknown_scenario_and_lists_valid_names(self):
        code, _stdout, stderr = self.run_cli("--scenarios", "typo_scenario")
        self.assertEqual(code, 2)
        for name in module.SCENARIOS:
            self.assertIn(name, stderr)
        self.assertFalse(self.out_dir.exists())

    def test_cli_scenario_subset_is_written(self):
        code, stdout, stderr = self.run_cli("--scenarios", "localization_uncertain_not_excluded")
        self.assertEqual(code, 0, stderr)
        self.assertEqual(
            list(json.loads(stdout)["scenarios"]),
            ["localization_uncertain_not_excluded"],
        )

    def test_cli_list_scenarios_prints_the_preregistered_names(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = module.main(["--list-scenarios"])
        self.assertEqual(code, 0, stderr)
        payload = json.loads(stdout.getvalue())
        self.assertEqual(set(payload), set(module.SCENARIOS))
        for name, description in payload.items():
            self.assertTrue(description.strip(), name)
        self.assertFalse(self.out_dir.exists())

    def test_cli_confounder_file_is_optional_and_recorded(self):
        confounders = write_tsv(self.tmp / "confounders.tsv", [{"accession": "T1_A"}])
        code, stdout, stderr = self.run_cli("--confounders", str(confounders))
        self.assertEqual(code, 0, stderr)
        report = json.loads(stdout)
        self.assertEqual(report["inputs"]["confounders"]["status"], "verified")
        strict = report["scenarios"]["current_v1_rules"]
        self.assertEqual(strict["first_failing_criteria"]["T1_A"], "v1_not_confounder")
        self.assertEqual(
            strict["excluded_reason_counts"],
            {"common_criteria": 1, "localization_conflict": 1},
        )

    def test_confounder_file_is_optional(self):
        code, stdout, stderr = self.run_cli()
        self.assertEqual(code, 0, stderr)
        report = json.loads(stdout)
        self.assertEqual(report["inputs"]["confounders"]["status"], "not_supplied")
        strict = report["scenarios"]["current_v1_rules"]
        self.assertNotIn("v1_not_confounder", strict["pending_criteria"])
        self.assertNotIn("v1_not_confounder", strict["first_failing_criteria"].values())
        self.assertEqual(strict["excluded_reason_counts"], {"localization_conflict": 1})

    def test_cli_records_input_provenance(self):
        code, stdout, stderr = self.run_cli()
        self.assertEqual(code, 0, stderr)
        inputs = json.loads(stdout)["inputs"]
        for key in ("v1", "v2", "evidence"):
            self.assertEqual(inputs[key]["status"], "verified", key)
            self.assertEqual(len(inputs[key]["sha256"]), 64, key)


if __name__ == "__main__":
    unittest.main()
