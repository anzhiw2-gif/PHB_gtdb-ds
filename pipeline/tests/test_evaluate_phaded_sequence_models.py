"""Tests for the independent sequence-family model evaluator (Task 6).

Written before the implementation (test-first).  The evaluator decides whether a
*sequence* family model earns the ``sequence_family_hmm_validated`` layer from
held-out recovery and nearest-family/competing-fold confusability alone.

Load-bearing invariants (never relax):

* an experimental positive is **not** a precondition for building or validating a
  discovery HMM or a sequence family HMM (approved design assumption 1): a family
  with ``experimental_anchor_count=0`` that passes every sequence gate is
  validated, and a family with three experimental anchors that fails a sequence
  gate is *not*;
* ``functional_calibration_status`` is copied from the input row verbatim and is
  never upgraded here -- this module emits no ``calibrated_candidate_model`` at
  all, on any layer, for any input;
* a held-out panel smaller than five must recover **every** held-out sequence;
  for panels of five or more the configured minimum recall applies, and the two
  rules meet exactly at ``heldout_count=5, heldout_recovered=4`` (4/5 = 0.8);
* ``EFFN`` is a diagnostic.  ``eff_nseq`` below one is reported and never blocks.
"""

import ast
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


TESTS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = TESTS_DIR.parent / "scripts"
SCRIPT = SCRIPTS_DIR / "evaluate_phaded_sequence_models.py"
SCHEMA_SCRIPT = SCRIPTS_DIR / "phaded_evidence_schema.py"
GATES = TESTS_DIR.parent / "config" / "phaded_sequence_model_gates.tsv"

#: The preregistered ``sequence-gate-v1`` row, as ``load_sequence_model_gates``
#: returns it (``dict[str, str]``).
GATE_FIELDS = [
    "gate_version", "min_training_sequences", "min_heldout_sequences", "min_heldout_recall",
    "max_unexplained_confounder_hits", "require_unique_best_model", "require_alignment_hash",
]
GATE_ROW = dict(zip(
    GATE_FIELDS,
    ["sequence-gate-v1", "3", "1", "0.8", "0", "true", "true"],
))

HEX64 = "a" * 64

PROFILE_ONLY_FIELDS = [
    "profile_id", "profile_kind", "model_layer",
]


def load_module():
    """Load the evaluator from disk; it does not exist yet, so this fails first."""
    spec = importlib.util.spec_from_file_location(
        "evaluate_phaded_sequence_models", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_schema():
    spec = importlib.util.spec_from_file_location("phaded_evidence_schema", SCHEMA_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gates(**overrides):
    """The preregistered gate row (strings, exactly as loaded from the TSV)."""
    row = dict(GATE_ROW)
    row.update({key: str(value) for key, value in overrides.items()})
    return row


def record(**overrides):
    """A manifest row that passes every sequence gate unless overridden.

    Only ``overrides`` decide the outcome of a test; every other field is set to
    a value that satisfies the ``sequence-gate-v1`` legs.
    """
    row = {
        "profile_id": "family_DED_hfam_70_9818be7f78e3",
        "profile_kind": "family",
        "model_layer": "discovery_hmm_uncalibrated",
        "functional_calibration_status": "not_function_calibrated",
        "training_sequence_count": "8",
        "experimental_anchor_count": "0",
        "annotation_only_count": "8",
        "training_accessions": "AAB02914.1;AAA;BBB",
        "training_set_sha256": "b" * 64,
        "alignment_sha256": HEX64,
        "hmm_sha256": "d" * 64,
        "mafft_version": "7.525",
        "hmmer_version": "3.4",
        "bit_reproducible": "false",
        # evaluation inputs supplied by the CLI
        "heldout_count": "5",
        "heldout_recovered": "5",
        "unexplained_confounder_hits": "0",
        "unique_best_model": "true",
        "eff_nseq": "0.748535",
        "coverage_fraction": "0.97",
        "family_assignment_stable": "true",
        "nearest_family_score": "210.4",
        "competing_fold_score": "31.2",
        "sequence_family_call": "",
    }
    row.update({key: ("" if value is None else str(value)) for key, value in overrides.items()})
    return row


class FailingTestFirstCheck(unittest.TestCase):
    """The evaluator is written after this module (test-first)."""

    def test_module_exists(self):
        self.assertTrue(SCRIPT.is_file(), "missing %s" % SCRIPT)


class SequenceGateTests(unittest.TestCase):
    """The gate is earned by sequence performance, never by experimental count."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_sequence_gate_does_not_require_experimental_positive(self):
        result = self.module.evaluate_model(record(
            training_sequence_count=8,
            heldout_count=5,
            heldout_recovered=5,
            unexplained_confounder_hits=0,
            unique_best_model=True,
            alignment_sha256="a" * 64,
            experimental_anchor_count=0,
        ), gates())
        self.assertEqual(result["model_layer"], "sequence_family_hmm_validated")
        self.assertEqual(result["functional_calibration_status"], "not_function_calibrated")

    def test_effn_below_one_is_diagnostic_not_an_automatic_failure(self):
        item = record(eff_nseq=0.7)
        self.assertNotIn("eff_nseq_below_one", self.module.blocking_reasons(item, gates()))

    def test_effn_is_reported_as_a_diagnostic_field(self):
        result = self.module.evaluate_model(record(eff_nseq=0.7), gates())
        self.assertEqual(result["metrics"]["eff_nseq"], "0.7")
        self.assertEqual(result["blocking_reasons"], [])
        self.assertEqual(result["model_layer"], "sequence_family_hmm_validated")

    def test_a_passing_row_has_no_blocking_reason(self):
        self.assertEqual(self.module.blocking_reasons(record(), gates()), [])

    def test_three_experimental_anchors_do_not_earn_the_layer(self):
        """A real experimental anchor count >= 3 never bypasses a failed sequence leg."""
        result = self.module.evaluate_model(
            record(experimental_anchor_count=3, unexplained_confounder_hits=2), gates()
        )
        self.assertEqual(result["model_layer"], "discovery_hmm_uncalibrated")
        self.assertIn("unexplained_confounder_hits_present", result["blocking_reasons"])
        self.assertFalse(result["promoted"])

    def test_experimental_count_is_recorded_but_never_gates(self):
        zero = self.module.evaluate_model(record(experimental_anchor_count=0), gates())
        many = self.module.evaluate_model(record(experimental_anchor_count=9), gates())
        self.assertEqual(zero["promoted"], many["promoted"])
        self.assertEqual(zero["metrics"]["experimental_anchor_count"], 0)
        self.assertEqual(many["metrics"]["experimental_anchor_count"], 9)
        self.assertIsNone(
            self.module.compute_metrics(record(experimental_anchor_count=""))[
                "experimental_anchor_count"
            ],
            "an absent experimental count stays pending, never 0",
        )

    def test_min_training_sequences_blocks_at_two(self):
        reasons = self.module.blocking_reasons(record(training_sequence_count=2), gates())
        self.assertIn("training_sequences_below_minimum", reasons)
        self.assertNotIn("training_sequences_below_minimum", self.module.blocking_reasons(
            record(training_sequence_count=3), gates()
        ))

    def test_training_count_is_read_not_assumed_from_accession_string(self):
        """Fail closed: the count column decides, an absent count is ``pending``."""
        reasons = self.module.blocking_reasons(
            record(training_sequence_count="", training_accessions="A;B;C"), gates()
        )
        self.assertIn("training_sequence_count_pending", reasons)


class HeldoutRecallTests(unittest.TestCase):
    """Held-out rules, including the explicit smaller-than-five panel rule."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_small_panel_requires_every_sequence_recovered(self):
        for size in range(1, 5):
            with self.subTest(heldout_count=size):
                self.assertEqual(
                    self.module.blocking_reasons(
                        record(heldout_count=size, heldout_recovered=size), gates()
                    ),
                    [],
                )
                self.assertIn(
                    "heldout_recall_below_minimum",
                    self.module.blocking_reasons(
                        record(heldout_count=size, heldout_recovered=size - 1), gates()
                    ),
                )

    def test_four_of_four_passes_and_four_of_three_is_impossible(self):
        self.assertEqual(
            self.module.blocking_reasons(record(heldout_count=4, heldout_recovered=4), gates()),
            [],
        )
        self.assertIn(
            "heldout_recovered_exceeds_count",
            self.module.blocking_reasons(record(heldout_count=4, heldout_recovered=5), gates()),
        )

    def test_three_of_four_fails_the_small_panel_rule(self):
        self.assertIn(
            "heldout_recall_below_minimum",
            self.module.blocking_reasons(record(heldout_count=4, heldout_recovered=3), gates()),
        )

    def test_five_of_four_passes_at_the_configured_recall(self):
        """4/5 = 0.8 meets the configured minimum: the two rules meet exactly here."""
        detail = {}
        reasons = self.module.blocking_reasons(
            record(heldout_count=5, heldout_recovered=4), gates(), detail=detail
        )
        self.assertEqual(reasons, [])
        self.assertEqual(detail["heldout_rule"], "min_recall_at_or_above_five")
        self.assertEqual(detail["heldout_recall"], "0.8")

    def test_small_panel_rule_is_encoded_explicitly(self):
        detail = {}
        self.module.blocking_reasons(
            record(heldout_count=4, heldout_recovered=4), gates(), detail=detail
        )
        self.assertEqual(detail["heldout_rule"], "all_recovered_below_five")
        self.assertEqual(detail["small_heldout_panel_threshold"], 5)

    def test_panel_below_the_minimum_heldout_count_blocks(self):
        self.assertIn(
            "heldout_panel_below_minimum",
            self.module.blocking_reasons(record(heldout_count=0, heldout_recovered=0), gates()),
        )

    def test_pending_heldout_values_fail_closed(self):
        reasons = self.module.blocking_reasons(
            record(heldout_count="pending", heldout_recovered="pending"), gates()
        )
        self.assertIn("heldout_count_pending", reasons)
        self.assertIn("heldout_recovered_pending", reasons)

    def test_missing_heldout_columns_fail_closed(self):
        item = record()
        item.pop("heldout_count")
        item.pop("heldout_recovered")
        reasons = self.module.blocking_reasons(item, gates())
        self.assertIn("heldout_count_pending", reasons)
        self.assertIn("heldout_recovered_pending", reasons)

    def test_min_heldout_recall_override_is_honoured_for_large_panels(self):
        """Panels of five or more use the configured recall, not the all-recovered rule."""
        self.assertEqual(
            self.module.blocking_reasons(record(heldout_count=10, heldout_recovered=9), gates()),
            [],
        )
        self.assertEqual(
            self.module.blocking_reasons(record(heldout_count=10, heldout_recovered=8), gates()),
            [],
        )
        detail = {}
        self.assertIn(
            "heldout_recall_below_minimum",
            self.module.blocking_reasons(
                record(heldout_count=10, heldout_recovered=7), gates(), detail=detail
            ),
        )
        self.assertEqual(detail["heldout_rule"], "min_recall_at_or_above_five")
        self.assertEqual(detail["heldout_recall"], "0.7")


class ConfounderAndUniquenessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_confounder_hits_are_counted_exactly(self):
        gate_row = gates(max_unexplained_confounder_hits=1)
        self.assertEqual(
            self.module.blocking_reasons(record(unexplained_confounder_hits=1), gate_row), []
        )
        reasons = self.module.blocking_reasons(record(unexplained_confounder_hits=2), gate_row)
        self.assertIn("unexplained_confounder_hits_present", reasons)
        self.assertNotIn("unexplained_confounder_hits_present", self.module.blocking_reasons(
            record(unexplained_confounder_hits=0), gates()
        ))

    def test_confounder_count_is_exact_and_pending_fails_closed(self):
        result = self.module.evaluate_model(
            record(unexplained_confounder_hits=7), gates(max_unexplained_confounder_hits=7)
        )
        self.assertEqual(result["metrics"]["unexplained_confounder_hits"], 7)
        self.assertEqual(result["blocking_reasons"], [])
        self.assertIn(
            "unexplained_confounder_hits_pending",
            self.module.blocking_reasons(record(unexplained_confounder_hits=""), gates()),
        )

    def test_require_unique_best_model_blocks_a_non_unique_winner(self):
        self.assertIn(
            "not_unique_best_model",
            self.module.blocking_reasons(record(unique_best_model="false"), gates()),
        )
        self.assertNotIn(
            "not_unique_best_model",
            self.module.blocking_reasons(record(unique_best_model="true"), gates()),
        )

    def test_require_unique_best_model_pending_fails_closed(self):
        reasons = self.module.blocking_reasons(record(unique_best_model="pending"), gates())
        self.assertIn("unique_best_model_pending", reasons)

    def test_unique_best_model_is_not_required_when_the_gate_says_so(self):
        self.assertEqual(
            self.module.blocking_reasons(
                record(unique_best_model="false"), gates(require_unique_best_model="false")
            ),
            [],
        )

    def test_final_score_margin_is_a_metric_and_a_missing_margin_stays_pending(self):
        result = self.module.evaluate_model(record(), gates())
        self.assertEqual(result["metrics"]["score_margin"], "179.2")
        blank = self.module.evaluate_model(
            record(nearest_family_score="", competing_fold_score=""), gates()
        )
        self.assertEqual(blank["metrics"]["score_margin"], "pending")
        self.assertEqual(blank["blocking_reasons"], [])


class AlignmentHashTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_missing_pending_and_empty_alignment_hash_block(self):
        for value in ("", "pending", "PENDING", "not_available"):
            with self.subTest(alignment_sha256=value):
                self.assertIn(
                    "alignment_sha256_pending",
                    self.module.blocking_reasons(record(alignment_sha256=value), gates()),
                )

    def test_missing_alignment_hash_column_blocks(self):
        item = record()
        item.pop("alignment_sha256")
        self.assertIn(
            "alignment_sha256_pending", self.module.blocking_reasons(item, gates())
        )

    def test_a_64_hex_alignment_hash_passes(self):
        self.assertNotIn(
            "alignment_sha256_pending",
            self.module.blocking_reasons(record(alignment_sha256=HEX64), gates()),
        )
        self.assertNotIn(
            "alignment_sha256_invalid",
            self.module.blocking_reasons(record(alignment_sha256="0123456789abcdef" * 4), gates()),
        )

    def test_a_malformed_short_alignment_hash_blocks(self):
        reasons = self.module.blocking_reasons(record(alignment_sha256="abc123"), gates())
        self.assertIn("alignment_sha256_invalid", reasons)

    def test_alignment_hash_is_not_required_when_the_gate_says_so(self):
        self.assertEqual(
            self.module.blocking_reasons(
                record(alignment_sha256=""), gates(require_alignment_hash="false")
            ),
            [],
        )
        self.assertIn(
            "alignment_sha256_invalid",
            self.module.blocking_reasons(
                record(alignment_sha256="zz"), gates(require_alignment_hash="false")
            ),
        )

    def test_bit_non_reproducibility_is_recorded_and_never_blocks(self):
        """MAFFT 7.525 is not bit-reproducible; the alignment hash is mandatory instead."""
        result = self.module.evaluate_model(record(bit_reproducible="false"), gates())
        self.assertEqual(result["blocking_reasons"], [])
        self.assertEqual(result["metrics"]["bit_reproducible"], "false")


class ModelLayerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_reference_query_only_is_never_promoted(self):
        result = self.module.evaluate_model(
            record(model_layer="reference_query_only", training_sequence_count=2), gates()
        )
        self.assertEqual(result["model_layer"], "reference_query_only")
        self.assertIn("model_layer_reference_query_only_unpromotable", result["blocking_reasons"])

    def test_discovery_layer_without_a_family_call_is_promotable(self):
        result = self.module.evaluate_model(record(model_layer="discovery_hmm_uncalibrated"), gates())
        self.assertEqual(result["model_layer"], "sequence_family_hmm_validated")
        self.assertEqual(result["input_model_layer"], "discovery_hmm_uncalibrated")

    def test_discovery_layer_may_not_carry_a_family_call(self):
        reasons = self.module.blocking_reasons(
            record(model_layer="discovery_hmm_uncalibrated", sequence_family_call="DED_hfam_70"),
            gates(),
        )
        self.assertIn("discovery_family_call_present", reasons)

    def test_family_call_aliases_are_checked_too(self):
        self.assertIn(
            "discovery_family_call_present",
            self.module.blocking_reasons(record(family_call="DED_hfam_52"), gates()),
        )

    def test_an_already_validated_layer_keeps_its_family_call(self):
        result = self.module.evaluate_model(
            record(model_layer="sequence_family_hmm_validated", sequence_family_call="DED_hfam_70"),
            gates(),
        )
        self.assertEqual(result["model_layer"], "sequence_family_hmm_validated")
        self.assertEqual(result["blocking_reasons"], [])

    def test_unknown_layer_fails_closed(self):
        reasons = self.module.blocking_reasons(record(model_layer="whatever"), gates())
        self.assertIn("model_layer_unpromotable", reasons)
        self.assertIn("model_layer_unknown", reasons)

    def test_an_unknown_layer_is_echoed_verbatim_and_never_renamed(self):
        result = self.module.evaluate_model(record(model_layer="whatever"), gates())
        self.assertEqual(result["model_layer"], "whatever")
        self.assertEqual(result["input_model_layer"], "pending")
        self.assertFalse(result["promoted"])

    def test_a_row_without_a_layer_cell_reports_pending(self):
        item = record()
        item.pop("model_layer")
        result = self.module.evaluate_model(item, gates())
        self.assertEqual(result["model_layer"], "pending")
        self.assertEqual(result["input_model_layer"], "pending")
        self.assertIn("model_layer_unknown", result["blocking_reasons"])


    def test_the_module_never_invents_the_calibrated_layer(self):
        """A calibrated input layer is echoed verbatim but never produced here."""
        result = self.module.evaluate_model(
            record(model_layer="calibrated_candidate_model"), gates()
        )
        self.assertEqual(result["model_layer"], "calibrated_candidate_model")
        self.assertFalse(result["promoted"])
        self.assertIn("model_layer_unpromotable", result["blocking_reasons"])
        module = load_module()
        for layer in (module.MODEL_LAYER, *sorted(module.PROMOTABLE_MODEL_LAYERS)):
            self.assertNotEqual(layer, "calibrated_candidate_model")

    def test_rows_never_report_a_layer_the_module_did_not_earn(self):
        """No code path promotes a blocked row, whichever single leg fails."""
        matrix = [
            ("confounders", record(unexplained_confounder_hits=1)),
            ("heldout", record(heldout_count=4, heldout_recovered=3)),
            ("unique best", record(unique_best_model="false")),
            ("alignment hash", record(alignment_sha256="zz")),
            ("training count", record(training_sequence_count=2)),
            ("pending values", record(heldout_count="", heldout_recovered="")),
            ("unpromotable layer", record(model_layer="reference_query_only")),
        ]
        for label, item in matrix:
            with self.subTest(leg=label):
                result = self.module.evaluate_model(item, gates())
                self.assertFalse(result["promoted"])
                self.assertEqual(result["model_layer"], item["model_layer"])
                self.assertNotEqual(result["model_layer"], "sequence_family_hmm_validated")


class OverclaimingBoundaryTests(unittest.TestCase):
    """Every output carries the sequence-claim boundary; calibration never moves."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_every_output_carries_the_three_boundary_keys(self):
        matrix = [
            record(),
            record(model_layer="reference_query_only"),
            record(model_layer="discovery_hmm_uncalibrated", sequence_family_call="DED_hfam_4"),
            record(heldout_count=4, heldout_recovered=3),
            record(unexplained_confounder_hits=5),
            record(alignment_sha256="pending"),
            record(training_sequence_count=2),
            record(functional_calibration_status="candidate_gate_passed_not_promoted"),
            record(functional_calibration_status="not_function_calibrated"),
        ]
        for item in matrix:
            with self.subTest(profile_id=item["profile_id"], layer=item["model_layer"]):
                result = self.module.evaluate_model(item, gates())
                self.assertEqual(result["sequence_claim"], "sequence-defined family homology")
                self.assertEqual(
                    result["phenotype_boundary"], "does not establish PHB/PHA degradation"
                )
                self.assertEqual(
                    result["functional_calibration_status"],
                    item["functional_calibration_status"],
                )
                self.assertNotEqual(result["model_layer"], "calibrated_candidate_model")

    def test_functional_status_is_copied_verbatim_when_every_sequence_gate_passes(self):
        result = self.module.evaluate_model(
            record(
                functional_calibration_status="not_function_calibrated",
                experimental_anchor_count=0,
            ),
            gates(),
        )
        self.assertEqual(result["model_layer"], "sequence_family_hmm_validated")
        self.assertEqual(result["functional_calibration_status"], "not_function_calibrated")

    def test_a_passed_gate_is_not_a_promotion_to_calibration(self):
        result = self.module.evaluate_model(record(), gates())
        self.assertTrue(result["promoted"])
        self.assertEqual(
            result["functional_calibration_status"], "not_function_calibrated"
        )
        self.assertFalse(result["functional_calibration_changed"])

    def test_a_missing_functional_status_column_defaults_to_not_calibrated(self):
        item = record()
        item.pop("functional_calibration_status")
        result = self.module.evaluate_model(item, gates())
        self.assertEqual(result["functional_calibration_status"], "not_function_calibrated")

    def test_calibrated_status_in_the_input_is_never_invented_by_this_module(self):
        """An unrecognized/calibrated input status is copied verbatim, never promoted."""
        result = self.module.evaluate_model(
            record(functional_calibration_status="calibrated_candidate_model"), gates()
        )
        self.assertEqual(
            result["functional_calibration_status"], "calibrated_candidate_model"
        )
        self.assertEqual(result["model_layer"], "discovery_hmm_uncalibrated")
        self.assertIn(
            "input_functional_calibration_status_not_allowed", result["blocking_reasons"]
        )


class SourceBoundaryTests(unittest.TestCase):
    """Static proof that no code path can emit a calibrated model."""

    @classmethod
    def setUpClass(cls):
        cls.source = SCRIPT.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)

    def _literal_strings_outside_docstrings(self):
        docstrings = set()
        for node in ast.walk(self.tree):
            if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                body = getattr(node, "body", [])
                if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
                    if isinstance(body[0].value.value, str):
                        docstrings.add(id(body[0].value))
        return [
            (node.lineno, node.value)
            for node in ast.walk(self.tree)
            if isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and id(node) not in docstrings
        ]

    def test_calibrated_candidate_model_appears_in_no_executable_string(self):
        literals = self._literal_strings_outside_docstrings()
        offenders = [item for item in literals if "calibrated_candidate_model" in item[1]]
        self.assertEqual(offenders, [], "the evaluator must never name the calibrated layer")

    def test_calibrated_candidate_model_is_only_explained_in_the_docstring(self):
        self.assertIn("calibrated_candidate_model", self.source)
        for lineno, value in self._literal_strings_outside_docstrings():
            self.assertNotIn("calibrated_candidate_model", value, "line %d" % lineno)

    def _pinned_literal(self, name):
        """Return the literal string set assigned to a module-level constant.

        Accepts both ``X = {"a", "b"}`` and ``X = frozenset({"a", "b"})``.
        """
        for node in self.tree.body:
            if not isinstance(node, ast.Assign):
                continue
            if not any(
                isinstance(target, ast.Name) and target.id == name for target in node.targets
            ):
                continue
            value = node.value
            if isinstance(value, ast.Call):  # frozenset({...})
                value = value.args[0]
            return {
                element.value
                for element in getattr(value, "elts", [])
                if isinstance(element, ast.Constant)
            }
        return None

    def test_promotable_layers_are_pinned_to_the_two_expected_layers(self):
        pinned = self._pinned_literal("PROMOTABLE_MODEL_LAYERS")
        self.assertIsNotNone(pinned, "PROMOTABLE_MODEL_LAYERS must be a module-level constant")
        self.assertEqual(
            pinned, {"discovery_hmm_uncalibrated", "sequence_family_hmm_validated"}
        )
        module = load_module()
        self.assertNotIn("calibrated_candidate_model", module.PROMOTABLE_MODEL_LAYERS)
        self.assertEqual(module.MODEL_LAYER, "sequence_family_hmm_validated")

    def test_the_calibrated_status_is_not_in_the_allowlist(self):
        """Static proof: the statuses this module accepts are read from the schema."""
        assignments = [
            node
            for node in self.tree.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id == "_ALLOWED_FUNCTIONAL_STATUSES"
                for target in node.targets
            )
        ]
        self.assertEqual(len(assignments), 1, "_ALLOWED_FUNCTIONAL_STATUSES must be pinned once")
        module = load_module()
        schema = load_schema()
        self.assertEqual(
            set(module._ALLOWED_FUNCTIONAL_STATUSES),
            set(schema.ALLOWED_MODEL_STATES["sequence_family_hmm_validated"]),
        )

    def test_no_assignment_target_or_keyword_is_the_functional_calibration_status(self):
        """The status may only be *read* from the row, never assigned or constructed."""
        for node in ast.walk(self.tree):
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        self.assertNotEqual(target.id, "functional_calibration_status")
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for argument in list(node.args.args) + list(node.args.kwonlyargs):
                    self.assertNotEqual(argument.arg, "functional_calibration_status")
            if isinstance(node, ast.Dict):
                values = list(node.values)
                for index, key in enumerate(node.keys):
                    if isinstance(key, ast.Constant) and key.value == "functional_calibration_status":
                        # The only sanctioned write is the copied input value.
                        self.assertIsInstance(
                            values[index], (ast.Name, ast.Call, ast.Subscript),
                            "line %d writes functional_calibration_status from a literal"
                            % getattr(key, "lineno", -1),
                        )


class OutputFileTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "inputs").mkdir()
        self.profiles = self.root / "inputs" / "profile_manifest.tsv"
        self.gates_path = self.root / "inputs" / "phaded_sequence_model_gates.tsv"
        self.out = self.root / "results"

    def tearDown(self):
        self.tmp.cleanup()

    def _write(self, path, fieldnames, rows):
        with Path(path).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fieldnames})

    def _seed(self, profile_rows, evaluation_rows=None):
        fields = sorted(set().union(*[set(row) for row in profile_rows]))
        self._write(self.profiles, fields, profile_rows)
        self._write(self.gates_path, GATE_FIELDS, [GATE_ROW])
        evaluation = None
        if evaluation_rows is not None:
            evaluation = self.root / "inputs" / "evaluation_inputs.tsv"
            eval_fields = sorted(set().union(*[set(row) for row in evaluation_rows]))
            self._write(evaluation, eval_fields, evaluation_rows)
        return evaluation

    def _manifest_rows(self):
        base = record()
        pending = record(
            profile_id="family_DED_hfam_4_pnone",
            model_layer="discovery_hmm_uncalibrated",
            heldout_count="4",
            heldout_recovered="3",
        )
        return [base, pending]

    def test_cli_writes_the_three_required_outputs_and_prints_json(self):
        self._seed(self._manifest_rows())
        payload = self.module.run(
            self.module.build_arg_namespace([
                str(self.profiles), str(self.out), "--gates", str(self.gates_path),
            ])
        )
        for name in (
            "sequence_model_evaluation.tsv",
            "sequence_model_gate_detail.tsv",
            "sequence_model_evaluation_manifest.json",
        ):
            self.assertTrue((self.out / name).is_file(), name)
        self.assertEqual(payload["counts"]["models_evaluated"], 2)
        self.assertEqual(payload["counts"]["models_validated"], 1)
        self.assertEqual(payload["model_layer"], "sequence_family_hmm_validated")
        self.assertIn("gates", payload)

    def test_cli_refuses_to_overwrite_an_existing_output_directory(self):
        self._seed(self._manifest_rows())
        (self.out).mkdir(parents=True)
        (self.out / "sequence_model_evaluation.tsv").write_text("keep\n", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            self.module.run(self.module.build_arg_namespace([
                str(self.profiles), str(self.out), "--gates", str(self.gates_path),
            ]))
        self.assertEqual(
            (self.out / "sequence_model_evaluation.tsv").read_text(encoding="utf-8"), "keep\n"
        )

    def test_cli_validates_that_every_input_path_exists(self):
        self._seed(self._manifest_rows())
        missing = self.root / "inputs" / "nope.tsv"
        with self.assertRaises(FileNotFoundError):
            self.module.run(self.module.build_arg_namespace([
                str(self.profiles), str(self.out), "--gates", str(missing),
            ]))
        with self.assertRaises(FileNotFoundError):
            self.module.run(self.module.build_arg_namespace([
                str(missing), str(self.out), "--gates", str(self.gates_path),
            ]))
        with self.assertRaises(FileNotFoundError):
            self.module.run(self.module.build_arg_namespace([
                str(self.profiles), str(self.out), "--gates", str(self.gates_path),
                "--evaluation", str(missing),
            ]))

    def test_main_returns_one_for_a_missing_input_path(self):
        self._seed(self._manifest_rows())
        code = self.module.main([
            str(self.profiles), str(self.out), "--gates", str(self.root / "inputs" / "nope.tsv"),
        ])
        self.assertEqual(code, 1)

    def test_main_returns_zero_and_prints_a_json_summary(self):
        import contextlib
        import io

        self._seed(self._manifest_rows())
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            code = self.module.main([
                str(self.profiles), str(self.out), "--gates", str(self.gates_path),
            ])
        self.assertEqual(code, 0)
        payload = json.loads(stream.getvalue().strip().splitlines()[-1])
        self.assertEqual(payload["counts"]["models_evaluated"], 2)

    def test_manifest_records_input_hashes_and_the_gate_version(self):
        self._seed(self._manifest_rows())
        payload = self.module.run(self.module.build_arg_namespace([
            str(self.profiles), str(self.out), "--gates", str(self.gates_path),
        ]))
        manifest = payload["manifest"]
        self.assertEqual(manifest["gate_version"], "sequence-gate-v1")
        self.assertEqual(manifest["model_layer"], "sequence_family_hmm_validated")
        self.assertEqual(manifest["inputs"]["profile_manifest"]["sha256"], self.module.sha256_file(
            self.profiles
        ))
        self.assertEqual(manifest["inputs"]["sequence_model_gates"]["status"], "verified")
        self.assertFalse(manifest["functional_calibration_status_changed"])
        self.assertFalse(manifest["functional_calibration_status_upgraded"])
        self.assertEqual(manifest["threads"]["max_single_task_threads_allowed"], 40)

    def test_evaluation_tsv_carries_the_boundary_keys_for_every_row(self):
        self._seed(self._manifest_rows())
        self.module.run(self.module.build_arg_namespace([
            str(self.profiles), str(self.out), "--gates", str(self.gates_path),
        ]))
        with (self.out / "sequence_model_evaluation.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        for row in rows:
            self.assertEqual(row["sequence_claim"], "sequence-defined family homology")
            self.assertEqual(row["phenotype_boundary"], "does not establish PHB/PHA degradation")
            self.assertEqual(row["functional_calibration_status"], "not_function_calibrated")
            self.assertNotEqual(row["model_layer"], "calibrated_candidate_model")

    def test_gate_detail_tsv_has_one_row_per_model_with_stable_reason_codes(self):
        self._seed(self._manifest_rows())
        self.module.run(self.module.build_arg_namespace([
            str(self.profiles), str(self.out), "--gates", str(self.gates_path),
        ]))
        with (self.out / "sequence_model_gate_detail.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        by_id = {row["profile_id"]: row for row in rows}
        self.assertEqual(by_id[record()["profile_id"]]["blocking_reasons"], "")
        self.assertEqual(by_id[record()["profile_id"]]["heldout_rule"], "min_recall_at_or_above_five")
        self.assertEqual(
            by_id["family_DED_hfam_4_pnone"]["blocking_reasons"], "heldout_recall_below_minimum"
        )
        self.assertEqual(
            by_id["family_DED_hfam_4_pnone"]["heldout_rule"], "all_recovered_below_five"
        )
        self.assertEqual(by_id["family_DED_hfam_4_pnone"]["eff_nseq"], "0.748535")

    def test_cli_merges_external_evaluation_inputs(self):
        profile_rows = self._manifest_rows()
        for row in profile_rows:
            # Blank in the manifest: still ``pending`` until the CLI merges the
            # evaluation-input TSV in.
            for field in (
                "heldout_count", "heldout_recovered", "unexplained_confounder_hits",
                "unique_best_model", "eff_nseq",
            ):
                row[field] = ""
        evaluation_rows = [
            {
                "profile_id": record()["profile_id"],
                "heldout_count": "5",
                "heldout_recovered": "5",
                "unexplained_confounder_hits": "0",
                "unique_best_model": "true",
                "eff_nseq": "0.5",
            },
            {
                "profile_id": "family_DED_hfam_4_pnone",
                "heldout_count": "4",
                "heldout_recovered": "3",
                "unexplained_confounder_hits": "0",
                "unique_best_model": "true",
                "eff_nseq": "1.4",
            },
        ]
        evaluation = self._seed(profile_rows, evaluation_rows)
        payload = self.module.run(self.module.build_arg_namespace([
            str(self.profiles), str(self.out), "--gates", str(self.gates_path),
            "--evaluation", str(evaluation),
        ]))
        self.assertEqual(payload["counts"]["models_validated"], 1)
        by_id = {row["profile_id"]: row for row in payload["rows"]}
        self.assertEqual(by_id["family_DED_hfam_4_pnone"]["metrics"]["eff_nseq"], "1.4")
        self.assertEqual(
            by_id["family_DED_hfam_4_pnone"]["blocking_reasons"],
            ["heldout_recall_below_minimum"],
        )

    def test_evaluation_inputs_may_not_override_manifest_identity(self):
        profile_rows = self._manifest_rows()
        evaluation_rows = [
            {
                "profile_id": record()["profile_id"],
                "model_layer": "sequence_family_hmm_validated",
                "functional_calibration_status": "candidate_gate_passed_not_promoted",
                "heldout_count": "5",
                "heldout_recovered": "5",
            }
        ]
        evaluation = self._seed(profile_rows, evaluation_rows)
        payload = self.module.run(self.module.build_arg_namespace([
            str(self.profiles), str(self.out), "--gates", str(self.gates_path),
            "--evaluation", str(evaluation),
        ]))
        first = payload["rows"][0]
        self.assertEqual(first["input_model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(first["functional_calibration_status"], "not_function_calibrated")

    def test_an_unlisted_evaluation_row_fails_closed(self):
        profile_rows = self._manifest_rows()
        evaluation_rows = [{"profile_id": "not_in_the_manifest", "heldout_count": "5"}]
        evaluation = self._seed(profile_rows, evaluation_rows)
        with self.assertRaises(ValueError):
            self.module.run(self.module.build_arg_namespace([
                str(self.profiles), str(self.out), "--gates", str(self.gates_path),
                "--evaluation", str(evaluation),
            ]))


if __name__ == "__main__":
    unittest.main()
