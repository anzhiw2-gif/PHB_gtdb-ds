"""Tests for the PhaDED candidate-flow ledger and reconciliation checks (Task 8).

Written before the implementation (test-first).  The module under test closes the
candidate flow ledger for the redesign: every accession in the frozen universe
gets exactly one primary disposition, every hold/demotion/deferred reason is
recorded non-exclusively, and each documented count discrepancy of the
2026-09-21 review is explained by an explicit accession-level set difference
instead of by arithmetic on aggregates.

Load-bearing invariants:

* ``unaccounted`` and ``multiply_disposed`` are reported separately and both
  must be zero for the flow to close;
* overlapping reason sets are counted once (union), never by summing the
  per-reason counts;
* every named reconciliation check reports ``status`` and can be escalated to a
  hard failure through :func:`assert_check_ok`;
* no input path is ever treated as empty: a missing or empty file is an error.

All fixtures here are synthetic.  Nothing in this module reads or writes
``runs/``, ``results/`` or ``deploy/``.
"""

import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "reconcile_phaded_candidate_flow.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "reconcile_phaded_candidate_flow", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


module = load_module()


def record(accession, disposition, **extra):
    """Synthetic per-accession record in the module's own record shape."""
    return module.to_record(accession, disposition, **extra)


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join("" if c is None else str(c) for c in row))
    write_text(path, "\n".join(lines) + "\n")


def write_list(path, identifiers):
    write_text(path, "".join("%s\n" % i for i in identifiers))


def authority_row(layer_id, layer_kind, **extra):
    row = {
        "layer_id": layer_id,
        "layer_kind": layer_kind,
        "source_evidence_type": (
            "functional_prior"
            if layer_kind == "superfamily"
            else "sequence_clustering_2009"
        ),
        "source_version": "Knoll2009_v1.1_frozen",
        "registry_eligible": "true" if layer_kind == "superfamily" else "false",
        "gate_profile_id": "%s_%s" % (layer_kind, layer_id),
    }
    row.update(extra)
    return row


def family_row(superfamily, family_id, **extra):
    row = {"phaded_superfamily": superfamily, "phaded_family_id": family_id}
    row.update(extra)
    return row


def profile_row(profile_id, status, **extra):
    row = {"profile_id": profile_id, "model_status": status}
    row.update(extra)
    return row


def profile_id_for(layer_row):
    """The profile id the authority table binds for a layer.

    In the frozen evidence ``profile_manifest.profile_id`` equals
    ``phaded_classification_authority.gate_profile_id``, so the fixture keeps the
    same relation instead of inventing a second naming scheme.
    """
    return layer_row["gate_profile_id"]


def synthetic_profile_inventory():
    """8 superfamilies + 38 families (46 definitions) and 9 + 37 profiles."""
    superfamilies = ["SF_%d" % i for i in range(1, 9)]
    authority = [authority_row(sf, "superfamily") for sf in superfamilies]
    families = []
    for index in range(38):
        family_id = "DED_hfam_%d" % (index + 2)
        families.append(family_row(superfamilies[index % 8], family_id))
        authority.append(authority_row(family_id, "family"))
    trained_ids = {profile_id_for(row) for row in authority[:9]}
    manifest = [
        profile_row(
            profile_id_for(row),
            "trained" if profile_id_for(row) in trained_ids else "reference_only",
        )
        for row in authority
    ]
    return authority, families, manifest


class ReconcileTests(unittest.TestCase):
    def test_every_candidate_has_one_primary_disposition(self):
        rows, summary = module.reconcile(universe={"A", "B", "C"}, records=[
            record("A", "core_sequence_homolog"),
            record("B", "function_unresolved"),
            record("C", "remote_homolog_candidate"),
        ])
        self.assertEqual(summary["unaccounted"], 0)
        self.assertEqual(summary["multiply_disposed"], 0)

        self.assertEqual(summary["universe_count"], 3)
        self.assertEqual(len(rows), 3)
        self.assertEqual([row["accession"] for row in rows], ["A", "B", "C"])
        self.assertEqual(
            summary["primary_disposition_counts"],
            {
                "core_sequence_homolog": 1,
                "function_unresolved": 1,
                "remote_homolog_candidate": 1,
            },
        )

    def test_summary_reports_the_required_ledger_fields(self):
        rows, summary = module.reconcile(
            universe={"A", "B", "C", "D"},
            records=[
                record("A", "core_sequence_homolog"),
                record("B", "function_unresolved"),
                record(
                    "C",
                    "deferred_structure_review",
                    hold_reasons=["localization_conflict"],
                ),
                record(
                    "D",
                    "function_unresolved",
                    demote_reasons=["localization_conflict_signalp_predicted_export"],
                ),
            ],
        )
        for key in (
            "universe_count",
            "primary_disposition_counts",
            "unaccounted",
            "multiply_disposed",
            "demoted_from_v1",
            "hold_reason_counts_nonexclusive",
            "hold_reason_union_count",
        ):
            self.assertIn(key, summary)

        self.assertEqual(summary["hold_reason_counts_nonexclusive"],
                         {"localization_conflict": 1})
        self.assertEqual(summary["hold_reason_union_count"], 1)
        self.assertEqual(summary["demoted_from_v1"], 1)
        self.assertEqual(summary["demote_reason_counts"], {
            "localization_conflict_signalp_predicted_export": 1,
        })
        self.assertEqual(summary["disposition_partition_closes"], True)
        self.assertEqual(
            sum(summary["primary_disposition_counts"].values())
            + summary["unaccounted"],
            summary["universe_count"],
        )

    def test_hold_reasons_are_not_exclusive_but_the_union_is_counted_once(self):
        rows, summary = module.reconcile(
            universe={"A", "B"},
            records=[
                record("A", "function_unresolved",
                       hold_reasons=["localization_conflict", "architecture_criteria"]),
                record("B", "function_unresolved",
                       hold_reasons=["common_criteria"]),
            ],
        )
        self.assertEqual(
            summary["hold_reason_counts_nonexclusive"],
            {
                "architecture_criteria": 1,
                "common_criteria": 1,
                "localization_conflict": 1,
            },
        )
        self.assertEqual(summary["sum_of_hold_reason_counts"], 3)
        self.assertEqual(summary["hold_reason_union_count"], 2)
        self.assertEqual(summary["unaccounted"], 0)
        self.assertEqual(summary["multiply_disposed"], 0)
        row_a = {row["accession"]: row for row in rows}["A"]
        self.assertEqual(
            row_a["hold_reasons"], "architecture_criteria;localization_conflict"
        )

    def test_missing_disposition_is_unaccounted(self):
        rows, summary = module.reconcile(
            universe={"A", "B"},
            records=[record("A", "core_sequence_homolog")],
        )
        self.assertEqual(summary["unaccounted"], 1)
        self.assertEqual(summary["unaccounted_accessions"], ["B"])
        self.assertEqual(summary["multiply_disposed"], 0)
        self.assertEqual(summary["disposition_partition_closes"], True)
        row_b = {row["accession"]: row for row in rows}["B"]
        self.assertEqual(row_b["primary_disposition"], "")

    def test_conflicting_dispositions_are_multiply_disposed(self):
        rows, summary = module.reconcile(
            universe={"A", "B"},
            records=[
                record("A", "core_sequence_homolog"),
                record("A", "remote_homolog_candidate"),
                record("B", "function_unresolved"),
            ],
        )
        self.assertEqual(summary["multiply_disposed"], 1)
        self.assertEqual(summary["multiply_disposed_accessions"], ["A"])
        self.assertEqual(summary["unaccounted"], 0)
        row_a = {row["accession"]: row for row in rows}["A"]
        self.assertEqual(row_a["primary_disposition"], "core_sequence_homolog")
        self.assertEqual(
            row_a["disposition_conflict"], "core_sequence_homolog;remote_homolog_candidate"
        )

    def test_repeated_identical_disposition_is_not_multiply_disposed(self):
        rows, summary = module.reconcile(
            universe={"A"},
            records=[
                record("A", "core_sequence_homolog"),
                record("A", "core_sequence_homolog", hold_reasons=["common_criteria"]),
            ],
        )
        self.assertEqual(summary["multiply_disposed"], 0)
        self.assertEqual(summary["unaccounted"], 0)
        self.assertEqual(summary["hold_reason_union_count"], 1)
        row_a = rows[0]
        self.assertEqual(row_a["record_count"], 2)

    def test_records_outside_the_universe_are_rejected(self):
        with self.assertRaises(ValueError):
            module.reconcile(
                universe={"A"},
                records=[record("B", "core_sequence_homolog")],
            )

    def test_unknown_disposition_is_rejected(self):
        with self.assertRaises(ValueError):
            module.reconcile(
                universe={"A"},
                records=[record("A", "high_confidence")],
            )

    def test_empty_universe_is_rejected(self):
        with self.assertRaises(ValueError):
            module.reconcile(universe=set(), records=[])

    def test_duplicate_universe_accessions_are_rejected(self):
        with self.assertRaises(ValueError):
            module.reconcile(universe=["A", "A"], records=[
                record("A", "core_sequence_homolog"),
            ])

    def test_records_accept_plain_mappings(self):
        rows, summary = module.reconcile(
            universe={"A"},
            records=[{"accession": "A", "primary_disposition": "core_sequence_homolog"}],
        )
        self.assertEqual(summary["unaccounted"], 0)
        self.assertEqual(rows[0]["primary_disposition"], "core_sequence_homolog")


class ReasonSetSummaryTests(unittest.TestCase):
    def test_overlapping_reason_sets_are_not_added_as_unique_proteins(self):
        summary = module.summarize_reason_sets({"PhaC": {"A", "B"}, "iPhaZ": {"B", "C"}})
        self.assertEqual(summary["union_count"], 3)
        self.assertEqual(summary["sum_of_reason_counts"], 4)

    def test_overlap_details_are_machine_readable(self):
        summary = module.summarize_reason_sets({"PhaC": {"A", "B"}, "iPhaZ": {"B", "C"}})
        self.assertEqual(summary["reason_counts"], {"PhaC": 2, "iPhaZ": 2})
        self.assertEqual(summary["overlap_count"], 1)
        self.assertEqual(summary["overlap_accessions"], ["B"])
        self.assertEqual(set(summary["overlap_accessions"]), {"B"})
        self.assertEqual(
            summary["pairwise_overlaps"],
            [{"reason_a": "PhaC", "reason_b": "iPhaZ", "overlap_count": 1,
              "union_count": 3}],
        )

    def test_disjoint_reason_sets_have_zero_overlap(self):
        summary = module.summarize_reason_sets({"A": {"x"}, "B": {"y"}})
        self.assertEqual(summary["union_count"], 2)
        self.assertEqual(summary["sum_of_reason_counts"], 2)
        self.assertEqual(summary["overlap_count"], 0)
        self.assertEqual(summary["overlap_accessions"], [])

    def test_empty_reason_sets_mapping_is_rejected(self):
        with self.assertRaises(ValueError):
            module.summarize_reason_sets({})

    def test_identical_reason_sets_are_counted_once(self):
        summary = module.summarize_reason_sets({"A": {"x", "y"}, "B": {"x", "y"}})
        self.assertEqual(summary["union_count"], 2)
        self.assertEqual(summary["sum_of_reason_counts"], 4)
        self.assertEqual(summary["overlap_count"], 2)


class OverlapMatrixTests(unittest.TestCase):
    def test_overlap_matrix_rows_use_unions_not_sums(self):
        rows = module.overlap_matrix({"PhaC": {"A", "B"}, "iPhaZ": {"B", "C"}})
        self.assertEqual(len(rows), 1)
        pair = rows[0]
        self.assertEqual((pair["reason_a"], pair["reason_b"]), ("PhaC", "iPhaZ"))
        self.assertEqual(pair["overlap_count"], 1)
        self.assertEqual(pair["union_count"], 3)
        self.assertEqual(pair["sum_of_counts"], 4)
        self.assertEqual(pair["overlap_accessions"], "B")
        self.assertEqual(pair["accessions_only_in_a"], "A")
        self.assertEqual(pair["accessions_only_in_b"], "C")

    def test_every_unordered_reason_pair_is_reported_once(self):
        sets = {"a": {"A"}, "b": {"B"}, "c": {"C"}}
        rows = module.overlap_matrix(sets)
        self.assertEqual(
            [(row["reason_a"], row["reason_b"]) for row in rows],
            [("a", "b"), ("a", "c"), ("b", "c")],
        )

    def test_self_pair_is_not_emitted(self):
        rows = module.overlap_matrix({"PhaC": {"A"}, "iPhaZ": {"A"}})
        self.assertNotIn(("PhaC", "PhaC"), [(r["reason_a"], r["reason_b"]) for r in rows])


class HoldResidualTests(unittest.TestCase):
    def test_residual_is_an_explicit_set_difference(self):
        universe = {"A", "B", "C", "D", "E"}
        high_confidence = {"A", "B"}
        hold = {"C", "D"}
        other = {"E"}
        payload = module.check_hold_residual(
            universe=universe,
            high_confidence=high_confidence,
            hold=hold,
            other_dispositions={"excluded_input_quality": other},
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["universe_count"], 5)
        self.assertEqual(payload["high_confidence_count"], 2)
        self.assertEqual(payload["hold_count"], 2)
        self.assertEqual(payload["residual_count"], 1)
        self.assertEqual(payload["residual_accessions"], ["E"])
        self.assertEqual(payload["residual_explained_count"], 1)
        self.assertEqual(payload["residual_unexplained_count"], 0)
        self.assertIn("109087", module.HOLD_RESIDUAL_EXPRESSION)

    def test_residual_bucket_counts_sum_to_the_residual(self):
        payload = module.check_hold_residual(
            universe={"A", "B", "C", "D"},
            high_confidence={"A"},
            hold={"B"},
            other_dispositions={"excluded_input_quality": {"C"}, "deferred": {"D"}},
        )
        self.assertEqual(
            sum(payload["residual_bucket_counts"].values()), payload["residual_count"]
        )

    def test_unexplained_residual_is_a_mismatch(self):
        payload = module.check_hold_residual(
            universe={"A", "B", "C", "D"},
            high_confidence={"A"},
            hold={"B"},
            other_dispositions={"excluded_input_quality": {"C"}},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["residual_unexplained_count"], 1)
        self.assertEqual(payload["residual_unexplained_accessions"], ["D"])
        with self.assertRaises(ValueError):
            module.assert_check_ok(payload)

    def test_overlapping_high_confidence_and_hold_is_a_mismatch(self):
        payload = module.check_hold_residual(
            universe={"A", "B"},
            high_confidence={"A", "B"},
            hold={"B"},
            other_dispositions={},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["overlap_high_confidence_hold_count"], 1)
        self.assertEqual(payload["overlap_high_confidence_hold_accessions"], ["B"])

    def test_residual_below_the_expected_count_is_a_mismatch(self):
        payload = module.check_hold_residual(
            universe={"A", "B", "C"},
            high_confidence={"A"},
            hold={"B"},
            other_dispositions={"unexpected": {"C"}},
            expected_residual_count=616,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["expected_residual_count"], 616)
        self.assertEqual(payload["residual_count"], 1)

    def test_accessions_outside_the_universe_are_a_mismatch(self):
        payload = module.check_hold_residual(
            universe={"A"},
            high_confidence={"A", "Z"},
            hold=set(),
            other_dispositions={},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["outside_universe_accessions"], ["Z"])


class CysDemotionTests(unittest.TestCase):
    def test_demotions_joining_hold_are_stated_as_such(self):
        payload = module.check_cys_demotion(
            demoted={"A", "B"},
            hold={"A", "B", "C"},
            high_confidence_after={"C"},
            high_confidence_before={"A", "B", "C"},
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["demoted_count"], 2)
        self.assertEqual(payload["join_hold"], True)
        self.assertEqual(payload["separate_disposition"], False)
        self.assertEqual(payload["join_hold_count"], 2)
        self.assertEqual(payload["separate_disposition_count"], 0)
        self.assertEqual(payload["hold_reason"], "localization_conflict")
        self.assertEqual(payload["demoted_accessions"], ["A", "B"])
        self.assertEqual(payload["unaccounted_accessions"], [])
        self.assertIn("join the hold set", payload["explanation"])

    def test_demotions_forming_a_separate_disposition_are_stated_as_such(self):
        payload = module.check_cys_demotion(
            demoted={"A", "B"},
            hold={"C"},
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["join_hold"], False)
        self.assertEqual(payload["separate_disposition"], True)
        self.assertEqual(payload["separate_disposition_accessions"], ["A", "B"])
        self.assertIn("separate disposition", payload["explanation"])

    def test_demotions_held_outside_the_high_confidence_layer_state_a_separate_disposition(self):
        payload = module.check_cys_demotion(
            demoted={"A", "B"},
            hold={"C", "D"},
            high_confidence_after={"C", "D"},
            high_confidence_before={"A", "B", "C", "D"},
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["join_hold"], False)
        self.assertEqual(payload["separate_disposition"], True)
        self.assertEqual(payload["separate_disposition_accessions"], ["A", "B"])
        self.assertIn("separate disposition", payload["explanation"])
    def test_a_demoted_accession_missing_from_every_layer_is_reported(self):
        payload = module.check_cys_demotion(
            demoted=["A"],
            hold=set(),
            high_confidence_after=set(),
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["separate_disposition_accessions"], ["A"])
        self.assertEqual(payload["unaccounted_count"], 1)
        self.assertEqual(payload["unaccounted_accessions"], ["A"])

    def test_a_demoted_accession_still_in_high_confidence_is_a_mismatch(self):
        payload = module.check_cys_demotion(
            demoted=["A"], hold={"A"}, high_confidence_after={"A"}
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["still_high_confidence_accessions"], ["A"])

    def test_each_demoted_accession_appears_exactly_once(self):
        payload = module.check_cys_demotion(demoted=["A", "A"], hold={"A"})
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["duplicate_demoted_accessions"], ["A"])
        self.assertEqual(payload["demoted_accessions"], ["A"])
        self.assertEqual(payload["demoted_count"], 1)

    def test_expected_demotion_count_must_match(self):
        payload = module.check_cys_demotion(
            demoted=["A"], hold={"A"}, expected_count=52
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["expected_count"], 52)
        self.assertEqual(payload["demoted_count"], 1)

    def test_empty_demotion_layer_is_a_mismatch_when_required(self):
        payload = module.check_cys_demotion(demoted=[], hold=set(), required=True)
        self.assertEqual(payload["status"], "mismatch")


class CysFunnelUnionTests(unittest.TestCase):
    def test_overlap_explains_the_union_difference(self):
        phac = {"A", "B", "C"}
        iphaz = {"C", "D"}
        payload = module.check_cys_funnel_union(phac=phac, iphaz=iphaz)
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["union_count"], 4)
        self.assertEqual(payload["sum_of_reason_counts"], 5)
        self.assertEqual(payload["overlap_count"], 1)
        self.assertEqual(payload["overlap_accessions"], ["C"])
        self.assertEqual(
            payload["sum_of_reason_counts"] - payload["overlap_count"],
            payload["union_count"],
        )

    def test_union_count_is_required_and_checked(self):
        payload = module.check_cys_funnel_union(
            phac={"A"}, iphaz={"B"}, expected_union_count=4301
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["expected_union_count"], 4301)
        self.assertEqual(payload["union_count"], 2)

    def test_funnel_payload_reports_the_universe_intersection(self):
        payload = module.check_cys_funnel_union(
            phac={"A", "B"}, iphaz={"B", "C"}, universe={"B", "C", "D"}
        )
        self.assertEqual(payload["universe_count"], 3)
        self.assertEqual(payload["funnel_in_universe_count"], 2)
        self.assertEqual(payload["funnel_outside_universe_count"], 1)
        self.assertEqual(payload["funnel_outside_universe_accessions"], ["A"])

    def test_in_universe_funnel_hits_must_be_merged_high_confidence(self):
        payload = module.check_cys_funnel_union(
            phac={"A", "B"},
            iphaz={"B"},
            universe={"A", "B"},
            merged_high_confidence={"A", "B"},
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["merged_high_confidence_count"], 2)
        self.assertEqual(payload["in_universe_not_merged_count"], 0)

    def test_an_in_universe_funnel_hit_missing_from_the_merge_is_a_mismatch(self):
        payload = module.check_cys_funnel_union(
            phac={"A", "B"},
            iphaz={"B"},
            universe={"A", "B"},
            merged_high_confidence={"A"},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["in_universe_not_merged_accessions"], ["B"])

    def test_merge_cross_check_needs_a_universe(self):
        with self.assertRaises(ValueError):
            module.check_cys_funnel_union(
                phac={"A"}, iphaz={"B"}, merged_high_confidence={"A", "B"}
            )

    def test_empty_funnel_hit_sets_are_rejected(self):
        with self.assertRaises(ValueError):
            module.check_cys_funnel_union(phac=set(), iphaz=set())


class GrodonManifestDeltaTests(unittest.TestCase):
    def test_delta_buckets_sum_to_the_difference(self):
        eligible = {"G1", "G2", "G3", "G4", "G5"}
        manifest = {"G1", "G2"}
        payload = module.check_grodon_manifest_delta(
            eligible_positives=eligible,
            manifest_positives=manifest,
            exclusion_reasons={"G3": "no_genus_control", "G4": "missing_fasta"},
            prediction_failures={"G5": "too_few_ribosomal_hits"},
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["eligible_positive_count"], 5)
        self.assertEqual(payload["manifest_positive_count"], 2)
        self.assertEqual(payload["difference"], 3)
        self.assertEqual(payload["delta_explained_count"], 3)
        self.assertEqual(payload["delta_unexplained_count"], 0)
        self.assertEqual(payload["no_genus_control_count"], 1)
        self.assertEqual(payload["missing_fasta_count"], 1)
        self.assertEqual(payload["failed_prediction_count"], 1)
        self.assertEqual(
            payload["reason_bucket_counts"],
            {
                "missing_fasta": 1,
                "no_genus_control": 1,
                "too_few_ribosomal_hits": 1,
            },
        )
        self.assertEqual(
            sum(payload["reason_bucket_counts"].values()), payload["difference"]
        )

    def test_manifest_positives_outside_the_eligible_set_are_a_mismatch(self):
        payload = module.check_grodon_manifest_delta(
            eligible_positives={"G1"},
            manifest_positives={"G1", "G2"},
            exclusion_reasons={},
            prediction_failures={},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["manifest_not_eligible_count"], 1)
        self.assertEqual(payload["manifest_not_eligible_accessions"], ["G2"])

    def test_unexplained_delta_is_a_mismatch(self):
        payload = module.check_grodon_manifest_delta(
            eligible_positives={"G1", "G2", "G3"},
            manifest_positives={"G1"},
            exclusion_reasons={"G2": "no_genus_control"},
            prediction_failures={},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["delta_unexplained_count"], 1)
        self.assertEqual(payload["delta_unexplained_accessions"], ["G3"])

    def test_double_counted_reason_is_a_mismatch(self):
        payload = module.check_grodon_manifest_delta(
            eligible_positives={"G1", "G2"},
            manifest_positives={"G1"},
            exclusion_reasons={"G2": "no_genus_control"},
            prediction_failures={"G2": "hmmsearch_failed"},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["multiply_explained_accessions"], ["G2"])

    def test_expected_counts_pin_the_plan_numbers(self):
        payload = module.check_grodon_manifest_delta(
            eligible_positives={"G1", "G2", "G3"},
            manifest_positives={"G1"},
            exclusion_reasons={"G2": "no_genus_control"},
            prediction_failures={"G3": "too_few_ribosomal_hits"},
            expected_eligible_count=4507,
            expected_manifest_count=4441,
            expected_difference=66,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["count_expectation_failures"],
                         ["eligible_positive_count", "manifest_positive_count",
                          "difference"])


class MotifCoverageTests(unittest.TestCase):
    def test_coverage_and_resolvability_are_two_distinct_numbers(self):
        payload = module.check_motif_coverage_vs_resolvability(
            covered={"A", "B", "C", "D"},
            resolvable={"A", "B"},
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["covered_count"], 4)
        self.assertEqual(payload["resolvable_count"], 2)
        self.assertNotEqual(payload["covered_count"], payload["resolvable_count"])
        self.assertEqual(payload["covered_but_unresolvable_accessions"], ["C", "D"])
        self.assertEqual(payload["covered_but_unresolvable_count"], 2)
        self.assertEqual(payload["resolvable_but_uncovered_count"], 0)

    def test_resolvable_but_uncovered_is_a_mismatch(self):
        payload = module.check_motif_coverage_vs_resolvability(
            covered={"A"},
            resolvable={"A", "B"},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["resolvable_but_uncovered_accessions"], ["B"])

    def test_covered_but_unresolvable_can_be_tolerated_explicitly(self):
        payload = module.check_motif_coverage_vs_resolvability(
            covered={"A", "B"},
            resolvable={"A"},
            covered_but_unresolvable_is_expected=True,
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["covered_but_unresolvable_count"], 1)

    def test_per_criterion_rows_are_reported(self):
        payload = module.check_motif_coverage_vs_resolvability(
            covered={"A", "B"},
            resolvable={"A"},
            criteria_resolved={"lipase_box": {"A"}, "his": set()},
            covered_but_unresolvable_is_expected=True,
        )
        self.assertEqual(payload["criteria"]["lipase_box"]["resolvable_count"], 1)
        self.assertEqual(payload["criteria"]["his"]["resolvable_count"], 0)
        self.assertEqual(payload["criteria"]["his"]["unresolvable_count"], 2)

    def test_a_criterion_resolved_outside_the_coverage_layer_is_a_mismatch(self):
        payload = module.check_motif_coverage_vs_resolvability(
            covered={"A"},
            resolvable={"A"},
            criteria_resolved={"lipase_box": {"A", "B"}},
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["criteria"]["lipase_box"]["uncovered_count"], 1)

    def test_empty_covered_layer_is_rejected(self):
        with self.assertRaises(ValueError):
            module.check_motif_coverage_vs_resolvability(covered=set(), resolvable=set())


class ProfileInventoryTests(unittest.TestCase):
    def test_inventory_totals_reconcile(self):
        authority, families, manifest = synthetic_profile_inventory()
        payload = module.check_profile_inventory(
            authority_rows=authority,
            family_rows=families,
            profile_rows=manifest,
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["superfamily_count"], 8)
        self.assertEqual(payload["family_count"], 38)
        self.assertEqual(payload["definition_count"], 46)
        self.assertEqual(payload["profile_count"], 46)
        self.assertEqual(payload["trained_count"], 9)
        self.assertEqual(payload["reference_only_count"], 37)
        self.assertEqual(payload["trained_plus_reference_only"], 46)
        self.assertEqual(payload["families_bound_to_a_superfamily_count"], 38)

    def test_wrong_superfamily_count_is_a_mismatch(self):
        authority, families, manifest = synthetic_profile_inventory()
        without_one_superfamily = [
            row
            for row in authority
            if not (row["layer_kind"] == "superfamily" and row["layer_id"] == "SF_8")
        ]
        payload = module.check_profile_inventory(
            authority_rows=without_one_superfamily,
            family_rows=families,
            profile_rows=manifest,
            expected_superfamilies=8,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["superfamily_count"], 7)
        self.assertEqual(
            payload["families_with_unknown_superfamily"],
            ["DED_hfam_17", "DED_hfam_25", "DED_hfam_33", "DED_hfam_9"],
        )

    def test_wrong_family_count_is_a_mismatch(self):
        authority, families, manifest = synthetic_profile_inventory()
        payload = module.check_profile_inventory(
            authority_rows=authority,
            family_rows=families[:-1],
            profile_rows=manifest,
            expected_families=38,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["family_count"], 38)
        self.assertEqual(payload["family_definition_row_count"], 37)
        self.assertEqual(
            sum(payload["definition_counts"].values()), payload["definition_count"]
        )

    def test_authority_family_rows_must_match_the_family_definition_table(self):
        authority, families, manifest = synthetic_profile_inventory()
        payload = module.check_profile_inventory(
            authority_rows=[
                row
                for row in authority
                if not (row["layer_kind"] == "family" and row["layer_id"] == "DED_hfam_39")
            ],
            family_rows=families,
            profile_rows=manifest,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["family_count"], 37)
        self.assertEqual(payload["family_definition_row_count"], 38)
        self.assertEqual(payload["trained_count"], 9)
        self.assertEqual(payload["reference_only_count"], 37)

    def test_profile_and_definition_totals_must_be_equal(self):
        authority, families, manifest = synthetic_profile_inventory()
        payload = module.check_profile_inventory(
            authority_rows=authority,
            family_rows=families,
            profile_rows=manifest[:-1],
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["profile_count"], 45)
        self.assertEqual(payload["definition_count"], 46)
        self.assertEqual(payload["profile_definition_mismatch"], 1)

    def test_unknown_model_status_is_a_mismatch(self):
        authority, families, manifest = synthetic_profile_inventory()
        manifest[0] = profile_row(profile_id_for(authority[0]), "candidate_gate_passed")
        payload = module.check_profile_inventory(
            authority_rows=authority,
            family_rows=families,
            profile_rows=manifest,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["unknown_model_statuses"], ["candidate_gate_passed"])

    def test_family_bound_to_an_unknown_superfamily_is_a_mismatch(self):
        authority, families, manifest = synthetic_profile_inventory()
        families[0] = family_row("superfamily_that_does_not_exist", "DED_hfam_2")
        payload = module.check_profile_inventory(
            authority_rows=authority,
            family_rows=families,
            profile_rows=manifest,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(
            payload["families_with_unknown_superfamily"], ["DED_hfam_2"]
        )

    def test_superfamilies_and_families_both_have_to_be_present(self):
        authority, families, manifest = synthetic_profile_inventory()
        payload = module.check_profile_inventory(
            authority_rows=[row for row in authority if row["layer_kind"] == "family"],
            family_rows=families,
            profile_rows=manifest,
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["superfamily_count"], 0)

    def test_empty_family_table_is_rejected(self):
        authority, families, manifest = synthetic_profile_inventory()
        with self.assertRaises(ValueError):
            module.check_profile_inventory(
                authority_rows=authority,
                family_rows=[],
                profile_rows=manifest,
            )

    def test_missing_required_columns_are_rejected(self):
        with self.assertRaises(ValueError):
            module.check_profile_inventory(
                authority_rows=[{"layer_id": "x"}],
                family_rows=[],
                profile_rows=[],
            )


class CustomCountCheckTests(unittest.TestCase):
    def test_custom_counts_match(self):
        payload = module.check_custom_counts({"his_positive": 7317}, {"his_positive": 7317})
        self.assertEqual(payload["status"], "ok")

    def test_custom_counts_mismatch(self):
        payload = module.check_custom_counts({"held": 71860}, {"held": 71850})
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["mismatches"][0]["name"], "held")
        self.assertEqual(payload["mismatches"][0]["difference"], 10)

    def test_unpinned_name_is_a_mismatch(self):
        payload = module.check_custom_counts({"held": 71860}, {})
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["unpinned_names"], ["held"])


class AssertCheckOkTests(unittest.TestCase):
    def test_mismatch_raises_with_the_check_name(self):
        payload = module.check_motif_coverage_vs_resolvability(
            covered={"A"}, resolvable={"A", "B"}
        )
        with self.assertRaises(ValueError) as ctx:
            module.assert_check_ok(payload)
        self.assertIn("check_motif_coverage_vs_resolvability", str(ctx.exception))

    def test_ok_payload_does_not_raise(self):
        payload = module.check_custom_counts({"x": 1}, {"x": 1})
        module.assert_check_ok(payload)

    def test_every_named_check_carries_a_status_and_an_explanation(self):
        payloads = [
            module.check_hold_residual(
                universe={"A"}, high_confidence={"A"}, hold=set(), other_dispositions={}
            ),
            module.check_cys_demotion(demoted=["A"], hold={"A"}),
            module.check_cys_funnel_union(phac={"A"}, iphaz={"A"}),
            module.check_grodon_manifest_delta(
                eligible_positives={"A"},
                manifest_positives={"A"},
                exclusion_reasons={},
                prediction_failures={},
            ),
            module.check_motif_coverage_vs_resolvability(
                covered={"A"}, resolvable={"A"}
            ),
            module.check_profile_inventory(
                authority_rows=synthetic_profile_inventory()[0],
                family_rows=synthetic_profile_inventory()[1],
                profile_rows=synthetic_profile_inventory()[2],
            ),
        ]
        for payload in payloads:
            self.assertIn(payload["status"], {"ok", "mismatch"}, payload["check"])
            self.assertTrue(payload["explanation"], payload["check"])
            self.assertNotIn("problems", payload.keys() - {"problems"})

    def test_every_check_payload_is_json_serialisable(self):
        payload = module.check_cys_funnel_union(phac={"A"}, iphaz={"B"})
        self.assertIn('"union_count"', json.dumps(payload))

    def test_every_named_check_reports_mismatch_rather_than_raising(self):
        payloads = [
            module.check_hold_residual(
                universe={"A", "B"},
                high_confidence={"A"},
                hold=set(),
                other_dispositions={},
            ),
            module.check_cys_demotion(demoted=["A"], hold=set(), high_confidence_after={"A"}),
            module.check_cys_funnel_union(phac={"A"}, iphaz={"A"}, expected_union_count=9),
            module.check_grodon_manifest_delta(
                eligible_positives={"A", "B"},
                manifest_positives={"A"},
                exclusion_reasons={},
                prediction_failures={},
            ),
            module.check_motif_coverage_vs_resolvability(covered={"A"}, resolvable={"B"}),
            module.check_profile_inventory(
                authority_rows=[authority_row("SF_1", "superfamily")],
                family_rows=[family_row("SF_1", "DED_hfam_2")],
                profile_rows=[profile_row("superfamily_SF_1", "trained")],
            ),
        ]
        for payload in payloads:
            self.assertEqual(payload["status"], "mismatch", payload["check"])
            with self.assertRaises(ValueError):
                module.assert_check_ok(payload)


class ExpectationRoutingTests(unittest.TestCase):
    def test_expectations_parse(self):
        pins = module.parse_expectations(["held=71860", " merged = 40690 "])
        self.assertEqual(pins, {"held": 71860, "merged": 40690})

    def test_expectation_without_a_value_is_rejected(self):
        with self.assertRaises(ValueError):
            module.parse_expectations(["held"])

    def test_non_integer_expectation_is_rejected(self):
        with self.assertRaises(ValueError):
            module.parse_expectations(["held=many"])

    def test_routing_names_the_owning_check(self):
        routed = module.expectation_key_map({"cys_demoted_count": 52})
        self.assertEqual(routed, {"check_cys_demotion": {"cys_demoted_count": 52}})

    def test_unknown_expectation_name_is_rejected(self):
        with self.assertRaises(ValueError):
            module.expectation_key_map({"not_a_real_count": 1})

    def test_observed_counts_are_namespaced_per_check(self):
        payload = module.check_cys_funnel_union(phac={"A", "B"}, iphaz={"B"})
        observed = module.observed_check_counts(payload)
        self.assertEqual(
            observed,
            {"cys_funnel_union_count": 2, "cys_funnel_sum_of_reason_counts": 3},
        )


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def test_missing_path_is_an_error(self):
        with self.assertRaises(ValueError):
            module.read_accession_list(self.tmp / "absent.tsv")

    def test_directory_is_not_a_regular_file(self):
        with self.assertRaises(ValueError):
            module.read_accession_list(self.tmp)

    def test_empty_file_is_an_error(self):
        path = self.tmp / "empty.tsv"
        write_text(path, "")
        with self.assertRaises(ValueError):
            module.read_accession_list(path)

    def test_header_only_file_is_an_error(self):
        path = self.tmp / "header_only.tsv"
        write_text(path, "accession\n")
        with self.assertRaises(ValueError):
            module.read_accession_list(path)

    def test_headerless_single_column_list(self):
        path = self.tmp / "list.txt"
        write_list(path, ["A", "B", "B"])
        self.assertEqual(module.read_accession_list(path), ["A", "B"])

    def test_headerless_tsv_is_detected(self):
        path = self.tmp / "universe.tsv"
        write_text(path, "A\nB\n")
        self.assertEqual(module.read_accession_list(path), ["A", "B"])

    def test_named_column_is_used(self):
        path = self.tmp / "universe.tsv"
        write_tsv(path, ["accession", "genome"], [["A", "G1"], ["B", "G2"]])
        self.assertEqual(module.read_accession_list(path), ["A", "B"])

    def test_missing_required_column_is_an_error(self):
        path = self.tmp / "bad.tsv"
        write_tsv(path, ["protein_name"], [["A"]])
        with self.assertRaises(ValueError):
            module.read_records_table(path)

    def test_records_table_reads_the_full_ledger_schema(self):
        path = self.tmp / "records.tsv"
        write_tsv(
            path,
            ["accession", "primary_disposition", "hold_reasons", "demote_reasons",
             "deferred", "disposition_basis"],
            [["A", "function_unresolved", "common_criteria;architecture_criteria",
              "localization_conflict_signalp_predicted_export", "false", "test"]],
        )
        records = module.read_records_table(path)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["accession"], "A")
        self.assertEqual(records[0]["hold_reasons"],
                         ["common_criteria", "architecture_criteria"])
        self.assertEqual(records[0]["defer_reasons"], [])

    def test_reason_table_groups_accessions_by_reason(self):
        path = self.tmp / "holds.tsv"
        write_tsv(
            path,
            ["accession", "hold_reason"],
            [["A", "common_criteria"], ["B", "common_criteria"], ["C", "architecture"]],
        )
        grouped = module.read_reason_table(path, "hold_reason")
        self.assertEqual(grouped, {"common_criteria": ["A", "B"], "architecture": ["C"]})
        self.assertEqual(module.read_accession_list(path), ["A", "B", "C"])


class CliTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def build_inputs(self, with_optional=True):
        inputs = self.tmp / "inputs"
        universe = inputs / "universe.tsv"
        write_tsv(universe, ["accession"], [["A"], ["B"], ["C"], ["D"]])
        records = inputs / "records.tsv"
        write_tsv(
            records,
            ["accession", "primary_disposition", "hold_reasons", "demote_reasons"],
            [
                ["A", "core_sequence_homolog", "", ""],
                ["B", "remote_homolog_candidate", "common_criteria", ""],
                ["C", "deferred_structure_review",
                 "localization_conflict", "localization_conflict_signalp_predicted_export"],
                ["D", "function_unresolved", "architecture_criteria", ""],
            ],
        )
        paths = {"universe": universe, "records": records}
        if with_optional:
            hold = inputs / "hold.tsv"
            write_tsv(
                hold,
                ["accession", "hold_reason"],
                [["B", "common_criteria"], ["C", "localization_conflict"],
                 ["D", "architecture_criteria"]],
            )
            demoted = inputs / "demoted.tsv"
            write_tsv(
                demoted,
                ["accession", "demote_reason"],
                [["C", "localization_conflict_signalp_predicted_export"]],
            )
            high_confidence = inputs / "high_confidence.tsv"
            write_tsv(high_confidence, ["accession"], [["A"]])
            v1 = inputs / "v1.tsv"
            write_tsv(v1, ["accession"], [["A"], ["C"]])
            paths.update(
                hold=hold,
                demoted=demoted,
                high_confidence=high_confidence,
                v1=v1,
            )
        return paths

    def run_cli(self, argv):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(io.StringIO()):
            code = module.main(argv)
        return code, buffer.getvalue()

    def test_cli_writes_the_three_required_outputs(self):
        paths = self.build_inputs()
        out_dir = self.tmp / "out"
        code, printed = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(paths["records"]),
            "--output-dir", str(out_dir),
        ])
        self.assertEqual(code, 0)
        payload = json.loads(printed)
        self.assertEqual(payload["unaccounted"], 0)
        self.assertEqual(payload["multiply_disposed"], 0)
        self.assertTrue((out_dir / "candidate_flow.tsv").is_file())
        self.assertTrue((out_dir / "candidate_flow_summary.json").is_file())
        self.assertTrue((out_dir / "overlap_matrix.tsv").is_file())
        summary = json.loads((out_dir / "candidate_flow_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["universe_count"], 4)
        self.assertEqual(summary["primary_disposition_counts"], {
            "core_sequence_homolog": 1,
            "deferred_structure_review": 1,
            "function_unresolved": 1,
            "remote_homolog_candidate": 1,
        })
        flow_lines = (out_dir / "candidate_flow.tsv").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(flow_lines), 5)

    def test_cli_exits_non_zero_when_the_partition_does_not_close(self):
        paths = self.build_inputs(with_optional=False)
        records = paths["records"]
        write_tsv(
            records,
            ["accession", "primary_disposition"],
            [["A", "core_sequence_homolog"], ["A", "remote_homolog_candidate"]],
        )
        out_dir = self.tmp / "out_bad"
        code, printed = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(records),
            "--output-dir", str(out_dir),
        ])
        self.assertEqual(code, 2)
        payload = json.loads(printed)
        self.assertEqual(payload["multiply_disposed"], 1)
        self.assertGreater(payload["unaccounted"], 0)

    def test_cli_runs_the_named_checks_that_were_given_inputs(self):
        paths = self.build_inputs()
        out_dir = self.tmp / "out_checks"
        code, printed = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(paths["records"]),
            "--output-dir", str(out_dir),
            "--hold-reasons", str(paths["hold"]),
            "--cys-demoted", str(paths["demoted"]),
            "--high-confidence", str(paths["high_confidence"]),
            "--v1-high-confidence", str(paths["v1"]),
        ])
        self.assertEqual(code, 0)
        payload = json.loads(printed)
        self.assertIn("check_hold_residual", payload["checks"])
        self.assertIn("check_cys_demotion", payload["checks"])
        self.assertEqual(payload["checks"]["check_hold_residual"]["status"], "ok")
        self.assertTrue((out_dir / "check_hold_residual.json").is_file())
        self.assertTrue((out_dir / "checks_summary.json").is_file())
        residual = json.loads(
            (out_dir / "check_hold_residual.json").read_text(encoding="utf-8")
        )
        self.assertEqual(residual["residual_count"], 0)

    def test_cli_reports_skipped_checks_and_never_invents_empty_inputs(self):
        paths = self.build_inputs(with_optional=False)
        out_dir = self.tmp / "out_skipped"
        code, printed = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(paths["records"]),
            "--output-dir", str(out_dir),
        ])
        self.assertEqual(code, 0)
        payload = json.loads(printed)
        self.assertIn("check_profile_inventory", payload["checks_skipped"])
        for reason in payload["checks_skipped"].values():
            self.assertTrue(reason)

    def test_cli_fails_on_a_mismatching_check(self):
        paths = self.build_inputs()
        out_dir = self.tmp / "out_mismatch"
        code, printed = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(paths["records"]),
            "--output-dir", str(out_dir),
            "--hold-reasons", str(paths["hold"]),
            "--cys-demoted", str(paths["demoted"]),
            "--high-confidence", str(paths["high_confidence"]),
            "--expect", "hold_residual_count=616",
        ])
        self.assertEqual(code, 1)
        payload = json.loads(printed)
        self.assertEqual(
            payload["checks"]["check_hold_residual"]["status"], "mismatch"
        )

    def test_cli_fails_when_an_expectation_names_a_check_that_did_not_run(self):
        paths = self.build_inputs(with_optional=False)
        out_dir = self.tmp / "out_unpinned"
        code, printed = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(paths["records"]),
            "--output-dir", str(out_dir),
            "--expect", "cys_funnel_union_count=4301",
        ])
        self.assertEqual(code, 1)
        payload = json.loads(printed)
        self.assertEqual(payload["checks"]["check_custom_counts"]["status"], "mismatch")
        self.assertIn("check_cys_funnel_union", payload["checks_skipped"])

    def test_cli_refuses_a_missing_input_path(self):
        paths = self.build_inputs()
        with self.assertRaises(ValueError):
            self.run_cli([
                "--universe", str(paths["universe"]),
                "--records", str(self.tmp / "inputs" / "absent.tsv"),
                "--output-dir", str(self.tmp / "out_absent"),
            ])

    def test_cli_refuses_an_existing_non_empty_output_dir(self):
        paths = self.build_inputs()
        out_dir = self.tmp / "out_full"
        write_text(out_dir / "keep.txt", "evidence\n")
        with self.assertRaises(ValueError):
            self.run_cli([
                "--universe", str(paths["universe"]),
                "--records", str(paths["records"]),
                "--output-dir", str(out_dir),
            ])

    def test_cli_reuses_an_existing_empty_output_dir(self):
        paths = self.build_inputs()
        out_dir = self.tmp / "out_empty"
        out_dir.mkdir(parents=True)
        code, _ = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(paths["records"]),
            "--output-dir", str(out_dir),
        ])
        self.assertEqual(code, 0)

    def test_cli_writes_the_overlap_matrix_from_real_reason_sets(self):
        paths = self.build_inputs()
        out_dir = self.tmp / "out_overlap"
        write_list(self.tmp / "inputs" / "phac.txt", ["A", "B"])
        write_list(self.tmp / "inputs" / "iphaz.txt", ["B", "C"])
        code, printed = self.run_cli([
            "--universe", str(paths["universe"]),
            "--records", str(paths["records"]),
            "--output-dir", str(out_dir),
            "--phac-hits", str(self.tmp / "inputs" / "phac.txt"),
            "--iphaz-hits", str(self.tmp / "inputs" / "iphaz.txt"),
            "--expect", "cys_funnel_union_count=3",
        ])
        self.assertEqual(code, 0)
        payload = json.loads(printed)
        self.assertEqual(
            payload["checks"]["check_cys_funnel_union"]["union_count"], 3
        )
        lines = (out_dir / "overlap_matrix.tsv").read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        self.assertEqual(header[0], "reason_a")
        funnel_rows = [
            line.split("\t")
            for line in lines[1:]
            if line.startswith("funnel_PhaC\tfunnel_iPhaZ\t")
        ]
        self.assertEqual(len(funnel_rows), 1)
        self.assertEqual(funnel_rows[0][header.index("overlap_count")], "1")
        self.assertEqual(funnel_rows[0][header.index("union_count")], "3")
        self.assertEqual(funnel_rows[0][header.index("overlap_accessions")], "B")

    def test_cli_requires_the_core_arguments(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit):
                module.main(["--universe", "x"])
        self.assertIn("usage", stderr.getvalue())

    def test_cli_preflight_rejects_a_directory_as_universe(self):
        with self.assertRaises(ValueError):
            module.main([
                "--universe", str(self.tmp),
                "--records", str(self.tmp / "nope.tsv"),
                "--output-dir", str(self.tmp / "out_dir_universe"),
            ])


class UniverseCoverageTests(unittest.TestCase):
    def test_universe_flow_rows_cover_every_accession(self):
        universe = {"A", "B", "C", "D"}
        records = [
            record("A", "core_sequence_homolog"),
            record("B", "function_unresolved", hold_reasons=["common_criteria"]),
            record("C", "deferred_structure_review", defer_reasons=["with_lipase_deferred"]),
            record("D", "excluded_input_quality"),
        ]
        rows, summary = module.reconcile(universe=universe, records=records)
        self.assertEqual({row["accession"] for row in rows}, universe)
        self.assertEqual(summary["unaccounted"], 0)
        self.assertEqual(summary["multiply_disposed"], 0)
        self.assertEqual(summary["deferred_reason_counts"], {"with_lipase_deferred": 1})
        self.assertEqual(summary["deferred_reason_union_count"], 1)

    def test_reconcile_reports_hold_and_demote_overlap(self):
        rows, summary = module.reconcile(
            universe={"A"},
            records=[
                record(
                    "A",
                    "function_unresolved",
                    hold_reasons=["localization_conflict"],
                    demote_reasons=["localization_conflict_signalp_predicted_export"],
                )
            ],
        )
        self.assertEqual(summary["demoted_from_v1"], 1)
        self.assertEqual(summary["hold_reason_union_count"], 1)
        self.assertEqual(summary["hold_and_demote_overlap_accessions"], ["A"])
        self.assertEqual(rows[0]["hold_or_demote"], "hold+demote")


if __name__ == "__main__":
    unittest.main()
