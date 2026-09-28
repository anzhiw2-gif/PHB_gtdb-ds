"""Superfamily resolution must be evaluable, and the family gate must not move (rule #4).

``finalize_phaded_reference_panel_calibration.py`` L75 gates on ``kind == "family"``,
so every ``profile_kind == "superfamily"`` row is structurally excluded even when its
evidence is complete -- contradicting the layering authority
(``20260917_phaded_classification_authority_01``: 8 superfamilies are the functional
prior frame with ``registry_eligible=true``; 38 families are the refinement layer).

These tests pin the corrected convention before it is implemented:

* ``profile_kind == "superfamily"`` gets its **own** explicit gate with the same
  numeric thresholds (3 / 1 / 1 / 1) plus the zero-unexplained-hit condition, read
  from explicitly named columns whose absence **fails closed** (blocks, never passes);
* the family gate keeps its exact expression, so re-running it on the frozen readiness
  input must still reproduce the frozen decision file byte-for-byte;
* axis-anchored counts without an audited anchor basis are refused outright;
* ``new_family_call_created`` other than false still raises.
"""

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "pipeline" / "scripts" / "finalize_phaded_reference_panel_calibration.py"

READINESS = (
    REPO / "runs" / "20260915_phaded_reference_panel_acquisition_02" / "results"
    / "amendment" / "profile_calibration_readiness.tsv"
)
FROZEN_DECISIONS = (
    REPO / "runs" / "20260915_phaded_reference_panel_acquisition_02" / "results"
    / "calibration" / "leaveout_calibration_decisions.tsv"
)
FROZEN_DECISIONS_SHA256 = "4a04d5eabb9c1ed285349cfe669745f160102ba247c8ff97b620c5aabb228c93"

BASE_FIELDS = [
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "existing_experimental_positive", "external_bound_positive", "external_bound_negative",
    "external_bound_challenge", "heldout_positive", "family_resolved_negative_status",
    "challenge_status", "calibration_decision", "new_family_call_created",
]


def load_module():
    spec = importlib.util.spec_from_file_location("finalize_phaded_reference_panel_calibration", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def family_row(**overrides) -> dict:
    row = {
        "profile_id": "family_DED_hfam_X_deadbeef",
        "profile_kind": "family",
        "phaded_superfamily": "sf_x",
        "phaded_family_id": "DED_hfam_X",
        "existing_experimental_positive": "0",
        "external_bound_positive": "0",
        "external_bound_negative": "0",
        "external_bound_challenge": "0",
        "heldout_positive": "0",
        "family_resolved_negative_status": "missing",
        "challenge_status": "missing",
        "calibration_decision": "reference_only_insufficient_panel",
        "new_family_call_created": "false",
    }
    row.update(overrides)
    return row


def superfamily_row(**overrides) -> dict:
    row = {
        "profile_id": "superfamily_sf_x",
        "profile_kind": "superfamily",
        "phaded_superfamily": "sf_x",
        "phaded_family_id": "",
        "existing_experimental_positive": "0",
        "external_bound_positive": "0",
        "external_bound_negative": "0",
        "external_bound_challenge": "0",
        "heldout_positive": "0",
        "family_resolved_negative_status": "missing",
        "challenge_status": "missing",
        "calibration_decision": "planned_not_run_superfamily_requires_family_resolution",
        "new_family_call_created": "false",
        # explicitly named superfamily-resolution columns
        "superfamily_positive_count": "0",
        "superfamily_heldout_positive": "0",
        "axis_anchored_negative_count": "0",
        "axis_anchored_challenge_count": "0",
        "axis_anchor_basis_status": "not_audited",
        "zero_unexplained_hits_status": "not_assessable",
    }
    row.update(overrides)
    return row


def run_finalize(rows, tmp: Path, *, extra_fields=(), promote_calibrated_model=False):
    module = load_module()
    fields = BASE_FIELDS + [field for field in extra_fields if field not in BASE_FIELDS]
    readiness = tmp / "readiness.tsv"
    with readiness.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result = module.finalize(
        readiness, tmp / "out", run_id="test",
        promote_calibrated_model=promote_calibrated_model,
    )
    return module, result, tmp / "out"


SUPERFAMILY_EXTRA_FIELDS = (
    "superfamily_positive_count", "superfamily_heldout_positive",
    "axis_anchored_negative_count", "axis_anchored_challenge_count",
    "axis_anchor_basis_status", "zero_unexplained_hits_status",
)

#: Added by the 2026-09-28 evidence-model redesign (plan Task 5): the family
#: gate counts qualified independent evidence, and the two columns must be
#: present in the readiness row to be read at all.
FAMILY_QUALIFICATION_FIELDS = (
    "qualified_e2_e3_positive_count", "independent_genus_count",
)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SuperfamilyGateTests(unittest.TestCase):
    """profile_kind == superfamily must be evaluable and must fail closed."""

    def _decision(self, rows, tmp, profile_id="superfamily_sf_x"):
        _module, result, _out = run_finalize(rows, tmp, extra_fields=SUPERFAMILY_EXTRA_FIELDS)
        return {row["profile_id"]: row for row in result["decisions"]}[profile_id]

    def test_superfamily_gate_blocked_when_positive_insufficient(self):
        with tempfile.TemporaryDirectory() as temporary:
            row = self._decision(
                [
                    superfamily_row(
                        superfamily_positive_count="2",
                        superfamily_heldout_positive="1",
                        axis_anchored_negative_count="1",
                        axis_anchored_challenge_count="1",
                        axis_anchor_basis_status="audited",
                        zero_unexplained_hits_status="sufficient",
                    )
                ],
                Path(temporary),
            )
            self.assertEqual(row["decision"], "planned_not_run_superfamily_requires_family_resolution")

    def test_superfamily_gate_blocked_when_heldout_missing(self):
        with tempfile.TemporaryDirectory() as temporary:
            row = self._decision(
                [
                    superfamily_row(
                        superfamily_positive_count="3",
                        superfamily_heldout_positive="0",
                        axis_anchored_negative_count="1",
                        axis_anchored_challenge_count="1",
                        axis_anchor_basis_status="audited",
                        zero_unexplained_hits_status="sufficient",
                    )
                ],
                Path(temporary),
            )
            self.assertEqual(row["decision"], "planned_not_run_superfamily_requires_family_resolution")

    def test_superfamily_gate_blocked_when_axis_credit_missing(self):
        with tempfile.TemporaryDirectory() as temporary:
            row = self._decision(
                [
                    superfamily_row(
                        superfamily_positive_count="3",
                        superfamily_heldout_positive="1",
                        axis_anchored_negative_count="0",
                        axis_anchored_challenge_count="1",
                        axis_anchor_basis_status="audited",
                        zero_unexplained_hits_status="sufficient",
                    )
                ],
                Path(temporary),
            )
            self.assertEqual(row["decision"], "planned_not_run_superfamily_requires_family_resolution")

    def test_superfamily_gate_blocked_when_zero_unexplained_hits_not_assessable(self):
        with tempfile.TemporaryDirectory() as temporary:
            row = self._decision(
                [
                    superfamily_row(
                        superfamily_positive_count="3",
                        superfamily_heldout_positive="1",
                        axis_anchored_negative_count="1",
                        axis_anchored_challenge_count="1",
                        axis_anchor_basis_status="audited",
                        zero_unexplained_hits_status="not_assessable",
                    )
                ],
                Path(temporary),
            )
            self.assertEqual(row["decision"], "planned_not_run_superfamily_requires_family_resolution")

    def test_superfamily_gate_blocked_when_axis_basis_is_not_audited(self):
        with tempfile.TemporaryDirectory() as temporary:
            row = self._decision(
                [
                    superfamily_row(
                        superfamily_positive_count="3",
                        superfamily_heldout_positive="1",
                        axis_anchored_negative_count="0",
                        axis_anchored_challenge_count="0",
                        axis_anchor_basis_status="not_audited",
                        zero_unexplained_hits_status="sufficient",
                    )
                ],
                Path(temporary),
            )
            self.assertEqual(row["decision"], "planned_not_run_superfamily_requires_family_resolution")

    def test_superfamily_gate_fails_closed_when_columns_are_absent(self):
        """Absent new columns must block, never pass."""
        with tempfile.TemporaryDirectory() as temporary:
            row = self._decision(
                [
                    superfamily_row(
                        existing_experimental_positive="30",
                        external_bound_positive="30",
                        heldout_positive="5",
                        external_bound_negative="5",
                        external_bound_challenge="5",
                        family_resolved_negative_status="sufficient",
                        challenge_status="sufficient",
                    )
                ],
                Path(temporary),
            )
            self.assertEqual(row["decision"], "planned_not_run_superfamily_requires_family_resolution")
            self.assertEqual(row["calibration_status"], "not_run")

    def test_superfamily_gate_is_reachable_when_every_slot_is_met(self):
        """The structural exclusion is removed: complete evidence now passes.

        Updated by the 2026-09-28 evidence-model redesign (plan Task 5): a
        superfamily row is still *evaluable* and its own gate still reads
        ``true``, but promotion is now a separate authorized action, so the
        default finalize reports ``candidate_gate_passed_not_promoted`` instead
        of ``calibrated_candidate_model``.  The structural exclusion this test
        originally pinned (superfamilies were previously blocked outright) is
        still gone -- the gate is evaluated and passes.
        """
        with tempfile.TemporaryDirectory() as temporary:
            module, result, out = run_finalize(
                [
                    superfamily_row(
                        superfamily_positive_count="3",
                        superfamily_heldout_positive="1",
                        axis_anchored_negative_count="1",
                        axis_anchored_challenge_count="1",
                        axis_anchor_basis_status="audited",
                        zero_unexplained_hits_status="sufficient",
                    )
                ],
                Path(temporary),
                extra_fields=SUPERFAMILY_EXTRA_FIELDS,
            )
            row = {item["profile_id"]: item for item in result["decisions"]}["superfamily_sf_x"]
            self.assertEqual(row["decision"], "candidate_gate_passed_not_promoted")
            self.assertEqual(row["calibration_status"], "candidate_gate_passed_not_promoted")
            detail = read_tsv(out / "leaveout_calibration_gate_detail.tsv")[0]
            self.assertEqual(detail["gate_resolution"], "superfamily")
            self.assertEqual(detail["gate_passed"], "true")
            self.assertEqual(detail["gate_kind"], "superfamily")
            self.assertEqual(detail["promotion_authorized"], "false")
            self.assertIn("superfamily", result["report"]["gate"])
            # The superfamily gate still passes; with promotion authorized the
            # promotion path is exercised (authorization is orthogonal to kind).
            promoted = module.finalize(
                Path(temporary) / "readiness.tsv", Path(temporary) / "out_promoted",
                run_id="test", promote_calibrated_model=True,
            )
            promoted_row = {
                item["profile_id"]: item for item in promoted["decisions"]
            }["superfamily_sf_x"]
            self.assertEqual(promoted_row["decision"], "calibrated_candidate_model")

    def test_axis_anchored_counts_without_audited_basis_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                run_finalize(
                    [
                        superfamily_row(
                            axis_anchored_negative_count="1",
                            axis_anchor_basis_status="not_audited",
                        )
                    ],
                    Path(temporary),
                    extra_fields=SUPERFAMILY_EXTRA_FIELDS,
                )

    def test_new_family_call_created_true_raises_for_both_kinds(self):
        for row in (
            family_row(new_family_call_created="true"),
            superfamily_row(new_family_call_created="true"),
        ):
            with tempfile.TemporaryDirectory() as temporary:
                with self.assertRaises(ValueError):
                    run_finalize([row], Path(temporary), extra_fields=SUPERFAMILY_EXTRA_FIELDS)


class FamilyGateRegressionTests(unittest.TestCase):
    """The family gate must not move by one byte."""

    def test_family_gate_is_byte_identical_on_frozen_readiness(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            module.finalize(READINESS, Path(temporary), run_id="regression")
            regenerated = Path(temporary) / "leaveout_calibration_decisions.tsv"
            self.assertEqual(
                hashlib.sha256(regenerated.read_bytes()).hexdigest(), FROZEN_DECISIONS_SHA256
            )
            self.assertEqual(regenerated.read_bytes(), FROZEN_DECISIONS.read_bytes())

    def test_frozen_rows_keep_their_historical_decisions(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            result = module.finalize(READINESS, Path(temporary), run_id="regression")
            decisions = {row["profile_id"]: row["decision"] for row in result["decisions"]}
            self.assertEqual(len(decisions), 37)
            self.assertEqual(sum(value == "calibrated_candidate_model" for value in decisions.values()), 0)
            self.assertEqual(
                sum(value == "planned_not_run_superfamily_requires_family_resolution" for value in decisions.values()),
                4,
            )

    def test_family_rows_ignore_axis_anchored_credit(self):
        """Axis credit is reported for family rows but must never satisfy the family gate."""
        with tempfile.TemporaryDirectory() as temporary:
            module, result, out = run_finalize(
                [
                    family_row(
                        existing_experimental_positive="3",
                        heldout_positive="1",
                        external_bound_negative="0",
                        external_bound_challenge="1",
                        family_resolved_negative_status="missing",
                        challenge_status="sufficient",
                        axis_anchored_negative_count="5",
                        axis_anchored_challenge_count="5",
                        axis_anchor_basis_status="audited",
                    )
                ],
                Path(temporary),
                extra_fields=("axis_anchored_negative_count", "axis_anchored_challenge_count", "axis_anchor_basis_status"),
            )
            row = {item["profile_id"]: item for item in result["decisions"]}[family_row()["profile_id"]]
            self.assertEqual(row["decision"], "reference_only_insufficient_panel")
            detail = read_tsv(out / "leaveout_calibration_gate_detail.tsv")[0]
            self.assertEqual(detail["gate_resolution"], "family")
            self.assertEqual(detail["gate_passed"], "false")
            self.assertEqual(detail["axis_anchored_negative_count"], "5")

    def test_family_gate_still_passes_on_the_legacy_slots(self):
        """Superseded by the 2026-09-28 evidence-model redesign (plan Task 5).

        The legacy slots alone no longer reach ``calibrated_candidate_model``:

        * the family gate now counts *qualified* independent evidence
          (``qualified_e2_e3_positive_count`` / ``independent_genus_count``)
          rather than raw ``existing_experimental_positive``, so a readiness row
          that declares only legacy slots fails closed to 0;
        * a passed-but-unpromoted gate reports
          ``candidate_gate_passed_not_promoted``; ``calibrated_candidate_model``
          is emitted only by an explicitly authorized promotion
          (``--promote-calibrated-model``).

        The frozen-file regression
        (``test_family_gate_is_byte_identical_on_frozen_readiness``) is
        unaffected and still pins the historical decisions byte-for-byte.
        """
        with tempfile.TemporaryDirectory() as temporary:
            _module, result, out = run_finalize(
                [
                    family_row(
                        existing_experimental_positive="3",
                        heldout_positive="1",
                        external_bound_negative="1",
                        external_bound_challenge="1",
                        family_resolved_negative_status="sufficient",
                        challenge_status="sufficient",
                    )
                ],
                Path(temporary),
                extra_fields=FAMILY_QUALIFICATION_FIELDS,
            )
            row = {item["profile_id"]: item for item in result["decisions"]}[family_row()["profile_id"]]
            self.assertEqual(row["decision"], "reference_only_insufficient_panel")
            detail = read_tsv(out / "leaveout_calibration_gate_detail.tsv")[0]
            self.assertEqual(detail["qualified_e2_e3_positive_count"], "0")
            self.assertEqual(detail["independent_genus_count"], "0")
            self.assertEqual(detail["promotion_authorized"], "false")

    def test_family_gate_passes_but_is_not_promoted_without_authorization(self):
        """A fully qualified family passes the gate yet stays unpromoted."""
        with tempfile.TemporaryDirectory() as temporary:
            _module, result, out = run_finalize(
                [
                    family_row(
                        heldout_positive="1",
                        external_bound_negative="1",
                        external_bound_challenge="1",
                        family_resolved_negative_status="sufficient",
                        challenge_status="sufficient",
                        qualified_e2_e3_positive_count="3",
                        independent_genus_count="3",
                    )
                ],
                Path(temporary),
                extra_fields=FAMILY_QUALIFICATION_FIELDS,
            )
            row = {item["profile_id"]: item for item in result["decisions"]}[family_row()["profile_id"]]
            self.assertEqual(row["decision"], "candidate_gate_passed_not_promoted")
            detail = read_tsv(out / "leaveout_calibration_gate_detail.tsv")[0]
            self.assertEqual(detail["gate_passed"], "true")
            self.assertEqual(detail["promotion_authorized"], "false")
            self.assertEqual(result["summary"]["calibrated_candidate_model"], 0)
            self.assertEqual(result["summary"]["candidate_gate_passed_not_promoted"], 1)

    def test_authorized_promotion_is_the_only_path_to_calibrated_candidate_model(self):
        """``calibrated_candidate_model`` requires the explicit finalize action."""
        with tempfile.TemporaryDirectory() as temporary:
            module, result, _out = run_finalize(
                [
                    family_row(
                        heldout_positive="1",
                        external_bound_negative="1",
                        external_bound_challenge="1",
                        family_resolved_negative_status="sufficient",
                        challenge_status="sufficient",
                        qualified_e2_e3_positive_count="3",
                        independent_genus_count="3",
                    )
                ],
                Path(temporary),
                extra_fields=FAMILY_QUALIFICATION_FIELDS,
                promote_calibrated_model=True,
            )
            row = {item["profile_id"]: item for item in result["decisions"]}[family_row()["profile_id"]]
            self.assertEqual(row["decision"], "calibrated_candidate_model")
            self.assertTrue(result["report"]["gate"]["promotion"]["authorized"])

    def test_gate_detail_labels_resolution_for_every_row(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            result = module.finalize(READINESS, Path(temporary), run_id="regression")
            rows = read_tsv(Path(temporary) / "leaveout_calibration_gate_detail.tsv")
            self.assertEqual(len(rows), 37)
            for row in rows:
                self.assertEqual(row["gate_resolution"], row["profile_kind"])
                self.assertIn(row["gate_resolution"], {"family", "superfamily"})
                self.assertEqual(row["gate_evaluated"], "true")
            superfamilies = [row for row in rows if row["profile_kind"] == "superfamily"]
            self.assertEqual(len(superfamilies), 4)
            for row in superfamilies:
                self.assertEqual(row["gate_passed"], "false")
            self.assertEqual(result["summary"]["calibrated_candidate_model"], 0)

    def test_gate_report_declares_both_resolutions_and_unchanged_thresholds(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            result = module.finalize(READINESS, Path(temporary), run_id="regression")
            gate = result["report"]["gate"]
            self.assertEqual(gate["family"]["resolution"], "family")
            self.assertEqual(gate["superfamily"]["resolution"], "superfamily")
            for block in (gate["family"], gate["superfamily"]):
                self.assertEqual(block["minimum_positive_count"], 3)
                self.assertEqual(block["minimum_heldout_positive_count"], 1)
                self.assertEqual(block["minimum_negative_count"], 1)
                self.assertEqual(block["minimum_challenge_count"], 1)
                self.assertTrue(block["zero_unexplained_negative_or_challenge_hits"])
                self.assertFalse(block["new_family_calls_allowed"])
            self.assertFalse(result["report"]["summary"]["new_family_call_created"])
            self.assertFalse(result["summary"]["new_family_call_created"])


if __name__ == "__main__":
    unittest.main()
