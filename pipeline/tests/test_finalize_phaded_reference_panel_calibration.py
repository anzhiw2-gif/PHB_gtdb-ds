import contextlib
import csv
import importlib.util
import io
import itertools
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "finalize_phaded_reference_panel_calibration.py"

BASE_FIELDS = [
    "profile_id", "profile_kind", "external_bound_positive", "external_bound_negative",
    "external_bound_challenge", "heldout_positive", "family_resolved_negative_status",
    "challenge_status", "calibration_decision", "new_family_call_created",
]

#: Columns added by the 2026-09-28 evidence-model redesign (Task 5).  They are absent
#: from the frozen readiness input, so absence must fail closed (counted as 0).
QUALIFIED_FIELDS = ("qualified_e2_e3_positive_count", "independent_genus_count")

SUPPLIED_DECISION = "reference_only_insufficient_panel"
NOT_PROMOTED = "candidate_gate_passed_not_promoted"
PROMOTED = "calibrated_candidate_model"

_OUTPUTS = itertools.count(1)


def _load_module():
    spec = importlib.util.spec_from_file_location("finalize_phaded_reference_panel_calibration", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def calibration_row(**overrides):
    row = {
        "profile_id": "family_DED_hfam_test_deadbeef",
        "profile_kind": "family",
        "phaded_superfamily": "extracellular dPHASCL type 1",
        "phaded_family_id": "DED_hfam_test",
        "existing_experimental_positive": "0",
        "external_bound_positive": "0",
        "external_bound_negative": "0",
        "external_bound_challenge": "0",
        "heldout_positive": "0",
        "family_resolved_negative_status": "missing",
        "challenge_status": "missing",
        "calibration_decision": SUPPLIED_DECISION,
        "new_family_call_created": "false",
    }
    row.update(overrides)
    return row


def passing_row(**overrides):
    """A family row whose functional slots are complete under the new gate."""
    row = calibration_row(
        heldout_positive="1",
        external_bound_negative="1",
        external_bound_challenge="1",
        family_resolved_negative_status="sufficient",
        challenge_status="sufficient",
        qualified_e2_e3_positive_count="3",
        independent_genus_count="3",
    )
    row.update(overrides)
    return row


def write_readiness(path, rows):
    fields = list(BASE_FIELDS)
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = {field: "" for field in fields}
            output.update(row)
            writer.writerow(output)
    return path


def run_finalize(rows, root, *, promote=None):
    module = _load_module()
    readiness = write_readiness(Path(root) / "readiness.tsv", rows)
    output_dir = Path(root) / f"out-{next(_OUTPUTS)}"
    kwargs = {} if promote is None else {"promote_calibrated_model": promote}
    result = module.finalize(readiness, output_dir, run_id="test", **kwargs)
    return module, result, output_dir


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class FinalizeReferencePanelCalibrationTests(unittest.TestCase):
    def test_all_profiles_without_panel_gate_are_planned_not_run(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            readiness = root / "readiness.tsv"
            with readiness.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["profile_id", "profile_kind", "external_bound_positive", "external_bound_negative", "external_bound_challenge", "heldout_positive", "family_resolved_negative_status", "challenge_status", "calibration_decision", "new_family_call_created"], delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows([
                    {"profile_id": "p1", "profile_kind": "family", "external_bound_positive": "0", "external_bound_negative": "0", "external_bound_challenge": "0", "heldout_positive": "0", "family_resolved_negative_status": "missing", "challenge_status": "missing", "calibration_decision": "reference_only_insufficient_panel", "new_family_call_created": "false"},
                    {"profile_id": "p2", "profile_kind": "superfamily", "external_bound_positive": "0", "external_bound_negative": "0", "external_bound_challenge": "0", "heldout_positive": "0", "family_resolved_negative_status": "missing", "challenge_status": "missing", "calibration_decision": "planned_not_run_superfamily_requires_family_resolution", "new_family_call_created": "false"},
                ])
            result = module.finalize(readiness, root / "out", run_id="test")
            self.assertEqual(result["summary"]["calibrated_candidate_model"], 0)
            self.assertEqual(result["summary"]["new_family_call_created"], False)
            rows = list(csv.DictReader((root / "out" / "leaveout_calibration_decisions.tsv").open(encoding="utf-8"), delimiter="\t"))
            self.assertEqual({row["decision"] for row in rows}, {"reference_only_insufficient_panel", "planned_not_run_superfamily_requires_family_resolution"})


class QualifiedEvidenceFamilyGateTests(unittest.TestCase):
    """The functional gate keeps its thresholds but counts qualified independent evidence."""

    def test_family_gate_requires_qualified_e2_e3_positives_and_independent_genera(self):
        with tempfile.TemporaryDirectory() as temporary:
            _module, result, _out = run_finalize([passing_row()], Path(temporary))
            decision = result["decisions"][0]
            self.assertEqual(decision["decision"], NOT_PROMOTED)
            self.assertEqual(decision["calibration_status"], NOT_PROMOTED)
            self.assertEqual(result["summary"][NOT_PROMOTED], 1)
            self.assertEqual(result["summary"][PROMOTED], 0)

    def test_two_qualified_positives_do_not_pass_the_family_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            _module, result, _out = run_finalize(
                [passing_row(qualified_e2_e3_positive_count="2")], Path(temporary)
            )
            self.assertEqual(result["decisions"][0]["decision"], SUPPLIED_DECISION)
            self.assertEqual(result["decisions"][0]["calibration_status"], "not_run")

    def test_two_independent_genera_do_not_pass_the_family_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            _module, result, _out = run_finalize(
                [passing_row(independent_genus_count="2")], Path(temporary)
            )
            self.assertEqual(result["decisions"][0]["decision"], SUPPLIED_DECISION)
            self.assertEqual(result["decisions"][0]["calibration_status"], "not_run")

    def test_legacy_positive_count_no_longer_satisfies_the_family_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            _module, result, out = run_finalize(
                [
                    calibration_row(
                        existing_experimental_positive="3",
                        external_bound_positive="1",
                        heldout_positive="1",
                        external_bound_negative="1",
                        external_bound_challenge="1",
                        family_resolved_negative_status="sufficient",
                        challenge_status="sufficient",
                    )
                ],
                Path(temporary),
            )
            decision = result["decisions"][0]
            self.assertEqual(decision["decision"], SUPPLIED_DECISION)
            self.assertEqual(decision["positive_count_for_gate"], "4")
            detail = read_tsv(out / "leaveout_calibration_gate_detail.tsv")[0]
            self.assertEqual(detail["gate_passed"], "false")
            self.assertEqual(detail["qualified_e2_e3_positive_count"], "0")
            self.assertEqual(detail["independent_genus_count"], "0")
            self.assertIn("qualified_e2_e3_positive_below_3", detail["blocking_slots"])
            self.assertIn("independent_genus_below_3", detail["blocking_slots"])

    def test_gate_detail_reports_the_qualified_slots_of_a_passing_profile(self):
        with tempfile.TemporaryDirectory() as temporary:
            _module, _result, out = run_finalize([passing_row()], Path(temporary))
            detail = read_tsv(out / "leaveout_calibration_gate_detail.tsv")[0]
            self.assertEqual(detail["gate_kind"], "family")
            self.assertEqual(detail["gate_resolution"], "family")
            self.assertEqual(detail["gate_passed"], "true")
            self.assertEqual(detail["gate_evaluated"], "true")
            self.assertEqual(detail["qualified_e2_e3_positive_count"], "3")
            self.assertEqual(detail["independent_genus_count"], "3")
            self.assertEqual(detail["required_positive_count"], "3")
            self.assertEqual(detail["required_independent_genus_count"], "3")
            self.assertEqual(detail["promotion_authorized"], "false")

    def test_negative_qualified_counts_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                run_finalize([passing_row(qualified_e2_e3_positive_count="-1")], Path(temporary))

    def test_promotion_is_reported_even_when_no_gate_passes(self):
        with tempfile.TemporaryDirectory() as temporary:
            _module, result, out = run_finalize([calibration_row()], Path(temporary), promote=True)
            detail = read_tsv(out / "leaveout_calibration_gate_detail.tsv")[0]
            self.assertEqual(detail["promotion_authorized"], "true")
            self.assertEqual(detail["gate_passed"], "false")
            self.assertEqual(result["decisions"][0]["decision"], SUPPLIED_DECISION)


class CalibratedModelPromotionTests(unittest.TestCase):
    """calibrated_candidate_model is emitted only by an explicitly authorized action."""

    def test_calibrated_candidate_model_requires_the_promotion_flag(self):
        with tempfile.TemporaryDirectory() as temporary:
            _module, unpromoted, _out = run_finalize([passing_row()], Path(temporary))
            self.assertEqual(unpromoted["summary"][PROMOTED], 0)
            self.assertEqual(unpromoted["summary"][NOT_PROMOTED], 1)
        with tempfile.TemporaryDirectory() as temporary:
            _module, promoted, out = run_finalize([passing_row()], Path(temporary), promote=True)
            decision = promoted["decisions"][0]
            self.assertEqual(decision["decision"], PROMOTED)
            self.assertEqual(decision["calibration_status"], "passed")
            self.assertEqual(promoted["summary"][PROMOTED], 1)
            self.assertEqual(promoted["summary"][NOT_PROMOTED], 0)
            self.assertTrue(promoted["report"]["gate"]["promotion"]["authorized"])
            self.assertEqual(read_tsv(out / "leaveout_calibration_decisions.tsv")[0]["decision"], PROMOTED)

    def test_promotion_flag_cannot_bypass_a_failing_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                run_finalize(
                    [passing_row(qualified_e2_e3_positive_count="2", calibration_decision=PROMOTED)],
                    Path(temporary),
                    promote=True,
                )

    def test_calibrated_decision_without_a_gate_still_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                run_finalize([calibration_row(calibration_decision=PROMOTED)], Path(temporary))

    def test_cli_defaults_to_not_promoting(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            readiness = write_readiness(root / "readiness.tsv", [passing_row()])
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = module.main([
                    "--readiness", str(readiness),
                    "--output-dir", str(root / "cli-default"),
                    "--run-id", "cli-default",
                ])
            self.assertEqual(exit_code, 0)
            summary = json.loads(stdout.getvalue())
            self.assertEqual(summary[PROMOTED], 0)
            self.assertEqual(summary[NOT_PROMOTED], 1)

    def test_cli_promotes_only_with_the_explicit_flag(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            readiness = write_readiness(root / "readiness.tsv", [passing_row()])
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = module.main([
                    "--readiness", str(readiness),
                    "--output-dir", str(root / "cli-promote"),
                    "--run-id", "cli-promote",
                    "--promote-calibrated-model",
                ])
            self.assertEqual(exit_code, 0)
            summary = json.loads(stdout.getvalue())
            self.assertEqual(summary[PROMOTED], 1)
            self.assertEqual(summary[NOT_PROMOTED], 0)


if __name__ == "__main__":
    unittest.main()
