"""The promotion run must be mechanical, scoped and auditable.

``promote_phaded_calibrated_candidate.py`` consumes the unmodified
``finalize_phaded_reference_panel_calibration.py`` for a run whose authorization is
*scoped*:

* every profile in the ``authorized_promotion_frame`` scope is handed to ``finalize`` with
  ``promote_calibrated_model=True`` -- exactly the profiles whose gate returns true become
  ``calibrated_candidate_model``;
* a profile marked ``deferred_design_spec_special_family`` is handed to ``finalize`` in its
  own invocation with the flag **off**, which is the module's own
  ``candidate_gate_passed_not_promoted`` path -- the deferral is therefore visible in the
  gate detail (``gate_passed=true``) instead of being hidden or silently promoted;
* a third, full-frame invocation with the flag off is the control that proves no profile can
  reach ``calibrated_candidate_model`` without the explicit flag;
* the governing ``leaveout_calibration_decisions.tsv`` / ``leaveout_calibration_gate_detail.tsv``
  are the row-for-row union of the scope invocations, in readiness order, and every merged
  row must equal its source row exactly.

These tests pin all of that before the runner exists.
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "pipeline" / "scripts" / "promote_phaded_calibrated_candidate.py"

READINESS_FIELDS = (
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "existing_experimental_positive", "external_bound_positive", "external_bound_negative",
    "external_bound_challenge", "heldout_positive", "family_resolved_negative_status",
    "challenge_status", "calibration_decision", "new_family_call_created",
    "qualified_e2_e3_positive_count", "independent_genus_count", "qualification_source",
    "readiness_row_origin", "readiness_notes", "promotion_scope",
)

AUTHORIZED_SCOPE = "authorized_promotion_frame"
DEFERRED_SCOPE = "deferred_design_spec_special_family"


def load_module():
    spec = importlib.util.spec_from_file_location("promote_phaded_calibrated_candidate", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(READINESS_FIELDS), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in READINESS_FIELDS})
    return path


def readiness_row(profile_id, scope, *, kind="family", qualified="4", genera="4",
                  heldout="1", negative="1", challenge="2", negative_status="sufficient",
                  challenge_status="sufficient", decision="candidate_gate_passed_not_promoted"):
    return {
        "profile_id": profile_id, "profile_kind": kind, "phaded_superfamily": "sf",
        "phaded_family_id": "DED_hfam_900" if kind == "family" else "",
        "existing_experimental_positive": "0", "external_bound_positive": "0",
        "external_bound_negative": negative, "external_bound_challenge": challenge,
        "heldout_positive": heldout, "family_resolved_negative_status": negative_status,
        "challenge_status": challenge_status, "calibration_decision": decision,
        "new_family_call_created": "false",
        "qualified_e2_e3_positive_count": qualified, "independent_genus_count": genera,
        "qualification_source": "ledger.tsv::independence_group",
        "readiness_row_origin": "frozen_readiness_verbatim", "readiness_notes": "",
        "promotion_scope": scope,
    }


class Fixture:
    def __init__(self, root: Path):
        self.root = root
        self.rows = [
            readiness_row("family_PASS_aaaa", AUTHORIZED_SCOPE),
            readiness_row("family_DEFER_bbbb", DEFERRED_SCOPE),
            readiness_row("family_FAIL_cccc", AUTHORIZED_SCOPE, qualified="1", genera="1",
                          decision="reference_only_insufficient_panel"),
            readiness_row("superfamily_SF_dddd", AUTHORIZED_SCOPE, kind="superfamily",
                          qualified="0", genera="0", heldout="0", negative="0", challenge="0",
                          negative_status="missing", challenge_status="missing",
                          decision="planned_not_run_superfamily_requires_family_resolution"),
            readiness_row("family_FROZEN_eeee", AUTHORIZED_SCOPE, qualified="0", genera="0",
                          heldout="0", negative="0", challenge="0", negative_status="missing",
                          challenge_status="missing",
                          decision="reference_only_insufficient_panel"),
        ]
        self.readiness = write_tsv(root / "calibration_readiness_v2.tsv", self.rows)
        self.out = root / "results"

    def promote(self, module, **kwargs):
        options = dict(readiness=self.readiness, output_dir=self.out, run_id="test",
                       defer_profiles=("family_DEFER_bbbb",))
        options.update(kwargs)
        output_dir = options.pop("output_dir")
        return module.promote(output_dir=output_dir, **options)


class PromotionTests(unittest.TestCase):
    def test_only_the_authorized_passing_profile_is_promoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=True)
            fields, decisions = read_tsv(fixture.out / "leaveout_calibration_decisions.tsv")
            promoted = [row["profile_id"] for row in decisions
                        if row["decision"] == "calibrated_candidate_model"]
            self.assertEqual(promoted, ["family_PASS_aaaa"])

    def test_deferred_profile_is_reported_as_gate_passed_not_promoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=True)
            _, decisions = read_tsv(fixture.out / "leaveout_calibration_decisions.tsv")
            deferred = {row["profile_id"]: row for row in decisions}["family_DEFER_bbbb"]
            self.assertEqual(deferred["decision"], "candidate_gate_passed_not_promoted")
            self.assertEqual(deferred["calibration_status"], "candidate_gate_passed_not_promoted")
            _, gate_detail = read_tsv(fixture.out / "leaveout_calibration_gate_detail.tsv")
            detail = {row["profile_id"]: row for row in gate_detail}["family_DEFER_bbbb"]
            self.assertEqual(detail["gate_passed"], "true")
            self.assertEqual(detail["blocking_slots"], "")
            self.assertEqual(detail["promotion_authorized"], "false")

    def test_promotion_authorized_is_true_only_inside_the_authorized_frame(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=True)
            _, gate_detail = read_tsv(fixture.out / "leaveout_calibration_gate_detail.tsv")
            authorized = {row["profile_id"] for row in gate_detail
                          if row["promotion_authorized"] == "true"}
            self.assertEqual(authorized, {"family_PASS_aaaa", "family_FAIL_cccc",
                                          "superfamily_SF_dddd", "family_FROZEN_eeee"})

    def test_merged_files_cover_the_frame_in_readiness_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=True)
            _, decisions = read_tsv(fixture.out / "leaveout_calibration_decisions.tsv")
            self.assertEqual([row["profile_id"] for row in decisions],
                             [row["profile_id"] for row in fixture.rows])
            _, gate_detail = read_tsv(fixture.out / "leaveout_calibration_gate_detail.tsv")
            self.assertEqual([row["profile_id"] for row in gate_detail],
                             [row["profile_id"] for row in fixture.rows])

    def test_merged_rows_equal_their_scope_source_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            result = fixture.promote(module, promote_calibrated_model=True)
            _, merged = read_tsv(fixture.out / "leaveout_calibration_decisions.tsv")
            merged_by_id = {row["profile_id"]: row for row in merged}
            for scope in (AUTHORIZED_SCOPE, DEFERRED_SCOPE):
                path = fixture.out / "scopes" / scope / "leaveout_calibration_decisions.tsv"
                _, scope_rows = read_tsv(path)
                for row in scope_rows:
                    self.assertEqual(merged_by_id[row["profile_id"]], row)
            verification = json.loads(
                (fixture.out / "promotion_verification.json").read_text(encoding="utf-8"))
            self.assertTrue(all(verification["checks"].values()), verification["checks"])
            self.assertEqual(result["promoted"], ["family_PASS_aaaa"])

    def test_flag_off_promotes_nothing_and_is_recorded_as_the_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=False)
            _, decisions = read_tsv(fixture.out / "leaveout_calibration_decisions.tsv")
            self.assertEqual([row["profile_id"] for row in decisions
                              if row["decision"] == "calibrated_candidate_model"], [])
            by_id = {row["profile_id"]: row["decision"] for row in decisions}
            self.assertEqual(by_id["family_PASS_aaaa"], "candidate_gate_passed_not_promoted")
            _, control = read_tsv(fixture.out / "control_flag_off" / "leaveout_calibration_decisions.tsv")
            self.assertEqual([row["profile_id"] for row in control
                              if row["decision"] == "calibrated_candidate_model"], [])
            verification = json.loads(
                (fixture.out / "promotion_verification.json").read_text(encoding="utf-8"))
            self.assertFalse(verification["promotion_authorized"])
            self.assertTrue(verification["checks"]["flag_off_control_has_no_promotion"])

    def test_gate_detail_flags_match_the_module_output_per_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=True)
            _, merged = read_tsv(fixture.out / "leaveout_calibration_gate_detail.tsv")
            by_id = {row["profile_id"]: row for row in merged}
            self.assertEqual(by_id["family_PASS_aaaa"]["promotion_authorized"], "true")
            self.assertEqual(by_id["family_DEFER_bbbb"]["promotion_authorized"], "false")
            self.assertEqual(by_id["family_FAIL_cccc"]["gate_passed"], "false")
            self.assertEqual(by_id["superfamily_SF_dddd"]["gate_kind"], "none")


class GateTableTests(unittest.TestCase):
    def test_gate_table_covers_every_profile_with_its_scope_and_decision(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=True)
            fields, table = read_tsv(fixture.out / "per_family_gate_table.tsv")
            self.assertEqual(len(table), len(fixture.rows))
            by_id = {row["profile_id"]: row for row in table}
            for required in ("profile_kind", "phaded_superfamily", "phaded_family_id",
                             "qualified_e2_e3_positive_count", "independent_genus_count",
                             "required_qualified_e2_e3_positive_count",
                             "required_independent_genus_count", "gate_passed", "blocking_slots",
                             "promotion_scope", "final_decision", "calibration_status", "promoted"):
                self.assertIn(required, fields)
            self.assertEqual(by_id["family_PASS_aaaa"]["final_decision"], "calibrated_candidate_model")
            self.assertEqual(by_id["family_PASS_aaaa"]["promoted"], "true")
            self.assertEqual(by_id["family_DEFER_bbbb"]["final_decision"],
                             "candidate_gate_passed_not_promoted")
            self.assertEqual(by_id["family_DEFER_bbbb"]["promoted"], "false")
            self.assertEqual(by_id["family_DEFER_bbbb"]["gate_passed"], "true")
            self.assertEqual(by_id["family_DEFER_bbbb"]["promotion_scope"], DEFERRED_SCOPE)
            self.assertEqual(by_id["family_FAIL_cccc"]["blocking_slots"],
                             "qualified_e2_e3_positive_below_3;independent_genus_below_3")

    def test_scope_map_names_the_invocation_that_decided_each_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fixture.promote(module, promote_calibrated_model=True)
            _, mapping = read_tsv(fixture.out / "promotion_scope_map.tsv")
            by_id = {row["profile_id"]: row for row in mapping}
            self.assertEqual(len(mapping), len(fixture.rows))
            self.assertEqual(by_id["family_DEFER_bbbb"]["source_invocation"],
                             "deferred_design_spec_special_family")
            self.assertEqual(by_id["family_PASS_aaaa"]["source_invocation"],
                             "authorized_promotion_frame")


class FailClosedTests(unittest.TestCase):
    def test_unacknowledged_deferral_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            with self.assertRaises(module.PromotionError):
                fixture.promote(module, promote_calibrated_model=True, defer_profiles=())

    def test_deferral_naming_a_non_deferred_profile_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            with self.assertRaises(module.PromotionError):
                fixture.promote(module, promote_calibrated_model=True,
                                defer_profiles=("family_DEFER_bbbb", "family_FAIL_cccc"))

    def test_unknown_scope_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            rows = [dict(row) for row in fixture.rows]
            rows[1]["promotion_scope"] = "promote_everything"
            write_tsv(fixture.readiness, rows)
            module = load_module()
            with self.assertRaises(module.PromotionError):
                fixture.promote(module, promote_calibrated_model=True,
                                defer_profiles=("family_DEFER_bbbb",))

    def test_empty_authorized_frame_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            rows = [dict(row) for row in fixture.rows]
            for row in rows:
                row["promotion_scope"] = DEFERRED_SCOPE
            write_tsv(fixture.readiness, rows)
            module = load_module()
            with self.assertRaises(module.PromotionError):
                fixture.promote(module, promote_calibrated_model=True,
                                defer_profiles=tuple(row["profile_id"] for row in rows))

    def test_existing_output_files_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            fixture.out.mkdir(parents=True, exist_ok=True)
            (fixture.out / "leaveout_calibration_decisions.tsv").write_text("x", encoding="utf-8")
            module = load_module()
            with self.assertRaises(FileExistsError):
                fixture.promote(module, promote_calibrated_model=True)

    def test_missing_promotion_scope_column_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = Fixture(Path(tmp))
            module = load_module()
            fields, rows = read_tsv(fixture.readiness)
            write_tsv(fixture.readiness, [{k: v for k, v in row.items() if k != "promotion_scope"}
                                          for row in rows])
            with self.assertRaises(module.PromotionError):
                fixture.promote(module, promote_calibrated_model=True,
                                defer_profiles=("family_DEFER_bbbb",))


if __name__ == "__main__":
    unittest.main()
