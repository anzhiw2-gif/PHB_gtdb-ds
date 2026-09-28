"""Axis-anchored credit conventions for the PhaDED reference-only gate (rule #3).

The frozen per-profile gate counters only sum records bound to the profile's **own**
family (``reconcile_phaded_reference_panel.py`` L135-147 and
``acquire_phaded_reference_panel_amendment.py`` L184-191).  Every control that the
discrimination panel anchors to a **discrimination axis** instead of to a family
therefore contributes a constant zero, and the 960-record panel is invisible to
calibration.

These tests pin the corrected convention *before* it is implemented:

* the legacy counters (``external_bound_positive`` / ``external_bound_negative`` /
  ``external_bound_challenge`` / ``heldout_positive``) keep their exact meaning and
  are never fed by the new fields;
* the new ``axis_anchored_*`` fields only count a record when its anchor basis is
  auditable (axis + ``anchor_justification`` + ``source_doi_or_pmid``), when its
  ``anchor_evidence_kind`` is in the axis' pre-registered qualifying list, and when
  the credited profile holds a family-resolved anchor on the **opposite** end of the
  same axis.  Anything else is refused, fail-closed.
"""

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "pipeline" / "scripts" / "credit_phaded_axis_anchored_controls.py"
RECONCILE = REPO / "pipeline" / "scripts" / "reconcile_phaded_reference_panel.py"

PANEL = REPO / "runs" / "20260917_phaded_discrimination_panel_01" / "results" / "challenge_panel_by_axis.tsv"
AXIS_REACHABILITY = REPO / "runs" / "20260917_phaded_discrimination_panel_01" / "results" / "axis_reachability.json"
READINESS = REPO / "runs" / "20260915_phaded_reference_panel_acquisition_02" / "results" / "amendment" / "profile_calibration_readiness.tsv"
PROFILE_MANIFEST = REPO / "runs" / "20260917_phaded_classification_authority_01" / "inputs" / "profile_manifest.tsv"

AXIS_1 = "axis_1_substrate_chain_length"
AXIS_3 = "axis_3_localization"
AXIS_4 = "axis_4_fold_neighborhood"

CHANNEL_NEGATIVE = "axis_anchored_negative"
CHANNEL_CHALLENGE = "axis_anchored_challenge"
CHANNEL_NONE = "no_axis_credit"


def load_module():
    spec = importlib.util.spec_from_file_location("credit_phaded_axis_anchored_controls", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_reconcile():
    spec = importlib.util.spec_from_file_location("reconcile_phaded_reference_panel", RECONCILE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, fields, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# A hermetic two-axis definition set: one axis whose negative end qualifies an
# experimental negative assay, and one axis that has no qualifying negative kind.
SYNTHETIC_AXES = {
    "axis_x": {
        "axis_id": "axis_x",
        "positive_end_state": "target_state",
        "negative_end_state": "contrast_state",
        "positive_end_qualifying_evidence_kinds": "experimental_positive_assay",
        "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
        "positive_end_challenge_evidence_kinds": "reported_label",
        "negative_end_challenge_evidence_kinds": "reported_label;supported_domain_state;candidate_polar_x1_pattern_state",
        "negative_end_anchor_status": "anchored_experimental_family_unresolved",
        "source_doi_or_pmid": "10.0000/synthetic",
    },
    "axis_y": {
        "axis_id": "axis_y",
        "positive_end_state": "target_state_y",
        "negative_end_state": "contrast_state_y",
        "positive_end_qualifying_evidence_kinds": "experimental_positive_assay",
        "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
        "positive_end_challenge_evidence_kinds": "reported_label",
        "negative_end_challenge_evidence_kinds": "reported_label",
        "negative_end_anchor_status": "no_formal_negative_record_available",
        "source_doi_or_pmid": "10.0000/synthetic",
    },
}

TARGET_FAMILY = {
    "profile_id": "family_DED_hfam_X_deadbeef",
    "profile_kind": "family",
    "phaded_superfamily": "sf_x",
    "phaded_family_id": "DED_hfam_X",
}


def record(**overrides) -> dict:
    base = {
        "axis_id": "axis_x",
        "axis_end": "negative",
        "role": "negative_anchor",
        "accession": "NEG1",
        "evidence_status": "experimental_negative",
        "anchor_evidence_kind": "experimental_negative_assay",
        "family_binding_status": "unresolved",
        "bound_family": "",
        "bound_reference_only_profile_id": "",
        "anchor_justification": "explicitly reported as not hydrolysing PHB",
        "source_doi_or_pmid": "10.0000/synthetic",
    }
    base.update(overrides)
    return base


def anchor(**overrides) -> dict:
    base = {
        "axis_id": "axis_x",
        "axis_end": "positive",
        "role": "positive_anchor",
        "accession": "POS1",
        "evidence_status": "experimental_positive",
        "anchor_evidence_kind": "experimental_positive_assay",
        "family_binding_status": "resolved_existing_family",
        "bound_family": "DED_hfam_X",
        "bound_reference_only_profile_id": "family_DED_hfam_X_deadbeef",
        "anchor_justification": "ledger experimental positive at the target end",
        "source_doi_or_pmid": "10.0000/synthetic",
    }
    base.update(overrides)
    return base


def credit_for(module, records, profiles=(TARGET_FAMILY,), family_superfamily=None, axes=None):
    return module.attribute_axis_credit(
        records,
        axes if axes is not None else SYNTHETIC_AXES,
        list(profiles),
        family_superfamily if family_superfamily is not None else {"DED_hfam_X": "sf_x"},
    )


class AxisAnchorAuditTests(unittest.TestCase):
    """The new axis fields must refuse credit without an auditable anchor basis."""

    def test_missing_anchor_justification_is_rejected_fail_closed(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.audit_anchor_basis(record(anchor_justification=""))
        with self.assertRaises(ValueError):
            module.audit_anchor_basis(record(anchor_justification="   "))

    def test_missing_source_doi_or_pmid_is_rejected_fail_closed(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.audit_anchor_basis(record(source_doi_or_pmid=""))

    def test_missing_axis_id_is_rejected_fail_closed(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.audit_anchor_basis(record(axis_id=""))

    def test_unknown_axis_is_rejected_fail_closed(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.classify_record(record(axis_id="axis_does_not_exist"), SYNTHETIC_AXES)

    def test_audited_record_returns_the_audited_status(self):
        module = load_module()
        self.assertEqual(module.audit_anchor_basis(record()), "audited")

    def test_loaded_frozen_axis_definitions_cover_every_frozen_panel_axis(self):
        module = load_module()
        axes = module.load_axis_definitions(AXIS_REACHABILITY)
        observed = {row["axis_id"] for row in read_tsv(PANEL)}
        self.assertEqual(set(axes), observed)
        self.assertEqual(axes[AXIS_1]["negative_end_qualifying_evidence_kinds"], "experimental_negative_assay")


class AxisCreditChannelTests(unittest.TestCase):
    """Only pre-registered qualifying evidence kinds may open a credit channel."""

    def test_negative_channel_requires_a_qualifying_experimental_negative_assay(self):
        module = load_module()
        channel, _reason = module.classify_record(record(), SYNTHETIC_AXES)
        self.assertEqual(channel, CHANNEL_NEGATIVE)
        channel, reason = module.classify_record(
            record(anchor_evidence_kind="reported_label"), SYNTHETIC_AXES
        )
        self.assertEqual(channel, CHANNEL_NONE)
        self.assertIn("reported_label", reason)

    def test_polar_x1_pattern_never_counts_as_negative(self):
        module = load_module()
        channel, reason = module.classify_record(
            record(
                axis_id=AXIS_4,
                anchor_evidence_kind="candidate_polar_x1_pattern_state",
                evidence_status="annotation_only",
                role="challenge",
            ),
            module.load_axis_definitions(AXIS_REACHABILITY),
        )
        self.assertNotEqual(channel, CHANNEL_NEGATIVE)
        self.assertEqual(channel, CHANNEL_CHALLENGE)
        # and a record that claims the negative end but carries only a pattern state
        channel, reason = module.classify_record(
            record(
                axis_id=AXIS_4,
                anchor_evidence_kind="candidate_polar_x1_pattern_state",
                evidence_status="experimental_negative",
            ),
            module.load_axis_definitions(AXIS_REACHABILITY),
        )
        self.assertEqual(channel, CHANNEL_NONE)
        self.assertIn("candidate_polar_x1_pattern_state", reason)

    def test_frozen_confounder_records_are_never_negative_credit(self):
        module = load_module()
        result = credit_for(
            module,
            read_tsv(PANEL),
            profiles=read_tsv(READINESS),
            family_superfamily=module.load_family_superfamily(PROFILE_MANIFEST),
            axes=module.load_axis_definitions(AXIS_REACHABILITY),
        )
        offenders = [
            row["accession"]
            for row in result["records"]
            if row["anchor_evidence_kind"] == "candidate_polar_x1_pattern_state"
            and row["credit_channel"] == CHANNEL_NEGATIVE
        ]
        self.assertEqual(offenders, [])


class AxisCreditAttributionTests(unittest.TestCase):
    """Credit needs a family-resolved anchor on the opposite end of the same axis."""

    def test_axis_credit_requires_an_opposing_end_resolved_anchor(self):
        module = load_module()
        without_anchor = credit_for(module, [record()])
        self.assertEqual(without_anchor["profiles"][TARGET_FAMILY["profile_id"]]["axis_anchored_negative_count"], 0)
        self.assertEqual(
            without_anchor["records"][0]["credit_reason"], "no_opposing_end_resolved_anchor"
        )

        with_anchor = credit_for(module, [record(), anchor()])
        self.assertEqual(with_anchor["profiles"][TARGET_FAMILY["profile_id"]]["axis_anchored_negative_count"], 1)
        self.assertEqual(with_anchor["records"][0]["credited_profile_ids"], [TARGET_FAMILY["profile_id"]])

    def test_same_end_anchor_does_not_create_discrimination_credit(self):
        module = load_module()
        result = credit_for(
            module,
            [record(), anchor(axis_end="negative", role="challenge", anchor_evidence_kind="reported_label")],
        )
        self.assertEqual(result["profiles"][TARGET_FAMILY["profile_id"]]["axis_anchored_negative_count"], 0)

    def test_unresolved_anchors_do_not_credit_anything(self):
        module = load_module()
        result = credit_for(
            module,
            [record(), anchor(family_binding_status="unresolved", bound_family="", bound_reference_only_profile_id="")],
        )
        self.assertEqual(result["profiles"][TARGET_FAMILY["profile_id"]]["axis_anchored_negative_count"], 0)

    def test_axis_anchor_basis_status_is_not_audited_without_credit_records(self):
        module = load_module()
        result = credit_for(module, [anchor()])
        self.assertEqual(result["profiles"][TARGET_FAMILY["profile_id"]]["axis_anchor_basis_status"], "not_audited")
        self.assertEqual(result["profiles"][TARGET_FAMILY["profile_id"]]["axis_anchored_negative_count"], 0)

    def test_frozen_panel_credits_q84c08_only(self):
        module = load_module()
        result = credit_for(
            module,
            read_tsv(PANEL),
            profiles=read_tsv(READINESS),
            family_superfamily=module.load_family_superfamily(PROFILE_MANIFEST),
            axes=module.load_axis_definitions(AXIS_REACHABILITY),
        )
        by_accession = {row["accession"]: row for row in result["records"]}
        for accession in ("Q84C08", "O87189", "Q7WT48", "Q7WT49"):
            self.assertIn(accession, by_accession)
            self.assertEqual(by_accession[accession]["family_resolved_negative_slot_credited"], False)
        self.assertTrue(by_accession["Q84C08"]["axis_anchored_slot_credited"])
        for accession in ("O87189", "Q7WT48", "Q7WT49"):
            self.assertFalse(by_accession[accession]["axis_anchored_slot_credited"])
            self.assertEqual(by_accession[accession]["credit_reason"], "no_opposing_end_resolved_anchor")

    def test_frozen_panel_superfamily_credit_stays_below_every_threshold(self):
        module = load_module()
        result = credit_for(
            module,
            read_tsv(PANEL),
            profiles=read_tsv(READINESS),
            family_superfamily=module.load_family_superfamily(PROFILE_MANIFEST),
            axes=module.load_axis_definitions(AXIS_REACHABILITY),
        )
        superfamilies = {
            profile_id: payload
            for profile_id, payload in result["profiles"].items()
            if profile_id.startswith("superfamily_")
        }
        self.assertEqual(len(superfamilies), 4)
        for profile_id, payload in superfamilies.items():
            self.assertLessEqual(payload["axis_anchored_negative_count"], 1, profile_id)


class LegacyCounterIsolationTests(unittest.TestCase):
    """The new fields must not feed, alter or dilute any legacy counter."""

    def _reconcile_inputs(self, root: Path, panel_axes: Path | None):
        profiles = root / "profiles.tsv"
        ledger = root / "ledger.tsv"
        panel = root / "panel.tsv"
        _write(
            profiles,
            ["profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id", "model_status"],
            [
                {"profile_id": "family_DED_hfam_X_deadbeef", "profile_kind": "family", "phaded_superfamily": "sf_x", "phaded_family_id": "DED_hfam_X", "model_status": "reference_only"},
            ],
        )
        _write(
            ledger,
            ["accession", "phaded_superfamily", "phaded_family_id", "evidence_status"],
            [],
        )
        _write(
            panel,
            ["accession", "panel", "pmid", "doi", "sequence_length", "normalized_sequence_sha256", "source_sha256", "independence_status", "decision"],
            [
                {"accession": "NEG1", "panel": "mcl_pha_non_phb_negative", "pmid": "PMID:1", "doi": "", "sequence_length": "100", "normalized_sequence_sha256": "a" * 64, "source_sha256": "b" * 64, "independence_status": "independent", "decision": "accept_panel"},
            ],
        )
        return profiles, ledger, panel

    def test_reconcile_without_axis_panel_emits_legacy_columns_only(self):
        module = load_reconcile()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles, ledger, panel = self._reconcile_inputs(root, None)
            module.reconcile(profiles, ledger, panel, root / "out", expected_reference_only_count=1)
            header = (root / "out" / "profile_calibration_readiness.tsv").read_text(encoding="utf-8").splitlines()[0]
            self.assertNotIn("axis_anchored", header)
            self.assertEqual(
                header,
                "\t".join(
                    [
                        "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
                        "existing_experimental_positive", "external_bound_positive",
                        "external_bound_negative", "external_bound_challenge", "heldout_positive",
                        "family_resolved_negative_status", "challenge_status", "calibration_decision",
                        "new_family_call_created",
                    ]
                ),
            )

    def test_legacy_counters_are_not_changed_by_axis_credit(self):
        module = load_reconcile()
        credit = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles, ledger, panel = self._reconcile_inputs(root, None)
            baseline = module.reconcile(profiles, ledger, panel, root / "base", expected_reference_only_count=1)
            baseline_row = baseline["readiness"][0]

            axis_panel = root / "axis_panel.tsv"
            _write(
                axis_panel,
                sorted(
                    {
                        "axis_id", "axis_end", "role", "accession", "evidence_status",
                        "anchor_evidence_kind", "family_binding_status", "bound_family",
                        "bound_reference_only_profile_id", "anchor_justification", "source_doi_or_pmid",
                        "axis_name", "record_origin", "anchor_evidence_basis", "formal_negative_eligible",
                        "gate_credit_slot", "evidence_source_path", "secondary_axis_capability",
                    }
                ),
                [
                    {"axis_id": "axis_x", "axis_end": "negative", "role": "negative_anchor", "accession": "NEG1", "evidence_status": "experimental_negative", "anchor_evidence_kind": "experimental_negative_assay", "family_binding_status": "unresolved", "bound_family": "", "bound_reference_only_profile_id": "", "anchor_justification": "explicitly not hydrolysing PHB", "source_doi_or_pmid": "10.0000/synthetic", "axis_name": "x", "record_origin": "external_acquisition", "anchor_evidence_basis": "b", "formal_negative_eligible": "false", "gate_credit_slot": "none", "evidence_source_path": "p", "secondary_axis_capability": ""},
                    {"axis_id": "axis_x", "axis_end": "positive", "role": "positive_anchor", "accession": "POS1", "evidence_status": "experimental_positive", "anchor_evidence_kind": "experimental_positive_assay", "family_binding_status": "resolved_existing_family", "bound_family": "DED_hfam_X", "bound_reference_only_profile_id": "family_DED_hfam_X_deadbeef", "anchor_justification": "target end anchor", "source_doi_or_pmid": "10.0000/synthetic", "axis_name": "x", "record_origin": "reference_ledger", "anchor_evidence_basis": "b", "formal_negative_eligible": "false", "gate_credit_slot": "none", "evidence_source_path": "p", "secondary_axis_capability": ""},
                ],
            )
            axes = root / "axes.json"
            axes.write_text(json.dumps({"axis_definitions": [SYNTHETIC_AXES["axis_x"]]}), encoding="utf-8")
            enriched = module.reconcile(
                profiles, ledger, panel, root / "enriched",
                expected_reference_only_count=1,
                axis_panel=axis_panel,
                axis_definitions=axes,
            )
            enriched_row = enriched["readiness"][0]
            for field in (
                "existing_experimental_positive", "external_bound_positive",
                "external_bound_negative", "external_bound_challenge", "heldout_positive",
                "family_resolved_negative_status", "challenge_status", "calibration_decision",
            ):
                self.assertEqual(enriched_row[field], baseline_row[field], field)
            self.assertEqual(enriched_row["external_bound_negative"], 0)
            self.assertEqual(enriched_row["axis_anchored_negative_count"], 1)
            self.assertEqual(enriched_row["axis_anchor_basis_status"], "audited")

    def test_enriched_reconcile_is_refused_without_an_audited_axis_panel(self):
        module = load_reconcile()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles, ledger, panel = self._reconcile_inputs(root, None)
            with self.assertRaises(ValueError):
                module.reconcile(
                    profiles, ledger, panel, root / "out",
                    expected_reference_only_count=1,
                    axis_panel=root / "missing_panel.tsv",
                )


class SuperfamilyRollupTests(unittest.TestCase):
    """The superfamily-resolution roll-up sums member families and never invents evidence."""

    def test_superfamily_rollup_sums_member_families(self):
        module = load_module()
        readiness = [
            {"profile_id": "family_DED_hfam_A", "profile_kind": "family", "phaded_superfamily": "sf_x", "phaded_family_id": "DED_hfam_A", "existing_experimental_positive": "2", "external_bound_positive": "1", "heldout_positive": "0"},
            {"profile_id": "family_DED_hfam_B", "profile_kind": "family", "phaded_superfamily": "sf_x", "phaded_family_id": "DED_hfam_B", "existing_experimental_positive": "0", "external_bound_positive": "0", "heldout_positive": "0"},
            {"profile_id": "superfamily_sf_x", "profile_kind": "superfamily", "phaded_superfamily": "sf_x", "phaded_family_id": "", "existing_experimental_positive": "0", "external_bound_positive": "0", "heldout_positive": "0"},
        ]
        rollup = module.rollup_superfamily_evidence(
            readiness, {"DED_hfam_A": "sf_x", "DED_hfam_B": "sf_x"}
        )
        payload = rollup["superfamily_sf_x"]
        self.assertEqual(payload["member_family_ids"], ["DED_hfam_A", "DED_hfam_B"])
        self.assertEqual(payload["superfamily_positive_count"], 3)
        self.assertEqual(payload["superfamily_heldout_positive"], 0)

    def test_rollup_reports_member_families_without_a_readiness_row(self):
        module = load_module()
        readiness = [
            {"profile_id": "superfamily_sf_y", "profile_kind": "superfamily", "phaded_superfamily": "sf_y", "phaded_family_id": "", "existing_experimental_positive": "0", "external_bound_positive": "0", "heldout_positive": "0"},
        ]
        rollup = module.rollup_superfamily_evidence(
            readiness, {"DED_hfam_Z": "sf_y"}
        )
        self.assertEqual(rollup["superfamily_sf_y"]["superfamily_positive_count"], 0)
        self.assertEqual(rollup["superfamily_sf_y"]["member_families_without_readiness_row"], ["DED_hfam_Z"])


class FrozenPanelIntegrityTests(unittest.TestCase):
    """This stage is read-only; the frozen panel must stay byte-identical."""

    def test_frozen_panel_and_readiness_are_unmodified(self):
        self.assertEqual(
            hashlib.sha256(PANEL.read_bytes()).hexdigest(),
            "0431a9edf17685b353ac835e45469bed9f5f3bba995c47fd895332b7237d0755",
        )
        self.assertEqual(
            hashlib.sha256(READINESS.read_bytes()).hexdigest(),
            "5a88a4c71e53f4b72e8473c4e5e79b6ff348f66f728095abf99301c58038c7d6",
        )


if __name__ == "__main__":
    unittest.main()
