"""Fail-closed invariants for the axis-anchored PhaDED discrimination panel.

These tests encode the Task 5 requirements:

* every panel record binds to exactly one discrimination axis and records its
  anchoring basis (accession, evidence source, why the record sits at that end);
* ``annotation_only`` records are never written as formal negatives;
* records that are not bound to an existing family stay ``unresolved``;
* no gate-lowering field or logic exists anywhere in the design stage.

The tests run against small fixtures so they do not depend on the dated run, plus
a few read-only assertions that only execute when the dated run outputs exist.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "pipeline" / "scripts" / "design_phaded_discrimination_panel.py"
FINALIZE = ROOT / "pipeline" / "scripts" / "finalize_phaded_reference_panel_calibration.py"
RUN_DIR = ROOT / "runs" / "20260917_phaded_discrimination_panel_01"
PANEL_OUT = RUN_DIR / "results" / "challenge_panel_by_axis.tsv"
REACHABILITY_OUT = RUN_DIR / "results" / "axis_reachability.json"

AXIS_1 = "axis_1_substrate_chain_length"
AXIS_2 = "axis_2_particle_state"
AXIS_3 = "axis_3_localization"
AXIS_4 = "axis_4_fold_neighborhood"

AXIS_DEFINITIONS_HEADER = [
    "axis_id", "axis_name", "positive_end_state", "negative_end_state",
    "positive_end_qualifying_evidence_kinds", "negative_end_qualifying_evidence_kinds",
    "positive_end_challenge_evidence_kinds", "negative_end_challenge_evidence_kinds",
    "positive_end_anchor_status", "negative_end_anchor_status", "axis_boundary_note",
    "source_doi_or_pmid",
]
ASSIGNMENT_RULES_HEADER = [
    "rule_id", "priority", "source_layer", "axis_id", "axis_end",
    "qualifying_evidence_kind", "condition_field", "condition_operator",
    "condition_value", "rule_description",
]
DECLARATIONS_HEADER = [
    "accession", "include", "axis_id", "axis_end", "anchor_evidence_kind",
    "evidence_status", "family_binding_status_expected", "anchor_justification",
    "source_doi_or_pmid", "evidence_source_path", "exclusion_reason",
]


def load_module():
    spec = importlib.util.spec_from_file_location("design_phaded_discrimination_panel", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    lines.extend("\t".join(str(row.get(field, "")) for field in header) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def panel_row(module, **overrides):
    """Build a minimally valid panel row, then apply overrides."""
    row = {field: "" for field in module.PANEL_FIELDS}
    row.update(
        axis_id=AXIS_1,
        axis_name="substrate chain length",
        axis_end="negative",
        role="challenge",
        accession="TEST0001",
        record_origin="external_acquisition",
        evidence_status="annotation_only",
        anchor_evidence_kind="reported_label",
        anchor_evidence_basis="ledger_substrate_class_MCL",
        family_binding_status="unresolved",
        bound_family="",
        bound_reference_only_profile_id="",
        formal_negative_eligible="false",
        gate_credit_slot="none",
        anchor_justification="synthetic fixture row: label-only MCL evidence, negative end",
        source_doi_or_pmid="10.1186/1471-2105-10-89 (PMID 19296857)",
        evidence_source_path="fixture",
        secondary_axis_capability="",
    )
    row.update(overrides)
    return row


class AxisDefinitionTests(unittest.TestCase):
    def test_declared_axis_ids_are_exactly_the_four_planned_axes(self):
        module = load_module()
        self.assertEqual(
            set(module.AXIS_IDS), {AXIS_1, AXIS_2, AXIS_3, AXIS_4}
        )
        self.assertEqual(len(module.AXIS_IDS), 4)

    def test_axis_definition_requires_both_ends_and_qualifying_evidence_kinds(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "axes.tsv"
            write_tsv(
                path,
                AXIS_DEFINITIONS_HEADER,
                [
                    {
                        "axis_id": AXIS_1,
                        "axis_name": "substrate",
                        "positive_end_state": "SCL",
                        "negative_end_state": "",
                        "positive_end_qualifying_evidence_kinds": "experimental_positive_assay",
                        "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
                        "source_doi_or_pmid": "x",
                    }
                ],
            )
            with self.assertRaises(ValueError):
                module.load_axis_definitions(path)

    def test_axis_definition_table_must_declare_every_known_end(self):
        module = load_module()
        self.assertEqual(set(module.AXIS_ENDS), {"positive", "negative"})
        self.assertEqual(set(module.ROLES), {"positive_anchor", "negative_anchor", "challenge"})
        self.assertEqual(
            set(module.EVIDENCE_STATUSES),
            {"experimental_positive", "experimental_negative", "annotation_only"},
        )
        for axis_id, definition in module.builtin_axis_definitions().items():
            for end in module.AXIS_ENDS:
                self.assertTrue(definition[f"{end}_end_qualifying_evidence_kinds"].strip(), (axis_id, end))

    def test_valid_single_axis_fixture_loads(self):
        module = load_module()
        definition = module.builtin_axis_definitions()[AXIS_1]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "axes.tsv"
            write_tsv(path, AXIS_DEFINITIONS_HEADER, [definition])
            loaded = module.load_axis_definitions(path)
        self.assertEqual(set(loaded), {AXIS_1})
        self.assertEqual(
            loaded[AXIS_1]["negative_end_qualifying_evidence_kinds"], "experimental_negative_assay"
        )

    def test_unknown_evidence_kind_in_axis_definition_is_rejected(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "axes.tsv"
            write_tsv(
                path,
                AXIS_DEFINITIONS_HEADER,
                [
                    {
                        "axis_id": AXIS_1,
                        "axis_name": "substrate",
                        "positive_end_state": "SCL",
                        "negative_end_state": "MCL",
                        "positive_end_qualifying_evidence_kinds": "made_up_kind",
                        "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
                        "source_doi_or_pmid": "x",
                    }
                ],
            )
            with self.assertRaises(ValueError):
                module.load_axis_definitions(path)

    def test_assignment_rules_must_match_implemented_rule_ids(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rules.tsv"
            write_tsv(
                path,
                ASSIGNMENT_RULES_HEADER,
                [
                    {
                        "rule_id": "R4242_invented_rule",
                        "priority": "1",
                        "source_layer": "reference_ledger",
                        "axis_id": AXIS_1,
                        "axis_end": "negative",
                        "qualifying_evidence_kind": "reported_label",
                        "condition_field": "substrate_class",
                        "condition_operator": "equals",
                        "condition_value": "MCL",
                        "rule_description": "invented",
                    }
                ],
            )
            with self.assertRaises(ValueError):
                module.load_assignment_rules(path)


class RoleAndEligibilityTests(unittest.TestCase):
    def test_negative_anchor_requires_formal_experimental_negative(self):
        module = load_module()
        axis = module.builtin_axis_definitions()[AXIS_1]
        self.assertEqual(
            module.classify_role("negative", "experimental_negative_assay", axis),
            "negative_anchor",
        )
        for kind in ("reported_label", "supported_domain_state", "candidate_polar_x1_pattern_state"):
            self.assertEqual(
                module.classify_role("negative", kind, axis), "challenge", kind
            )

    def test_annotation_only_negative_end_record_is_a_challenge(self):
        module = load_module()
        axis = module.builtin_axis_definitions()[AXIS_1]
        role = module.classify_role("negative", "reported_label", axis)
        self.assertEqual(role, "challenge")
        row = panel_row(module, role=role)
        module.validate_panel_invariants([row])

    def test_annotation_only_record_can_never_be_a_formal_negative(self):
        module = load_module()
        row = panel_row(
            module,
            role="negative_anchor",
            evidence_status="annotation_only",
            anchor_evidence_kind="experimental_negative_assay",
            formal_negative_eligible="false",
        )
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([row])

    def test_formal_negative_requires_resolved_binding_to_reference_only_profile(self):
        module = load_module()
        row = panel_row(
            module,
            axis_id=AXIS_1,
            role="negative_anchor",
            evidence_status="experimental_negative",
            anchor_evidence_kind="experimental_negative_assay",
            family_binding_status="unresolved",
            formal_negative_eligible="true",
            gate_credit_slot="family_resolved_negative",
        )
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([row])

    def test_resolved_binding_without_reference_only_profile_is_not_eligible(self):
        module = load_module()
        row = panel_row(
            module,
            axis_id=AXIS_1,
            role="negative_anchor",
            evidence_status="experimental_negative",
            anchor_evidence_kind="experimental_negative_assay",
            family_binding_status="resolved_existing_family",
            bound_family="DED_hfam_4",
            bound_reference_only_profile_id="",
            formal_negative_eligible="false",
            gate_credit_slot="none",
        )
        module.validate_panel_invariants([row])
        self.assertFalse(module.is_formal_negative_eligible(row))

    def test_unresolved_record_must_not_carry_a_bound_family(self):
        module = load_module()
        row = panel_row(
            module,
            family_binding_status="unresolved",
            bound_family="DED_hfam_44",
        )
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([row])

    def test_ambiguous_binding_is_kept_ambiguous_and_never_upgraded(self):
        module = load_module()
        index = module.build_ledger_index(
            [
                {"accession": "AAA111.1", "phaded_family_id": "DED_hfam_44"},
                {"accession": "AAA111.2", "phaded_family_id": "DED_hfam_45"},
                {"accession": "BBB222.1", "phaded_family_id": "DED_hfam_46"},
            ]
        )
        self.assertEqual(module.bind_family("AAA111.1", index)[0], "resolved_existing_family")
        self.assertEqual(module.bind_family("AAA111.3", index)[0], "ambiguous_existing_family")
        self.assertEqual(module.bind_family("BBB222", index)[0], "resolved_existing_family_versionless")
        self.assertEqual(module.bind_family("CCC333.1", index)[0], "unresolved")
        self.assertEqual(module.bind_family("CCC333.1", index)[1], "")


class PanelInvariantTests(unittest.TestCase):
    def test_every_record_binds_exactly_one_axis(self):
        module = load_module()
        rows = [
            panel_row(module, accession="A1", axis_id=AXIS_1),
            panel_row(module, accession="A2", axis_id=AXIS_4, axis_end="positive",
                      anchor_evidence_kind="supported_domain_state", role="positive_anchor"),
        ]
        module.validate_panel_invariants(rows)
        broken = panel_row(module, accession="A3", axis_id="axis_1_substrate_chain_length;axis_2_particle_state")
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([broken])

    def test_duplicate_accession_is_rejected(self):
        module = load_module()
        rows = [panel_row(module, accession="DUP"), panel_row(module, accession="DUP")]
        with self.assertRaises(ValueError):
            module.validate_panel_invariants(rows)

    def test_missing_justification_or_source_is_rejected(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([panel_row(module, accession="A1", anchor_justification="")])
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([panel_row(module, accession="A2", source_doi_or_pmid="")])
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([panel_row(module, accession="A3", anchor_evidence_basis="")])

    def test_unknown_role_or_axis_end_is_rejected(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([panel_row(module, accession="A1", role="formal_negative")])
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([panel_row(module, accession="A2", axis_end="neutral")])

    def test_lid_pattern_state_only_anchors_the_intracellular_end_of_axis_3(self):
        module = load_module()
        good = panel_row(
            module,
            accession="LID1",
            axis_id=AXIS_3,
            axis_end="negative",
            anchor_evidence_kind="supported_pattern_state",
            anchor_evidence_basis="lid_loop_detected_supported",
        )
        module.validate_panel_invariants([good])
        wrong_axis = panel_row(
            module,
            accession="LID2",
            axis_id=AXIS_4,
            axis_end="negative",
            anchor_evidence_kind="supported_pattern_state",
            anchor_evidence_basis="lid_loop_detected_supported",
        )
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([wrong_axis])

    def test_pf06850_is_never_used_as_an_extracellular_sbd_anchor(self):
        module = load_module()
        row = panel_row(
            module,
            accession="P1",
            axis_id=AXIS_3,
            axis_end="positive",
            role="positive_anchor",
            anchor_evidence_kind="experimental_positive_assay",
            evidence_status="experimental_positive",
            anchor_evidence_basis="pf06850_extracellular_sbd",
        )
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([row])

    def test_candidate_universe_rows_cannot_be_lid_anchored(self):
        module = load_module()
        row = panel_row(
            module,
            accession="GCA_000000000.1|c1_1",
            record_origin="candidate_universe",
            axis_id=AXIS_3,
            axis_end="negative",
            anchor_evidence_kind="supported_pattern_state",
            anchor_evidence_basis="lid_loop_detected_supported",
        )
        with self.assertRaises(ValueError):
            module.validate_panel_invariants([row])


class DeclarationTests(unittest.TestCase):
    def test_declaration_requires_reason_and_no_axis_when_excluded(self):
        module = load_module()
        rows = [
            {
                "accession": "X1", "include": "false", "axis_id": "", "axis_end": "",
                "anchor_evidence_kind": "", "evidence_status": "experimental_negative",
                "family_binding_status_expected": "unresolved", "anchor_justification": "",
                "source_doi_or_pmid": "PMID:1", "evidence_source_path": "fixture",
                "exclusion_reason": "",
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "decl.tsv"
            write_tsv(path, DECLARATIONS_HEADER, rows)
            with self.assertRaises(ValueError):
                module.load_external_declarations(path)

    def test_declaration_evidence_kind_must_be_allowed_for_the_axis_end(self):
        module = load_module()
        rows = [
            {
                "accession": "X2", "include": "true", "axis_id": AXIS_1, "axis_end": "negative",
                "anchor_evidence_kind": "supported_domain_state",
                "evidence_status": "annotation_only",
                "family_binding_status_expected": "unresolved",
                "anchor_justification": "a fold domain state cannot anchor the substrate-chain-length axis",
                "source_doi_or_pmid": "PMID:1", "evidence_source_path": "fixture",
                "exclusion_reason": "",
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "decl.tsv"
            write_tsv(path, DECLARATIONS_HEADER, rows)
            with self.assertRaises(ValueError):
                module.load_external_declarations(path)

    def test_fold_domain_state_is_a_challenge_at_the_fold_axis_negative_end(self):
        """α/β-hydrolase domain rows sit at axis 4's negative end, but never as formal negatives."""
        module = load_module()
        row = panel_row(
            module,
            accession="FOLD1",
            axis_id=AXIS_4,
            axis_end="negative",
            anchor_evidence_kind="supported_domain_state",
            anchor_evidence_basis="pfam_alpha_beta_hydrolase_without_PF10503",
        )
        module.validate_panel_invariants([row])
        self.assertEqual(row["role"], "challenge")

    def test_declaration_evidence_status_must_match_frozen_acquisition_row(self):
        module = load_module()
        frozen = {"Y1": {"evidence_status": "experimental_negative"}}
        rows = [
            {
                "accession": "Y1", "include": "true", "axis_id": AXIS_1, "axis_end": "negative",
                "anchor_evidence_kind": "experimental_negative_assay",
                "evidence_status": "annotation_only",
                "family_binding_status_expected": "unresolved",
                "anchor_justification": "upgraded by mistake",
                "source_doi_or_pmid": "PMID:1", "evidence_source_path": "fixture",
                "exclusion_reason": "",
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "decl.tsv"
            write_tsv(path, DECLARATIONS_HEADER, rows)
            with self.assertRaises(ValueError):
                module.load_external_declarations(path, frozen_status=frozen)
            declarations = [module.normalize_declaration(rows[0])]
            with self.assertRaises(ValueError):
                module.validate_external_declarations(declarations, frozen_status=frozen)

    def test_declared_binding_must_match_frozen_acquisition_binding(self):
        module = load_module()
        frozen = {"Z1": {"family_binding_status": "unresolved", "evidence_status": "experimental_negative"}}
        row = {
            "accession": "Z1", "include": "true", "axis_id": AXIS_3, "axis_end": "negative",
            "anchor_evidence_kind": "experimental_negative_assay",
            "evidence_status": "experimental_negative",
            "family_binding_status_expected": "resolved_existing_family",
            "anchor_justification": "invented binding",
            "source_doi_or_pmid": "PMID:1", "evidence_source_path": "fixture",
            "exclusion_reason": "",
        }
        declarations = [module.normalize_declaration(row)]
        with self.assertRaises(ValueError):
            module.validate_external_declarations(declarations, frozen_status=frozen)

    def test_prior_panel_name_maps_to_expected_evidence_status(self):
        module = load_module()
        cases = {
            "independent_experimental_positive": "experimental_positive",
            "mcl_pha_experimental_positive": "experimental_positive",
            "mcl_pha_non_phb_negative": "experimental_negative",
            "intracellular_non_ephaz_negative": "experimental_negative",
            "fragment_or_incomplete_negative": "experimental_negative",
            "annotation_only_near_neighbor_negative": "annotation_only",
        }
        for panel, expected in cases.items():
            self.assertEqual(module.prior_panel_evidence_status(panel), expected, panel)


class ReachabilityTests(unittest.TestCase):
    def axis_rows(self, module):
        return [
            panel_row(module, accession="N1", axis_id=AXIS_1, axis_end="negative",
                      role="challenge", anchor_evidence_kind="reported_label"),
            panel_row(module, accession="N2", axis_id=AXIS_1, axis_end="negative",
                      role="negative_anchor", evidence_status="experimental_negative",
                      anchor_evidence_kind="experimental_negative_assay",
                      formal_negative_eligible="false"),
            panel_row(module, accession="P1", axis_id=AXIS_1, axis_end="positive",
                      role="positive_anchor", evidence_status="experimental_positive",
                      anchor_evidence_kind="experimental_positive_assay"),
        ]

    def test_family_resolved_negative_count_and_gap(self):
        module = load_module()
        rows = self.axis_rows(module)
        report = module.summarize_axis_reachability(rows, module.builtin_axis_definitions())
        axis_1 = report[AXIS_1]
        self.assertEqual(axis_1["anchored_negative_records"], 2)
        self.assertEqual(axis_1["formal_experimental_negatives"], 1)
        self.assertEqual(axis_1["family_resolved_negatives"], 0)
        self.assertEqual(axis_1["gate_credit_negatives"], 0)
        self.assertEqual(axis_1["gap_to_one_family_resolved_negative"], 1)
        self.assertEqual(axis_1["positive_anchors"], 1)

    def test_family_resolved_negative_closes_the_gap_only_when_fully_bound(self):
        module = load_module()
        rows = self.axis_rows(module) + [
            panel_row(
                module,
                accession="N3",
                axis_id=AXIS_1,
                axis_end="negative",
                role="negative_anchor",
                evidence_status="experimental_negative",
                anchor_evidence_kind="experimental_negative_assay",
                family_binding_status="resolved_existing_family",
                bound_family="DED_hfam_44",
                bound_reference_only_profile_id="family_DED_hfam_44_5f6fdcb6699d",
                formal_negative_eligible="true",
                gate_credit_slot="family_resolved_negative",
            )
        ]
        report = module.summarize_axis_reachability(rows, module.builtin_axis_definitions())
        self.assertEqual(report[AXIS_1]["family_resolved_negatives"], 1)
        self.assertEqual(report[AXIS_1]["gate_credit_negatives"], 1)
        self.assertEqual(report[AXIS_1]["gap_to_one_family_resolved_negative"], 0)

    def test_gate_constants_mirror_the_authoritative_finalize_gate(self):
        module = load_module()
        self.assertEqual(module.GATE["minimum_positive_count"], 3)
        self.assertEqual(module.GATE["minimum_heldout_positive_count"], 1)
        self.assertEqual(module.GATE["minimum_family_resolved_negative_count"], 1)
        self.assertEqual(module.GATE["minimum_challenge_count"], 1)
        self.assertTrue(module.GATE["zero_unexplained_negative_or_challenge_hits"])
        source = FINALIZE.read_text(encoding="utf-8")
        for fragment in ("positive_count >= 3", "heldout >= 1", "negative >= 1", "challenge >= 1"):
            self.assertIn(fragment, source)

    def test_profile_gap_reports_every_missing_gate_slot(self):
        module = load_module()
        readiness = [
            {
                "profile_id": "family_DED_hfam_44_5f6fdcb6699d",
                "profile_kind": "family",
                "phaded_superfamily": "extracellular dPHASCL type 1",
                "phaded_family_id": "DED_hfam_44",
                "existing_experimental_positive": "0",
                "external_bound_positive": "0",
                "external_bound_negative": "0",
                "external_bound_challenge": "1",
                "heldout_positive": "0",
                "family_resolved_negative_status": "missing",
                "challenge_status": "sufficient",
                "calibration_decision": "reference_only_insufficient_panel",
            },
            {
                "profile_id": "family_DED_hfam_53_0750b4051ce6",
                "profile_kind": "family",
                "phaded_superfamily": "extracellular dPHASCL type 1",
                "phaded_family_id": "DED_hfam_53",
                "existing_experimental_positive": "2",
                "external_bound_positive": "0",
                "external_bound_negative": "0",
                "external_bound_challenge": "7",
                "heldout_positive": "0",
                "family_resolved_negative_status": "missing",
                "challenge_status": "sufficient",
                "calibration_decision": "reference_only_insufficient_panel",
            },
        ]
        profiles = module.summarize_profile_gaps(readiness, {"family_DED_hfam_44_5f6fdcb6699d": False})
        first = profiles[0]
        self.assertEqual(first["current_positive_count_for_gate"], 0)
        self.assertEqual(first["positive_gap"], 3)
        self.assertEqual(first["heldout_gap"], 1)
        self.assertEqual(first["negative_gap"], 1)
        self.assertEqual(first["challenge_gap"], 0)
        self.assertFalse(first["gate_passed"])
        second = profiles[1]
        self.assertEqual(second["current_positive_count_for_gate"], 2)
        self.assertEqual(second["positive_gap"], 1)
        self.assertIn("positive_gap", second["blocking_reasons"])

    def test_axis_reachability_json_is_gate_preserving(self):
        module = load_module()
        rows = self.axis_rows(module)
        document = module.build_reachability_document(
            rows,
            module.builtin_axis_definitions(),
            [],
            [],
            [],
            run_id="fixture_run",
        )
        self.assertFalse(document["gate"]["gate_relaxation_applied"])
        self.assertFalse(document["gate"]["gate_relaxation_supported"])
        self.assertEqual(document["gate"]["minimum_positive_count"], 3)
        self.assertEqual(document["summary"]["profiles_meeting_full_gate"], 0)
        self.assertNotIn("relax", json.dumps(document).lower().replace("gate_relaxation", ""))


class NoGateLoweringTests(unittest.TestCase):
    def test_panel_field_names_contain_no_gate_lowering_vocabulary(self):
        module = load_module()
        forbidden = ("relax", "lower", "loosen", "discount", "waive", "override_threshold")
        for field in module.PANEL_FIELDS:
            for token in forbidden:
                self.assertNotIn(token, field.lower(), field)

    def test_module_exposes_no_gate_lowering_attribute_or_function(self):
        module = load_module()
        forbidden = ("relax", "lower", "loosen", "waive", "downgrade")
        for name in dir(module):
            if name.startswith("__"):
                continue
            for token in forbidden:
                self.assertNotIn(token, name.lower(), name)

    def test_output_writer_refuses_to_overwrite_existing_results(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "challenge_panel_by_axis.tsv"
            out.write_text("existing\n", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                module.write_panel_tsv(out, [])


@unittest.skipUnless(PANEL_OUT.is_file() and REACHABILITY_OUT.is_file(), "dated run outputs not present")
class DatedRunTests(unittest.TestCase):
    def test_generated_panel_satisfies_all_invariants(self):
        module = load_module()
        rows = module.read_tsv(PANEL_OUT)
        self.assertEqual(len(rows), len({row["accession"] for row in rows}))
        self.assertTrue(rows)
        module.validate_panel_invariants(rows)
        self.assertEqual({row["axis_id"] for row in rows} <= set(module.AXIS_IDS), True)

    def test_generated_panel_has_no_formal_negative_from_annotation_only(self):
        module = load_module()
        rows = module.read_tsv(PANEL_OUT)
        for row in rows:
            if row["evidence_status"] == "annotation_only":
                self.assertEqual(row["formal_negative_eligible"], "false", row["accession"])
                self.assertNotEqual(row["role"], "negative_anchor", row["accession"])

    def test_generated_reachability_reports_all_profiles_blocked(self):
        module = load_module()
        document = json.loads(REACHABILITY_OUT.read_text(encoding="utf-8"))
        self.assertEqual(document["summary"]["reference_only_profiles"], 37)
        self.assertEqual(document["summary"]["profiles_meeting_full_gate"], 0)
        self.assertFalse(document["gate"]["gate_relaxation_applied"])

    def test_q84c08_is_a_single_axis1_negative_anchor_with_unresolved_binding(self):
        module = load_module()
        rows = [row for row in module.read_tsv(PANEL_OUT) if row["accession"] == "Q84C08"]
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["axis_id"], AXIS_1)
        self.assertEqual(row["axis_end"], "negative")
        self.assertEqual(row["role"], "negative_anchor")
        self.assertEqual(row["family_binding_status"], "unresolved")
        self.assertEqual(row["formal_negative_eligible"], "false")


if __name__ == "__main__":
    unittest.main()
