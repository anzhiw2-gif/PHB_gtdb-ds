"""Tests for the v2 PhaDED candidate catalog builder (evidence-separated semantics).

Written before the implementation (test-first).  The v1 "high-confidence" filter
(``filter_phaded_high_confidence.filter_candidates``) collapses nucleophile type,
transport prediction, localization and architecture into one binary label.  The
v2 catalog keeps them separate and gives every protein exactly one
``primary_disposition``.

Load-bearing invariants pinned here:

* AHSMG is a *Ser motif class*, never a third nucleophile identity
  (``nucleophile_identity`` x ``motif_class`` are two columns);
* SBD is an accessory-domain annotation and never participates in the type 1 /
  type 2 catalytic-domain definition;
* a SignalP class is a transport *prediction* and is never promoted to
  experimental localization evidence;
* ``excluded_input_quality`` is reachable only through a documented
  input-quality reason: scientific uncertainty never routes there, in either
  direction, including a row that carries both kinds of cue;
* localization discordance is an annotation and a functional-view gate, never a
  silent deletion or demotion of a sequence homolog;
* the functional view is derived only and can never change the master catalog
  disposition;
* the explicit 52-record Cys demotion input accounts for every input identifier
  exactly once, and a count that contradicts the declared expectation raises
  instead of proceeding silently.

All fixtures are synthetic and in-memory / tempfile: no project data is read or
recomputed (plan Authorization Checkpoint 5 is not authorized).
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import build_phaded_candidate_catalog_v2 as module  # noqa: E402
import phaded_evidence_schema as schema  # noqa: E402


# --- synthetic fixture vocabulary (mirrors, never reads, project values) ------

PHAZ7 = "extracellular native-SCL/PhaZ7-like"
CYS_NO_LIPASE_BOX = "intracellular nPHASCL without lipase box"
NPHASCL_WITH_LIPASE_BOX = "intracellular nPHASCL with lipase box"
NPHAMCL = "intracellular nPHAMCL"
D_TYPE1 = "extracellular dPHASCL type 1"

DEMOTION_COLUMNS = [
    "accession", "genome", "superfamily", "signalp_class", "profile_best_evalue",
    "profile_evidence_status", "lipase_box_state", "lipase_box_x1",
    "catalytic_domain_type", "catalytic_residue_verification", "lid_state",
    "sbd_pf06850_binding_state", "ahsmg_state", "interpro_status",
    "high_confidence", "hold_reason", "signalp_class_prediction", "demote_basis",
]


def build_row(**overrides) -> dict[str, str]:
    """Build one raw input row with the documented v2 evidence columns.

    Defaults describe an unremarkable validated-lipase-box candidate; each test
    overrides only the field under test.
    """
    row = {
        "accession": "ACC_0001",
        "genome": "GCA_000000001.1",
        "superfamily": D_TYPE1,
        "model_layer": "sequence_family_hmm_validated",
        "functional_calibration_status": "not_function_calibrated",
        "sequence_family_call": "hfam_test",
        "assignment_unique": "true",
        "lipase_box_state": "supported",
        "ahsmg_state": "not_detected_pattern",
        "motif_state": "",
        "sbd_state": "not_detected",
        "lid_state": "not_detected_pattern",
        "linker_state": "not_detected_pattern",
        "with_lipase_state": "",
        "catalytic_domain_type": "type1_verified",
        "signalp_class": "SP",
        "localization_evidence": "unknown",
        "experimental_evidence_grade": "",
    }
    row.update(overrides)
    return row


def demotion_row(index: int, **overrides) -> dict[str, str]:
    """Build one synthetic row in the frozen demotion-input shape."""
    row = {
        "accession": "DEMO_%04d" % index,
        "genome": "GCA_0000000%02d.1" % (index % 100),
        "superfamily": CYS_NO_LIPASE_BOX,
        "signalp_class": "pending",
        "profile_best_evalue": "1e-100",
        "profile_evidence_status": "discovery_hit",
        "lipase_box_state": "not_detected_pattern",
        "lipase_box_x1": "",
        "catalytic_domain_type": "undetermined_no_oxyanion",
        "catalytic_residue_verification": "inferred_cys_not_verified",
        "lid_state": "not_detected_pattern",
        "sbd_pf06850_binding_state": "detected",
        "ahsmg_state": "not_detected_pattern",
        "interpro_status": "interpro_supported",
        "high_confidence": "1",
        "hold_reason": "localization_conflict_signalp_predicted_export",
        "signalp_class_prediction": "SP",
        "demote_basis": "SignalP 6.0h prediction (not experimental); synthetic fixture",
    }
    row.update(overrides)
    return row


def write_tsv(path: Path, columns, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t",
            lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    return path


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


class CatalogTestCase(unittest.TestCase):
    """Base fixture: a temp working directory plus catalog helpers."""

    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_catalog_v2_"))
        self.candidates = self.tmp / "candidates.tsv"
        self.demotion = self.tmp / "catalog_universe_demotion.tsv"
        self.out_dir = self.tmp / "out"

    def build_catalog(self, rows, demotion_rows=(), **kwargs):
        write_tsv(self.candidates, module.CANDIDATE_COLUMNS, rows)
        write_tsv(self.demotion, DEMOTION_COLUMNS, demotion_rows)
        kwargs.setdefault("expect_demoted", len(list(demotion_rows)))
        return module.build_catalog(
            candidates_path=self.candidates,
            demotion_path=self.demotion,
            output_dir=self.out_dir,
            **kwargs,
        )


# ---------------------------------------------------------------------------
# Step 1: corrected field semantics
# ---------------------------------------------------------------------------

class FieldSemanticsTests(unittest.TestCase):
    """Plan Step 1: the three semantic corrections, pinned verbatim."""

    def test_ahsmg_candidate_is_ser_with_ahsmg_motif(self):
        row = build_row(superfamily=PHAZ7, ahsmg_state="supported")
        out = module.classify_row(row)
        self.assertEqual(out["nucleophile_identity"], "ser")
        self.assertEqual(out["motif_class"], "AHSMG")

    def test_sbd_does_not_define_type_one_or_type_two(self):
        row = build_row(catalytic_domain_type="type2_verified", sbd_state="supported")
        out = module.classify_row(row)
        self.assertEqual(out["catalytic_domain_type"], "type2_verified")
        self.assertEqual(out["accessory_domain_architecture"], "SBD:supported")

    def test_signalp_other_is_not_equivalent_to_experimental_intracellular(self):
        out = module.classify_row(
            build_row(signalp_class="OTHER", localization_evidence="unknown")
        )
        self.assertEqual(out["transport_signal_prediction"], "OTHER")
        self.assertEqual(out["localization_evidence"], "unknown")

    def test_ahsmg_is_never_a_nucleophile_identity(self):
        """AHSMG must not survive as a third nucleophile value anywhere."""
        out = module.classify_row(build_row(superfamily=PHAZ7, ahsmg_state="supported"))
        self.assertIn(out["nucleophile_identity"], schema.NUCLEOPHILE_IDENTITIES)
        self.assertNotEqual(out["nucleophile_identity"], "ahsmg")
        with self.assertRaises(ValueError):
            schema.validate_candidate_evidence(
                {**out, "nucleophile_identity": "ahsmg"}
            )

    def test_lipase_box_ser_candidate_keeps_gxsxg_motif_class(self):
        out = module.classify_row(build_row(superfamily=D_TYPE1, lipase_box_state="supported"))
        self.assertEqual(out["nucleophile_identity"], "ser")
        self.assertEqual(out["motif_class"], "GxSxG")

    def test_cys_superfamily_is_cys_with_cys_associated_motif(self):
        out = module.classify_row(
            build_row(
                superfamily=CYS_NO_LIPASE_BOX,
                lipase_box_state="not_detected_pattern",
                catalytic_domain_type="unresolved",
            )
        )
        self.assertEqual(out["nucleophile_identity"], "cys")
        self.assertEqual(out["motif_class"], "Cys-associated")

    def test_accessory_architecture_records_sbd_and_lid_independently(self):
        out = module.classify_row(
            build_row(sbd_state="supported", lid_state="supported")
        )
        self.assertEqual(
            out["accessory_domain_architecture"], "SBD:supported|lid:supported"
        )

    def test_unknown_evidence_columns_emit_documented_defaults_never_invented(self):
        """Missing values keep a documented default; nothing is inferred.

        An absent text column is the empty string (so "no call was made" is
        never confused with a family named ``pending``), while columns whose
        vocabulary has a weakest member use it (``not_tested``, ``unknown``,
        ``unresolved``).
        """
        out = module.classify_row({"accession": "ACC_X"})
        self.assertEqual(out["sequence_family_call"], "")
        self.assertEqual(out["sequence_family_confidence"], "")
        self.assertEqual(out["experimental_evidence_grade"], "")
        self.assertEqual(out["input_quality_reason"], "")
        self.assertEqual(out["catalytic_domain_type"], "unresolved")
        self.assertEqual(out["nucleophile_identity"], "unresolved")
        self.assertEqual(out["motif_class"], "unresolved")
        self.assertEqual(out["localization_evidence"], "unknown")
        self.assertEqual(out["transport_signal_prediction"], "not_tested")
        self.assertEqual(out["accessory_domain_architecture"], "")
        schema.validate_candidate_evidence(out)

    def test_discovery_layer_rows_carry_no_family_call_at_all(self):
        """The discovery layer never makes a family call (AGENTS.md boundary)."""
        out = module.classify_row(
            build_row(model_layer="discovery_hmm_uncalibrated",
                      sequence_family_call="hfam_must_be_dropped")
        )
        self.assertEqual(out["sequence_family_call"], "")
        schema.validate_candidate_evidence(out)
        # A validated layer does carry the call it was given.
        self.assertEqual(
            module.classify_row(build_row())["sequence_family_call"], "hfam_test"
        )

    def test_discovery_layer_row_never_claims_a_core_homolog(self):
        out = module.classify_row(
            build_row(model_layer="discovery_hmm_uncalibrated", assignment_unique="true")
        )
        self.assertEqual(out["primary_disposition"], "probable_sequence_homolog")
        schema.validate_candidate_evidence(out)


# ---------------------------------------------------------------------------
# Step 2: special-family dispositions and the input-quality guard
# ---------------------------------------------------------------------------

class DispositionTests(unittest.TestCase):
    """Plan Step 2: special families plus input-quality exclusion."""

    def test_with_lipase_box_goes_to_deferred_structure_review(self):
        row = build_row(superfamily=NPHASCL_WITH_LIPASE_BOX)
        self.assertEqual(
            module.primary_disposition(row), "deferred_structure_review"
        )
        self.assertEqual(
            module.classify_row(row)["primary_disposition"],
            "deferred_structure_review",
        )

    def test_nphamcl_goes_to_function_unresolved(self):
        row = build_row(superfamily=NPHAMCL)
        self.assertEqual(module.primary_disposition(row), "function_unresolved")

    def test_special_family_precedes_a_validated_model_layer(self):
        """A validated layer never overrides an unresolved competing function."""
        row = build_row(
            superfamily=NPHAMCL,
            model_layer="sequence_family_hmm_validated",
            assignment_unique="true",
        )
        self.assertEqual(module.primary_disposition(row), "function_unresolved")

    def test_validated_unique_assignment_is_core_sequence_homolog(self):
        row = build_row(
            model_layer="sequence_family_hmm_validated", assignment_unique="true"
        )
        self.assertEqual(module.primary_disposition(row), "core_sequence_homolog")

    def test_validated_non_unique_assignment_is_probable(self):
        row = build_row(
            model_layer="sequence_family_hmm_validated", assignment_unique="false"
        )
        self.assertEqual(module.primary_disposition(row), "probable_sequence_homolog")

    def test_discovery_layer_is_probable_at_most(self):
        row = build_row(model_layer="discovery_hmm_uncalibrated", assignment_unique="true")
        self.assertEqual(module.primary_disposition(row), "probable_sequence_homolog")

    def test_reference_query_only_is_remote_homolog_candidate(self):
        row = build_row(model_layer="reference_query_only")
        self.assertEqual(module.primary_disposition(row), "remote_homolog_candidate")

    def test_calibrated_layer_is_never_core_without_unique_assignment(self):
        row = build_row(
            model_layer="calibrated_candidate_model",
            functional_calibration_status="calibrated_candidate_model",
            assignment_unique="false",
        )
        self.assertEqual(module.primary_disposition(row), "probable_sequence_homolog")

    def test_every_disposition_is_a_declared_vocabulary_value(self):
        rows = [
            build_row(superfamily=NPHASCL_WITH_LIPASE_BOX),
            build_row(superfamily=NPHAMCL),
            build_row(model_layer="sequence_family_hmm_validated", assignment_unique="true"),
            build_row(model_layer="discovery_hmm_uncalibrated"),
            build_row(model_layer="reference_query_only"),
            build_row(input_quality_reason="sequence_missing"),
        ]
        for row in rows:
            self.assertIn(module.primary_disposition(row), schema.PRIMARY_DISPOSITIONS)


class InputQualityGuardTests(unittest.TestCase):
    """``excluded_input_quality`` is reachable only through a documented reason."""

    SCIENTIFIC_UNCERTAINTY_ROWS = (
        {"superfamily": "intracellular nPHAMCL"},          # competing function
        {"model_layer": "reference_query_only"},           # ambiguous family
        {"assignment_unique": "false"},                    # not uniquely assigned
        {"superfamily": "not a known superfamily"},        # ambiguous family
        {"signalp_class": "not_tested",
         "localization_evidence": "unknown"},              # unresolved localization
        {"experimental_evidence_grade": "A"},              # low confidence
    )

    def test_documented_input_quality_reason_may_exclude(self):
        for reason in ("sequence_missing", "sequence_truncated", "input_error",
                       "empty_sequence", "malformed_record"):
            out = module.classify_row(build_row(input_quality_reason=reason))
            self.assertEqual(out["primary_disposition"], "excluded_input_quality")
            self.assertEqual(out["input_quality_reason"], reason)

    def test_input_quality_reason_aliases_normalize(self):
        out = module.classify_row(build_row(input_quality_reason="truncated"))
        self.assertEqual(out["input_quality_reason"], "sequence_truncated")
        self.assertEqual(out["primary_disposition"], "excluded_input_quality")

    def test_scientific_uncertainty_never_routes_to_excluded_input_quality(self):
        """The direction that matters most: uncertainty is not an input error."""
        for overrides in self.SCIENTIFIC_UNCERTAINTY_ROWS:
            with self.subTest(overrides=overrides):
                out = module.classify_row(build_row(**overrides))
                self.assertNotEqual(
                    out["primary_disposition"], "excluded_input_quality"
                )
                self.assertIn(out["primary_disposition"], schema.PRIMARY_DISPOSITIONS)
                self.assertNotIn(
                    "excluded_input_quality", out["primary_disposition"]
                )

    def test_scientific_uncertainty_wording_in_the_reason_field_is_rejected(self):
        """A scientific excuse written into an input-quality field fails closed."""
        for wording in ("function_unresolved", "low_confidence", "ambiguous_family",
                        "scientific_uncertainty", "remote_homolog"):
            with self.subTest(wording=wording):
                with self.assertRaises(module.EvidenceError) as caught:
                    module.classify_row(build_row(input_quality_reason=wording))
                self.assertIn(
                    module.RULE_SCIENTIFIC_UNCERTAINTY_NOT_INPUT_ERROR,
                    str(caught.exception),
                )

    def test_adversarial_row_with_both_cues_is_excluded_by_documented_precedence(self):
        """Adjudicated precedence, stated and proved.

        A broken input cannot be classified at all, so the input-quality
        verdict wins -- but only because a *documented* input-quality reason is
        present.  The scientific-uncertainty cues are still recorded in
        ``evidence_flags`` and in ``scientific_uncertainty_present``, so the
        exclusion is auditable and reversible; the row is never *added* to the
        exclusion bucket merely for being uncertain.
        """
        row = build_row(
            superfamily=NPHAMCL,                 # competing function (scientific)
            assignment_unique="false",           # scientific
            input_quality_reason="sequence_truncated",  # documented input quality
        )
        out = module.classify_row(row)
        self.assertEqual(out["primary_disposition"], "excluded_input_quality")
        self.assertEqual(out["input_quality_reason"], "sequence_truncated")
        self.assertEqual(out["scientific_uncertainty_present"], "true")
        self.assertIn(module.FLAG_FUNCTION_COMPETING, out["evidence_flags"])

        # Reverse direction: the same scientific cues without the documented
        # input-quality reason must not reach the exclusion disposition.
        without_input_error = dict(row)
        without_input_error.pop("input_quality_reason")
        self.assertNotEqual(
            module.classify_row(without_input_error)["primary_disposition"],
            "excluded_input_quality",
        )

    def test_excluded_row_is_not_a_sequence_homolog(self):
        out = module.classify_row(
            build_row(
                model_layer="sequence_family_hmm_validated",
                assignment_unique="true",
                input_quality_reason="sequence_missing",
            )
        )
        self.assertEqual(out["primary_disposition"], "excluded_input_quality")


# ---------------------------------------------------------------------------
# Step 3: localization is a compatibility annotation
# ---------------------------------------------------------------------------

class LocalizationCompatibilityTests(unittest.TestCase):
    """Plan Step 3: discordance annotates, it never deletes or demotes."""

    def test_signalp_export_on_an_intracellular_family_is_only_a_flag(self):
        row = build_row(
            superfamily=CYS_NO_LIPASE_BOX,
            lipase_box_state="not_detected_pattern",
            signalp_class="SP",
            localization_evidence="unknown",
            model_layer="sequence_family_hmm_validated",
            assignment_unique="true",
        )
        out = module.classify_row(row)
        self.assertIn(module.FLAG_LOCALIZATION_DISCORDANT, out["evidence_flags"])
        self.assertEqual(out["transport_signal_prediction"], "SP")
        self.assertEqual(out["localization_evidence"], "unknown")
        # Annotation only: the sequence homolog is neither deleted nor demoted.
        self.assertEqual(out["primary_disposition"], "core_sequence_homolog")
        # The master verdict is identical with and without the discordance: the
        # flag annotates a resolved homolog, it does not resolve it differently.
        concordant = module.classify_row(
            {**row, "signalp_class": "OTHER"}
        )
        self.assertEqual(concordant["primary_disposition"],
                         out["primary_disposition"])
        self.assertNotIn(module.FLAG_LOCALIZATION_DISCORDANT,
                         concordant["evidence_flags"])

    def test_concordant_prediction_carries_no_discordance_flag(self):
        out = module.classify_row(build_row(superfamily=D_TYPE1, signalp_class="SP"))
        self.assertNotIn(module.FLAG_LOCALIZATION_DISCORDANT, out["evidence_flags"])

    def test_functional_view_computation_does_not_mutate_primary_disposition(self):
        """Proof: computing the derived view cannot change the master catalog."""
        rows = [
            build_row(accession="ACC_1", signalp_class="SP",
                      experimental_evidence_grade="E3",
                      localization_evidence="other_evidence"),
            build_row(accession="ACC_2", superfamily=CYS_NO_LIPASE_BOX,
                      lipase_box_state="not_detected_pattern", signalp_class="SP"),
            build_row(accession="ACC_3", superfamily=NPHAMCL),
            build_row(accession="ACC_4", model_layer="reference_query_only"),
            build_row(accession="ACC_5", input_quality_reason="sequence_missing"),
        ]
        classified = [module.classify_row(row) for row in rows]
        before = [row["primary_disposition"] for row in classified]
        before_rows = [dict(row) for row in classified]

        eligible = [module.functional_view_eligible(row) for row in classified]
        view = module.build_functional_view(classified)

        after = [row["primary_disposition"] for row in classified]
        self.assertEqual(before, after)
        self.assertEqual(before_rows, [dict(row) for row in classified])
        self.assertEqual(len(eligible), len(classified))
        # The view is a strict, accessioned subset of the master catalog.
        master_dispositions = {row["accession"]: row["primary_disposition"]
                               for row in classified}
        for view_row in view:
            self.assertIn(view_row["accession"], master_dispositions)
            self.assertEqual(
                view_row["primary_disposition"],
                master_dispositions[view_row["accession"]],
            )
        self.assertTrue(view, "fixture must contain at least one eligible row")
        self.assertLess(len(view), len(classified))

    def test_functional_view_excludes_localization_discordant_rows(self):
        row = build_row(
            superfamily=CYS_NO_LIPASE_BOX,
            lipase_box_state="not_detected_pattern",
            signalp_class="SP",
            model_layer="sequence_family_hmm_validated",
            assignment_unique="true",
        )
        classified = module.classify_row(row)
        self.assertIn(module.FLAG_LOCALIZATION_DISCORDANT, classified["evidence_flags"])
        self.assertFalse(module.functional_view_eligible(classified))
        self.assertEqual(module.build_functional_view([classified]), [])

    def test_functional_view_requires_a_real_experimental_grade(self):
        """Annotation-only grades (A) and absent grades never open the view."""
        for grade in ("", "A", "pending", "E1"):
            with self.subTest(grade=grade):
                classified = module.classify_row(
                    build_row(experimental_evidence_grade=grade,
                              localization_evidence="other_evidence")
                )
                self.assertFalse(module.functional_view_eligible(classified))

    def test_functional_view_admits_a_fully_resolved_row(self):
        classified = module.classify_row(
            build_row(experimental_evidence_grade="E3", assignment_unique="true",
                      localization_evidence="other_evidence")
        )
        self.assertTrue(module.functional_view_eligible(classified))

    def test_functional_view_requires_explicit_localization_support(self):
        """The view may demand localization support the master catalog does not."""
        for evidence in ("unknown", "historical_family_label_only"):
            with self.subTest(localization_evidence=evidence):
                classified = module.classify_row(
                    build_row(experimental_evidence_grade="E3",
                              localization_evidence=evidence)
                )
                # The master verdict is untouched by the stricter view.
                self.assertEqual(classified["primary_disposition"],
                                 "core_sequence_homolog")
                self.assertFalse(module.functional_view_eligible(classified))

    def test_functional_view_eligibility_is_not_a_calibration_claim(self):
        """Sequence qualification and functional calibration stay separate."""
        classified = module.classify_row(
            build_row(experimental_evidence_grade="E3",
                      localization_evidence="other_evidence",
                      functional_calibration_status="not_function_calibrated")
        )
        self.assertTrue(module.functional_view_eligible(classified))
        self.assertEqual(classified["functional_calibration_status"],
                         "not_function_calibrated")

    def test_primary_disposition_does_not_depend_on_the_functional_view(self):
        """Source-level guard: the master classifier never consults the view."""
        source = Path(module.__file__).read_text(encoding="utf-8")

        def body_of(name):
            return source.split("def %s(" % name, 1)[1].split("\ndef ", 1)[0]

        for name in ("primary_disposition", "classify_row"):
            with self.subTest(function=name):
                self.assertNotIn("functional_view_eligible", body_of(name))
        # The derived view is computed only after every master verdict is final.
        catalog_body = body_of("build_catalog_rows")
        self.assertLess(
            catalog_body.index("rows.append(classified)"),
            catalog_body.index('row["functional_view_eligible"]'),
        )


# ---------------------------------------------------------------------------
# Step 4: the explicit 52-record demotion input
# ---------------------------------------------------------------------------

class DemotionWiringTests(CatalogTestCase):
    """Plan Step 4: every demotion accession is accounted for exactly once."""

    DEMOTION_COUNT = 52

    def demotion_rows_52(self, count=DEMOTION_COUNT):
        return [demotion_row(index) for index in range(1, count + 1)]

    def test_all_52_demotion_ids_are_accounted_for_exactly_once(self):
        rows = self.demotion_rows_52()
        result = self.build_catalog(rows, rows, expect_demoted=self.DEMOTION_COUNT)

        accounting = result["summary"]["demotion_accounting"]
        self.assertEqual(accounting["declared_expected"], self.DEMOTION_COUNT)
        self.assertEqual(accounting["input_ids"], self.DEMOTION_COUNT)
        self.assertEqual(accounting["accounted_ids"], self.DEMOTION_COUNT)
        self.assertEqual(accounting["missing_ids"], [])
        self.assertEqual(accounting["duplicate_ids"], [])
        self.assertTrue(accounting["exactly_once"])
        self.assertEqual(accounting["counts_match_declared"], True)

        # Every input accession appears exactly once in the written catalog.
        columns, catalog_rows = read_tsv(self.out_dir / "phaded_candidate_catalog_v2.tsv")
        accessions = [row["accession"] for row in catalog_rows]
        self.assertEqual(len(accessions), len(set(accessions)))
        expected = {row["accession"] for row in rows}
        self.assertEqual(set(accessions), expected)
        for accession in expected:
            self.assertEqual(accessions.count(accession), 1)

    def test_every_demoted_row_gets_the_discordant_flag_and_leaves_the_view(self):
        rows = self.demotion_rows_52()
        result = self.build_catalog(rows, rows, expect_demoted=self.DEMOTION_COUNT)
        _, catalog_rows = read_tsv(self.out_dir / "phaded_candidate_catalog_v2.tsv")

        self.assertEqual(len(catalog_rows), self.DEMOTION_COUNT)
        for row in catalog_rows:
            self.assertIn(module.FLAG_LOCALIZATION_DISCORDANT,
                          row["evidence_flags"].split(module.EVIDENCE_FLAG_SEPARATOR))
            self.assertEqual(row["functional_view_eligible"], "false")
        self.assertEqual(module.build_functional_view(catalog_rows), [])

    def test_demotion_count_mismatch_raises_instead_of_proceeding(self):
        rows = self.demotion_rows_52(count=self.DEMOTION_COUNT - 1)
        with self.assertRaises(module.EvidenceError) as caught:
            self.build_catalog(rows, rows, expect_demoted=self.DEMOTION_COUNT)
        message = str(caught.exception)
        self.assertIn(module.RULE_DEMOTION_COUNT_MISMATCH, message)
        self.assertIn(str(self.DEMOTION_COUNT - 1), message)
        self.assertIn(str(self.DEMOTION_COUNT), message)
        # Fail closed: nothing was written.
        self.assertFalse(self.out_dir.exists())

    def test_demotion_count_mismatch_raises_even_in_dry_run(self):
        rows = self.demotion_rows_52(count=self.DEMOTION_COUNT + 1)
        with self.assertRaises(module.EvidenceError):
            self.build_catalog(rows, rows, expect_demoted=self.DEMOTION_COUNT,
                               dry_run=True)
        self.assertFalse(self.out_dir.exists())

    def test_declared_expectation_is_overridable_with_a_stated_basis(self):
        rows = self.demotion_rows_52(count=3)
        result = self.build_catalog(rows, rows, expect_demoted=3)
        self.assertEqual(
            result["summary"]["demotion_accounting"]["declared_expected"], 3
        )
        self.assertEqual(
            result["summary"]["demotion_accounting"]["accounted_ids"], 3
        )
        self.assertEqual(
            result["summary"]["demotion_accounting"]["basis"], "explicit_expectation"
        )

    def test_module_declares_the_real_52_with_its_frozen_basis(self):
        """The real value is declared, documented and auditable."""
        self.assertEqual(module.EXPECTED_CYS_EXPORT_DEMOTIONS, 52)
        self.assertIn("20260920_phaded_three_gaps_01",
                      module.EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS)
        self.assertIn("hold_candidates_addendum_52.tsv",
                      module.EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS)
        self.assertIn("52", module.EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS)

    def test_duplicate_demotion_accessions_are_rejected(self):
        rows = self.demotion_rows_52()
        rows.append(demotion_row(1))
        with self.assertRaises(module.EvidenceError) as caught:
            self.build_catalog(rows, rows, expect_demoted=len(rows))
        self.assertIn(module.RULE_DEMOTION_ACCOUNTING, str(caught.exception))

    def test_demotion_only_accession_is_added_once_and_still_demoted(self):
        rows = [build_row(accession="ACC_CAND")]
        demotions = [demotion_row(1), demotion_row(2)]
        result = self.build_catalog(rows, demotions, expect_demoted=2)
        _, catalog_rows = read_tsv(self.out_dir / "phaded_candidate_catalog_v2.tsv")
        self.assertEqual(len(catalog_rows), 3)
        self.assertEqual(
            result["summary"]["demotion_accounting"]["accounted_ids"], 2
        )
        self.assertEqual(
            result["summary"]["demotion_accounting"]["demotion_only_ids"], 2
        )
        self.assertEqual(
            result["summary"]["universe"]["demotion_only_count"], 2
        )
        for row in catalog_rows:
            if row["accession"].startswith("DEMO_"):
                self.assertIn(module.FLAG_LOCALIZATION_DISCORDANT,
                              row["evidence_flags"])

    def test_demotion_input_supplies_the_frozen_signalp_export_prediction(self):
        """The demotion row's transport prediction is honoured, not overwritten."""
        demotions = [demotion_row(1, signalp_class="pending",
                                  signalp_class_prediction="TAT")]
        result = self.build_catalog([], demotions, expect_demoted=1)
        _, catalog_rows = read_tsv(self.out_dir / "phaded_candidate_catalog_v2.tsv")
        self.assertEqual(catalog_rows[0]["transport_signal_prediction"], "TAT")
        self.assertEqual(result["summary"]["demotion_accounting"]["accounted_ids"], 1)


# ---------------------------------------------------------------------------
# Catalog shape, validation and fail-closed behaviour
# ---------------------------------------------------------------------------

class CatalogArtifactTests(CatalogTestCase):
    def mixed_rows(self):
        return [
            build_row(accession="ACC_CORE"),
            build_row(accession="ACC_PROBABLE", assignment_unique="false"),
            build_row(accession="ACC_REMOTE", model_layer="reference_query_only"),
            build_row(accession="ACC_FUNC", superfamily=NPHAMCL),
            build_row(accession="ACC_DEFERRED", superfamily=NPHASCL_WITH_LIPASE_BOX),
            build_row(accession="ACC_EXCLUDED", input_quality_reason="sequence_missing"),
        ]

    def test_one_row_per_protein_with_exactly_one_disposition(self):
        rows = self.mixed_rows()
        self.build_catalog(rows, [])
        columns, catalog_rows = read_tsv(
            self.out_dir / "phaded_candidate_catalog_v2.tsv"
        )
        accessions = [row["accession"] for row in catalog_rows]
        self.assertEqual(len(accessions), len(set(accessions)))
        self.assertEqual(len(catalog_rows), len(rows))
        for row in catalog_rows:
            self.assertIn(row["primary_disposition"], schema.PRIMARY_DISPOSITIONS)

    def test_every_emitted_row_carries_the_required_evidence_columns(self):
        self.build_catalog(self.mixed_rows(), [])
        columns, catalog_rows = read_tsv(
            self.out_dir / "phaded_candidate_catalog_v2.tsv"
        )
        required = {
            "historical_superfamily", "sequence_family_call",
            "sequence_family_confidence", "catalytic_domain_type",
            "accessory_domain_architecture", "transport_signal_prediction",
            "localization_evidence", "nucleophile_identity", "motif_class",
            "experimental_evidence_grade", "model_layer",
            "functional_calibration_status", "primary_disposition",
            "evidence_flags", "functional_view_eligible",
        }
        self.assertTrue(required.issubset(set(columns)), sorted(required - set(columns)))
        for row in catalog_rows:
            # ``evidence_flags`` is multi-valued and non-exclusive; it is never
            # a single collapsed label.
            flags = [flag for flag in
                     row["evidence_flags"].split(module.EVIDENCE_FLAG_SEPARATOR) if flag]
            self.assertTrue(flags, "every classified row carries at least one flag")
            for flag in flags:
                self.assertIn(flag, module.EVIDENCE_FLAGS)

    def test_every_emitted_row_passes_the_shared_validator(self):
        self.build_catalog(self.mixed_rows(), [])
        _, catalog_rows = read_tsv(self.out_dir / "phaded_candidate_catalog_v2.tsv")
        for row in catalog_rows:
            schema.validate_candidate_evidence(row)

    def test_writer_fails_closed_on_an_unvalidatable_row(self):
        """A row the shared ontology rejects is never written."""
        bad = module.classify_row(build_row(accession="ACC_BAD"))
        bad["nucleophile_identity"] = "cys"
        bad["motif_class"] = "AHSMG"
        with self.assertRaises(module.EvidenceError) as caught:
            module.write_catalog(self.out_dir, [bad], {"run_id": "x"}, {"run_id": "x"})
        self.assertIn(module.RULE_SCHEMA_VALIDATION_FAILED, str(caught.exception))
        self.assertFalse(self.out_dir.exists())

    def test_outputs_are_the_three_declared_artifacts(self):
        self.build_catalog(self.mixed_rows(), [])
        written = sorted(path.name for path in self.out_dir.iterdir())
        self.assertEqual(
            written,
            ["catalog_manifest.json", "catalog_summary.json",
             "phaded_candidate_catalog_v2.tsv"],
        )

    def test_summary_counts_agree_with_the_written_catalog(self):
        self.build_catalog(self.mixed_rows(), [])
        _, catalog_rows = read_tsv(
            self.out_dir / "phaded_candidate_catalog_v2.tsv"
        )
        summary = json.loads(
            (self.out_dir / "catalog_summary.json").read_text(encoding="utf-8")
        )
        counts = summary["primary_disposition_counts"]
        self.assertEqual(sum(counts.values()), len(catalog_rows))
        for row in catalog_rows:
            self.assertGreaterEqual(counts[row["primary_disposition"]], 1)
        self.assertIn("excluded_input_quality", counts)

    def test_manifest_binds_inputs_and_declares_the_boundary(self):
        result = self.build_catalog(self.mixed_rows(), [], dry_run=True)
        manifest = result["manifest"]
        self.assertEqual(manifest["counts"]["catalog_rows"], len(self.mixed_rows()))
        self.assertIn("inputs", manifest)
        self.assertIn("candidates", manifest["inputs"])
        self.assertIn("demotion", manifest["inputs"])
        self.assertFalse(manifest["boundaries"]["phenotype_claimed"])
        self.assertFalse(manifest["boundaries"]["functional_calibration_performed"])
        self.assertTrue(manifest["boundaries"]["v1_frozen_results_untouched"])
        self.assertTrue(manifest["boundaries"]["one_disposition_per_protein"])
        self.assertIn("counts", manifest)

    def test_dry_run_validates_without_writing(self):
        result = self.build_catalog(self.mixed_rows(), [], dry_run=True)
        self.assertTrue(result["dry_run"])
        self.assertFalse(self.out_dir.exists())
        self.assertEqual(result["summary"]["universe"]["catalog_rows"],
                         len(self.mixed_rows()))

    def test_refuses_an_existing_non_empty_output_dir(self):
        self.out_dir.mkdir(parents=True)
        (self.out_dir / "existing.tsv").write_text("x\n", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            self.build_catalog(self.mixed_rows(), [])
        self.assertTrue((self.out_dir / "existing.tsv").is_file())
        self.assertEqual(sorted(p.name for p in self.out_dir.iterdir()), ["existing.tsv"])

    def test_refuses_missing_input_paths(self):
        with self.assertRaises(FileNotFoundError):
            module.build_catalog(
                candidates_path=self.tmp / "nope.tsv",
                demotion_path=self.demotion,
                output_dir=self.out_dir,
                dry_run=True,
            )

    def test_deterministic_row_order_and_repeatability(self):
        rows = self.mixed_rows()
        result = module.build_catalog_rows(rows, [], expect_demoted=0)
        self.assertEqual(
            [row["accession"] for row in result["rows"]],
            sorted(row["accession"] for row in rows),
        )
        again = module.build_catalog_rows(list(reversed(rows)), [], expect_demoted=0)
        self.assertEqual(result["rows"], again["rows"])
        self.assertEqual(result["summary"], again["summary"])


class FlagVocabularyTests(unittest.TestCase):
    def test_unknown_flag_is_rejected(self):
        with self.assertRaises(module.EvidenceError) as caught:
            module.normalize_evidence_flags("not_a_real_flag")
        self.assertIn(module.RULE_UNKNOWN_EVIDENCE_FLAG, str(caught.exception))

    def test_flags_are_deduplicated_and_canonically_ordered(self):
        flags = module.normalize_evidence_flags(
            "localization_discordant|ser_motif_lipase_box|localization_discordant"
        )
        self.assertEqual(flags.split(module.EVIDENCE_FLAG_SEPARATOR),
                         ["ser_motif_lipase_box", "localization_discordant"])

    def test_declared_vocabulary_has_no_separator_or_duplicate_hazards(self):
        self.assertEqual(len(module.EVIDENCE_FLAGS), len(set(module.EVIDENCE_FLAGS)))
        for flag in module.EVIDENCE_FLAGS:
            self.assertNotIn(module.EVIDENCE_FLAG_SEPARATOR, flag)
            self.assertNotIn(" ", flag)

    def test_input_quality_flags_are_never_scientific(self):
        for flag in module.INPUT_QUALITY_FLAGS:
            self.assertIn(flag, module.EVIDENCE_FLAGS)
        for flag in module.SCIENTIFIC_UNCERTAINTY_FLAGS:
            self.assertIn(flag, module.EVIDENCE_FLAGS)
        self.assertFalse(set(module.INPUT_QUALITY_FLAGS)
                         & set(module.SCIENTIFIC_UNCERTAINTY_FLAGS))


class SchemaVocabularyReuseTests(unittest.TestCase):
    """The catalog must reuse the shared ontology, never restate it."""

    def test_dispositions_come_from_the_schema_module(self):
        for disposition in schema.PRIMARY_DISPOSITIONS:
            self.assertIn(disposition, module.PRIMARY_DISPOSITIONS)

    def test_motif_and_nucleophile_columns_use_schema_vocabularies(self):
        out = module.classify_row(build_row(ahsmg_state="supported", superfamily=PHAZ7))
        self.assertIn(out["motif_class"], schema.MOTIF_CLASSES)
        self.assertIn(out["nucleophile_identity"], schema.NUCLEOPHILE_IDENTITIES)
        self.assertIn(out["transport_signal_prediction"],
                      schema.TRANSPORT_SIGNAL_PREDICTIONS)
        self.assertIn(out["localization_evidence"], schema.LOCALIZATION_EVIDENCE)
        self.assertIn(out["catalytic_domain_type"], schema.CATALYTIC_DOMAIN_TYPES)
        self.assertIn(out["model_layer"], schema.MODEL_LAYERS)

    def test_schema_module_is_imported_not_duplicated(self):
        self.assertTrue(hasattr(module, "SCHEMA"))
        self.assertIs(module.SCHEMA.PRIMARY_DISPOSITIONS, schema.PRIMARY_DISPOSITIONS)


class CliTests(CatalogTestCase):
    def test_main_prints_a_json_summary_and_writes_the_catalog(self):
        import contextlib
        import io

        rows = [build_row(accession="ACC_CORE")]
        write_tsv(self.candidates, module.CANDIDATE_COLUMNS, rows)
        write_tsv(self.demotion, DEMOTION_COLUMNS, [])
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = module.main([
                "--candidates", str(self.candidates),
                "--demotion", str(self.demotion),
                "--output-dir", str(self.out_dir),
                "--expect-demoted", "0",
            ])
        self.assertEqual(code, 0)
        payload = json.loads(buffer.getvalue())
        self.assertEqual(payload["counts"]["catalog_rows"], 1)
        self.assertTrue((self.out_dir / "phaded_candidate_catalog_v2.tsv").is_file())

    def test_main_states_the_basis_for_the_declared_demotion_expectation(self):
        """The CLI default records the frozen 52 with its documented basis."""
        import contextlib
        import io

        rows = [demotion_row(index) for index in range(1, 53)]
        write_tsv(self.candidates, module.CANDIDATE_COLUMNS, rows)
        write_tsv(self.demotion, DEMOTION_COLUMNS, rows)
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = module.main([
                "--candidates", str(self.candidates),
                "--demotion", str(self.demotion),
                "--output-dir", str(self.out_dir),
            ])
        self.assertEqual(code, 0)
        accounting = json.loads(buffer.getvalue())["demotion_accounting"]
        self.assertEqual(accounting["basis"], module.BASIS_DECLARED_DEFAULT)
        self.assertEqual(accounting["basis_detail"],
                         module.EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS)
        self.assertEqual(accounting["declared_expected"], 52)
        self.assertTrue(accounting["exactly_once"])

    def test_a_stated_expectation_records_its_own_basis(self):
        """Overriding the count is allowed, but never silently."""
        rows = [demotion_row(1)]
        write_tsv(self.demotion, DEMOTION_COLUMNS, rows)
        result = module.build_catalog_rows(
            [build_row(accession="ACC_CORE")], rows, expect_demoted=1,
        )
        accounting = result["summary"]["demotion_accounting"]
        self.assertEqual(accounting["basis"], module.BASIS_EXPLICIT_EXPECTATION)
        self.assertEqual(accounting["declared_expected"], 1)
        self.assertNotEqual(accounting["basis_detail"],
                            module.EXPECTED_CYS_EXPORT_DEMOTIONS_BASIS)

    def test_main_returns_nonzero_on_a_failing_input(self):
        write_tsv(self.candidates, module.CANDIDATE_COLUMNS, [build_row()])
        write_tsv(self.demotion, DEMOTION_COLUMNS, [demotion_row(1)])
        code = module.main([
            "--candidates", str(self.candidates),
            "--demotion", str(self.demotion),
            "--output-dir", str(self.out_dir),
            "--expect-demoted", "52",
        ])
        self.assertEqual(code, 1)
        self.assertFalse(self.out_dir.exists())

    def test_main_reports_missing_input_paths(self):
        code = module.main([
            "--candidates", str(self.tmp / "missing.tsv"),
            "--demotion", str(self.demotion),
            "--output-dir", str(self.out_dir),
        ])
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
