"""Tests for the shared PhaDED evidence ontology (controlled vocabulary + validators)."""

import contextlib
import csv
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "phaded_evidence_schema.py"
GATES = Path(__file__).resolve().parents[1] / "config" / "phaded_sequence_model_gates.tsv"

GATE_FIELDS = [
    "gate_version", "min_training_sequences", "min_heldout_sequences", "min_heldout_recall",
    "max_unexplained_confounder_hits", "require_unique_best_model", "require_alignment_hash",
]
GATE_ROW = dict(zip(
    GATE_FIELDS,
    ["sequence-gate-v1", "3", "1", "0.8", "0", "true", "true"],
))

REFERENCE_FIELDS = [
    "reference_id", "accession", "sequence_sha256", "phaded_superfamily", "phaded_family_id",
    "reported_localization", "catalytic_residues", "lipase_box_or_ahsmg", "evidence_status",
    "primary_doi",
]


def _load_module():
    spec = importlib.util.spec_from_file_location("phaded_evidence_schema", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_tsv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def _candidate_row(**overrides):
    row = {
        "accession": "P12345",
        "model_layer": "sequence_family_hmm_validated",
        "functional_calibration_status": "not_function_calibrated",
        "primary_disposition": "core_sequence_homolog",
        "sequence_family_call": "hfam_70",
        "experimental_evidence_grade": "A",
        "catalytic_domain_type": "type1_verified",
        "transport_signal_prediction": "SP",
        "localization_evidence": "experimental_confirmed",
        "localization_evidence_source": "experimental_assay",
        "nucleophile_identity": "ser",
        "motif_class": "GxSxG",
        "evidence_flags": "",
    }
    row.update(overrides)
    return row


def _reference_row(**overrides):
    row = {
        "reference_id": "ref-1", "accession": "P12345",
        "sequence_sha256": "a" * 64,
        "phaded_superfamily": "extracellular native-SCL/PhaZ7-like",
        "phaded_family_id": "family-test-01",
        "reported_localization": "extracellular",
        "catalytic_residues": "Ser190, Asp224, His266",
        "lipase_box_or_ahsmg": "AHSMG",
        "evidence_status": "annotation_only",
        "primary_doi": "10.1186/1471-2105-10-89",
    }
    row.update(overrides)
    return row


def _write_gates(root, **overrides):
    row = dict(GATE_ROW)
    row.update(overrides)
    fieldnames = [field for field in GATE_FIELDS if field in row]
    fieldnames += [field for field in row if field not in GATE_FIELDS]
    return _write_tsv(root / "gates.tsv", fieldnames, [row])


def _run_main(argv):
    module = _load_module()
    stdout, stderr = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        code = module.main(argv)
    return code, stdout.getvalue(), stderr.getvalue()


class PhaDEDEvidenceVocabularyTests(unittest.TestCase):
    def test_controlled_vocabularies_match_the_preregistered_design(self):
        module = _load_module()
        self.assertEqual(module.MODEL_LAYERS, {
            "reference_query_only",
            "discovery_hmm_uncalibrated",
            "sequence_family_hmm_validated",
            "calibrated_candidate_model",
        })
        self.assertEqual(module.NUCLEOPHILE_IDENTITIES, {"ser", "cys", "unresolved"})
        self.assertEqual(module.MOTIF_CLASSES, {"GxSxG", "AHSMG", "Cys-associated", "other", "unresolved"})
        self.assertEqual(module.EVIDENCE_GRADES, {"E3", "E2", "E1", "A"})
        self.assertEqual(module.PRIMARY_DISPOSITIONS, {
            "core_sequence_homolog",
            "probable_sequence_homolog",
            "remote_homolog_candidate",
            "function_unresolved",
            "deferred_structure_review",
            "excluded_input_quality",
        })
        self.assertEqual(module.FUNCTIONAL_CALIBRATION_STATUSES, {
            "not_function_calibrated",
            "candidate_gate_passed_not_promoted",
            "calibrated_candidate_model",
            "blocked_contradictory_experimental_negative",
        })
        self.assertEqual(module.TRANSPORT_SIGNAL_PREDICTIONS, {
            "SP", "LIPO", "TAT", "TATLIPO", "OTHER", "not_tested", "invalid_input",
        })
        self.assertEqual(module.LOCALIZATION_EVIDENCE, {
            "experimental_confirmed", "other_evidence", "historical_family_label_only", "unknown",
        })
        self.assertEqual(module.CATALYTIC_DOMAIN_TYPES, {
            "type1_verified", "type2_verified", "other", "unresolved",
        })

    def test_ahsmg_is_a_motif_class_and_never_a_nucleophile_identity(self):
        module = _load_module()
        self.assertIn("AHSMG", module.MOTIF_CLASSES)
        self.assertNotIn("ahsmg", module.NUCLEOPHILE_IDENTITIES)

    def test_small_heldout_panel_rule_is_recorded_for_the_evaluator(self):
        module = _load_module()
        self.assertEqual(module.SMALL_HELDOUT_PANEL_THRESHOLD, 5)
        self.assertIn("smaller than five", module.__doc__)
        self.assertIn("Task 6", module.__doc__)


class PhaDEDNucleophileTests(unittest.TestCase):
    def test_ahsmg_is_a_ser_motif_not_a_nucleophile_identity(self):
        module = _load_module()
        row = module.normalize_nucleophile("extracellular native-SCL/PhaZ7-like", "", "supported")
        self.assertEqual(row, {"nucleophile_identity": "ser", "motif_class": "AHSMG"})

    def test_motif_state_may_name_ahsmg_directly(self):
        module = _load_module()
        row = module.normalize_nucleophile("extracellular native-SCL/PhaZ7-like", "AHSMG", "")
        self.assertEqual(row, {"nucleophile_identity": "ser", "motif_class": "AHSMG"})

    def test_canonical_ser_motif_is_ser(self):
        module = _load_module()
        row = module.normalize_nucleophile("extracellular dPHASCL type 1", "GxSxG", "")
        self.assertEqual(row, {"nucleophile_identity": "ser", "motif_class": "GxSxG"})

    def test_cys_associated_motif_is_cys_even_in_an_scl_superfamily(self):
        module = _load_module()
        row = module.normalize_nucleophile("intracellular nPHASCL without lipase box", "Cys-associated", "")
        self.assertEqual(row, {"nucleophile_identity": "cys", "motif_class": "Cys-associated"})

    def test_superfamily_prior_applies_only_without_informative_motif_evidence(self):
        module = _load_module()
        self.assertEqual(
            module.normalize_nucleophile("intracellular nPHASCL with lipase box", "", ""),
            {"nucleophile_identity": "ser", "motif_class": "unresolved"},
        )
        self.assertEqual(
            module.normalize_nucleophile("intracellular nPHAMCL", "", ""),
            {"nucleophile_identity": "unresolved", "motif_class": "unresolved"},
        )
        self.assertEqual(
            module.normalize_nucleophile("", "", ""),
            {"nucleophile_identity": "unresolved", "motif_class": "unresolved"},
        )

    def test_unknown_motif_state_is_rejected_naming_field_and_value(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "motif_state='banana'"):
            module.normalize_nucleophile("extracellular dPHASCL type 1", "banana", "")

    def test_unknown_ahsmg_state_is_rejected_naming_field_and_value(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "ahsmg_state='maybe'"):
            module.normalize_nucleophile("extracellular dPHASCL type 1", "", "maybe")


class PhaDEDLocalizationTests(unittest.TestCase):
    def test_signalp_does_not_become_localization_truth(self):
        module = _load_module()
        row = module.normalize_localization("OTHER", "historical_family_label_only")
        self.assertEqual(row["transport_signal_prediction"], "OTHER")
        self.assertEqual(row["localization_evidence"], "historical_family_label_only")

    def test_normalize_localization_never_infers_localization_from_signalp(self):
        module = _load_module()
        row = module.normalize_localization("SP", "")
        self.assertEqual(row, {"transport_signal_prediction": "SP", "localization_evidence": "unknown"})

    def test_normalize_localization_defaults_without_inputs(self):
        module = _load_module()
        self.assertEqual(
            module.normalize_localization("", ""),
            {"transport_signal_prediction": "not_tested", "localization_evidence": "unknown"},
        )

    def test_normalize_localization_accepts_every_controlled_class(self):
        module = _load_module()
        for prediction in sorted(module.TRANSPORT_SIGNAL_PREDICTIONS):
            with self.subTest(prediction=prediction):
                row = module.normalize_localization(prediction.lower(), "unknown")
                self.assertEqual(row["transport_signal_prediction"], prediction)
        for evidence in sorted(module.LOCALIZATION_EVIDENCE):
            with self.subTest(evidence=evidence):
                row = module.normalize_localization("not_tested", evidence.upper())
                self.assertEqual(row["localization_evidence"], evidence)

    def test_normalize_localization_rejects_unknown_classes(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "signalp_class='SECRETED'"):
            module.normalize_localization("SECRETED", "unknown")
        with self.assertRaisesRegex(ValueError, "localization_evidence='cytoplasmic_guess'"):
            module.normalize_localization("SP", "cytoplasmic_guess")


class PhaDEDModelStateTests(unittest.TestCase):
    def test_model_layers_are_disjoint_from_functional_status(self):
        module = _load_module()
        module.validate_model_state("sequence_family_hmm_validated", "not_function_calibrated")

    def test_unknown_model_layer_is_rejected_naming_field_and_value(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "model_layer='discovery_hmm'"):
            module.validate_model_state("discovery_hmm", "not_function_calibrated")

    def test_unknown_functional_status_is_rejected_naming_field_and_value(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "functional_calibration_status='calibrated'"):
            module.validate_model_state("sequence_family_hmm_validated", "calibrated")

    def test_calibrated_status_requires_the_calibrated_layer(self):
        module = _load_module()
        for layer in ("discovery_hmm_uncalibrated", "reference_query_only", "sequence_family_hmm_validated"):
            with self.subTest(layer=layer):
                with self.assertRaisesRegex(ValueError, "calibrated_candidate_model"):
                    module.validate_model_state(layer, "calibrated_candidate_model")

    def test_calibrated_layer_requires_the_calibrated_status(self):
        module = _load_module()
        module.validate_model_state("calibrated_candidate_model", "calibrated_candidate_model")
        with self.assertRaisesRegex(ValueError, "functional_calibration_status"):
            module.validate_model_state("calibrated_candidate_model", "not_function_calibrated")

    def test_discovery_layer_may_pass_the_candidate_gate_without_being_calibrated(self):
        module = _load_module()
        module.validate_model_state("discovery_hmm_uncalibrated", "candidate_gate_passed_not_promoted")


class PhaDEDCandidateEvidenceTests(unittest.TestCase):
    def test_candidate_row_read_from_tsv_is_accepted_and_normalized(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = _write_tsv(Path(temporary_directory) / "candidates.tsv", list(_candidate_row()), [_candidate_row()])
            with path.open(encoding="utf-8-sig", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            normalized = module.validate_candidate_evidence(row)
        self.assertEqual(normalized["model_layer"], "sequence_family_hmm_validated")
        self.assertEqual(normalized["functional_calibration_status"], "not_function_calibrated")
        self.assertEqual(normalized["primary_disposition"], "core_sequence_homolog")
        self.assertEqual(normalized["nucleophile_identity"], "ser")
        self.assertEqual(normalized["motif_class"], "GxSxG")
        self.assertEqual(normalized["transport_signal_prediction"], "SP")
        self.assertEqual(normalized["localization_evidence"], "experimental_confirmed")
        self.assertEqual(normalized["localization_evidence_source"], "experimental_assay")

    def test_candidate_defaults_are_the_weakest_layer(self):
        module = _load_module()
        normalized = module.validate_candidate_evidence({})
        self.assertEqual(normalized["model_layer"], "reference_query_only")
        self.assertEqual(normalized["functional_calibration_status"], "not_function_calibrated")
        self.assertEqual(normalized["transport_signal_prediction"], "not_tested")
        self.assertEqual(normalized["localization_evidence"], "unknown")
        self.assertEqual(normalized["nucleophile_identity"], "unresolved")
        self.assertEqual(normalized["motif_class"], "unresolved")
        self.assertEqual(normalized["catalytic_domain_type"], "unresolved")
        self.assertEqual(normalized["primary_disposition"], "")
        self.assertEqual(normalized["experimental_evidence_grade"], "")

    def test_candidate_validator_rejects_ahsmg_nucleophile_identity(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "AHSMG is a Ser motif class"):
            module.validate_candidate_evidence(_candidate_row(nucleophile_identity="ahsmg"))

    def test_candidate_validator_rejects_unknown_nucleophile_identity(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "nucleophile_identity='asp'"):
            module.validate_candidate_evidence(_candidate_row(nucleophile_identity="asp"))

    def test_candidate_validator_rejects_cys_identity_with_ahsmg_motif_class(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "motif_class='AHSMG'"):
            module.validate_candidate_evidence(_candidate_row(nucleophile_identity="cys", motif_class="AHSMG"))

    def test_signalp_only_provenance_cannot_claim_experimental_localization(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "localization_evidence_source='signalp_only'"):
            module.validate_candidate_evidence(_candidate_row(localization_evidence_source="signalp_only"))
        module.validate_candidate_evidence(_candidate_row(localization_evidence_source="experimental_assay"))
        module.validate_candidate_evidence(_candidate_row(
            localization_evidence="other_evidence", localization_evidence_source="signalp_only",
        ))

    def test_experimental_localization_requires_a_stated_source(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "localization_evidence_source"):
            module.validate_candidate_evidence(_candidate_row(localization_evidence_source=""))

    def test_discovery_layer_row_with_a_family_call_is_rejected(self):
        module = _load_module()
        for field in ("sequence_family_call", "family_call"):
            row = _candidate_row(
                model_layer="discovery_hmm_uncalibrated",
                primary_disposition="probable_sequence_homolog",
                sequence_family_call="",
            )
            row[field] = "hfam_70"
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, field):
                    module.validate_candidate_evidence(row)

    def test_discovery_layer_row_may_score_without_a_family_call(self):
        module = _load_module()
        normalized = module.validate_candidate_evidence(_candidate_row(
            model_layer="discovery_hmm_uncalibrated",
            primary_disposition="probable_sequence_homolog",
            sequence_family_call="",
        ))
        self.assertEqual(normalized["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(normalized["primary_disposition"], "probable_sequence_homolog")
        self.assertEqual(normalized["sequence_family_call"], "")

    def test_unexplained_discovery_exclusion_is_rejected(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "excluded_input_quality"):
            module.validate_candidate_evidence(_candidate_row(
                model_layer="discovery_hmm_uncalibrated",
                primary_disposition="excluded_input_quality",
                sequence_family_call="",
            ))

    def test_documented_input_quality_exclusion_is_accepted(self):
        module = _load_module()
        normalized = module.validate_candidate_evidence(_candidate_row(
            primary_disposition="excluded_input_quality",
            sequence_family_call="",
            input_quality_reason="sequence_truncated",
        ))
        self.assertEqual(normalized["primary_disposition"], "excluded_input_quality")
        self.assertEqual(normalized["input_quality_reason"], "sequence_truncated")

    def test_scientific_uncertainty_is_not_an_input_error(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "scientific uncertainty"):
            module.validate_candidate_evidence(_candidate_row(
                primary_disposition="excluded_input_quality",
                sequence_family_call="",
                input_quality_reason="function_unresolved",
            ))

    def test_discovery_layer_cannot_claim_a_core_sequence_homolog(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "core_sequence_homolog"):
            module.validate_candidate_evidence(_candidate_row(
                model_layer="discovery_hmm_uncalibrated",
                primary_disposition="core_sequence_homolog",
                sequence_family_call="",
            ))

    def test_candidate_validator_rejects_unknown_disposition_and_transport_class(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "primary_disposition='high_confidence'"):
            module.validate_candidate_evidence(_candidate_row(primary_disposition="high_confidence"))
        with self.assertRaisesRegex(ValueError, "transport_signal_prediction='SECRETED'"):
            module.validate_candidate_evidence(_candidate_row(transport_signal_prediction="SECRETED"))

    def test_candidate_validator_ignores_unrelated_and_missing_columns(self):
        module = _load_module()
        row = _candidate_row(notes=None)
        row["future_column"] = "whatever"
        module.validate_candidate_evidence(row)


class PhaDEDReferenceEvidenceTests(unittest.TestCase):
    def test_reference_row_without_extended_columns_is_accepted(self):
        module = _load_module()
        normalized = module.validate_reference_evidence(_reference_row())
        self.assertEqual(normalized["reference_id"], "ref-1")
        self.assertEqual(normalized["model_layer"], "reference_query_only")
        self.assertEqual(normalized["functional_calibration_status"], "not_function_calibrated")
        self.assertEqual(normalized["experimental_evidence_grade"], "")
        self.assertEqual(normalized["nucleophile_identity"], "unresolved")
        self.assertEqual(normalized["motif_class"], "unresolved")
        self.assertEqual(normalized["localization_evidence"], "unknown")

    def test_reference_validator_rejects_ahsmg_nucleophile_identity(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "AHSMG is a Ser motif class"):
            module.validate_reference_evidence(_reference_row(nucleophile_identity="ahsmg"))

    def test_reference_validator_rejects_cys_identity_against_ahsmg_ledger_column(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "lipase_box_or_ahsmg"):
            module.validate_reference_evidence(_reference_row(nucleophile_identity="cys"))
        normalized = module.validate_reference_evidence(_reference_row(nucleophile_identity="ser"))
        self.assertEqual(normalized["nucleophile_identity"], "ser")

    def test_reference_signalp_only_localization_claim_is_rejected(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "localization_evidence_source='signalp_only'"):
            module.validate_reference_evidence(_reference_row(
                localization_evidence="experimental_confirmed",
                localization_evidence_source="signalp_only",
            ))
        normalized = module.validate_reference_evidence(_reference_row(
            localization_evidence="experimental_confirmed",
            localization_evidence_source="experimental_assay",
        ))
        self.assertEqual(normalized["localization_evidence"], "experimental_confirmed")

    def test_reference_validator_rejects_invalid_evidence_grade(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "experimental_evidence_grade='E4'"):
            module.validate_reference_evidence(_reference_row(experimental_evidence_grade="E4"))

    def test_reference_functional_calibration_requires_e2_or_e3(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "functional_calibration_eligible"):
            module.validate_reference_evidence(_reference_row(
                experimental_evidence_grade="E1", functional_calibration_eligible="true",
            ))
        module.validate_reference_evidence(_reference_row(
            experimental_evidence_grade="A", discovery_training_eligible="true",
        ))
        normalized = module.validate_reference_evidence(_reference_row(
            experimental_evidence_grade="E3", functional_calibration_eligible="true",
        ))
        self.assertEqual(normalized["functional_calibration_eligible"], "true")

    def test_reference_validator_checks_model_state_when_present(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "model_layer='discovery_hmm'"):
            module.validate_reference_evidence(_reference_row(model_layer="discovery_hmm"))
        with self.assertRaisesRegex(ValueError, "calibrated_candidate_model"):
            module.validate_reference_evidence(_reference_row(
                model_layer="reference_query_only", functional_calibration_status="calibrated_candidate_model",
            ))


class PhaDEDSequenceModelGateTests(unittest.TestCase):
    def test_gate_loader_reads_the_preregistered_row(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            gates = module.load_sequence_model_gates(_write_gates(Path(temporary_directory)))
        self.assertEqual(gates, GATE_ROW)

    def test_real_gate_config_is_preregistered(self):
        module = _load_module()
        self.assertEqual(module.SEQUENCE_MODEL_GATES_PATH, GATES)
        self.assertTrue(GATES.is_file())
        self.assertEqual(module.load_sequence_model_gates(GATES), GATE_ROW)
        self.assertEqual(module.SEQUENCE_MODEL_GATE_VERSION, "sequence-gate-v1")

    def test_gate_loader_rejects_missing_or_extra_columns(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            missing = _write_gates(root, gate_version=None)
            text = missing.read_text(encoding="utf-8").replace("gate_version\t", "", 1)
            missing.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "column"):
                module.load_sequence_model_gates(missing)
            extra = _write_gates(root, extra_gate="true")
            with self.assertRaisesRegex(ValueError, "column"):
                module.load_sequence_model_gates(extra)

    def test_gate_loader_rejects_malformed_files(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            cases = {
                "wrong gate version": {"gate_version": "sequence-gate-v2"},
                "non-integer minimum": {"min_training_sequences": "three"},
                "zero training sequences": {"min_training_sequences": "0"},
                "zero held-out sequences": {"min_heldout_sequences": "0"},
                "recall above one": {"min_heldout_recall": "1.5"},
                "recall of zero": {"min_heldout_recall": "0"},
                "non-numeric recall": {"min_heldout_recall": "high"},
                "non-boolean flag": {"require_alignment_hash": "yes"},
                "negative confounder budget": {"max_unexplained_confounder_hits": "-1"},
            }
            for label, overrides in cases.items():
                with self.subTest(case=label):
                    with self.assertRaises(ValueError):
                        module.load_sequence_model_gates(_write_gates(root, **overrides))

            empty = root / "empty.tsv"
            empty.write_text("", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "column"):
                module.load_sequence_model_gates(empty)

            duplicate = _write_gates(root)
            _write_tsv(duplicate, GATE_FIELDS, [GATE_ROW, GATE_ROW])
            with self.assertRaisesRegex(ValueError, "exactly one"):
                module.load_sequence_model_gates(duplicate)

    def test_gate_loader_rejects_a_missing_or_symlinked_file(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            with self.assertRaisesRegex(ValueError, "not a regular file"):
                module.load_sequence_model_gates(root / "absent.tsv")


class PhaDEDEvidenceSchemaCliTests(unittest.TestCase):
    def test_cli_validates_candidate_rows_and_prints_a_json_summary(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = _write_tsv(
                Path(temporary_directory) / "candidates.tsv",
                list(_candidate_row()),
                [_candidate_row(), _candidate_row(accession="Q99999")],
            )
            code, stdout, stderr = _run_main([str(path)])
        self.assertEqual(code, 0)
        self.assertEqual(stderr, "")
        self.assertEqual(json.loads(stdout), {"rows": 2, "valid": 2, "invalid": 0})

    def test_cli_stops_at_the_first_violation_with_a_readable_message(self):
        rows = [
            _candidate_row(),
            _candidate_row(nucleophile_identity="ahsmg"),
            _candidate_row(primary_disposition="high_confidence"),
        ]
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = _write_tsv(Path(temporary_directory) / "candidates.tsv", list(rows[0]), rows)
            code, stdout, stderr = _run_main([str(path)])
        self.assertEqual(code, 1)
        self.assertIn("nucleophile_identity", stderr)
        self.assertIn(":3:", stderr)
        self.assertNotIn(":4:", stderr)
        self.assertEqual(json.loads(stdout), {"rows": 3, "valid": 1, "invalid": 1})

    def test_cli_keep_going_reports_every_violation(self):
        rows = [
            _candidate_row(),
            _candidate_row(nucleophile_identity="ahsmg"),
            _candidate_row(primary_disposition="high_confidence"),
        ]
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = _write_tsv(Path(temporary_directory) / "candidates.tsv", list(rows[0]), rows)
            code, stdout, stderr = _run_main([str(path), "--keep-going"])
        self.assertEqual(code, 1)
        self.assertIn(":3:", stderr)
        self.assertIn(":4:", stderr)
        self.assertEqual(json.loads(stdout), {"rows": 3, "valid": 1, "invalid": 2})

    def test_cli_validates_reference_rows(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = _write_tsv(Path(temporary_directory) / "ledger.tsv", REFERENCE_FIELDS, [_reference_row()])
            code, stdout, stderr = _run_main([str(path), "--kind", "reference"])
        self.assertEqual(code, 0)
        self.assertEqual(stderr, "")
        self.assertEqual(json.loads(stdout), {"rows": 1, "valid": 1, "invalid": 0})

    def test_cli_can_validate_the_gate_configuration(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            good = _write_gates(root)
            code, stdout, stderr = _run_main([str(good), "--kind", "candidate", "--gates", str(good)])
            self.assertEqual(code, 0)
            self.assertEqual(stderr, "")
            self.assertEqual(json.loads(stdout), {"rows": 1, "valid": 1, "invalid": 0})
            bad = _write_gates(root, gate_version="sequence-gate-v2")
            code, _, stderr = _run_main([str(good), "--gates", str(bad)])
        self.assertEqual(code, 1)
        self.assertIn("gate_version", stderr)

    def test_cli_rejects_unreadable_input(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            code, _, stderr = _run_main([str(Path(temporary_directory) / "absent.tsv")])
        self.assertEqual(code, 2)
        self.assertIn("not a regular file", stderr)


if __name__ == "__main__":
    unittest.main()
