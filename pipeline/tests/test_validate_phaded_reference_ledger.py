import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_phaded_reference_ledger.py"
EVIDENCE_COLUMNS_PATH = Path(__file__).resolve().parents[1] / "config" / "phaded_reference_evidence_columns.tsv"


FIELDS = [
    "reference_id", "accession", "sequence_sha256", "sequence_source", "source_database",
    "source_database_version", "retrieval_date", "organism", "taxonomy_id",
    "phaded_superfamily", "phaded_family_id", "reported_localization", "substrate_class",
    "substrate_detail", "experimental_assay", "experimental_result", "evidence_status",
    "positive_negative_control", "catalytic_residues", "lipase_box_or_ahsmg",
    "oxyanion_hole_evidence", "domain_architecture", "primary_doi", "pmid", "pmcid", "notes",
]

# The v2 evidence-quality columns (Task 4). The frozen 723-row ledger predates them, so
# every fixture that exercises them is synthetic and in-memory/tempfile only.
EVIDENCE_FIELDS = [
    "accession_version", "gene_id", "study_id", "experiment_unit_id", "independence_group",
    "experimental_evidence_grade", "assay_directness", "experimental_substrate_class",
    "experimental_system", "negative_control_description", "catalytic_mutant_evidence",
    "localization_evidence", "sequence_integrity", "discovery_training_eligible",
    "functional_calibration_eligible",
]

ALL_FIELDS = FIELDS + EVIDENCE_FIELDS
AMINO_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# Grade-aware defaults so a fixture that only names the grade is still self-consistent.
_GRADE_DEFAULTS = {
    "E3": {"assay_directness": "direct"},
    "E2": {"assay_directness": "direct"},
    "E1": {"assay_directness": "indirect"},
    "A": {
        "assay_directness": "annotation", "evidence_status": "annotation_only",
        "experimental_assay": "pending", "experimental_result": "pending",
        "study_id": "pending", "experiment_unit_id": "pending", "independence_group": "pending",
        "negative_control_description": "not_assessed", "catalytic_mutant_evidence": "not_assessed",
        "localization_evidence": "not_assessed",
    },
}


def reference_sequence(accession, length=21):
    digest = hashlib.sha256(accession.encode("utf-8")).hexdigest()
    body = "".join(
        AMINO_ALPHABET[int(digest[index] + digest[index + 1], 16) % len(AMINO_ALPHABET)]
        for index in range(length - 1)
    )
    return "M" + body


def valid_row(**overrides):
    """A synthetic row that satisfies every declared value unless overridden."""
    grade = overrides.get("experimental_evidence_grade", "E3")
    row = {field: "" for field in ALL_FIELDS}
    row.update({
        "reference_id": "ref-1", "accession": "P12345",
        "sequence_source": "local-test", "source_database": "UniProtKB",
        "source_database_version": "2026-09-09", "retrieval_date": "2026-09-09",
        "organism": "Escherichia coli", "taxonomy_id": "562",
        "phaded_superfamily": "extracellular dPHASCL type 1", "phaded_family_id": "family-test-01",
        "reported_localization": "extracellular", "substrate_class": "PHB",
        "substrate_detail": "native PHB granules",
        "experimental_assay": "purified enzyme PHB hydrolysis", "experimental_result": "positive",
        "evidence_status": "experimental_positive",
        "positive_negative_control": "boiled enzyme control",
        "catalytic_residues": "Ser139-Asp191-His228", "lipase_box_or_ahsmg": "GxSxG",
        "oxyanion_hole_evidence": "backbone amides confirmed",
        "domain_architecture": "catalytic domain + SBD",
        "primary_doi": "10.1186/1471-2105-10-89", "pmid": "19296857", "pmcid": "PMC2666664",
        "notes": "synthetic test row",
        "accession_version": "P12345.1", "gene_id": "phaZ",
        "study_id": "STUDY-001", "experiment_unit_id": "UNIT-001",
        "independence_group": "INDEP-001", "experimental_evidence_grade": "E3",
        "assay_directness": "direct", "experimental_substrate_class": "PHB",
        "experimental_system": "purified enzyme, in vitro",
        "negative_control_description": "boiled enzyme control",
        "catalytic_mutant_evidence": "S139A abolished hydrolysis",
        "localization_evidence": "experimental", "sequence_integrity": "complete",
        "discovery_training_eligible": "false", "functional_calibration_eligible": "false",
    })
    row.update(_GRADE_DEFAULTS.get(grade, {}))
    row.update(overrides)
    if "accession_version" not in overrides:
        row["accession_version"] = row["accession"] + ".1"
    if "sequence_sha256" not in overrides:
        row["sequence_sha256"] = hashlib.sha256(
            reference_sequence(row["accession"]).encode("ascii")
        ).hexdigest()
    return row


def write_evidence_ledger(root, rows, name="ledger_evidence.tsv", fasta_name="reference_evidence.faa"):
    ledger = root / name
    with ledger.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ALL_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in ALL_FIELDS})
    fasta = root / fasta_name
    lines = []
    for row in rows:
        lines.append(">" + row["reference_id"] + "|" + row["accession"])
        lines.append(reference_sequence(row["accession"]))
    fasta.write_text("\n".join(lines) + "\n", encoding="ascii")
    return ledger, fasta


def validate_rows(root, module, rows, **kwargs):
    ledger, fasta = write_evidence_ledger(root, rows)
    return module.validate_ledger(
        ledger, write_family_definitions(root), fasta, expected_family_count=1, **kwargs
    )


def _load_module():
    spec = importlib.util.spec_from_file_location("validate_phaded_reference_ledger", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_family_definitions(root):
    # Test-only one-family fixture; production validation defaults to 38.
    path = root / "families.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["phaded_superfamily", "phaded_family_id", "source_reference", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({
            "phaded_superfamily": "extracellular dPHASCL type 1",
            "phaded_family_id": "family-test-01",
            "source_reference": "Knoll2009:PMC2666664",
            "status": "source_bound",
        })
    return path


def write_reference_fasta(root, sequence="M" + "A" * 19):
    path = root / "reference.faa"
    path.write_text(">ref-1|P12345\n" + sequence + "\n", encoding="ascii")
    return path, sequence


def write_ledger(root, accession="P12345", primary_doi="10.1186/1471-2105-10-89", evidence_status="experimental_positive", strict_training="no"):
    path = root / "ledger.tsv"
    fasta, sequence = write_reference_fasta(root)
    row = {field: "" for field in FIELDS}
    row.update({
        "reference_id": "ref-1", "accession": accession,
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "sequence_source": "local-test", "source_database": "UniProtKB",
        "source_database_version": "2026-09-09", "retrieval_date": "2026-09-09",
        "organism": "Test organism", "taxonomy_id": "123",
        "phaded_superfamily": "extracellular dPHASCL type 1", "phaded_family_id": "family-test-01",
        "experimental_assay": "purified enzyme PHB hydrolysis", "experimental_result": "positive",
        "evidence_status": evidence_status, "primary_doi": primary_doi,
        "notes": "strict_training=" + strict_training,
    })
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    return path, fasta


class PhaDEDLedgerTests(unittest.TestCase):
    def test_rejects_reference_without_accession_or_primary_identifier(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root, accession="", primary_doi="")
            with self.assertRaisesRegex(ValueError, "accession"):
                module.validate_ledger(ledger, write_family_definitions(root), fasta, expected_family_count=1)

    def test_rejects_annotation_only_record_from_strict_training(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root, evidence_status="annotation_only", strict_training="yes")
            with self.assertRaisesRegex(ValueError, "strict-training"):
                module.validate_ledger(
                    ledger, write_family_definitions(root), fasta,
                    strict_training=True, expected_family_count=1,
                )

    def test_rejects_fasta_to_ledger_mismatch(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root)
            fasta.write_text(">other|P12345\nMAAAAAAAAAAAAAAAAAAA\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "FASTA-to-ledger mismatch"):
                module.validate_ledger(ledger, write_family_definitions(root), fasta, expected_family_count=1)

    def test_rejects_unbound_family_assignment(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root)
            text = ledger.read_text(encoding="utf-8").replace("family-test-01", "family-not-bound")
            ledger.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not source-bound"):
                module.validate_ledger(ledger, write_family_definitions(root), fasta, expected_family_count=1)

    def test_rejects_invalid_family_definition_or_unexpected_cardinality(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root)
            definitions = write_family_definitions(root)
            with self.assertRaisesRegex(ValueError, "expected 38 unique family IDs"):
                module.validate_ledger(ledger, definitions, fasta)

            definitions.write_text(
                definitions.read_text(encoding="utf-8").replace(
                    "extracellular dPHASCL type 1", "unsupported superfamily"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "invalid PhaDED superfamily"):
                module.validate_ledger(ledger, definitions, fasta, expected_family_count=1)

    def test_rejects_strict_training_row_with_missing_source_provenance(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root)
            text = ledger.read_text(encoding="utf-8").replace("local-test", "")
            ledger.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "sequence_source"):
                module.validate_ledger(
                    ledger, write_family_definitions(root), fasta,
                    strict_training=True, expected_family_count=1,
                )

    def test_rejects_fasta_header_with_mismatched_accession(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root)
            fasta.write_text(">ref-1|Q99999\nMAAAAAAAAAAAAAAAAAAA\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "FASTA accession mismatch"):
                module.validate_ledger(ledger, write_family_definitions(root), fasta, expected_family_count=1)

    # --- Task 4: evidence quality, grade provenance, and independence -----------------

    def test_direct_experimental_record_requires_complete_provenance(self):
        module = _load_module()
        row = valid_row(evidence_status="experimental_positive", experimental_evidence_grade="E3")
        row["study_id"] = ""
        with self.assertRaisesRegex(ValueError, "E3.*study_id"):
            module.validate_reference_row(row)

    def test_annotation_only_record_can_be_discovery_training_eligible(self):
        module = _load_module()
        row = valid_row(evidence_status="annotation_only", experimental_evidence_grade="A")
        row["discovery_training_eligible"] = "true"
        module.validate_reference_row(row)

    def test_rejects_unknown_evidence_grade(self):
        module = _load_module()
        row = valid_row(experimental_evidence_grade="E4")
        with self.assertRaisesRegex(ValueError, "experimental_evidence_grade"):
            module.validate_reference_row(row)

    def test_rejects_unknown_assay_directness_and_sequence_integrity(self):
        module = _load_module()
        with self.assertRaisesRegex(ValueError, "assay_directness"):
            module.validate_reference_row(valid_row(assay_directness="presumed"))
        with self.assertRaisesRegex(ValueError, "sequence_integrity"):
            module.validate_reference_row(valid_row(sequence_integrity="mostly_complete"))

    def test_rejects_e2_record_without_experiment_unit(self):
        module = _load_module()
        row = valid_row(experimental_evidence_grade="E2")
        row["experiment_unit_id"] = "pending"
        with self.assertRaisesRegex(ValueError, "E2.*experiment_unit_id"):
            module.validate_reference_row(row)

    def test_rejects_direct_grade_without_direct_assay_or_evidence_source(self):
        module = _load_module()
        row = valid_row(experimental_evidence_grade="E3", assay_directness="indirect")
        with self.assertRaisesRegex(ValueError, "assay_directness"):
            module.validate_reference_row(row)
        row = valid_row(experimental_evidence_grade="E3")
        for column in ("primary_doi", "pmid", "pmcid"):
            row[column] = ""
        with self.assertRaisesRegex(ValueError, "evidence source"):
            module.validate_reference_row(row)

    def test_e1_record_is_audited_but_never_a_direct_functional_positive(self):
        module = _load_module()
        row = valid_row(
            experimental_evidence_grade="E1", assay_directness="direct", independence_group="INDEP-E1"
        )
        module.validate_reference_row(row)
        summary = module.summarize_independence([row])
        self.assertEqual(summary["functional_positive_count"], 0)
        self.assertEqual(summary["functional_positive_row_count"], 0)
        self.assertEqual(summary["grade_counts"]["E1"], 1)
        self.assertEqual(summary["e1_rows_with_direct_assay"], ["P12345"])
        self.assertEqual(summary["audited_row_count"], 1)

    def test_discovery_eligibility_requires_complete_sequence_and_family_consistency(self):
        module = _load_module()
        row = valid_row(
            evidence_status="annotation_only", experimental_evidence_grade="A",
            sequence_integrity="truncated", discovery_training_eligible="true",
        )
        with self.assertRaisesRegex(ValueError, "sequence_integrity"):
            module.validate_reference_row(row)
        row = valid_row(
            evidence_status="annotation_only", experimental_evidence_grade="A",
            phaded_family_id="pending", discovery_training_eligible="true",
        )
        with self.assertRaisesRegex(ValueError, "phaded_family_id"):
            module.validate_reference_row(row)

    def test_functional_calibration_eligibility_is_forbidden_for_grade_a(self):
        module = _load_module()
        row = valid_row(
            evidence_status="annotation_only", experimental_evidence_grade="A",
            functional_calibration_eligible="true",
        )
        with self.assertRaisesRegex(ValueError, "functional_calibration_eligible"):
            module.validate_reference_row(row)

    def test_contradictory_direct_negative_blocks_calibration_but_keeps_the_row_audited(self):
        module = _load_module()
        row = valid_row(
            evidence_status="experimental_negative", experimental_evidence_grade="E3",
            experimental_result="no PHB hydrolysis detected", functional_calibration_eligible="true",
        )
        with self.assertRaisesRegex(ValueError, "functional_calibration_eligible"):
            module.validate_reference_row(row)
        row["functional_calibration_eligible"] = "false"
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            report = validate_rows(root, module, [row])
        summary = report["independence_summary"]
        self.assertEqual(report["status"], "valid")
        self.assertEqual(summary["audited_row_count"], 1)
        self.assertIn("P12345", summary["audited_accessions"])
        self.assertEqual(summary["contradictory_direct_negative_accessions"], ["P12345"])
        self.assertEqual(summary["functional_calibration_blocked_accessions"], ["P12345"])
        self.assertEqual(summary["functional_positive_count"], 0)

    def test_same_independence_group_counts_as_one_functional_positive(self):
        module = _load_module()
        first = valid_row(accession="A00001", reference_id="ref-a", independence_group="INDEP-X")
        second = valid_row(accession="A00002", reference_id="ref-b", independence_group="INDEP-X")
        summary = module.summarize_independence([first, second])
        self.assertEqual(summary["functional_positive_row_count"], 2)
        self.assertEqual(summary["functional_positive_count"], 1)
        self.assertEqual(summary["functional_positive_independence_groups"], ["INDEP-X"])
        self.assertEqual(summary["functional_positive_accessions"], ["A00001", "A00002"])
        self.assertEqual(summary["audited_row_count"], 2)
        self.assertEqual(summary["counts_by_independence_group"]["INDEP-X"], 2)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            report = validate_rows(root, module, [first, second])
        self.assertEqual(report["independence_summary"]["functional_positive_count"], 1)
        self.assertEqual(report["functional_positive_count"], 1)
        self.assertEqual(report["functional_positive_row_count"], 2)

    def test_summarize_independence_reports_every_required_dimension(self):
        module = _load_module()
        rows = [
            valid_row(accession="A00001", reference_id="ref-a", gene_id="g1", study_id="S1",
                      experiment_unit_id="U1", independence_group="G1", organism="Escherichia coli"),
            valid_row(accession="A00002", reference_id="ref-b", gene_id="g1", study_id="S1",
                      experiment_unit_id="U1", independence_group="G1", organism="Escherichia coli"),
            valid_row(accession="A00003", reference_id="ref-c", gene_id="g2", study_id="S2",
                      experiment_unit_id="U2", independence_group="G2", organism="Bacillus subtilis",
                      experimental_evidence_grade="E2"),
        ]
        summary = module.summarize_independence(rows)
        for key in (
            "counts_by_accession", "counts_by_accession_version", "counts_by_sequence_sha256",
            "counts_by_gene", "counts_by_experiment_unit", "counts_by_study", "counts_by_genus",
            "counts_by_independence_group",
        ):
            self.assertIn(key, summary)
        self.assertEqual(summary["counts_by_gene"]["g1"], 2)
        self.assertEqual(summary["counts_by_experiment_unit"]["U1"], 2)
        self.assertEqual(summary["counts_by_study"]["S1"], 2)
        self.assertEqual(summary["counts_by_genus"]["Escherichia"], 2)
        self.assertEqual(summary["counts_by_independence_group"]["G1"], 2)
        self.assertEqual(len(summary["counts_by_sequence_sha256"]), 3)
        self.assertEqual(summary["unique_accessions"], 3)
        self.assertEqual(summary["unique_genes"], 2)
        self.assertEqual(summary["functional_positive_count"], 2)
        self.assertEqual(summary["duplicate_counts"]["gene"], {"g1": 2})

    def test_report_records_the_functional_calibration_gate_inputs(self):
        module = _load_module()
        rows = [
            valid_row(accession="A00001", reference_id="ref-a", independence_group="G1",
                      study_id="S1", experiment_unit_id="U1"),
            valid_row(accession="A00002", reference_id="ref-b", independence_group="G2",
                      study_id="S1", experiment_unit_id="U2"),
            valid_row(accession="A00003", reference_id="ref-c", independence_group="G3",
                      study_id="S2", experiment_unit_id="U3"),
        ]
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            report = validate_rows(root, module, rows)
        gate = report["functional_calibration_gate"]
        self.assertEqual(gate["minimum_independent_positive_count"], 3)
        self.assertEqual(gate["positive_count_basis"], "unique_independence_group")
        self.assertEqual(gate["family_functional_positive_counts"]["family-test-01"], 3)
        self.assertEqual(gate["family_distinct_genus_counts"]["family-test-01"], 1)
        self.assertEqual(gate["families_meeting_minimum"], ["family-test-01"])
        self.assertEqual(gate["families_below_minimum"], [])

    def test_legacy_ledger_without_evidence_columns_still_validates_in_compat_mode(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root)
            definitions = write_family_definitions(root)
            report = module.validate_ledger(ledger, definitions, fasta, expected_family_count=1)
            self.assertEqual(report["status"], "valid")
            self.assertEqual(report["row_count"], 1)
            self.assertEqual(report["evidence_column_mode"], "compat")
            self.assertEqual(report["evidence_columns_declared"], [])
            self.assertEqual(report["pending_evidence_columns"], sorted(module.EVIDENCE_COLUMNS))
            self.assertEqual(report["independence_summary"]["functional_positive_count"], 0)
            self.assertEqual(
                report["independence_summary"]["uncurated_experimental_positive_accessions"], ["P12345"]
            )
            with self.assertRaisesRegex(ValueError, "missing required evidence columns"):
                module.validate_ledger(
                    ledger, definitions, fasta, expected_family_count=1, evidence_column_mode="strict"
                )
            with self.assertRaisesRegex(ValueError, "missing required evidence columns"):
                module.validate_ledger(
                    ledger, definitions, fasta, expected_family_count=1, strict_evidence_columns=True
                )

    def test_strict_evidence_column_mode_requires_curated_values(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            row = valid_row(
                experimental_evidence_grade="pending", assay_directness="pending",
                sequence_integrity="pending",
            )
            ledger, fasta = write_evidence_ledger(root, [row])
            definitions = write_family_definitions(root)
            with self.assertRaisesRegex(ValueError, "experimental_evidence_grade"):
                module.validate_ledger(
                    ledger, definitions, fasta, expected_family_count=1, evidence_column_mode="strict"
                )
            report = module.validate_ledger(ledger, definitions, fasta, expected_family_count=1)
            self.assertEqual(report["status"], "valid")

    def test_strict_evidence_column_mode_accepts_a_fully_curated_ledger(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            report = validate_rows(
                root, module, [valid_row(experimental_evidence_grade="E3")],
                evidence_column_mode="strict",
            )
            self.assertEqual(report["status"], "valid")
            self.assertEqual(report["functional_positive_count"], 1)

    def test_rejects_unknown_evidence_column_mode(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_ledger(root)
            with self.assertRaisesRegex(ValueError, "evidence_column_mode"):
                module.validate_ledger(
                    ledger, write_family_definitions(root), fasta,
                    expected_family_count=1, evidence_column_mode="whatever",
                )

    def test_evidence_columns_config_matches_the_validator_vocabulary(self):
        module = _load_module()
        declared = module.load_evidence_columns()
        self.assertEqual(list(declared), list(module.EVIDENCE_COLUMNS))
        self.assertTrue(EVIDENCE_COLUMNS_PATH.is_file())
        config_lines = [
            line.strip()
            for line in EVIDENCE_COLUMNS_PATH.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        self.assertEqual(config_lines, list(module.EVIDENCE_COLUMNS))

    def test_load_evidence_columns_rejects_a_drifted_config(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            config = root / "columns.tsv"
            config.write_text("\n".join(sorted(module.EVIDENCE_COLUMNS)[:-1]) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "evidence columns config"):
                module.load_evidence_columns(config)

    def test_cli_prints_the_independence_summary_as_json(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_evidence_ledger(root, [valid_row()])
            definitions = write_family_definitions(root)
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = module.main([
                    str(ledger), str(definitions), str(fasta),
                    "--expected-family-count", "1", "--independence-summary-json",
                ])
        self.assertEqual(exit_code, 0)
        payload = json.loads(stdout.getvalue())
        self.assertEqual(payload["functional_positive_count"], 1)
        self.assertEqual(payload["counts_by_independence_group"]["INDEP-001"], 1)

    def test_cli_writes_the_independence_summary_to_a_file(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            ledger, fasta = write_evidence_ledger(root, [valid_row()])
            definitions = write_family_definitions(root)
            target = root / "summary" / "independence.json"
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                exit_code = module.main([
                    str(ledger), str(definitions), str(fasta),
                    "--expected-family-count", "1",
                    "--independence-summary-json", str(target),
                ])
            self.assertEqual(exit_code, 0)
            self.assertTrue(target.is_file())
            payload = json.loads(target.read_text(encoding="utf-8"))
        self.assertEqual(payload["functional_positive_count"], 1)


if __name__ == "__main__":
    unittest.main()
