import csv
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_phb_literature_ledger.py"
FIELDS = [
    "record_id", "organism_name", "strain_name", "protein_accession", "gene_accession",
    "genome_accession", "ncbi_taxid", "gtdb_taxonomy", "source_database",
    "source_database_version", "retrieval_date", "primary_doi", "pmid", "pmcid",
    "experimental_substrate", "experimental_result", "evidence_type", "evidence_status",
    "exact_identifier_scope", "sequence_sha256", "notes",
]


def load_module():
    spec = importlib.util.spec_from_file_location("validate_phb_literature_ledger", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_fixture(root, **changes):
    sequence = "M" + "A" * 9
    row = {field: "pending" for field in FIELDS}
    row.update({
        "record_id": "phb-001", "organism_name": "Test bacterium", "strain_name": "T1",
        "protein_accession": "P12345", "gene_accession": "NC_000001.1",
        "genome_accession": "GCA_000001.1", "ncbi_taxid": "123",
        "source_database": "UniProtKB", "source_database_version": "pending",
        "retrieval_date": "2026-09-13", "pmid": "12345678",
        "experimental_substrate": "PHB", "experimental_result": "hydrolysis observed",
        "evidence_type": "direct_experiment", "evidence_status": "experimental_positive",
        "exact_identifier_scope": "protein_accession", "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
    })
    row.update(changes)
    path = root / "ledger.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    return path


class ValidatePhbLiteratureLedgerTests(unittest.TestCase):
    def test_validates_fixture_and_counts_status(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as directory:
            report = module.validate_ledger(write_fixture(Path(directory)))
            self.assertEqual(report["row_count"], 1)
            self.assertEqual(report["evidence_status_counts"]["experimental_positive"], 1)

    def test_rejects_duplicate_accession(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = write_fixture(root)
            text = path.read_text(encoding="utf-8")
            duplicate = text.splitlines()[1].replace("phb-001", "phb-002", 1)
            path.write_text(text + duplicate + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate protein_accession"):
                module.validate_ledger(path)

    def test_rejects_missing_citation_for_non_pending_record(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as directory:
            path = write_fixture(Path(directory), pmid="", pmcid="", primary_doi="")
            with self.assertRaisesRegex(ValueError, "citation"):
                module.validate_ledger(path)

    def test_rejects_unsupported_status(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as directory:
            path = write_fixture(Path(directory), evidence_status="reported")
            with self.assertRaisesRegex(ValueError, "evidence_status"):
                module.validate_ledger(path)

    def test_rejects_invalid_sequence_hash(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as directory:
            path = write_fixture(Path(directory), sequence_sha256="not-a-hash")
            with self.assertRaisesRegex(ValueError, "sequence_sha256"):
                module.validate_ledger(path)


if __name__ == "__main__":
    unittest.main()
