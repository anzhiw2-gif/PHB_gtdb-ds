import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_phaded_panel_acquisition_03.py"


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_phaded_panel_acquisition_03", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PanelAcquisition03Tests(unittest.TestCase):
    def test_exact_ledger_binding_is_allowed_but_does_not_create_family(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            ledger = Path(temporary) / "ledger.tsv"
            ledger.write_text("accession\tphaded_family_id\nAAF86381.1\tDED_hfam_55\n", encoding="utf-8")
            bound = module.bind_existing_family("AAF86381.1", ledger)
            self.assertEqual(bound, ("resolved_existing_family", "DED_hfam_55"))

    def test_unresolved_external_candidate_cannot_be_formal_negative(self):
        module = load_module()
        row = module.classify_candidate("family_resolved_negative", "unresolved", "experimental_negative")
        self.assertEqual(row["role"], "family_resolved_negative")
        self.assertEqual(row["eligibility"], "not_eligible_unresolved_family")
        self.assertEqual(row["new_family_call_created"], "false")

    def test_training_accession_cannot_be_heldout_positive(self):
        module = load_module()
        row = module.classify_candidate("independent_positive", "resolved_existing_family", "experimental_positive", training_accession=True)
        self.assertEqual(row["eligibility"], "not_eligible_training_accession")

    def test_write_panel_marks_training_accession(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "panel.tsv"
            module.write_panel(out, [{
                "accession": "X1.1", "role": "independent_positive", "evidence_status": "experimental_positive",
                "family_binding_status": "resolved_existing_family", "bound_family": "DED_hfam_3",
                "training_accession": "true",
            }])
            row = next(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
            self.assertEqual(row["training_accession"], "true")
            self.assertEqual(row["eligibility"], "not_eligible_training_accession")

    def test_manifest_preserves_independence_and_provenance(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "panel.tsv"
            rows = [{
                "accession": "X1.1", "role": "independent_positive", "evidence_status": "experimental_positive",
                "family_binding_status": "unresolved", "pmid": "123", "source_database": "NCBI Protein",
                "sequence_sha256": "a" * 64,
            }]
            module.write_panel(out, rows)
            read = list(csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
            self.assertEqual(read[0]["independence_status"], "pending_family_resolution")
            self.assertEqual(read[0]["new_family_call_created"], "false")


if __name__ == "__main__":
    unittest.main()
