import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_ephaz_bridge_ledger.py"


class BridgeLedgerValidatorTests(unittest.TestCase):
    def _load_module(self):
        spec = importlib.util.spec_from_file_location("validate_ephaz_bridge_ledger", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def _write_ledger(self, root, row):
        path = root / "bridge.tsv"
        fields = [
            "candidate_id", "accession", "protein_name", "organism", "sequence_source",
            "sequence_url", "sequence_length", "primary_reference", "pmid", "doi",
            "experimental_evidence_type", "substrate", "extracellular_evidence",
            "architecture_evidence", "completeness_status", "core_overlap_status",
            "iPhaZ_challenge_status", "decision", "decision_reason", "review_date", "notes",
        ]
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(row)
        return path

    def test_accepts_only_traceable_experimental_complete_bridge_rows(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            row = {
                "candidate_id": "bridge_1", "accession": "Q51871", "protein_name": "PHB depolymerase A",
                "organism": "Paucimonas lemoignei", "sequence_source": "UniProt",
                "sequence_url": "https://rest.uniprot.org/uniprotkb/Q51871.fasta", "sequence_length": "433",
                "primary_reference": "J Bacteriol 1995", "pmid": "7836292", "doi": "10.1128/jb.177.3.596-607.1995",
                "experimental_evidence_type": "cloned_gene_and_extracellular_enzyme_characterization",
                "substrate": "PHB", "extracellular_evidence": "paper_reports_extracellular_localization; UniProt_signal_1_26",
                "architecture_evidence": "Esterase_PHB; alpha_beta_hydrolase", "completeness_status": "complete",
                "core_overlap_status": "not_in_current_core", "iPhaZ_challenge_status": "not_challenge",
                "decision": "accept_bridge_candidate", "decision_reason": "all_inclusion_criteria_met",
                "review_date": "2026-08-30", "notes": "",
            }
            ledger = self._write_ledger(Path(tmp), row)
            accepted = module.validate_ledger(ledger)
            self.assertEqual(accepted, ["Q51871"])

    def test_rejects_accepted_row_without_primary_experimental_evidence(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            row = {
                "candidate_id": "bridge_2", "accession": "Q5SLU4", "protein_name": "Carboxylesterase",
                "organism": "Thermus thermophilus", "sequence_source": "UniProt",
                "sequence_url": "https://rest.uniprot.org/uniprotkb/Q5SLU4.fasta", "sequence_length": "515",
                "primary_reference": "", "pmid": "", "doi": "",
                "experimental_evidence_type": "annotation_only", "substrate": "PHB",
                "extracellular_evidence": "signal", "architecture_evidence": "carboxylesterase",
                "completeness_status": "complete", "core_overlap_status": "not_in_current_core",
                "iPhaZ_challenge_status": "not_challenge", "decision": "accept_bridge_candidate",
                "decision_reason": "incorrect", "review_date": "2026-08-30", "notes": "",
            }
            ledger = self._write_ledger(Path(tmp), row)
            with self.assertRaises(module.LedgerError):
                module.validate_ledger(ledger)


if __name__ == "__main__":
    unittest.main()
