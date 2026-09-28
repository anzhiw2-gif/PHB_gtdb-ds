import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "assemble_gtdb_phb_evidence.py"
SPEC = importlib.util.spec_from_file_location("assemble_gtdb_phb_evidence", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
try:
    SPEC.loader.exec_module(module)
except FileNotFoundError:
    module = None


def write_tsv(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8")


class AssembleGtdbPhbEvidenceTests(unittest.TestCase):
    def test_genome_tier_uses_exact_literature_bridge_and_keeps_holds(self):
        if module is None:
            self.fail("assemble_gtdb_phb_evidence.py is not implemented")
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            assignment = root / "assignment.tsv"
            feature = root / "feature.tsv"
            literature = root / "literature.tsv"
            manual = root / "manual.tsv"
            output = root / "out.tsv"
            write_tsv(assignment, ["accession", "genome", "phaded_superfamily_best", "phaded_family_best", "assignment_status"], [
                ["pA", "G_A", "extracellular dPHASCL type 1", "family_1", "assigned"],
                ["pB", "G_B", "extracellular dPHASCL type 1", "family_1", "assigned"],
                ["pC", "G_C", "", "", "unassigned_PhaDED_like"],
            ])
            write_tsv(feature, ["accession", "architecture_consistency", "signalp_class", "sequence_integrity"], [
                ["pA", "consistent", "SP", "valid"],
                ["pB", "partial", "OTHER", "valid"],
                ["pC", "pending", "OTHER", "valid"],
            ])
            write_tsv(literature, ["genome_accession", "evidence_status", "exact_identifier_scope"], [["G_A", "experimental_positive", "exact_genome"]])
            write_tsv(manual, ["candidate_id", "manual_review_decision"], [["pB", "hold_architecture_conflict"]])
            rows = module.assemble_genome_evidence(assignment, feature, literature, manual, output)
            self.assertEqual(rows["G_A"]["evidence_tier"], "A")
            self.assertEqual(rows["G_B"]["evidence_tier"], "H")
            self.assertEqual(rows["G_C"]["evidence_tier"], "C")
            self.assertEqual(output.exists(), True)

    def test_missing_literature_is_not_a_negative(self):
        if module is None:
            self.fail("assemble_gtdb_phb_evidence.py is not implemented")
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            assignment = root / "a.tsv"
            feature = root / "f.tsv"
            literature = root / "l.tsv"
            manual = root / "m.tsv"
            output = root / "o.tsv"
            write_tsv(assignment, ["accession", "genome", "phaded_superfamily_best", "phaded_family_best", "assignment_status"], [["p", "G", "sf", "fam", "assigned"]])
            write_tsv(feature, ["accession", "architecture_consistency", "signalp_class", "sequence_integrity"], [["p", "consistent", "OTHER", "valid"]])
            write_tsv(literature, ["genome_accession", "evidence_status", "exact_identifier_scope"], [])
            write_tsv(manual, ["candidate_id", "manual_review_decision"], [])
            rows = module.assemble_genome_evidence(assignment, feature, literature, manual, output)
            self.assertEqual(rows["G"]["evidence_tier"], "B")
            self.assertEqual(rows["G"]["literature_status"], "not_found")


if __name__ == "__main__":
    unittest.main()
