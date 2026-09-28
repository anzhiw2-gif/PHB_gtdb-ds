import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.audit_phaded_dphASCL2_structure import summarize_pae, summarize_pdb


class AuditDphASCL2StructureTests(unittest.TestCase):
    def test_domain_and_nterm_bands_are_separated(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "x.pdb"
            lines = ["MODEL     1"]
            for i, score in enumerate([40.0, 50.0, 80.0, 90.0], 1):
                lines.append(f"ATOM      1  CA  ALA A{i:4d}      0.000   0.000   0.000  1.00{score:6.2f}           C")
            pdb.write_text("\n".join(lines) + "\nENDMDL\n", encoding="ascii")
            result = summarize_pdb(pdb, sequence_length=4, domain_start=3, domain_end=4, nterm_end=2)
            self.assertEqual(result["domain_residue_count"], 2)
            self.assertEqual(result["domain_plddt_mean"], 85.0)
            self.assertEqual(result["nterm_plddt_mean"], 45.0)

    def test_low_confidence_runs_report_residue_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "x.pdb"
            lines = ["MODEL     1"]
            for i, score in enumerate([40.0, 50.0, 80.0, 60.0], 1):
                lines.append(f"ATOM      1  CA  ALA A{i:4d}      0.000   0.000   0.000  1.00{score:6.2f}           C")
            pdb.write_text("\n".join(lines) + "\nENDMDL\n", encoding="ascii")
            result = summarize_pdb(pdb, sequence_length=4, domain_start=3, domain_end=4, nterm_end=2)
            self.assertEqual(result["low_confidence_runs"], "1-2;4-4")

    def test_pae_reports_domain_boundary_and_internal_error(self):
        matrix = [[0.0, 4.0, 20.0, 21.0], [4.0, 0.0, 18.0, 19.0], [20.0, 18.0, 0.0, 5.0], [21.0, 19.0, 5.0, 0.0]]
        result = summarize_pae(matrix, sequence_length=4, domain_start=3, domain_end=4, nterm_end=2)
        self.assertEqual(result["pae_domain_internal_mean"], 5.0)
        self.assertEqual(result["pae_nterm_domain_mean"], 19.5)


if __name__ == "__main__":
    unittest.main()
