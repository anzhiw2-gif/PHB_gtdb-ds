import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_phaded_interpro_result.py"


def load_module():
    spec = importlib.util.spec_from_file_location("validate_phaded_interpro_result", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def interpro_row(accession: str, sequence: str) -> list[str]:
    return [
        accession,
        hashlib.md5(sequence.encode("ascii")).hexdigest(),
        str(len(sequence)),
        "Pfam",
        "PF10503",
        "Esterase PHB depolymerase",
        "1",
        str(len(sequence)),
        "1.0E-20",
        "T",
        "14-09-2026",
        "IPR010126",
        "Esterase, PHB depolymerase",
        "-",
        "-",
    ]


class ValidateInterProResultTests(unittest.TestCase):
    def test_validates_sequence_binding_and_excluded_accounting(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            accessions = ["candidate_1", "candidate_2", "excluded_1"]
            (root / "accessions.txt").write_text("\n".join(accessions) + "\n", encoding="ascii")
            (root / "input.faa").write_text(">candidate_1\nMAAA\n>candidate_2\nMCCC\n", encoding="ascii")
            with (root / "excluded.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["accession", "reason"])
                writer.writerow(["excluded_1", "invalid_amino_acid"])
            with (root / "result.tsv").open("w", encoding="ascii", newline="") as handle:
                handle.write("\t".join(interpro_row("candidate_1", "MAAA")) + "\n")
                handle.write("\t".join(interpro_row("candidate_2", "MCCC")) + "\n")
            report = module.validate_paths(root / "result.tsv", root / "input.faa", root / "accessions.txt", root / "excluded.tsv")
            self.assertEqual(report["status"], "validated_candidate_only")
            self.assertEqual(report["result_accession_count"], 2)
            self.assertEqual(report["missing_excluded_count"], 1)
            self.assertEqual(report["field_error_count"], 0)

    def test_detects_md5_and_coordinate_errors(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "accessions.txt").write_text("candidate_1\n", encoding="ascii")
            (root / "input.faa").write_text(">candidate_1\nMAAA\n", encoding="ascii")
            (root / "excluded.tsv").write_text("accession\treason\n", encoding="ascii")
            row = interpro_row("candidate_1", "MAAA")
            row[1] = "0" * 32
            row[6] = "0"
            row[7] = "9"
            with (root / "result.tsv").open("w", encoding="ascii", newline="") as handle:
                handle.write("\t".join(row) + "\n")
            report = module.validate_paths(root / "result.tsv", root / "input.faa", root / "accessions.txt", root / "excluded.tsv")
            self.assertGreaterEqual(report["md5_mismatch_count"], 1)
            self.assertGreaterEqual(report["coordinate_error_count"], 1)
            self.assertEqual(report["status"], "validation_failed")


if __name__ == "__main__":
    unittest.main()
