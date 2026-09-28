import importlib.util
import tempfile
import unittest

from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "extract_priority_fasta.py"


def load_module():
    spec = importlib.util.spec_from_file_location("extract_priority_fasta", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ExtractPriorityFastaTests(unittest.TestCase):
    def test_extract_strips_only_terminal_stop_and_rejects_internal_invalid_characters(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta = root / "all.faa"
            fasta.write_text(">A\nMAA*\n>B\nMA*AA\n", encoding="ascii")
            output = root / "selected.faa"
            with self.assertRaisesRegex(ValueError, "invalid amino acid"):
                module.extract(fasta, ["B"], output, sanitize_terminal_stops=True)
            fasta.write_text(">A\nMAA*\n", encoding="ascii")
            count = module.extract(fasta, ["A"], output, sanitize_terminal_stops=True)
            self.assertEqual(count, 1)
            self.assertEqual(output.read_text(encoding="ascii"), ">A\nMAA\n")

    def test_extract_can_record_and_exclude_non_iupac_records(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta = root / "all.faa"
            fasta.write_text(">A\nMAA*\n>B\nMAX\n", encoding="ascii")
            output = root / "selected.faa"
            excluded = root / "excluded.tsv"
            count = module.extract(fasta, ["A", "B"], output, sanitize_terminal_stops=True, exclude_invalid=excluded)
            self.assertEqual(count, 1)
            self.assertIn("B\tinvalid_amino_acid", excluded.read_text(encoding="utf-8"))

    def test_extract_preserves_queue_order_and_rejects_missing_accession(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta = root / "all.faa"
            fasta.write_text(">B annotation\nMCCC\n>A\nMAAA\n", encoding="ascii")
            output = root / "selected.faa"
            count = module.extract(fasta, ["A", "B"], output)
            self.assertEqual(count, 2)
            self.assertEqual(output.read_text(encoding="ascii"), ">A\nMAAA\n>B annotation\nMCCC\n")
            with self.assertRaisesRegex(ValueError, "missing from FASTA"):
                module.extract(fasta, ["MISSING"], output)


if __name__ == "__main__":
    unittest.main()
