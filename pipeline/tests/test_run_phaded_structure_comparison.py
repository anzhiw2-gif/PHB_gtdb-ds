import tempfile
import unittest
from pathlib import Path

from pipeline.scripts import run_phaded_structure_comparison as module


class StructureComparisonTests(unittest.TestCase):
    def test_foldseek_alignment_fields_match_requested_columns(self):
        self.assertEqual(len(module.ALIGN_FIELDS), 17)

    def test_normalizes_utf8_bom_on_first_header(self):
        self.assertEqual(module.normalize_fieldnames(["\ufeffcandidate_id", "reference_pdb"])[0], "candidate_id")

    def test_rejects_server_threads_above_forty(self):
        with self.assertRaisesRegex(ValueError, "40"):
            module.validate_threads(41, server=True)

    def test_requires_existing_pair_files(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(FileNotFoundError, "candidate"):
                module.validate_pair({"candidate_pdb": str(Path(temp) / "missing"), "reference_pdb": str(Path(temp) / "ref")})


if __name__ == "__main__":
    unittest.main()
