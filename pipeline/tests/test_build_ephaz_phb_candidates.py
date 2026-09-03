import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.build_ephaz_phb_candidates import build_candidates


class BuildEphazPhbCandidatesTests(unittest.TestCase):
    def test_writes_exclusion_and_independent_variants(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            core = root / "core.faa"
            independent = root / "independent.faa"
            out = root / "out"
            core.write_text(">B2NHN2|core\nAAAA\n>Q51718|mcl\nCCCC\n", encoding="utf-8")
            independent.write_text(">AAB40611.1|phb\nGGGG\n>O24719|phb\nTTTT\n", encoding="utf-8")
            result = build_candidates(core, independent, out, ["AAB40611.1", "O24719"])
            self.assertEqual(result["PHB-focused_no_Q51718"], ["B2NHN2"])
            self.assertEqual(result["PHB-focused_plus_independent"], ["AAB40611.1", "B2NHN2", "O24719"])
            self.assertNotIn("Q51718", out.joinpath("PHB-focused_no_Q51718.faa").read_text())

    def test_rejects_missing_independent_accession(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            core = root / "core.faa"
            independent = root / "independent.faa"
            core.write_text(">B2NHN2|core\nAAAA\n", encoding="utf-8")
            independent.write_text(">AAB40611.1|phb\nGGGG\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                build_candidates(core, independent, root / "out", ["AAB40611.1", "O24719"])


if __name__ == "__main__":
    unittest.main()
