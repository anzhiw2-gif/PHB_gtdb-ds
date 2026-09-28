import tempfile
import unittest
from pathlib import Path

from pipeline.scripts import materialize_phaded_ded_references as module


class MaterializePhaDEDTests(unittest.TestCase):
    def test_extracts_exact_family_count_and_deduplicates_sequences(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "family_index.html").write_text(
                '<a href="index.pl?page=hfam&id=61">i-nPHAscl (no lipase box) homologous family 1</a>\n'
                '<a href="index.pl?page=hfam&id=62">i-nPHAscl (no lipase box) homologous family 2</a>\n',
                encoding="latin1",
            )
            (root / "seq61.txt").write_text(">A1\nMKK\n>A2\nMKK\n", encoding="ascii")
            (root / "seq62.txt").write_text(">B1\nMNN\n", encoding="ascii")
            (root / "family_61.html").write_text("<html></html>", encoding="latin1")
            (root / "family_62.html").write_text("<html></html>", encoding="latin1")
            output = root / "out"
            result = module.materialize(root, output, expected_family_count=2, seed_gis=set())
            self.assertEqual(result["family_count"], 2)
            self.assertEqual(result["reference_count"], 2)
            self.assertTrue((output / "phaded_reference.faa").is_file())


if __name__ == "__main__":
    unittest.main()
