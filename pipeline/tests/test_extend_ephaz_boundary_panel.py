import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.extend_ephaz_boundary_panel import extend_panel


class ExtendEphazBoundaryPanelTests(unittest.TestCase):
    def test_appends_unique_boundary_records_and_writes_ledger(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel = root / "intracellular_non_ephaz_negative.faa"
            additions = root / "boundary.faa"
            ledger = root / "boundary.tsv"
            panel.write_text(">O87189|existing\nAAAA\n", encoding="utf-8")
            additions.write_text(">Q5Y152|mcl_boundary\nCCCC\n>Q4W8C9|phb_state_boundary\nGGGG\n", encoding="utf-8")
            rows = extend_panel(panel, additions, ledger, "intracellular_non_ephaz_negative")
            self.assertEqual(rows, ["Q4W8C9", "Q5Y152"])
            text = panel.read_text(encoding="utf-8")
            self.assertIn(">Q4W8C9|phb_state_boundary", text)
            self.assertIn("Q5Y152\tintracellular_non_ephaz_negative", ledger.read_text(encoding="utf-8"))

    def test_rejects_duplicate_accession(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel = root / "panel.faa"
            additions = root / "boundary.faa"
            panel.write_text(">Q5Y152|existing\nAAAA\n", encoding="utf-8")
            additions.write_text(">Q5Y152|duplicate\nCCCC\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                extend_panel(panel, additions, root / "ledger.tsv", "panel")

    def test_parses_uniprot_prefixed_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel = root / "panel.faa"
            additions = root / "boundary.faa"
            panel.write_text(">O87189|existing\nAAAA\n", encoding="utf-8")
            additions.write_text(">sp|Q4W8C9|3HBOH\nCCCC\n", encoding="utf-8")
            self.assertEqual(extend_panel(panel, additions, root / "ledger.tsv", "panel"), ["Q4W8C9"])


if __name__ == "__main__":
    unittest.main()
