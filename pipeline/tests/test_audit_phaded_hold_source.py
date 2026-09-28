import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.audit_phaded_hold_source import summarize_gbff


class AuditPhadedHoldSourceTests(unittest.TestCase):
    def test_exact_translation_match_and_partial_location(self):
        with tempfile.TemporaryDirectory() as tmp:
            gbff = Path(tmp) / "x.gbff"
            gbff.write_text(
                "LOCUS       TEST                 12 bp    DNA     linear   01-JAN-2020\n"
                "FEATURES             Location/Qualifiers\n"
                "     CDS             complement(1..>12)\n"
                "                     /locus_tag=\"L1\"\n"
                "                     /translation=\"MPEP\"\n"
                "ORIGIN\n"
                "        1 atgcctgaaccc\n"
                "//\n",
                encoding="ascii",
            )
            rows = summarize_gbff(gbff, {"candidate": "MPEP"})
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["locus_tag"], "L1")
            self.assertEqual(rows[0]["strand"], -1)
            self.assertTrue(rows[0]["partial_location"])


if __name__ == "__main__":
    unittest.main()
