import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "04b_build_hmms_v2.sh"


class HMMBuilderInterfaceTests(unittest.TestCase):
    def test_builder_supports_explicit_layer_family_list(self):
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("--families", text)
        self.assertIn("--family-list", text)
        self.assertIn("ePhaZ_curated_core", text)
        self.assertIn("ePhaZ_broad_discovery", text)

    def test_builder_keeps_the_legacy_nine_family_default(self):
        text = SCRIPT.read_text(encoding="utf-8")
        for family in ("ePhaZ", "iPhaZ", "OH", "BdhA", "ArchPhaZ_patatin", "ArchPhaZ_hydrolase", "PhaJ", "PhaC", "phasin"):
            self.assertIn(family, text)


if __name__ == "__main__":
    unittest.main()
