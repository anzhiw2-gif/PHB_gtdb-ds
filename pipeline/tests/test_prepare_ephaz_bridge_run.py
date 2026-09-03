import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_ephaz_bridge_run.py"


class PrepareBridgeRunTests(unittest.TestCase):
    def test_writes_isolated_auditable_bridge_run(self):
        spec = importlib.util.spec_from_file_location("prepare_ephaz_bridge_run", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.faa"
            source.write_text(">B2NHN2\nAAAA\n>O05527\nCCCC\n>P12625\nGGGG\n>Q51718\nTTTT\n", encoding="utf-8")
            control = root / "controls.faa"
            control.write_text(">A0ABY7N197\nAAAA\n>A0A3Q9BTL6\nCCCC\n>A0ABY9RN41\nGGGG\n>A0A3N5XWX8\nTTTT\n>A0ABT5L4W6\nAAAA\n>A0ABT7T1G3\nCCCC\n", encoding="utf-8")
            negative = root / "negative.faa"
            negative.write_text(">N1\nAAAA\n", encoding="utf-8")
            run = module.prepare(root, "20260830_test", source, control, negative)
            self.assertTrue((run / "input_contract.json").is_file())
            self.assertTrue((run / "results" / "ephaz_bridge_candidate.faa").is_file())
            self.assertIn("Q51871", (run / "results" / "ephaz_bridge_candidate.faa").read_text(encoding="utf-8"))
            self.assertTrue((run / "inputs" / "iPhaZ_like_challenge.faa").is_file())


if __name__ == "__main__":
    unittest.main()
