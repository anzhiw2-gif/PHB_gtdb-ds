import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_phaded_mapping_inputs.py"


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_phaded_mapping_inputs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PreparePhaDEDInputsTests(unittest.TestCase):
    def test_identical_cross_tier_duplicate_is_allowed(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "tier.faa"
            path.write_text(">candidate\nMABC\n", encoding="ascii")
            records = {}
            module.read_fasta(path, records)
            module.read_fasta(path, records)
            self.assertEqual(records, {"candidate": "MABC"})

    def test_conflicting_cross_tier_duplicate_fails_closed(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            first = Path(temporary_directory) / "first.faa"
            second = Path(temporary_directory) / "second.faa"
            first.write_text(">candidate\nMABC\n", encoding="ascii")
            second.write_text(">candidate\nMABD\n", encoding="ascii")
            records = {}
            module.read_fasta(first, records)
            with self.assertRaisesRegex(ValueError, "conflicting FASTA sequence"):
                module.read_fasta(second, records)


if __name__ == "__main__":
    unittest.main()
