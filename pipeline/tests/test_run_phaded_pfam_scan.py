import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_phaded_pfam_scan.py"


class PhaDEDPfamScanTests(unittest.TestCase):
    def test_rejects_cpu_above_server_limit_before_running_tool(self):
        spec = importlib.util.spec_from_file_location("run_phaded_pfam_scan", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            database = root / "database.hmm"; database.write_text("x", encoding="ascii")
            fasta = root / "input.faa"; fasta.write_text(">a\nM\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "between 1 and 40"):
                module.run_scan("not-run", database, fasta, root / "out.domtblout", 41)

    def test_scan_command_suppresses_alignment_stdout(self):
        spec = importlib.util.spec_from_file_location("run_phaded_pfam_scan", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertIn("--noali", module.scan_command("hmmscan", Path("db.hmm"), Path("input.faa"), Path("out.dom"), 40))
