import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from pipeline.scripts.run_phaded_full_signalp import build_signalp_command


class FullSignalPCommandTests(unittest.TestCase):
    def test_builds_pinned_signalp6_command_with_at_most_forty_threads(self):
        with TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "candidates.faa"
            fasta.write_text(">x\nM\n", encoding="ascii")
            output = Path(tmp) / "output"
            command = build_signalp_command(
                fasta,
                output,
                torch_threads=32,
                write_processes=8,
            )
        self.assertIn("signalp6", command)
        self.assertIn("--fastafile", command)
        self.assertIn(str(fasta), command)
        self.assertIn("--output_dir", command)
        self.assertIn(str(output), command)
        self.assertIn("--torch_num_threads", command)
        self.assertIn("32", command)
        self.assertIn("--write_procs", command)

    def test_rejects_thread_count_above_project_limit(self):
        with TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "candidates.faa"
            fasta.write_text(">x\nM\n", encoding="ascii")
            with self.assertRaises(ValueError):
                build_signalp_command(
                    fasta,
                    Path(tmp) / "output",
                    torch_threads=41,
                    write_processes=8,
                )
