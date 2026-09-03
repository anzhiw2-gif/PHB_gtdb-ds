import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "ephaz_bridge_loo_calibration.py"


class BridgeLooCalibrationTests(unittest.TestCase):
    def _load(self):
        spec = importlib.util.spec_from_file_location("ephaz_bridge_loo_calibration", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_loo_calibration_requires_disjoint_holdout_and_preserves_challenge_roles(self):
        module = self._load()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            core = root / "core.faa"
            bridge = root / "bridge.faa"
            challenge = root / "challenge.faa"
            negative = root / "negative.faa"
            core.write_text(">C1\nAAAA\n>C2\nCCCC\n>C3\nGGGG\n>C4\nTTTT\n", encoding="utf-8")
            bridge.write_text(">B1\nACGT\n", encoding="utf-8")
            challenge.write_text(">X1\nACGT\n", encoding="utf-8")
            negative.write_text(">N1\nACGT\n", encoding="utf-8")

            calls = []
            def fake_run(command, **kwargs):
                calls.append(command)
                if command[0] == "mafft":
                    return type("R", (), {"stdout": ">aln\nACGT\n", "stderr": ""})()
                if command[0] == "hmmbuild":
                    Path(command[1]).write_text("HMMER3/f\n", encoding="utf-8")
                    return type("R", (), {"stdout": "", "stderr": ""})()
                if command[0] == "hmmsearch":
                    tbl = Path(command[command.index("--tblout") + 1])
                    dom = Path(command[command.index("--domtblout") + 1])
                    tbl.write_text("C4 - q - 1e-30 90 0 1e-30 90 0\n", encoding="utf-8")
                    dom.write_text("C4 - 4 q - 100 1e-30 90 0 1 1 1e-30 1e-30 90 0 1 80 1 80 1 80 0.9 d\n", encoding="utf-8")
                    return type("R", (), {"stdout": "", "stderr": ""})()
                raise AssertionError(command)

            with patch.object(module.subprocess, "run", side_effect=fake_run):
                result = module.run_loo_calibration(core, bridge, challenge, negative, root / "out", thresholds=(1e-5,))
            self.assertEqual(result["fold_count"], 4)
            rows = (root / "out" / "loo_summary.tsv").read_text(encoding="utf-8")
            self.assertIn("baseline_without_C4", rows)
            self.assertIn("bridge_augmented_without_C4", rows)
            self.assertTrue(any(call[0] == "hmmsearch" for call in calls))
            fold = root / "out" / "folds" / "baseline_without_C4"
            for name in ("training.faa", "alignment.faa", "model.hmm", "probe.faa", "hits.tblout", "hits.domtblout"):
                self.assertTrue((fold / name).is_file(), name)
            metadata = json.loads((root / "out" / "loo_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(len(metadata["models"]), 8)
            self.assertIn("sha256", metadata["models"][0]["artifacts"]["model_hmm"])
            self.assertEqual(set(metadata["software"]), {"mafft", "hmmbuild", "hmmsearch"})

    def test_main_binds_explicit_paths_and_tool_versions(self):
        module = self._load()
        with patch.object(module, "run_loo_calibration", return_value={"fold_count": 4}) as calibration:
            with patch.object(module, "tool_record", side_effect=lambda tool: {"path": tool, "version": f"{tool} version", "size": 1, "sha256": "a" * 64}):
                status = module.main([
                    "--core", "core.faa", "--bridge", "bridge.faa", "--challenge", "challenge.faa",
                    "--negative", "negative.faa", "--outdir", "out", "--mafft", "m", "--hmmbuild", "b", "--hmmsearch", "s",
                ])
        self.assertEqual(status, 0)
        self.assertEqual(calibration.call_args.kwargs["software_records"]["mafft"]["version"], "m version")
        self.assertEqual(calibration.call_args.kwargs["mafft_bin"], "m")

    def test_tool_version_uses_mafft_version_flag(self):
        module = self._load()
        with patch.object(module.subprocess, "run", return_value=type("R", (), {"returncode": 0, "stdout": "", "stderr": "MAFFT v7.520\n"})()) as run:
            self.assertEqual(module.tool_version("/opt/bin/mafft"), "MAFFT v7.520")
        self.assertEqual(run.call_args.args[0], ["/opt/bin/mafft", "--version"])

    def test_tool_record_captures_executable_hash(self):
        module = self._load()
        with tempfile.TemporaryDirectory() as tmp:
            executable = Path(tmp) / "hmmbuild"
            executable.write_text("tool", encoding="utf-8")
            with patch.object(module, "tool_version", return_value="HMMER 3.4"):
                record = module.tool_record(str(executable))
        self.assertEqual(record["path"], str(executable.resolve()))
        self.assertEqual(record["size"], 4)
        self.assertEqual(len(record["sha256"]), 64)

    def test_fasta_and_tblout_normalize_sp_and_tr_accessions(self):
        module = self._load()
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "probe.faa"
            fasta.write_text(">sp|P12345|PHAZ example\nAAAA\n>tr|A0A999|PHAZ example\nCCCC\n", encoding="utf-8")
            records = module.read_fasta(fasta)
            self.assertEqual(set(records), {"P12345", "A0A999"})
            tbl = Path(tmp) / "hits.tblout"
            tbl.write_text("sp|P12345|PHAZ - 1e-30 100 2e-30 90 0\ntr|A0A999|PHAZ - 1e-20 80 3e-20 70 0\n", encoding="utf-8")
            self.assertEqual(module.parse_hits(tbl), {"P12345": 2e-30, "A0A999": 3e-20})


if __name__ == "__main__":
    unittest.main()
