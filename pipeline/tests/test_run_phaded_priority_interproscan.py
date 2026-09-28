import hashlib
import importlib.util
import unittest
from unittest.mock import patch

from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_phaded_priority_interproscan.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_phaded_priority_interproscan", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PriorityInterProScanTests(unittest.TestCase):
    def test_run_resolves_output_and_temp_paths_before_launch(self):
        module = load_module()
        root = Path(self.id().replace('.', '_'))
        fasta = root / 'input.faa'
        output = root / 'nested' / 'out.tsv'
        temp = root / 'nested' / 'tmp'
        try:
            fasta.parent.mkdir(parents=True, exist_ok=True)
            fasta.write_text('>x\nMSEQ\n', encoding='ascii')
            output.parent.mkdir(parents=True, exist_ok=True)
            def fake_run(command, **kwargs):
                Path(command[4]).write_text('result\n', encoding='ascii')
                return type('Completed', (), {'stdout': '', 'stderr': ''})()
            with patch.object(module.subprocess, 'run', side_effect=fake_run) as mocked:
                result = module.run('interproscan.sh', fasta, output, temp, 1, hashlib.sha256(fasta.read_bytes()).hexdigest())
            command = mocked.call_args.args[0]
            self.assertTrue(Path(command[2]).is_absolute())
            self.assertTrue(Path(command[4]).is_absolute())
            self.assertTrue(Path(command[11]).is_absolute())
            self.assertEqual(result['output']['path'], str(output.resolve()))
        finally:
            import shutil
            shutil.rmtree(root, ignore_errors=True)

    def test_command_is_tsv_iprlookup_and_capped_at_40_cpus(self):
        module = load_module()
        command = module.interpro_command("interproscan.sh", Path("input.faa"), Path("out.tsv"), Path("tmp"), 40)
        self.assertEqual(command, ["interproscan.sh", "-i", str(Path("input.faa").resolve()), "-o", str(Path("out.tsv").resolve()), "-f", "TSV", "--iprlookup", "-cpu", "40", "-T", str(Path("tmp").resolve())])
        with self.assertRaisesRegex(ValueError, "between 1 and 40"):
            module.interpro_command("interproscan.sh", Path("input.faa"), Path("out.tsv"), Path("tmp"), 41)

    def test_input_verification_rejects_hash_mismatch(self):
        module = load_module()
        path = Path(self.id().replace(".", "_"))
        try:
            path.write_text("MSEQ\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                module.verify_input(path, "0" * 64)
            self.assertEqual(module.verify_input(path, hashlib.sha256(path.read_bytes()).hexdigest()), path.stat().st_size)
        finally:
            path.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
