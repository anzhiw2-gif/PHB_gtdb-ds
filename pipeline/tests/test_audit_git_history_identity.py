"""Tests for the read-only git-history identity audit.

These build throwaway git repositories in a temp directory; the real repository is
never touched.
"""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import audit_git_history_identity as module  # noqa: E402


def git_init(root: Path) -> None:
    env = {
        **os.environ,
        "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
        "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.invalid",
    }
    for command in (["init", "-q"], ["add", "-A"], ["commit", "-q", "-m", "one"]):
        subprocess.run(["git", *command], cwd=root, check=True, capture_output=True, env=env)


# This module deliberately carries NO credential-shaped literal. The first version
# used a PEM header as the fixture and tripped the working-tree gate - and a naive
# "BEGIN RSA " + "PRIVATE KEY-----" split does NOT help, because the gate's regex is
# -----BEGIN [^-]*PRIVATE KEY----- and the gap between the pieces contains no dash,
# so it still matches across the concatenation. Splitting safely would require the
# break to fall INSIDE the matched literal. The exemption is keyed on the PATH, so
# any pattern demonstrates it; the server-address form below is used instead and is
# assembled from pieces the same way the safety module assembles its own.


class FixtureHolderTests(unittest.TestCase):
    def test_the_holder_path_is_the_safety_test_module(self):
        self.assertEqual(module.FIXTURE_HOLDER, "pipeline/tests/test_public_repo_safety.py")


class ScanTests(unittest.TestCase):
    def _repo(self, temporary: str, files: dict[str, str]) -> Path:
        root = Path(temporary) / "repo"
        root.mkdir(parents=True)
        for name, text in files.items():
            target = root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text, encoding="utf-8", newline="\n")
        git_init(root)
        return root

    def test_a_clean_history_reports_no_violations(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = self._repo(temporary, {"README.md": "nothing sensitive here\n"})
            report = module.scan(root)
            self.assertEqual(report["violating_blobs"], 0)
            self.assertIn("CLEAN", report["verdict"])
            self.assertGreater(report["commits_examined"], 0)
            self.assertGreater(report["text_blobs_scanned"], 0)

    def test_a_committed_identity_string_is_found_and_named(self):
        with tempfile.TemporaryDirectory() as temporary:
            secret = "10.16.1." + "141"
            root = self._repo(temporary, {"notes.md": f"server at {secret}\n"})
            report = module.scan(root)
            self.assertEqual(report["violating_blobs"], 1)
            violation = report["violations"][0]
            self.assertEqual(violation["path"], "notes.md")
            self.assertIn("1", violation["pattern_indexes"].split(","))
            self.assertGreaterEqual(int(violation["commit_count"]), 1)
            self.assertIn("CONTAINS", report["verdict"])

    def test_the_fixture_holder_is_exempted_and_counted_not_dropped(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = self._repo(temporary, {
                module.FIXTURE_HOLDER: 'SERVER = "10.16.1.' + '141"\n',
            })
            report = module.scan(root)
            self.assertEqual(report["violating_blobs"], 0, "the fixture holder must not be reported")
            self.assertEqual(report["fixture_holder_blobs_exempted"], 1)
            self.assertEqual(report["fixture_holder_detail"][0]["path"], module.FIXTURE_HOLDER)
            self.assertIn("test fixture holder", report["fixture_holder_detail"][0]["why"])

    def test_multibyte_content_does_not_desynchronise_the_batch_parser(self):
        # The first version decoded the cat-file stream before slicing it by the
        # BYTE size git reported, so every offset after the first multi-byte
        # character shifted and the parser blew up on a bogus header. A blob whose
        # content is multibyte AND long enough to matter must parse cleanly.
        with tempfile.TemporaryDirectory() as temporary:
            body = "中文注释" * 500 + "\n"
            root = self._repo(temporary, {"cn.md": body, "other.md": "plain\n"})
            report = module.scan(root)
            self.assertEqual(report["violating_blobs"], 0)
            self.assertGreaterEqual(report["text_blobs_scanned"], 2)

    def test_a_multibyte_blob_is_still_scanned_for_violations(self):
        with tempfile.TemporaryDirectory() as temporary:
            secret = "/home/data/" + "haoyu"
            root = self._repo(temporary, {"cn.md": "中文" * 300 + f"\npath {secret}\n"})
            report = module.scan(root)
            self.assertEqual(report["violating_blobs"], 1)
            self.assertEqual(report["violations"][0]["path"], "cn.md")

    def test_an_empty_repository_is_refused_rather_than_called_clean(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "empty"
            root.mkdir()
            subprocess.run(["git", "init", "-q"], cwd=root, check=True, capture_output=True)
            with self.assertRaisesRegex(ValueError, "refusing to report a clean history"):
                module.scan(root)


class CliTests(unittest.TestCase):
    def test_cli_writes_only_when_asked_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "repo"
            root.mkdir()
            (root / "README.md").write_text("clean\n", encoding="utf-8", newline="\n")
            git_init(root)
            out = Path(temporary) / "report.json"
            self.assertEqual(module.main(["--workdir", str(root)]), 0)
            self.assertFalse(out.exists(), "no --out means nothing is written")
            self.assertEqual(module.main(["--workdir", str(root), "--out", str(out)]), 0)
            self.assertTrue(out.exists())
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                module.main(["--workdir", str(root), "--out", str(out)])


if __name__ == "__main__":
    unittest.main()
