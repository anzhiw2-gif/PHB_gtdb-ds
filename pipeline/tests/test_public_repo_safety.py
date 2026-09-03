"""Public-repository safety checks for tracked text and configuration files."""
import os
import re
import subprocess
import unittest


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


class PublicRepoSafetyTests(unittest.TestCase):
    def test_tracked_files_do_not_expose_machine_identity_or_credentials(self):
        tracked = subprocess.run(
            ["git", "ls-files", "-z"],
            cwd=ROOT,
            check=True,
            capture_output=True,
        ).stdout.decode().split("\0")
        patterns = (
            re.escape("10.16.1." + "141"),
            r"(?i)(?<!<SERVER_USER>)" + re.escape("hao" + "yu@"),
            r"(?i)(?<!\$\{PHB_REMOTE_ROOT\})" + re.escape("/home/data/" + "haoyu"),
            r"(?i)" + re.escape("C:" + "\\\\Users\\\\HUAWEI"),
            r"(?i)" + re.escape("D:" + "\\\\PHB_gtdb-ds"),
            r"-----BEGIN [^-]*PRIVATE KEY-----",
            r"\b(?:ghp|github_pat)_[A-Za-z0-9_]+\b",
        )
        violations = []
        for rel in tracked:
            if rel.replace("\\", "/") == "pipeline/tests/test_public_repo_safety.py":
                continue
            if not rel or os.path.isdir(os.path.join(ROOT, rel)):
                continue
            try:
                with open(os.path.join(ROOT, rel), "rb") as handle:
                    data = handle.read()
            except OSError:
                continue
            if b"\0" in data:
                continue
            text = data.decode("utf-8", errors="ignore")
            for pattern in patterns:
                if re.search(pattern, text):
                    violations.append(f"{rel}: {pattern}")
        self.assertEqual([], violations, "public safety violations:\n" + "\n".join(violations))


if __name__ == "__main__":
    unittest.main()
