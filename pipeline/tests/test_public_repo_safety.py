"""Public-repository safety checks for tracked text and configuration files."""
import os
import re
import subprocess
import unittest


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

SELF_RELPATH = "pipeline/tests/test_public_repo_safety.py"

# ---------------------------------------------------------------------------
# Local-machine identity patterns.
#
# IMPORTANT (Windows path escaping): every Windows-path pattern below uses
# ``[\\/]{1,2}`` -- a regex character class matching ONE backslash, TWO
# backslashes or ONE forward slash. Writing ``re.escape("C:" + "\\\\Users")``
# instead would build the two-character sequence ``\\`` (escaped backslash),
# which matches only a *doubled* backslash and therefore silently misses every
# real single-backslash path such as ``D:\PHB_gtdb-ds``. Regression tests for
# this live in :class:`IdentityPatternRegressionTests`.
#
# The same rationale applies to the negative lookbehinds: they must tolerate
# both the single- and doubled-backslash spellings of an already-sanitized
# placeholder (``${PHB_REMOTE_ROOT}`` / ``\<REPO_ROOT\>``) so sanitized text is
# not re-flagged.
# ---------------------------------------------------------------------------
IDENTITY_PATTERNS = (
    re.escape("10.16.1." + "141"),
    r"(?i)(?<!<SERVER_USER>)" + re.escape("hao" + "yu@"),
    r"(?i)(?<![\\/])(?<!\$\{PHB_REMOTE_ROOT\})" + re.escape("/home/data/" + "haoyu")
    + r"(?![0-9A-Za-z_])",
    r"(?i)" + re.escape("C:") + r"[\\/]{1,2}" + re.escape("Users") + r"[\\/]{1,2}"
    + re.escape("HUAWEI") + r"(?![0-9A-Za-z_])",
    r"(?i)" + re.escape("D:") + r"[\\/]{1,2}" + re.escape("PHB_gtdb-ds"),
    r"-----BEGIN [^-]*PRIVATE KEY-----",
    r"\b(?:ghp|github_pat)_[A-Za-z0-9_]+\b",
)


def identity_patterns():
    """Return the public-repository identity/credential pattern tuple."""
    return IDENTITY_PATTERNS


def identity_violations(text):
    """Return the pattern strings that ``text`` violates (empty when clean)."""
    return [pattern for pattern in IDENTITY_PATTERNS if re.search(pattern, text)]


class PublicRepoSafetyTests(unittest.TestCase):
    def test_tracked_files_do_not_expose_machine_identity_or_credentials(self):
        tracked = subprocess.run(
            ["git", "-c", "core.quotePath=false", "ls-files", "-z"],
            cwd=ROOT,
            check=True,
            capture_output=True,
        ).stdout.decode().split("\0")
        patterns = identity_patterns()
        violations = []
        for rel in tracked:
            if rel.replace("\\", "/") == SELF_RELPATH:
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


class IdentityPatternRegressionTests(unittest.TestCase):
    """Pins the Windows-path escaping fix and the no-over-matching guarantee."""

    FLAGGED = (
        # single backslash -- the spelling that the pre-fix pattern missed
        r"D:\PHB_gtdb-ds",
        r"D:\PHB_gtdb-ds\pipeline\scripts\06_screen.sh",
        r"C:\Users\HUAWEI",
        r"C:\Users\HUAWEI\AppData\Local\Temp",
        # doubled backslash
        r"D:\\PHB_gtdb-ds",
        r"C:\\Users\\HUAWEI",
        # forward slash
        "D:/PHB_gtdb-ds",
        "C:/Users/HUAWEI",
        # other identity patterns of the tuple
        "haoyu@10.16.1.141",
        "logger@10.16.1.141",
        "ssh haoyu@10.16.1.141",
        "/home/data/haoyu/PHB_gtdb-ds",
        # boundary tightening must not lose a real hit at these spellings
        "/home/data/haoyu`",
        "/home/data/haoyu)",
        r"C:\Users\HUAWEI\AppData",
        r"C:\Users\HUAWEI`",
        "-----BEGIN RSA PRIVATE KEY-----",
        "-----BEGIN OPENSSH PRIVATE KEY-----",
        "token ghp_" + "A" * 36,
        "token github_pat_" + "b" * 40,
    )

    NOT_FLAGGED = (
        "<REPO_ROOT>",
        "<LOCAL_HOME>",
        "<SERVER_USER>",
        "<SERVER_HOST>",
        "${PHB_REMOTE_ROOT}",
        "${PHB_REMOTE_ROOT}/PHB_gtdb-ds",
        "<SERVER_USER>@<SERVER_HOST>",
        r"<REPO_ROOT>\runs\<run_id>",
        # no over-matching on unrelated drives/dirs/users
        r"D:\projects\other",
        r"C:\Users\PUBLIC",
        r"C:\Users\HUAWEI2",
        "/home/data/haoyu2/other",
        "10.16.1.142",
        "haoyu2@example.org",
        "user haoyu-other",
    )

    def test_single_backslash_windows_paths_are_flagged(self):
        for sample in (r"D:\PHB_gtdb-ds", r"D:\PHB_gtdb-ds\runs",
                       r"D:\PHB_gtdb-ds\pipeline", r"C:\Users\HUAWEI"):
            with self.subTest(sample=sample):
                self.assertNotEqual([], identity_violations(sample))

    def test_doubled_and_forward_slash_spellings_are_flagged(self):
        for sample in (r"D:\\PHB_gtdb-ds", r"C:\\Users\\HUAWEI",
                       "D:/PHB_gtdb-ds", "C:/Users/HUAWEI"):
            with self.subTest(sample=sample):
                self.assertNotEqual([], identity_violations(sample))

    def test_all_forbidden_samples_are_flagged(self):
        for sample in self.FLAGGED:
            with self.subTest(sample=sample):
                self.assertNotEqual(
                    [], identity_violations(sample), f"not flagged: {sample!r}")

    def test_placeholders_and_benign_text_are_not_flagged(self):
        for sample in self.NOT_FLAGGED:
            with self.subTest(sample=sample):
                self.assertEqual(
                    [], identity_violations(sample), f"over-matched: {sample!r}")

    def test_forward_slash_remote_root_placeholder_is_not_flagged(self):
        # `${PHB_REMOTE_ROOT}/...` must survive; a bare /home/data/haoyu must not.
        self.assertEqual([], identity_violations("${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs"))
        self.assertNotEqual([], identity_violations("/home/data/haoyu/PHB_gtdb-ds/runs"))

    def test_predicate_ignores_already_sanitized_placeholder_forms(self):
        sanitized = (
            "工作区：`<REPO_ROOT>`（分支 `main`）\n"
            "服务器 ${PHB_REMOTE_ROOT}/PHB_gtdb-ds，用户 <SERVER_USER>@<SERVER_HOST>\n"
        )
        self.assertEqual([], identity_violations(sanitized))


if __name__ == "__main__":
    unittest.main()
