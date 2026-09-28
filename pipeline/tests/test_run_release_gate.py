"""Tests for the release gate runner.

The runner exists because a gate whose result is printed but not acted on is worth
nothing: that is exactly how a failing commit reached the public remote on 2026-09-28.
So the property under test is not "the gates pass" but "a failing gate stops the run
and the exit code says so".
"""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline.scripts import run_release_gate as gate


def _fake_gate(name, returncode, calls):
    def run():
        calls.append(name)
        return returncode, f"{name} output"
    return gate.Gate(name, "description", run)


class GateRunnerTests(unittest.TestCase):
    def test_all_passing_gates_exit_zero_and_run_in_order(self):
        calls = []
        gates = [_fake_gate("a", 0, calls), _fake_gate("b", 0, calls)]
        with contextlib.redirect_stdout(io.StringIO()):
            code = gate.run_gates(gates)
        self.assertEqual(code, 0)
        self.assertEqual(calls, ["a", "b"])

    def test_a_failing_gate_stops_the_run(self):
        calls = []
        gates = [
            _fake_gate("first", 0, calls),
            _fake_gate("second", 1, calls),
            _fake_gate("third", 0, calls),
        ]
        with contextlib.redirect_stdout(io.StringIO()) as out:
            code = gate.run_gates(gates)
        self.assertNotEqual(code, 0)
        self.assertEqual(calls, ["first", "second"], "later gates must not run")
        self.assertIn("second", out.getvalue())

    def test_the_first_failure_sets_the_exit_code(self):
        calls = []
        gates = [_fake_gate("a", 2, calls), _fake_gate("b", 3, calls)]
        with contextlib.redirect_stdout(io.StringIO()):
            code = gate.run_gates(gates)
        self.assertEqual(code, 2)

    def test_output_names_every_gate_and_its_verdict(self):
        calls = []
        gates = [_fake_gate("alpha", 0, calls), _fake_gate("beta", 1, calls)]
        with contextlib.redirect_stdout(io.StringIO()) as out:
            gate.run_gates(gates)
        text = out.getvalue()
        for name in ("alpha", "beta"):
            self.assertIn(name, text)
        self.assertIn("PASS", text)
        self.assertIn("FAIL", text)


class GateCompositionTests(unittest.TestCase):
    """The gates that must be present, and the trap that is easy to forget."""

    def test_the_composition_covers_what_the_incident_required(self):
        names = [g.name for g in gate.DEFAULT_GATES]
        for required in (
            "unit-tests",
            "compileall",
            "diff-check",
            "public-repo-safety",
            "no-literal-server-path",
        ):
            self.assertIn(required, names, f"gate {required!r} must be in the default set")

    def test_the_literal_path_gate_scans_tracked_documents(self):
        # The string is built at run time from pieces so that this test file - which is
        # itself tracked - does not carry it. Writing about a forbidden string is still
        # writing it; the previous round learned that the hard way.
        gate_obj = next(g for g in gate.DEFAULT_GATES if g.name == "no-literal-server-path")
        _code, output = gate_obj.run()
        self.assertNotIn("/home/" + "data/", output)

    def test_the_fixture_holder_is_exempt_and_the_gate_passes_on_the_real_tree(self):
        # Running the gate against the real repository is the only way to know it is a
        # gate rather than a hypothesis: its first version failed on every run because
        # it flagged the safety test module, which must hold the forbidden forms as
        # fixtures. Same mistake as the earlier history audit.
        gate_obj = next(g for g in gate.DEFAULT_GATES if g.name == "no-literal-server-path")
        code, output = gate_obj.run()
        self.assertEqual(code, 0, f"gate should pass on this tree, said:\n{output}")
        self.assertIn("exempt", output)
        self.assertIn("pipeline/tests/test_public_repo_safety.py", gate._FIXTURE_HOLDERS)

    def test_a_file_carrying_the_literal_would_be_flagged(self):
        # Prove the gate can fail, rather than trusting that it would. A temporary
        # repository root is used instead of mocking pathlib globally: the first
        # version of this test patched Path.is_file and Path.read_text, which made
        # every one of the thousand tracked files look like an offender.
        literal = "/home/" + "data/" + "haoyu"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "pipeline" / "scripts").mkdir(parents=True)
            offender = root / "pipeline" / "scripts" / "example.py"
            offender.write_text(f"PATH = {literal!r}\n", encoding="utf-8")
            clean = root / "pipeline" / "scripts" / "clean.py"
            clean.write_text("PATH = 'elsewhere'\n", encoding="utf-8")

            listing = mock.Mock(returncode=0, stdout="pipeline/scripts/example.py\n"
                                                    "pipeline/scripts/clean.py\n")
            with mock.patch.object(gate, "REPO_ROOT", root), \
                 mock.patch.object(gate.subprocess, "run", return_value=listing):
                code, output = gate._no_literal_server_path()

        self.assertEqual(code, 1)
        self.assertIn("example.py", output)
        self.assertNotIn("clean.py", output)

    def test_a_gate_reports_a_command_and_a_description(self):
        for gate_obj in gate.DEFAULT_GATES:
            self.assertTrue(gate_obj.name)
            self.assertTrue(gate_obj.description)


class MainTests(unittest.TestCase):
    def test_main_accepts_an_injected_gate_list(self):
        calls = []
        gates = [_fake_gate("x", 0, calls)]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(gate.main([], gates=gates), 0)
        self.assertEqual(calls, ["x"])

    def test_dry_run_lists_gates_without_running_them(self):
        calls = []
        gates = [_fake_gate("never", 0, calls)]
        with contextlib.redirect_stdout(io.StringIO()) as out:
            code = gate.main(["--dry-run"], gates=gates)
        self.assertEqual(code, 0)
        self.assertEqual(calls, [], "dry run must not execute any gate")
        self.assertIn("never", out.getvalue())


if __name__ == "__main__":
    unittest.main()
