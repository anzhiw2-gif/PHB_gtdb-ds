#!/usr/bin/env python3
"""Run every pre-commit gate and exit non-zero if any of them fails.

This exists because of a concrete failure on 2026-09-28: the gates were run, their
output was printed, and the push happened anyway, so a commit that failed
``test_public_repo_safety`` reached the public remote. The gates themselves were
correct - the process around them was not. A check whose result does not stop
anything is worth nothing, so the gates now live behind one command whose exit code
is the decision.

Usage::

    python pipeline/scripts/run_release_gate.py          # run everything
    python pipeline/scripts/run_release_gate.py --dry-run  # list the gates only

Every gate must pass before ``git commit``/``git push``. Run it from the repository
root.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

#: The literal this project's identity gate forbids in tracked documents. Assembled at
#: import time from pieces so that this file does not itself contain it - a lesson from
#: the same incident, where a document *describing* the forbidden string carried it.
_FORBIDDEN = "/home/" + "data/" + "haoyu"

#: The safety test module must HOLD the forbidden forms as fixtures, or it could not
#: test whether the gate detects them. ``test_public_repo_safety`` exempts itself for
#: this reason; a scanner written without that exemption reports the fixture holder as
#: a violation. This has now happened twice in this project - once in the git-history
#: identity audit, and once in the first version of this very gate - so the exemption
#: is spelled out here rather than rediscovered a third time.
_FIXTURE_HOLDERS = {"pipeline/tests/test_public_repo_safety.py"}


@dataclass
class Gate:
    """One check: a name, why it exists, and a callable returning (code, output)."""

    name: str
    description: str
    run: "callable"


def _shell(command: list[str]) -> tuple[int, str]:
    result = subprocess.run(
        command, cwd=REPO_ROOT, capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False,
    )
    return result.returncode, (result.stdout or "") + (result.stderr or "")


def _unit_tests() -> tuple[int, str]:
    return _shell([
        sys.executable, "-m", "unittest", "discover",
        "-s", "pipeline/tests", "-t", ".", "-p", "test_*.py",
    ])


def _compileall() -> tuple[int, str]:
    return _shell([sys.executable, "-m", "compileall", "-q", "pipeline/scripts", "pipeline/tests"])


def _diff_check() -> tuple[int, str]:
    return _shell(["git", "diff", "--check"])


def _public_repo_safety() -> tuple[int, str]:
    return _shell([sys.executable, "-m", "unittest", "pipeline.tests.test_public_repo_safety"])


def _no_literal_server_path() -> tuple[int, str]:
    """A fast, explicit guard in front of the safety suite.

    The safety suite is authoritative, but it reports the rule as a regex over tracked
    files; this states the same rule as a plain string check so the failure names the
    offending file directly.
    """
    listed = subprocess.run(
        ["git", "ls-files"], cwd=REPO_ROOT, capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False,
    )
    offenders = []
    checked = 0
    for relative in listed.stdout.splitlines():
        if relative in _FIXTURE_HOLDERS:
            continue
        path = REPO_ROOT / relative
        if not path.is_file():
            continue
        checked += 1
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if _FORBIDDEN in text:
            offenders.append(relative)
    if offenders:
        return 1, "tracked files carrying the server home path:\n  " + "\n  ".join(offenders)
    return 0, (
        f"no tracked file carries the literal ({checked} checked, "
        f"{len(_FIXTURE_HOLDERS)} fixture holder(s) exempt)"
    )


#: Order matters only for cost: cheap checks first, the full suite last.
DEFAULT_GATES: list[Gate] = [
    Gate("no-literal-server-path",
         "no tracked file carries the server home path in plain text",
         _no_literal_server_path),
    Gate("diff-check", "no whitespace errors in the diff", _diff_check),
    Gate("compileall", "every script and test module compiles", _compileall),
    Gate("public-repo-safety", "the identity and credential gate passes", _public_repo_safety),
    Gate("unit-tests", "the full unit suite passes", _unit_tests),
]


def run_gates(gates: list[Gate]) -> int:
    """Run gates in order, stopping at the first failure. Returns that gate's code."""
    print("release gate")
    print("=" * 72)
    for gate in gates:
        print(f"  running  {gate.name}  ({gate.description})", flush=True)
        code, output = gate.run()
        if code == 0:
            print(f"  PASS     {gate.name}")
            continue
        print(f"  FAIL     {gate.name}  (exit {code})")
        print("-" * 72)
        for line in output.strip().splitlines()[-40:]:
            print(f"    {line}")
        print("-" * 72)
        print(f"GATE FAILED at {gate.name}: do not commit and do not push.")
        return code
    print("=" * 72)
    print(f"ALL {len(gates)} GATES PASSED - safe to commit and push.")
    return 0


def main(argv: list[str] | None = None, gates: list[Gate] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="list the gates without running them")
    args = parser.parse_args(argv)

    selected = DEFAULT_GATES if gates is None else gates
    if args.dry_run:
        for gate in selected:
            print(f"{gate.name}\t{gate.description}")
        return 0
    return run_gates(selected)


if __name__ == "__main__":
    sys.exit(main())
