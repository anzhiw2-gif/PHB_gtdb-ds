"""Task 9 — one full-database ``-Z`` for every HMMER shard search.

Requirement (docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md,
"HMMER 分片尺度"): 正式入口必须记录完整目标序列数 `Z`; if shard search is used,
every shard must be searched with the *same* full-database ``-Z``, and the
manifest must store ``database_size_Z`` plus the per-shard sequence counts and
the command.  A shard E-value computed against HMMER's default ``Z`` (the
shard's own sequence count) is a *per-shard* scale and must never be reported as
a single full-library scale.

These tests are the failing test for Task 9.  They exercise:

* ``pipeline/scripts/hmmer_command.py`` — the dependency-free command builder
  shared with ``06_screen.sh``;
* the shell entrypoint text — proving one ``-Z`` per per-shard HMMER call;
* ``06_validate_screen_manifest.py`` — the manifest scale contract.

Task F4 extends the same contract to the *tier* path: the tier entrypoint
(``formal_scan13_tier_processing.sh``) and its rescoring step
(``08c_tier_rescore.py``) re-run HMMER, so they must also require one
full-library ``-Z``, pass it to every call, and record it in
``tier_processing_manifest.json``.  Missing scale evidence must make the
manifest invalid instead of being defaulted to a shard-local count.
"""

from __future__ import annotations

import ast
import importlib.util
import io
import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
MODULE_PATH = SCRIPTS / "hmmer_command.py"
SCREEN_SH = SCRIPTS / "06_screen.sh"
VALIDATOR = SCRIPTS / "06_validate_screen_manifest.py"

# Task F4 — the tier entrypoint and the rescore step it re-runs HMMER through.
TIER_SH = SCRIPTS / "formal_scan13_tier_processing.sh"
TIER_RESCORE = SCRIPTS / "08c_tier_rescore.py"

FULL_DATABASE_Z = 292_000_000


def load_hmmer_command():
    """Import ``pipeline/scripts/hmmer_command.py`` by path.

    The module lives under ``pipeline/scripts`` (which is not a package), so it
    is loaded by file path.  This mirrors how the shell entrypoint and the
    validator treat it: a self-contained file with no third-party imports.
    """
    if not MODULE_PATH.is_file():
        raise AssertionError(f"missing command builder: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location("hmmer_command", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("hmmer_command", module)
    spec.loader.exec_module(module)
    return module


def load_validator():
    """Import ``pipeline/scripts/06_validate_screen_manifest.py`` by path."""
    spec = importlib.util.spec_from_file_location(
        "screen_manifest_validator_task9", VALIDATOR
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_one_call_body(script_text: str) -> str:
    """Return the body of the per-shard ``hmmsearch`` worker function."""
    match = re.search(r"^run_one\(\)\s*\{(.*?)^\}\s*$", script_text, re.M | re.S)
    if match is None:
        match = re.search(
            r"^\s*(?:hmmsearch_shard|screen_shard)\(\)\s*\{(.*?)^\}\s*$",
            script_text,
            re.M | re.S,
        )
    if match is None:
        raise AssertionError("06_screen.sh has no per-shard HMMER worker function")
    return match.group(1)


def per_shard_loop_body(script_text: str) -> str:
    """Return the per-family/per-shard dispatch loop (from the first shard loop).

    The region deliberately starts at the first ``run_one`` definition, so an
    unrelated ``-Z`` mention elsewhere in the script (for example the run
    contract heredoc) cannot be mistaken for a HMMER flag.
    """
    start = script_text.index("run_one() {")
    end = script_text.index("parallel -j")
    return script_text[start:end]


def strip_shell_comments(text: str) -> str:
    """Drop whole-line shell comments so prose cannot satisfy a flag assertion."""
    return "\n".join(
        line for line in text.splitlines() if not line.lstrip().startswith("#")
    )


#: Shell source with comments removed: the counterpart of
#: :func:`python_executable_text` for the entrypoint scripts.
shell_code = strip_shell_comments


def load_rescore_module():
    """Import the tier rescoring step by path (it lives outside any package)."""
    if not TIER_RESCORE.is_file():
        raise AssertionError(f"missing tier rescoring script: {TIER_RESCORE}")
    spec = importlib.util.spec_from_file_location("tier_rescore_task_f4", TIER_RESCORE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tier_heredoc_blocks() -> list[str]:
    """Return every embedded Python program of the tier entrypoint."""
    return re.findall(r"<<'PY'\n(.*?)\nPY\n", TIER_SH.read_text(encoding="utf-8"), re.S)


def tier_heredoc_block(*markers: str) -> str:
    """Return the one embedded program containing all ``markers``."""
    matches = [
        block for block in tier_heredoc_blocks() if all(m in block for m in markers)
    ]
    if len(matches) != 1:
        raise AssertionError(
            f"expected exactly one tier-driver heredoc matching {markers}, "
            f"found {len(matches)}"
        )
    return matches[0]


def load_tier_manifest_programs(namespace=None):
    """Execute the entrypoint's embedded scale code and return its namespace.

    The database-size checks that make ``tier_processing_manifest.json`` valid
    live inside the driver as embedded Python.  Executing that code directly is
    the only way to test the *real* artefacts — a look-alike copy in the test
    would pass while the shipped driver stayed unscaled.  Only the
    manifest-writing program is executed (the others need real arguments).
    """
    namespace = {} if namespace is None else dict(namespace)
    block = tier_heredoc_block("tier_processing_manifest.json")
    code = compile(block, "<tier_driver_manifest_heredoc>", "exec")
    exec(code, namespace)
    return namespace


def dict_keys_from_node(node):
    """Literal string keys of a dict literal (nested dicts are not descended)."""
    if not isinstance(node, ast.Dict):
        raise AssertionError("expected a dict literal")
    return {
        key.value
        for key in node.keys
        if isinstance(key, ast.Constant) and isinstance(key.value, str)
    }


def flag_arguments(source: str) -> set[str]:
    """Every CLI flag the source passes as a *literal argument* string.

    Docstrings and comments are removed first, so prose cannot satisfy a flag
    assertion — but prose cannot break one either, because the word boundary
    makes ``--domZ`` a distinct flag from ``--domtblout``.
    """
    code = python_code_text(source)
    return set(re.findall(r"(?<![\w-])--[A-Za-z][\w-]*", code))


def manifest_payload_keys(block: str) -> set[str]:
    """Keys of the manifest dict the driver's builder returns.

    The payload is read out of the embedded program's own
    ``build_tier_processing_manifest`` function, so the assertion tracks what
    the driver really writes rather than a re-typed look-alike.
    """
    tree = ast.parse(block)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "build_tier_processing_manifest":
            for statement in node.body:
                if isinstance(statement, ast.Assign) and isinstance(statement.value, ast.Dict):
                    return dict_keys_from_node(statement.value)
            raise AssertionError("no manifest dict literal in the tier builder")
    raise AssertionError("the tier driver has no manifest builder function")


def called_name(node) -> str | None:
    """The function name a ``Call`` node invokes (plain name or attribute)."""
    func = getattr(node, "func", None)
    if isinstance(func, ast.Attribute):
        return func.attr
    if isinstance(func, ast.Name):
        return func.id
    return None


def _unescape_source(text: str) -> str:
    """Undo ``ast.unparse``'s source-level escapes."""

    def replace(match):
        body = match.group(1)
        simple = {"n": "\n", "t": "\t", "\\": "\\", "'": "'", '"': '"', "\n": ""}
        if body in simple:
            return simple[body]
        return chr(int(body[2:], 16)) if body.startswith(("u", "U", "x")) else match.group(0)

    return re.sub(r"\\(u[0-9a-fA-F]{4}|U[0-9a-fA-F]{8}|x[0-9a-fA-F]{2}|.)", replace, text)


def python_code_text(source: str) -> str:
    """Python source with comments and docstrings removed.

    Prose may legitimately explain that ``--domZ`` is never emitted; only real
    code and real flag strings may satisfy (or break) an assertion.
    """
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(
            node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ):
            body = getattr(node, "body", [])
            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)
            ):
                body.pop(0)
    return _unescape_source(ast.unparse(tree))


def environ_names_read_in_main(block: str) -> set[str]:
    """The environment lookups the rescore step performs anywhere in its source.

    This is the dynamic counterpart of the module-constant scan: the full-library
    ``Z`` may come from the command line or from the documented
    ``PHB_DATABASE_SIZE_Z`` override, never from a shard-local variable.  Reading
    *any* other variable would be an undocumented second source of scale, so the
    result must be exactly the two documented lookups.
    """
    names: set[str] = set()

    def record(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            names.add(node.value)
        elif isinstance(node, ast.Name):
            names.add(node.id)

    for child in ast.walk(ast.parse(block)):
        # ``os.environ["X"]`` and ``os.environ.get("X")`` both count.
        if isinstance(child, ast.Subscript):
            target, index = child.value, child.slice
            if isinstance(target, ast.Attribute) and target.attr == "environ":
                record(index)
        elif isinstance(child, ast.Call):
            func = child.func
            if (
                isinstance(func, ast.Attribute)
                and func.attr == "get"
                and isinstance(func.value, ast.Attribute)
                and func.value.attr == "environ"
                and child.args
            ):
                record(child.args[0])
    return names


class HmmerCommandBuilderTests(unittest.TestCase):
    """The shared builder must always pin the full-database ``-Z``."""

    def setUp(self):
        self.module = load_hmmer_command()

    def test_every_shard_uses_the_same_full_database_z(self):
        # The value here is deliberately arbitrary: this test is about pass-through,
        # not about any particular library size. It previously used 292000000, which
        # was the long-standing *approximation* for scan-13 and is now known to be
        # about 2.1x too low - using a refuted constant as an example invites copying
        # it. The measured scan-13 full-library Z is 615,969,589
        # (runs/20260928_phaded_scan13_z_scale_reconciliation_01/results/
        # shard_record_counts.tsv).
        command = self.module.build_hmmsearch_command(
            "model.hmm", "shard.faa", total_targets=123456789
        )
        self.assertIn("-Z", command)
        self.assertEqual(command[command.index("-Z") + 1], "123456789")

    def test_z_is_the_same_for_two_different_shards_of_one_library(self):
        first = self.module.build_hmmsearch_command(
            "model.hmm", "shard_0001.faa", total_targets=FULL_DATABASE_Z
        )
        second = self.module.build_hmmsearch_command(
            "model.hmm", "shard_0100.faa", total_targets=FULL_DATABASE_Z
        )
        self.assertEqual(
            first[first.index("-Z") + 1], second[second.index("-Z") + 1]
        )
        self.assertEqual(first[first.index("-Z") + 1], str(FULL_DATABASE_Z))

    def test_z_is_not_the_shard_sequence_count(self):
        """A shard-local default Z must never leak into the shared scale."""
        command = self.module.build_hmmsearch_command(
            "model.hmm", "shard_0051.faa", total_targets=FULL_DATABASE_Z
        )
        self.assertNotEqual(command[command.index("-Z") + 1], "2923820")

    def test_negative_zero_and_empty_targets_are_rejected(self):
        for bad in (0, -1, -292_000_000, ""):
            with self.subTest(total_targets=bad):
                with self.assertRaises(ValueError):
                    self.module.build_hmmsearch_command(
                        "model.hmm", "shard.faa", total_targets=bad
                    )

    def test_non_integer_targets_are_rejected(self):
        for bad in ("292000000.0", 292000000.5, None, [292000000], True):
            with self.subTest(total_targets=bad):
                with self.assertRaises(ValueError):
                    self.module.build_hmmsearch_command(
                        "model.hmm", "shard.faa", total_targets=bad
                    )

    def test_missing_or_empty_paths_are_rejected(self):
        for hmm, shard in (("", "shard.faa"), ("model.hmm", ""), (None, "s.faa")):
            with self.subTest(hmm=hmm, shard=shard):
                with self.assertRaises(ValueError):
                    self.module.build_hmmsearch_command(
                        hmm, shard, total_targets=FULL_DATABASE_Z
                    )

    def test_domz_is_never_added(self):
        command = self.module.build_hmmsearch_command(
            "model.hmm", "shard.faa", total_targets=FULL_DATABASE_Z
        )
        self.assertNotIn("--domZ", command)
        self.assertNotIn("--domE", command)
        self.assertFalse([item for item in command if item.startswith("--domZ")])

    def test_evalue_is_forwarded_as_an_argument(self):
        command = self.module.build_hmmsearch_command(
            "model.hmm", "shard.faa", total_targets=FULL_DATABASE_Z, evalue="1e-5"
        )
        self.assertEqual(command[command.index("-E") + 1], "1e-5")

    def test_tblout_and_domtblout_are_separate_flags(self):
        command = self.module.build_hmmsearch_command(
            "model.hmm",
            "shard.faa",
            total_targets=FULL_DATABASE_Z,
            tblout="fam__shard_0001.tbl",
            domtblout="fam__shard_0001.dom",
        )
        self.assertEqual(command[command.index("--tblout") + 1], "fam__shard_0001.tbl")
        self.assertEqual(
            command[command.index("--domtblout") + 1], "fam__shard_0001.dom"
        )

    def test_cpu_is_a_positive_integer_argument(self):
        command = self.module.build_hmmsearch_command(
            "model.hmm", "shard.faa", total_targets=FULL_DATABASE_Z, cpu=1
        )
        self.assertEqual(command[command.index("--cpu") + 1], "1")
        with self.assertRaises(ValueError):
            self.module.build_hmmsearch_command(
                "model.hmm", "shard.faa", total_targets=FULL_DATABASE_Z, cpu=0
            )

    def test_program_is_hmmsearch_and_model_precedes_the_shard(self):
        command = self.module.build_hmmsearch_command(
            "model.hmm", "shard.faa", total_targets=FULL_DATABASE_Z
        )
        self.assertEqual(command[0], "hmmsearch")
        self.assertLess(command.index("model.hmm"), command.index("shard.faa"))

    def test_command_is_subprocess_safe(self):
        """Every element is a literal str; no shell metacharacter interpolation."""
        hostile = "; rm -rf /tmp/x $(id) `id` | tee > out && echo pwned"
        command = self.module.build_hmmsearch_command(
            hostile, hostile, total_targets=FULL_DATABASE_Z, tblout=hostile
        )
        self.assertTrue(all(isinstance(item, str) for item in command))
        self.assertIn(hostile, command)
        # The hostile text survives verbatim as one argv element: nothing was
        # split, quoted and re-joined, or handed to a shell.
        self.assertEqual(command.count(hostile), 3)
        self.assertEqual(self.module.build_hmmsearch_command.__module__, "hmmer_command")

    def test_module_has_no_third_party_imports(self):
        source = MODULE_PATH.read_text(encoding="utf-8")
        imported = set(re.findall(r"^\s*(?:from|import)\s+([A-Za-z_][\w.]*)", source, re.M))
        allowed = {
            "argparse",
            "json",
            "os",
            "sys",
            "pathlib",
            "shlex",
            "typing",
            "__future__",
            "collections",
            "dataclasses",
            "re",
        }
        self.assertFalse(imported - allowed, f"unexpected imports: {imported - allowed}")


class ScreeningEntrypointScaleTests(unittest.TestCase):
    """``06_screen.sh`` must bind every shard call to one full-database Z."""

    def setUp(self):
        self.script = SCREEN_SH.read_text(encoding="utf-8")

    def test_shell_always_passes_the_database_size_to_hmmsearch(self):
        body = run_one_call_body(self.script)
        self.assertIn('hmmsearch --tblout "$out" --domtblout "${out%.tbl}.dom"', body)
        self.assertIn('-Z "$DATABASE_SIZE_Z"', body)
        self.assertIn('-E "$EVAL" --cpu 1 "$hmm" "$shard"', body)

    def test_exactly_one_z_flag_inside_the_per_shard_call(self):
        scope = strip_shell_comments(per_shard_loop_body(self.script))
        self.assertEqual(scope.count('-Z "$DATABASE_SIZE_Z"'), 1)
        body = strip_shell_comments(run_one_call_body(self.script))
        self.assertEqual(body.count('-Z "$DATABASE_SIZE_Z"'), 1)
        self.assertNotIn("--domZ", body)
        self.assertNotIn("-Z 292000000", self.script)

    def test_z_is_scoped_to_the_per_shard_call_not_a_family_constant(self):
        body = run_one_call_body(self.script)
        scope = per_shard_loop_body(self.script)
        self.assertIn("-Z", body)
        # The flag may only appear inside the per-shard call inside the loop,
        # never in the per-family dispatch or run contract.
        without_call = scope.replace(body, "")
        self.assertNotIn("-Z", without_call)

    def test_database_size_flag_is_parsed_and_exported(self):
        self.assertIn("--database-size-z)", self.script)
        self.assertIn('--database-size-z) DATABASE_SIZE_Z="$2"; shift 2 ;;', self.script)
        self.assertIn("export BUILD_DIR EVAL FAILED DATABASE_SIZE_Z", self.script)

    def test_database_size_is_validated_as_a_positive_integer(self):
        self.assertIn("DATABASE_SIZE_Z", self.script)
        self.assertRegex(
            self.script,
            r'\[\[ "\$DATABASE_SIZE_Z" =~ \^\[1-9\]\[0-9\]\*\$\s*\]\]',
        )
        self.assertRegex(self.script, r"\[ERROR\].*--database-size-z")

    def test_absent_database_size_refuses_instead_of_searching_unscaled(self):
        self.assertIn('DATABASE_SIZE_Z="${PHB_DATABASE_SIZE_Z:-}"', self.script)
        self.assertRegex(
            self.script,
            r'if \[ -z "\$DATABASE_SIZE_Z" \]; then[\s\S]{0,900}?exit 1',
        )
        body = run_one_call_body(self.script)
        self.assertNotIn("${DATABASE_SIZE_Z:-", body)
        self.assertNotIn('${DATABASE_SIZE_Z:-1}', self.script)

    def test_run_contract_records_the_database_size_and_its_basis(self):
        self.assertIn("DATABASE_SIZE_BASIS", self.script)
        self.assertIn("database_size_Z", self.script)
        self.assertIn("run_contract", self.script.lower())

    def test_manifest_call_forwards_the_scale(self):
        self.assertIn("--database-size-z \"$DATABASE_SIZE_Z\"", self.script)
        self.assertIn("--database-size-basis \"$DATABASE_SIZE_BASIS\"", self.script)

    def test_server_guard_is_preserved_untouched(self):
        self.assertIn("--server) SERVER=1; shift ;;", self.script)
        self.assertIn('python "$SCRIPT_DIR/server_resources.py"', self.script)
        self.assertIn("MEASURED_LIMIT", self.script)
        self.assertIn("THREADS", self.script)


class TierDatabaseSizeScaleTests(unittest.TestCase):
    """Task F4 — the tier entrypoint must bind one full-database Z (fail-closed).

    ``formal_scan13_tier_processing.sh`` re-runs HMMER through
    ``08c_tier_rescore.py``, so it needs exactly the same scale contract as
    ``06_screen.sh``: an explicit ``--database-size-z``, a positive-integer
    check, a refusal instead of a silent shard-local default, a ``-Z`` reaching
    every HMMER call, and the scale recorded in the tier manifest.
    """

    def setUp(self):
        self.script = TIER_SH.read_text(encoding="utf-8")

    def test_entrypoint_requires_an_explicit_database_size_flag(self):
        self.assertIn("--database-size-z", self.script)
        self.assertIn('DATABASE_SIZE_Z="${PHB_DATABASE_SIZE_Z:-}"', self.script)
        self.assertIn("--database-size-basis", self.script)
        self.assertIn('DATABASE_SIZE_BASIS="${PHB_DATABASE_SIZE_BASIS:-}"', self.script)

    def test_missing_database_size_is_fail_closed(self):
        self.assertRegex(
            self.script,
            r'if \[ -z "\$DATABASE_SIZE_Z" \]; then[\s\S]{0,1200}?exit 1',
        )
        self.assertRegex(self.script, r"\[ERROR\][^\n]*--database-size-z")
        # No fallback value and no derivation from the evidence in the run: the
        # scale is a supplied fact, never a counted default.  (The pre-existing
        # extraction check may still count lines; argparse/tests do.)
        self.assertNotIn("${DATABASE_SIZE_Z:-", self.script)
        for derived in ("grep -c '^>'", "seq_count", "Z_shard"):
            with self.subTest(derived=derived):
                self.assertNotIn(derived, self.script)
        # The one Z is forwarded to the rescoring step, never re-derived there.
        self.assertNotIn('Z="$(', self.script)

    def test_database_size_is_validated_as_a_positive_integer(self):
        self.assertRegex(
            self.script,
            r'\[\[ "\$DATABASE_SIZE_Z" =~ \^\[1-9\]\[0-9\]\*\$\s*\]\]',
        )
        self.assertRegex(
            self.script,
            r"\[ERROR\][^\n]*--database-size-z[^\n]*positive integer",
        )

    def test_scale_gate_precedes_every_side_effect(self):
        """Refusing must happen before the run directory or any copy exists."""
        gate = self.script.index('if [ -z "$DATABASE_SIZE_Z" ]; then')
        for effect in ("create_run_layout", "mkdir -p", 'cp "$SCRIPT_DIR/', "<<'PY'"):
            with self.subTest(effect=effect):
                self.assertLess(gate, self.script.index(effect))

    def test_one_z_is_threaded_to_every_hmm_rerun(self):
        self.assertIn(
            '--database-size-z "$DATABASE_SIZE_Z"', self.script
        )
        self.assertIn(
            '--database-size-basis "$DATABASE_SIZE_BASIS"', self.script
        )
        # The rescoring call itself carries the scale, next to the CPU.
        self.assertRegex(
            self.script,
            r'08c_tier_rescore\.py" --database-size-z "\$DATABASE_SIZE_Z" '
            r'--database-size-basis "\$DATABASE_SIZE_BASIS" --cpu "\$HMM_CPU"',
        )
        # The entrypoint never invokes hmmsearch itself: HMMER flags live only in
        # the rescoring step, which builds them through the shared builder.  The
        # recorded template may only appear as the manifest's literal string.
        shell = strip_shell_comments(self.script)
        self.assertEqual(
            [
                line.strip()
                for line in shell.splitlines()
                if line.strip().startswith("hmmsearch")
            ],
            [],
        )
        self.assertIn("hmmsearch --tblout", shell)

    def test_hmmer_command_module_is_shipped_with_the_rescoring_step(self):
        """The copied rescore step needs the shared builder beside it."""
        self.assertIn("hmmer_command.py", self.script)

    def test_entrypoint_offers_no_hmm_or_domain_flags(self):
        """Only run/scale flags reach the driver; HMMER flags come from the builder."""
        shell = strip_shell_comments(self.script)
        driver_flags = {
            match.group(1)
            for match in re.finditer(r"^\s*(--[a-z][\w-]*)\)", shell, re.M)
        }
        self.assertEqual(
            driver_flags,
            {
                "--run-id",
                "--parent-run",
                "--hmm-cpu",
                "--database-size-z",
                "--database-size-basis",
            },
        )

    def test_tier_hmm_cpu_is_single_threaded_by_default(self):
        values = [
            line.strip()
            for line in self.script.splitlines()
            if line.startswith("HMM_CPU=")
        ]
        self.assertEqual(values, ["HMM_CPU=1"], "per-task HMMER must be single-threaded")
        self.assertIn("--hmm-cpu must be 1..60", self.script)

    def test_manifest_records_the_scale_fields(self):
        block = next(
            blob
            for blob in tier_heredoc_blocks()
            if "tier_processing_manifest.json" in blob
        )
        payload_keys = manifest_payload_keys(block)
        for field in (
            "database_size_Z",
            "database_size_basis",
            "shards",
            "shard_sequence_total",
            "hmmsearch_command_template",
            "domz_used",
        ):
            with self.subTest(field=field):
                self.assertIn(field, payload_keys)
        # The pre-existing manifest fields must survive the change.
        for field in (
            "schema_version",
            "status",
            "run_id",
            "created_utc",
            "parent_run",
            "hmm_cpu",
            "input_contract_sha256",
            "tier_processing_summary_sha256",
            "tier1_counts",
        ):
            with self.subTest(preserved=field):
                self.assertIn(field, payload_keys)
        self.assertIn("count_fasta_records", block)
        self.assertIn("shard_*.faa", block)
        # The recorded command template carries the one -Z and never --domZ.
        self.assertIn("'-Z \"$DATABASE_SIZE_Z\"", block)
        template = block.split("command_template = (", 1)[1].split(")", 1)[0]
        self.assertNotIn("--domZ", template)
    def test_manifest_scale_gate_fails_closed_on_real_shard_counts(self):
        """The shipped gate — not a look-alike — must check real record counts."""
        programs = load_tier_manifest_programs()
        gate = programs.get("require_database_size_scale")
        self.assertTrue(callable(gate), "the tier driver has no scale gate")
        shards = Path(tempfile.mkdtemp(prefix="task_f4_shards_"))
        self.addCleanup(shutil.rmtree, shards, True)
        records = "".join(
            f">protein_{index}\n{'M' * 60}\n" for index in range(3)
        )
        (shards / "shard_0001.faa").write_text(records, encoding="utf-8")
        (shards / "shard_0002.faa").write_text(records, encoding="utf-8")
        with self.assertRaises(SystemExit) as caught:
            gate(str(shards), 6 + 1, "measured parent scan-13 shard counts")
        self.assertNotEqual(caught.exception.code, 0)
        counted, total = gate(str(shards), 6, "measured parent scan-13 shard counts")
        self.assertEqual(total, 6)
        self.assertEqual(
            {item["name"]: item["sequence_count"] for item in counted},
            {"shard_0001.faa": 3, "shard_0002.faa": 3},
        )
        observed = gate(str(shards), "6", "measured parent scan-13 shard counts")
        self.assertEqual(observed[1], 6)

    def test_manifest_scale_gate_rejects_missing_scale_evidence(self):
        """A missing/blank Z or basis makes the manifest invalid, not defaulted."""
        programs = load_tier_manifest_programs()
        gate = programs.get("require_database_size_scale")
        self.assertTrue(callable(gate), "the tier driver has no scale gate")
        shards = Path(tempfile.mkdtemp(prefix="task_f4_shards_"))
        self.addCleanup(shutil.rmtree, shards, True)
        records = "".join(
            f">protein_{index}\n{'M' * 60}\n" for index in range(3)
        )
        (shards / "shard_0001.faa").write_text(records, encoding="utf-8")
        for z, basis in (
            (None, "measured"),
            ("", "measured"),
            ("0", "measured"),
            ("-3", "measured"),
            ("2.5", "measured"),
            (3, ""),
            (3, "   "),
            (3, None),
        ):
            with self.subTest(z=z, basis=basis):
                with self.assertRaises(SystemExit):
                    gate(str(shards), z, basis)

    def test_manifest_scale_gate_rejects_an_empty_shard_dir(self):
        """Without the per-shard counts the scale cannot be proven."""
        programs = load_tier_manifest_programs()
        gate = programs.get("require_database_size_scale")
        self.assertTrue(callable(gate), "the tier driver has no scale gate")
        empty = Path(tempfile.mkdtemp(prefix="task_f4_empty_"))
        self.addCleanup(shutil.rmtree, empty, True)
        with self.assertRaises(SystemExit):
            gate(str(empty), 5, "measured")

    def test_manifest_validator_rejects_missing_scale_fields(self):
        """A manifest missing any scale field must be invalid, never defaulted."""
        programs = load_tier_manifest_programs()
        validate = programs.get("validate_tier_manifest")
        self.assertTrue(callable(validate), "the tier driver has no manifest validator")
        complete = {
            "database_size_Z": FULL_DATABASE_Z,
            "database_size_basis": "measured parent scan-13 shard counts",
            "shards": [{"name": "shard_0001.faa", "sequence_count": FULL_DATABASE_Z}],
            "shard_sequence_total": FULL_DATABASE_Z,
            "hmmsearch_command_template": programs["command_template"],
            "domz_used": False,
        }
        self.assertEqual(validate(dict(complete)), complete)
        for field in (
            "database_size_Z",
            "database_size_basis",
            "shards",
            "shard_sequence_total",
            "hmmsearch_command_template",
            "domz_used",
        ):
            with self.subTest(missing=field):
                broken = dict(complete)
                del broken[field]
                with self.assertRaises(SystemExit):
                    validate(broken)
        for override in (
            {"shard_sequence_total": FULL_DATABASE_Z // 2},
            {"domz_used": True},
            {"hmmsearch_command_template": "hmmsearch -E 1e-5 model.hmm shard.faa"},
            {"hmmsearch_command_template": '-Z "$DATABASE_SIZE_Z" --domZ 1e-5'},
        ):
            with self.subTest(override=override):
                broken = dict(complete)
                broken.update(override)
                with self.assertRaises(SystemExit):
                    validate(broken)


class TierManifestEndToEndTests(unittest.TestCase):
    """Run the driver's own manifest program against a synthetic parent run."""

    def setUp(self):
        self.workdir = Path(tempfile.mkdtemp(prefix="task_f4_manifest_"))
        self.addCleanup(shutil.rmtree, self.workdir, True)
        self.run = self.workdir / "runs" / "20260928_tier_f4_01"
        self.parent = self.workdir / "runs" / "20260901_formal_frozen_scan_13"
        self.shards = self.parent / "inputs" / "scan_shards"
        for path in (
            self.run / "inputs",
            self.run / "results" / "tier_processing" / "screen",
            self.run / "results" / "tier_processing" / "data" / "screen" / "tiers",
            self.shards,
        ):
            path.mkdir(parents=True, exist_ok=True)
        records = "".join(f">protein_{index}\n{'M' * 60}\n" for index in range(3))
        for name in ("shard_0001.faa", "shard_0002.faa"):
            (self.shards / name).write_text(records, encoding="utf-8")
        (self.run / "inputs" / "hmmer_command.py").write_text(
            "PROGRAM = 'hmmsearch'\n", encoding="utf-8"
        )
        (self.run / "input_contract.json").write_text("{}\n", encoding="utf-8")
        (self.run / "results" / "tier_processing" / "screen" / "summary.txt").write_text(
            "summary\n", encoding="utf-8"
        )
        (self.run / "results" / "tier_processing" / "data" / "screen" / "tiers"
         / "ePhaZ_tier1.faa").write_text(">p1\nMMM\n", encoding="utf-8")
        self.saved_argv = list(sys.argv)
        self.addCleanup(setattr, sys, "argv", self.saved_argv)

    def run_manifest_program(self, database_size_z, basis):
        """Run the driver's manifest builder and write it exactly as main() does."""
        programs = load_tier_manifest_programs()
        with redirect_stdout(io.StringIO()):
            manifest = programs["build_tier_processing_manifest"](
                self.run, self.parent, 1, str(database_size_z), basis
            )
            (self.run / "results" / "tier_processing_manifest.json").write_text(
                json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )

    def test_manifest_is_written_with_the_full_scale(self):
        self.run_manifest_program(6, "measured parent scan-13 shard counts")
        manifest = json.loads(
            (self.run / "results" / "tier_processing_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(manifest["database_size_Z"], 6)
        self.assertEqual(manifest["database_size_basis"],
                         "measured parent scan-13 shard counts")
        self.assertEqual(manifest["shard_sequence_total"], 6)
        self.assertEqual(
            {item["name"]: item["sequence_count"] for item in manifest["shards"]},
            {"shard_0001.faa": 3, "shard_0002.faa": 3},
        )
        self.assertIn('-Z "$DATABASE_SIZE_Z"', manifest["hmmsearch_command_template"])
        self.assertNotIn("--domZ", manifest["hmmsearch_command_template"])
        self.assertIs(manifest["domz_used"], False)
        # Pre-existing manifest evidence survives.
        self.assertEqual(manifest["run_id"], self.run.name)
        self.assertEqual(manifest["hmm_cpu"], 1)
        self.assertEqual(manifest["tier1_counts"], {"ePhaZ": 1})
        self.assertIn("input_contract_sha256", manifest)
        self.assertIn("hmmer_command_sha256", manifest)

    def test_manifest_is_not_written_when_the_scale_contradicts_the_counts(self):
        with self.assertRaises(SystemExit) as caught:
            self.run_manifest_program(7, "measured parent scan-13 shard counts")
        self.assertNotEqual(caught.exception.code, 0)
        self.assertFalse(
            (self.run / "results" / "tier_processing_manifest.json").exists(),
            "an invalid manifest must not be published",
        )

    def test_manifest_is_not_written_without_a_basis(self):
        with self.assertRaises(SystemExit):
            self.run_manifest_program(6, "   ")
        self.assertFalse(
            (self.run / "results" / "tier_processing_manifest.json").exists()
        )

    def test_manifest_is_not_written_when_the_library_shards_are_absent(self):
        shutil.rmtree(self.shards)
        with self.assertRaises(SystemExit):
            self.run_manifest_program(6, "measured parent scan-13 shard counts")
        self.assertFalse(
            (self.run / "results" / "tier_processing_manifest.json").exists()
        )


class TierRescoreScaleTests(unittest.TestCase):
    """Task F4 — every rescoring HMMER call shares one full-database Z."""

    def setUp(self):
        self.module = load_rescore_module()
        self.source = TIER_RESCORE.read_text(encoding="utf-8")
        self.script = TIER_SH.read_text(encoding="utf-8")
        self.calls = []
        self.workdir = None
        self._original_run = self.module.subprocess.run

    def tearDown(self):
        self.module.subprocess.run = self._original_run

    def run_rescore(self, *argv, expect_families=None):
        """Run the real ``main()`` with HMMER replaced by a recorder.

        Returns the argv of every HMMER call the step would have executed, in
        order.  Nothing is executed: ``subprocess.run`` is replaced.
        """
        module = self.module
        workdir = Path(tempfile.mkdtemp(prefix="task_f4_rescore_"))
        self.workdir = workdir
        self.addCleanup(shutil.rmtree, workdir, True)
        seqdir = workdir / "data" / "screen" / "family_seqs"
        seqdir.mkdir(parents=True, exist_ok=True)
        families = module.CURATED if expect_families is None else expect_families
        for name in families:
            (seqdir / f"{name}_validated.faa").write_text(
                f">protein_1\n{'M' * 60}\n", encoding="utf-8"
            )
            hmm = workdir / module.CURATED[name]
            hmm.parent.mkdir(parents=True, exist_ok=True)
            hmm.write_text("HMMER3/f\n", encoding="utf-8")

        def fake_run(command, **kwargs):
            self.calls.append(command)
            Path(command[command.index("--tblout") + 1]).write_text(
                "# tbl\n", encoding="utf-8"
            )

        module.subprocess.run = fake_run
        previous_cwd = os.getcwd()
        try:
            os.chdir(workdir)
            with redirect_stdout(io.StringIO()):
                module.main(list(argv))
        finally:
            os.chdir(previous_cwd)
        return self.calls

    def test_z_reaches_the_command_from_the_documented_override(self):
        previous = os.environ.get("PHB_DATABASE_SIZE_Z")
        os.environ["PHB_DATABASE_SIZE_Z"] = str(FULL_DATABASE_Z)
        try:
            # No --database-size-z on the command line: the documented
            # PHB_DATABASE_SIZE_Z override is the only source left.
            commands = self.run_rescore("--cpu", "1")
        finally:
            if previous is None:
                os.environ.pop("PHB_DATABASE_SIZE_Z", None)
            else:
                os.environ["PHB_DATABASE_SIZE_Z"] = previous
        self.assertTrue(commands)
        for command in commands:
            with self.subTest(command=command[:3]):
                self.assertEqual(
                    command[command.index("-Z") + 1], str(FULL_DATABASE_Z)
                )
        self.assertEqual(
            {command[command.index("-Z") + 1] for command in commands},
            {str(FULL_DATABASE_Z)},
        )

    def test_rescore_module_never_defines_a_shard_local_z(self):
        tree = ast.parse(self.source)
        constants = set()
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        constants.add(target.id)
        self.assertIn("DATABASE_SIZE_ENV_VAR", constants)
        # Only the environment-variable *names* may be module constants; a
        # module-level Z value can silently stand in for supplied evidence.
        z_value_constants = {
            name
            for name in constants
            if name not in {"DATABASE_SIZE_ENV_VAR", "DATABASE_SIZE_BASIS_ENV_VAR"}
            and ("Z" in name.upper() and ("SIZE" in name.upper() or name.upper().endswith("_Z")))
        }
        self.assertFalse(
            z_value_constants,
            f"a module-level Z constant can silently stand in for evidence: "
            f"{z_value_constants}",
        )

    def test_rescore_module_does_not_spawn_hmmsearch_outside_the_builder(self):
        tree = ast.parse(self.source)
        built = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and called_name(node) == "build_hmmsearch_command"
        ]
        self.assertTrue(built, "the rescore step must build every call via the builder")
        # Whatever runs a process must be handed an argv built by the shared
        # builder (`hmmsearch()`), not a hand-written literal command line.
        argv_names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and isinstance(node.value, ast.BinOp):
                for name in node.targets:
                    if isinstance(name, ast.Name):
                        argv_names.add(name.id)
        spawned = 0
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            if not isinstance(node.func, ast.Attribute) or node.func.attr != "run":
                continue
            spawned += 1
            first = node.args[0]
            if isinstance(first, ast.Name):
                self.assertIn(
                    first.id, argv_names,
                    "subprocess.run must receive the builder's argv, never a "
                    "hand-written literal command line",
                )
            else:
                self.assertEqual(called_name(first), "build_hmmsearch_command")
        self.assertEqual(spawned, 1, "the rescore step must have exactly one spawn site")

    def test_rescore_has_a_required_database_size_flag_and_no_domz(self):
        code = python_code_text(self.source)
        flags = flag_arguments(self.source)
        self.assertIn("--database-size-z", flags)
        self.assertIn("--cpu", flags)
        self.assertIn("PHB_DATABASE_SIZE_Z", code)
        self.assertIn("build_hmmsearch_command", code)
        # --domZ may appear in prose that says it is never emitted, but it is
        # never a real flag of the rescoring step.
        self.assertNotIn("--domZ", flags)
        # No raw '-Z' may be added by hand: exactly one quoted flag token, so the
        # scale cannot be duplicated or defaulted.
        self.assertEqual(re.findall(r"""["']\s*-Z["']""", code), ["'-Z'"])

    def test_rescore_reads_no_other_environment_variable(self):
        # The rescore step may read the full-library Z and its stated basis from
        # the documented environment overrides — and nothing else that could
        # smuggle in a second scale (for example a shard count).
        self.assertEqual(
            environ_names_read_in_main(self.source),
            {"DATABASE_SIZE_ENV_VAR", "DATABASE_SIZE_BASIS_ENV_VAR"},
            "only the two documented environment constants may be consulted",
        )
        self.assertEqual(
            re.findall(r'^(DATABASE_SIZE\w*_ENV_VAR) = "([^"]+)"', self.source, re.M),
            [
                ("DATABASE_SIZE_ENV_VAR", "PHB_DATABASE_SIZE_Z"),
                ("DATABASE_SIZE_BASIS_ENV_VAR", "PHB_DATABASE_SIZE_BASIS"),
            ],
            "those two constants must name the documented variables",
        )

    def test_every_requested_family_and_tier_is_hmm_scaled(self):
        """Two HMMER calls per family, all with the same full-database Z."""
        calls = self.run_rescore(
            "--database-size-z", str(FULL_DATABASE_Z), "--cpu", "1"
        )
        self.assertEqual(len(calls), 2 * len(self.module.CURATED))
        self.assertEqual(
            {command[command.index("-Z") + 1] for command in calls},
            {str(FULL_DATABASE_Z)},
        )
        for command in calls:
            with self.subTest(tbl=command[command.index("--tblout") + 1]):
                self.assertNotIn("--domZ", command)
                self.assertEqual(command[command.index("--cpu") + 1], "1")
        # Both tiers of every family were searched, not just the first.
        tblouts = " ".join(command[command.index("--tblout") + 1] for command in calls)
        for family in self.module.CURATED:
            with self.subTest(family=family):
                self.assertIn(f"{family}_tier1.tbl", tblouts)
                self.assertIn(f"{family}_tier2.tbl", tblouts)

    def test_rescore_refuses_without_the_database_size(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as caught:
                self.module.main([])
        self.assertNotEqual(caught.exception.code, 0)
        self.assertIn("--database-size-z", stderr.getvalue())
        self.assertIn("positive integer", stderr.getvalue())
        self.assertIn("no default", stderr.getvalue())

    def test_rescore_rejects_non_positive_or_non_numeric_z(self):
        for bad in ("0", "-1", "2.5", "292000000.0", "abc", "1e5", "", " "):
            with self.subTest(database_size_z=bad):
                stderr = io.StringIO()
                with redirect_stderr(stderr):
                    with self.assertRaises(SystemExit) as caught:
                        self.module.main(["--database-size-z", bad])
                self.assertNotEqual(caught.exception.code, 0)
                self.assertIn("positive integer", stderr.getvalue())

    def test_rescore_accepts_a_positive_integer_z(self):
        self.assertEqual(self.module.validate_database_size_z("292000000"), 292000000)
        self.assertEqual(self.module.validate_database_size_z(292000000), 292000000)
        self.assertEqual(self.module.validate_database_size_z(" 7 "), 7)
        for bad in (0, -1, None, True, False, 2.5, "", "   ", "-1", "0", "1.5"):
            with self.subTest(database_size_z=bad):
                with self.assertRaises(ValueError):
                    self.module.validate_database_size_z(bad)


class ScreenManifestScaleTests(unittest.TestCase):
    """The manifest must prove that all shard E-values share one scale."""

    def setUp(self):
        self.validator = load_validator()
        self.tmp = Path(tempfile.mkdtemp(prefix="task9_manifest_"))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.shard_dir = self.tmp / "shards"
        self.hmm_dir = self.tmp / "hmms"
        self.hmmout = self.tmp / "hmmsearch"
        for path in (self.shard_dir, self.hmm_dir, self.hmmout):
            path.mkdir(parents=True, exist_ok=True)
        self.counts = {"shard_0001": 3, "shard_0002": 2}
        self.total = 5
        for name, count in self.counts.items():
            residues = "M" * 60
            records = "".join(
                f">protein_{index}\n{residues}\n" for index in range(count)
            )
            (self.shard_dir / f"{name}.faa").write_text(records, encoding="utf-8")
        (self.hmm_dir / "familyA.hmm").write_text("HMMER3/f\n", encoding="utf-8")
        for name in self.counts:
            (self.hmmout / f"familyA__{name}.tbl").write_text("# tbl\n", encoding="utf-8")
            (self.hmmout / f"familyA__{name}.dom").write_text("# dom\n", encoding="utf-8")

    def build(self, **overrides):
        kwargs = dict(
            shard_dir=str(self.shard_dir),
            hmm_dir=str(self.hmm_dir),
            hmmout=str(self.hmmout),
            families=["familyA"],
            evalue="1e-5",
            database_size_z=self.total,
            database_size_basis="sum of per-shard sequence counts",
        )
        kwargs.update(overrides)
        return self.validator.build_manifest(**kwargs)

    def test_happy_path_records_the_full_scale(self):
        manifest = self.build()
        self.assertEqual(manifest["database_size_Z"], self.total)
        self.assertEqual(
            manifest["database_size_basis"], "sum of per-shard sequence counts"
        )
        self.assertEqual(manifest["shard_sequence_total"], self.total)
        self.assertEqual(
            {shard["name"]: shard["sequence_count"] for shard in manifest["shards"]},
            {f"{name}.faa": count for name, count in self.counts.items()},
        )
        self.assertEqual(
            manifest["hmmsearch_command_template"],
            self.validator.COMMAND_TEMPLATE,
        )
        self.assertIn('-Z "$DATABASE_SIZE_Z"', self.validator.COMMAND_TEMPLATE)
        self.assertNotIn("--domZ", self.validator.COMMAND_TEMPLATE)
        self.validator.validate_manifest(manifest)

    def test_every_task_carries_the_same_z(self):
        manifest = self.build()
        observed = {task["database_size_Z"] for task in manifest["tasks"]}
        self.assertEqual(observed, {self.total})

    def test_shard_runs_with_different_z_values_fail_validation(self):
        manifest = self.build()
        manifest["tasks"][0]["database_size_Z"] = self.total // 2
        with self.assertRaises(self.validator.ScreenManifestError) as caught:
            self.validator.validate_manifest(manifest)
        self.assertIn("database_size_Z", str(caught.exception))

    def test_shard_z_that_contradicts_the_declared_scale_fails(self):
        manifest = self.build()
        manifest["tasks"][0]["database_size_Z"] = 2923820  # shard-local default Z
        with self.assertRaises(self.validator.ScreenManifestError):
            self.validator.validate_manifest(manifest)

    def test_missing_database_size_z_is_not_defaulted(self):
        manifest = self.build()
        del manifest["database_size_Z"]
        with self.assertRaises(self.validator.ScreenManifestError) as caught:
            self.validator.validate_manifest(manifest)
        self.assertIn("database_size_Z", str(caught.exception))

    def test_legacy_manifest_without_database_size_is_rejected_by_default(self):
        """A missing scale must never be silently defaulted to the shard count."""
        legacy = {
            "schema_version": 1,
            "families": ["familyA"],
            "shards": [{"name": "shard_0001", "sha256": "a" * 64}],
            "tasks": [],
        }
        with self.assertRaises(self.validator.ScreenManifestError):
            self.validator.validate_manifest(legacy)

    def test_legacy_manifest_can_only_be_accepted_explicitly(self):
        legacy = {
            "schema_version": 1,
            "families": ["familyA"],
            "shards": [{"name": "shard_0001", "sha256": "a" * 64}],
            "tasks": [
                {
                    "family": "familyA",
                    "shard": "shard_0001",
                    "input_sha256": "a" * 64,
                    "hmm_sha256": "b" * 64,
                    "tbl_sha256": "c" * 64,
                    "dom_sha256": "d" * 64,
                    "evalue": "1e-5",
                }
            ],
        }
        accepted = self.validator.validate_manifest(
            legacy, require_database_size=False
        )
        self.assertNotIn("database_size_Z", accepted)

    def test_build_requires_an_explicit_positive_z(self):
        for bad in (None, 0, -5, "292000000", 1.5):
            with self.subTest(database_size_z=bad):
                with self.assertRaises(self.validator.ScreenManifestError):
                    self.build(database_size_z=bad)

    def test_build_requires_a_documented_basis(self):
        for bad in (None, "", "   "):
            with self.subTest(database_size_basis=bad):
                with self.assertRaises(self.validator.ScreenManifestError):
                    self.build(database_size_basis=bad)

    def test_shard_sequence_total_must_equal_the_declared_z(self):
        with self.assertRaises(self.validator.ScreenManifestError) as caught:
            self.build(database_size_z=self.total + 1)
        self.assertIn("database_size_Z", str(caught.exception))

    def test_shard_sequence_count_is_the_real_record_count(self):
        manifest = self.build()
        counts = {shard["name"]: shard["sequence_count"] for shard in manifest["shards"]}
        self.assertEqual(
            counts, {f"{name}.faa": count for name, count in self.counts.items()}
        )
        self.assertEqual(sum(counts.values()), manifest["shard_sequence_total"])

    def test_shard_sequence_count_mismatch_fails_validation(self):
        manifest = self.build()
        manifest["shards"][0]["sequence_count"] = 99
        with self.assertRaises(self.validator.ScreenManifestError) as caught:
            self.validator.validate_manifest(manifest)
        self.assertIn("sequence_count", str(caught.exception))

    def test_command_template_is_recorded_and_has_no_domz(self):
        manifest = self.build()
        template = manifest["hmmsearch_command_template"]
        self.assertIn("-Z", template)
        self.assertIn("-E", template)
        self.assertIn("--tblout", template)
        self.assertIn("--domtblout", template)
        self.assertNotIn("--domZ", template)


class RescaleConventionTests(unittest.TestCase):
    """Task 9 must reuse the project's established rescaling convention."""

    def setUp(self):
        spec = importlib.util.spec_from_file_location(
            "parse_phaded_cys_targeted_recall_task9",
            SCRIPTS / "parse_phaded_cys_targeted_recall.py",
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.rescale_evalue = module.rescale_evalue

    def test_rescale_reuses_the_established_convention(self):
        # shard-local default Z -> full-library Z is a pure multiplication.
        scaled = self.rescale_evalue(1.4e-08, 2923820, 292000000)
        self.assertEqual(scaled, 1.4e-08 * (292000000 / 2923820))
        # Cross-check against the frozen scan-13 evidence: the same target
        # scored at Z=2,923,820 (shard_0051 default Z) gave E=1.4e-08, and at
        # Z=109,087 gave E=5.1e-10.  The ratio is 26.80 predicted vs 27.45
        # measured, i.e. HMMER's two-significant-figure output precision.
        self.assertAlmostEqual(
            self.rescale_evalue(1.4e-08, 2923820, 109087) / 5.1e-10, 1.0, places=1
        )

    def test_rescale_rejects_non_positive_database_sizes(self):
        for from_z, to_z in ((0, 10), (10, 0), (-1, 10), (10, -1)):
            with self.subTest(from_z=from_z, to_z=to_z):
                with self.assertRaises(ValueError):
                    self.rescale_evalue(1e-5, from_z, to_z)

    def test_rescaling_to_a_smaller_database_shrinks_the_evalue(self):
        self.assertLess(self.rescale_evalue(1e-5, 292000000, 2923820), 1e-5)


if __name__ == "__main__":
    unittest.main()
