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
"""

from __future__ import annotations

import importlib.util
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
MODULE_PATH = SCRIPTS / "hmmer_command.py"
SCREEN_SH = SCRIPTS / "06_screen.sh"
VALIDATOR = SCRIPTS / "06_validate_screen_manifest.py"

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


class HmmerCommandBuilderTests(unittest.TestCase):
    """The shared builder must always pin the full-database ``-Z``."""

    def setUp(self):
        self.module = load_hmmer_command()

    def test_every_shard_uses_the_same_full_database_z(self):
        command = self.module.build_hmmsearch_command(
            "model.hmm", "shard.faa", total_targets=292000000
        )
        self.assertIn("-Z", command)
        self.assertEqual(command[command.index("-Z") + 1], "292000000")

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
