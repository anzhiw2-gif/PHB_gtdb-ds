"""Tests for the candidate-only structure-review preparation script.

Two measured defects of
``pipeline/scripts/prepare_phaded_structure_review.py`` are pinned here:

* the hard-coded candidate expectation that aborted every target set other than
  the original pilot, and
* a ``structure_predictor`` manifest block that asserted "no predictor
  available" without ever probing.

The removed literal is never written in this module either: it is assembled at
runtime from two one-character strings, so the static guard below can scan this
test module as well as the script under test.  The guard exists so the defect
cannot silently return.
"""

import json
import re
import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.prepare_phaded_structure_review import (
    DEFERRED_DISPOSITION,
    KEEP,
    build_parser,
    prepare,
    probe_structure_tools,
    summarize_tool_probe,
)

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "prepare_phaded_structure_review.py"

#: The removed hard-coded expectation, assembled at runtime on purpose.
OLD_COUNT_LITERAL = "".join(("1", "9"))
FORBIDDEN_LITERAL = re.compile(r"(?<!\d)" + re.escape(OLD_COUNT_LITERAL) + r"(?!\d)")

DECISION_1 = "retain_dPHASCL1_like_candidate"
DECISION_2 = "retain_extracellular_dPHASCL2_like_candidate"
BASE_HEADER = "source_accession\tphylo_decision\tnearest_reference_accession"


def write_inputs(root: Path, rows, *, header: str = BASE_HEADER):
    """Write a synthetic decisions TSV plus the matching candidate FASTA."""
    decisions = root / "decisions.tsv"
    candidates = root / "candidates.faa"
    decisions.write_text(
        header + "\n" + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )
    candidates.write_text("".join(f">{row[0]}\nMKT\n" for row in rows), encoding="ascii")
    return decisions, candidates


def keep_rows(count: int, *, decision: str = DECISION_1, prefix: str = "A"):
    return [(f"{prefix}{index}", decision, f"R{index}") for index in range(count)]


def read_manifest(output_dir: Path) -> dict:
    return json.loads((output_dir / "structure_review_manifest.json").read_text(encoding="utf-8"))


def stub_probe(*, available: bool):
    """Build a probe report through the real probe with injected collaborators."""
    if available:
        return probe_structure_tools(
            which=lambda name: "/opt/tools/foldseek" if name == "foldseek" else None,
            path_exists=lambda path: str(path) == "/opt/tools/foldseek",
            runner=lambda command: (0, "foldseek version 463739e0", ""),
            sha256_file=lambda path: "a" * 64,
        )
    return probe_structure_tools(
        which=lambda name: None,
        path_exists=lambda path: False,
        runner=lambda command: (1, "", "should not run"),
        sha256_file=lambda path: (_ for _ in ()).throw(AssertionError("no hash")),
    )


class ExpectedCandidateCountTests(unittest.TestCase):
    def test_builds_target_sets_of_any_size(self):
        for size in (0, 1, 987, 5000):
            with self.subTest(size=size), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                decisions, candidates = write_inputs(root, keep_rows(size))
                output_dir = root / "out"
                manifest = prepare(decisions, candidates, output_dir)
                self.assertEqual(manifest["candidate_count"], size)
                contract = manifest["candidate_count_contract"]
                self.assertEqual(contract["observed_candidates"], size)
                self.assertIsNone(contract["expected_candidates"])
                self.assertFalse(contract["enforced"])
                self.assertEqual(contract["status"], "no_fixed_expectation")
                written = (output_dir / "structure_candidates.tsv").read_text(encoding="utf-8").splitlines()
                self.assertEqual(len(written) - 1, size)
                fasta = (output_dir / "structure_candidates.faa").read_text(encoding="ascii")
                self.assertEqual(fasta.count(">"), size)

    def test_declared_expectation_can_be_satisfied_and_is_recorded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, keep_rows(987))
            manifest = prepare(
                decisions,
                candidates,
                root / "out",
                expected_candidates=987,
                expected_candidates_source="p4_input_contract",
            )
            contract = manifest["candidate_count_contract"]
            self.assertEqual(contract["expected_candidates"], 987)
            self.assertEqual(contract["expected_candidates_source"], "p4_input_contract")
            self.assertTrue(contract["enforced"])
            self.assertEqual(contract["status"], "satisfied")

    def test_declared_expectation_mismatch_names_the_declared_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, keep_rows(3))
            with self.assertRaises(ValueError) as caught:
                prepare(
                    decisions,
                    candidates,
                    root / "out",
                    expected_candidates=987,
                    expected_candidates_source="--expected-candidates",
                )
            message = str(caught.exception)
            self.assertIn("987", message)
            self.assertIn("--expected-candidates", message)
            self.assertIn("observed 3", message)
            self.assertNotIn(OLD_COUNT_LITERAL, message)

    def test_no_fixed_expectation_by_default(self):
        # The old defect: this single non-selected row used to abort the build.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(
                root, [("A", "PhaC_vs_PhaZ_unresolved", "R")]
            )
            manifest = prepare(decisions, candidates, root / "out")
            self.assertEqual(manifest["candidate_count"], 0)
            self.assertEqual(manifest["candidate_count_contract"]["status"], "no_fixed_expectation")

    def test_rejects_a_negative_declared_expectation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, keep_rows(1))
            with self.assertRaises(ValueError):
                prepare(decisions, candidates, root / "out", expected_candidates=-1)

    def test_cli_exposes_the_expected_candidate_arguments(self):
        parser = build_parser()
        args = parser.parse_args(
            [
                "--decisions", "d.tsv",
                "--candidates", "c.faa",
                "--output-dir", "out",
                "--expected-candidates", "987",
            ]
        )
        self.assertEqual(args.expected_candidates, 987)
        args = parser.parse_args(["--decisions", "d", "--candidates", "c", "--output-dir", "o"])
        self.assertIsNone(args.expected_candidates)


class StructureToolProbeTests(unittest.TestCase):
    def test_available_records_path_version_and_sha256(self):
        report = probe_structure_tools(
            which=lambda name: "/opt/tools/foldseek" if name == "foldseek" else None,
            path_exists=lambda path: str(path) == "/opt/tools/foldseek",
            runner=lambda command: (0, "foldseek version 463739e0\nsecond line", ""),
            sha256_file=lambda path: "b" * 64,
        )
        record = report["Foldseek"]
        self.assertEqual(record["status"], "available")
        self.assertEqual(record["path"], "/opt/tools/foldseek")
        self.assertEqual(record["version"], "foldseek version 463739e0")
        self.assertEqual(record["sha256"], "b" * 64)
        self.assertEqual(record["sha256_status"], "computed")
        self.assertTrue(record["measurement_complete"])
        self.assertIn("--version", record["probe_command"])
        self.assertEqual(record["probe_exit_code"], 0)
        self.assertIn("which(foldseek)", record["probe_lookup"])
        for key in ("status", "path", "version", "sha256", "probe_command", "reason"):
            self.assertIn(key, record)

    def test_unavailable_records_probe_command_and_reason(self):
        calls = []
        report = probe_structure_tools(
            which=lambda name: None,
            path_exists=lambda path: False,
            runner=lambda command: calls.append(list(command)) or (1, "", ""),
            sha256_file=lambda path: (_ for _ in ()).throw(AssertionError("no hash")),
        )
        record = report["Foldseek"]
        self.assertEqual(record["status"], "unavailable")
        self.assertIsNone(record["path"])
        self.assertIsNone(record["version"])
        self.assertEqual(record["sha256"], "pending")
        self.assertEqual(calls, [])
        self.assertIn("foldseek", record["probe_command"])
        self.assertTrue(record["reason"])
        self.assertIn("no executable found", record["reason"])
        self.assertFalse(record["measurement_complete"])

    def test_never_fabricates_a_hash(self):
        report = probe_structure_tools(
            which=lambda name: "/opt/tools/foldseek" if name == "foldseek" else None,
            path_exists=lambda path: True,
            runner=lambda command: (0, "foldseek 1.0", ""),
            sha256_file=lambda path: (_ for _ in ()).throw(PermissionError("denied")),
        )
        record = report["Foldseek"]
        self.assertEqual(record["sha256"], "pending")
        self.assertEqual(record["sha256_status"], "pending")
        self.assertFalse(record["measurement_complete"])
        self.assertEqual(record["status"], "available")
        self.assertIn("sha256", record["reason"])

    def test_rejects_a_digest_that_is_not_a_sha256(self):
        report = probe_structure_tools(
            which=lambda name: "/opt/tools/foldseek" if name == "foldseek" else None,
            path_exists=lambda path: True,
            runner=lambda command: (0, "foldseek 1.0", ""),
            sha256_file=lambda path: "not-a-digest",
        )
        record = report["Foldseek"]
        self.assertEqual(record["sha256"], "pending")
        self.assertEqual(record["sha256_status"], "pending")

    def test_declared_tool_path_is_probed_before_path_lookup(self):
        declared = "/srv/tools/foldseek_20260912/foldseek/bin/foldseek"
        report = probe_structure_tools(
            which=lambda name: None,
            path_exists=lambda path: Path(path).as_posix() == declared,
            runner=lambda command: (0, "foldseek 463739e0", ""),
            sha256_file=lambda path: "c" * 64,
            tool_paths={"Foldseek": declared},
        )
        record = report["Foldseek"]
        self.assertEqual(record["status"], "available")
        self.assertEqual(record["path"], declared)
        self.assertEqual(record["sha256"], "c" * 64)
        self.assertIn("declared path", record["probe_lookup"])

    def test_version_probe_failure_is_pending_not_available(self):
        report = probe_structure_tools(
            which=lambda name: "/opt/tools/foldseek" if name == "foldseek" else None,
            path_exists=lambda path: True,
            runner=lambda command: (2, "", "unknown option --version"),
            sha256_file=lambda path: "d" * 64,
        )
        record = report["Foldseek"]
        self.assertEqual(record["status"], "pending")
        self.assertIsNone(record["version"])
        self.assertIn("exited 2", record["reason"])
        summary = summarize_tool_probe(report)
        self.assertEqual(summary["status"], "pending_tool")
        self.assertIn("--version", summary["reason"])
        self.assertIn("exited 2", summary["reason"])

    def test_lookup_failure_is_pending_and_never_raises(self):
        def which(name):
            raise OSError("which is unavailable")

        report = probe_structure_tools(
            which=which,
            path_exists=lambda path: False,
            runner=lambda command: (0, "", ""),
            sha256_file=lambda path: "e" * 64,
        )
        self.assertTrue(report)
        for record in report.values():
            self.assertEqual(record["status"], "pending")
            self.assertIn("OSError", record["reason"])
        summary = summarize_tool_probe(report)
        self.assertEqual(summary["status"], "pending_tool")
        self.assertIn("OSError", summary["reason"])


class ManifestProvenanceTests(unittest.TestCase):
    def test_manifest_separates_tool_availability_from_prediction_performed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, keep_rows(4))
            manifest = prepare(
                decisions,
                candidates,
                root / "out",
                tool_probe=lambda: stub_probe(available=True),
            )
            self.assertEqual(manifest["status"], "planned_not_run")
            predictor = manifest["structure_predictor"]
            self.assertEqual(predictor["status"], "available")
            self.assertTrue(predictor["tool_available"])
            self.assertFalse(predictor["prediction_performed"])
            self.assertIsNone(predictor["prediction_evidence"])
            self.assertEqual(predictor["available_tools"], ["Foldseek"])
            self.assertIn("Foldseek", predictor["candidates_checked"])
            self.assertEqual(predictor["tools"]["Foldseek"]["path"], "/opt/tools/foldseek")
            self.assertEqual(manifest["tool_available"], predictor["tool_available"])
            self.assertEqual(manifest["prediction_performed"], predictor["prediction_performed"])
            self.assertFalse(manifest["prediction_performed"])

    def test_pending_tool_only_when_the_probe_genuinely_failed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, keep_rows(4))
            output_dir = root / "out"
            manifest = prepare(
                decisions, candidates, output_dir, tool_probe=lambda: stub_probe(available=False)
            )
            predictor = manifest["structure_predictor"]
            self.assertEqual(predictor["status"], "unavailable")
            self.assertFalse(predictor["tool_available"])
            self.assertFalse(predictor["prediction_performed"])
            self.assertIn("probe ran to completion", predictor["reason"])

    def test_pending_tool_when_the_probe_itself_raises(self):
        def broken_probe():
            raise RuntimeError("probe backend unavailable")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, keep_rows(4))
            manifest = prepare(decisions, candidates, root / "out", tool_probe=broken_probe)
            predictor = manifest["structure_predictor"]
            self.assertEqual(predictor["status"], "pending_tool")
            self.assertFalse(predictor["tool_available"])
            self.assertFalse(predictor["prediction_performed"])
            self.assertIn("RuntimeError", predictor["reason"])
            self.assertIn("probe backend unavailable", predictor["reason"])
            self.assertEqual(predictor["probe_command"], "not executed")

    def test_pending_tool_when_the_caller_skips_the_probe(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, keep_rows(4))
            manifest = prepare(decisions, candidates, root / "out", probe_tools=False)
            predictor = manifest["structure_predictor"]
            self.assertEqual(predictor["status"], "pending_tool")
            self.assertFalse(predictor["tool_available"])
            self.assertIn("skip-tool-probe", predictor["reason"])

    def test_manifest_keeps_every_pre_existing_key_and_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = keep_rows(16) + keep_rows(3, decision=DECISION_2, prefix="B")
            decisions, candidates = write_inputs(root, rows)
            reference = root / "reference.faa"
            reference.write_text(">sp|REF1|ref\nMKT\n", encoding="ascii")
            output_dir = root / "out"
            manifest = prepare(decisions, candidates, output_dir, reference)
            for key in (
                "status",
                "candidate_count",
                "decision_counts",
                "nearest_reference_controls_written",
                "inputs",
                "outputs",
                "structure_predictor",
                "phenotype_boundary",
            ):
                self.assertIn(key, manifest)
            self.assertEqual(manifest["candidate_count"], 16 + 3)
            self.assertEqual(manifest["decision_counts"][DECISION_1], 16)
            self.assertEqual(manifest["decision_counts"][DECISION_2], 3)
            for name in (
                "structure_candidates.faa",
                "structure_candidates.tsv",
                "nearest_reference_controls.faa",
                "structure_review_manifest.json",
            ):
                self.assertTrue((output_dir / name).is_file(), name)
            self.assertEqual(manifest, read_manifest(output_dir))


class CandidateOnlySemanticsTests(unittest.TestCase):
    def test_deferred_layer_never_enters_counts_and_nothing_is_dropped(self):
        header = BASE_HEADER + "\tprimary_disposition"
        rows = [
            (f"K{index}", DECISION_1, f"R{index}", "function_unresolved") for index in range(5)
        ] + [
            (f"D{index}", "with_lipase_deferred", f"S{index}", DEFERRED_DISPOSITION)
            for index in range(3)
        ]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, rows, header=header)
            output_dir = root / "out"
            manifest = prepare(decisions, candidates, output_dir)
            self.assertEqual(manifest["candidate_count"], 5)
            boundary = manifest["candidate_only_semantics"]
            self.assertEqual(boundary["deferred_disposition"], DEFERRED_DISPOSITION)
            self.assertEqual(boundary["deferred_rows_observed"], 3)
            self.assertEqual(boundary["deferred_rows_selected"], 0)
            self.assertFalse(boundary["deferred_layer_included_in_counts"])
            self.assertFalse(boundary["deletes_candidates"])
            self.assertFalse(boundary["demotes_candidates"])
            self.assertFalse(boundary["excludes_candidates"])
            selection = (output_dir / "structure_candidates.tsv").read_text(encoding="utf-8")
            self.assertNotIn(DEFERRED_DISPOSITION, selection)
            self.assertEqual(selection.count("K"), 5)
            self.assertEqual(
                (output_dir / "structure_candidates.faa").read_text(encoding="ascii").count(">"),
                5,
            )

    def test_non_selected_rows_are_kept_out_of_the_selection_but_nothing_is_written_back(self):
        rows = [("A", DECISION_1, "R"), ("B", "PhaC_vs_PhaZ_unresolved", "S")]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decisions, candidates = write_inputs(root, rows)
            before = decisions.read_text(encoding="utf-8")
            manifest = prepare(decisions, candidates, root / "out")
            self.assertEqual(manifest["candidate_count"], 1)
            self.assertEqual(decisions.read_text(encoding="utf-8"), before)


class NoHardCodedExpectationTests(unittest.TestCase):
    def test_the_removed_literal_no_longer_appears_in_the_sources(self):
        for path in (SCRIPT_PATH, Path(__file__).resolve()):
            text = path.read_text(encoding="utf-8")
            offenders = [
                f"{path.name}:{text[: match.start()].count(chr(10)) + 1}: {match.group(0)}"
                for match in FORBIDDEN_LITERAL.finditer(text)
            ]
            self.assertEqual(offenders, [], f"{OLD_COUNT_LITERAL!r} as an expectation is back: {offenders}")

    def test_the_old_comparison_and_message_are_gone(self):
        text = SCRIPT_PATH.read_text(encoding="utf-8")
        self.assertNotIn(f"!= {OLD_COUNT_LITERAL}", text)
        self.assertNotIn(f"== {OLD_COUNT_LITERAL}", text)
        self.assertNotIn(f"expected {OLD_COUNT_LITERAL}", text)
        self.assertIn("expected_candidates", text)
        self.assertIn("--expected-candidates", text)

    def test_keep_set_is_unchanged_and_never_contains_the_deferred_disposition(self):
        self.assertEqual(
            KEEP,
            {"retain_dPHASCL1_like_candidate", "retain_extracellular_dPHASCL2_like_candidate"},
        )
        self.assertNotIn(DEFERRED_DISPOSITION, KEEP)


if __name__ == "__main__":
    unittest.main()
