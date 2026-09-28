"""Tests for the discovery-layer / deferred-layer downstream misuse audit (F6, test-first).

Governance rules under audit (AGENTS.md, PhaDED evidence-model redesign):

* ``R1`` the discovery layer (``model_layer="discovery_hmm_uncalibrated"``) is for
  **recall only**; it must never filter, delete or demote a candidate;
* ``R2`` the discovery layer must never produce a family call;
* ``R3`` the discovery layer must never be written into
  ``pipeline/config/formal_scan_models.tsv``;
* ``R4`` the with-lipase (``DED_hfam_2``) deferred layer must be preserved, never deleted;
* ``R5`` the deferred layer must never enter a candidate / high-confidence count;
* ``R6`` ``reference_query_only`` profiles have no HMM at all and must never be scored
  as if they had one;
* ``R7`` ``model_layer`` and ``functional_calibration_status`` are independent and
  neither may be derived from the other.

These tests are synthetic-source tests: they pin the *classification contract* of
:func:`classify_usage` (an AST-based classifier) and the CLI contract of
``pipeline/scripts/audit_discovery_layer_usage.py``.  A final class pins the
governance invariant that the frozen formal-scan registry does not list the
discovery layer; it never pins the current defect inventory, because the
concurrent F3 task is fixing some of those call points.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "pipeline" / "scripts" / "audit_discovery_layer_usage.py"
REGISTRY_PATH = REPO_ROOT / "pipeline" / "config" / "formal_scan_models.tsv"

REQUIRED_COLUMNS = (
    "script", "line", "evidence_kind", "signal", "use_class", "violation",
    "rule_cited", "note",
)
USE_CLASSES = {
    "recall_only", "thresholding", "family_call", "counting",
    "registry_write", "deferred_layer_touch", "unknown",
}
VIOLATION_VALUES = {"true", "false", "unreviewed"}


def _load_module():
    spec = importlib.util.spec_from_file_location("audit_discovery_layer_usage", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _findings(source):
    module = _load_module()
    return module.classify_usage(source, "synthetic.py")


def _violations(findings):
    return [row for row in findings if row["violation"] == "true"]


def _by_signal(findings, signal):
    return [row for row in findings if row["signal"] == signal]


def _call_main(module, argv):
    stdout, stderr = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        code = module.main(argv)
    return code, stdout.getvalue(), stderr.getvalue()


# ---------------------------------------------------------------------------
# R1: the discovery layer recalls and ranks; it never filters
# ---------------------------------------------------------------------------

DISCOVERY_FILTER_SOURCE = '''
import csv


def keep_candidates(rows):
    kept = []
    for row in rows:
        if float(row["discovery_best_evalue"]) > 1e-5:
            continue
        kept.append(row)
    return kept
'''


RECALL_ONLY_SOURCE = '''
def rank_candidates(rows):
    scored = [
        {"accession": row["accession"], "discovery_score": float(row["discovery_score"])}
        for row in rows
    ]
    return sorted(scored, key=lambda row: -row["discovery_score"])
'''


class DiscoveryScoreUseTests(unittest.TestCase):
    def test_discovery_score_filter_is_a_violation(self):
        findings = _findings(DISCOVERY_FILTER_SOURCE)
        violations = _violations(findings)
        self.assertTrue(violations, f"a discovery-score filter must be flagged: {findings}")
        self.assertTrue(
            any(row["use_class"] == "thresholding" for row in violations),
            f"expected a thresholding violation: {violations}",
        )
        self.assertTrue(
            all(row["rule_cited"] == "R1_recall_only_only" for row in violations),
            f"expected rule R1: {violations}",
        )
        self.assertTrue(all(row["evidence_kind"] == "ast" for row in violations))

    def test_discovery_score_used_only_for_recall_is_not_a_violation(self):
        findings = _findings(RECALL_ONLY_SOURCE)
        self.assertEqual(_violations(findings), [], f"ranking is recall-only: {findings}")
        self.assertTrue(
            any(row["use_class"] == "recall_only" for row in findings),
            f"the audit must still report the recall-only use: {findings}",
        )
        self.assertTrue(
            any(row["signal"] == "discovery_score_sort_or_rank" for row in findings),
            f"expected the ranking signal: {findings}",
        )

    def test_discovery_score_count_into_a_candidate_count_is_a_violation(self):
        source = '''
def summarize(rows):
    high_confidence_count = 0
    for row in rows:
        if float(row["discovery_score"]) >= 20.0:
            high_confidence_count += 1
    return {"high_confidence_candidate_count": high_confidence_count}
'''
        findings = _findings(source)
        violations = _violations(findings)
        self.assertTrue(violations, f"a discovery score feeding a candidate count is a violation: {findings}")
        self.assertTrue(
            any(row["rule_cited"] == "R1_recall_only_only" for row in violations),
            f"expected rule R1: {violations}",
        )

    def test_commented_out_filter_is_not_a_violation_although_a_text_scan_flags_it(self):
        source = '''
def load_rows(path):
    """Rows are only scored; candidates are NEVER filtered on the discovery_score."""
    # historical draft, deliberately disabled:
    #     if row["discovery_score"] > 0.5: drop(row)
    return sorted(rows, key=lambda row: -row["discovery_score"])
'''
        module = _load_module()
        findings = module.classify_usage(source, "synthetic.py")
        self.assertEqual(
            _violations(findings), [],
            f"a comment/docstring mention must not be a violation: {findings}",
        )
        text_rows = module.scan_text_signals(source, "synthetic.py")
        self.assertTrue(
            any(row["signal"] == "discovery_score_comparison" for row in text_rows),
            f"the documented text scan is expected to produce this false positive: {text_rows}",
        )
        for row in text_rows:
            self.assertEqual(row["evidence_kind"], "text")

    def test_docstring_governance_claim_is_unreviewed_not_a_violation(self):
        source = '''
def build(rows):
    """The discovery layer never deletes a candidate; the deferred tier is untouched."""
    return list(rows)
'''
        findings = _findings(source)
        self.assertEqual(_violations(findings), [], f"prose is not code: {findings}")
        prose = _by_signal(findings, "prose_only_governance_claim")
        self.assertTrue(prose, f"the unverifiable prose claim must be recorded: {findings}")
        self.assertTrue(all(row["violation"] == "unreviewed" for row in prose))
        self.assertTrue(all(row["use_class"] == "unknown" for row in prose))


# ---------------------------------------------------------------------------
# R2: the discovery layer never produces a family call
# ---------------------------------------------------------------------------

class DiscoveryFamilyCallTests(unittest.TestCase):
    def test_discovery_family_assignment_is_a_violation(self):
        source = '''
def call_family(row):
    sequence_family_call = ""
    if row["discovery_best_evalue"] < 1e-5:
        sequence_family_call = row["discovery_best_family"]
    return {"sequence_family_call": sequence_family_call}
'''
        findings = _findings(source)
        violations = _violations(findings)
        self.assertTrue(violations, f"a discovery-derived family call must be flagged: {findings}")
        self.assertTrue(
            any(row["use_class"] == "family_call" and row["rule_cited"] == "R2_no_family_call"
                for row in violations),
            f"expected a family_call violation citing R2: {violations}",
        )


# ---------------------------------------------------------------------------
# R3: the frozen formal-scan registry
# ---------------------------------------------------------------------------

REGISTRY_WRITE_SOURCE = '''
from pathlib import Path

REGISTRY = Path("pipeline/config/formal_scan_models.tsv")


def add_discovery_model(model_name):
    with REGISTRY.open("w", encoding="utf-8") as handle:
        handle.write(model_name + "\\n")
'''


REGISTRY_READ_SOURCE = '''
import hashlib
from pathlib import Path

REGISTRY = Path("pipeline/config/formal_scan_models.tsv")


def registry_sha256():
    return hashlib.sha256(REGISTRY.read_bytes()).hexdigest()
'''


class RegistryTests(unittest.TestCase):
    def test_registry_write_is_flagged_as_registry_write_violation(self):
        findings = _findings(REGISTRY_WRITE_SOURCE)
        rows = _by_signal(findings, "formal_scan_models_reference")
        self.assertTrue(rows, f"the registry reference must be recorded: {findings}")
        self.assertTrue(all(row["use_class"] == "registry_write" for row in rows))
        self.assertTrue(all(row["evidence_kind"] == "text" for row in rows))
        violations = _violations(findings)
        self.assertTrue(violations, f"writing the registry is a violation: {findings}")
        self.assertTrue(all(row["rule_cited"] == "R3_registry_exclusion" for row in violations))

    def test_registry_read_only_reference_is_not_a_violation(self):
        findings = _findings(REGISTRY_READ_SOURCE)
        rows = _by_signal(findings, "formal_scan_models_reference")
        self.assertTrue(rows, f"the registry reference must still be recorded: {findings}")
        self.assertTrue(all(row["use_class"] == "registry_write" for row in rows))
        self.assertEqual(_violations(findings), [], f"a read-only hash is compliant: {findings}")
        self.assertTrue(all("read_only=true" in row["note"] for row in rows))


# ---------------------------------------------------------------------------
# R4/R5: the with-lipase deferred layer
# ---------------------------------------------------------------------------

DEFERRED_READ_SOURCE = '''
from pathlib import Path

ARCHIVE = Path(
    "runs/20260919_phaded_with_lipase_deferred_archive_01/results/"
    "pool_external_with_lipase_deferred.tsv"
)


def read_deferred_rows():
    with ARCHIVE.open(encoding="utf-8") as handle:
        return handle.read().splitlines()
'''


DEFERRED_WRITE_SOURCE = '''
from pathlib import Path


def rebuild(out_dir):
    path = Path(out_dir) / "pool_external_with_lipase_deferred.tsv"
    with path.open("w", encoding="utf-8") as handle:
        handle.write("protein_id\\n")
    return path
'''


DEFERRED_DELETE_SOURCE = '''
import shutil
from pathlib import Path

ARCHIVE = Path("runs/20260919_phaded_with_lipase_deferred_archive_01")


def prune():
    shutil.rmtree(ARCHIVE)
'''


class DeferredLayerTests(unittest.TestCase):
    def test_deferred_archive_read_is_a_read_only_touch(self):
        findings = _findings(DEFERRED_READ_SOURCE)
        rows = _by_signal(findings, "with_lipase_deferred_archive_reference")
        self.assertTrue(rows, f"the deferred archive touch must be recorded: {findings}")
        self.assertTrue(all(row["use_class"] == "deferred_layer_touch" for row in rows))
        self.assertTrue(all(row["evidence_kind"] == "text" for row in rows))
        self.assertTrue(all("read_only=true" in row["note"] for row in rows))
        self.assertEqual(_violations(findings), [], f"reading the layer is legitimate: {findings}")

    def test_deferred_archive_write_is_a_touch_but_not_a_violation(self):
        findings = _findings(DEFERRED_WRITE_SOURCE)
        rows = _by_signal(findings, "with_lipase_deferred_archive_reference")
        self.assertTrue(rows, f"the deferred archive touch must be recorded: {findings}")
        self.assertTrue(all("read_only=false" in row["note"] for row in rows))
        self.assertEqual(
            _violations(findings), [],
            f"rebuilding the deferred tier preserves it (R4): {findings}",
        )

    def test_deferred_archive_delete_is_a_violation(self):
        findings = _findings(DEFERRED_DELETE_SOURCE)
        violations = _violations(findings)
        self.assertTrue(violations, f"deleting the deferred layer breaks R4: {findings}")
        self.assertTrue(
            any(row["rule_cited"] == "R4_deferred_preserved" for row in violations),
            f"expected rule R4: {violations}",
        )

    def test_deferred_layer_entering_a_candidate_count_is_a_violation(self):
        source = '''
def build_summary(deferred_rows, candidates):
    high_confidence_candidates = list(candidates) + list(deferred_rows)
    return {"high_confidence_count": len(high_confidence_candidates)}
'''
        findings = _findings(source)
        # No deferred-layer token literal in this file: the construct is only
        # decidable from the caller's data, so the audit must not invent one.
        self.assertEqual(_violations(findings), [], f"no discovery/deferred token here: {findings}")


# ---------------------------------------------------------------------------
# R6/R7 and the legacy model_status layer proxy
# ---------------------------------------------------------------------------

class ReferenceQueryOnlyTests(unittest.TestCase):
    def test_scoring_reference_query_only_profile_is_a_violation(self):
        source = '''
def score_manifest(manifest, hmmsearch, fasta):
    for row in manifest:
        if row["model_status"] == "reference_only":
            hmmsearch(fasta, row["hmm_path"])
    return "done"
'''
        findings = _findings(source)
        violations = _violations(findings)
        self.assertTrue(violations, f"reference_query_only has no HMM to score: {findings}")
        self.assertTrue(
            any(row["rule_cited"] == "R6_reference_query_only_unscoreable" for row in violations),
            f"expected rule R6: {violations}",
        )

    def test_excluding_reference_query_only_from_scoring_is_compliant(self):
        source = '''
def unambiguous(value):
    """Discriminating layers only; reference_query_only has no HMM."""
    return value not in {"reference_query_only"}
'''
        findings = _findings(source)
        rows = _by_signal(findings, "reference_query_only_model_layer")
        self.assertTrue(rows, f"the reference_query_only mention must be recorded: {findings}")
        self.assertEqual(_violations(findings), [], f"excluding it is compliant: {findings}")


class LayerStatusIndependenceTests(unittest.TestCase):
    def test_calibration_status_derived_from_layer_is_a_violation(self):
        source = '''
def functional_calibration_status(row):
    functional_calibration_status = "not_function_calibrated"
    if row["model_layer"] == "calibrated_candidate_model":
        functional_calibration_status = "calibrated_candidate_model"
    return functional_calibration_status
'''
        findings = _findings(source)
        violations = _violations(findings)
        self.assertTrue(violations, f"deriving the status from the layer breaks R7: {findings}")
        self.assertTrue(
            any(row["rule_cited"] == "R7_layer_status_independent" for row in violations),
            f"expected rule R7: {violations}",
        )

    def test_layer_derived_from_calibration_status_is_unreviewed(self):
        source = '''
def resolve_model_layer(row):
    status = row["functional_calibration_status"]
    if status == "calibrated_candidate_model":
        return "calibrated_candidate_model"
    return "discovery_hmm_uncalibrated"
'''
        findings = _findings(source)
        rows = _by_signal(findings, "model_layer_derived_from_calibration_status")
        self.assertTrue(rows, f"the reverse derivation must be named: {findings}")
        self.assertTrue(all(row["rule_cited"] == "R7_layer_status_independent" for row in rows))
        self.assertTrue(all(row["violation"] == "unreviewed" for row in rows))
        self.assertEqual(
            _violations(findings), [],
            f"reverse derivation is a review item, not a proven violation: {findings}",
        )


MODEL_STATUS_SCORING_SOURCE = '''
def select_scoring_profiles(manifest):
    trained = [row for row in manifest if row.get("model_status") == "trained"]
    return [row["profile_id"] for row in trained]
'''


MODEL_STATUS_REPORT_SOURCE = '''
from collections import Counter


def status_counts(manifest):
    counts = Counter(row.get("model_status", "") for row in manifest)
    return {"model_status_counts": dict(counts), "trained_count": counts.get("trained", 0)}
'''


class ModelStatusProxyTests(unittest.TestCase):
    def test_model_status_trained_criterion_is_unreviewed_or_a_violation(self):
        findings = _findings(MODEL_STATUS_SCORING_SOURCE)
        rows = _by_signal(findings, "model_status_used_as_layer_criterion")
        self.assertTrue(rows, f"the legacy-column criterion must be named: {findings}")
        self.assertTrue(
            all(row["violation"] in {"true", "unreviewed"} for row in rows),
            f"never silently compliant: {rows}",
        )
        self.assertTrue(
            all("model_layer" in row["note"] for row in rows),
            f"the note must point at the model_layer criterion: {rows}",
        )

    def test_pre_f3_defect_shape_is_a_violation(self):
        """The shape F3 is fixing: discovery-layer profiles read as ``trained`` and scored."""
        source = '''
def run(manifest):
    trained = [row for row in manifest if row.get("model_status") == "trained"]
    score_rows = []
    for row in trained:
        command = ["hmmsearch", "--tblout", "out.tblout", row["hmm"], "candidates.faa"]
        run_command(command)
        score_rows.append({"accession": row["profile_id"]})
    return score_rows
'''
        findings = _findings(source)
        violations = _violations(findings)
        self.assertTrue(
            violations,
            f"selecting profiles for scoring by model_status must be flagged: {findings}",
        )
        self.assertTrue(
            any(row["rule_cited"] == "R1_recall_only_only" for row in violations),
            f"expected rule R1: {violations}",
        )

    def test_model_status_report_count_is_not_a_violation(self):
        findings = _findings(MODEL_STATUS_REPORT_SOURCE)
        self.assertEqual(_violations(findings), [], f"counting statuses reports only: {findings}")


# ---------------------------------------------------------------------------
# Uncertainty, schema, CLI
# ---------------------------------------------------------------------------

class UncertaintyTests(unittest.TestCase):
    def test_unclassifiable_discovery_construct_is_unreviewed(self):
        source = '''
def dispatch(row, ctx):
    if ctx["discovery_mode"]:
        return transform(row)
    return row
'''
        findings = _findings(source)
        unreviewed = [row for row in findings if row["violation"] == "unreviewed"]
        self.assertTrue(unreviewed, f"an undecidable construct must be unreviewed: {findings}")
        self.assertTrue(any(row["use_class"] == "unknown" for row in unreviewed))

    def test_findings_are_never_all_false(self):
        corpus = [DISCOVERY_FILTER_SOURCE, RECALL_ONLY_SOURCE, REGISTRY_WRITE_SOURCE]
        values = set()
        for source in corpus:
            values.update(row["violation"] for row in _findings(source))
        self.assertTrue(values - {"false"}, f"the audit must not answer false to everything: {values}")

    def test_syntax_error_is_reported_as_unreviewed(self):
        findings = _findings("def broken(:\n    pass\n")
        self.assertTrue(findings, "unparseable source must not be silently skipped")
        self.assertTrue(all(row["violation"] == "unreviewed" for row in findings))
        self.assertTrue(all(row["use_class"] == "unknown" for row in findings))


class SchemaTests(unittest.TestCase):
    def test_findings_have_the_required_columns_and_vocabulary(self):
        findings = _findings(DISCOVERY_FILTER_SOURCE + REGISTRY_WRITE_SOURCE)
        self.assertTrue(findings)
        for row in findings:
            self.assertEqual(set(row), set(REQUIRED_COLUMNS), row)
            self.assertIn(row["evidence_kind"], {"ast", "text"})
            self.assertIn(row["use_class"], USE_CLASSES)
            self.assertIn(row["violation"], VIOLATION_VALUES)
            self.assertIsInstance(row["line"], int)
            self.assertGreaterEqual(row["line"], 1)
            self.assertTrue(row["rule_cited"])
            self.assertTrue(row["note"])

    def test_findings_are_returned_in_a_deterministic_order(self):
        module = _load_module()
        first = module.classify_usage(DISCOVERY_FILTER_SOURCE, "z.py")
        second = module.classify_usage(DISCOVERY_FILTER_SOURCE, "z.py")
        self.assertEqual(first, second)
        keys = [(row["script"], row["line"], row["signal"], row["use_class"]) for row in first]
        self.assertEqual(keys, sorted(keys))


class CliTests(unittest.TestCase):
    def _write_script(self, root, name, source):
        path = root / name
        path.write_text(source, encoding="utf-8")
        return path

    def test_cli_writes_the_audit_table_and_summary(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            scripts = root / "scripts"
            scripts.mkdir()
            self._write_script(scripts, "clean.py", RECALL_ONLY_SOURCE)
            self._write_script(scripts, "filtering.py", DISCOVERY_FILTER_SOURCE)
            out = root / "audit"
            code, stdout, _ = _call_main(
                module,
                ["--scripts-dir", str(scripts), "--out-dir", str(out), "--exclude", ""],
            )
            self.assertEqual(code, 1, stdout)
            table = out / "discovery_layer_usage_audit.tsv"
            summary_path = out / "discovery_layer_usage_audit_summary.json"
            self.assertTrue(table.is_file())
            self.assertTrue(summary_path.is_file())
            header = table.read_text(encoding="utf-8").splitlines()[0].split("\t")
            self.assertEqual(header, list(REQUIRED_COLUMNS))
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["scripts_scanned"], 2)
            self.assertGreaterEqual(summary["violations"], 1)
            self.assertIn("unreviewed", summary)
            self.assertEqual(json.loads(stdout)["violations"], summary["violations"])

    def test_cli_exits_zero_when_no_violation_is_found(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            scripts = root / "scripts"
            scripts.mkdir()
            self._write_script(scripts, "clean.py", RECALL_ONLY_SOURCE)
            code, _, _ = _call_main(
                module,
                ["--scripts-dir", str(scripts), "--out-dir", str(root / "audit"), "--exclude", ""],
            )
            self.assertEqual(code, 0)

    def test_cli_refuses_to_write_into_a_non_empty_output_directory(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            scripts = root / "scripts"
            scripts.mkdir()
            self._write_script(scripts, "filtering.py", DISCOVERY_FILTER_SOURCE)
            out = root / "audit"
            out.mkdir()
            sentinel = out / "keep.txt"
            sentinel.write_text("frozen evidence\n", encoding="utf-8")
            code, _, stderr = _call_main(
                module,
                ["--scripts-dir", str(scripts), "--out-dir", str(out), "--exclude", ""],
            )
            self.assertEqual(code, 2)
            self.assertIn("non-empty", stderr)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "frozen evidence\n")
            self.assertFalse((out / "discovery_layer_usage_audit.tsv").exists())

    def test_cli_rejects_a_missing_scripts_directory(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            code, _, stderr = _call_main(
                module,
                ["--scripts-dir", str(root / "absent"), "--out-dir", str(root / "audit")],
            )
            self.assertEqual(code, 2)
            self.assertIn("--scripts-dir", stderr)

    def test_cli_rejects_a_missing_registry_file(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            scripts = root / "scripts"
            scripts.mkdir()
            self._write_script(scripts, "clean.py", RECALL_ONLY_SOURCE)
            code, _, stderr = _call_main(
                module,
                ["--scripts-dir", str(scripts), "--out-dir", str(root / "audit"),
                 "--registry", str(root / "absent.tsv")],
            )
            self.assertEqual(code, 2)
            self.assertIn("--registry", stderr)


# ---------------------------------------------------------------------------
# Registry absence check (R3) and real-repository invariants
# ---------------------------------------------------------------------------

class RegistryAbsenceTests(unittest.TestCase):
    def test_registry_absence_check_detects_a_discovery_layer_entry(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            registry = Path(temporary_directory) / "formal_scan_models.tsv"
            registry.write_text(
                "model\thmm_source\tthreshold\n"
                "ePhaZ_curated_core\tsplit_run\te-5\n"
                "cys_discovery_uncalibrated\tdiscovery_run\te-5\n",
                encoding="utf-8",
            )
            report = module.check_registry_absence(registry)
            self.assertFalse(report["discovery_layer_absent"])
            self.assertEqual(report["discovery_layer_models"], ["cys_discovery_uncalibrated"])
            self.assertEqual(len(report["sha256"]), 64)

    def test_registry_absence_check_accepts_the_ephaZ_broad_discovery_model(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            registry = Path(temporary_directory) / "formal_scan_models.tsv"
            registry.write_text(
                "model\thmm_source\tthreshold\n"
                "ePhaZ_curated_core\tsplit_run\te-5\n"
                "ePhaZ_broad_discovery\tsplit_run\te-5\n",
                encoding="utf-8",
            )
            report = module.check_registry_absence(registry)
            self.assertTrue(report["discovery_layer_absent"])
            self.assertEqual(report["discovery_layer_models"], [])


class RealRepositoryTests(unittest.TestCase):
    """Stable governance pins. They never pin the current defect inventory."""

    def test_the_frozen_registry_does_not_list_the_discovery_layer(self):
        module = _load_module()
        report = module.check_registry_absence(REGISTRY_PATH)
        self.assertTrue(
            report["discovery_layer_absent"],
            f"discovery layer found in formal_scan_models.tsv: {report['discovery_layer_models']}",
        )
        self.assertEqual(len(report["sha256"]), 64)

    def test_real_repository_scan_is_schema_valid_and_not_vacuous(self):
        module = _load_module()
        findings = module.classify_repository(
            REPO_ROOT / "pipeline" / "scripts",
            exclude={SCRIPT_PATH.name},
        )
        self.assertTrue(findings, "the real repository scan must report something")
        for row in findings:
            self.assertEqual(set(row), set(REQUIRED_COLUMNS), row)
            self.assertIn(row["use_class"], USE_CLASSES)
            self.assertIn(row["violation"], VIOLATION_VALUES)
            self.assertTrue(row["script"].startswith("pipeline/scripts/"), row["script"])

    def test_the_audit_declares_its_own_limits(self):
        module = _load_module()
        limits = " ".join(module.AUDIT_LIMITS).lower()
        for needle in ("dynamic dispatch", "data-dependent", "no execution"):
            self.assertIn(needle, limits, f"missing stated limit: {needle}")


if __name__ == "__main__":
    unittest.main()
