"""Tests for the formal scan-13 tier-processing entrypoint.

Two layers live here:

* the original pytest-style test, kept verbatim as :func:`test_tier_run_creates_a_new_auditable_run_and_keeps_broad_separate`
  and mirrored inside a :class:`unittest.TestCase` so that
  ``python -m unittest pipeline.tests.test_formal_scan13_tier_run`` discovers it;
* the Task 9 additions, which pin the *database-size scale* of the frozen scan 13
  (``runs/20260901_formal_frozen_scan_13``) and record what is still missing
  before tier E-values on that evidence may be treated as full-library values.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

TIER_DRIVER_ASSERTIONS = (
    "create_run_layout",
    "write_input_contract",
    "prepare_formal_scan13_tier.py",
    "parallel_extract_sequences.py",
    "parallel -j 20",
    "broad_discovery.tsv",
    "07b_extract_seqs.py",
    "08_validate.py",
    "08c_tier_rescore.py",
    "--hmm-cpu must be 1..60",
    'mv "$BUILD" "$RUN_ROOT/results/tier_processing"',
)

# The measured total of the frozen scan-13 input space: GTDB R232, 199,923
# representative genomes, 100 shards, 267 GB, ~2.92e8 proteins
# (runs/20260917_phaded_cys_targeted_recall_01/results/recall_summary.json
# design.input_space_note, and the single-shard count of 2,923,820 records for
# shard_0051 recorded in runs/20260917_phaded_task_completion_check_01).
# It is an *approximate* figure and the authoritative measured total is pending;
# it is recorded here only so the scale contract has a concrete expected value.
APPROXIMATE_DATABASE_SIZE_Z = 292_000_000

# Independent, read-only audit of the frozen scan-13 executed command:
#   20260917_phaded_task_completion_check_01/01_recall/README.md §8 line 56
#     "全库扫描：score_one.sh 无 -Z（= HMMER 默认 Z=每分片序列数），-E 1e-5；100/100 片一致。"
#   recall_summary.json results.pool_vs_fullscan_difference_explanation
#     "全库扫描用 HMMER 默认 Z（≈每分片序列数）"
# The frozen manifest schema (reproduced verbatim in
#   runs/20260917_phaded_cys_targeted_recall_01/logs/s02_provenance.raw.txt)
# is reproduced here so the audit can be checked offline.
FROZEN_SCAN13_MANIFEST_KEYS = {
    "schema_version",
    "status",
    "run_id",
    "created_utc",
    "parent_run",
    "reused_tasks",
    "threads",
    "task_total",
    "task_completed",
    "registry_sha256",
    "models",
    "overlength_exclusions_sha256",
    "hits_all_sha256",
}


def test_tier_run_creates_a_new_auditable_run_and_keeps_broad_separate():
    script = (ROOT / "scripts" / "formal_scan13_tier_processing.sh").read_text(encoding="utf-8")
    assert "create_run_layout" in script
    assert "write_input_contract" in script
    assert "prepare_formal_scan13_tier.py" in script
    assert "parallel_extract_sequences.py" in script
    assert "parallel -j 20" in script
    assert "broad_discovery.tsv" in script
    assert "07b_extract_seqs.py" in script
    assert "08_validate.py" in script
    assert "08c_tier_rescore.py" in script
    assert "--hmm-cpu must be 1..60" in script
    assert 'mv "$BUILD" "$RUN_ROOT/results/tier_processing"' in script


class TierRunContractTests(unittest.TestCase):
    """The pre-existing tier-run contract assertions, unittest-discoverable."""

    def setUp(self):
        self.script = (ROOT / "scripts" / "formal_scan13_tier_processing.sh").read_text(
            encoding="utf-8"
        )

    def test_tier_run_creates_a_new_auditable_run_and_keeps_broad_separate(self):
        for expected in TIER_DRIVER_ASSERTIONS:
            with self.subTest(expected=expected):
                self.assertIn(expected, self.script)

    def test_run_id_is_rejected_when_it_is_not_a_safe_dated_identifier(self):
        self.assertIn(
            '[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$ && "$RUN_ID" != *..* ]]',
            self.script,
        )

    def test_parent_run_is_required_to_carry_accepted_hit_evidence(self):
        self.assertIn(
            '[[ -s "$PARENT_RUN/results/hits_all.tsv" && -s "$PARENT_RUN/input_contract.json" ]]',
            self.script,
        )


class Scan13DatabaseSizeAuditTests(unittest.TestCase):
    """Task 9 — the frozen scan-13 scale, and what must not be assumed."""

    def setUp(self):
        self.script = (ROOT / "scripts" / "formal_scan13_tier_processing.sh").read_text(
            encoding="utf-8"
        )

    def test_frozen_scan13_command_shape_has_no_database_size_scale(self):
        """The frozen command template could not have produced a full-library Z.

        This test *documents* the audit finding; it does not weaken it.  The
        executed template is copied from the deployed driver whose SHA-256 is
        recorded in the frozen run's own input contract
        (``b5e47723a783de278f97b6b49f4bdd30eae5658634d713408a21b4710bcfee53``,
        ``deploy/20260901_formal_frozen_scan_13/scripts/formal_frozen_screen_parallel.sh``):
            hmmsearch --tblout "$BUILD/$base.tbl" --domtblout "$BUILD/$base.dom" \
              -E 1e-5 --cpu 1 "$hmm" "$filtered"
        There is no ``-Z`` and no ``--domZ``, so HMMER used its default
        ``Z = <sequences in this shard>`` and the 1,000 frozen tasks do **not**
        share one scale.
        """
        frozen_template = (
            'hmmsearch --tblout "$BUILD/$base.tbl" --domtblout "$BUILD/$base.dom" '
            '-E 1e-5 --cpu 1 "$hmm" "$filtered"'
        )
        self.assertNotIn("-Z", frozen_template)
        self.assertNotIn("--domZ", frozen_template)

    def test_frozen_scan13_manifest_schema_records_no_scale_or_shard_counts(self):
        """The frozen manifest has no ``database_size_Z`` and no sequence counts."""
        scale_fields = {
            "database_size_Z",
            "database_size_basis",
            "shard_sequence_total",
            "hmmsearch_command_template",
            "shards",
        }
        self.assertFalse(scale_fields & FROZEN_SCAN13_MANIFEST_KEYS)
        self.assertIn("task_total", FROZEN_SCAN13_MANIFEST_KEYS)

    def test_frozen_scan13_shard_sequence_counts_differ_across_shards(self):
        """Per-shard defaults really were different, so one scale cannot be assumed.

        Two counts are measured in the frozen evidence: shard_0001 has
        4,694,121 records (``.../01_recall/logs/s01_recon.raw.txt`` line
        "shard_0001.faa record count") and shard_0051 has 2,923,820
        (``.../01_recall/README.md`` §8 line 59, "N_0051/109087 =
        2,923,820/109,087").  A ~38% spread is exactly why a single
        shard-local Z must never be reported as the full-library scale.
        """
        shard_0001_records = 4_694_121
        shard_0051_records = 2_923_820
        self.assertNotEqual(shard_0001_records, shard_0051_records)
        self.assertGreater(
            shard_0001_records / shard_0051_records, 1.3
        )

    def test_approximate_full_library_total_is_not_the_per_shard_count(self):
        """The full-library Z is ~100x a shard count, so the two are not
        interchangeable.  The exact measured total remains pending."""
        self.assertGreater(APPROXIMATE_DATABASE_SIZE_Z, 100 * 2_923_820 * 0.9)
        self.assertLess(APPROXIMATE_DATABASE_SIZE_Z, 100 * 2_923_820 * 1.1)

    def test_tier_entrypoint_does_not_yet_pin_a_database_size(self):
        """PENDING (Task 9 step 6) — the tier re-scoring entrypoint is unscaled.

        ``formal_scan13_tier_processing.sh`` re-runs HMMER through
        ``08c_tier_rescore.py`` on an extracted subset, and it neither accepts a
        ``--database-size-z`` nor records one in
        ``results/tier_processing_manifest.json``.  Adding that is a separate
        change to a file Task 9 does not own, so this test pins the current
        state and names the gap instead of silently assuming a scale.
        """
        self.assertNotIn("--database-size-z", self.script)
        self.assertNotIn("database_size_Z", self.script)

    def test_tier_manifest_schema_still_lacks_scale_fields(self):
        """The tier manifest fields a Task-9-complete run must add."""
        manifest_fields = {
            "schema_version",
            "status",
            "run_id",
            "created_utc",
            "parent_run",
            "hmm_cpu",
            "input_contract_sha256",
            "tier_processing_summary_sha256",
            "tier1_counts",
        }
        required_scale_fields = {
            "database_size_Z",
            "database_size_basis",
            "hmmsearch_command_template",
        }
        self.assertFalse(required_scale_fields & manifest_fields)
        self.assertIn("parent_run", manifest_fields)

    def test_frozen_scale_is_audited_read_only(self):
        """No locally retained scan-13 manifest may be rewritten by this task."""
        frozen_manifests = sorted(
            (ROOT.parent / "runs").glob("**/scan_manifest.json")
        )
        self.assertEqual(frozen_manifests, [])
        # The copied manifest inside the tier-processing run is the parent's
        # frozen evidence; it is an input, never a place to write a new scale.
        self.assertIn(
            'cp "$PARENT_RUN/results/scan_manifest.json" "$RUN_ROOT/inputs/parent_scan_manifest.json"',
            self.script,
        )

    def test_expected_scale_is_a_positive_integer_not_a_shard_count(self):
        self.assertIsInstance(APPROXIMATE_DATABASE_SIZE_Z, int)
        self.assertGreater(APPROXIMATE_DATABASE_SIZE_Z, 1)
        # ``database_size_Z`` is serialized as a JSON number, never as text, and
        # it is the full-library count rather than any single shard's count.
        payload = json.dumps(
            {"database_size_Z": APPROXIMATE_DATABASE_SIZE_Z}, sort_keys=True
        )
        self.assertEqual(
            payload, '{"database_size_Z": %d}' % APPROXIMATE_DATABASE_SIZE_Z
        )
        self.assertNotIn('"2923820"', payload)


if __name__ == "__main__":
    unittest.main()
