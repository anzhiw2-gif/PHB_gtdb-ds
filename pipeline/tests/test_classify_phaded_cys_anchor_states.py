"""Tests for the P3 Cys anchor-state classification.

Synthetic fixtures only; the real frozen motif layer and candidate sequences are
read by the run.
"""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import classify_phaded_cys_anchor_states as module  # noqa: E402


EVIDENCE_FIELDS = [
    "accession", "pool", "n_cysteine", "cys_nucleophile_motif",
    "cys_nucleophile_motif_positions", "cys_nucleophile_motif_tripeptides",
    "evidence_basis", "catalytic_activity_verified", "independent_validation",
]
MERGED_FIELDS = ["accession", "superfamily", "nucleophile_type"]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


def evidence_row(accession: str, tier: str, n_cysteine: str = "3") -> dict[str, str]:
    return {
        "accession": accession, "pool": "pool_in", "n_cysteine": n_cysteine,
        "cys_nucleophile_motif": tier, "cys_nucleophile_motif_positions": "100",
        "cys_nucleophile_motif_tripeptides": "VCQ",
        "evidence_basis": "sequence_context_match_to_single_literature_anchor_phaZ1_C183_PMID16233560",
        "catalytic_activity_verified": "false", "independent_validation": "pending",
    }


class AnchorConstantTests(unittest.TestCase):
    def test_anchor_constants_match_the_frozen_calibration(self):
        self.assertEqual(module.ANCHOR_POSITION, 183)
        self.assertEqual(module.ANCHOR_CONTEXT, "VCQ")
        self.assertEqual(module.ANCHOR_PMID, "16233560")

    def test_caveat_forbids_reading_uncertain_as_a_negative(self):
        self.assertIn("C183S retains activity", module.CAVEAT)
        self.assertIn("never be read as evidence against", module.CAVEAT)


class StateRuleTests(unittest.TestCase):
    def test_pattern_tiers_become_pattern_present(self):
        for tier in module.TIER_PATTERN:
            state, reason = module.anchor_state(tier, 3, 400)
            self.assertEqual(state, module.STATE_PATTERN, tier)
            self.assertIn("anchor context reproduced", reason)

    def test_hydrophobic_only_is_a_substitution(self):
        state, reason = module.anchor_state("minus1_hydrophobic_only", 3, 400)
        self.assertEqual(state, module.STATE_SUBSTITUTION)
        self.assertIn("+1 is not Q", reason)

    def test_short_sequence_is_truncation_before_any_tier_check(self):
        state, reason = module.anchor_state("anchor_exact_vcq", 3, 182)
        self.assertEqual(state, module.STATE_TRUNCATION)
        self.assertIn("anchor cannot be present", reason)

    def test_full_length_with_no_cysteine_is_uncertain_not_substitution(self):
        state, reason = module.anchor_state("no_cysteine", 0, 400)
        self.assertEqual(state, module.STATE_UNCERTAIN)
        self.assertIn("neither a clean substitution call nor decidable", reason)

    def test_not_detected_is_uncertain(self):
        state, _reason = module.anchor_state("not_detected", 2, 400)
        self.assertEqual(state, module.STATE_UNCERTAIN)

    def test_unrecognised_tier_is_uncertain_and_says_so(self):
        state, reason = module.anchor_state("some_new_tier", 1, 400)
        self.assertEqual(state, module.STATE_UNCERTAIN)
        self.assertIn("unrecognised tier", reason)

    def test_exact_boundary_length_is_not_truncation(self):
        state, _reason = module.anchor_state("anchor_exact_vcq", 3, module.ANCHOR_POSITION)
        self.assertEqual(state, module.STATE_PATTERN)


class FastaLengthesTests(unittest.TestCase):
    def test_lengths_are_summed_over_wrapped_lines(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = tmp / "seqs.faa"
            path.write_text(">a\nACDEF\nGHIKL\n>b\nMNPQ\n", encoding="utf-8")
            self.assertEqual(module.read_fasta_lengths(path), {"a": 10, "b": 4})

    def test_duplicate_header_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "dup.faa"
            path.write_text(">a\nAA\n>a\nBB\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate FASTA header"):
                module.read_fasta_lengths(path)


class ClassifyAndCliTests(unittest.TestCase):
    def _fixture(self, tmp: Path):
        evidence = write_tsv(tmp / "evidence.tsv", EVIDENCE_FIELDS, [
            evidence_row("c1", "anchor_exact_vcq"),
            evidence_row("c2", "minus1_hydrophobic_only"),
            evidence_row("c3", "not_detected"),
            evidence_row("c4", "anchor_exact_vcq"),
            evidence_row("nc1", "anchor_exact_vcq"),   # only in the merged table as non-cys
        ])
        merged = write_tsv(tmp / "merged.tsv", MERGED_FIELDS, [
            {"accession": "c1", "superfamily": "x", "nucleophile_type": "cys"},
            {"accession": "c2", "superfamily": "x", "nucleophile_type": "cys"},
            {"accession": "c3", "superfamily": "x", "nucleophile_type": "cys"},
            {"accession": "c4", "superfamily": "x", "nucleophile_type": "cys"},
            {"accession": "nc1", "superfamily": "x", "nucleophile_type": "ser"},
        ])
        union = tmp / "union.faa"
        union.write_text(
            ">c1\n" + "A" * 400 + "\n>c2\n" + "A" * 400 + "\n>c3\n" + "A" * 400 + "\n"
            + ">c4\n" + "A" * 100 + "\n>nc1\n" + "A" * 400 + "\n",
            encoding="utf-8",
        )
        return evidence, merged, union

    def test_only_cys_candidates_are_scored(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            evidence, merged, union = self._fixture(tmp)
            lengths = module.read_fasta_lengths(union)
            targets = {"c1", "c2", "c3", "c4"}
            rows, summary = module.classify(module.read_rows(evidence), lengths, targets)
            self.assertEqual(len(rows), 4)
            self.assertNotIn("nc1", {row["accession"] for row in rows})
            self.assertEqual(summary["state_counts"][module.STATE_PATTERN], 1)
            self.assertEqual(summary["state_counts"][module.STATE_SUBSTITUTION], 1)
            self.assertEqual(summary["state_counts"][module.STATE_UNCERTAIN], 1)
            self.assertEqual(summary["state_counts"][module.STATE_TRUNCATION], 1)

    def test_missing_sequence_length_is_pending_and_uncertain(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            evidence, merged, union = self._fixture(tmp)
            rows, summary = module.classify(
                module.read_rows(evidence), {}, {"c1"},
            )
            self.assertEqual(rows[0]["anchor_state"], module.STATE_UNCERTAIN)
            self.assertEqual(rows[0]["length"], "pending")
            self.assertEqual(summary["accessions_missing_sequence_length"], ["c1"])

    def test_cli_writes_artifacts_and_reports_the_expectation(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            evidence, merged, union = self._fixture(tmp)
            out = tmp / "out"
            self.assertEqual(
                module.main(["--motif-evidence", str(evidence), "--candidate-union", str(union),
                             "--merged", str(merged), "--out-dir", str(out),
                             "--expect-candidates", "4"]), 0,
            )
            summary = json.loads((out / "cys_anchor_states_summary.json").read_text("utf-8"))
            self.assertEqual(summary["cys_candidates_from_merged"], 4)
            self.assertTrue(summary["expectation_matches"])
            self.assertIn("rule_order", summary)
            self.assertIn("NOT catalytic-activity verification", summary["caveat"])
            self.assertIn("deletes, demotes and promotes nothing", summary["boundary"])
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--motif-evidence", str(evidence), "--candidate-union", str(union),
                             "--merged", str(merged), "--out-dir", str(out)])

    def test_expectation_mismatch_is_reported_not_hidden(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            evidence, merged, union = self._fixture(tmp)
            out = tmp / "out"
            self.assertEqual(
                module.main(["--motif-evidence", str(evidence), "--candidate-union", str(union),
                             "--merged", str(merged), "--out-dir", str(out),
                             "--expect-candidates", "29974"]), 0,
            )
            summary = json.loads((out / "cys_anchor_states_summary.json").read_text("utf-8"))
            self.assertFalse(summary["expectation_matches"])
            self.assertEqual(summary["declared_expectation"], 29974)


if __name__ == "__main__":
    unittest.main()
