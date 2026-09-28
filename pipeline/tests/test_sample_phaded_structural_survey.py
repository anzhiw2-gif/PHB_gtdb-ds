"""Tests for the P4 structural survey sampler.

Synthetic fixtures only; the real P4 score table is read by the run.
"""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import sample_phaded_structural_survey as module  # noqa: E402


SCORE_FIELDS = [
    "accession", "anchor_bitscore", "anchor_evalue", "competitor_subject",
    "competitor_bitscore", "competitor_evalue", "disposition",
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


def score_row(accession: str, margin: float) -> dict[str, str]:
    """Build a row whose competitor-minus-anchor margin is exactly ``margin``."""
    anchor = 50.0
    return {
        "accession": accession, "anchor_bitscore": f"{anchor:.1f}", "anchor_evalue": "1e-20",
        "competitor_subject": "Q88N36", "competitor_bitscore": f"{anchor + margin:.1f}",
        "competitor_evalue": "1e-20", "disposition": "nphamcl_like_supported",
    }


class SanitizeTests(unittest.TestCase):
    def test_stop_codon_and_whitespace_are_stripped(self):
        self.assertEqual(module.sanitize_sequence("a", "ACDE*"), "ACDE")
        self.assertEqual(module.sanitize_sequence("a", "AC DE\nFG"), "ACDEFG")

    def test_ambiguous_is_kept_and_unsupported_fails_closed(self):
        self.assertEqual(module.sanitize_sequence("a", "ACDXE"), "ACDXE")
        with self.assertRaisesRegex(ValueError, "unsupported residue"):
            module.sanitize_sequence("a", "ACDE!")


class ScoringTests(unittest.TestCase):
    def test_margin_is_competitor_minus_anchor_and_order_is_ascending(self):
        rows = [score_row("hi", 3.0), score_row("lo", -12.0)]
        scored = module.scored_candidates(rows)
        self.assertEqual([row["accession"] for row in scored], ["lo", "hi"])
        self.assertAlmostEqual(scored[0]["margin"], -12.0)

    def test_unscored_rows_are_skipped(self):
        rows = [score_row("a", 1.0), {"accession": "b", "anchor_bitscore": "pending",
                                      "competitor_bitscore": "pending"}]
        scored = module.scored_candidates(rows)
        self.assertEqual([row["accession"] for row in scored], ["a"])

    def test_no_scored_rows_is_refused(self):
        with self.assertRaisesRegex(ValueError, "no candidate carries a usable panel margin"):
            module.scored_candidates([{"accession": "b", "anchor_bitscore": "pending",
                                       "competitor_bitscore": "pending"}])

    def test_row_without_accession_is_refused(self):
        with self.assertRaisesRegex(ValueError, "score row without an accession"):
            module.scored_candidates([{"anchor_bitscore": "1", "competitor_bitscore": "2"}])


class StratumTests(unittest.TestCase):
    def test_deciles_partition_every_candidate(self):
        rows = [score_row(f"c{i:02d}", float(i)) for i in range(25)]
        strata = module.decile_strata(module.scored_candidates(rows))
        self.assertEqual(len(strata), module.DECILES)
        self.assertEqual(sum(len(stratum) for stratum in strata), 25)
        self.assertEqual(strata[0][0]["accession"], "c00")
        self.assertEqual(strata[-1][-1]["accession"], "c24")

    def test_too_few_candidates_is_refused(self):
        rows = [score_row(f"c{i}", float(i)) for i in range(5)]
        with self.assertRaisesRegex(ValueError, "to form deciles"):
            module.decile_strata(module.scored_candidates(rows))

    def test_rule_of_one_is_deterministic_and_covers_every_decile(self):
        rows = [score_row(f"c{i:02d}", float(i)) for i in range(100)]
        strata = module.decile_strata(module.scored_candidates(rows))
        first = module.draw(strata, 3, 7)
        second = module.draw(strata, 3, 7)
        self.assertEqual([row["accession"] for row in first],
                         [row["accession"] for row in second])
        self.assertEqual(len(first), 30)
        self.assertEqual({row["margin_decile"] for row in first}, set(range(module.DECILES)))
        self.assertEqual(len({row["accession"] for row in first}), 30)

    def test_draw_is_within_stratum(self):
        rows = [score_row(f"c{i:02d}", float(i)) for i in range(100)]
        strata = module.decile_strata(module.scored_candidates(rows))
        drawn = module.draw(strata, 2, 11)
        for row in drawn:
            stratum = strata[row["margin_decile"]]
            self.assertIn(row["accession"], [member["accession"] for member in stratum])

    def test_asking_for_more_than_a_decile_holds_is_refused(self):
        rows = [score_row(f"c{i}", float(i)) for i in range(10)]
        strata = module.decile_strata(module.scored_candidates(rows))
        with self.assertRaisesRegex(ValueError, "fewer than the"):
            module.draw(strata, 5, 1)

    def test_zero_per_decile_is_refused(self):
        rows = [score_row(f"c{i}", float(i)) for i in range(10)]
        strata = module.decile_strata(module.scored_candidates(rows))
        with self.assertRaisesRegex(ValueError, "at least 1"):
            module.draw(strata, 0, 1)


class CliTests(unittest.TestCase):
    def _fixture(self, tmp: Path, n: int = 100, family: str = module.NPHAMCL_SUPERFAMILY):
        scores = write_tsv(tmp / "scores.tsv", SCORE_FIELDS,
                           [score_row(f"c{i:03d}", float(i) - 50.0) for i in range(n)])
        with (tmp / "union.faa").open("w", encoding="ascii", newline="\n") as handle:
            for i in range(n):
                handle.write(f">c{i:03d}\n" + "ACDEFGHIKL" * 30 + "*\n")
        merged = write_tsv(tmp / "merged.tsv", ["accession", "superfamily"],
                           [{"accession": f"c{i:03d}", "superfamily": family} for i in range(n)])
        return scores, merged, tmp / "union.faa"

    def test_cli_writes_a_sanitised_fasta_and_a_composition_table(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, merged, union = self._fixture(tmp)
            out = tmp / "out"
            self.assertEqual(
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(merged), "--out-dir", str(out),
                             "--per-decile", "3", "--seed", "5", "--expect-scored", "100"]), 0,
            )
            records = module.read_fasta(out / "survey_targets.faa")
            self.assertEqual(len(records), 30)
            self.assertFalse(any("*" in sequence for sequence in records.values()))
            summary = json.loads((out / "survey_summary.json").read_text("utf-8"))
            self.assertEqual(summary["sample_size"], 30)
            self.assertEqual(summary["seed"], 5)
            self.assertEqual(len(summary["stratum_sizes"]), 10)
            self.assertEqual(len(summary["draw_per_decile"]), 10)
            self.assertIn("proportional draw would barely touch", summary["why_deciles_not_proportional"])
            self.assertIn("TRANCHE", summary["scope"])
            self.assertIn("not a phenotype claim", summary["boundary"])
            rows = list(csv.DictReader((out / "survey_composition.tsv").open(encoding="utf-8"),
                                       delimiter="\t"))
            self.assertEqual(len(rows), 30)

    def test_a_changed_universe_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, merged, union = self._fixture(tmp)
            with self.assertRaisesRegex(ValueError, "does not match the declared"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(merged), "--out-dir", str(tmp / "out"),
                             "--expect-scored", "987"])

    def test_a_candidate_outside_the_superfamily_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, merged, union = self._fixture(tmp, family="extracellular dPHASCL type 1")
            with self.assertRaisesRegex(ValueError, "not 'intracellular nPHAMCL'"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(merged), "--out-dir", str(tmp / "out"),
                             "--expect-scored", "100"])

    def test_a_missing_sequence_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, merged, _union = self._fixture(tmp)
            short = tmp / "short.faa"
            short.write_text(">c000\nACDE\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "has no sequence in"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(short),
                             "--merged", str(merged), "--out-dir", str(tmp / "out"),
                             "--expect-scored", "100"])

    def test_a_pool_with_no_merged_table_skips_the_superfamily_check(self):
        # The deferred with-lipase pool's identity is fixed by the score table it came
        # from, so requiring a merged table there would be a check against nothing.
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, _merged, union = self._fixture(tmp)
            out = tmp / "out"
            self.assertEqual(
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--out-dir", str(out), "--per-decile", "3", "--seed", "5",
                             "--expect-scored", "100"]), 0,
            )
            summary = json.loads((out / "survey_summary.json").read_text("utf-8"))
            self.assertIn("skipped (no --merged table given", summary["superfamily_check"])

    def test_a_separate_sequence_source_is_used_and_reported(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, merged, union = self._fixture(tmp)
            other = tmp / "other.faa"
            with other.open("w", encoding="ascii", newline="\n") as handle:
                for i in range(100):
                    handle.write(f">c{i:03d}\n" + "MKLW" * 40 + "*\n")
            out = tmp / "out"
            self.assertEqual(
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--sequence-source", str(other), "--merged", str(merged),
                             "--out-dir", str(out), "--per-decile", "2", "--seed", "9",
                             "--expect-scored", "100"]), 0,
            )
            records = module.read_fasta(out / "survey_targets.faa")
            self.assertTrue(all(sequence.startswith("MKLW") for sequence in records.values()),
                            "sequences must come from --sequence-source")
            self.assertFalse(any("*" in sequence for sequence in records.values()))
            summary = json.loads((out / "survey_summary.json").read_text("utf-8"))
            self.assertIn("other.faa", summary["sequence_source"])

    def test_non_empty_output_dir_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, merged, union = self._fixture(tmp)
            out = tmp / "out"
            out.mkdir()
            (out / "existing").write_text("x", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(merged), "--out-dir", str(out)])


if __name__ == "__main__":
    unittest.main()
