"""Tests for the P4 structural-pilot assembly.

Synthetic fixtures only; the real frozen tables are read by the run.
"""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import build_phaded_structural_pilot as module  # noqa: E402


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


SCORE_FIELDS = [
    "accession", "anchor_bitscore", "anchor_evalue", "competitor_subject",
    "competitor_bitscore", "competitor_evalue", "disposition",
]


def score_row(accession: str, anchor: str, competitor: str, **overrides) -> dict[str, str]:
    row = {
        "accession": accession, "anchor_bitscore": anchor, "anchor_evalue": "1e-20",
        "competitor_subject": "Q88N36", "competitor_bitscore": competitor,
        "competitor_evalue": "1e-20", "disposition": "nphamcl_like_supported",
    }
    row.update(overrides)
    return row


class SanitizeTests(unittest.TestCase):
    def test_stop_codon_and_whitespace_are_stripped(self):
        cleaned, ambiguous = module.sanitize_sequence("a", "ACDE*")
        self.assertEqual(cleaned, "ACDE")
        self.assertEqual(ambiguous, [])
        cleaned, _ = module.sanitize_sequence("a", "AC DE\nFG")
        self.assertEqual(cleaned, "ACDEFG")

    def test_ambiguous_residues_are_kept_and_reported(self):
        cleaned, ambiguous = module.sanitize_sequence("a", "ACDXE")
        self.assertEqual(cleaned, "ACDXE")
        self.assertEqual(ambiguous, ["X"])

    def test_an_unsupported_residue_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "unsupported residue"):
            module.sanitize_sequence("a", "ACDE!")

    def test_a_sequence_that_is_only_stop_codons_is_refused(self):
        with self.assertRaisesRegex(ValueError, "empty after sanitisation"):
            module.sanitize_sequence("a", "**")


class ExtremeSelectionTests(unittest.TestCase):
    def test_extremes_are_the_min_and_max_margin(self):
        rows = [
            score_row("anchor_side", "60.0", "40.0"),
            score_row("competitor_side", "40.0", "50.0"),
            score_row("middle", "50.0", "48.0"),
        ]
        anchor_favoured, competitor_favoured = module.pilot_extremes_from_scores(rows)
        self.assertEqual(anchor_favoured["accession"], "anchor_side")
        self.assertEqual(competitor_favoured["accession"], "competitor_side")
        self.assertEqual(anchor_favoured["margin"], "-20.0000")
        self.assertEqual(competitor_favoured["margin"], "10.0000")

    def test_unscored_rows_are_skipped(self):
        rows = [
            score_row("a", "60.0", "40.0"),
            score_row("b", "40.0", "50.0"),
            score_row("pending", "pending", "pending"),
        ]
        anchor_favoured, competitor_favoured = module.pilot_extremes_from_scores(rows)
        self.assertEqual({anchor_favoured["accession"], competitor_favoured["accession"]},
                         {"a", "b"})

    def test_single_scored_candidate_is_refused(self):
        with self.assertRaisesRegex(ValueError, "fewer than two scored candidates"):
            module.pilot_extremes_from_scores([score_row("only", "60.0", "40.0")])

    def test_row_without_accession_is_refused(self):
        rows = [score_row("a", "60.0", "40.0"), {"anchor_bitscore": "1", "competitor_bitscore": "2"}]
        with self.assertRaisesRegex(ValueError, "score row without an accession"):
            module.pilot_extremes_from_scores(rows)


class FastaTests(unittest.TestCase):
    def test_duplicate_header_and_orphan_sequence_are_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            with self.assertRaisesRegex(ValueError, "duplicate FASTA header"):
                module.read_fasta(write(tmp / "dup.faa", ">a\nAA\n>a\nBB\n"))
            with self.assertRaisesRegex(ValueError, "before the first header"):
                module.read_fasta(write(tmp / "orphan.faa", "AA\n>a\nBB\n"))


class AssemblyTests(unittest.TestCase):
    def _fixture(self, tmp: Path):
        scores = write_tsv(tmp / "scores.tsv", SCORE_FIELDS, [
            score_row("cand_anchor", "60.0", "40.0"),
            score_row("cand_competitor", "40.0", "50.0",
                      disposition="competed_by_measured_confounder"),
        ])
        union = write(tmp / "union.faa",
                      ">cand_anchor\n" + "A" * 300 + "\n>cand_competitor\n" + "B" * 150 + "\n")
        merged = write_tsv(tmp / "merged.tsv", ["accession", "superfamily"], [
            {"accession": "cand_anchor", "superfamily": module.NPHAMCL_SUPERFAMILY},
            {"accession": "cand_competitor", "superfamily": module.NPHAMCL_SUPERFAMILY},
        ])
        challenge = write(tmp / "challenge.faa", ">P24640\nCCCC\n>Q02104\nDDDD\n>Q88N36\nEEEE\n")
        return scores, union, merged, challenge

    def test_pilot_has_the_two_margin_extremes_plus_the_competitors(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, union, merged, challenge = self._fixture(tmp)
            out = tmp / "out"
            self.assertEqual(
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(merged), "--challenge-fasta", str(challenge),
                             "--out-dir", str(out),
                             "--competitors", "P24640", "Q02104", "Q88N36"]), 0,
            )
            composition = json.loads((out / "pilot_composition.json").read_text("utf-8"))
            self.assertEqual(composition["member_count"], 5)
            self.assertEqual(composition["members"]["cand_anchor"]["role"], "most_anchor_favoured")
            self.assertEqual(composition["members"]["cand_competitor"]["role"],
                             "most_competitor_favoured")
            self.assertEqual(composition["members"]["P24640"]["role"], "hard_competitor")
            self.assertIn("experimental PDB structure", composition["anchor_source"])
            self.assertIn("cannot separate anything", composition["why_these_candidates"])
            self.assertIn("boundary", composition)
            records = module.read_fasta(out / "pilot_prediction_targets.faa")
            self.assertEqual(len(records), 5)

    def test_a_pilot_candidate_outside_the_superfamily_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, union, _merged, challenge = self._fixture(tmp)
            wrong = write_tsv(tmp / "wrong.tsv", ["accession", "superfamily"], [
                {"accession": "cand_anchor", "superfamily": "extracellular dPHASCL type 1"},
                {"accession": "cand_competitor", "superfamily": module.NPHAMCL_SUPERFAMILY},
            ])
            with self.assertRaisesRegex(ValueError, "not 'intracellular nPHAMCL'"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(wrong), "--challenge-fasta", str(challenge),
                             "--out-dir", str(tmp / "out"), "--competitors", "P24640"])

    def test_missing_candidate_sequence_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, _union, merged, challenge = self._fixture(tmp)
            short = write(tmp / "short.faa", ">cand_anchor\n" + "A" * 300 + "\n")
            with self.assertRaisesRegex(ValueError, "no sequence in the candidate union"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(short),
                             "--merged", str(merged), "--challenge-fasta", str(challenge),
                             "--out-dir", str(tmp / "out"), "--competitors", "P24640"])

    def test_missing_competitor_sequence_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, union, merged, challenge = self._fixture(tmp)
            with self.assertRaisesRegex(ValueError, "no sequence in the challenge FASTA"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(merged), "--challenge-fasta", str(challenge),
                             "--out-dir", str(tmp / "out"), "--competitors", "NOT_THERE"])

    def test_non_empty_output_dir_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            scores, union, merged, challenge = self._fixture(tmp)
            out = tmp / "out"
            out.mkdir()
            (out / "existing").write_text("x", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--competition-scores", str(scores), "--candidate-union", str(union),
                             "--merged", str(merged), "--challenge-fasta", str(challenge),
                             "--out-dir", str(out), "--competitors", "P24640"])


if __name__ == "__main__":
    unittest.main()
