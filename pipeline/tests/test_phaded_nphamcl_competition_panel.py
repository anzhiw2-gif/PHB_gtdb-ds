"""Tests for the nPHAMCL competition panel (P4) builder and scorer.

Synthetic fixtures only; the frozen panel evidence and the real BLAST run live
outside these tests.
"""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import analyze_phaded_nphamcl_competition as scorer  # noqa: E402
import build_phaded_nphamcl_competition_panel as builder  # noqa: E402


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


class FastaTests(unittest.TestCase):
    def test_round_trip_and_duplicate_header_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = write(tmp / "a.faa", ">x\nACDE\n>y\nFGHI\n")
            self.assertEqual(builder.read_fasta(path), {"x": "ACDE", "y": "FGHI"})
            dup = write(tmp / "dup.faa", ">x\nAA\n>x\nBB\n")
            with self.assertRaisesRegex(ValueError, "duplicate FASTA header"):
                builder.read_fasta(dup)

    def test_sequence_before_header_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = write(Path(temporary) / "bad.faa", "ACDE\n>x\nAA\n")
            with self.assertRaisesRegex(ValueError, "before the first header"):
                builder.read_fasta(path)


class SelectionTests(unittest.TestCase):
    def test_only_the_nphamcl_superfamily_is_selected(self):
        rows = [
            {"accession": "a", "superfamily": builder.NPHAMCL_SUPERFAMILY},
            {"accession": "b", "superfamily": "extracellular dPHASCL type 1"},
        ]
        self.assertEqual(builder.select_nphamcl_accessions(rows), ["a"])

    def test_duplicate_accession_rejected(self):
        rows = [
            {"accession": "a", "superfamily": builder.NPHAMCL_SUPERFAMILY},
            {"accession": "a", "superfamily": builder.NPHAMCL_SUPERFAMILY},
        ]
        with self.assertRaisesRegex(ValueError, "not accession-unique"):
            builder.select_nphamcl_accessions(rows)


class PanelMeasurementTests(unittest.TestCase):
    def _rows(self, accepted: list[str], rejected: list[str]) -> list[dict[str, str]]:
        return [
            {"family": "hfam_4", "challenge": name, "evalue": "5e-12", "rejected": "no"}
            for name in accepted
        ] + [
            {"family": "hfam_4", "challenge": name, "evalue": "999.0", "rejected": "yes"}
            for name in rejected
        ]

    def test_panel_refuses_when_measurement_contradicts_the_declaration(self):
        rows = self._rows(["P24640"], ["P25275"])
        with self.assertRaisesRegex(ValueError, "disagrees with the declared hard competitors"):
            builder.verify_panel_against_measurement(rows)

    def test_panel_accepts_the_measured_three(self):
        rows = self._rows(["P24640", "Q02104", "Q88N36"], ["P25275", "P08658"])
        hard, controls, rejected = builder.verify_panel_against_measurement(rows)
        self.assertEqual(hard, ["P24640", "Q02104", "Q88N36"])
        self.assertEqual(controls, ["P08658", "P25275"])
        self.assertEqual(rejected, controls)

    def test_missing_family_raises(self):
        with self.assertRaisesRegex(ValueError, "no rows for hfam_4"):
            builder.verify_panel_against_measurement(
                [{"family": "hfam_52", "challenge": "x", "rejected": "yes"}],
            )

    def test_build_panel_fails_closed_on_a_missing_member(self):
        with self.assertRaisesRegex(ValueError, "no sequence in the challenge FASTA"):
            builder.build_panel({"P24640": "AA"}, {"8YNV_A": "CC"}, ["P24640", "Q02104"], [])

    def test_build_panel_requires_exactly_one_anchor(self):
        with self.assertRaisesRegex(ValueError, "exactly one record"):
            builder.build_panel({"P24640": "AA"}, {"a": "CC", "b": "DD"}, ["P24640"], [])


class ScoringTests(unittest.TestCase):
    def _blast(self, tmp: Path, rows: list[tuple[str, str, str]]) -> Path:
        body = "".join(
            f"{q}\t{s}\t30.0\t200\t{e}\t100.0\n" for q, s, e in rows
        )
        return write(tmp / "blast.tsv", body)

    def test_competitor_winning_makes_the_candidate_competed(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [("q1", "8YNV_A", "1e-20"), ("q1", "P24640", "1e-30")])
            records, summary = scorer.classify(
                scorer.read_blast(path), competitor_ids=["P24640"], control_ids=[],
            )
            self.assertEqual(records[0]["disposition"], scorer.DISPOSITION_COMPETED)
            self.assertEqual(summary["disposition_counts"][scorer.DISPOSITION_COMPETED], 1)

    def test_anchor_winning_supports_the_candidate(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [("q1", "8YNV_A", "1e-30"), ("q1", "P24640", "1e-10")])
            records, _summary = scorer.classify(
                scorer.read_blast(path), competitor_ids=["P24640"], control_ids=[],
            )
            self.assertEqual(records[0]["disposition"], scorer.DISPOSITION_SUPPORTED)
            self.assertEqual(records[0]["competitor_over_anchor_evalue_ratio"], "1e+20")

    def test_competitor_only_hit_is_competed_and_anchor_only_is_supported(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [("q1", "Q02104", "1e-9"), ("q2", "8YNV_A", "1e-9")])
            records, _ = scorer.classify(
                scorer.read_blast(path), competitor_ids=["Q02104"], control_ids=[],
            )
            by_query = {row["accession"]: row for row in records}
            self.assertEqual(by_query["q1"]["disposition"], scorer.DISPOSITION_COMPETED)
            self.assertEqual(by_query["q1"]["anchor_evalue"], "pending")
            self.assertEqual(by_query["q2"]["disposition"], scorer.DISPOSITION_SUPPORTED)

    def test_no_hit_to_any_panel_member_is_unresolved_not_excluded(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [("q1", "P25275", "999.0")])
            records, summary = scorer.classify(
                scorer.read_blast(path), competitor_ids=["P24640"], control_ids=["P25275"],
            )
            self.assertEqual(records[0]["disposition"], scorer.DISPOSITION_UNRESOLVED)
            self.assertIn("NOT", summary["interpretation"])
            self.assertIn("function_unresolved", summary["interpretation"])

    def test_malformed_and_non_numeric_rows_raise(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            short = write(tmp / "short.tsv", "q1\t8YNV_A\t30.0\n")
            with self.assertRaisesRegex(ValueError, "malformed BLAST row"):
                scorer.read_blast(short)
            bad = write(tmp / "bad.tsv", "q1\t8YNV_A\t30.0\t200\tnotanumber\t100.0\n")
            with self.assertRaisesRegex(ValueError, "non-numeric BLAST score"):
                scorer.read_blast(bad)

    def test_summary_states_the_sequence_only_boundary(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [("q1", "8YNV_A", "1e-9")])
            _records, summary = scorer.classify(scorer.read_blast(path), competitor_ids=[])
            self.assertIn("Sequence-level competition only", summary["boundary"])
            self.assertIn("Foldseek", summary["boundary"])

    def test_a_zero_evalue_does_not_divide_by_zero(self):
        # Strong hits underflow to E = 0.0, and the ratio must not blow up. The
        # three cases are named rather than turned into a number that looks like data.
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [
                ("both_zero", "8YNV_A", "0.0"), ("both_zero", "P24640", "0.0"),
                ("anchor_zero", "8YNV_A", "0.0"), ("anchor_zero", "P24640", "1e-5"),
                ("competitor_zero", "8YNV_A", "1e-5"), ("competitor_zero", "P24640", "0.0"),
            ])
            records, _summary = scorer.classify(
                scorer.read_blast(path), competitor_ids=["P24640"], control_ids=[],
            )
            by_query = {row["accession"]: row for row in records}
            self.assertEqual(by_query["both_zero"]["competitor_over_anchor_evalue_ratio"],
                             "both_below_precision")
            self.assertEqual(by_query["anchor_zero"]["competitor_over_anchor_evalue_ratio"], "inf")
            self.assertEqual(by_query["competitor_zero"]["competitor_over_anchor_evalue_ratio"], "0")
            self.assertEqual(by_query["anchor_zero"]["disposition"], scorer.DISPOSITION_SUPPORTED)
            self.assertEqual(by_query["competitor_zero"]["disposition"],
                             scorer.DISPOSITION_COMPETED)

    def test_disposition_labels_can_be_renamed_for_another_pool(self):
        # The same panel logic scores the with-lipase deferred pool, whose members
        # are not nPHAMCL candidates; the labels must not claim otherwise.
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [("q1", "8YNV_A", "1e-30"), ("q1", "P24640", "1e-10"),
                                     ("q2", "P24640", "1e-9"), ("q3", "P08658", "999.0")])
            records, summary = scorer.classify(
                scorer.read_blast(path), competitor_ids=["P24640"], control_ids=["P08658"],
                supported_label="anchor_nearest", competed_label="competitor_nearest",
                unresolved_label="no_panel_hit",
            )
            by_query = {row["accession"]: row for row in records}
            self.assertEqual(by_query["q1"]["disposition"], "anchor_nearest")
            self.assertEqual(by_query["q2"]["disposition"], "competitor_nearest")
            self.assertEqual(by_query["q3"]["disposition"], "no_panel_hit")
            self.assertEqual(summary["disposition_counts"],
                             {"anchor_nearest": 1, "competitor_nearest": 1, "no_panel_hit": 1})
            self.assertAlmostEqual(summary["supported_share"], 1 / 3)
            self.assertAlmostEqual(summary["competed_share"], 1 / 3)

    def test_cli_writes_artifacts_and_refuses_a_non_empty_dir(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = self._blast(tmp, [("q1", "8YNV_A", "1e-9")])
            out = tmp / "out"
            self.assertEqual(
                scorer.main(["--blast", str(path), "--competitors", "P24640",
                             "--controls", "P25275", "--out-dir", str(out)]), 0,
            )
            payload = json.loads((out / "nphamcl_competition_summary.json").read_text("utf-8"))
            self.assertEqual(payload["candidates_scored"], 1)
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                scorer.main(["--blast", str(path), "--competitors", "P24640",
                             "--controls", "P25275", "--out-dir", str(out)])


class PanelBuilderCliTests(unittest.TestCase):
    def test_cli_builds_targets_and_panel(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            merged = write_tsv(tmp / "merged.tsv", ["accession", "superfamily"], [
                {"accession": "a", "superfamily": builder.NPHAMCL_SUPERFAMILY},
                {"accession": "b", "superfamily": "extracellular dPHASCL type 1"},
            ])
            union = write(tmp / "union.faa", ">a\nACDE\n>b\nFFFF\n")
            challenge = write(tmp / "challenge.faa", ">P24640\nAAAA\n>Q02104\nBBBB\n>Q88N36\nCCCC\n>P25275\nDDDD\n")
            validation = write_tsv(tmp / "cv.tsv", ["family", "challenge", "evalue", "rejected"], [
                {"family": "hfam_4", "challenge": "P24640", "evalue": "5.4e-12", "rejected": "no"},
                {"family": "hfam_4", "challenge": "Q02104", "evalue": "4.1e-11", "rejected": "no"},
                {"family": "hfam_4", "challenge": "Q88N36", "evalue": "6.6e-11", "rejected": "no"},
                {"family": "hfam_4", "challenge": "P25275", "evalue": "999.0", "rejected": "yes"},
            ])
            anchor = write(tmp / "8ynv.faa", ">8YNV_A\nEEEE\n")
            out = tmp / "out"
            self.assertEqual(
                builder.main(["--merged", str(merged), "--candidate-union", str(union),
                              "--challenge-fasta", str(challenge),
                              "--challenge-validation", str(validation),
                              "--anchor-fasta", str(anchor), "--out-dir", str(out),
                              "--expect-targets", "1"]), 0,
            )
            provenance = json.loads((out / "panel_provenance.json").read_text("utf-8"))
            self.assertEqual(provenance["target_count"], 1)
            self.assertEqual(provenance["panel_size"], 5)
            self.assertTrue(provenance["target_count_matches_declared"])
            self.assertIn("forced by frozen measurement", provenance["derivation"])
            panel = builder.read_fasta(out / "competition_panel.faa")
            self.assertEqual(sorted(panel), ["8YNV_A", "P24640", "P25275", "Q02104", "Q88N36"])


if __name__ == "__main__":
    unittest.main()
