"""Tests for the adapter additivity verifier.

Synthetic fixtures only; the real 75 MB outputs are read by the run.
"""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import verify_phaded_catalog_input_is_additive as module  # noqa: E402


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    return path


class AdditivityTests(unittest.TestCase):
    def test_an_added_populated_column_is_additive(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["accession", "x"],
                                [{"accession": "p1", "x": "1"}, {"accession": "p2", "x": "2"}])
            revision = write_tsv(tmp / "b.tsv", ["accession", "x", "gap"],
                                 [{"accession": "p1", "x": "1", "gap": "1.3"},
                                  {"accession": "p2", "x": "2", "gap": "0.4"}])
            report = module.compare(shipped, revision)
            self.assertEqual(report["added_columns"], ["gap"])
            self.assertEqual(report["removed_columns"], [])
            self.assertEqual(report["mismatched_value_count"], 0)
            self.assertIn("ADDITIVE", report["verdict"])
            self.assertEqual(report["added_column_fill"]["gap"],
                             {"rows": 2, "non_empty": 2})

    def test_a_changed_value_is_not_additive(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["accession", "x"], [{"accession": "p1", "x": "1"}])
            revision = write_tsv(tmp / "b.tsv", ["accession", "x"], [{"accession": "p1", "x": "9"}])
            report = module.compare(shipped, revision)
            self.assertEqual(report["mismatched_value_count"], 1)
            self.assertIn("NOT ADDITIVE", report["verdict"])
            self.assertIn("x:", report["mismatch_examples"][0])

    def test_a_reordered_or_renamed_row_is_not_additive(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["accession", "x"],
                                [{"accession": "p1", "x": "1"}, {"accession": "p2", "x": "2"}])
            revision = write_tsv(tmp / "b.tsv", ["accession", "x"],
                                 [{"accession": "p2", "x": "2"}, {"accession": "p1", "x": "1"}])
            report = module.compare(shipped, revision)
            self.assertIn("NOT ADDITIVE", report["verdict"])
            self.assertIn("accession", report["mismatch_examples"][0])

    def test_a_row_count_difference_is_not_additive(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["accession", "x"],
                                [{"accession": "p1", "x": "1"}, {"accession": "p2", "x": "2"}])
            revision = write_tsv(tmp / "b.tsv", ["accession", "x"], [{"accession": "p1", "x": "1"}])
            report = module.compare(shipped, revision)
            self.assertIn("NOT ADDITIVE", report["verdict"])
            self.assertEqual(report["shipped_rows"], 2)
            self.assertEqual(report["revision_rows"], 1)

    def test_a_removed_column_is_not_additive(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["accession", "x", "y"],
                                [{"accession": "p1", "x": "1", "y": "2"}])
            revision = write_tsv(tmp / "b.tsv", ["accession", "x"], [{"accession": "p1", "x": "1"}])
            report = module.compare(shipped, revision)
            self.assertEqual(report["removed_columns"], ["y"])
            self.assertIn("NOT ADDITIVE", report["verdict"])

    def test_an_empty_added_column_is_reported_as_empty(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["accession"], [{"accession": "p1"}])
            revision = write_tsv(tmp / "b.tsv", ["accession", "gap"], [{"accession": "p1"}])
            report = module.compare(shipped, revision)
            self.assertEqual(report["added_column_fill"]["gap"],
                             {"rows": 1, "non_empty": 0})

    def test_missing_accession_column_is_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["other"], [{"other": "1"}])
            revision = write_tsv(tmp / "b.tsv", ["accession"], [{"accession": "p1"}])
            with self.assertRaisesRegex(ValueError, "must carry an accession column"):
                module.compare(shipped, revision)


class CliTests(unittest.TestCase):
    def test_cli_writes_only_when_asked_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            shipped = write_tsv(tmp / "a.tsv", ["accession"], [{"accession": "p1"}])
            revision = write_tsv(tmp / "b.tsv", ["accession", "gap"],
                                 [{"accession": "p1", "gap": "1.0"}])
            out = tmp / "report.json"
            self.assertEqual(module.main(["--shipped", str(shipped),
                                          "--revision", str(revision)]), 0)
            self.assertFalse(out.exists(), "no --out means nothing is written")
            self.assertEqual(module.main(["--shipped", str(shipped), "--revision", str(revision),
                                          "--out", str(out)]), 0)
            payload = json.loads(out.read_text("utf-8"))
            self.assertIn("boundary", payload)
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                module.main(["--shipped", str(shipped), "--revision", str(revision),
                             "--out", str(out)])


if __name__ == "__main__":
    unittest.main()
