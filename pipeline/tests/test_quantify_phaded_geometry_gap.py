"""Tests for the F14 geometry-gap quantification.

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

import quantify_phaded_geometry_gap as module  # noqa: E402


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


class VerdictTests(unittest.TestCase):
    def test_verified_types_pass_through(self):
        self.assertEqual(module.normalize_verdict("type1_verified"), "type1_verified")
        self.assertEqual(module.normalize_verdict("type2_verified"), "type2_verified")

    def test_undetermined_verdicts_stay_unresolved(self):
        for value in ("undetermined_no_oxyanion", "undetermined_position_conflict", "", "other"):
            self.assertEqual(module.normalize_verdict(value), module.UNRESOLVED, value)


class QuantifyTests(unittest.TestCase):
    def test_only_unresolved_rows_are_examined(self):
        adapted = [
            {"accession": "a", "catalytic_domain_type": "type1_verified"},
            {"accession": "b", "catalytic_domain_type": module.UNRESOLVED},
        ]
        geometry = [{"accession": "a", "catalytic_domain_type": "type1_verified"},
                    {"accession": "b", "catalytic_domain_type": "type2_verified"}]
        rows, summary = module.quantify(adapted, geometry)
        self.assertEqual([row["accession"] for row in rows], ["b"])
        self.assertEqual(summary["unresolved_before"], 1)
        self.assertEqual(summary["closed_by_the_geometry_table"], 1)
        self.assertEqual(summary["closed_breakdown"], {"type2_verified": 1})

    def test_an_undetermined_frozen_verdict_does_not_close_the_gap(self):
        adapted = [{"accession": "b", "catalytic_domain_type": module.UNRESOLVED}]
        geometry = [{"accession": "b", "catalytic_domain_type": "undetermined_no_oxyanion"}]
        rows, summary = module.quantify(adapted, geometry)
        self.assertEqual(rows[0]["would_become"], module.UNRESOLVED)
        self.assertEqual(rows[0]["closes_the_gap"], "no")
        self.assertEqual(summary["closed_by_the_geometry_table"], 0)
        self.assertEqual(summary["still_unresolved"], 1)
        self.assertIn("undetermined_no_oxyanion", summary["still_unresolved_breakdown"])

    def test_an_accession_absent_from_the_geometry_table_is_reported(self):
        adapted = [{"accession": "gone", "catalytic_domain_type": module.UNRESOLVED}]
        rows, summary = module.quantify(adapted, [{"accession": "other",
                                                   "catalytic_domain_type": "type1_verified"}])
        self.assertEqual(rows[0]["geometry_verdict"], "absent_from_the_geometry_table")
        self.assertEqual(summary["still_unresolved_breakdown"],
                         {"absent_from_the_geometry_table": 1})

    def test_duplicate_geometry_accessions_are_refused(self):
        with self.assertRaisesRegex(ValueError, "duplicate accessions"):
            module.quantify(
                [{"accession": "b", "catalytic_domain_type": module.UNRESOLVED}],
                [{"accession": "b", "catalytic_domain_type": "type1_verified"},
                 {"accession": "b", "catalytic_domain_type": "type2_verified"}],
            )

    def test_row_without_accession_is_refused(self):
        with self.assertRaisesRegex(ValueError, "adapted row without an accession"):
            module.quantify([{"catalytic_domain_type": module.UNRESOLVED}], [])

    def test_closure_rate_and_boundary_are_reported(self):
        adapted = [
            {"accession": "a", "catalytic_domain_type": module.UNRESOLVED},
            {"accession": "b", "catalytic_domain_type": module.UNRESOLVED},
            {"accession": "c", "catalytic_domain_type": module.UNRESOLVED},
            {"accession": "d", "catalytic_domain_type": module.UNRESOLVED},
        ]
        geometry = [
            {"accession": "a", "catalytic_domain_type": "type1_verified"},
            {"accession": "b", "catalytic_domain_type": "type2_verified"},
            {"accession": "c", "catalytic_domain_type": "undetermined_no_oxyanion"},
            {"accession": "d", "catalytic_domain_type": "undetermined_position_conflict"},
        ]
        _rows, summary = module.quantify(adapted, geometry)
        self.assertAlmostEqual(summary["closure_rate"], 0.5)
        self.assertIn("third source", summary["interpretation"])
        self.assertIn("no disposition changes", summary["boundary"])


class CliTests(unittest.TestCase):
    def test_cli_writes_both_artifacts_and_refuses_a_non_empty_dir(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            adapted = write_tsv(tmp / "adapted.tsv", ["accession", "catalytic_domain_type"], [
                {"accession": "a", "catalytic_domain_type": module.UNRESOLVED},
            ])
            geometry = write_tsv(tmp / "geometry.tsv", ["accession", "catalytic_domain_type"], [
                {"accession": "a", "catalytic_domain_type": "type1_verified"},
            ])
            out = tmp / "out"
            self.assertEqual(
                module.main(["--adapted-input", str(adapted), "--geometry", str(geometry),
                             "--out-dir", str(out)]), 0,
            )
            summary = json.loads((out / "geometry_gap_summary.json").read_text("utf-8"))
            self.assertEqual(summary["closed_by_the_geometry_table"], 1)
            self.assertEqual(summary["geometry_rows"], 1)
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--adapted-input", str(adapted), "--geometry", str(geometry),
                             "--out-dir", str(out)])

    def test_missing_input_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            geometry = write_tsv(tmp / "geometry.tsv", ["accession", "catalytic_domain_type"], [
                {"accession": "a", "catalytic_domain_type": "type1_verified"},
            ])
            with self.assertRaises(ValueError):
                module.main(["--adapted-input", "/nonexistent/adapted.tsv",
                             "--geometry", str(geometry), "--out-dir", str(tmp / "out")])


if __name__ == "__main__":
    unittest.main()
