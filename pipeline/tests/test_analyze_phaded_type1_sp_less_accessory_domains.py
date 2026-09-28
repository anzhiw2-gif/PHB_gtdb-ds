"""Tests for the P2 SP-less type-1 accessory-domain analysis.

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

import analyze_phaded_type1_sp_less_accessory_domains as module  # noqa: E402


HOLD_FIELDS = [
    "accession", "genome", "superfamily", "hold_reason", "catalytic_domain_type",
    "signalp_class",
]
MATRIX_FIELDS = ["accession", "pfam_accessions", "interpro_signatures"]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


def hold_row(accession: str, **overrides) -> dict[str, str]:
    row = {
        "accession": accession,
        "genome": "G1",
        "superfamily": module.TYPE1_SUPERFAMILY,
        "hold_reason": module.LOCALIZATION_HOLD_REASON,
        "catalytic_domain_type": module.TYPE1_GEOMETRY,
        "signalp_class": "OTHER",
    }
    row.update(overrides)
    return row


class MarkerSetTests(unittest.TestCase):
    def test_marker_set_matches_the_frozen_layer(self):
        self.assertEqual(
            module.MARKERS["fn3"],
            ("PF00041", "SM00060", "cd00063", "SSF49265", "PS50853", "G3DSA:2.60.40.10"),
        )
        self.assertEqual(module.MARKERS["ig_like"], ("PF16403", "PF17957"))
        self.assertEqual(module.MARKERS["tsp3"], ("PF02412",))
        self.assertEqual(module.MARKERS["chb"], ("PF13290",))

    def test_tokens_splits_every_separator(self):
        self.assertEqual(
            module.tokens("PF00041;PF02412", "PF13290,PF17957", "SM00060|cd00063"),
            {"PF00041", "PF02412", "PF13290", "PF17957", "SM00060", "cd00063"},
        )


class SelectionTests(unittest.TestCase):
    def test_selects_only_the_three_key_intersection(self):
        rows = [
            hold_row("keep"),
            hold_row("other-reason", hold_reason="architecture_criteria"),
            hold_row("other-geometry", catalytic_domain_type="type2_verified"),
            hold_row("other-superfamily", superfamily="extracellular dPHASCL type 2"),
        ]
        selected = module.select_sp_less_type1(rows)
        self.assertEqual([row["accession"] for row in selected], ["keep"])

    def test_duplicate_accession_raises(self):
        with self.assertRaisesRegex(ValueError, "not accession-unique"):
            module.select_sp_less_type1([hold_row("dup"), hold_row("dup")])

    def test_missing_column_raises(self):
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            module.select_sp_less_type1([{"accession": "a"}])


class PresenceTests(unittest.TestCase):
    def test_presence_is_read_from_pfam_and_interpro(self):
        hits = module.accessory_hit("", "PF00041;PF10503", "")
        self.assertEqual(hits["fn3"], "yes")
        self.assertEqual(hits["ig_like"], "no")

    def test_interpro_only_marker_is_detected(self):
        hits = module.accessory_hit("G3DSA:2.60.40.10", "", "")
        self.assertEqual(hits["fn3"], "yes")

    def test_no_signature_means_no_not_pending(self):
        hits = module.accessory_hit("PF10503", "PF10503", "")
        self.assertEqual(hits, {name: "no" for name in module.MARKER_ORDER})


class JoinTests(unittest.TestCase):
    def test_row_absent_from_matrix_is_pending_not_absent(self):
        target = [hold_row("a1")]
        joined, summary = module.classify_rows(target, {})
        self.assertEqual(joined[0]["in_matrix"], "no")
        self.assertEqual(joined[0]["has_fn3"], "pending")
        self.assertEqual(summary["target_rows_absent_from_matrix"], ["a1"])

    def test_counts_and_combinations(self):
        target = [hold_row("a1"), hold_row("a2"), hold_row("a3")]
        matrix = {
            "a1": {"accession": "a1", "pfam_accessions": "PF00041;PF10503",
                   "interpro_signatures": ""},
            "a2": {"accession": "a2", "pfam_accessions": "PF10503",
                   "interpro_signatures": ""},
            "a3": {"accession": "a3", "pfam_accessions": "PF00041;PF02412;PF13290",
                   "interpro_signatures": ""},
        }
        _joined, summary = module.classify_rows(target, matrix)
        self.assertEqual(summary["target_count"], 3)
        self.assertEqual(summary["with_any_accessory_domain"], 2)
        self.assertEqual(summary["per_domain_counts"]["fn3"], 2)
        self.assertEqual(summary["per_domain_counts"]["tsp3"], 1)
        self.assertIn("OPTIONAL-module", summary["interpretation_note"])
        self.assertIn("NEVER a pass/hold criterion", summary["interpretation_note"])


class CliTests(unittest.TestCase):
    def _fixture(self, tmp: Path) -> tuple[Path, Path]:
        hold = write_tsv(tmp / "hold.tsv", HOLD_FIELDS, [hold_row("a1")])
        matrix = write_tsv(tmp / "matrix.tsv", MATRIX_FIELDS, [
            {"accession": "a1", "pfam_accessions": "PF00041", "interpro_signatures": ""},
        ])
        return hold, matrix

    def test_cli_writes_artifacts_and_reports_a_target_mismatch(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            hold, matrix = self._fixture(tmp)
            out = tmp / "out"
            self.assertEqual(
                module.main(["--hold", str(hold), "--matrix", str(matrix),
                             "--out-dir", str(out), "--expect-targets", "11610"]), 0,
            )
            summary = json.loads((out / "type1_sp_less_accessory_summary.json").read_text("utf-8"))
            self.assertEqual(summary["target_count"], 1)
            self.assertEqual(summary["target_count_mismatch"]["declared"], 11610)
            self.assertEqual(summary["target_count_mismatch"]["observed"], 1)

    def test_cli_refuses_a_non_empty_output_dir(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            hold, matrix = self._fixture(tmp)
            out = tmp / "out"
            out.mkdir()
            (out / "existing").write_text("x", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--hold", str(hold), "--matrix", str(matrix), "--out-dir", str(out)])

    def test_missing_input_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            _hold, matrix = self._fixture(tmp)
            with self.assertRaises(ValueError):
                module.main(["--hold", "/nonexistent/hold.tsv", "--matrix", str(matrix),
                             "--out-dir", str(tmp / "out")])


if __name__ == "__main__":
    unittest.main()
