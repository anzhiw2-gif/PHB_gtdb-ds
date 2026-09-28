"""Tests for the reference-panel SignalP gate calibration (packet P1).

Synthetic fixtures only: the real 723-reference panel and the real SignalP run
live under runs/ and are read by the run, not by these tests.
"""
import csv
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import analyze_phaded_reference_signalp as module  # noqa: E402


SIGNALP_HEADER = (
    "# SignalP-6.0\tOrganism: Other\tTimestamp: 20260928221440\n"
    "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tLIPO(Sec/SPII)\tTAT(Tat/SPI)\t"
    "TATLIPO(Tat/SPII)\tPILIN(Sec/SPIII)\tCS Position\n"
)


def write_signalp(tmp: Path, rows: list[tuple[str, str]]) -> Path:
    path = tmp / "prediction_results.txt"
    body = "".join(
        f"{sid}\t{pred}\t1.000000\t0.000000\t0.000000\t0.000000\t0.000000\t0.000000\t\n"
        for sid, pred in rows
    )
    path.write_text(SIGNALP_HEADER + body, encoding="utf-8")
    return path


def write_ledger(tmp: Path, rows: list[dict[str, str]]) -> Path:
    path = tmp / "ledger.tsv"
    fields = [
        "reference_id", "accession", "phaded_superfamily", "phaded_family_id",
        "experimental_evidence_grade", "sequence_integrity",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


def ledger_row(reference_id: str, superfamily: str, **overrides) -> dict[str, str]:
    row = {
        "reference_id": reference_id,
        "accession": overrides.pop("accession", reference_id.split("_")[-1]),
        "phaded_superfamily": superfamily,
        "phaded_family_id": "DED_hfam_x",
        "experimental_evidence_grade": overrides.pop("grade", "A"),
        "sequence_integrity": overrides.pop("integrity", "complete"),
    }
    row.update(overrides)
    return row


class ParsingTests(unittest.TestCase):
    def test_parses_ids_and_strips_the_accession(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            path = write_signalp(tmp, [("DED_hfam_1_0001|AAA111.1", "SP")])
            records = module.parse_prediction_results(path)
            self.assertEqual(list(records), ["DED_hfam_1_0001"])
            self.assertEqual(records["DED_hfam_1_0001"]["Prediction"], "SP")
            self.assertEqual(records["DED_hfam_1_0001"]["signalp_id"], "DED_hfam_1_0001|AAA111.1")

    def test_missing_file_raises(self):
        with self.assertRaises(ValueError):
            module.parse_prediction_results(Path("/nonexistent/prediction_results.txt"))

    def test_duplicate_prediction_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = write_signalp(Path(temporary), [("DED_1|A", "SP"), ("DED_1|A", "OTHER")])
            with self.assertRaisesRegex(ValueError, "duplicate SignalP prediction"):
                module.parse_prediction_results(path)

    def test_empty_file_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "prediction_results.txt"
            path.write_text("# only comments\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no SignalP predictions"):
                module.parse_prediction_results(path)


class LocalizationMappingTests(unittest.TestCase):
    def test_inherited_labels_map_from_the_historical_names(self):
        self.assertEqual(
            module.historical_localization("extracellular dPHASCL type 1"), "extracellular",
        )
        self.assertEqual(
            module.historical_localization("intracellular nPHAMCL"), "intracellular",
        )
        self.assertEqual(
            module.historical_localization("periplasmic PHA depolymerases"), "periplasmic",
        )

    def test_unknown_superfamily_is_unmapped_not_guessed(self):
        self.assertEqual(module.historical_localization("some new family"), "unmapped")


class CalibrationTests(unittest.TestCase):
    def test_export_prediction_classes_are_the_four_export_types(self):
        self.assertEqual(
            module.EXPORT_PREDICTIONS, {"SP", "LIPO", "TAT", "TATLIPO"},
        )

    def test_rates_are_computed_per_inherited_stratum(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            ledger = write_ledger(tmp, [
                ledger_row("r1", "extracellular dPHASCL type 1"),
                ledger_row("r2", "extracellular dPHASCL type 1"),
                ledger_row("r3", "intracellular nPHAMCL"),
                ledger_row("r4", "intracellular nPHAMCL"),
            ])
            signalp = write_signalp(tmp, [
                ("r1|A", "SP"), ("r2|A", "OTHER"), ("r3|A", "OTHER"), ("r4|A", "LIPO"),
            ])
            joined, summary = module.classify_rows(
                module.read_reference_rows(ledger), module.parse_prediction_results(signalp),
            )
            self.assertEqual(len(joined), 4)
            self.assertAlmostEqual(summary["extracellular_export_signal_rate"], 0.5)
            self.assertAlmostEqual(summary["intracellular_export_signal_rate"], 0.5)
            self.assertEqual(summary["prediction_class_counts"]["OTHER"], 2)

    def test_absent_prediction_is_not_counted_as_intracellular_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            ledger = write_ledger(tmp, [ledger_row("r1", "extracellular dPHASCL type 1")])
            signalp = write_signalp(tmp, [("other|A", "SP")])
            joined, summary = module.classify_rows(
                module.read_reference_rows(ledger), module.parse_prediction_results(signalp),
            )
            self.assertEqual(joined[0]["has_export_signal"], "")
            self.assertEqual(summary["rows_without_prediction"], ["r1"])
            bucket = summary["by_inherited_localization"]["extracellular"]
            self.assertEqual(bucket["export_signal"], 0)
            self.assertEqual(bucket["no_export_signal"], 0)

    def test_unmapped_superfamily_is_reported_not_bucketed(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            ledger = write_ledger(tmp, [ledger_row("r1", "mystery family")])
            signalp = write_signalp(tmp, [("r1|A", "SP")])
            _joined, summary = module.classify_rows(
                module.read_reference_rows(ledger), module.parse_prediction_results(signalp),
            )
            self.assertEqual(summary["rows_with_unmapped_superfamily"], ["r1"])
            self.assertIn("unmapped", summary["by_inherited_localization"])

    def test_documented_localization_is_carried_through(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            ledger = write_ledger(tmp, [ledger_row("r1", "extracellular dPHASCL type 1")])
            signalp = write_signalp(tmp, [("r1|A", "OTHER")])
            joined, _summary = module.classify_rows(
                module.read_reference_rows(ledger),
                module.parse_prediction_results(signalp),
                experimental_localization={"r1": "extracellular"},
            )
            self.assertEqual(joined[0]["documented_localization"], "extracellular")
            self.assertEqual(joined[0]["inherited_localization"], "extracellular")
            self.assertEqual(joined[0]["has_export_signal"], "false")

    def test_summary_states_the_specificity_boundary(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            ledger = write_ledger(tmp, [ledger_row("r1", "intracellular nPHAMCL")])
            signalp = write_signalp(tmp, [("r1|A", "OTHER")])
            _joined, summary = module.classify_rows(
                module.read_reference_rows(ledger), module.parse_prediction_results(signalp),
            )
            self.assertIn("not experimental", summary["specificity_note"])
            self.assertIn("NOT evidence of an intracellular", summary["specificity_note"])


class CliTests(unittest.TestCase):
    def test_writes_both_artifacts_and_refuses_a_non_empty_dir(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            ledger = write_ledger(tmp, [ledger_row("r1", "intracellular nPHAMCL")])
            signalp = write_signalp(tmp, [("r1|A", "OTHER")])
            out = tmp / "out"
            self.assertEqual(
                module.main(["--ledger", str(ledger), "--signalp", str(signalp),
                             "--out-dir", str(out)]), 0,
            )
            self.assertTrue((out / "reference_signalp_stratification.tsv").is_file())
            payload = json.loads((out / "reference_signalp_calibration.json").read_text("utf-8"))
            self.assertEqual(payload["total_references"], 1)
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--ledger", str(ledger), "--signalp", str(signalp),
                             "--out-dir", str(out)])

    def test_missing_ledger_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            signalp = write_signalp(tmp, [("r1|A", "SP")])
            with self.assertRaises(ValueError):
                module.main(["--ledger", "/nonexistent/ledger.tsv", "--signalp", str(signalp),
                             "--out-dir", str(tmp / "out")])


if __name__ == "__main__":
    unittest.main()
