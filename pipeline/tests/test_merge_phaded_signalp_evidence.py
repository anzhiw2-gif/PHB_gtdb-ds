import csv
import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.merge_phaded_signalp_evidence import (
    merge_signalp_into_features,
    parse_signalp_results,
)


class SignalPMergeTests(unittest.TestCase):
    def test_parses_prediction_types_and_marks_missing_accessions(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "prediction_results.txt"
            path.write_text(
                "# header\n"
                "# ID\tPrediction\tOTHER\tSP(Sec/SPI)\n"
                "A\tSP\t0.1\t0.9\n"
                "B\tOTHER\t0.9\t0.1\n",
                encoding="utf-8",
            )
            rows = parse_signalp_results(path, ["A", "B", "C"], set())
        self.assertEqual(rows["A"]["signalp_class"], "SP")
        self.assertEqual(rows["B"]["signalp_class"], "OTHER")
        self.assertEqual(rows["C"]["signalp_status"], "signalp_not_reported")

    def test_merges_signalp_without_dropping_feature_rows(self):
        features = [
            {"accession": "A", "signalp_evidence": "pending", "signalp_class": "pending"},
            {"accession": "B", "signalp_evidence": "pending", "signalp_class": "pending"},
        ]
        signalp = {
            "A": {"signalp_class": "SP", "signalp_status": "signalp_supported"},
            "B": {"signalp_class": "OTHER", "signalp_status": "signalp_supported"},
        }
        merged = merge_signalp_into_features(features, signalp)
        self.assertEqual(len(merged), 2)
        self.assertEqual(merged[0]["signalp_evidence"], "signalp_supported")
        self.assertEqual(merged[1]["signalp_class"], "OTHER")
