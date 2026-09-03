import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "stratify_ephaz_ambiguous.py"


class StratifiedReviewTests(unittest.TestCase):
    def _load(self):
        self.assertTrue(SCRIPT.is_file(), SCRIPT)
        spec = importlib.util.spec_from_file_location("stratify_ephaz_ambiguous", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module

    def test_bins_and_pending_evidence_are_explicit(self):
        module = self._load()
        row = module.enrich_row(
            {"accession": "GCA_1|contig_1", "length": "230", "ephaz_bitscore": "40", "iphaz_bitscore": "43", "assignment": "ambiguous"},
            {"GCA_1|contig_1": {"ephaz": 0.42, "iphaz": 0.77}},
            {"GCA_1|contig_1": "marker_present"},
        )
        self.assertEqual(row["delta_abs_bin"], "2-<5")
        self.assertEqual(row["length_bin"], "<250")
        self.assertEqual(row["domain_coverage_bin"], "0.5-<0.8")
        self.assertEqual(row["neighborhood_bin"], "marker_present")
        pending = module.enrich_row(
            {"accession": "GCA_2|contig_2", "length": "700", "ephaz_bitscore": "50", "iphaz_bitscore": "51", "assignment": "ambiguous"},
            {}, {},
        )
        self.assertEqual(pending["domain_coverage_bin"], "pending")
        self.assertEqual(pending["neighborhood_bin"], "no_record")

    def test_sampling_is_deterministic_and_capped_per_stratum(self):
        module = self._load()
        rows = [{"accession": f"A{i}", "stratum": "same"} for i in range(20)]
        first = module.sample_rows(rows, per_stratum=3, seed="20260829")
        second = module.sample_rows(list(reversed(rows)), per_stratum=3, seed="20260829")
        self.assertEqual([row["accession"] for row in first], [row["accession"] for row in second])
        self.assertEqual(len(first), 3)

    def test_domtblout_uses_hmm_coordinates_and_full_candidate_id(self):
        module = self._load()
        # domtblout columns 15-16 are HMM coordinates; 19-20 are envelope
        # coordinates and must not be used for model coverage.
        fields = ["GCA_1|contig_1", "-", "500", "query", "-", "400", "1e-20", "80", "0", "1", "1", "1", "1", "1", "1", "101", "300", "1", "200", "1", "200", "0.9", "desc"]
        with tempfile.TemporaryDirectory() as tmp:
            dom = Path(tmp) / "hits.dom"
            dom.write_text(" ".join(fields) + "\n", encoding="utf-8")
            coverage = module.parse_domtblout([dom])
        self.assertAlmostEqual(coverage["GCA_1|contig_1"], 0.5)

    def test_domtblout_merges_non_overlapping_hmm_intervals(self):
        module = self._load()
        base = ["GCA_1|contig_1", "-", "500", "query", "-", "400", "1e-20", "80", "0", "1", "1", "1", "1", "1", "1"]
        first = base + ["1", "100", "1", "100", "1", "100", "0.9", "desc"]
        second = base + ["201", "300", "101", "200", "101", "200", "0.9", "desc"]
        with tempfile.TemporaryDirectory() as tmp:
            dom = Path(tmp) / "hits.dom"
            dom.write_text(" ".join(first) + "\n" + " ".join(second) + "\n", encoding="utf-8")
            coverage = module.parse_domtblout([dom])
        self.assertAlmostEqual(coverage["GCA_1|contig_1"], 0.5)

    def test_missing_required_competition_column_fails_closed(self):
        module = self._load()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.tsv"
            path.write_text("accession\tlength\nA\t100\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                module.read_competition(path)

    def test_sample_fasta_keeps_only_requested_full_ids(self):
        module = self._load()
        records = {"GCA_1|contig_1": ("GCA_1|contig_1 description", "MMMM"), "GCA_1|contig_2": ("GCA_1|contig_2", "AAAA")}
        output = module.fasta_text(records, ["GCA_1|contig_2"])
        self.assertEqual(output, ">GCA_1|contig_2\nAAAA\n")


if __name__ == "__main__":
    unittest.main()
