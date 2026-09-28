import tempfile
import unittest
from pathlib import Path

from pipeline.scripts.merge_phaded_foldseek_evidence import (
    merge_foldseek_into_matrix,
    read_foldseek_summary,
    merge_files,
)


class FoldseekMergeTests(unittest.TestCase):
    def test_reads_accession_bound_foldseek_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "structural_evidence_summary.tsv"
            path.write_text(
                "candidate_id\tclass\treference_model_id\tfident\talnlen\talntmscore\tqtmscore\tttmscore\trmsd\tstructural_evidence_band\tinterpretation\n"
                "A\tdPHASCL1\tAF-X-F1\t0.4\t300\t0.69\t0.64\t0.68\t10.0\tmoderate\tcandidate-only support\n",
                encoding="utf-8",
            )
            rows = read_foldseek_summary(path)
        self.assertEqual(rows["A"]["foldseek_status"], "foldseek_support_candidate_only")
        self.assertEqual(rows["A"]["foldseek_reference_model_id"], "AF-X-F1")

    def test_merge_preserves_all_matrix_rows_and_marks_unmatched_pending(self):
        matrix = [{"accession": "A"}, {"accession": "B"}]
        foldseek = {
            "A": {
                "foldseek_status": "foldseek_support_candidate_only",
                "foldseek_reference_model_id": "AF-X-F1",
                "foldseek_alntmscore": "0.69",
                "foldseek_interpretation": "candidate-only support",
            }
        }
        merged = merge_foldseek_into_matrix(matrix, foldseek)
        self.assertEqual(len(merged), 2)
        self.assertEqual(merged[0]["foldseek_status"], "foldseek_support_candidate_only")
        self.assertEqual(merged[1]["foldseek_status"], "foldseek_not_tested_full_library")

    def test_merge_files_writes_run_contract_with_bound_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "run"
            inputs = root / "inputs"
            results = root / "results"
            inputs.mkdir(parents=True)
            results.mkdir(parents=True)
            matrix = inputs / "matrix.tsv"
            summary = inputs / "summary.tsv"
            output = results / "merged.tsv"
            matrix.write_text("accession\nA\n", encoding="utf-8")
            summary.write_text(
                "candidate_id\tclass\treference_model_id\tstructural_evidence_band\n"
                "A\tdPHASCL1\tAF-X-F1\tmoderate\n",
                encoding="utf-8",
            )
            merge_files(matrix, summary, output)
            contract = root / "input_contract.json"
            self.assertTrue(contract.exists())
            self.assertIn("foldseek_summary", contract.read_text(encoding="utf-8"))
