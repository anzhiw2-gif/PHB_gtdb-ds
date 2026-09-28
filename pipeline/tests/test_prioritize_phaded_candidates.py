import importlib.util
import csv
import unittest
import tempfile
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prioritize_phaded_candidates.py"


def load_module():
    spec = importlib.util.spec_from_file_location("prioritize_phaded_candidates", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PriorityCandidateTests(unittest.TestCase):
    def test_rank_exposes_components_and_keeps_ambiguous_review(self):
        module = load_module()
        rows = [
            {
                "accession": "A", "genome": "G1", "phaded_superfamily_best": "sf",
                "phaded_family_best": "fam", "assignment_status": "assigned",
                "architecture_consistency": "consistent", "sequence_integrity": "valid",
                "signalp_evidence": "accepted", "superfamily_score_gap": "10",
                "reference_novelty": "high",
            },
            {
                "accession": "B", "genome": "G2", "phaded_superfamily_best": "",
                "phaded_family_best": "", "assignment_status": "ambiguous_superfamily",
                "architecture_consistency": "pending", "sequence_integrity": "valid",
                "signalp_evidence": "no_record", "superfamily_score_gap": "0.2",
                "reference_novelty": "low",
            },
        ]
        ranked = module.rank_candidates(rows)
        self.assertEqual(ranked[0]["accession"], "A")
        self.assertIn("priority_score", ranked[0])
        self.assertEqual(ranked[1]["review_bucket"], "ambiguous")

    def test_terminal_stop_has_valid_score_while_invalid_states_keep_penalty(self):
        module = load_module()
        base = {
            "genome": "G1", "phaded_superfamily_best": "sf", "phaded_family_best": "fam",
            "assignment_status": "assigned", "architecture_consistency": "consistent",
            "signalp_evidence": "no_record", "superfamily_score_gap": "0",
        }
        rows = [
            dict(base, accession="valid", sequence_integrity="valid"),
            dict(base, accession="terminal", sequence_integrity="terminal_stop_only"),
            dict(base, accession="internal", sequence_integrity="invalid_internal_character"),
            dict(base, accession="legacy", sequence_integrity="invalid_character"),
        ]
        scores = {row["accession"]: row["priority_score"] for row in module.rank_candidates(rows)}
        self.assertEqual(scores["terminal"], scores["valid"])
        self.assertEqual(scores["internal"], scores["valid"] - 5.0)
        self.assertEqual(scores["legacy"], scores["valid"] - 5.0)

    def test_summary_keeps_protein_and_genome_denominators_separate(self):
        module = load_module()
        rows = [
            {"accession": "A", "genome": "G1", "phaded_superfamily_best": "sf", "phaded_family_best": "fam"},
            {"accession": "B", "genome": "G1", "phaded_superfamily_best": "sf", "phaded_family_best": "fam2"},
        ]
        summary = module.summarize_candidates(rows)
        self.assertEqual(summary["protein_count"], 2)
        self.assertEqual(summary["genome_superfamily_count"], 1)
        self.assertEqual(summary["genome_family_count"], 2)

    def test_representative_selection_deduplicates_exact_sequences(self):
        module = load_module()
        rows = [
            {"accession": "A", "sequence_sha256": "same", "phaded_superfamily_best": "sf", "assignment_status": "assigned"},
            {"accession": "B", "sequence_sha256": "same", "phaded_superfamily_best": "sf", "assignment_status": "ambiguous_family"},
            {"accession": "C", "sequence_sha256": "other", "phaded_superfamily_best": "sf2", "assignment_status": "assigned"},
        ]
        selected = module.select_representatives(rows, max_per_class=2)
        self.assertEqual(len(selected), 2)
        self.assertEqual({row["accession"] for row in selected}, {"A", "C"})

    def test_downstream_outputs_have_fixed_headers_and_boundary(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            ranked = [{"accession": "A", "priority_score": 1.0}]
            summary = {"protein_count": 1, "genome_superfamily_count": 1, "genome_family_count": 1}
            representatives = [{"accession": "A", "selection_reason": "validated_unique_sequence"}]
            module.write_downstream_outputs(ranked, summary, representatives, output_dir)
            with (output_dir / "taxonomy_summary.tsv").open(encoding="utf-8", newline="") as handle:
                self.assertEqual(next(csv.reader(handle, delimiter="\t")), module.SUMMARY_FIELDS)
            with (output_dir / "representative_selection.tsv").open(encoding="utf-8", newline="") as handle:
                self.assertEqual(next(csv.reader(handle, delimiter="\t")), ["accession", "selection_reason"])
            self.assertIn("candidate-only", (output_dir / "interpretation_boundary.md").read_text(encoding="utf-8").lower())

    def test_representatives_without_sequence_hash_are_not_marked_validated(self):
        module = load_module()
        rows = [{"accession": "A", "phaded_superfamily_best": "sf", "assignment_status": "assigned", "sequence_sha256": ""}]
        selected = module.select_representatives(rows)
        self.assertEqual(selected, [])


if __name__ == "__main__":
    unittest.main()
