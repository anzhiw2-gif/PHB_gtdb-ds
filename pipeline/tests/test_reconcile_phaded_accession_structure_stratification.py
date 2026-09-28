import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "reconcile_phaded_accession_structure_stratification.py"


def load_module():
    spec = importlib.util.spec_from_file_location("reconcile_phaded_accession_structure_stratification", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path, header, rows):
    path.write_text("\t".join(header) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8")


class ReconcileTests(unittest.TestCase):
    def test_stratify_and_bridge_are_explicit_and_fail_closed(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            assignment = root / "assignment.tsv"
            evidence = root / "evidence.tsv"
            manual = root / "manual.tsv"
            bridge = root / "bridge.tsv"
            stats = root / "stats.tsv"
            holds = root / "holds.tsv"
            matches = root / "matches.tsv"
            write_tsv(assignment, ["accession", "genome", "layer", "phaded_superfamily_best", "phaded_family_best", "assignment_status"], [
                ["p1", "G1", "L1", "sf1", "fam1", "assigned"],
                ["p2", "G2", "L2", "", "", "unassigned_PhaDED_like"],
            ])
            write_tsv(evidence, ["accession", "architecture_consistency", "feature_evidence_status"], [
                ["p1", "consistent", "supported"], ["p2", "conflicting", "conflict"],
            ])
            write_tsv(manual, ["candidate_id", "manual_review_decision", "superfamily"], [["p2", "hold_architecture_conflict", "sf2"]])
            write_tsv(bridge, ["record_id", "gtdb_genome_accession", "evidence_tier", "match_class", "pmid"], [["r1", "G1", "direct_literature_supported", "exact_genome_accession", "1"], ["r2", "G9", "sequence_only_candidate", "taxonomy_only_match", "2"]])
            result = module.run_reconciliation(assignment, evidence, manual, bridge, stats, holds, matches)
            self.assertEqual(result["candidate_count"], 2)
            self.assertEqual(result["exact_candidate_genome_count"], 1)
            self.assertEqual(result["manual_review_count"], 1)
            self.assertIn("unassigned_PhaDED_like", stats.read_text(encoding="utf-8"))
            self.assertIn("hold_architecture_conflict", holds.read_text(encoding="utf-8"))
            self.assertIn("direct_literature_supported", matches.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
