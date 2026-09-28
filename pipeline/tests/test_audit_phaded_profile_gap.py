import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_phaded_profile_gap.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("audit_phaded_profile_gap", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ProfileGapAuditTests(unittest.TestCase):
    def test_reference_only_profile_has_zero_direct_hits_and_parent_pool(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            manifest = root / "profile_manifest.tsv"
            ledger = root / "ledger.tsv"
            assignments = root / "assignments.tsv"
            scores = root / "scores.tsv"
            panel = root / "panel.tsv"
            output = root / "gap.tsv"
            _write_tsv(
                manifest,
                ["profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id", "training_count", "model_status", "model_reason"],
                [{"profile_id": "family_DED_hfam_2_hash", "profile_kind": "family", "phaded_superfamily": "sf", "phaded_family_id": "DED_hfam_2", "training_count": "1", "model_status": "reference_only", "model_reason": "fewer_than_three_independent_experimental_positive_records"}],
            )
            _write_tsv(
                ledger,
                ["accession", "phaded_superfamily", "phaded_family_id", "evidence_status", "organism", "positive_negative_control"],
                [{"accession": "P1", "phaded_superfamily": "sf", "phaded_family_id": "DED_hfam_2", "evidence_status": "experimental_positive", "organism": "Org1", "positive_negative_control": "not_assessed"}],
            )
            _write_tsv(
                assignments,
                ["accession", "phaded_superfamily_best", "phaded_family_best", "assignment_status"],
                [{"accession": "C1", "phaded_superfamily_best": "sf", "phaded_family_best": "", "assignment_status": "unassigned_PhaDED_like"}],
            )
            _write_tsv(scores, ["accession", "profile_id", "score"], [])
            _write_tsv(panel, ["accession", "panel"], [])
            report = module.audit(
                manifest, ledger, assignments, scores, panel, output,
            )
            row = report["profiles"][0]
            self.assertEqual(row["direct_score_row_count"], 0)
            self.assertEqual(row["direct_best_candidate_count"], 0)
            self.assertEqual(row["parent_superfamily_candidate_count"], 1)
            self.assertEqual(row["reference_experimental_positive_count"], 1)
            self.assertEqual(row["independent_positive_sufficiency"], "insufficient")

    def test_panel_evidence_is_not_family_resolved_without_ledger_overlap(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            manifest = root / "profile_manifest.tsv"
            ledger = root / "ledger.tsv"
            assignments = root / "assignments.tsv"
            scores = root / "scores.tsv"
            panel = root / "panel.tsv"
            output = root / "gap.tsv"
            _write_tsv(manifest, ["profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id", "training_count", "model_status", "model_reason"], [{"profile_id": "superfamily_x", "profile_kind": "superfamily", "phaded_superfamily": "sf", "phaded_family_id": "", "training_count": "1", "model_status": "reference_only", "model_reason": "insufficient"}])
            _write_tsv(ledger, ["accession", "phaded_superfamily", "phaded_family_id", "evidence_status", "organism", "positive_negative_control"], [])
            _write_tsv(assignments, ["accession", "phaded_superfamily_best", "phaded_family_best", "assignment_status"], [])
            _write_tsv(scores, ["accession", "profile_id", "score"], [])
            _write_tsv(panel, ["accession", "panel"], [{"accession": "P1", "panel": "fragment_challenge"}])
            report = module.audit(manifest, ledger, assignments, scores, panel, output)
            row = report["profiles"][0]
            self.assertEqual(row["external_panel_exact_profile_count"], 0)
            self.assertEqual(row["external_panel_unresolved_count"], 1)


if __name__ == "__main__":
    unittest.main()
