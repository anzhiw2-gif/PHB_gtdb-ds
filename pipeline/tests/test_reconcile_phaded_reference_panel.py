import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "reconcile_phaded_reference_panel.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("reconcile_phaded_reference_panel", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ReconcileReferencePanelTests(unittest.TestCase):
    def test_versionless_accession_binding_does_not_assign_unresolved_reference_profile(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles = root / "profiles.tsv"
            ledger = root / "ledger.tsv"
            panel = root / "panel.tsv"
            _write(
                profiles,
                ["profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id", "model_status"],
                [{"profile_id": "p1", "profile_kind": "family", "phaded_superfamily": "sf", "phaded_family_id": "DED_hfam_1", "model_status": "reference_only"}],
            )
            _write(
                ledger,
                ["accession", "phaded_superfamily", "phaded_family_id", "evidence_status"],
                [{"accession": "P26495.1", "phaded_superfamily": "other", "phaded_family_id": "DED_hfam_4", "evidence_status": "experimental_positive"}],
            )
            _write(
                panel,
                ["accession", "panel", "pmid", "doi", "sequence_length", "normalized_sequence_sha256", "source_sha256", "independence_status", "decision"],
                [{"accession": "P26495", "panel": "annotation_only_near_neighbor_negative", "pmid": "PMID:1", "doi": "", "sequence_length": "100", "normalized_sequence_sha256": "a" * 64, "source_sha256": "b" * 64, "independence_status": "independent", "decision": "accept_panel"}],
            )
            result = module.reconcile(profiles, ledger, panel, root / "out", expected_reference_only_count=1)
            row = result["records"][0]
            self.assertEqual(row["family_binding_status"], "resolved_existing_family_versionless")
            self.assertEqual(row["bound_family"], "DED_hfam_4")
            self.assertEqual(row["reference_only_profile_id"], "")
            self.assertEqual(row["calibration_eligibility"], "challenge_unresolved_family")

    def test_unresolved_panel_record_stays_unresolved_and_summary_is_planned(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _write(root / "profiles.tsv", ["profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id", "model_status"], [{"profile_id": "p1", "profile_kind": "family", "phaded_superfamily": "sf", "phaded_family_id": "DED_hfam_1", "model_status": "reference_only"}])
            _write(root / "ledger.tsv", ["accession", "phaded_superfamily", "phaded_family_id", "evidence_status"], [])
            _write(root / "panel.tsv", ["accession", "panel", "pmid", "doi", "sequence_length", "normalized_sequence_sha256", "source_sha256", "independence_status", "decision"], [{"accession": "X1", "panel": "independent_experimental_positive", "pmid": "PMID:2", "doi": "", "sequence_length": "100", "normalized_sequence_sha256": "c" * 64, "source_sha256": "d" * 64, "independence_status": "independent", "decision": "accept_panel"}])
            result = module.reconcile(root / "profiles.tsv", root / "ledger.tsv", root / "panel.tsv", root / "out", expected_reference_only_count=1)
            row = result["records"][0]
            self.assertEqual(row["family_binding_status"], "unresolved")
            self.assertEqual(row["calibration_eligibility"], "not_eligible_unresolved_family")
            self.assertEqual(result["summary"]["calibrated_candidate_model"], 0)
            self.assertEqual(result["summary"]["reference_only_insufficient_panel"], 1)


if __name__ == "__main__":
    unittest.main()
