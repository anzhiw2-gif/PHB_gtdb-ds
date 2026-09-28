import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_phaded_reference_panel.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("prepare_phaded_reference_panel", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReferencePanelPreparationTests(unittest.TestCase):
    def test_creates_three_planned_roles_per_reference_only_profile(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "profile_manifest.tsv"
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "profile_id",
                        "profile_kind",
                        "phaded_superfamily",
                        "phaded_family_id",
                        "model_status",
                        "model_reason",
                    ],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "profile_id": "family_DED_hfam_44_hash",
                            "profile_kind": "family",
                            "phaded_superfamily": "extracellular dPHASCL type 1",
                            "phaded_family_id": "DED_hfam_44",
                            "model_status": "reference_only",
                            "model_reason": "fewer_than_three_independent_experimental_positive_records",
                        },
                        {
                            "profile_id": "superfamily_intracellular_hash",
                            "profile_kind": "superfamily",
                            "phaded_superfamily": "intracellular nPHAMCL",
                            "phaded_family_id": "",
                            "model_status": "reference_only",
                            "model_reason": "fewer_than_three_independent_experimental_positive_records",
                        },
                        {
                            "profile_id": "family_DED_hfam_52_hash",
                            "profile_kind": "family",
                            "phaded_superfamily": "extracellular dPHASCL type 1",
                            "phaded_family_id": "DED_hfam_52",
                            "model_status": "trained",
                            "model_reason": "sufficient_independent_experimental_positive_records",
                        },
                    ]
                )

            result = module.prepare(manifest, root / "out", expected_reference_only_count=2)

            self.assertEqual(result["status"], "planned_not_run")
            self.assertEqual(result["reference_only_profile_count"], 2)
            self.assertEqual(result["record_count"], 6)
            with (root / "out" / "reference_panel_acquisition.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["role"] for row in rows}, {"independent_positive", "family_resolved_negative", "challenge_control"})
            self.assertEqual(set(rows[0]), set(module.OUTPUT_FIELDS))
            self.assertEqual({row["evidence_status"] for row in rows}, {"planned_not_run"})
            self.assertEqual({row["retrieval_status"] for row in rows}, {"pending"})
            self.assertEqual({row["independence_status"] for row in rows}, {"pending"})
            self.assertTrue(all(row["accession"] == "" for row in rows))
            self.assertTrue(all(row["sequence_sha256"] == "pending" for row in rows))
            self.assertTrue(all(row["source"] for row in rows))
            report = (root / "out" / "reference_panel_acquisition.json").read_text(encoding="utf-8")
            self.assertIn("cannot produce new family assignments", report)

    def test_rejects_duplicate_profile_and_non_reference_status(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "profile_manifest.tsv"
            manifest.write_text(
                "profile_id\tprofile_kind\tphaded_superfamily\tphaded_family_id\tmodel_status\tmodel_reason\n"
                "p1\tfamily\tsf\tDED_hfam_1\treference_only\tinsufficient\n"
                "p1\tfamily\tsf\tDED_hfam_2\treference_only\tinsufficient\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate profile_id"):
                module.prepare(manifest, root / "out", expected_reference_only_count=2)

    def test_rejects_reference_only_count_mismatch_and_does_not_assign_family(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "profile_manifest.tsv"
            manifest.write_text(
                "profile_id\tprofile_kind\tphaded_superfamily\tphaded_family_id\tmodel_status\tmodel_reason\n"
                "p1\tfamily\tsf\tDED_hfam_1\treference_only\tinsufficient\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "expected 2 reference-only profiles"):
                module.prepare(manifest, root / "out", expected_reference_only_count=2)


if __name__ == "__main__":
    unittest.main()
