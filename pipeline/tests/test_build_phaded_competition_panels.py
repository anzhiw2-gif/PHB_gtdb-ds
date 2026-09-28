import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_phaded_competition_panels.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_phaded_competition_panels", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CompetitionPanelTests(unittest.TestCase):
    def test_builds_four_disjoint_candidate_panels_with_reference_and_challenge_roles(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            priority = root / "priority.tsv"
            dossier = root / "dossier.tsv"
            interpro = root / "interpro.tsv"
            candidates = root / "candidates.faa"
            references = root / "references.faa"
            ledger = root / "ledger.tsv"
            rows = [
                {"accession": "A", "phaded_family_best": "family_DED_hfam_70_hash", "phaded_superfamily_best": "extracellular dPHASCL type 2"},
                {"accession": "B", "phaded_family_best": "family_DED_hfam_4_hash", "phaded_superfamily_best": "intracellular nPHAMCL"},
                {"accession": "C", "phaded_family_best": "family_DED_hfam_55_hash", "phaded_superfamily_best": "extracellular dPHASCL type 1"},
                {"accession": "D", "phaded_family_best": "family_DED_hfam_70_hash", "phaded_superfamily_best": "extracellular dPHASCL type 2"},
            ]
            with priority.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader(); writer.writerows(rows)
            dossier_fields = ["accession", "provisional_label", "review_order", "phaded_family_best", "phaded_superfamily_best"]
            with dossier.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=dossier_fields, delimiter="\t")
                writer.writeheader()
                writer.writerows([
                    {"accession": "A", "provisional_label": "A_retain_high_priority_architecture", "review_order": "1", "phaded_family_best": rows[0]["phaded_family_best"], "phaded_superfamily_best": rows[0]["phaded_superfamily_best"]},
                    {"accession": "B", "provisional_label": "E_hold_synthase_like_structure", "review_order": "2", "phaded_family_best": rows[1]["phaded_family_best"], "phaded_superfamily_best": rows[1]["phaded_superfamily_best"]},
                    {"accession": "C", "provisional_label": "E_hold_structure_assignment_review", "review_order": "3", "phaded_family_best": rows[2]["phaded_family_best"], "phaded_superfamily_best": rows[2]["phaded_superfamily_best"]},
                ])
            with interpro.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["accession", "review_stage", "review_order", "phaded_family_best", "phaded_superfamily_best"], delimiter="\t")
                writer.writeheader(); writer.writerow({"accession": "D", "review_stage": "ePhaZ_competition_review", "review_order": "4", "phaded_family_best": rows[3]["phaded_family_best"], "phaded_superfamily_best": rows[3]["phaded_superfamily_best"]})
            candidates.write_text(">A\nMAA*\n>B\nMYYY\n>C\nMCCC\n>D\nMDDD\n", encoding="ascii")
            references.write_text(">DED_hfam_70_0001|R70\nMREF\n>DED_hfam_4_0001|R4\nMREF\n>DED_hfam_55_0001|R55\nMREF\n", encoding="ascii")
            ledger_fields = ["reference_id", "accession", "phaded_family_id", "phaded_superfamily"]
            with ledger.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=ledger_fields, delimiter="\t")
                writer.writeheader(); writer.writerows([
                    {"reference_id": "DED_hfam_70_0001", "accession": "R70", "phaded_family_id": "DED_hfam_70", "phaded_superfamily": "extracellular dPHASCL type 2"},
                    {"reference_id": "DED_hfam_4_0001", "accession": "R4", "phaded_family_id": "DED_hfam_4", "phaded_superfamily": "intracellular nPHAMCL"},
                    {"reference_id": "DED_hfam_55_0001", "accession": "R55", "phaded_family_id": "DED_hfam_55", "phaded_superfamily": "extracellular dPHASCL type 1"},
                ])
            output = root / "out"
            summary = module.build_panels(priority, dossier, interpro, candidates, references, ledger, output)
            self.assertEqual(summary["candidate_counts"], {"extracellular_explicit": 1, "phaC_vs_phaZ": 1, "structure_anomaly": 1, "ephaz_competition": 1})
            self.assertEqual((output / "panel_extracellular_explicit.phylo.faa").read_text(encoding="ascii").splitlines()[0], ">candidate|extracellular_explicit|A")
            self.assertIn(">reference_target|extracellular_explicit|DED_hfam_70_0001|R70", (output / "panel_extracellular_explicit.phylo.faa").read_text(encoding="ascii"))
            manifest = json.loads((output / "competition_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "candidate-only")
            self.assertEqual(set(manifest["panels"]), {"extracellular_explicit", "phaC_vs_phaZ", "structure_anomaly", "ephaz_competition"})

    def test_rejects_overlapping_buckets(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            common = "accession\tphaded_family_best\tphaded_superfamily_best\nA\tfamily_DED_hfam_1_x\tintracellular nPHAMCL\n"
            for name in ("priority.tsv", "dossier.tsv", "interpro.tsv"):
                (root / name).write_text(common if name == "priority.tsv" else "accession\tprovisional_label\treview_order\tphaded_family_best\tphaded_superfamily_best\nA\tA_retain_high_priority_architecture\t1\tfamily_DED_hfam_1_x\tintracellular nPHAMCL\n", encoding="utf-8")
            (root / "candidates.faa").write_text(">A\nMAA\n", encoding="ascii")
            (root / "references.faa").write_text(">DED_hfam_1_0001|R\nMREF\n", encoding="ascii")
            (root / "ledger.tsv").write_text("reference_id\taccession\tphaded_family_id\tphaded_superfamily\nDED_hfam_1_0001\tR\tintracellular nPHAMCL\tintracellular nPHAMCL\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing review_stage"):
                module.build_panels(root / "priority.tsv", root / "dossier.tsv", root / "interpro.tsv", root / "candidates.faa", root / "references.faa", root / "ledger.tsv", root / "out")


if __name__ == "__main__":
    unittest.main()
