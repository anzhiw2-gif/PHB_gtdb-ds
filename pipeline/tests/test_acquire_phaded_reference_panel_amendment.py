import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "acquire_phaded_reference_panel_amendment.py"

PROFILE_FIELDS = ["profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id", "training_accessions", "model_status"]
LEDGER_FIELDS = ["accession", "phaded_family_id", "phaded_superfamily", "evidence_status", "pmid", "primary_doi", "sequence_sha256", "source_database"]
AXIS_FIELDS = [
    "axis_id", "axis_end", "role", "accession", "evidence_status", "anchor_evidence_kind",
    "family_binding_status", "bound_family", "bound_reference_only_profile_id",
    "anchor_justification", "source_doi_or_pmid",
]
AXIS_DEFINITIONS = {
    "axis_definitions": [
        {
            "axis_id": "axis_x",
            "positive_end_qualifying_evidence_kinds": "experimental_positive_assay",
            "negative_end_qualifying_evidence_kinds": "experimental_negative_assay",
            "positive_end_challenge_evidence_kinds": "reported_label",
            "negative_end_challenge_evidence_kinds": "reported_label",
        }
    ]
}


def _load_module():
    spec = importlib.util.spec_from_file_location("acquire_phaded_reference_panel_amendment", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AmendmentTests(unittest.TestCase):
    def test_infer_run_id_uses_dated_run_parent(self):
        module = _load_module()
        self.assertEqual(
            module.infer_run_id("runs/20260915_phaded_reference_panel_acquisition_02/results/amendment"),
            "20260915_phaded_reference_panel_acquisition_02",
        )

    def test_parse_genpept_keeps_accession_and_metadata(self):
        module = _load_module()
        text = (
            "LOCUS       X1          300 aa\n"
            "DEFINITION  hypothetical protein [Example bacterium].\n"
            "ACCESSION   X1\n"
            "VERSION     X1.1\n"
            "REFERENCE   1\n"
            "  PUBMED      123456\n"
            "//\n"
        )
        record = module.parse_genpept(text.splitlines())[0]
        self.assertEqual(record["accession"], "X1.1")
        self.assertIn("hypothetical protein", record["definition"])
        self.assertEqual(record["pmids"], "123456")

    def test_annotation_only_is_never_formal_negative(self):
        module = _load_module()
        self.assertEqual(module.role_for("annotation_only", "same_family"), "challenge_control")
        self.assertEqual(module.eligibility_for("annotation_only", "same_family"), "not_eligible_annotation_only")

    def test_existing_family_binding_is_not_a_new_family_call(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "ledger.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["accession", "phaded_family_id"], delimiter="\t")
                writer.writeheader()
                writer.writerow({"accession": "X1.1", "phaded_family_id": "DED_hfam_2"})
            exact, bases = module.load_ledger(path)
            status, family = module.bind_accession("X1", exact, bases)
            self.assertEqual(status, "resolved_existing_family_versionless")
            self.assertEqual(family, "DED_hfam_2")


class AmendmentAxisCreditTests(unittest.TestCase):
    """Rule #3 at the amendment counter: new axis fields, legacy counters frozen."""

    def _fixture(self, root: Path):
        profiles = root / "profiles.tsv"
        ledger = root / "ledger.tsv"
        panel = root / "axis_panel.tsv"
        definitions = root / "axes.json"
        with profiles.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=PROFILE_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow({"profile_id": "family_DED_hfam_X_deadbeef", "profile_kind": "family", "phaded_superfamily": "sf_x", "phaded_family_id": "DED_hfam_X", "training_accessions": "", "model_status": "reference_only"})
            for index in range(1, 37):
                writer.writerow({"profile_id": f"family_DED_hfam_{index}_pad{index:06d}", "profile_kind": "family", "phaded_superfamily": "sf_other", "phaded_family_id": f"DED_hfam_{index}", "training_accessions": "", "model_status": "reference_only"})
        with ledger.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
        with panel.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=AXIS_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow({"axis_id": "axis_x", "axis_end": "negative", "role": "negative_anchor", "accession": "NEG1", "evidence_status": "experimental_negative", "anchor_evidence_kind": "experimental_negative_assay", "family_binding_status": "unresolved", "bound_family": "", "bound_reference_only_profile_id": "", "anchor_justification": "explicitly not hydrolysing PHB", "source_doi_or_pmid": "10.0000/synthetic"})
            writer.writerow({"axis_id": "axis_x", "axis_end": "positive", "role": "positive_anchor", "accession": "POS1", "evidence_status": "experimental_positive", "anchor_evidence_kind": "experimental_positive_assay", "family_binding_status": "resolved_existing_family", "bound_family": "DED_hfam_X", "bound_reference_only_profile_id": "family_DED_hfam_X_deadbeef", "anchor_justification": "target end anchor", "source_doi_or_pmid": "10.0000/synthetic"})
        definitions.write_text(json.dumps(AXIS_DEFINITIONS), encoding="utf-8")
        return profiles, ledger, panel, definitions

    def _out(self, root: Path, run_id: str) -> Path:
        return root / "runs" / run_id / "results" / "amendment"

    def test_axis_enrichment_is_optional_and_keeps_legacy_columns(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles, ledger, panel, definitions = self._fixture(root)
            base_out = self._out(root, "20260101_test_acquire_base")
            baseline = module.build_amendment(profiles, ledger, [], base_out)
            with (base_out / "profile_calibration_readiness.tsv").open(encoding="utf-8", newline="") as handle:
                baseline_header = handle.readline().strip()
            self.assertNotIn("axis_anchored", baseline_header)

            enriched_out = self._out(root, "20260101_test_acquire_enriched")
            report = module.build_amendment(
                profiles, ledger, [], enriched_out,
                axis_panel=panel, axis_definitions=definitions,
            )
            with (enriched_out / "profile_calibration_readiness.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_id = {row["profile_id"]: row for row in rows}
            target = by_id["family_DED_hfam_X_deadbeef"]
            self.assertEqual(target["external_bound_negative"], "0")
            self.assertEqual(target["family_resolved_negative_status"], "missing")
            self.assertEqual(target["axis_anchored_negative_count"], "1")
            self.assertEqual(target["axis_anchor_basis_status"], "audited")
            others = {row["axis_anchored_negative_count"] for pid, row in by_id.items() if pid != "family_DED_hfam_X_deadbeef"}
            self.assertEqual(others, {"0"})
            self.assertEqual(report["axis_anchored_credit"]["legacy_counters_unchanged"], True)
            self.assertEqual(baseline["summary"]["formal_negative_records"], 0)

    def test_axis_panel_without_audit_columns_is_refused(self):
        module = _load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles, ledger, panel, definitions = self._fixture(root)
            broken = root / "broken_panel.tsv"
            with panel.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            rows[0]["anchor_justification"] = ""
            with broken.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=AXIS_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaises(ValueError):
                module.build_amendment(
                    profiles, ledger, [], self._out(root, "20260101_test_acquire_broken"),
                    axis_panel=broken, axis_definitions=definitions,
                )


if __name__ == "__main__":
    unittest.main()
