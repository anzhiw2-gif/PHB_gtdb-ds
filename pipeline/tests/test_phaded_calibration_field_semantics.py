"""Field-name semantics for the Task 2 implementation calibration report (P4).

The Phase 1 verification report (P4) found that
``criteria.<criterion>.after.independent_detection_supported`` carried the *total*
supported count, even though part of that count is not an independent detection at
all: it is inferred by transferring a literature-anchored DED alignment column from
an anchor row to another row (conservation gate >= 0.5).

These tests pin the corrected contract of the report:

* ``resolved_supported`` is the total supported count for a criterion;
* ``independent_detection_supported`` excludes rows that rest only on an
  anchor-column transfer;
* ``transferred_from_anchor_column`` is reported as its own field and never
  silently folded into the total under an "independent" name.

The tests are written before the implementation change on purpose: they fail against
the Phase 1 report schema.
"""

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "pipeline" / "scripts" / "complete_phaded_reference_residue_mapping.py"
RUN = REPO / "runs" / "20260917_phaded_reference_residue_mapping_01"
LITERATURE = RUN / "inputs" / "literature_residue_reference.tsv"
ANCHOR_FASTA = RUN / "inputs" / "literature_anchor_sequences.faa"

REFERENCE_FASTA = REPO / "runs" / "20260911_phaded_pfam_architecture_01" / "inputs" / "phaded_reference.faa"
REFERENCE_LEDGER = REPO / "runs" / "20260911_phaded_pfam_architecture_01" / "inputs" / "phaded_reference_ledger.tsv"
REFERENCE_DOMTBLOUT = REPO / "runs" / "20260911_phaded_pfam_architecture_01" / "results" / "raw" / "reference_pfam.domtblout"
ALIGNMENT_DIR = REPO / "runs" / "20260909_phaded_reference_evidence_01" / "inputs" / "raw" / "ded"
ALIGNMENT_MANIFEST = REPO / "runs" / "20260915_phaded_motif_reconciliation_01" / "results" / "reference_alignment_manifest.tsv"
LEGACY_PANEL = REPO / "runs" / "20260915_phaded_motif_reconciliation_01" / "results" / "reference_motif_panel.tsv"

TRANSFER_PREFIX = "anchor_column_transfer:"
ANCHOR_SELF_PREFIX = "literature_anchor:"

# criterion -> manifest column that records how the state was resolved
BASIS_COLUMNS = {
    "his_state": "his_mapping_basis",
    "asp_state": "asp_mapping_basis",
    "oxyanion_hole_state": "oxyanion_hole_mapping_basis",
}

CRITERIA = (
    "lipase_box_state", "catalytic_ser_cys_state", "his_state", "asp_state",
    "oxyanion_hole_state", "ahsmg_state", "sbd_state", "linker_state", "lid_state",
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "complete_phaded_reference_residue_mapping", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CalibrationFieldSemanticsTests(unittest.TestCase):
    """The report must not label transferred evidence as independent detection."""

    @classmethod
    def setUpClass(cls):
        module = load_module()
        cls.module = module
        cls.tmp = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmp.name) / "results"
        cls.result = module.run_completion(
            reference_fasta=REFERENCE_FASTA,
            reference_ledger=REFERENCE_LEDGER,
            reference_domtblout=REFERENCE_DOMTBLOUT,
            alignment_dir=ALIGNMENT_DIR,
            alignment_manifest=ALIGNMENT_MANIFEST,
            literature_reference=LITERATURE,
            anchor_fasta=ANCHOR_FASTA,
            legacy_panel=LEGACY_PANEL,
            output_dir=cls.output_dir,
            run_id="20260917_phaded_reference_residue_mapping_01",
        )
        cls.report = cls.result["calibration_report"]
        cls.rows = read_tsv(cls.output_dir / "reference_mapping_manifest.tsv")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _expected_transfer_split(self, criterion: str) -> dict:
        basis_column = BASIS_COLUMNS.get(criterion)
        if basis_column is None:
            return {"transferred": 0, "anchor_self": 0, "resolved": 0}
        support = 0
        transferred = 0
        anchor_self = 0
        for row in self.rows:
            if not row[criterion].startswith("supported"):
                continue
            support += 1
            basis = row[basis_column]
            if basis.startswith(TRANSFER_PREFIX):
                transferred += 1
            elif basis.startswith(ANCHOR_SELF_PREFIX):
                anchor_self += 1
        return {"transferred": transferred, "anchor_self": anchor_self, "resolved": support}

    def test_every_criterion_reports_resolved_supported_as_the_total(self):
        for criterion in CRITERIA:
            after = self.report["criteria"][criterion]["after"]
            self.assertIn("resolved_supported", after, criterion)
            self.assertIn("transferred_from_anchor_column", after, criterion)
            self.assertEqual(after["resolved_supported"], after["supported"], criterion)
            self.assertEqual(
                after["resolved_supported"],
                sum(
                    row[criterion].startswith("supported") for row in self.rows
                ),
                criterion,
            )

    def test_independent_detection_no_longer_carries_the_total(self):
        """Where anchored columns were transferred, the name must not mean the total."""
        for criterion in BASIS_COLUMNS:
            after = self.report["criteria"][criterion]["after"]
            expected = self._expected_transfer_split(criterion)
            self.assertEqual(
                after["transferred_from_anchor_column"],
                expected["transferred"],
                criterion,
            )
            self.assertEqual(
                after["independent_detection_supported"],
                after["resolved_supported"] - expected["transferred"],
                criterion,
            )
            if expected["transferred"]:
                self.assertNotEqual(
                    after["independent_detection_supported"],
                    after["resolved_supported"],
                    f"{criterion}: transferred evidence is still reported as independent detection",
                )

    def test_his_and_asp_do_rest_on_anchor_column_transfers(self):
        """The two criteria that actually carry transferred support must say so."""
        for criterion in ("his_state", "asp_state"):
            after = self.report["criteria"][criterion]["after"]
            self.assertGreater(after["transferred_from_anchor_column"], 0, criterion)
            self.assertLess(
                after["independent_detection_supported"],
                after["resolved_supported"],
                criterion,
            )

    def test_oxyanion_carries_no_transferred_coordinate(self):
        """PDB 2D80's oxyanion column is not conserved, so nothing is transferred for it.

        The single non-pattern oxyanion support is the 2D80 anchor row carrying its own
        recorded coordinate, which is why it counts as independent of any transfer.
        """
        after = self.report["criteria"]["oxyanion_hole_state"]["after"]
        self.assertEqual(after["transferred_from_anchor_column"], 0)
        self.assertEqual(after["anchor_self_literature_coordinate"], 1)
        self.assertEqual(
            after["independent_detection_supported"], after["resolved_supported"]
        )

    def test_transfer_split_is_broken_out_by_evidence_source(self):
        for criterion in BASIS_COLUMNS:
            after = self.report["criteria"][criterion]["after"]
            breakdown = after["supported_evidence_breakdown"]
            expected = self._expected_transfer_split(criterion)
            self.assertEqual(breakdown["anchor_column_transfer"], expected["transferred"])
            self.assertEqual(
                breakdown["anchor_self_literature_coordinate"], expected["anchor_self"]
            )
            self.assertEqual(
                after["supported_evidence_breakdown"]["direct_row_evidence"],
                expected["resolved"] - expected["transferred"] - expected["anchor_self"],
            )
            self.assertEqual(
                sum(breakdown.values()),
                after["resolved_supported"],
                f"{criterion}: evidence breakdown does not sum to the resolved total",
            )

    def test_non_transfer_criteria_report_a_zero_transfer_split(self):
        for criterion in CRITERIA:
            if criterion in BASIS_COLUMNS:
                continue
            after = self.report["criteria"][criterion]["after"]
            self.assertEqual(after["transferred_from_anchor_column"], 0, criterion)
            self.assertEqual(
                after["independent_detection_supported"],
                after["resolved_supported"],
                criterion,
            )

    def test_report_documents_the_field_semantics(self):
        semantics = self.report["field_semantics"]
        self.assertIn("resolved_supported", semantics["fields"])
        self.assertIn("independent_detection_supported", semantics["fields"])
        self.assertIn("transferred_from_anchor_column", semantics["fields"])
        self.assertIn("transfer", semantics["fields"]["independent_detection_supported"].lower())

    def test_frozen_phase1_report_is_untouched_by_the_rename(self):
        """The Phase 1 report file stays as evidence; this run only changes future output."""
        frozen = RUN / "results" / "implementation_calibration_report.json"
        self.assertTrue(frozen.is_file())
        import json

        payload = json.loads(frozen.read_text(encoding="utf-8"))
        self.assertIn(
            "independent_detection_supported",
            payload["criteria"]["his_state"]["after"],
        )


if __name__ == "__main__":
    unittest.main()
