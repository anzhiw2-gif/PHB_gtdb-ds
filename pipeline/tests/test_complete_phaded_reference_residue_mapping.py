"""Invariants for the PhaDED reference residue mapping completion (Task 2).

The tests are written before the implementation on purpose: the His/Asp/lid
invariants below are exactly the ones the 2026-09-15 motif panel violated
(His 0/723, Asp restricted to one superfamily, lid derived from a subtype
prior).
"""

import csv
import hashlib
import importlib.util
import inspect
import re
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
LEGACY_AUDIT = REPO / "pipeline" / "scripts" / "audit_phaded_motifs.py"

EXPECTED_RESIDUE = {
    "catalytic_serine": "S",
    "catalytic_histidine": "H",
    "catalytic_aspartate": "D",
    "oxyanion_hole_cysteine": "C",
    "lid_loop_1": None,
    "lid_loop_2": None,
}

PHAZKT_ACCESSION = "AAM63408.1"
PHAZGK13_ACCESSION = "Q51718.1"
TYPE2_ACCESSION = "2D80A"


def _load(name: str, path: Path):
    if not path.is_file():
        raise AssertionError(f"required file is missing: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_module():
    return _load("complete_phaded_reference_residue_mapping", SCRIPT)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SharedFixture(unittest.TestCase):
    """Run the completion once for the whole class."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.tmp = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmp.name) / "results"
        cls.result = cls.module.run_completion(
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
        cls.rows = read_tsv(cls.output_dir / "reference_mapping_manifest.tsv")
        cls.by_accession = {row["accession"]: row for row in cls.rows}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()


class LegacyDefectTests(unittest.TestCase):
    """The pre-existing detector must be demonstrably blind to His, and fixed here."""

    def test_legacy_his_pattern_never_fires_on_any_reference(self):
        legacy = _load("audit_phaded_motifs", LEGACY_AUDIT)
        sequences = legacy.read_fasta(REFERENCE_FASTA)
        cleaned = {key.split("|", 1)[1]: value for key, value in sequences.items()}
        fired = [a for a, s in cleaned.items() if legacy.HIS_PATTERN.search(s.rstrip("*"))]
        self.assertEqual(fired, [], "legacy His pattern was expected to be blind across all 723 references")

    def test_calibrated_his_and_asp_patterns_fire_on_the_8daj_anchor(self):
        module = load_module()
        anchor = module.read_fasta(ANCHOR_FASTA)["8DAJ_1"]
        his = module.detect_literature_pattern(anchor, "catalytic_histidine")
        self.assertEqual(his["state"], "supported")
        self.assertIn(270, module.coordinate_positions(his["coordinates"]))
        asp = module.detect_literature_pattern(anchor, "catalytic_aspartate")
        self.assertEqual(asp["state"], "supported")
        self.assertIn(197, module.coordinate_positions(asp["coordinates"]))


class LiteratureReferenceTests(unittest.TestCase):
    def test_literature_table_records_full_provenance(self):
        rows = read_tsv(LITERATURE)
        self.assertGreaterEqual(len(rows), 10)
        for row in rows:
            self.assertTrue(row["entry_id"], row)
            self.assertTrue(row["residue_role"], row)
            self.assertTrue(row["coordinate_source"], row)
            self.assertTrue(row["source_database"], row)
            self.assertTrue(row["database_version"] if "database_version" in row else row["source_version"], row)
            self.assertTrue(row["retrieval_date"], row)
            self.assertTrue(row["evidence_type"], row)
            for field in ("doi", "pmid", "pmcid"):
                self.assertNotEqual(row[field].strip(), "", f"{field} must be pending at worst: {row['entry_id']}")

    def test_anchor_sequence_hashes_match_recorded_values(self):
        module = load_module()
        anchors = module.read_fasta(ANCHOR_FASTA)
        recorded = {
            row["accession"]: row["anchor_sequence_sha256"]
            for row in read_tsv(LITERATURE)
            if row["anchor_sequence_sha256"] not in ("", "pending")
        }
        panel = {
            key.split("|", 1)[1]: value
            for key, value in module.read_fasta(REFERENCE_FASTA).items()
        }
        checked = 0
        for accession, digest in recorded.items():
            sequence = anchors.get(accession, panel.get(accession))
            self.assertIsNotNone(sequence, f"anchor sequence for {accession} not found")
            self.assertEqual(hashlib.sha256(sequence.encode()).hexdigest(), digest, accession)
            checked += 1
        self.assertGreaterEqual(checked, 4)

    def test_every_anchored_reference_row_repeats_the_recorded_hash(self):
        recorded = {
            row["accession"]: row["anchor_sequence_sha256"]
            for row in read_tsv(LITERATURE)
            if row["anchor_sequence_sha256"] not in ("", "pending")
        }
        self.assertEqual(sorted(recorded), ["2D80A", "8DAJ_1", "AAM63408.1", "Q51718.1"])

    def test_every_literature_coordinate_matches_its_own_anchor_sequence(self):
        module = load_module()
        anchors = module.read_fasta(ANCHOR_FASTA)
        panel = {
            key.split("|", 1)[1]: value
            for key, value in module.read_fasta(REFERENCE_FASTA).items()
        }
        verified = 0
        for row in module.load_literature_reference(LITERATURE):
            coordinate = row["alignment_coordinate"]
            if not coordinate.isdigit():
                continue
            sequence = anchors.get(row["accession"], panel.get(row["accession"]))
            self.assertIsNotNone(sequence, row["accession"])
            expected = EXPECTED_RESIDUE[row["residue_role"]]
            self.assertEqual(sequence[int(coordinate) - 1], expected, row["entry_id"])
            verified += 1
        self.assertGreaterEqual(verified, 8)


class ReferenceMappingTests(SharedFixture):
    def test_manifest_covers_all_723_references(self):
        self.assertEqual(len(self.rows), 723)
        ledger = read_tsv(REFERENCE_LEDGER)
        self.assertEqual({row["accession"] for row in ledger}, set(self.by_accession))

    def test_his_state_is_supported_for_the_8daj_calibrated_superfamily(self):
        supported = [row for row in self.rows if row["his_state"].startswith("supported")]
        self.assertGreater(len(supported), 0, "His must no longer be 0/723")
        for row in supported:
            self.assertTrue(row["his_coordinates"], row["accession"])
            self.assertTrue(row["his_mapping_basis"], row["accession"])

    def test_asp_state_is_supported_for_the_type2_two_d80_reference(self):
        row = self.by_accession[TYPE2_ACCESSION]
        self.assertTrue(row["asp_state"].startswith("supported"), row["asp_state"])
        self.assertIn(101, self.module.coordinate_positions(row["asp_coordinates"]))
        self.assertIn("2D80", row["asp_mapping_basis"])
        self.assertIn(135, self.module.coordinate_positions(row["his_coordinates"]))
        self.assertIn(230, self.module.coordinate_positions(row["oxyanion_hole_coordinates"]))

    def test_anchor_column_transfer_resolves_the_intracellular_mcl_family(self):
        row = self.by_accession[PHAZKT_ACCESSION]
        self.assertTrue(row["his_state"].startswith("supported"), row["his_state"])
        self.assertIn(248, self.module.coordinate_positions(row["his_coordinates"]))
        self.assertIn(221, self.module.coordinate_positions(row["asp_coordinates"]))
        self.assertIn("DED_hfam_4", row["his_mapping_basis"])

    def test_anchor_column_transfer_resolves_the_extracellular_mcl_family(self):
        row = self.by_accession[PHAZGK13_ACCESSION]
        self.assertIn(260, self.module.coordinate_positions(row["his_coordinates"]))
        self.assertIn(228, self.module.coordinate_positions(row["asp_coordinates"]))

    def test_supported_states_always_carry_a_coordinate_and_a_basis(self):
        for row in self.rows:
            for role in ("lipase_box", "his", "asp", "oxyanion_hole", "lid"):
                state = row[f"{role}_state"]
                if state.startswith("supported"):
                    self.assertTrue(row[f"{role}_coordinates"], (row["accession"], role))
                    self.assertTrue(row[f"{role}_mapping_basis"], (row["accession"], role))

    def test_missing_evidence_stays_pending_and_is_never_an_negative(self):
        pending = [row for row in self.rows if row["his_state"].startswith("pending")]
        self.assertGreater(len(pending), 0, "families without a literature anchor must stay pending")
        for row in pending:
            self.assertTrue(
                row["his_mapping_basis"] in {
                    "pending_no_anchor_for_this_role_in_alignment",
                    "pending_no_alignment_row",
                }
                or row["his_mapping_basis"].startswith("pending_gap_at_mapped_column")
                or row["his_mapping_basis"].startswith("pending_anchor_column_not_conserved"),
                row["his_mapping_basis"],
            )
            self.assertNotIn("not_detected", row["his_state"])

    def test_unconserved_transferred_columns_fail_closed_to_pending(self):
        """The type 2 fold is circularly permuted: its columns do not transfer."""
        type2 = [row for row in self.rows if row["phaded_superfamily"] == "extracellular dPHASCL type 2"]
        transferred = [
            row for row in type2
            if row["his_mapping_basis"].startswith("pending_anchor_column_not_conserved")
        ]
        self.assertEqual(len(transferred), 12)
        for row in transferred:
            self.assertEqual(row["his_state"], "pending_reference_annotation")
            self.assertEqual(row["asp_state"], "pending_reference_annotation")
            self.assertIn("conserved_fraction=0.0769", row["his_mapping_basis"])

    def test_sbd_is_bound_from_the_existing_pf06850_domtblout(self):
        expected = set()
        for raw in REFERENCE_DOMTBLOUT.read_text(encoding="utf-8", errors="replace").splitlines():
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) >= 22 and fields[1].split(".", 1)[0] == "PF06850":
                expected.add(fields[3].split("|", 1)[-1])
        supported = {row["accession"] for row in self.rows if row["sbd_state"] == "supported"}
        self.assertEqual(supported, expected)
        self.assertGreater(len(supported), 200)

    def test_linker_is_delimited_from_flanking_domain_coordinates_or_pending(self):
        delimited = [
            row for row in self.rows
            if row["linker_state"] == "delimited_candidate_region_not_profiled"
        ]
        self.assertEqual(len(delimited), 41)
        for row in delimited:
            self.assertTrue(row["linker_coordinates"])
            self.assertTrue(row["linker_mapping_basis"])
            self.assertIn("no_linker_profile_bound", row["linker_mapping_basis"])
        for row in self.rows:
            self.assertIn(row["linker_state"], {
                "delimited_candidate_region_not_profiled",
                "pending_reference_annotation",
            })
        self.assertGreater(
            sum(row["linker_state"] == "pending_reference_annotation" for row in self.rows), 0
        )

    def test_ahsmg_positive_controls_still_pass(self):
        for accession in ("2VTVA", "AAK07742.1"):
            self.assertEqual(self.by_accession[accession]["ahsmg_state"], "supported")

    def test_lipase_box_is_now_detected_independently_of_the_family_prior(self):
        """The 2026-09-15 panel reported 441 by exempting two superfamilies from the scan.

        Independent detection finds 444: three members of the
        "without lipase box" superfamily do carry a chance Gx1Sx2G match.
        """
        supported = [row for row in self.rows if row["lipase_box_state"] == "supported"]
        self.assertEqual(len(supported), 444)
        self.assertEqual(
            sorted({row["lipase_box_state"] for row in self.rows}),
            ["not_detected_pattern", "supported"],
        )


class LidDecircularisationTests(SharedFixture):
    def test_lid_detector_input_excludes_subtype_and_superfamily_labels(self):
        source = inspect.getsource(self.module.detect_lid)
        for forbidden in ("subtype", "superfamily", "phaded_family", "family_id"):
            self.assertNotIn(forbidden, source, f"lid detector must not read {forbidden}")
        parameters = list(inspect.signature(self.module.detect_lid).parameters)
        self.assertIn("sequence", parameters)
        self.assertLessEqual(set(parameters), {"sequence", "signature"})

    def test_lid_detection_call_sites_never_pass_a_label(self):
        source = inspect.getsource(self.module)
        lines = [
            line for line in source.splitlines()
            if "detect_lid(" in line and "def detect_lid" not in line
        ]
        self.assertTrue(lines)
        for line in lines:
            for forbidden in ("subtype", "superfamily", "phaded_family"):
                self.assertNotIn(forbidden, line, line)

    def test_lid_is_detected_on_phazkt_and_not_on_the_extracellular_control(self):
        phazkt = self.by_accession[PHAZKT_ACCESSION]
        self.assertEqual(phazkt["lid_state"], "supported")
        self.assertEqual(phazkt["lid_coordinates"], "34-38;190-195")
        self.assertEqual(phazkt["lid_loop_spacing"], "156")
        gk13 = self.by_accession[PHAZGK13_ACCESSION]
        self.assertEqual(gk13["lid_state"], "not_detected_pattern")
        self.assertEqual(gk13["lid_coordinates"], "")

    def test_lid_state_is_never_derived_from_a_subtype_prior(self):
        allowed = {
            "supported",
            "partial_lid_loop1_only",
            "partial_lid_loop2_only",
            "conflict_lid_loop_spacing",
            "not_detected_pattern",
        }
        for row in self.rows:
            self.assertIn(row["lid_state"], allowed, row["accession"])
        self.assertEqual(
            sorted({row["lid_expectation_from_superfamily"] for row in self.rows}),
            ["expected_for_superfamily", "not_expected_for_superfamily"],
        )

    def test_lid_detector_flags_broken_loop_spacing(self):
        module = load_module()
        sequence = "M" + "A" * 5 + "FNGIG" + "A" * 400 + "YYWQLF" + "A" * 10
        result = module.detect_lid(sequence)
        self.assertEqual(result["state"], "conflict_lid_loop_spacing")
        self.assertEqual(result["loop_spacing"], "405")

    def test_lid_detected_set_is_reported_independently_of_the_prior(self):
        detected = {row["accession"] for row in self.rows if row["lid_state"] == "supported"}
        self.assertEqual(len(detected), 39)
        for row in self.rows:
            if row["accession"] in detected:
                self.assertEqual(row["lid_expectation_from_superfamily"], "expected_for_superfamily")


class CalibrationReportTests(SharedFixture):
    def test_calibration_report_compares_before_and_after(self):
        report_path = self.output_dir / "implementation_calibration_report.json"
        self.assertTrue(report_path.is_file())
        report = self.module.json.loads(report_path.read_text(encoding="utf-8"))
        before = report["criteria"]["his_state"]["before"]
        after = report["criteria"]["his_state"]["after"]
        self.assertEqual(before["supported"], 0)
        self.assertGreater(after["supported"], 0)
        self.assertEqual(report["criteria"]["asp_state"]["before"]["supported"], 25)
        self.assertGreater(report["criteria"]["asp_state"]["after"]["supported"], 25)
        self.assertEqual(report["criteria"]["lid_state"]["before"]["independent_detection_supported"], 0)
        self.assertEqual(report["criteria"]["lid_state"]["after"]["supported"], 39)
        self.assertEqual(report["criteria"]["ahsmg_state"]["after"]["supported"], 2)
        self.assertEqual(report["reference_count"], 723)

    def test_calibration_report_records_the_defect_and_remaining_pending(self):
        report = self.module.json.loads(
            (self.output_dir / "implementation_calibration_report.json").read_text(encoding="utf-8")
        )
        defect = report["defects_corrected"]
        self.assertTrue(any(item["defect_id"] == "DEFECT_A_LID_CIRCULARITY" for item in defect))
        self.assertTrue(any(item["defect_id"] == "DEFECT_B_HIS_PATTERN_NEVER_FIRED" for item in defect))
        self.assertIn("linker_state", report["criteria"])
        self.assertGreater(report["criteria"]["linker_state"]["after"]["pending_reference_annotation"], 0)
        self.assertEqual(report["phenotype_boundary_holds"], True)


if __name__ == "__main__":
    unittest.main()
