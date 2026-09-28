"""Tests for the Cys discovery-layer training set builder (2026-09-17, Step 3).

Written before the implementation (test-first).  The builder selects the Cys
superfamily (``intracellular nPHASCL without lipase box``) from the *frozen*
723-row reference ledger, applies an explicit integrity filter, and emits a
training FASTA plus a per-record manifest that records keep/drop and the exact
reason for every dropped record.  It never edits the frozen ledger.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_phaded_cys_discovery_training_set.py"

CYS = "intracellular nPHASCL without lipase box"

LEDGER_HEADER = [
    "reference_id", "accession", "sequence_sha256", "sequence_source",
    "source_database", "source_database_version", "retrieval_date", "organism",
    "taxonomy_id", "phaded_superfamily", "phaded_family_id",
    "reported_localization", "substrate_class", "substrate_detail",
    "experimental_assay", "experimental_result", "evidence_status",
    "positive_negative_control", "catalytic_residues", "lipase_box_or_ahsmg",
    "oxyanion_hole_evidence", "domain_architecture", "primary_doi", "pmid",
    "pmcid", "notes",
]

MAPPING_HEADER = [
    "reference_id", "accession", "phaded_superfamily", "phaded_family_id",
    "evidence_status", "lipase_box_state", "lipase_box_coordinates",
    "lipase_box_x1", "lipase_box_mapping_basis",
    "lipase_box_expectation_from_superfamily", "pfam_accessions",
    "motif_panel_status",
]


def load_module():
    spec = importlib.util.spec_from_file_location(
        "build_phaded_cys_discovery_training_set", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(str(row.get(col, "")) for col in header))
    write_text(path, "\n".join(lines) + "\n")


LONG = "M" + "A" * 219          # 220 aa, passes the >= 200 core minimum
SHORT = "M" + "A" * 99          # 100 aa, below the core minimum
# canonical lipase box G-x1-S-x2-G (Knoll 2009 pattern) embedded in a long chain
WITH_LIPASE_BOX = "M" + "A" * 40 + "GLSHG" + "A" * 175


class TrainingSetBuilderTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.ledger = self.root / "ledger.tsv"
        write_tsv(self.ledger, LEDGER_HEADER, [
            {"reference_id": "DED_hfam_65_0001", "accession": "CAJ92291.1",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_65",
             "reported_localization": "intracellular", "evidence_status": "experimental_positive"},
            {"reference_id": "DED_hfam_61_0001", "accession": "BAF86293.1",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_61",
             "reported_localization": "intracellular", "evidence_status": "annotation_only"},
            {"reference_id": "DED_hfam_61_0043", "accession": "94498690",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_61",
             "reported_localization": "intracellular", "evidence_status": "annotation_only"},
            {"reference_id": "DED_hfam_61_0004", "accession": "ABS67380.1",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_61",
             "reported_localization": "intracellular", "evidence_status": "annotation_only"},
            {"reference_id": "DED_hfam_61_0002", "accession": "SECRETED1.1",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_61",
             "reported_localization": "extracellular", "evidence_status": "annotation_only"},
            {"reference_id": "DED_hfam_57_0038", "accession": "CAJ96855.1",
             "phaded_superfamily": "extracellular dPHASCL type 1", "phaded_family_id": "DED_hfam_57",
             "reported_localization": "extracellular", "evidence_status": "annotation_only"},
        ])
        self.faa = self.root / "ref.faa"
        write_text(self.faa,
                   ">DED_hfam_65_0001|CAJ92291.1\n%s\n" % LONG +
                   ">DED_hfam_61_0001|BAF86293.1\n%s\n" % LONG +
                   ">DED_hfam_61_0043|94498690\n%s\n" % SHORT +
                   ">DED_hfam_61_0004|ABS67380.1\n%s\n" % WITH_LIPASE_BOX +
                   ">DED_hfam_61_0002|SECRETED1.1\n%s\n" % LONG +
                   ">DED_hfam_57_0038|CAJ96855.1\n%s\n" % LONG)
        self.mapping = self.root / "mapping.tsv"
        write_tsv(self.mapping, MAPPING_HEADER, [
            {"reference_id": "DED_hfam_65_0001", "accession": "CAJ92291.1", "phaded_superfamily": CYS,
             "lipase_box_state": "not_detected_pattern", "lipase_box_expectation_from_superfamily":
             "not_expected_for_superfamily", "pfam_accessions": "PF06850"},
            {"reference_id": "DED_hfam_61_0001", "accession": "BAF86293.1", "phaded_superfamily": CYS,
             "lipase_box_state": "not_detected_pattern", "lipase_box_expectation_from_superfamily":
             "not_expected_for_superfamily", "pfam_accessions": "PF06850"},
            {"reference_id": "DED_hfam_61_0043", "accession": "94498690", "phaded_superfamily": CYS,
             "lipase_box_state": "not_detected_pattern", "lipase_box_expectation_from_superfamily":
             "not_expected_for_superfamily", "pfam_accessions": "PF06850"},
            {"reference_id": "DED_hfam_61_0004", "accession": "ABS67380.1", "phaded_superfamily": CYS,
             "lipase_box_state": "supported", "lipase_box_coordinates": "43-43",
             "lipase_box_expectation_from_superfamily": "not_expected_for_superfamily",
             "pfam_accessions": "PF06850"},
            {"reference_id": "DED_hfam_61_0002", "accession": "SECRETED1.1", "phaded_superfamily": CYS,
             "lipase_box_state": "not_detected_pattern", "lipase_box_expectation_from_superfamily":
             "not_expected_for_superfamily", "pfam_accessions": "PF06850"},
        ])

    def tearDown(self):
        self.tmp.cleanup()

    # ---- canonical lipase box detector -----------------------------------
    def test_detect_canonical_lipase_box_finds_gxsxg(self):
        self.assertEqual(self.module.detect_canonical_lipase_box(WITH_LIPASE_BOX), [42])
        self.assertEqual(self.module.detect_canonical_lipase_box(LONG), [])

    # ---- decisions -------------------------------------------------------
    def test_build_manifest_keeps_core_and_records_every_exclusion_reason(self):
        rows, fasta_rows = self.module.build_training_set(
            ledger_rows=self.module.read_tsv(self.ledger),
            mapping_rows=self.module.read_tsv(self.mapping),
            reference_by_id=self.module.load_reference_sequences(self.faa),
            min_core_length=200,
        )
        by_id = {r["reference_id"]: r for r in rows}
        self.assertEqual(len(rows), 5, "only Cys superfamily rows are considered")
        self.assertEqual(by_id["DED_hfam_65_0001"]["include_in_training"], "true")
        self.assertEqual(by_id["DED_hfam_65_0001"]["decision"], "training_core")
        self.assertEqual(by_id["DED_hfam_61_0043"]["include_in_training"], "false")
        self.assertEqual(by_id["DED_hfam_61_0043"]["exclusion_reason"], "short_sequence_below_core_minimum")
        self.assertEqual(by_id["DED_hfam_61_0004"]["include_in_training"], "false")
        self.assertEqual(by_id["DED_hfam_61_0004"]["exclusion_reason"],
                         "lipase_box_annotation_conflict_with_superfamily")
        self.assertEqual(by_id["DED_hfam_61_0002"]["include_in_training"], "false")
        self.assertEqual(by_id["DED_hfam_61_0002"]["exclusion_reason"],
                         "reported_localization_not_intracellular")
        self.assertEqual([r["reference_id"] for r in fasta_rows],
                         ["DED_hfam_65_0001", "DED_hfam_61_0001"])

    def test_missing_signalp_evidence_is_pending_and_retains_sequence(self):
        rows, _ = self.module.build_training_set(
            ledger_rows=self.module.read_tsv(self.ledger),
            mapping_rows=self.module.read_tsv(self.mapping),
            reference_by_id=self.module.load_reference_sequences(self.faa),
            min_core_length=200,
        )
        for row in rows:
            self.assertEqual(row["signalp_evidence"], "pending_reference_layer_has_no_signalp_run")
        kept = [r for r in rows if r["include_in_training"] == "true"]
        self.assertTrue(kept)
        self.assertTrue(all(r["signalp_evidence"].startswith("pending") for r in kept))

    def test_manifest_records_per_record_sequence_sha256(self):
        rows, _ = self.module.build_training_set(
            ledger_rows=self.module.read_tsv(self.ledger),
            mapping_rows=self.module.read_tsv(self.mapping),
            reference_by_id=self.module.load_reference_sequences(self.faa),
            min_core_length=200,
        )
        for row in rows:
            self.assertTrue(row["sequence_sha256"])
            self.assertRegex(row["sequence_sha256"], r"^[0-9a-f]{64}$")
            self.assertTrue(row["sequence_sha256_is_of_protein_sequence"] == "true")

    def test_training_fasta_headers_are_labelled_uncalibrated_discovery_layer(self):
        _, fasta_rows = self.module.build_training_set(
            ledger_rows=self.module.read_tsv(self.ledger),
            mapping_rows=self.module.read_tsv(self.mapping),
            reference_by_id=self.module.load_reference_sequences(self.faa),
            min_core_length=200,
        )
        for row in fasta_rows:
            self.assertIn("discovery_hmm_uncalibrated", row["header"])

    def test_written_outputs_are_self_describing_and_frozen_inputs_untouched(self):
        out_faa = self.root / "out" / "training_set.faa"
        out_manifest = self.root / "out" / "training_set_manifest.tsv"
        out_summary = self.root / "out" / "training_set_summary.json"
        before = {p: p.read_bytes() for p in (self.ledger, self.faa, self.mapping)}
        summary = self.module.write_training_set(
            out_faa, out_manifest, out_summary,
            ledger_rows=self.module.read_tsv(self.ledger),
            mapping_rows=self.module.read_tsv(self.mapping),
            reference_by_id=self.module.load_reference_sequences(self.faa),
            min_core_length=200,
        )
        for path, payload in before.items():
            self.assertEqual(path.read_bytes(), payload, "frozen input was modified")
        self.assertTrue(out_faa.is_file() and out_manifest.is_file() and out_summary.is_file())
        self.assertEqual(summary["training_sequence_count"], 2)
        self.assertEqual(summary["excluded_count"], 3)
        self.assertEqual(summary["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(summary["frozen_ledger_modified"], False)
        payload = json.loads(out_summary.read_text(encoding="utf-8"))
        self.assertEqual(payload["criteria"]["min_core_length"], 200)


if __name__ == "__main__":
    unittest.main()
