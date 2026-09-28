"""Tests for the intracellular nPHASCL (no lipase box) Cys-type annotation layer.

Written before the implementation (test-first). The module under test builds a
*candidate-layer annotation evidence* table only: it fits no HMM, opens no
registry gate, changes no ``subtype_call`` and writes nothing outside its own
output directory.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_phaded_inphascl_cys_evidence.py"

CANDIDATE_HEADER = [
    "accession",
    "genome",
    "candidate_prior_superfamily",
    "sbd_pf06850_binding_state",
    "sbd_pf06850_coordinates",
    "lipase_box_state",
]
# the reference manifest is read for its Pfam/lipase-box evidence columns only
REFERENCE_MANIFEST_HEADER = [
    "accession",
    "reference_id",
    "phaded_superfamily",
    "phaded_family_id",
    "pfam_accessions",
    "lipase_box_state",
]


def load_module():
    spec = importlib.util.spec_from_file_location(
        "annotate_phaded_inphascl_cys_evidence", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_fasta(path, records):
    lines = []
    for name, sequence in records:
        lines.append(f">{name}")
        lines.append(sequence)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def synthetic_inputs(root, *, candidate_rows, reference_rows, candidate_seqs, reference_seqs):
    candidate_evidence = root / "candidate.tsv"
    reference_manifest = root / "reference.tsv"
    candidate_fasta = root / "candidate.faa"
    reference_fasta = root / "reference.faa"
    write_tsv(candidate_evidence, CANDIDATE_HEADER, candidate_rows)
    write_tsv(reference_manifest, REFERENCE_MANIFEST_HEADER, reference_rows)
    write_fasta(candidate_fasta, candidate_seqs)
    write_fasta(reference_fasta, reference_seqs)
    return candidate_evidence, reference_manifest, candidate_fasta, reference_fasta


class CysEvidenceClassTests(unittest.TestCase):
    """The three primitive states must combine into an exhaustive, disjoint set."""

    def test_class_is_derived_from_the_three_primitive_states(self):
        module = load_module()
        for pf_state in module.PF06850_BINDING_STATES:
            for lipase_state in module.LIPASE_BOX_STATES:
                for val_state in module.VAL_CYS_STATES:
                    label = module.cys_evidence_class(pf_state, lipase_state, val_state)
                    self.assertIn(label, module.CYS_EVIDENCE_CLASSES, label)
        self.assertEqual(
            len(module.CYS_EVIDENCE_CLASSES),
            len(module.PF06850_BINDING_STATES)
            * len(module.LIPASE_BOX_STATES)
            * len(module.VAL_CYS_STATES),
        )
        self.assertEqual(len(set(module.CYS_EVIDENCE_CLASSES)), len(module.CYS_EVIDENCE_CLASSES))

    def test_unknown_state_fails_closed(self):
        module = load_module()
        with self.assertRaises(ValueError):
            module.cys_evidence_class("maybe", "not_detected_pattern", "present")
        with self.assertRaises(ValueError):
            module.cys_evidence_class("detected", "not_detected_in_tested_pfam", "present")
        with self.assertRaises(ValueError):
            module.cys_evidence_class("detected", "not_detected_pattern", "yes")

    def test_primary_evidence_requires_pf06850_and_absent_lipase_box(self):
        module = load_module()
        self.assertTrue(
            module.is_cys_type_candidate_evidence("detected", "not_detected_pattern")
        )
        self.assertFalse(
            module.is_cys_type_candidate_evidence("detected", "supported")
        )
        self.assertFalse(
            module.is_cys_type_candidate_evidence("not_detected_in_tested_pfam", "not_detected_pattern")
        )
        self.assertFalse(
            module.is_cys_type_candidate_evidence("not_tested_no_pfam_evidence", "not_detected_pattern")
        )
        # an unknown state is never silently treated as "not a positive"
        with self.assertRaises(ValueError):
            module.is_cys_type_candidate_evidence("detected", "not_tested")


class ValCysDescriptorTests(unittest.TestCase):
    """Knoll 2009: the catalytic cysteine carries an (almost always valine) -1 residue."""

    def test_dyad_positions_are_one_based_and_valine_specific(self):
        module = load_module()
        # V at 2 -> C at 3 ; V at 8 -> C at 9 ; the C at 1 has no preceding residue
        self.assertEqual(module.val_cys_cysteine_positions("AVCAVVAVCA"), [3, 9])
        self.assertEqual(module.val_cys_cysteine_positions("VCA"), [2])
        self.assertEqual(module.val_cys_cysteine_positions("ACV"), [])

    def test_generic_hydrophobic_minus_one_is_tracked_separately(self):
        module = load_module()
        self.assertEqual(module.hydrophobic_minus1_cysteine_positions("ILC"), [3])
        self.assertEqual(module.hydrophobic_minus1_cysteine_positions("DEC"), [])
        self.assertEqual(module.hydrophobic_minus1_cysteine_positions("AVCILC"), [3, 6])
        # valine is a subset of the hydrophobic set, so a pure V-C dyad must agree
        sequence = "AVC"
        self.assertEqual(
            module.val_cys_cysteine_positions(sequence),
            module.hydrophobic_minus1_cysteine_positions(sequence),
        )

    def test_descriptor_takes_only_a_sequence_so_it_cannot_be_circular(self):
        import inspect

        module = load_module()
        for name in (
            "val_cys_cysteine_positions",
            "hydrophobic_minus1_cysteine_positions",
        ):
            parameters = list(inspect.signature(getattr(module, name)).parameters)
            self.assertEqual(len(parameters), 1, name)
            self.assertNotIn("superfamily", parameters[0], name)
            self.assertNotIn("family", parameters[0], name)
            self.assertNotIn("subtype", parameters[0], name)
        # the classification takes the three measured states, never a family label
        parameters = list(inspect.signature(module.cys_evidence_class).parameters)
        self.assertEqual(
            parameters, ["pf06850_binding_state", "lipase_box_state", "val_cys_state"]
        )

    def test_sequence_with_no_cysteine_yields_absent_state(self):
        module = load_module()
        self.assertEqual(module.val_cys_cysteine_positions("MKLPQ"), [])
        self.assertEqual(module.hydrophobic_minus1_cysteine_positions("MKLPQ"), [])

    def test_empty_sequence_fails_closed(self):
        module = load_module()
        for bad in ("", "   ", None):
            with self.assertRaises(ValueError):
                module.val_cys_cysteine_positions(bad)


class AnnotationRunTests(unittest.TestCase):
    def _fixtures(self, root):
        candidate_rows = [
            ["p1", "G1", "extracellular dPHASCL type 1", "not_detected_in_tested_pfam", "", "supported"],
            ["p2", "G1", "", "detected", "211-412", "not_detected_pattern"],
            ["p3", "G1", "", "detected", "150-330", "not_detected_pattern"],
            ["p4", "G2", "", "detected", "90-270", "supported"],
            ["p5", "G2", "intracellular nPHAMCL", "not_tested_no_pfam_evidence", "", "not_detected_pattern"],
        ]
        reference_rows = [
            ["RA1", "RA1_0001", "intracellular nPHASCL without lipase box", "DED_hfam_65", "PF06850", "not_detected_pattern"],
            ["RA2", "RA2_0001", "intracellular nPHASCL without lipase box", "DED_hfam_61", "PF06850", "not_detected_pattern"],
            ["RB1", "RB1_0001", "extracellular dPHASCL type 1", "DED_hfam_52", "PF10503", "supported"],
            ["RB2", "RB2_0001", "extracellular dPHASCL type 2", "DED_hfam_70", "", "supported"],
        ]
        candidate_seqs = [
            ("p1", "MAVCGT"),
            ("p2", "MKVCGT"),
            ("p3", "MKLAGT"),
            ("p4", "MKVCGT"),
            ("p5", "MKVCGT"),
        ]
        reference_seqs = [
            ("RA1_0001|RA1", "MAVCGT"),
            ("RA2_0001|RA2", "MIVCGT"),
            ("RB1_0001|RB1", "MAISGT"),
            ("RB2_0001|RB2", "MAISGT"),
        ]
        return synthetic_inputs(
            root,
            candidate_rows=candidate_rows,
            reference_rows=reference_rows,
            candidate_seqs=candidate_seqs,
            reference_seqs=reference_seqs,
        )

    def test_end_to_end_outputs_are_self_consistent_and_inputs_untouched(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            before = [path.read_bytes() for path in inputs]
            manifest = module.annotate(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=root / "results",
                run_id="test_run",
            )
            after = [path.read_bytes() for path in inputs]
            self.assertEqual(before, after)

            rows = (root / "results" / "cys_candidate_evidence.tsv").read_text(
                encoding="utf-8"
            ).splitlines()
            header = rows[0].split("\t")
            self.assertEqual(rows[0], "\t".join(module.CANDIDATE_COLUMNS))
            self.assertEqual([line.split("\t")[0] for line in rows[1:]], ["p1", "p2", "p3", "p4", "p5"])
            by_accession = {line.split("\t")[0]: dict(zip(header, line.split("\t"))) for line in rows[1:]}

            # p2 (detected + no lipase box + V-C) and p3 (same without V-C) are both declared
            # positives: the val_cys descriptor is declared non-load-bearing
            self.assertEqual(by_accession["p2"]["cys_type_candidate_evidence"], "true")
            self.assertEqual(by_accession["p3"]["cys_type_candidate_evidence"], "true")
            self.assertNotEqual(
                by_accession["p2"]["cys_evidence_class"], by_accession["p3"]["cys_evidence_class"]
            )
            self.assertEqual(by_accession["p4"]["cys_type_candidate_evidence"], "false")
            self.assertEqual(by_accession["p1"]["cys_type_candidate_evidence"], "false")
            self.assertEqual(by_accession["p1"]["val_cys_dyad_state"], "present")
            self.assertEqual(by_accession["p3"]["val_cys_dyad_state"], "absent")
            self.assertEqual(by_accession["p5"]["val_cys_dyad_state"], "present")
            self.assertEqual(by_accession["p5"]["cys_type_candidate_evidence"], "false")
            self.assertEqual(by_accession["p2"]["sbd_pf06850_coordinates"], "211-412")
            self.assertEqual(by_accession["p2"]["evidence_layer"], module.EVIDENCE_LAYER)

            for accession, row in by_accession.items():
                self.assertEqual(row["evidence_layer"], module.EVIDENCE_LAYER, accession)
                self.assertEqual(row["family_call_made"], "false", accession)
                self.assertEqual(row["new_family_call_made"], "false", accession)
                self.assertEqual(row["hmm_score_present"], "false", accession)
                self.assertEqual(row["subtype_call_impact"], "none", accession)
                self.assertEqual(row["evidence_boundary"], module.PHENOTYPE_BOUNDARY, accession)

            self.assertEqual(manifest["candidate_count"], 5)
            self.assertEqual(manifest["cys_type_candidate_evidence_count"], 2)
            self.assertEqual(manifest["subtype_call_rows_changed"], 0)
            self.assertEqual(manifest["hmm_fitted"], False)
            self.assertEqual(manifest["run_id"], "test_run")
            self.assertEqual(manifest["status"], "completed_candidate_only")
            self.assertEqual(
                sum(manifest["cys_evidence_class_counts"].values()),
                manifest["candidate_count"],
            )
            self.assertEqual(
                sum(manifest["pf06850_binding_state_counts"].values()),
                manifest["candidate_count"],
            )
            self.assertEqual(manifest["outputs"]["cys_candidate_evidence.tsv"]["size"],
                             (root / "results" / "cys_candidate_evidence.tsv").stat().st_size)

            summary = json.loads((root / "results" / "cys_decision_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["candidate_count"], 5)

    def test_reference_contrast_counts_and_rates(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            module.annotate(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=root / "results",
                run_id="test_run",
            )
            lines = (root / "results" / "cys_reference_contrast.tsv").read_text(
                encoding="utf-8"
            ).splitlines()
            header = lines[0].split("\t")
            self.assertEqual(lines[0], "\t".join(module.REFERENCE_CONTRAST_COLUMNS))
            by_superfamily = {
                dict(zip(header, line.split("\t")))["phaded_superfamily"]: dict(
                    zip(header, line.split("\t"))
                )
                for line in lines[1:]
            }
            cys = by_superfamily["intracellular nPHASCL without lipase box"]
            self.assertEqual(cys["reference_count"], "2")
            self.assertEqual(cys["pf06850_binding_count"], "2")
            self.assertEqual(cys["pf06850_binding_rate"], "1.000000")
            self.assertEqual(cys["cys_evidence_classifier_count"], "2")
            extra = by_superfamily["extracellular dPHASCL type 1"]
            self.assertEqual(extra["pf06850_binding_count"], "0")
            self.assertEqual(extra["cys_evidence_classifier_count"], "0")

    def test_falsification_gate_withholds_evidence_when_pf06850_is_not_specific(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            # an extracellular reference that carries PF06850 must trip the pre-registered gate
            write_tsv(
                inputs[1],
                REFERENCE_MANIFEST_HEADER,
                [
                    ["RA1", "RA1_0001", "intracellular nPHASCL without lipase box", "DED_hfam_65", "PF06850", "not_detected_pattern"],
                    ["RB1", "RB1_0001", "extracellular dPHASCL type 1", "DED_hfam_52", "PF06850", "not_detected_pattern"],
                ],
            )
            manifest = module.annotate(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=root / "results",
                run_id="test_run",
            )
            gates = manifest["falsification_gates"]
            self.assertEqual(gates["pf06850_reference_specificity"]["outcome"], "failed")
            self.assertEqual(manifest["status"], "falsification_triggered")
            self.assertEqual(manifest["cys_type_candidate_evidence_count"], 0)
            rows = (root / "results" / "cys_candidate_evidence.tsv").read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertIn(module.WITHHELD_CLASS, rows[1])

    def test_reference_sensitivity_gate_fails_below_the_declared_floor(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            # keep the Cys superfamily but strip PF06850 from one of its two references
            write_tsv(
                inputs[1],
                REFERENCE_MANIFEST_HEADER,
                [
                    ["RA1", "RA1_0001", "intracellular nPHASCL without lipase box", "DED_hfam_65", "PF06850", "not_detected_pattern"],
                    ["RA2", "RA2_0001", "intracellular nPHASCL without lipase box", "DED_hfam_61", "", "not_detected_pattern"],
                ],
            )
            manifest = module.annotate(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=root / "results",
                run_id="test_run",
            )
            gate = manifest["falsification_gates"]["pf06850_reference_sensitivity"]
            self.assertEqual(gate["outcome"], "failed")
            self.assertGreaterEqual(gate["floor"], 0.9)
            self.assertEqual(manifest["status"], "falsification_triggered")

    def test_fail_closed_on_missing_column_duplicate_accession_and_missing_sequence(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)

            broken = root / "broken.tsv"
            broken.write_text("accession\tgenome\np1\tG1\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                module.read_candidate_evidence(broken)

            duplicate = root / "dup.tsv"
            write_tsv(
                duplicate,
                CANDIDATE_HEADER,
                [
                    ["p1", "G1", "", "detected", "", "not_detected_pattern"],
                    ["p1", "G1", "", "detected", "", "not_detected_pattern"],
                ],
            )
            with self.assertRaises(ValueError):
                module.read_candidate_evidence(duplicate)

            unknown = root / "unknown.tsv"
            write_tsv(
                unknown,
                CANDIDATE_HEADER,
                [["p1", "G1", "", "probably", "", "not_detected_pattern"]],
            )
            with self.assertRaises(ValueError):
                module.read_candidate_evidence(unknown)

            short_fasta = root / "short.faa"
            write_fasta(short_fasta, [("p1", "MAVCGT")])
            with self.assertRaises(ValueError):
                module.annotate(
                    candidate_evidence=inputs[0],
                    reference_manifest=inputs[1],
                    candidate_fasta=short_fasta,
                    reference_fasta=inputs[3],
                    output_dir=root / "results_missing_seq",
                    run_id="test_run",
                )

    def test_module_never_touches_the_formal_registry_or_config(self):
        module = load_module()
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("formal_scan_models", source)
        self.assertNotIn("git commit", source)
        self.assertNotIn("git push", source)
        self.assertIn("subtype_call_impact", module.CANDIDATE_COLUMNS)
        self.assertNotIn("subtype_call", module.CANDIDATE_COLUMNS)

    def test_main_returns_non_zero_on_a_triggered_gate(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            write_tsv(
                inputs[1],
                REFERENCE_MANIFEST_HEADER,
                [
                    ["RA1", "RA1_0001", "intracellular nPHASCL without lipase box", "DED_hfam_65", "PF06850", "not_detected_pattern"],
                    ["RB1", "RB1_0001", "extracellular dPHASCL type 1", "DED_hfam_52", "PF06850", "not_detected_pattern"],
                ],
            )
            exit_code = module.main(
                [
                    "--candidate-evidence", str(inputs[0]),
                    "--reference-manifest", str(inputs[1]),
                    "--candidate-fasta", str(inputs[2]),
                    "--reference-fasta", str(inputs[3]),
                    "--output-dir", str(root / "results"),
                    "--run-id", "test_run",
                ]
            )
            self.assertNotEqual(exit_code, 0)


class TypingDescriptorPassTests(unittest.TestCase):
    """The exploratory VCQ descriptor pass must not touch the pre-registered artifact."""

    def _fixtures(self, root):
        candidate_rows = [
            ["p1", "G1", "", "detected", "211-412", "not_detected_pattern"],
            ["p2", "G1", "", "detected", "150-330", "not_detected_pattern"],
            ["p3", "G2", "", "not_detected_in_tested_pfam", "", "supported"],
        ]
        reference_rows = [
            ["RA1", "RA1_0001", "intracellular nPHASCL without lipase box", "DED_hfam_65", "PF06850", "not_detected_pattern"],
            ["RA2", "RA2_0001", "intracellular nPHASCL without lipase box", "DED_hfam_61", "PF06850", "not_detected_pattern"],
            ["RB1", "RB1_0001", "extracellular dPHASCL type 1", "DED_hfam_52", "PF10503", "supported"],
        ]
        candidate_seqs = [("p1", "MAVCQGT"), ("p2", "MAVCAGT"), ("p3", "MAISGT")]
        reference_seqs = [
            ("RA1_0001|RA1", "MAVCQGT"),
            ("RA2_0001|RA2", "MIVCQGT"),
            ("RB1_0001|RB1", "MAVCQGT"),
        ]
        return synthetic_inputs(
            root,
            candidate_rows=candidate_rows,
            reference_rows=reference_rows,
            candidate_seqs=candidate_seqs,
            reference_seqs=reference_seqs,
        )

    def test_vcq_positions_are_one_based_and_glutamine_specific(self):
        module = load_module()
        self.assertEqual(module.val_cys_gln_cysteine_positions("AVCQ"), [3])
        self.assertEqual(module.val_cys_gln_cysteine_positions("AVCQQ"), [3])
        self.assertEqual(module.val_cys_gln_cysteine_positions("AVC"), [])
        self.assertEqual(module.val_cys_gln_cysteine_positions("AVCA"), [])
        self.assertEqual(module.val_cys_gln_cysteine_positions("DKVCQ"), [4])
        self.assertEqual(module.val_cys_gln_cysteine_positions("VCQ"), [2])
        with self.assertRaises(ValueError):
            module.val_cys_gln_cysteine_positions("")

    def test_descriptor_pass_requires_the_preregistered_artifact_and_leaves_it_byte_identical(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            output = root / "results"
            preregistered_missing = module.annotate_typing_descriptor(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=output,
                run_id="test_run",
                _allow_without_preregistered_artifact=True,
            )
            self.assertEqual(preregistered_missing["preregistered_artifact_checked"], False)
            with self.assertRaises(FileNotFoundError):
                module.annotate_typing_descriptor(
                    candidate_evidence=inputs[0],
                    reference_manifest=inputs[1],
                    candidate_fasta=inputs[2],
                    reference_fasta=inputs[3],
                    output_dir=output,
                    run_id="test_run",
                )

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            output = root / "results"
            module.annotate(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=output,
                run_id="test_run",
            )
            artifact = output / "cys_candidate_evidence.tsv"
            before = artifact.read_bytes()
            manifest = module.annotate_typing_descriptor(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=output,
                run_id="test_run",
            )
            self.assertEqual(artifact.read_bytes(), before)
            self.assertEqual(manifest["preregistered_artifact_checked"], True)
            self.assertEqual(manifest["preregistered_artifact_sha256"], module._sha256(artifact))
            self.assertEqual(manifest["evidence_role"], module.EXPLORATORY_ROLE)
            self.assertEqual(manifest["held_out_set"], "none")

            lines = (output / "cys_typing_descriptor_evidence.tsv").read_text(
                encoding="utf-8"
            ).splitlines()
            header = lines[0].split("\t")
            self.assertEqual(lines[0], "\t".join(module.DESCRIPTOR_COLUMNS))
            by_accession = {
                line.split("\t")[0]: dict(zip(header, line.split("\t"))) for line in lines[1:]
            }
            self.assertEqual(sorted(by_accession), ["p1", "p2", "p3"])
            self.assertEqual(by_accession["p1"]["val_cys_gln_state"], "present")
            self.assertEqual(by_accession["p2"]["val_cys_gln_state"], "absent")
            self.assertEqual(by_accession["p2"]["val_cys_dyad_state"], "present")
            self.assertEqual(by_accession["p3"]["val_cys_gln_state"], "absent")
            self.assertEqual(by_accession["p1"]["family_call_made"], "false")
            self.assertEqual(by_accession["p1"]["evidence_role"], module.EXPLORATORY_ROLE)

            discovery = json.loads(
                (output / "cys_typing_descriptor_discovery.json").read_text(encoding="utf-8")
            )
            cys = discovery["descriptors"]["val_cys_gln_tripeptide"]["cys_superfamily"]
            other = discovery["descriptors"]["val_cys_gln_tripeptide"]["other_superfamilies"]
            self.assertEqual(cys["hits"], 2)
            self.assertEqual(cys["reference_count"], 2)
            self.assertEqual(other["hits"], 1)
            self.assertEqual(other["reference_count"], 1)
            self.assertEqual(other["hit_accessions"], ["RB1"])
            self.assertEqual(discovery["in_sample"], True)
            self.assertEqual(discovery["held_out_set"], "none")
            self.assertEqual(manifest["candidate_count"], 3)
            self.assertEqual(manifest["val_cys_gln_present_count"], 1)

    def test_descriptor_pass_fails_closed_on_missing_sequence(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            inputs = self._fixtures(root)
            module.annotate(
                candidate_evidence=inputs[0],
                reference_manifest=inputs[1],
                candidate_fasta=inputs[2],
                reference_fasta=inputs[3],
                output_dir=root / "results",
                run_id="test_run",
            )
            short = root / "short.faa"
            write_fasta(short, [("p1", "MAVCQGT")])
            with self.assertRaises(ValueError):
                module.annotate_typing_descriptor(
                    candidate_evidence=inputs[0],
                    reference_manifest=inputs[1],
                    candidate_fasta=short,
                    reference_fasta=inputs[3],
                    output_dir=root / "results",
                    run_id="test_run",
                )


if __name__ == "__main__":
    unittest.main()
