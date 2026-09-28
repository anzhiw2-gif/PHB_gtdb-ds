"""Failing-first tests for PhaDED catalytic-domain type 1/2 assignment.

Convention under test (Knoll M, et al. BMC Bioinformatics 2009;10:89, PMC2666664):
"Within the sequences of type 1 catalytic domains, the oxyanion hole can be found
N-terminal to the lipase box, similar to lipases. Within the sequences of type 2
catalytic domains, the oxyanion hole is found C-terminal to the catalytic triad."

That is the same relative-position convention already implemented in
``pipeline/scripts/audit_phaded_motifs.py`` (L177-180): a type-1-labelled sequence
whose oxyanion pattern starts *after* the lipase box is recorded there as
``conflict_relative_position``.  The geometric comparison is therefore
``oxyanion_start < lipase_start -> type 1`` and ``oxyanion_start > lipase_start
-> type 2``, with strict inequality and fail-closed ``undetermined_*`` otherwise.

Evidence boundary: a catalytic-domain type call is primary-sequence geometry of
candidate homology only.  It is not a validated PHB/PHA degradation phenotype.
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


PIPELINE = Path(__file__).resolve().parents[1]
SCRIPT = PIPELINE / "scripts" / "assign_phaded_catalytic_domain_type.py"
AUDIT_SCRIPT = PIPELINE / "scripts" / "audit_phaded_motifs.py"
PROJECT_ROOT = PIPELINE.parent
REAL_REFERENCE_FASTA = (
    PROJECT_ROOT
    / "runs"
    / "20260911_phaded_pfam_architecture_01"
    / "inputs"
    / "phaded_reference.faa"
)

ALLOWED_TYPES = {
    "type1_verified",
    "type2_verified",
    "undetermined_no_oxyanion",
    "undetermined_position_conflict",
}

MOTIF_HEADER = [
    "accession",
    "oxyanion_hole_state",
    "oxyanion_hole_coordinates",
    "lipase_box_coordinates",
    "ahsmg_state",
]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


class OrientationConventionTests(unittest.TestCase):
    """The single source of truth for which side of the lipase box is type 1."""

    def setUp(self):
        self.module = load_module("assign_phaded_catalytic_domain_type", SCRIPT)

    def test_default_orientation_is_knoll_2009_type1_oxyanion_n_terminal(self):
        self.assertIs(self.module.TYPE1_OXYANION_N_TERMINAL, True)
        docstring = self.module.__doc__ or ""
        self.assertIn("N-terminal to the lipase box", docstring)
        self.assertIn("Knoll", docstring)
        self.assertIn("PMC2666664", docstring)

    def test_type_vocabulary_is_closed_and_exactly_the_four_documented_values(self):
        self.assertEqual(
            set(self.module.ALLOWED_CATALYTIC_DOMAIN_TYPES), ALLOWED_TYPES
        )

    def test_relative_orientation_labels_are_documented_and_stable(self):
        self.assertEqual(
            self.module.RELATIVE_ORIENTATIONS,
            (
                "oxyanion_n_terminal_of_lipase_box",
                "oxyanion_c_terminal_of_lipase_box",
                "coordinate_start_positions_equal",
                "coordinates_unusable",
            ),
        )


class FailClosedTypeAssignmentTests(unittest.TestCase):

    def setUp(self):
        self.module = load_module("assign_phaded_catalytic_domain_type", SCRIPT)

    def assign(self, *args, **kwargs):
        return self.module.assign_catalytic_domain_type(*args, **kwargs)

    def test_type_requires_both_coordinate_fields(self):
        cases = [
            ("supported", "", ""),
            ("supported", "80-85", ""),
            ("supported", "", "121-125"),
        ]
        for state, oxyanion, lipase in cases:
            with self.subTest(state=state, oxyanion=oxyanion, lipase=lipase):
                result = self.assign(state, oxyanion, lipase)
                self.assertEqual(
                    result["catalytic_domain_type"], "undetermined_no_oxyanion"
                )
                self.assertEqual(result["assignment_basis"], "oxyanion_coordinates_missing")
                self.assertEqual(result["oxyanion_start"], "")
                self.assertEqual(result["lipase_start"], "")

    def test_not_detected_pattern_never_receives_a_type(self):
        for oxyanion, lipase in [("", ""), ("80-85", "121-125"), ("300-305", "121-125")]:
            with self.subTest(oxyanion=oxyanion, lipase=lipase):
                result = self.assign("not_detected_pattern", oxyanion, lipase)
                self.assertEqual(
                    result["catalytic_domain_type"], "undetermined_no_oxyanion"
                )
                self.assertEqual(
                    result["assignment_basis"], "oxyanion_signature_not_supported"
                )
                # No geometry is inferred for a signature that was never detected.
                self.assertEqual(result["relative_orientation"], "coordinates_unusable")

    def test_malformed_or_reversed_coordinates_fail_closed(self):
        for oxyanion, lipase in [
            ("abc", "121-125"),
            ("80", "121-125"),
            ("85-80", "121-125"),
            ("80-85", "125-121"),
            ("80-85", "not_a_range"),
            ("0-5", "121-125"),
            ("80..85", "121-125"),
        ]:
            with self.subTest(oxyanion=oxyanion, lipase=lipase):
                result = self.assign("supported", oxyanion, lipase)
                self.assertEqual(
                    result["catalytic_domain_type"], "undetermined_no_oxyanion"
                )
                self.assertEqual(
                    result["assignment_basis"], "invalid_coordinate_format"
                )

    def test_unknown_oxyanion_state_is_rejected_loudly(self):
        with self.assertRaises(ValueError):
            self.assign("unexpected_state", "80-85", "121-125")

    def test_every_basis_is_a_documented_basis(self):
        documented = set(self.module.ASSIGNMENT_BASES)
        self.assertTrue(documented)
        observations = [
            ("supported", "80-85", "121-125"),
            ("supported", "300-305", "121-125"),
            ("supported", "121-125", "121-125"),
            ("supported", "", ""),
            ("supported", "abc", "121-125"),
            ("not_detected_pattern", "", ""),
            ("conflict_relative_position", "300-305", "121-125"),
        ]
        for state, oxyanion, lipase in observations:
            with self.subTest(state=state, oxyanion=oxyanion, lipase=lipase):
                result = self.assign(state, oxyanion, lipase)
                self.assertIn(result["catalytic_domain_type"], ALLOWED_TYPES)
                self.assertIn(result["assignment_basis"], documented)


class GeometricTypeAssignmentTests(unittest.TestCase):

    def setUp(self):
        self.module = load_module("assign_phaded_catalytic_domain_type", SCRIPT)

    def assign(self, *args, **kwargs):
        return self.module.assign_catalytic_domain_type(*args, **kwargs)

    def test_type1_requires_oxyanion_strictly_n_terminal_of_lipase_box(self):
        result = self.assign("supported", "80-85", "121-125")
        self.assertEqual(result["catalytic_domain_type"], "type1_verified")
        self.assertEqual(
            result["relative_orientation"], "oxyanion_n_terminal_of_lipase_box"
        )
        self.assertEqual(result["oxyanion_start"], 80)
        self.assertEqual(result["lipase_start"], 121)
        self.assertEqual(result["assignment_basis"], "oxyanion_n_terminal_of_lipase_box")

    def test_type2_requires_oxyanion_strictly_c_terminal_of_lipase_box(self):
        result = self.assign("supported", "300-305", "121-125")
        self.assertEqual(result["catalytic_domain_type"], "type2_verified")
        self.assertEqual(
            result["relative_orientation"], "oxyanion_c_terminal_of_lipase_box"
        )
        self.assertEqual(result["assignment_basis"], "oxyanion_c_terminal_of_lipase_box")

    def test_comparison_uses_the_first_residue_of_each_coordinate_range(self):
        adjacent = self.assign("supported", "100-105", "101-105")
        self.assertEqual(adjacent["catalytic_domain_type"], "type1_verified")
        touching = self.assign("supported", "101-105", "101-105")
        self.assertEqual(touching["catalytic_domain_type"], "undetermined_position_conflict")
        self.assertEqual(touching["assignment_basis"], "coordinate_start_positions_equal")

    def test_relative_position_conflict_rows_are_never_promoted(self):
        for oxyanion, lipase in [("300-305", "121-125"), ("80-85", "121-125")]:
            with self.subTest(oxyanion=oxyanion, lipase=lipase):
                result = self.assign("conflict_relative_position", oxyanion, lipase)
                self.assertEqual(
                    result["catalytic_domain_type"], "undetermined_position_conflict"
                )
                self.assertEqual(
                    result["assignment_basis"], "oxyanion_relative_position_conflict"
                )

    def test_documented_orientation_switch_swaps_the_calls(self):
        """Guards the flag itself: flipping it inverts the geometry, so it stays explicit."""
        upstream = self.assign("supported", "80-85", "121-125")
        reversed_flag = self.assign(
            "supported", "80-85", "121-125", type1_oxyanion_n_terminal=False
        )
        self.assertEqual(upstream["catalytic_domain_type"], "type1_verified")
        self.assertEqual(reversed_flag["catalytic_domain_type"], "type2_verified")

    def test_convention_matches_audit_phaded_motifs_conflict_semantics(self):
        """Same geometry must give the same call in both modules (no second convention)."""
        audit = load_module("audit_phaded_motifs", AUDIT_SCRIPT)
        upstream_sequence = "M" + "A" * 20 + "HAGCQ" + "A" * 60 + "GLSAG" + "A" * 40
        downstream_sequence = "M" + "A" * 20 + "GLSAG" + "A" * 60 + "HAGCQ" + "A" * 40

        upstream = audit.audit_sequence(
            "UPSTREAM", upstream_sequence, "extracellular dPHASCL type 1"
        )
        self.assertEqual(upstream["oxyanion_hole_state"], "supported")
        self.assertEqual(
            self.assign(
                upstream["oxyanion_hole_state"],
                upstream["oxyanion_hole_coordinates"],
                upstream["lipase_box_coordinates"],
            )["catalytic_domain_type"],
            "type1_verified",
        )

        downstream = audit.audit_sequence(
            "DOWNSTREAM", downstream_sequence, "extracellular dPHASCL type 1"
        )
        self.assertEqual(downstream["oxyanion_hole_state"], "conflict_relative_position")
        self.assertEqual(
            self.assign(
                downstream["oxyanion_hole_state"],
                downstream["oxyanion_hole_coordinates"],
                downstream["lipase_box_coordinates"],
            )["catalytic_domain_type"],
            "undetermined_position_conflict",
        )

        type2 = audit.audit_sequence(
            "TYPE2", downstream_sequence, "extracellular dPHASCL type 2"
        )
        self.assertEqual(type2["oxyanion_hole_state"], "supported")
        self.assertEqual(
            self.assign(
                type2["oxyanion_hole_state"],
                type2["oxyanion_hole_coordinates"],
                type2["lipase_box_coordinates"],
            )["catalytic_domain_type"],
            "type2_verified",
        )


class AhsmgPositiveControlTests(unittest.TestCase):

    def setUp(self):
        self.module = load_module("assign_phaded_catalytic_domain_type", SCRIPT)

    def test_positive_control_is_reproduced_from_a_reference_fasta(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "panel.faa"
            fasta.write_text(
                ">DED_hfam_7_0001|AAK07742.1\n"
                + "M" * 20
                + "IVAHSMGVSMSL\n"
                + "M" * 20
                + "\n>DED_hfam_7_0003|ABV38134.1\n"
                + "M" * 60
                + "\n>DED_hfam_52_0001|P12625.1\n"
                + "M" * 60
                + "\n",
                encoding="ascii",
            )
            hits = self.module.reproduce_ahsmg_positive_controls(fasta)
            self.assertEqual(
                hits,
                [
                    {
                        "reference_id": "DED_hfam_7_0001",
                        "accession": "AAK07742.1",
                        # 20 M + IVAHSMGVSMSL (12) + 20 M = 52 residues; AHSMG at 23-27.
                        "sequence_length": 52,
                        "ahsmg_coordinates": "23-27",
                    }
                ],
            )

    def test_real_reference_panel_positive_controls_are_the_two_phaz7_references(self):
        if not REAL_REFERENCE_FASTA.is_file():
            self.skipTest(f"reference panel fasta not available: {REAL_REFERENCE_FASTA}")
        hits = self.module.reproduce_ahsmg_positive_controls(REAL_REFERENCE_FASTA)
        self.assertEqual(len(hits), 2)
        by_accession = {hit["accession"]: hit for hit in hits}
        self.assertEqual(set(by_accession), {"AAK07742.1", "2VTVA"})
        self.assertEqual(by_accession["AAK07742.1"]["ahsmg_coordinates"], "172-176")
        self.assertEqual(by_accession["2VTVA"]["ahsmg_coordinates"], "134-138")
        for hit in hits:
            self.assertTrue(hit["reference_id"].startswith("DED_hfam_7_"))

    def test_absence_statement_boundary_forbids_family_absence_claims(self):
        boundary = self.module.AHSMG_CONCLUSION_BOUNDARY
        self.assertIn("not detected", boundary)
        self.assertIn("not evidence", boundary)
        self.assertIn("may simply not contain", boundary)
        for forbidden in ("family does not exist", "enzyme does not exist", "该族不存在"):
            self.assertNotIn(forbidden, boundary)


class SummarySeparationTests(unittest.TestCase):

    def setUp(self):
        self.module = load_module("assign_phaded_catalytic_domain_type", SCRIPT)

    def build_rows(self):
        observations = [
            ("A1", "supported", "80-85", "121-125"),  # label type1 -> verified type1
            ("A2", "supported", "300-305", "121-125"),  # label type1 -> geometric type2
            ("A3", "conflict_relative_position", "300-305", "121-125"),
            ("A4", "not_detected_pattern", "", ""),
            ("A5", "supported", "", ""),  # fail-closed
            ("A6", "supported", "300-305", "121-125"),  # label type2 -> verified type2
            ("A7", "supported", "80-85", "121-125"),  # label type1 -> verified type1
            ("A8", "supported", "400-405", "121-125"),  # label type2 -> verified type2
        ]
        rows = []
        for accession, state, oxyanion, lipase in observations:
            rows.append(
                {
                    "accession": accession,
                    "oxyanion_hole_state": state,
                    "oxyanion_hole_coordinates": oxyanion,
                    "lipase_box_coordinates": lipase,
                    "ahsmg_state": "not_detected_pattern",
                }
            )
        return rows

    def build_calls(self):
        return {
            "A1": "dPHASCL1_like_candidate",
            "A2": "dPHASCL1_like_candidate",
            "A3": "dPHASCL1_like_candidate",
            "A4": "dPHASCL1_like_candidate",
            "A5": "dPHASCL1_like_candidate",
            "A6": "dPHASCL2_like_candidate",
            "A7": "dPHASCL1_like_candidate",
            "A8": "dPHASCL2_like_candidate",
        }

    def test_label_and_geometric_verification_are_reported_separately(self):
        rows = self.build_rows()
        calls = self.build_calls()
        assigned = self.module.build_assignments(rows, calls)
        summary = self.module.summarize_assignments(assigned)

        self.assertEqual(summary["total_candidates"], 8)
        self.assertEqual(summary["type_counts"]["type1_verified"], 2)
        self.assertEqual(summary["type_counts"]["type2_verified"], 3)
        self.assertEqual(summary["type_counts"]["undetermined_position_conflict"], 1)
        self.assertEqual(summary["type_counts"]["undetermined_no_oxyanion"], 2)
        self.assertTrue(summary["totals_self_consistent"])

        label_view = summary["label_vs_geometric_verification"]
        self.assertEqual(
            label_view["dPHASCL1_like_candidate"]["label_name"], "dPHASCL1_like_candidate"
        )
        self.assertEqual(label_view["dPHASCL1_like_candidate"]["labeled_count"], 6)
        self.assertEqual(label_view["dPHASCL1_like_candidate"]["type1_verified"], 2)
        self.assertEqual(label_view["dPHASCL1_like_candidate"]["type2_verified"], 1)
        self.assertEqual(label_view["dPHASCL1_like_candidate"]["undetermined_position_conflict"], 1)
        self.assertEqual(label_view["dPHASCL1_like_candidate"]["undetermined_no_oxyanion"], 2)
        self.assertEqual(label_view["dPHASCL2_like_candidate"]["labeled_count"], 2)
        self.assertEqual(label_view["dPHASCL2_like_candidate"]["type2_verified"], 2)
        self.assertEqual(label_view["dPHASCL2_like_candidate"]["type1_verified"], 0)
        self.assertNotEqual(
            label_view["dPHASCL1_like_candidate"]["labeled_count"],
            label_view["dPHASCL1_like_candidate"]["type1_verified"],
        )

    def test_label_is_never_used_to_assign_a_type(self):
        """The same geometry must produce the same call under either label."""
        for accession, state, oxyanion, lipase in [("A1", "supported", "80-85", "121-125")]:
            row = {
                "accession": accession,
                "oxyanion_hole_state": state,
                "oxyanion_hole_coordinates": oxyanion,
                "lipase_box_coordinates": lipase,
                "ahsmg_state": "not_detected_pattern",
            }
            for label in ("dPHASCL1_like_candidate", "dPHASCL2_like_candidate", ""):
                assigned = self.module.build_assignments([row], {accession: label})
                self.assertEqual(
                    assigned[0]["catalytic_domain_type"], "type1_verified"
                )

    def test_inverted_orientation_sensitivity_is_quantified(self):
        rows = self.build_rows()
        calls = {row["accession"]: "dPHASCL1_like_candidate" for row in rows}
        assigned = self.module.build_assignments(rows, calls)
        summary = self.module.summarize_assignments(assigned)
        inverted = summary["orientation_sensitivity_if_inverted"]
        # 3 supported rows carry the oxyanion hole C-terminal to the lipase box and
        # 2 carry it N-terminal, so flipping the single documented assumption swaps
        # the type1/type2 totals (2/3 -> 3/2) and moves nothing else.
        self.assertEqual(inverted["type1_verified"], 3)
        self.assertEqual(inverted["type2_verified"], 2)
        self.assertEqual(inverted["undetermined_position_conflict"], 1)
        self.assertEqual(inverted["undetermined_no_oxyanion"], 2)
        self.assertIn("Knoll", inverted["note"])


class LabelSourceConsistencyTests(unittest.TestCase):

    def setUp(self):
        self.module = load_module("assign_phaded_catalytic_domain_type", SCRIPT)

    def matrices(self, second_values):
        first = {"A1": "dPHASCL1_like_candidate", "A2": "dPHASCL2_like_candidate"}
        second = dict(first)
        second.update(second_values)
        return {"first.tsv": first, "second.tsv": second}

    def test_agreeing_matrices_merge_and_are_reported_as_consistent(self):
        matrices = self.matrices({})
        merged = self.module.merge_subtype_calls(matrices)
        self.assertEqual(merged["A1"], "dPHASCL1_like_candidate")
        report = self.module.describe_subtype_matrix_agreement(matrices)
        self.assertEqual(report["checked_matrices"], ["first.tsv", "second.tsv"])
        self.assertEqual(report["checked_accessions"], 2)
        self.assertTrue(report["accession_sets_identical"])
        self.assertTrue(report["subtype_call_identical"])
        self.assertEqual(report["disagreeing_accession_count"], 0)

    def test_disagreeing_matrices_are_rejected(self):
        with self.assertRaises(ValueError):
            self.module.merge_subtype_calls(
                self.matrices({"A2": "dPHASCL1_like_candidate"})
            )

    def test_accession_set_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):
            self.module.merge_subtype_calls(
                {"first.tsv": {"A1": "x"}, "second.tsv": {"A2": "x"}}
            )


class RunIntegrationTests(unittest.TestCase):

    def setUp(self):
        self.module = load_module("assign_phaded_catalytic_domain_type", SCRIPT)

    def test_cli_creates_run_layout_and_writes_closed_vocabulary_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            motif_rows = [
                {
                    "accession": "A1",
                    "oxyanion_hole_state": "supported",
                    "oxyanion_hole_coordinates": "80-85",
                    "lipase_box_coordinates": "121-125",
                    "ahsmg_state": "not_detected_pattern",
                },
                {
                    "accession": "A2",
                    "oxyanion_hole_state": "not_detected_pattern",
                    "oxyanion_hole_coordinates": "",
                    "lipase_box_coordinates": "",
                    "ahsmg_state": "not_detected_pattern",
                },
                {
                    "accession": "A3",
                    "oxyanion_hole_state": "conflict_relative_position",
                    "oxyanion_hole_coordinates": "300-305",
                    "lipase_box_coordinates": "121-125",
                    "ahsmg_state": "not_detected_pattern",
                },
                {
                    "accession": "A4",
                    "oxyanion_hole_state": "supported",
                    "oxyanion_hole_coordinates": "300-305",
                    "lipase_box_coordinates": "121-125",
                    "ahsmg_state": "not_detected_pattern",
                },
            ]
            motif_path = write_tsv(root / "inputs" / "motif.tsv", MOTIF_HEADER, motif_rows)
            matrix_path = write_tsv(
                root / "inputs" / "matrix.tsv",
                ["accession", "subtype_call"],
                [
                    {"accession": "A1", "subtype_call": "dPHASCL1_like_candidate"},
                    {"accession": "A2", "subtype_call": "dPHASCL1_like_candidate"},
                    {"accession": "A3", "subtype_call": "dPHASCL1_like_candidate"},
                    {"accession": "A4", "subtype_call": "dPHASCL2_like_candidate"},
                ],
            )
            panel_path = write_tsv(
                root / "inputs" / "panel.tsv",
                ["accession", "phaded_family_id", "ahsmg_state", "ahsmg_coordinates"],
                [
                    {
                        "accession": "AAK07742.1",
                        "phaded_family_id": "DED_hfam_7",
                        "ahsmg_state": "supported",
                        "ahsmg_coordinates": "172-176",
                    },
                    {
                        "accession": "P12625.1",
                        "phaded_family_id": "DED_hfam_52",
                        "ahsmg_state": "not_detected_pattern",
                        "ahsmg_coordinates": "",
                    },
                ],
            )
            fasta_path = root / "inputs" / "panel.faa"
            fasta_path.write_text(
                ">DED_hfam_7_0001|AAK07742.1\n" + "M" * 20 + "IVAHSMGVSMSL" + "M" * 20 + "\n",
                encoding="ascii",
            )
            run_dir = root / "runs" / "20260917_test_phaded_catalytic_domain_type_01"

            exit_code = self.module.main(
                [
                    "--motif-evidence",
                    str(motif_path),
                    "--subtype-matrix",
                    str(matrix_path),
                    "--reference-panel",
                    str(panel_path),
                    "--reference-fasta",
                    str(fasta_path),
                    "--run-dir",
                    str(run_dir),
                    "--create-run",
                    "--expected-total-candidates",
                    "4",
                    "--expected-dphascl1-labeled",
                    "3",
                ]
            )
            self.assertEqual(exit_code, 0)

            for name in ("logs", "inputs", "results"):
                self.assertTrue((run_dir / name).is_dir(), name)
            self.assertTrue((run_dir / "input_contract.json").is_file())

            table = run_dir / "results" / "catalytic_domain_type.tsv"
            with table.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 4)
            self.assertEqual(
                {row["catalytic_domain_type"] for row in rows},
                {
                    "type1_verified",
                    "type2_verified",
                    "undetermined_no_oxyanion",
                    "undetermined_position_conflict",
                },
            )
            self.assertEqual({row["accession"] for row in rows}, {"A1", "A2", "A3", "A4"})

            summary = json.loads(
                (run_dir / "results" / "type_summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(summary["total_candidates"], 4)
            self.assertEqual(summary["type_counts"], {
                "type1_verified": 1,
                "type2_verified": 1,
                "undetermined_no_oxyanion": 1,
                "undetermined_position_conflict": 1,
            })
            self.assertTrue(summary["totals_self_consistent"])
            self.assertEqual(
                summary["label_vs_geometric_verification"]["dPHASCL1_like_candidate"][
                    "labeled_count"
                ],
                3,
            )
            ahsmg = summary["ahsmg_evidence"]
            self.assertEqual(ahsmg["reference_positive_control_count"], 1)
            self.assertEqual(len(ahsmg["reference_positive_controls"]), 1)
            self.assertEqual(ahsmg["reference_positive_controls"][0]["accession"], "AAK07742.1")
            self.assertEqual(ahsmg["candidate_not_detected_pattern"], 4)
            self.assertEqual(len(ahsmg["reference_fasta_pattern_reproduction"]), 1)
            self.assertIn("not evidence", ahsmg["conclusion_boundary"])

            manifest = json.loads(
                (run_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertFalse(manifest["authorization"]["server_execution_started"])
            self.assertFalse(manifest["authorization"]["ssh_started"])
            self.assertFalse(manifest["authorization"]["historical_run_modified"])
            self.assertFalse(manifest["authorization"]["formal_scan_models_modified"])
            self.assertTrue(
                manifest["comparison_convention"][
                    "matches_audit_phaded_motifs_relative_position_semantics"
                ]
            )
            self.assertEqual(
                manifest["comparison_convention"]["type1_oxyanion_n_terminal_flag"], True
            )
            # The label source check is recorded, not merely assumed.
            consistency = manifest["label_source_consistency"]
            self.assertTrue(consistency["subtype_call_identical"])
            self.assertEqual(consistency["checked_accessions"], 4)
            self.assertEqual(consistency["disagreeing_accession_count"], 0)
            self.assertEqual(len(consistency["checked_matrices"]), 1)

            contract = json.loads(
                (run_dir / "input_contract.json").read_text(encoding="utf-8")
            )
            self.assertEqual(contract["run_id"], run_dir.name)
            # GTDB taxonomy/metadata/tree are not consumed by this analysis, so the
            # contract stays honestly "pending" while every declared input verifies.
            self.assertEqual(contract["status"], "pending")
            for name, record in contract["inputs"].items():
                self.assertEqual(record["status"], "verified", name)
                self.assertEqual(len(record["sha256"]), 64, name)
                self.assertGreater(record["size"], 0, name)
            for name, record in contract["gtdb"].items():
                self.assertEqual(record["status"], "pending", name)
                self.assertIsNone(record["sha256"], name)


if __name__ == "__main__":
    unittest.main()
