#!/usr/bin/env python3
"""Tests for ``build_phaded_v2_catalog_input.py`` (task F14).

Failing-test-first, synthetic fixtures only: the adapter is exercised over tiny
hand-written tables in a temporary directory, so no project data and no frozen
run is read or written.

The tests pin the *evidence semantics*, not the plumbing:

* a missing superfamily stays absent (``pending``) and is never defaulted to a
  family -- this is the 616-row case F1 documented;
* a SignalP class is a transport prediction and never becomes localization
  evidence;
* AHSMG is a Ser motif class (``motif_class=AHSMG`` with
  ``nucleophile_identity=ser``), a Cys-family call without a GxSxG is ``cys``,
  and a Ser family without a GxSxG is *not* silently reported as ``ser``;
* an unrecognized geometry verdict is ``unresolved``, never a verified type, and
  an SBD/lid/linker annotation never changes the catalytic type;
* ``sequence_family_hmm_validated`` can only be emitted when the profile passed
  the held-out/confounder evaluation, and a family below F13's corrected
  independent-positive bar is not calibrated;
* every universe accession appears exactly once, and a mismatch fails closed.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


def _load(name: str):
    path = SCRIPTS_DIR / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault(name, module)
    spec.loader.exec_module(module)
    return module


ADAPTER = _load("build_phaded_v2_catalog_input")
CATALOG = _load("build_phaded_candidate_catalog_v2")


def write_tsv(path: Path, columns, rows) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
    return path


UNIVERSE_COLUMNS = [
    "accession", "genome", "family", "source_layer", "signal_type",
    "assignment", "length", "layer",
]

SUBTYPE_COLUMNS = [
    "accession", "genome", "family", "layer", "signal_type", "pha_subtype_evidence",
    "phaded_superfamily_best", "phaded_superfamily_second", "superfamily_score_gap",
    "phaded_family_best", "phaded_family_second", "family_score_gap",
    "assignment_status", "evidence_status", "evidence_model", "evidence_threshold",
    "sequence_integrity", "signalp_evidence", "signalp_class",
    "lipase_box_state", "lipase_box_coordinates", "ahsmg_state", "ahsmg_coordinates",
    "catalytic_ser_cys_state", "catalytic_ser_cys_coordinates",
    "his_state", "his_coordinates", "asp_state", "asp_coordinates",
    "oxyanion_hole_state", "oxyanion_hole_coordinates",
    "sbd_state", "sbd_coordinates", "linker_state", "linker_coordinates",
    "lid_state", "lid_coordinates", "pfam_interpro_state", "architecture_consistency",
    "feature_evidence_status", "pfam_accessions", "pfam_coordinates",
    "assignment_review", "evidence_grade", "evidence_basis",
    "structure_evidence_status", "phylogeny_evidence_status",
    "interpro_evidence_status", "literature_evidence_tier",
    "motif_lipase_box", "motif_catalytic_ser_cys", "motif_his", "motif_asp",
    "motif_oxyanion_hole", "motif_sbd", "motif_linker", "motif_lid",
    "profile_model_status", "profile_evidence_status", "domain_evidence_status",
    "domain_pfam_state", "domain_interpro_state", "motif_evidence_status",
    "motif_completeness", "localization_evidence_status",
    "structure_evidence_status_v2", "phylogeny_evidence_status_v2",
    "manual_review_status", "subtype_call", "subtype_confidence",
    "subtype_evidence_boundary", "interpro_status", "interpro_secondary_state",
    "interpro_hit_count", "interpro_signatures", "interpro_ids",
]

MERGE_COLUMNS = [
    "accession", "pool_origin", "genome", "superfamily", "signalp_class",
    "profile_best_evalue", "profile_evidence_status", "lipase_box_state",
    "lipase_box_x1", "catalytic_domain_type", "catalytic_residue_verification",
    "nucleophile_type", "nucleophile_family", "lid_state",
    "sbd_pf06850_binding_state", "ahsmg_state", "interpro_status",
    "high_confidence",
]

HOLD_COLUMNS = [
    "accession", "genome", "superfamily", "nucleophile_type", "hold_reason",
    "lipase_box_state", "sbd_pf06850_binding_state", "signalp_class",
    "catalytic_domain_type", "interpro_status", "profile_evidence_status",
]

DEMOTION_COLUMNS = [
    "accession", "genome", "superfamily", "signalp_class", "profile_best_evalue",
    "profile_evidence_status", "lipase_box_state", "lipase_box_x1",
    "catalytic_domain_type", "catalytic_residue_verification", "lid_state",
    "sbd_pf06850_binding_state", "ahsmg_state", "interpro_status",
    "high_confidence", "hold_reason", "signalp_class_prediction", "demote_basis",
]

AUTHORITY_COLUMNS = [
    "layer_id", "layer_kind", "source_evidence_type", "source_version",
    "registry_eligible", "gate_profile_id", "notes",
]

PROFILE_COLUMNS = [
    "profile_id", "profile_kind", "phaded_superfamily", "phaded_family_id",
    "training_accessions", "training_count", "model_status", "model_reason",
    "alignment_sha256", "hmm_sha256", "mafft_version", "hmmer_version",
    "alignment_path", "hmm_path",
]

LEDGER_COLUMNS = ["reference_id", "accession", "phaded_family_id", "experimental_evidence_grade"]


class Fixture:
    """A tiny synthetic project laid out the way the real inputs are."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self.universe = {
            "accession": [
                "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8",
            ],
        }
        self.universe_rows = [
            {"accession": "A1", "genome": "G1", "family": "ePhaZ", "layer": "ePhaZ_curated_core_secreted"},
            {"accession": "A2", "genome": "G1", "family": "ePhaZ", "layer": "ePhaZ_curated_core_nonsecreted"},
            {"accession": "A3", "genome": "G2", "family": "iPhaZ", "layer": "iPhaZ_tier1"},
            {"accession": "A4", "genome": "G2", "family": "iPhaZ", "layer": "iPhaZ_tier1"},
            {"accession": "A5", "genome": "G3", "family": "iPhaZ", "layer": "iPhaZ_tier2_review"},
            {"accession": "A6", "genome": "G3", "family": "ePhaZ", "layer": "ePhaZ_tier2_review"},
            {"accession": "A7", "genome": "G4", "family": "ePhaZ", "layer": "ePhaZ_tier2_review"},
            {"accession": "A8", "genome": "G4", "family": "ePhaZ", "layer": "ePhaZ_curated_core_secreted"},
        ]
        # A1: trained profile hit on the calibrated family DED_hfam_52, GxSxG, SP
        # A2: calibrated family DED_hfam_70, AHSMG motif
        # A3: family DED_hfam_4 (below the corrected >=3 bar), no GxSxG, Cys family
        # A4: family DED_hfam_55 (below the bar), GxSxG
        # A5: DED_hfam_8 (below the bar), AHSMG
        # A6: no profile hit at all (the 472-row class)
        # A7: ambiguous superfamily, no profile call (no superfamily at all)
        # A8: profile call with an unrecognized geometry verdict + SBD/lid records
        self.subtype_rows = [
            self._subtype("A1", sf="extracellular dPHASCL type 1", fam="profile_hfam_52",
                          signalp="SP", lipase="supported", ser_cys="supported",
                          ahsmg="not_detected_pattern", oxy="supported",
                          profile_status="trained", assignment="assigned",
                          architecture="partial", integrity="valid"),
            self._subtype("A2", sf="extracellular native-SCL/PhaZ7-like", fam="profile_hfam_70",
                          signalp="TATLIPO", lipase="supported", ser_cys="supported",
                          ahsmg="supported", oxy="supported",
                          profile_status="trained", assignment="assigned",
                          architecture="partial", integrity="valid"),
            self._subtype("A3", sf="intracellular nPHASCL without lipase box", fam="profile_hfam_4",
                          signalp="OTHER", lipase="not_detected_pattern", ser_cys="pending",
                          ahsmg="not_detected_pattern", oxy="not_detected_pattern",
                          profile_status="trained", assignment="assigned",
                          architecture="partial", integrity="valid"),
            self._subtype("A4", sf="intracellular nPHAMCL", fam="profile_hfam_55",
                          signalp="pending", lipase="supported", ser_cys="supported",
                          ahsmg="not_detected_pattern", oxy="supported",
                          profile_status="trained", assignment="assigned",
                          architecture="partial", integrity="valid"),
            self._subtype("A5", sf="extracellular dPHASCL type 1", fam="profile_hfam_8",
                          signalp="LIPO", lipase="supported", ser_cys="supported",
                          ahsmg="supported", oxy="supported",
                          profile_status="trained", assignment="assigned",
                          architecture="partial", integrity="valid"),
            self._subtype("A6", sf="", fam="", signalp="pending", lipase="not_detected_pattern",
                          ser_cys="pending", ahsmg="not_detected_pattern",
                          oxy="not_detected_pattern", profile_status="unassigned",
                          assignment="unassigned_PhaDED_like", architecture="pending",
                          integrity="terminal_stop_only"),
            self._subtype("A7", sf="", fam="", signalp="OTHER", lipase="supported",
                          ser_cys="supported", ahsmg="not_detected_pattern", oxy="supported",
                          profile_status="ambiguous", assignment="ambiguous_superfamily",
                          architecture="partial", integrity="valid"),
            self._subtype("A8", sf="intracellular nPHASCL with lipase box", fam="profile_hfam_52",
                          signalp="SP", lipase="supported", ser_cys="supported",
                          ahsmg="not_detected_pattern", oxy="not_detected_pattern",
                          profile_status="trained", assignment="assigned",
                          architecture="conflicting", integrity="valid"),
        ]
        self.merge_rows = [
            {"accession": "A1", "pool_origin": "pool_in", "genome": "G1",
             "superfamily": "extracellular dPHASCL type 1", "signalp_class": "SP",
             "profile_best_evalue": "1e-200", "profile_evidence_status": "profile_trained_hit",
             "lipase_box_state": "supported", "lipase_box_x1": "A",
             "catalytic_domain_type": "type1_verified",
             "catalytic_residue_verification": "motif_level", "nucleophile_type": "ser",
             "nucleophile_family": "GxSxG", "lid_state": "not_detected_pattern",
             "sbd_pf06850_binding_state": "not_detected_in_tested_pfam",
             "ahsmg_state": "not_detected_pattern", "interpro_status": "interpro_supported",
             "high_confidence": "1"},
            {"accession": "A2", "pool_origin": "pool_in", "genome": "G1",
             "superfamily": "extracellular native-SCL/PhaZ7-like", "signalp_class": "TATLIPO",
             "profile_best_evalue": "1e-180", "profile_evidence_status": "profile_trained_hit",
             "lipase_box_state": "supported", "lipase_box_x1": "A",
             "catalytic_domain_type": "type1_verified",
             "catalytic_residue_verification": "motif_level", "nucleophile_type": "ahsmg",
             "nucleophile_family": "GxSxG", "lid_state": "not_detected_pattern",
             "sbd_pf06850_binding_state": "not_detected_in_tested_pfam",
             "ahsmg_state": "supported", "interpro_status": "interpro_supported",
             "high_confidence": "1"},
        ]
        self.hold_rows = [
            {"accession": "A3", "genome": "G2",
             "superfamily": "intracellular nPHASCL without lipase box",
             "nucleophile_type": "cys", "hold_reason": "localization_conflict",
             "lipase_box_state": "not_detected_pattern",
             "sbd_pf06850_binding_state": "detected", "signalp_class": "OTHER",
             "catalytic_domain_type": "undetermined_no_oxyanion",
             "interpro_status": "interpro_supported",
             "profile_evidence_status": "profile_trained_hit"},
            {"accession": "A4", "genome": "G2", "superfamily": "intracellular nPHAMCL",
             "nucleophile_type": "ser", "hold_reason": "architecture_criteria",
             "lipase_box_state": "supported",
             "sbd_pf06850_binding_state": "not_detected_in_tested_pfam",
             "signalp_class": "pending", "catalytic_domain_type": "type2_verified",
             "interpro_status": "interpro_supported",
             "profile_evidence_status": "profile_trained_hit"},
        ]
        self.demotion_rows = [
            {"accession": "A3", "genome": "G2",
             "superfamily": "intracellular nPHASCL without lipase box",
             "signalp_class": "pending", "profile_best_evalue": "1e-100",
             "profile_evidence_status": "profile_trained_hit",
             "lipase_box_state": "not_detected_pattern", "lipase_box_x1": "",
             "catalytic_domain_type": "undetermined_no_oxyanion",
             "catalytic_residue_verification": "inferred_cys_not_verified",
             "lid_state": "not_detected_pattern",
             "sbd_pf06850_binding_state": "detected", "ahsmg_state": "not_detected_pattern",
             "interpro_status": "interpro_supported", "high_confidence": "1",
             "hold_reason": "localization_conflict_signalp_predicted_export",
             "signalp_class_prediction": "SP",
             "demote_basis": "SignalP 6.0h prediction (not experimental)"},
        ]
        self.authority_rows = [
            {"layer_id": "extracellular dPHASCL type 1", "layer_kind": "superfamily",
             "registry_eligible": "true", "gate_profile_id": "profile_sf1"},
            {"layer_id": "superfamily_calibrated", "layer_kind": "superfamily",
             "registry_eligible": "true", "gate_profile_id": "profile_sf2"},
            {"layer_id": "intracellular nPHASCL without lipase box", "layer_kind": "superfamily",
             "registry_eligible": "true", "gate_profile_id": "profile_sf3"},
            {"layer_id": "DED_hfam_52", "layer_kind": "family", "registry_eligible": "false",
             "gate_profile_id": "profile_hfam_52"},
            {"layer_id": "DED_hfam_70", "layer_kind": "family", "registry_eligible": "false",
             "gate_profile_id": "profile_hfam_70"},
            {"layer_id": "DED_hfam_4", "layer_kind": "family", "registry_eligible": "false",
             "gate_profile_id": "profile_hfam_4"},
            {"layer_id": "DED_hfam_55", "layer_kind": "family", "registry_eligible": "false",
             "gate_profile_id": "profile_hfam_55"},
            {"layer_id": "DED_hfam_8", "layer_kind": "family", "registry_eligible": "false",
             "gate_profile_id": "profile_hfam_8"},
        ]
        self.profile_rows = [
            {"profile_id": "profile_sf1", "profile_kind": "superfamily",
             "phaded_superfamily": "extracellular dPHASCL type 1", "phaded_family_id": "",
             "model_status": "trained", "hmm_sha256": "aa11"},
            {"profile_id": "profile_sf2", "profile_kind": "superfamily",
             "phaded_superfamily": "superfamily_calibrated", "phaded_family_id": "",
             "model_status": "trained", "hmm_sha256": "bb22"},
            {"profile_id": "profile_sf3", "profile_kind": "superfamily",
             "phaded_superfamily": "intracellular nPHASCL without lipase box",
             "phaded_family_id": "", "model_status": "reference_only", "hmm_sha256": ""},
            {"profile_id": "profile_hfam_52", "profile_kind": "family",
             "phaded_superfamily": "extracellular dPHASCL type 1",
             "phaded_family_id": "DED_hfam_52", "model_status": "trained",
             "hmm_sha256": "cc33"},
            {"profile_id": "profile_hfam_70", "profile_kind": "family",
             "phaded_superfamily": "extracellular native-SCL/PhaZ7-like",
             "phaded_family_id": "DED_hfam_70", "model_status": "trained",
             "hmm_sha256": "dd44"},
            {"profile_id": "profile_hfam_4", "profile_kind": "family",
             "phaded_superfamily": "intracellular nPHAMCL", "phaded_family_id": "DED_hfam_4",
             "model_status": "trained", "hmm_sha256": "ee55"},
            {"profile_id": "profile_hfam_55", "profile_kind": "family",
             "phaded_superfamily": "extracellular dPHASCL type 1",
             "phaded_family_id": "DED_hfam_55", "model_status": "trained",
             "hmm_sha256": "ff66"},
            {"profile_id": "profile_hfam_8", "profile_kind": "family",
             "phaded_superfamily": "extracellular dPHAMCL", "phaded_family_id": "DED_hfam_8",
             "model_status": "trained", "hmm_sha256": "gg77"},
        ]
        # F13's corrected bar: only DED_hfam_52 and DED_hfam_70 reach >= 3 unique
        # independence groups.  DED_hfam_4 = 1, DED_hfam_55 = 2, DED_hfam_8 = 2.
        self.family_positive_counts = {
            "DED_hfam_52": 6, "DED_hfam_70": 7,
            "DED_hfam_4": 1, "DED_hfam_55": 2, "DED_hfam_8": 2,
        }
        self.independence_summary = {
            "minimum_independent_positive_count": 3,
            "positive_count_basis": "unique_independence_group",
            "family_functional_positive_counts": {
                "DED_hfam_52": 6, "DED_hfam_70": 7, "DED_hfam_4": 1,
                "DED_hfam_55": 2, "DED_hfam_8": 2,
            },
            "independence_summary": {
                "family_functional_positive_counts": dict(self.family_positive_counts),
                "minimum_independent_positive_count": 3,
            },
        }
        # The curated ledger is disjoint from the universe, exactly as in the
        # frozen evidence (0 of 723 reference accessions are GTDB candidates).
        self.ledger_rows = [
            {"reference_id": "R1", "accession": "REF1", "phaded_family_id": "DED_hfam_52",
             "experimental_evidence_grade": "E3"},
        ]

    @staticmethod
    def _subtype(accession, *, sf, fam, signalp, lipase, ser_cys, ahsmg, oxy,
                 profile_status, assignment, architecture, integrity):
        return {
            "accession": accession, "genome": "G1", "family": "ePhaZ", "layer": "x",
            "signal_type": signalp if signalp != "pending" else "", "pha_subtype_evidence": "x",
            "phaded_superfamily_best": sf, "phaded_superfamily_second": "",
            "superfamily_score_gap": "2.0" if sf else "0.0",
            "phaded_family_best": fam, "phaded_family_second": "", "family_score_gap": "2.0",
            "assignment_status": assignment, "evidence_status": "profile_score",
            "evidence_model": "PhaDED_profile_score",
            "evidence_threshold": "min_score_gap_1.0_bits_HMMER_E1e6",
            "sequence_integrity": integrity, "signalp_evidence": "accepted_layer_signalp",
            "signalp_class": signalp, "lipase_box_state": lipase, "lipase_box_coordinates": "10-14",
            "ahsmg_state": ahsmg, "ahsmg_coordinates": "",
            "catalytic_ser_cys_state": ser_cys, "catalytic_ser_cys_coordinates": "12",
            "his_state": "pending", "his_coordinates": "", "asp_state": "pending",
            "asp_coordinates": "", "oxyanion_hole_state": oxy, "oxyanion_hole_coordinates": "16",
            "sbd_state": "pending", "sbd_coordinates": "", "linker_state": "pending",
            "linker_coordinates": "", "lid_state": "pending", "lid_coordinates": "",
            "pfam_interpro_state": "supported", "architecture_consistency": architecture,
            "feature_evidence_status": "x", "pfam_accessions": "", "pfam_coordinates": "",
            "assignment_review": "no_conflict_recorded", "evidence_grade": "L0_profile_only",
            "evidence_basis": "x", "structure_evidence_status": "x",
            "phylogeny_evidence_status": "x", "interpro_evidence_status": "x",
            "literature_evidence_tier": "x", "motif_lipase_box": "", "motif_catalytic_ser_cys": "",
            "motif_his": "", "motif_asp": "", "motif_oxyanion_hole": "", "motif_sbd": "",
            "motif_linker": "", "motif_lid": "", "profile_model_status": profile_status,
            "profile_evidence_status": "profile_trained_hit", "domain_evidence_status": "x",
            "domain_pfam_state": "", "domain_interpro_state": "", "motif_evidence_status": "",
            "motif_completeness": "pattern_partial",
            "localization_evidence_status": "localization_unassigned",
            "structure_evidence_status_v2": "not_tested_full_library",
            "phylogeny_evidence_status_v2": "not_tested_full_library",
            "manual_review_status": "x", "subtype_call": "x", "subtype_confidence": "x",
            "subtype_evidence_boundary": "x", "interpro_status": "interpro_supported",
            "interpro_secondary_state": "", "interpro_hit_count": "1",
            "interpro_signatures": "", "interpro_ids": "",
        }

    def write(self):
        write_tsv(self.root / "protein_layers.tsv", UNIVERSE_COLUMNS, self.universe_rows)
        write_tsv(self.root / "subtype_matrix.tsv", SUBTYPE_COLUMNS, self.subtype_rows)
        write_tsv(self.root / "merged.tsv", MERGE_COLUMNS, self.merge_rows)
        write_tsv(self.root / "hold.tsv", HOLD_COLUMNS, self.hold_rows)
        write_tsv(self.root / "demotion.tsv", DEMOTION_COLUMNS, self.demotion_rows)
        write_tsv(self.root / "authority.tsv", AUTHORITY_COLUMNS, self.authority_rows)
        write_tsv(self.root / "profile_manifest.tsv", PROFILE_COLUMNS, self.profile_rows)
        write_tsv(self.root / "ledger.tsv", LEDGER_COLUMNS, self.ledger_rows)
        (self.root / "independence.json").write_text(
            json.dumps(self.independence_summary, indent=2), encoding="utf-8"
        )
        return self

    def build(self, **overrides):
        kwargs = dict(
            universe_path=self.root / "protein_layers.tsv",
            subtype_matrix_path=self.root / "subtype_matrix.tsv",
            v1_merge_path=self.root / "merged.tsv",
            v1_hold_path=self.root / "hold.tsv",
            demotion_path=self.root / "demotion.tsv",
            authority_path=self.root / "authority.tsv",
            profile_manifest_path=self.root / "profile_manifest.tsv",
            ledger_path=self.root / "ledger.tsv",
            independence_summary_path=self.root / "independence.json",
        )
        kwargs.update(overrides)
        return ADAPTER.build_adapted_input(**kwargs)


class ScoreGapColumnTests(unittest.TestCase):
    """The catalogue must carry the value its own threshold is compared against.

    It records ``evidence_threshold = min_score_gap_1.0_bits_HMMER_E1e6`` and an
    ``assignment_status`` of ``assigned`` versus ``ambiguous_*``, but until this
    change it did not carry the score gap those two are defined on. Measured on the
    1,087 rows whose best family is the promoted DED_hfam_70 profile, the gap
    separates them perfectly at the documented 1.0-bit threshold: every ``assigned``
    row is at 1.0-3.3 bits and every ``ambiguous_family`` row at 0.0-0.9, with no
    overlap. A status without its number is an assertion without its evidence.
    """

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.fixture = Fixture(self.root).write()
        self.result = self.fixture.build()
        self.rows = {row["accession"]: row for row in self.result["rows"]}

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def test_both_gap_columns_are_carried(self):
        for column in ("family_score_gap", "superfamily_score_gap"):
            self.assertIn(column, self.result["rows"][0], f"{column} missing from the output")

    def test_gap_values_come_from_the_subtype_matrix(self):
        # The fixture writes family_score_gap 2.0.
        self.assertEqual(self.rows["A1"]["family_score_gap"], "2.0")

    def test_the_threshold_is_carried_alongside_the_gap(self):
        # Carrying the gap is only meaningful next to the threshold it is compared to.
        self.assertEqual(self.rows["A1"]["evidence_threshold"],
                         "min_score_gap_1.0_bits_HMMER_E1e6")

    def test_derive_score_gaps_reads_the_subtype_row_and_never_invents_a_zero(self):
        self.assertEqual(ADAPTER.derive_score_gaps({"family_score_gap": "1.3",
                                            "superfamily_score_gap": "0.4"}),
                         ("1.3", "0.4"))
        # A subtype row without the columns reports absence rather than a zero.
        self.assertEqual(ADAPTER.derive_score_gaps({}), ("", ""))


class AdapterFixtureTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.fixture = Fixture(self.root).write()
        self.result = self.fixture.build()
        self.rows = {row["accession"]: row for row in self.result["rows"]}

    def tearDown(self) -> None:
        self._tmp.cleanup()

    # -- row identity ------------------------------------------------------

    def test_one_row_per_universe_accession_in_universe_order(self):
        self.assertEqual(len(self.result["rows"]), 8)
        self.assertEqual(
            [row["accession"] for row in self.result["rows"]],
            ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8"],
        )
        self.assertEqual(self.result["report"]["counts"]["universe_rows"], 8)
        self.assertEqual(self.result["report"]["counts"]["joined_rows"], 8)

    def test_duplicate_accession_fails_closed(self):
        rows = list(self.fixture.universe_rows) + [dict(self.fixture.universe_rows[0])]
        write_tsv(self.root / "universe_dup.tsv", UNIVERSE_COLUMNS, rows)
        with self.assertRaises(ADAPTER.InputError) as ctx:
            self.fixture.build(universe_path=self.root / "universe_dup.tsv")
        self.assertEqual(ctx.exception.rule, ADAPTER.RULE_UNIVERSE_DUPLICATE)

    def test_universe_row_without_subtype_matrix_row_fails_closed(self):
        subtype = [row for row in self.fixture.subtype_rows if row["accession"] != "A5"]
        write_tsv(self.root / "subtype_short.tsv", SUBTYPE_COLUMNS, subtype)
        with self.assertRaises(ADAPTER.InputError) as ctx:
            self.fixture.build(subtype_matrix_path=self.root / "subtype_short.tsv")
        self.assertEqual(ctx.exception.rule, ADAPTER.RULE_UNIVERSE_NOT_COVERED)
        self.assertIn("A5", ctx.exception.message)

    # -- historical superfamily -------------------------------------------

    def test_absent_superfamily_stays_pending_never_a_family(self):
        row = self.rows["A6"]
        self.assertEqual(row["historical_superfamily"], "")
        classified = CATALOG.classify_row(row)
        self.assertEqual(classified["historical_superfamily"], "pending")
        self.assertEqual(classified["primary_disposition"], "remote_homolog_candidate")
        self.assertNotIn("excluded_input_quality", classified["primary_disposition"])

    def test_ambiguous_superfamily_stays_pending(self):
        self.assertEqual(self.rows["A7"]["historical_superfamily"], "")

    def test_v1_superfamily_used_only_as_a_fallback(self):
        # A1's historical superfamily exists in the subtype matrix and in the v1
        # merge and agrees; the fallback only supplies a missing value.
        self.assertEqual(self.rows["A1"]["historical_superfamily"],
                         "extracellular dPHASCL type 1")
        report = self.result["report"]
        self.assertEqual(report["superfamily_source_counts"]["subtype_matrix"], 6)
        self.assertEqual(report["superfamily_source_counts"]["v1_fallback"], 0)
        self.assertEqual(report["superfamily_source_counts"]["absent"], 2)

    def test_v1_fallback_supplies_a_missing_superfamily(self):
        subtype = [
            dict(row, phaded_superfamily_best="" if row["accession"] == "A1" else row["phaded_superfamily_best"])
            for row in self.fixture.subtype_rows
        ]
        write_tsv(self.root / "subtype_nosf.tsv", SUBTYPE_COLUMNS, subtype)
        result = self.fixture.build(subtype_matrix_path=self.root / "subtype_nosf.tsv")
        rows = {row["accession"]: row for row in result["rows"]}
        self.assertEqual(rows["A1"]["historical_superfamily"], "extracellular dPHASCL type 1")
        self.assertEqual(
            result["report"]["superfamily_source_counts"]["v1_fallback"], 1
        )

    def test_v1_fallback_disagreement_is_reported_not_silent(self):
        subtype = [
            dict(row, phaded_superfamily_best="intracellular nPHAMCL"
                 if row["accession"] == "A1" else row["phaded_superfamily_best"])
            for row in self.fixture.subtype_rows
        ]
        write_tsv(self.root / "subtype_conflict.tsv", SUBTYPE_COLUMNS, subtype)
        result = self.fixture.build(subtype_matrix_path=self.root / "subtype_conflict.tsv")
        conflicts = result["report"]["superfamily_conflicts"]
        self.assertEqual(conflicts["count"], 1)
        self.assertEqual(conflicts["accessions"], ["A1"])
        rows = {row["accession"]: row for row in result["rows"]}
        # the subtype-matrix claim (the v2 authority) wins
        self.assertEqual(rows["A1"]["historical_superfamily"], "intracellular nPHAMCL")

    # -- transport signal --------------------------------------------------

    def test_signalp_class_maps_to_transport_signal_prediction(self):
        self.assertEqual(self.rows["A1"]["signalp_class"], "SP")
        self.assertEqual(self.rows["A2"]["signalp_class"], "TATLIPO")

    def test_untested_signalp_becomes_not_tested_never_a_prediction(self):
        # Both transport columns name the same state: the v2 column and its
        # v1-column-name alias both carry the ontology's not-tested token for an
        # untested row, so the catalog builder and the comparison module agree.
        self.assertEqual(self.rows["A6"]["transport_signal_prediction"], ADAPTER.NOT_TESTED)
        self.assertEqual(self.rows["A6"]["signalp_class"], ADAPTER.NOT_TESTED)
        classified = CATALOG.classify_row(self.rows["A6"])
        self.assertEqual(classified["transport_signal_prediction"], ADAPTER.NOT_TESTED)
        self.assertNotIn(classified["transport_signal_prediction"], CATALOG.EXPORT_PREDICTIONS)
        # a non-pending class is carried through unchanged on both columns
        self.assertEqual(self.rows["A1"]["transport_signal_prediction"], "SP")
        self.assertEqual(self.rows["A1"]["signalp_class"], "SP")

    def test_a_genuinely_blank_signalp_class_uses_the_documented_default(self):
        subtype = [
            dict(row, signalp_class="", signal_type="")
            if row["accession"] == "A6" else row
            for row in self.fixture.subtype_rows
        ]
        write_tsv(self.root / "subtype_blanksp.tsv", SUBTYPE_COLUMNS, subtype)
        result = self.fixture.build(subtype_matrix_path=self.root / "subtype_blanksp.tsv")
        rows = {row["accession"]: row for row in result["rows"]}
        self.assertEqual(rows["A6"]["transport_signal_prediction"], ADAPTER.NOT_TESTED)
        self.assertEqual(rows["A6"]["signalp_class"], ADAPTER.NOT_TESTED)

    def test_transport_prediction_never_becomes_localization_evidence(self):
        classified = CATALOG.classify_row(self.rows["A2"])
        self.assertEqual(classified["transport_signal_prediction"], "TATLIPO")
        self.assertNotEqual(classified["localization_evidence"], "experimental_confirmed")
        self.assertEqual(classified["localization_evidence_source"], "")

    def test_inherited_family_label_is_historical_family_label_only(self):
        from phaded_evidence_schema import LOCALIZATION_EVIDENCE
        for row in self.result["rows"]:
            self.assertIn(row["localization_evidence"], LOCALIZATION_EVIDENCE)
            self.assertNotEqual(row["localization_evidence"], "experimental_confirmed")

    # -- nucleophile identity and motif class ------------------------------

    def test_ahsmg_is_a_ser_motif_class(self):
        classified = CATALOG.classify_row(self.rows["A2"])
        self.assertEqual(classified["motif_class"], "AHSMG")
        self.assertEqual(classified["nucleophile_identity"], "ser")

    def test_ahsmg_row_never_carries_a_third_nucleophile_identity(self):
        from phaded_evidence_schema import NUCLEOPHILE_IDENTITIES
        for row in self.result["rows"]:
            classified = CATALOG.classify_row(row)
            self.assertIn(classified["nucleophile_identity"], NUCLEOPHILE_IDENTITIES)
            self.assertNotEqual(classified["nucleophile_identity"], "ahsmg")

    def test_missing_gxsxg_with_a_cys_family_call_is_cys(self):
        classified = CATALOG.classify_row(self.rows["A3"])
        self.assertEqual(classified["nucleophile_identity"], "cys")
        self.assertEqual(classified["motif_class"], "Cys-associated")

    def test_mcl_family_keeps_an_unresolved_motif_without_a_canonical_box(self):
        """An MCL family is not assumed to carry a canonical Ser lipase box."""
        from phaded_evidence_schema import normalize_nucleophile
        out = normalize_nucleophile("intracellular nPHAMCL", "", "")
        self.assertEqual(out["motif_class"], "unresolved")
        # A4's frozen state does support the box, so the box wins over the prior
        self.assertEqual(self.rows["A4"]["motif_class"], "GxSxG")
        # and an SCL family without the box is the Cys-associated case (A3)
        scl = normalize_nucleophile(
            "intracellular nPHASCL without lipase box", "Cys-associated", ""
        )
        self.assertEqual(scl, {"nucleophile_identity": "cys", "motif_class": "Cys-associated"})

    def test_catalytic_ser_cys_state_pending_is_not_promoted(self):
        # A3's serialized ligand state is pending; the Cys family call plus the
        # missing GxSxG is what decides the identity.
        self.assertEqual(self.rows["A3"]["catalytic_ser_cys_state"], "pending")

    def test_lipase_box_supported_becomes_a_gxsxg_motif(self):
        classified = CATALOG.classify_row(self.rows["A1"])
        self.assertEqual(classified["motif_class"], "GxSxG")
        self.assertEqual(classified["nucleophile_identity"], "ser")

    # -- catalytic domain type --------------------------------------------

    def test_unrecognized_geometry_verdict_becomes_unresolved(self):
        classified = CATALOG.classify_row(self.rows["A3"])
        self.assertEqual(self.rows["A3"]["catalytic_domain_type"], "unresolved")
        self.assertEqual(classified["catalytic_domain_type"], "unresolved")
        self.assertEqual(
            self.result["report"]["catalytic_domain_type_raw_counts"]["undetermined_no_oxyanion"], 1
        )

    def test_verified_domain_type_passes_through(self):
        self.assertEqual(self.rows["A1"]["catalytic_domain_type"], "type1_verified")
        self.assertEqual(
            CATALOG.classify_row(self.rows["A1"])["catalytic_domain_type"], "type1_verified"
        )

    def test_unrecognized_accessory_state_does_not_break_the_catalog(self):
        """Regression: a frozen PF06850 ``not_tested_*`` spelling must not raise.

        The catalog builder resolves a support flag from the SBD column and
        rejects a token outside the ontology's own sets, so the adapter only
        passes a normalized token through and the verbatim value stays in the
        dedicated comparison column.
        """
        hold = [
            dict(row, sbd_pf06850_binding_state="not_tested_no_pfam_evidence",
                 lid_state="partial_lid_loop1_only", linker_state="pending_reference_annotation")
            if row["accession"] == "A3" else row
            for row in self.fixture.hold_rows
        ]
        write_tsv(self.root / "hold_odd.tsv", HOLD_COLUMNS, hold)
        result = self.fixture.build(v1_hold_path=self.root / "hold_odd.tsv")
        rows = {row["accession"]: row for row in result["rows"]}
        # the verbatim frozen value survives for the comparison
        self.assertEqual(
            rows["A3"]["sbd_pf06850_binding_state_raw"], "not_tested_no_pfam_evidence"
        )
        # and the flag column carries an ontology token
        self.assertEqual(
            rows["A3"]["sbd_pf06850_binding_state"], "not_detected_in_tested_pfam"
        )
        # classifying never raises
        classified = {row["accession"]: row for row in CATALOG.classify_rows(result["rows"])}
        self.assertEqual(classified["A3"]["primary_disposition"],
                         "remote_homolog_candidate")

    def test_detected_binding_state_is_recorded_as_supported(self):
        hold = [
            dict(row, sbd_pf06850_binding_state="detected")
            if row["accession"] == "A4" else row
            for row in self.fixture.hold_rows
        ]
        write_tsv(self.root / "hold_detected.tsv", HOLD_COLUMNS, hold)
        result = self.fixture.build(v1_hold_path=self.root / "hold_detected.tsv")
        rows = {row["accession"]: row for row in result["rows"]}
        self.assertEqual(rows["A4"]["sbd_pf06850_binding_state_raw"], "detected")
        self.assertEqual(rows["A4"]["sbd_pf06850_binding_state"], "supported")
        self.assertEqual(rows["A4"]["sbd_state"], "supported")
        classified = CATALOG.classify_row(rows["A4"])
        self.assertIn(CATALOG.FLAG_SBD_ACCESSORY, classified["evidence_flags"])
        self.assertEqual(classified["catalytic_domain_type"], "type2_verified")

    def test_sbd_state_never_influences_catalytic_domain_type(self):
        """The subtype matrix's own SBD state is independent of the type call."""
        with_sbd = [
            dict(row, sbd_state="supported" if row["accession"] == "A7" else row["sbd_state"])
            for row in self.fixture.subtype_rows
        ]
        write_tsv(self.root / "subtype_sbd.tsv", SUBTYPE_COLUMNS, with_sbd)
        result = self.fixture.build(subtype_matrix_path=self.root / "subtype_sbd.tsv")
        rows = {row["accession"]: row for row in result["rows"]}
        # A7 has no v1 row, so the subtype matrix SBD state is what the adapter reads
        self.assertEqual(rows["A7"]["sbd_pf06850_binding_state_raw"], "supported")
        self.assertEqual(rows["A7"]["sbd_state"], "supported")
        classified = CATALOG.classify_row(rows["A7"])
        self.assertIn(CATALOG.FLAG_SBD_ACCESSORY, classified["evidence_flags"])
        # the accessory annotation is recorded and the type call stays untouched
        self.assertIn("SBD:supported", classified["accessory_domain_architecture"])
        self.assertEqual(classified["catalytic_domain_type"],
                         CATALOG.classify_row(self.rows["A7"])["catalytic_domain_type"])

    # -- model layer -------------------------------------------------------

    def test_calibrated_family_reaches_the_discovery_layer(self):
        self.assertEqual(self.rows["A1"]["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(self.rows["A2"]["model_layer"], "discovery_hmm_uncalibrated")

    def test_family_below_the_bar_does_not_reach_the_discovery_layer(self):
        for accession in ("A3", "A4", "A5"):
            self.assertEqual(
                self.rows[accession]["model_layer"], "reference_query_only", accession
            )

    def test_unassigned_profile_is_reference_query_only(self):
        self.assertEqual(self.rows["A6"]["model_layer"], "reference_query_only")
        self.assertEqual(self.rows["A7"]["model_layer"], "reference_query_only")

    def test_sequence_family_hmm_validated_is_never_emitted_without_the_evaluation(self):
        for row in self.result["rows"]:
            self.assertNotEqual(row["model_layer"], "sequence_family_hmm_validated")
        report = self.result["report"]
        self.assertFalse(report["sequence_model_evaluation"]["inputs_available"])
        self.assertEqual(report["model_layer_distribution"], {
            "discovery_hmm_uncalibrated": 3, "reference_query_only": 5,
        })

    def test_corrected_bar_is_recorded_with_the_below_bar_families(self):
        bar = self.result["report"]["calibration_bar"]
        self.assertEqual(bar["minimum_independent_positive_count"], 3)
        self.assertEqual(bar["positive_count_basis"], "unique_independence_group")
        self.assertEqual(bar["families_at_or_above_bar"], {"DED_hfam_52": 6, "DED_hfam_70": 7})
        self.assertEqual(
            bar["families_below_bar"],
            {"DED_hfam_4": 1, "DED_hfam_55": 2, "DED_hfam_8": 2},
        )
        self.assertFalse(bar["frozen_v1_trained_label_kept_for_below_bar_families"])

    def test_functional_calibration_status_is_independent_and_not_inferred(self):
        for row in self.result["rows"]:
            self.assertEqual(row["functional_calibration_status"], ADAPTER.NOT_FUNCTION_CALIBRATED)
        self.assertTrue(
            self.result["report"]["boundaries"]["functional_calibration_status_independent"]
        )

    # -- experimental grade ------------------------------------------------

    def test_no_grade_is_fabricated_for_gtdb_candidates(self):
        for row in self.result["rows"]:
            self.assertEqual(row["experimental_evidence_grade"], "")
        report = self.result["report"]
        self.assertEqual(report["counts"]["reference_rows_in_universe"], 0)
        self.assertTrue(report["boundaries"]["experimental_grade_fabricated"] is False)

    # -- mapping report ----------------------------------------------------

    def test_mapping_report_names_a_source_file_and_column_for_every_column(self):
        mapping = self.result["report"]["column_mapping"]
        for column in ADAPTER.OUTPUT_COLUMNS:
            self.assertIn(column, mapping, column)
            entry = mapping[column]
            self.assertTrue(entry["source"], column)
            self.assertTrue(entry["rule"], column)
            self.assertIn("pending_count", entry)
            self.assertIn("value_counts", entry)

    def test_mapping_report_counts_pending_values(self):
        pending = self.result["report"]["pending_value_classes"]
        self.assertIn("historical_superfamily", pending)
        self.assertEqual(pending["historical_superfamily"]["count"], 2)

    def test_unmappable_rows_are_reported_not_dropped(self):
        report = self.result["report"]
        self.assertIn("unmappable_rows", report)
        self.assertEqual(report["unmappable_rows"]["dropped_rows"], 0)

    # -- end-to-end through the catalog ------------------------------------

    def test_every_row_gets_exactly_one_primary_disposition(self):
        classified = CATALOG.classify_rows(self.result["rows"])
        self.assertEqual(len(classified), 8)
        dispositions = [row["primary_disposition"] for row in classified]
        from phaded_evidence_schema import PRIMARY_DISPOSITIONS
        for disposition in dispositions:
            self.assertIn(disposition, PRIMARY_DISPOSITIONS)
        self.assertEqual(len(dispositions), len(set(row["accession"] for row in classified)))

    def test_no_universe_row_is_excluded_for_input_quality(self):
        classified = CATALOG.classify_rows(self.result["rows"])
        self.assertEqual(
            [row["accession"] for row in classified
             if row["primary_disposition"] == "excluded_input_quality"],
            [],
        )
        for row in classified:
            self.assertEqual(row["input_quality_reason"], "")

    def test_scientific_uncertainty_rows_are_still_candidates(self):
        classified = {row["accession"]: row for row in CATALOG.classify_rows(self.result["rows"])}
        self.assertEqual(classified["A6"]["primary_disposition"], "remote_homolog_candidate")
        self.assertEqual(classified["A7"]["primary_disposition"], "remote_homolog_candidate")


class RealFrozenInputsTest(unittest.TestCase):
    """The frozen 109,087-row inputs, read-only (skipped when absent)."""

    UNIVERSE = Path("runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv")
    SUBTYPE = Path("runs/20260916_phaded_full_library_signalp_01/results/phaded_full_library_localization_subtype_matrix.tsv")

    def test_frozen_ahsmg_column_carries_no_positive_in_this_universe(self):
        """The AHSMG detector never fired for the frozen 109,087 candidate universe.

        This is the documented reason the adapted catalog emits zero AHSMG rows;
        it is asserted here so a future change of the frozen evidence cannot pass
        unnoticed.
        """
        if not self.UNIVERSE.is_file() or not self.SUBTYPE.is_file():
            self.skipTest("frozen inputs are not present in this checkout")
        with self.SUBTYPE.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 109087)
        positives = [row["accession"] for row in rows if row["ahsmg_state"] == "supported"]
        self.assertEqual(positives, [])


if __name__ == "__main__":
    unittest.main()
