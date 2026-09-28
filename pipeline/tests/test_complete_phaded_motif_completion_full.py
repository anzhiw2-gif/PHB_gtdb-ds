"""Invariants for the Task 7 full-library motif completion layer rerun.

Written before the implementation on purpose.  The invariants below are exactly
the ones that the 2026-09-15 candidate panel could not satisfy:

* every criterion must publish an ``assignable`` and a ``not_assignable`` count
  that add up to the candidate universe (109,087);
* a row that could not be evaluated must never be written as a negative
  (``not_detected_pattern`` / ``not_detected_in_tested_pfam`` /
  ``not_expected_for_subtype``);
* no detector may be fed a superfamily/subtype/family label, because the whole
  point of the candidate panel is that it can be measured against those priors;
* the reference-anchored residue transfer of Task 2 reaches 723 DED references
  and *no* candidate, so the candidate layer must say so instead of pretending
  the mapping carried over.
"""

import csv
import importlib.util
import inspect
import json
import re
import tempfile
import unittest
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "pipeline" / "scripts" / "complete_phaded_motif_completion_full.py"
TASK2_SCRIPT = REPO / "pipeline" / "scripts" / "complete_phaded_reference_residue_mapping.py"
RUN = REPO / "runs" / "20260917_phaded_motif_completion_full_01"
RESULTS = RUN / "results"

MOTIFS = REPO / "runs" / "20260915_phaded_motif_reconciliation_01"
LEGACY_EVIDENCE = MOTIFS / "results" / "motif_candidate_evidence.tsv"
CANDIDATE_FASTA = REPO / "runs" / "20260911_phaded_pfam_architecture_01" / "inputs" / "candidate_union.faa"
TASK2_RUN = REPO / "runs" / "20260917_phaded_reference_residue_mapping_01"
REFERENCE_MANIFEST = TASK2_RUN / "results" / "reference_mapping_manifest.tsv"
SUPERFAMILY_TABLE = TASK2_RUN / "results" / "superfamily_residue_mapping.tsv"
CALIBRATION_REPORT = TASK2_RUN / "results" / "implementation_calibration_report.json"
LITERATURE_REFERENCE = TASK2_RUN / "inputs" / "literature_residue_reference.tsv"
ANCHOR_FASTA = TASK2_RUN / "inputs" / "literature_anchor_sequences.faa"
REFERENCE_FASTA = REPO / "runs" / "20260911_phaded_pfam_architecture_01" / "inputs" / "phaded_reference.faa"
DED_ALIGNMENT_DIR = REPO / "runs" / "20260909_phaded_reference_evidence_01" / "inputs" / "raw" / "ded"

CANDIDATE_UNIVERSE = 109087
#: The five criteria this task was asked to recompute.
REQUIRED_CRITERIA = ("his", "asp", "lid", "linker", "sbd")
CRITERIA = ("lipase_box", "his", "asp", "oxyanion_hole", "ahsmg", "lid", "linker", "sbd")
#: Words that assert a result a row is not entitled to.  A ``not_assignable``
#: row carrying one of these would be exactly the failure mode this task exists
#: to prevent.
FORBIDDEN_STATE_WORDS = (
    "not_detected_pattern",
    "not_detected_in_tested_pfam",
    "not_detected_at_mapped_column",
    "not_expected_for_subtype",
    "not_expected_for_superfamily",
)

#: Measured 2026-09-15 candidate-layer distribution of the two criteria whose
#: coverage this task has to answer for.  The plan text quotes the *reference*
#: layer (His 0/723, Asp 25/723); these are the full-library numbers.
LEGACY_HIS_SUPPORTED = 265
LEGACY_ASP_SUPPORTED = 1290

LID_LOOP_1 = "FNGIG"
LID_LOOP_2 = "YYWQLF"


def _load(name: str, path: Path):
    if not path.is_file():
        raise AssertionError(f"required file is missing: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_module():
    return _load("complete_phaded_motif_completion_full", SCRIPT)


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


_SCAN: dict = {}


def full_library_scan() -> dict:
    """One cached pass over the produced layer, used by several assertions."""
    if not _SCAN:
        counters = {criterion: {"assignable": 0, "not_assignable": 0} for criterion in CRITERIA}
        violations: list[str] = []
        rows = 0
        with (RESULTS / "motif_completion_full.tsv").open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows += 1
                for criterion in CRITERIA:
                    bucket = row[f"{criterion}_assignability"]
                    if bucket not in counters[criterion]:
                        violations.append(f"{row['accession']}.{criterion}:unknown_bucket:{bucket}")
                        continue
                    counters[criterion][bucket] += 1
                    if bucket == "not_assignable":
                        state = row[f"{criterion}_state"]
                        if not state.startswith("not_assignable"):
                            violations.append(f"{row['accession']}.{criterion}:{state}")
                        for word in FORBIDDEN_STATE_WORDS:
                            if word in state:
                                violations.append(f"{row['accession']}.{criterion}:{state}")
                        if not row[f"{criterion}_assignability_reason"]:
                            violations.append(f"{row['accession']}.{criterion}:missing_reason")
        _SCAN.update(rows=rows, counters=counters, violations=violations)
    return _SCAN


def synthetic_sequences() -> dict[str, str]:
    """Six candidates covering every assignability branch of the layer."""
    lid_pair = "M" * 33 + LID_LOOP_1 + "A" * (190 - 39) + LID_LOOP_2 + "A" * 20
    return {
        # lid loops at 34-38 and 190-195 (spacing 156, the PhaZKT window)
        # plus a lipase box, a His motif, an Asp motif, a Cys oxyanion motif
        # and the AHSMG alternative motif
        "syn_lid_and_triad": "MADK" + lid_pair[4:] + "GLSAGG" + "GMTH" + "GATDYTV" + "HGCTQ" + "AHSMG",
        # nothing at all: every pattern must stay silent
        "syn_no_pattern": "M" + "P" * 80,
        # a single His-pattern hit in a circularly permuted (type 2) prior
        "syn_type2_his_only": "M" + "P" * 20 + "GMTH" + "P" * 20,
        # unassigned candidate that nonetheless carries PF06850
        "syn_unassigned_pf06850": "M" + "P" * 80,
        # unassigned candidate with no Pfam evidence at all
        "syn_no_pfam_evidence": "M" + "P" * 80,
        # lipase box only
        "syn_lipase_box_only": "M" + "K" * 10 + "GVSAG" + "K" * 10,
        # anchored prior (intracellular nPHAMCL) with no His/Asp pattern at all
        "syn_nphamcl_no_triad_pattern": "M" + "P" * 80,
    }


def synthetic_legacy_rows() -> list[dict]:
    return [
        {
            "accession": "syn_lid_and_triad",
            "genome": "GCF_syn_1",
            "motif_reference_subtype": "intracellular nPHAMCL",
            "assignment_status": "assigned",
            "pfam_accessions": "PF10503;PF06850",
            "pfam_coordinates": "PF10503:20-200;PF06850:220-300",
            "lipase_box_state": "supported",
            "his_state": "not_detected_pattern",
            "asp_state": "not_detected_pattern",
            "oxyanion_hole_state": "not_detected_pattern",
            "ahsmg_state": "not_detected_pattern",
            "lid_state": "pending_reference_annotation",
            "linker_state": "pending_reference_annotation",
            "sbd_state": "not_detected_in_tested_pfam",
        },
        {
            "accession": "syn_no_pattern",
            "genome": "GCF_syn_1",
            "motif_reference_subtype": "extracellular dPHASCL type 1",
            "assignment_status": "assigned",
            "pfam_accessions": "PF10503",
            "pfam_coordinates": "PF10503:5-120",
            "lipase_box_state": "not_detected_pattern",
            "his_state": "not_detected_pattern",
            "asp_state": "not_detected_pattern",
            "oxyanion_hole_state": "not_detected_pattern",
            "ahsmg_state": "not_detected_pattern",
            "lid_state": "not_expected_for_subtype",
            "linker_state": "pending_reference_annotation",
            "sbd_state": "not_detected_in_tested_pfam",
        },
        {
            "accession": "syn_type2_his_only",
            "genome": "GCF_syn_2",
            "motif_reference_subtype": "extracellular dPHASCL type 2",
            "assignment_status": "ambiguous_family",
            "pfam_accessions": "PF10503",
            "pfam_coordinates": "PF10503:2-90",
            "lipase_box_state": "not_detected_pattern",
            "his_state": "not_detected_pattern",
            "asp_state": "not_detected_pattern",
            "oxyanion_hole_state": "not_detected_pattern",
            "ahsmg_state": "not_detected_pattern",
            "lid_state": "not_expected_for_subtype",
            "linker_state": "pending_reference_annotation",
            "sbd_state": "not_detected_in_tested_pfam",
        },
        {
            "accession": "syn_unassigned_pf06850",
            "genome": "GCF_syn_3",
            "motif_reference_subtype": "",
            "assignment_status": "unassigned_PhaDED_like",
            "pfam_accessions": "PF06850",
            "pfam_coordinates": "PF06850:40-160",
            "lipase_box_state": "not_detected_pattern",
            "his_state": "not_detected_pattern",
            "asp_state": "not_detected_pattern",
            "oxyanion_hole_state": "not_detected_pattern",
            "ahsmg_state": "not_detected_pattern",
            "lid_state": "not_expected_for_subtype",
            "linker_state": "pending_reference_annotation",
            "sbd_state": "supported",
        },
        {
            "accession": "syn_no_pfam_evidence",
            "genome": "GCF_syn_4",
            "motif_reference_subtype": "",
            "assignment_status": "unassigned_PhaDED_like",
            "pfam_accessions": "",
            "pfam_coordinates": "",
            "lipase_box_state": "not_detected_pattern",
            "his_state": "not_detected_pattern",
            "asp_state": "not_detected_pattern",
            "oxyanion_hole_state": "not_detected_pattern",
            "ahsmg_state": "not_detected_pattern",
            "lid_state": "not_expected_for_subtype",
            "linker_state": "pending_reference_annotation",
            "sbd_state": "pending_tool_input",
        },
        {
            "accession": "syn_lipase_box_only",
            "genome": "GCF_syn_5",
            "motif_reference_subtype": "extracellular dPHAMCL",
            "assignment_status": "assigned",
            "pfam_accessions": "PF10503",
            "pfam_coordinates": "PF10503:1-60",
            "lipase_box_state": "supported",
            "his_state": "not_detected_pattern",
            "asp_state": "not_detected_pattern",
            "oxyanion_hole_state": "not_detected_pattern",
            "ahsmg_state": "not_detected_pattern",
            "lid_state": "not_expected_for_subtype",
            "linker_state": "not_expected_for_subtype",
            "sbd_state": "not_detected_in_tested_pfam",
        },
        {
            "accession": "syn_nphamcl_no_triad_pattern",
            "genome": "GCF_syn_6",
            "motif_reference_subtype": "intracellular nPHAMCL",
            "assignment_status": "ambiguous_family",
            "pfam_accessions": "PF10503",
            "pfam_coordinates": "PF10503:3-100",
            "lipase_box_state": "not_detected_pattern",
            "his_state": "not_detected_pattern",
            "asp_state": "not_detected_pattern",
            "oxyanion_hole_state": "not_detected_pattern",
            "ahsmg_state": "not_detected_pattern",
            "lid_state": "pending_reference_annotation",
            "linker_state": "pending_reference_annotation",
            "sbd_state": "not_detected_in_tested_pfam",
        },
    ]


class SyntheticFixture(unittest.TestCase):
    """Run the completion once on a six-candidate fixture."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.tmp = tempfile.TemporaryDirectory()
        root = Path(cls.tmp.name)
        cls.fasta = root / "candidates.faa"
        cls.fasta.write_text(
            "".join(f">{name}\n{sequence}\n" for name, sequence in synthetic_sequences().items()),
            encoding="ascii",
        )
        cls.evidence = root / "legacy_candidate_evidence.tsv"
        rows = synthetic_legacy_rows()
        fields: list[str] = []
        for row in rows:
            for field in row:
                if field not in fields:
                    fields.append(field)
        with cls.evidence.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        cls.output_dir = root / "results"
        cls.manifest = cls.module.run_motif_completion(
            candidate_fasta=cls.fasta,
            legacy_evidence=cls.evidence,
            reference_manifest=REFERENCE_MANIFEST,
            superfamily_table=SUPERFAMILY_TABLE,
            calibration_report=CALIBRATION_REPORT,
            literature_reference=LITERATURE_REFERENCE,
            reference_fasta=REFERENCE_FASTA,
            anchor_fasta=ANCHOR_FASTA,
            ded_alignment_dir=DED_ALIGNMENT_DIR,
            output_dir=cls.output_dir,
            run_id="synthetic_fixture",
        )
        cls.rows = read_tsv(cls.output_dir / "motif_completion_full.tsv")
        cls.by_accession = {row["accession"]: row for row in cls.rows}
        cls.summary = json.loads(
            (cls.output_dir / "motif_completion_summary.json").read_text(encoding="utf-8")
        )
        cls.coverage = read_tsv(cls.output_dir / "coverage_decomposition.tsv")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_outputs_are_produced_for_every_input_candidate(self):
        self.assertTrue((self.output_dir / "motif_completion_full.tsv").is_file())
        self.assertTrue((self.output_dir / "motif_completion_summary.json").is_file())
        self.assertTrue((self.output_dir / "coverage_decomposition.tsv").is_file())
        self.assertEqual(len(self.rows), len(synthetic_sequences()))
        self.assertEqual(
            sorted(self.by_accession), sorted(synthetic_sequences())
        )

    def test_assignable_plus_not_assignable_equals_the_candidate_count(self):
        totals = self.module.criterion_assignability_totals(self.rows)
        self.assertEqual(sorted(totals), sorted(CRITERIA))
        for criterion in CRITERIA:
            bucket = totals[criterion]
            self.assertEqual(
                bucket["assignable"] + bucket["not_assignable"],
                len(self.rows),
                f"{criterion} coverage must cover every candidate exactly once",
            )
            self.assertEqual(bucket["assignable"], bucket["assignable_positive"] + bucket["assignable_no_positive"])
        for criterion in REQUIRED_CRITERIA:
            self.assertIn(criterion, totals)

    def test_not_assignable_rows_never_carry_a_negative_state(self):
        for row in self.rows:
            for criterion in CRITERIA:
                if row[f"{criterion}_assignability"] == "not_assignable":
                    state = row[f"{criterion}_state"]
                    self.assertTrue(
                        state.startswith("not_assignable"),
                        f"{row['accession']}.{criterion} state {state!r}",
                    )
                    for word in FORBIDDEN_STATE_WORDS:
                        self.assertNotIn(word, state, f"{row['accession']}.{criterion}")
                    self.assertTrue(row[f"{criterion}_assignability_reason"])

    def test_assignable_rows_always_carry_a_reason(self):
        for row in self.rows:
            listed = [item for item in row["not_assignable_criteria"].split(";") if item]
            self.assertEqual(len(listed), len(set(listed)))
            for criterion in CRITERIA:
                self.assertIn(
                    row[f"{criterion}_assignability"], {"assignable", "not_assignable"}
                )
                self.assertTrue(row[f"{criterion}_assignability_reason"])
                self.assertEqual(
                    criterion in listed,
                    row[f"{criterion}_assignability"] == "not_assignable",
                    f"{row['accession']}.{criterion}",
                )
            self.assertEqual(int(row["criteria_not_assignable_count"]), len(listed))
            self.assertEqual(
                int(row["criteria_assignable_count"]) + len(listed), len(CRITERIA)
            )

    def test_detector_call_sites_never_receive_a_grouping_label(self):
        source = inspect.getsource(self.module)
        call_lines = [
            line for line in source.splitlines()
            if ("detect_literature_pattern(" in line or "detect_lid(" in line)
            and "def " not in line
        ]
        self.assertTrue(call_lines)
        for line in call_lines:
            for forbidden in ("subtype", "superfamily", "phaded_family", "prior", "family_id"):
                self.assertNotIn(forbidden, line, line)
        for function in (self.module.criterion_detection,):
            parameters = list(inspect.signature(function).parameters)
            self.assertLessEqual(set(parameters), {"sequence"})

    def test_his_and_asp_hits_are_assignable_and_misses_are_not(self):
        lid_row = self.by_accession["syn_lid_and_triad"]
        self.assertEqual(lid_row["his_state"], "supported")
        self.assertEqual(lid_row["his_assignability"], "assignable")
        self.assertEqual(lid_row["asp_state"], "supported")
        self.assertEqual(lid_row["asp_assignability"], "assignable")
        miss = self.by_accession["syn_no_pattern"]
        self.assertEqual(miss["his_assignability"], "not_assignable")
        self.assertEqual(miss["asp_assignability"], "not_assignable")
        self.assertTrue(miss["his_state"].startswith("not_assignable"))
        self.assertTrue(miss["his_assignability_basis"])

    def test_not_assignable_reason_depends_on_the_recorded_prior(self):
        # an anchor exists for the intracellular nPHAMCL prior, but no candidate
        # is a DED alignment member, so the transfer cannot reach it
        anchored = self.by_accession["syn_nphamcl_no_triad_pattern"]
        self.assertEqual(
            anchored["his_assignability_reason"],
            "not_assignable_candidate_not_in_any_ded_alignment",
        )
        # no literature anchor exists for the extracellular dPHASCL type 1 prior
        unanchored = self.by_accession["syn_no_pattern"]
        self.assertEqual(
            unanchored["his_assignability_reason"],
            "not_assignable_no_literature_anchor_for_prior_superfamily",
        )
        # the type 2 fold is circularly permuted: the anchored column is not transferable
        permuted = self.by_accession["syn_type2_his_only"]
        self.assertEqual(
            permuted["asp_assignability_reason"],
            "not_assignable_circularly_permuted_fold_column_not_transferable",
        )
        # no superfamily prior at all
        no_prior = self.by_accession["syn_no_pfam_evidence"]
        self.assertEqual(
            no_prior["his_assignability_reason"],
            "not_assignable_no_superfamily_prior_and_no_alignment_membership",
        )
        self.assertEqual(self.summary["unknown_prior_superfamilies"], [])

    def test_his_pattern_hit_stays_assignable_even_when_the_prior_is_permuted(self):
        row = self.by_accession["syn_type2_his_only"]
        self.assertEqual(row["his_state"], "supported")
        self.assertEqual(row["his_assignability"], "assignable")

    def test_lid_is_assignable_for_every_candidate_and_never_uses_the_prior(self):
        for row in self.rows:
            self.assertEqual(row["lid_assignability"], "assignable")
            self.assertIn(
                row["lid_state"],
                {
                    "supported",
                    "partial_lid_loop1_only",
                    "partial_lid_loop2_only",
                    "conflict_lid_loop_spacing",
                    "not_detected_pattern",
                },
            )
        row = self.by_accession["syn_lid_and_triad"]
        self.assertEqual(row["lid_state"], "supported")
        self.assertEqual(row["lid_coordinates"], "34-38;190-195")
        self.assertEqual(row["lid_loop_spacing"], "156")
        self.assertEqual(self.by_accession["syn_no_pattern"]["lid_state"], "not_detected_pattern")
        # the prior is reported in its own column and is never the detector input
        self.assertEqual(row["lid_expectation_from_prior"], "expected_for_prior")
        self.assertEqual(row["lid_prior_agrees_with_detection"], "true")

    def test_linker_and_sbd_cannot_be_assigned_at_the_candidate_layer(self):
        for row in self.rows:
            for criterion in ("linker", "sbd"):
                self.assertEqual(row[f"{criterion}_assignability"], "not_assignable")
                self.assertTrue(row[f"{criterion}_state"].startswith("not_assignable"))
        self.assertEqual(
            sorted({row["linker_assignability_reason"] for row in self.rows}),
            ["not_assignable_no_linker_profile_bound"],
        )
        self.assertEqual(
            sorted({row["sbd_assignability_reason"] for row in self.rows}),
            ["not_assignable_pf06850_is_not_the_extracellular_sbd_in_this_panel"],
        )

    def test_pf06850_raw_evidence_is_reported_separately_from_the_criterion(self):
        with_pf = self.by_accession["syn_unassigned_pf06850"]
        self.assertEqual(with_pf["sbd_pf06850_binding_state"], "detected")
        self.assertEqual(with_pf["sbd_pf06850_coordinates"], "40-160")
        self.assertEqual(with_pf["sbd_assignability"], "not_assignable")
        without = self.by_accession["syn_no_pfam_evidence"]
        self.assertEqual(without["sbd_pf06850_binding_state"], "not_tested_no_pfam_evidence")
        tested = self.by_accession["syn_no_pattern"]
        self.assertEqual(tested["sbd_pf06850_binding_state"], "not_detected_in_tested_pfam")
        for row in self.rows:
            self.assertTrue(row["sbd_pf06850_role_caveat"])

    def test_linker_region_delimitation_is_reported_but_not_promoted(self):
        row = self.by_accession["syn_lid_and_triad"]
        self.assertEqual(row["linker_region_delimitation"], "delimited_candidate_region_not_profiled")
        self.assertEqual(row["linker_region_coordinates"], "201-219")
        self.assertIn("no_linker_profile_bound", row["linker_region_basis"])
        self.assertEqual(row["linker_state"], "not_assignable_no_linker_profile_bound")
        self.assertEqual(
            self.by_accession["syn_no_pfam_evidence"]["linker_region_delimitation"],
            "no_pfam_annotation_available",
        )

    def test_lipase_box_and_oxyanion_misses_are_assignable_pattern_results(self):
        row = self.by_accession["syn_no_pattern"]
        self.assertEqual(row["lipase_box_state"], "not_detected_pattern")
        self.assertEqual(row["lipase_box_assignability"], "assignable")
        self.assertEqual(row["oxyanion_hole_state"], "not_detected_pattern")
        self.assertEqual(row["oxyanion_hole_assignability"], "assignable")
        self.assertEqual(self.by_accession["syn_lipase_box_only"]["lipase_box_state"], "supported")
        self.assertEqual(self.by_accession["syn_lid_and_triad"]["ahsmg_state"], "supported")

    def test_miss_assignability_is_derived_from_measured_pattern_coverage(self):
        basis = self.summary["criterion_miss_assignability_basis"]
        for criterion in CRITERIA:
            self.assertIn(criterion, basis)
        # the nucleophile pattern reproduces every literature anchor that has a
        # coordinate for it, the His and Asp proxies do not
        self.assertEqual(basis["lipase_box"]["anchors_reproduced"], basis["lipase_box"]["anchors_with_coordinate"])
        self.assertLess(basis["his"]["anchors_reproduced"], basis["his"]["anchors_with_coordinate"])
        self.assertLess(basis["asp"]["anchors_reproduced"], basis["asp"]["anchors_with_coordinate"])
        self.assertTrue(basis["his"]["miss_is_assignable"] is False)
        self.assertTrue(basis["asp"]["miss_is_assignable"] is False)
        self.assertTrue(basis["lipase_box"]["miss_is_assignable"] is True)
        self.assertTrue(basis["ahsmg"]["miss_is_assignable"] is True)
        self.assertGreater(basis["ahsmg"]["reference_panel_positive_controls"], 0)

    def test_summary_compares_before_and_after_for_every_criterion(self):
        criteria = self.summary["criteria"]
        self.assertEqual(sorted(criteria), sorted(CRITERIA))
        for criterion, payload in criteria.items():
            self.assertIn("before", payload)
            self.assertIn("after", payload)
            self.assertIn("assignable", payload)
            self.assertIn("not_assignable", payload)
            self.assertEqual(
                payload["assignable"] + payload["not_assignable"], len(self.rows), criterion
            )
            self.assertTrue(payload["before_source"])
            self.assertTrue(payload["after_source"])

    def test_coverage_decomposition_rows_cover_every_candidate_once(self):
        totals = self.module.criterion_assignability_totals(self.rows)
        counted: dict[str, int] = {}
        for row in self.coverage:
            counted.setdefault(row["criterion"], 0)
            counted[row["criterion"]] += int(row["candidate_count"])
        self.assertEqual(sorted(counted), sorted(CRITERIA))
        for criterion, total in counted.items():
            self.assertEqual(total, len(self.rows), criterion)
            self.assertEqual(
                totals[criterion]["assignable"], int(
                    next(
                        item["assignable"]
                        for item in self.coverage
                        if item["criterion"] == criterion and item["bucket"] == "assignable"
                    )
                )
            )

    def test_every_row_states_that_no_new_family_call_is_made(self):
        boundary = self.module.BOUNDARY
        for row in self.rows:
            self.assertEqual(row["motif_phenotype_boundary"], boundary)
            self.assertEqual(row["motif_evidence_level"], "full_library_sequence_pattern_layer_only")
            self.assertEqual(row["new_family_call_made"], "false")
        self.assertTrue(self.summary["phenotype_boundary_holds"])
        self.assertFalse(self.summary["new_family_or_superfamily_call_made"])


class FullLibraryArtifactTests(unittest.TestCase):
    """Check the produced 109,087-row layer itself."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.summary = json.loads(
            (RESULTS / "motif_completion_summary.json").read_text(encoding="utf-8")
        )
        cls.coverage = read_tsv(RESULTS / "coverage_decomposition.tsv")

    def test_the_run_outputs_exist_and_declare_the_candidate_universe(self):
        for name in ("motif_completion_full.tsv", "motif_completion_summary.json", "coverage_decomposition.tsv"):
            self.assertTrue((RESULTS / name).is_file(), name)
        self.assertEqual(self.summary["candidate_count"], CANDIDATE_UNIVERSE)
        self.assertEqual(self.summary["run_id"], RUN.name)

    def test_full_library_assignability_counts_sum_to_the_candidate_universe(self):
        for criterion, payload in self.summary["criteria"].items():
            self.assertEqual(
                payload["assignable"] + payload["not_assignable"],
                CANDIDATE_UNIVERSE,
                criterion,
            )
        counted = {}
        for row in self.coverage:
            counted.setdefault(row["criterion"], 0)
            counted[row["criterion"]] += int(row["candidate_count"])
        self.assertEqual(sorted(counted), sorted(CRITERIA))
        for criterion in CRITERIA:
            self.assertEqual(counted[criterion], CANDIDATE_UNIVERSE, criterion)

    def test_the_full_library_tsv_never_writes_a_not_assignable_row_as_a_negative(self):
        scan = full_library_scan()
        self.assertEqual(scan["rows"], CANDIDATE_UNIVERSE)
        self.assertEqual(scan["violations"], [])

    def test_full_library_counts_match_the_independent_recount(self):
        """Recount the TSV state columns without trusting the summary."""
        scan = full_library_scan()
        for criterion in CRITERIA:
            payload = self.summary["criteria"][criterion]
            self.assertEqual(payload["assignable"], scan["counters"][criterion]["assignable"], criterion)
            self.assertEqual(
                payload["not_assignable"], scan["counters"][criterion]["not_assignable"], criterion
            )

    def test_the_legacy_candidate_layer_is_quoted_as_measured(self):
        """The plan quotes the reference layer; the candidate layer is different."""
        criteria = self.summary["criteria"]
        self.assertEqual(criteria["his"]["before"]["supported"], LEGACY_HIS_SUPPORTED)
        self.assertEqual(criteria["asp"]["before"]["supported"], LEGACY_ASP_SUPPORTED)
        self.assertEqual(criteria["lid"]["before"]["independent_detection_supported"], 0)
        self.assertIn("motif_candidate_evidence.tsv", criteria["his"]["before_source"])
        corrections = self.summary["premise_corrections"]
        self.assertTrue(any("his" in note for note in corrections))
        self.assertTrue(any("asp" in note for note in corrections))

    def test_his_coverage_improved_and_asp_pattern_is_unchanged(self):
        criteria = self.summary["criteria"]
        self.assertGreater(criteria["his"]["after"]["supported"], criteria["his"]["before"]["supported"])
        self.assertEqual(
            criteria["asp"]["after"]["pattern_supported"], criteria["asp"]["before"]["supported"]
        )
        self.assertEqual(criteria["his"]["after"]["supported"], criteria["his"]["assignable"])
        self.assertEqual(criteria["asp"]["after"]["supported"], criteria["asp"]["assignable"])

    def test_lid_went_from_zero_independent_detections_to_full_coverage(self):
        lid = self.summary["criteria"]["lid"]
        self.assertEqual(lid["before"]["independent_detection_supported"], 0)
        self.assertEqual(lid["assignable"], CANDIDATE_UNIVERSE)
        self.assertGreater(lid["after"]["supported"], 0)
        self.assertEqual(
            lid["after"]["supported"] + lid["after"]["partial"], lid["after"]["any_loop_detected"]
        )

    def test_linker_and_sbd_remain_fully_unassignable(self):
        for criterion in ("linker", "sbd"):
            self.assertEqual(self.summary["criteria"][criterion]["assignable"], 0)
            self.assertEqual(
                self.summary["criteria"][criterion]["not_assignable"], CANDIDATE_UNIVERSE
            )
        self.assertEqual(self.summary["criteria"]["sbd"]["after"]["pf06850_detected"], 30621)

    def test_no_candidate_is_a_ded_alignment_member(self):
        self.assertEqual(self.summary["candidate_ded_alignment_members"], 0)
        self.assertEqual(
            self.summary["reference_alignment_column_transfer"]["references_with_anchor_column_transfer"],
            65,
        )

    def test_summary_records_the_censused_inputs_with_hashes(self):
        for name, record in self.summary["inputs"].items():
            self.assertTrue(record["path"], name)
            self.assertTrue(re.match(r"^[0-9a-f]{64}$", record["sha256"]), name)
            self.assertGreater(record["size"], 0, name)
        self.assertEqual(self.summary["candidate_alignment_layer"], "sequence_pattern_only_no_alignment_membership")
        self.assertEqual(self.summary["status"], "completed_candidate_only")

    def test_lid_prior_agreement_is_reported_without_using_the_prior_as_input(self):
        agreement = self.summary["lid_prior_agreement"]
        self.assertEqual(
            agreement["detected_and_expected"]
            + agreement["detected_but_not_expected"]
            + agreement["not_detected_but_expected"]
            + agreement["not_detected_and_not_expected"],
            CANDIDATE_UNIVERSE,
        )
        self.assertEqual(
            agreement["evidence_level"],
            "independent_sequence_detection_compared_to_a_recorded_prior",
        )

    def test_status_document_exists_and_states_the_boundary(self):
        document = REPO / "docs" / "T141_20260917_phaded_motif_completion_full_status.md"
        self.assertTrue(document.is_file())
        text = document.read_text(encoding="utf-8")
        for required in (
            "partial_catalytic_pattern_panel",
            "不产生任何新的 family 判定",
            "不可判定不等于阴性",
            "109,087",
        ):
            self.assertIn(required, text)


class BeforeBlockFieldNamingTests(unittest.TestCase):
    """V5: the ``before`` block must not look like an independent-detection count.

    The frozen ``motif_completion_summary.json`` of this run carries the total
    supported count of each criterion under the name
    ``before.independent_detection_supported``, and hard-codes ``0`` for ``lid``,
    ``linker`` and ``sbd``.  The name is over-strong (those numbers are the
    *resolved* total, whatever the evidence origin) and the literal ``0`` reads as
    "a detector ran and found nothing" for three criteria where no independent
    detector existed at all.  The frozen artifact cannot be rewritten, so this
    test pins the corrected vocabulary of any future run and the recommendation
    recorded beside the frozen file.
    """

    LEGACY_FIELD = "independent_detection_supported"
    RESOLVED_FIELD = "before_resolved_supported"
    INDEPENDENT_FIELD = "before_independent_detection_supported"
    NOT_AVAILABLE_FIELD = "before_not_available"
    NOT_AVAILABLE_REASON_FIELD = "before_not_available_reason"
    STATE_COUNTS_FIELD = "state_counts"
    NOT_AVAILABLE = "not_available"
    #: criteria whose 2026-09-15 column never came from an independent detector
    NO_INDEPENDENT_DETECTOR = ("lid", "linker", "sbd")

    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.before = cls.module.legacy_state_distributions(read_tsv(LEGACY_EVIDENCE))
        counted: dict[str, Counter] = {criterion: Counter() for criterion in CRITERIA}
        with LEGACY_EVIDENCE.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for criterion in CRITERIA:
                    counted[criterion][row[f"{criterion}_state"]] += 1
        cls.counted = counted

    def partition(self, criterion: str) -> dict:
        return dict(self.before[criterion][self.STATE_COUNTS_FIELD])

    def test_before_block_drops_the_over_strong_legacy_field_name(self):
        for criterion in CRITERIA:
            self.assertNotIn(self.LEGACY_FIELD, self.before[criterion], criterion)

    def test_partition_and_aggregates_are_kept_in_separate_namespaces(self):
        """A flat dict mixed state names with aggregates (V5's ambiguity)."""
        for criterion in CRITERIA:
            payload = self.before[criterion]
            self.assertEqual(
                sorted(payload),
                sorted(
                    [
                        self.STATE_COUNTS_FIELD,
                        self.RESOLVED_FIELD,
                        self.INDEPENDENT_FIELD,
                        self.NOT_AVAILABLE_FIELD,
                        self.NOT_AVAILABLE_REASON_FIELD,
                    ]
                ),
                criterion,
            )
            for field in (
                self.RESOLVED_FIELD,
                self.INDEPENDENT_FIELD,
                self.NOT_AVAILABLE_FIELD,
                self.NOT_AVAILABLE_REASON_FIELD,
            ):
                self.assertNotIn(field, payload[self.STATE_COUNTS_FIELD], field)

    def test_before_resolved_supported_is_the_independently_recounted_total(self):
        for criterion in CRITERIA:
            measured = sum(
                count for state, count in self.counted[criterion].items()
                if state.startswith("supported")
            )
            self.assertEqual(
                self.before[criterion][self.RESOLVED_FIELD], measured, criterion
            )
            self.assertEqual(
                measured,
                sum(
                    count for state, count in self.partition(criterion).items()
                    if state.startswith("supported")
                ),
                criterion,
            )

    def test_before_state_partition_is_unchanged_and_separately_summable(self):
        for criterion in CRITERIA:
            self.assertEqual(
                self.partition(criterion),
                dict(sorted(self.counted[criterion].items())),
                criterion,
            )
            self.assertEqual(sum(self.partition(criterion).values()), CANDIDATE_UNIVERSE, criterion)

    def test_criteria_without_an_independent_before_detector_are_not_available(self):
        for criterion in self.NO_INDEPENDENT_DETECTOR:
            payload = self.before[criterion]
            self.assertTrue(payload[self.NOT_AVAILABLE_FIELD], criterion)
            self.assertEqual(payload[self.INDEPENDENT_FIELD], self.NOT_AVAILABLE, criterion)
            self.assertTrue(payload[self.NOT_AVAILABLE_REASON_FIELD], criterion)
            # the whole point: no literal 0 pretending a detector ran
            self.assertNotEqual(payload[self.INDEPENDENT_FIELD], 0, criterion)

    def test_measured_criteria_report_an_integer_independent_before_count(self):
        measured_criteria = [
            criterion for criterion in CRITERIA
            if criterion not in self.NO_INDEPENDENT_DETECTOR
        ]
        self.assertEqual(len(measured_criteria), 5)
        for criterion in measured_criteria:
            payload = self.before[criterion]
            self.assertFalse(payload[self.NOT_AVAILABLE_FIELD], criterion)
            self.assertEqual(payload[self.NOT_AVAILABLE_REASON_FIELD], "", criterion)
            independence = payload[self.INDEPENDENT_FIELD]
            self.assertIsInstance(independence, int, criterion)
            self.assertEqual(independence, payload[self.RESOLVED_FIELD], criterion)

    def test_field_semantics_documents_the_before_vocabulary(self):
        semantics = self.module.BEFORE_FIELD_SEMANTICS
        for field in (
            self.RESOLVED_FIELD,
            self.INDEPENDENT_FIELD,
            self.NOT_AVAILABLE_FIELD,
            self.NOT_AVAILABLE_REASON_FIELD,
        ):
            self.assertIn(field, semantics, field)
            self.assertTrue(semantics[field], field)
        self.assertIn(self.LEGACY_FIELD, semantics)
        self.assertIn("removed", semantics[self.LEGACY_FIELD].lower())

    def test_before_distribution_docstring_names_the_new_fields(self):
        docstring = inspect.getdoc(self.module.legacy_state_distributions) or ""
        for field in (self.RESOLVED_FIELD, self.INDEPENDENT_FIELD, self.NOT_AVAILABLE_FIELD):
            self.assertIn(field, docstring, field)


if __name__ == "__main__":
    unittest.main()
