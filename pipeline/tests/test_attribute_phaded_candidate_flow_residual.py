"""Tests for the 616-accession candidate-flow residual attribution (Task F1).

Written before the implementation (test-first).  The module under test explains
the ``109,087 - 36,611 - 71,860 = 616`` residual of the frozen v1 candidate flow
by attributing every residual accession to the stage that provably dropped it,
and by proposing a disposition *rule* per attribution.

Load-bearing invariants asserted here:

* the residual is an explicit accession **set difference** -- never arithmetic on
  the three published aggregates;
* every residual accession carries exactly one ``attribution``, and the per
  attribution counts sum to the residual size;
* an accession whose frozen evidence does not prove a stage is reported as
  ``unexplained`` -- it is never silently pushed into a bucket;
* a duplicate accession in the residual is an error, not a de-duplication;
* a missing input path is an error, never an empty input;
* the script never writes outside ``--output-dir`` and refuses a non-empty one;
* the intentional ``check_hold_residual`` mismatch is reported, never silenced;
* the attribution logic reuses ``reconcile_phaded_candidate_flow`` rather than
  re-implementing its set logic.

All fixtures are synthetic.  Nothing in this module reads or writes ``runs/``,
``results/`` or ``deploy/``.
"""

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = TESTS_DIR.parent / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

SCRIPT = SCRIPTS_DIR / "attribute_phaded_candidate_flow_residual.py"

REGISTRY_TYPE1 = "extracellular dPHASCL type 1"
REGISTRY_TYPE2 = "extracellular dPHASCL type 2"
REGISTRY_CYS = "intracellular nPHASCL without lipase box"
REGISTRY_LIPASE = "intracellular nPHASCL with lipase box"
OUT_OF_REGISTRY = "some unregistered superfamily"


def load_module():
    """Load the script under test by path (it is not an importable package)."""
    if not SCRIPT.is_file():
        raise FileNotFoundError("script under test is missing: %s" % SCRIPT)
    spec = importlib.util.spec_from_file_location(
        "attribute_phaded_candidate_flow_residual", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in fields) + "\n")
    return path


def universe_row(entry):
    return {
        "accession": entry["accession"],
        "genome": entry["accession"].split("|")[0],
        "family": entry.get("family", "iPhaZ"),
        "source_layer": entry.get("source_layer", "tier2"),
        "signal_type": entry.get("signal_type", "OTHER"),
        "assignment": entry.get("assignment", ""),
        "length": entry.get("length", "300"),
        "layer": entry.get("layer", "iPhaZ_tier2_review"),
    }


def flow_rows(entries, kind):
    """Rows for the high-confidence (``kind='hc'``) or hold (``kind='hold'``) table."""
    if kind == "hc":
        fields = ["accession", "superfamily", "high_confidence"]
        return [
            {
                "accession": entry["accession"],
                "superfamily": entry.get("flow_superfamily", REGISTRY_TYPE1),
                "high_confidence": "1",
            }
            for entry in entries
            if entry.get("flow") == "hc"
        ], fields
    fields = ["accession", "superfamily", "hold_reason"]
    return [
        {
            "accession": entry["accession"],
            "superfamily": entry.get("flow_superfamily", REGISTRY_TYPE1),
            "hold_reason": entry.get("hold_reason", "common_criteria"),
        }
        for entry in entries
        if entry.get("flow") == "hold"
    ], fields


MATRIX_FIELDS = [
    "accession",
    "assignment_status",
    "evidence_status",
    "evidence_threshold",
    "phaded_superfamily_best",
    "superfamily_score_gap",
    "profile_evidence_status",
    "sequence_integrity",
    "evidence_grade",
]

CLASSIFICATION_FIELDS = [
    "protein_id",
    "superfamily_claim",
    "superfamily_confidence",
    "n_superfamilies_hit",
    "best_evalue",
    "model_layer",
]

MOTIF_FIELDS = [
    "accession",
    "candidate_prior_superfamily",
    "candidate_assignment_status",
    "prior_anchor_basis",
    "lipase_box_state",
    "sbd_pf06850_binding_state",
]

CONFOUNDER_FIELDS = ["accession", "confounder_flag", "recommended_confidence_cap"]


def build_fixture(root, entries, registry=(REGISTRY_TYPE1, REGISTRY_TYPE2, REGISTRY_CYS, REGISTRY_LIPASE),
                  pool_external_layers=None):
    """Write a complete synthetic frozen-input set and return the CLI paths."""
    root = Path(root)
    inputs = root / "frozen"
    inputs.mkdir(parents=True, exist_ok=True)

    universe = write_tsv(
        inputs / "protein_layers.tsv",
        ["accession", "genome", "family", "source_layer", "signal_type", "assignment", "length", "layer"],
        [universe_row(entry) for entry in entries if entry.get("in_universe", True)],
    )
    hc_rows, hc_fields = flow_rows(entries, "hc")
    high_confidence = write_tsv(inputs / "high_confidence_candidates.tsv", hc_fields, hc_rows)
    hold_rows, hold_fields = flow_rows(entries, "hold")
    hold = write_tsv(inputs / "hold_candidates.tsv", hold_fields, hold_rows)

    matrix_rows = []
    classification_rows = []
    motif_rows = []
    confounder_rows = []
    for entry in entries:
        matrix = entry.get("matrix")
        if matrix is None and entry.get("flow") in {"hc", "hold"}:
            # a flow row implies the frozen filter found an admissible assignment
            matrix = dict(ASSIGNED_MATRIX)
        if matrix is not None:
            row = {
                "accession": entry["accession"],
                "assignment_status": "unassigned_PhaDED_like",
                "evidence_status": "no_profile_score",
                "evidence_threshold": "min_score_gap_1.0_bits_HMMER_E1e6",
                "phaded_superfamily_best": "",
                "superfamily_score_gap": "0.0",
                "profile_evidence_status": "profile_unassigned",
                "sequence_integrity": "terminal_stop_only",
                "evidence_grade": "PENDING_profile_unassigned",
            }
            row.update(matrix)
            matrix_rows.append(row)
        classification = entry.get("classification")
        if classification is not None:
            row = {
                "protein_id": entry["accession"],
                "superfamily_claim": "",
                "superfamily_confidence": "unique",
                "n_superfamilies_hit": "1",
                "best_evalue": "1e-40",
                "model_layer": "discovery_hmm_uncalibrated",
            }
            row.update(classification)
            classification_rows.append(row)
        motif = entry.get("motif")
        if motif is not None:
            row = {
                "accession": entry["accession"],
                "candidate_prior_superfamily": "",
                "candidate_assignment_status": "unassigned_PhaDED_like",
                "prior_anchor_basis": "no_superfamily_prior",
                "lipase_box_state": "not_detected_pattern",
                "sbd_pf06850_binding_state": "not_detected_in_tested_pfam",
            }
            row.update(motif)
            motif_rows.append(row)
        confounder = entry.get("confounder")
        if confounder is not None:
            row = {
                "accession": entry["accession"],
                "confounder_flag": "suspected_lipase_esterase_confounder",
                "recommended_confidence_cap": "moderate_candidate_only",
            }
            row.update(confounder)
            confounder_rows.append(row)

    matrix = write_tsv(inputs / "subtype_matrix.tsv", MATRIX_FIELDS, matrix_rows) if matrix_rows else None
    classification = (
        write_tsv(inputs / "superfamily_classification.tsv", CLASSIFICATION_FIELDS, classification_rows)
        if classification_rows
        else None
    )
    motif = write_tsv(inputs / "motif_completion_full.tsv", MOTIF_FIELDS, motif_rows) if motif_rows else None
    confounders = (
        write_tsv(inputs / "confounder_candidates.tsv", CONFOUNDER_FIELDS, confounder_rows)
        if confounder_rows
        else None
    )

    per_superfamily = {}
    for name in registry:
        pass_rows = [e for e in entries if e.get("flow") == "hc" and e.get("flow_superfamily", REGISTRY_TYPE1) == name]
        hold_selected = [e for e in entries if e.get("flow") == "hold" and e.get("flow_superfamily", REGISTRY_TYPE1) == name]
        total = len(pass_rows) + len(hold_selected)
        if not total:
            continue
        per_superfamily[name] = {
            "total": total,
            "pass": len(pass_rows),
            "fail_localization": sum(1 for e in hold_selected if e.get("hold_reason") == "localization_conflict"),
            "fail_architecture": sum(1 for e in hold_selected if e.get("hold_reason") == "architecture_criteria"),
            "fail_common": sum(1 for e in hold_selected if e.get("hold_reason") == "common_criteria"),
        }
    hold_breakdown = {}
    for entry in entries:
        if entry.get("flow") == "hold":
            reason = entry.get("hold_reason", "common_criteria")
            hold_breakdown[reason] = hold_breakdown.get(reason, 0) + 1
    filter_summary = inputs / "filter_summary.json"
    filter_summary.write_text(
        json.dumps(
            {
                "total_high_confidence": len(hc_rows),
                "total_held": len(hold_rows),
                "hold_breakdown": hold_breakdown,
                "per_superfamily": per_superfamily,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    layer_paths = {}
    for name, rows in (pool_external_layers or {}).items():
        layer_paths[name] = write_tsv(
            inputs / ("pool_external_%s.tsv" % name),
            ["protein_id", "superfamily", "score_tier"],
            rows,
        )
    return {
        "universe": universe,
        "high_confidence": high_confidence,
        "hold": hold,
        "filter_summary": filter_summary,
        "matrix": matrix,
        "classification": classification,
        "motif": motif,
        "confounders": confounders,
        "pool_external_layers": layer_paths,
    }


UNASSIGNED = {
    "matrix": {},  # all matrix defaults: unassigned_PhaDED_like, no profile score
    "motif": {},
}

#: The frozen column values of a flow-row accession: the trained-profile branch of
#: the frozen filter admits it, which is why it reaches the high-confidence or hold
#: table at all.
ASSIGNED_MATRIX = {
    "assignment_status": "assigned",
    "evidence_status": "profile_score",
    "phaded_superfamily_best": REGISTRY_TYPE1,
    "profile_evidence_status": "profile_trained_hit",
    "superfamily_score_gap": "12.5",
    "evidence_grade": "L1_profile_plus_partial_architecture",
}


def residual_of(module, paths):
    """Read the frozen fixture inputs and return the residual accession set."""
    universe = set(module.reconcile.read_accession_list(paths["universe"], "universe"))
    high_confidence = set(
        module.reconcile.read_accession_list(paths["high_confidence"], "high confidence")
    )
    hold = set(module.read_hold_accessions(paths["hold"]))
    return sorted(universe - high_confidence - hold)


class AttributionVocabularyTests(unittest.TestCase):
    def test_registry_range_comes_from_the_frozen_v1_filter(self):
        module = load_module()
        frozen = module.load_frozen_filter_registry()
        self.assertEqual(len(frozen), 8)
        self.assertIn(REGISTRY_CYS, frozen)
        self.assertIn("periplasmic PHA depolymerases", frozen)
        self.assertEqual(frozen, module.frozen_filter_registry_for_testing())

    def test_vocabulary_keeps_the_requested_values_and_documents_extensions(self):
        module = load_module()
        for name in (
            "tier_below_threshold",
            "outside_registry_family_range",
            "pool_external_boundary",
            "score_tier_not_upgraded",
            "prefilter_excluded",
            "table_construction_gap",
            "unexplained",
        ):
            self.assertIn(name, module.ATTRIBUTION_ORDER, name)
            self.assertNotIn(name, module.ATTRIBUTION_EXTENSIONS, name)
        for name in ("profile_unassigned_no_discovery_claim", "discovery_claim_not_unique"):
            self.assertIn(name, module.ATTRIBUTION_ORDER, name)
            self.assertIn(name, module.ATTRIBUTION_EXTENSIONS, name)
            self.assertIn("EXTENSION", module.ATTRIBUTION_BASIS[name])
        self.assertEqual(
            len(module.ATTRIBUTION_ORDER), len(set(module.ATTRIBUTION_ORDER))
        )
        for name in module.ATTRIBUTION_ORDER:
            self.assertIn(name, module.ATTRIBUTION_STATUS, name)

    def test_every_attribution_has_a_disposition_rule_with_proving_evidence(self):
        module = load_module()
        for name in module.ATTRIBUTION_ORDER:
            rule = module.disposition_rule(name)
            self.assertIn(rule["primary_disposition"], module.reconcile.DISPOSITIONS, name)
            self.assertTrue(rule["justification"].strip(), name)
            self.assertTrue(rule["proving_evidence"].strip(), name)
            self.assertIn("::", rule["proving_evidence"], name)


class AttributeResidualTests(unittest.TestCase):
    #: Frozen column values of a residual accession that never received a call.
    UNASSIGNED_EVIDENCE = {
        "matrix_present": "true",
        "matrix_assignment_status": "unassigned_PhaDED_like",
        "matrix_evidence_status": "no_profile_score",
        "matrix_evidence_threshold": "min_score_gap_1.0_bits_HMMER_E1e6",
        "matrix_superfamily_score_gap": "0.0",
        "matrix_profile_evidence_status": "profile_unassigned",
    }

    def setUp(self):
        self.module = load_module()

    def evidence(self, accession, **overrides):
        values = dict(self.UNASSIGNED_EVIDENCE)
        values.update(overrides)
        return self.module.make_evidence_row(accession, **values)

    def test_duplicate_residual_accessions_raise(self):
        rows = [self.evidence("A_1"), self.evidence("A_2")]
        with self.assertRaises(ValueError):
            self.module.attribute_residual(["A_1", "A_1"], rows, {}, None)

    def test_empty_residual_is_reported_not_an_error(self):
        rows, summary = self.module.attribute_residual([], [], {}, None)
        self.assertEqual(rows, [])
        self.assertEqual(summary["residual_count"], 0)
        self.assertEqual(summary["sum_attribution_counts"], 0)
        self.assertTrue(summary["attribution_counts_sum_to_residual"])

    def test_attribution_counts_sum_to_residual_size(self):
        rows = [
            self.evidence("A_1"),  # unassigned, no discovery row
            self.evidence(
                "A_2",
                discovery_superfamily_confidence="ambiguous_2",
                discovery_n_superfamilies_hit="2",
            ),
            self.evidence(
                "A_3",
                matrix_assignment_status="ambiguous_superfamily",
                matrix_evidence_status="profile_score",
                matrix_phaded_superfamily_best=REGISTRY_TYPE1,
                matrix_profile_evidence_status="profile_ambiguous_superfamily",
                matrix_superfamily_score_gap="0.4",
            ),
        ]
        residual = ["A_1", "A_2", "A_3"]
        attributed, summary = self.module.attribute_residual(residual, rows, {}, None)
        self.assertEqual(len(attributed), 3)
        self.assertEqual(summary["sum_attribution_counts"], 3)
        self.assertEqual(summary["residual_count"], 3)
        self.assertTrue(summary["attribution_counts_sum_to_residual"])
        self.assertEqual(summary["attribution_counts"]["profile_unassigned_no_discovery_claim"], 1)
        self.assertEqual(summary["attribution_counts"]["discovery_claim_not_unique"], 1)
        self.assertEqual(summary["attribution_counts"]["score_tier_not_upgraded"], 1)
        self.assertEqual(summary["unexplained_count"], 0)

    def test_unattributable_accession_is_explicitly_unexplained(self):
        # a row whose frozen columns match no documented mechanism: a discovery
        # claim outside the registry range is attributed, so use a unique in-range
        # claim whose flow row is missing -> that is a table construction gap, and
        # a row with NO evidence at all -> unexplained.
        rows = [self.evidence("A_9", matrix_evidence_status="something_new")]
        attributed, summary = self.module.attribute_residual(["A_9"], rows, {}, None)
        self.assertEqual(attributed[0]["attribution"], "unexplained")
        self.assertEqual(summary["unexplained_count"], 1)
        self.assertIn("A_9", summary["unexplained_accessions"])
        self.assertEqual(summary["sum_attribution_counts"], 1)

    def test_missing_evidence_row_is_unexplained_never_bucketed(self):
        rows = [self.evidence("A_1")]
        attributed, summary = self.module.attribute_residual(["A_1", "A_2"], rows, {}, None)
        by_accession = {row["accession"]: row for row in attributed}
        self.assertEqual(by_accession["A_1"]["attribution"], "profile_unassigned_no_discovery_claim")
        self.assertEqual(by_accession["A_2"]["attribution"], "unexplained")
        self.assertIn("no evidence row", by_accession["A_2"]["attribution_basis"])
        self.assertEqual(summary["attribution_counts"]["unexplained"], 1)
        self.assertEqual(summary["sum_attribution_counts"], 2)

    def test_extra_evidence_rows_are_reported_not_silently_dropped(self):
        rows = [self.evidence("A_1"), self.evidence("A_7")]
        attributed, summary = self.module.attribute_residual(["A_1"], rows, {}, None)
        self.assertEqual(len(attributed), 1)
        self.assertEqual(summary["evidence_rows_outside_residual"], 1)
        self.assertEqual(summary["evidence_rows_total"], 2)

    def test_pool_external_layer_is_attributed_as_pool_external_boundary(self):
        rows = [self.evidence("A_1"), self.evidence("A_2")]
        external = {"pool_external_strong": ["A_2"]}
        attributed, summary = self.module.attribute_residual(
            ["A_1", "A_2"], rows, {}, external
        )
        by_accession = {row["accession"]: row for row in attributed}
        self.assertEqual(by_accession["A_2"]["attribution"], "pool_external_boundary")
        self.assertEqual(by_accession["A_2"]["pool_external_layers"], "pool_external_strong")
        self.assertEqual(summary["attribution_counts"]["pool_external_boundary"], 1)
        self.assertEqual(summary["pool_external_layers_supplied"], ["pool_external_strong"])
        attributed = {row["accession"]: row for row in attributed}
        self.assertEqual(attributed["A_1"]["pool_external_layers"], "")

    def test_absent_pool_external_layer_is_pending_not_a_proven_zero(self):
        rows = [self.evidence("A_1")]
        _, summary = self.module.attribute_residual(["A_1"], rows, {}, None)
        self.assertEqual(summary["pool_external_layers_supplied"], [])
        self.assertIn("pool_external_boundary", summary["pending_buckets"])

    def test_outside_registry_family_range_attribution(self):
        rows = [
            self.evidence("A_1", discovery_superfamily_claim=OUT_OF_REGISTRY),
            self.evidence("A_2", matrix_phaded_superfamily_best=OUT_OF_REGISTRY),
        ]
        attributed, summary = self.module.attribute_residual(["A_1", "A_2"], rows, {}, None)
        self.assertEqual(
            [row["attribution"] for row in attributed],
            ["outside_registry_family_range", "outside_registry_family_range"],
        )
        self.assertEqual(summary["attribution_counts"]["outside_registry_family_range"], 2)

    def test_unique_in_range_claim_is_a_table_construction_gap(self):
        rows = [
            self.evidence(
                "A_1",
                discovery_superfamily_claim=REGISTRY_TYPE1,
                discovery_superfamily_confidence="unique",
                discovery_n_superfamilies_hit="1",
            )
        ]
        attributed, _ = self.module.attribute_residual(["A_1"], rows, {}, None)
        self.assertEqual(attributed[0]["attribution"], "table_construction_gap")

    def test_confounder_is_a_flag_not_a_primary_attribution(self):
        rows = [self.evidence("A_1", confounder_flag="suspected_lipase_esterase_confounder")]
        attributed, summary = self.module.attribute_residual(["A_1"], rows, {}, None)
        self.assertEqual(attributed[0]["attribution"], "profile_unassigned_no_discovery_claim")
        self.assertEqual(attributed[0]["confounder_flag"], "suspected_lipase_esterase_confounder")
        self.assertEqual(summary["attribution_counts"]["prefilter_excluded"], 0)
        self.assertEqual(summary["confounder_flagged_count"], 1)

    def test_tier_below_threshold_is_documented_as_unreachable_for_a_residual(self):
        module = self.module
        self.assertEqual(module.ATTRIBUTION_STATUS["tier_below_threshold"], "structural_zero")
        self.assertIn(
            "min_score_gap_1.0_bits_HMMER_E1e6",
            module.ATTRIBUTION_BASIS["tier_below_threshold"],
        )
        rows = [self.evidence("A_1", discovery_best_evalue="1.0")]
        attributed, _ = module.attribute_residual(["A_1"], rows, {}, None)
        self.assertEqual(attributed[0]["attribution"], "profile_unassigned_no_discovery_claim")

    def test_evidence_row_without_accession_raises(self):
        with self.assertRaises(ValueError):
            self.module.attribute_residual(["A_1"], [{"matrix_evidence_status": "x"}], {}, None)

    def test_filter_summary_is_optional_for_the_pure_function_and_recorded(self):
        rows = [self.evidence("A_1")]
        _, without = self.module.attribute_residual(["A_1"], rows, {}, None)
        self.assertFalse(without["filter_summary_used"])
        _, with_summary = self.module.attribute_residual(
            ["A_1"], rows, {"per_superfamily": {REGISTRY_TYPE1: {"total": 1}}}, None
        )
        self.assertTrue(with_summary["filter_summary_used"])
        self.assertEqual(with_summary["registry_range_check"]["status"], "ok")


class AssignmentCompletenessTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_completeness_check_reports_ok_when_the_sets_are_equal(self):
        payload = self.module.check_assignment_completeness(
            admissible=["A_1", "A_2"],
            high_confidence=["A_1"],
            hold=["A_2"],
        )
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["admissible_count"], 2)
        self.assertEqual(payload["flow_union_count"], 2)

    def test_completeness_check_reports_the_two_directions_separately(self):
        payload = self.module.check_assignment_completeness(
            admissible=["A_1", "A_2", "A_3"],
            high_confidence=["A_1", "A_4"],
            hold=["A_2"],
        )
        self.assertEqual(payload["status"], "mismatch")
        self.assertEqual(payload["admissible_without_a_flow_row"], ["A_3"])
        self.assertEqual(payload["flow_row_without_an_admissible_assignment"], ["A_4"])


class EndToEndTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def entries(self):
        return [
            {"accession": "GCA_1|A_1", "flow": "hc"},
            {"accession": "GCA_1|A_2", "flow": "hold", "hold_reason": "localization_conflict"},
            # residual, mechanism 1: no profile score and no discovery row
            {"accession": "GCA_1|A_3", **UNASSIGNED},
            # residual, mechanism 2: non-unique discovery claim
            {
                "accession": "GCA_1|A_4",
                "matrix": {},
                "motif": {},
                "classification": {
                    "superfamily_claim": "ambiguous:%s|%s" % (REGISTRY_CYS, REGISTRY_LIPASE),
                    "superfamily_confidence": "ambiguous_2",
                    "n_superfamilies_hit": "2",
                    "best_evalue": "3e-4",
                },
            },
            # residual, mechanism 3: profile score that was not upgraded
            {
                "accession": "GCA_1|A_5",
                "matrix": {
                    "assignment_status": "ambiguous_superfamily",
                    "evidence_status": "profile_score",
                    "phaded_superfamily_best": REGISTRY_TYPE2,
                    "profile_evidence_status": "profile_ambiguous_superfamily",
                    "superfamily_score_gap": "0.9",
                },
                "motif": {},
            },
            # residual, recorded confounder: still attributed by its stage
            {
                "accession": "GCA_1|A_6",
                "matrix": {},
                "motif": {},
                "confounder": {},
            },
            # residual with a complete evidence row
            {"accession": "GCA_1|A_7", **UNASSIGNED},
            # residual with NO row in any frozen evidence table -> unexplained
            {"accession": "GCA_1|A_8"},
        ]

    def fixture(self, **kwargs):
        kwargs.setdefault(
            "pool_external_layers",
            {"strong": [{"protein_id": "GCA_1|A_1", "superfamily": REGISTRY_TYPE1, "score_tier": "strong"}]},
        )
        return build_fixture(self.root, self.entries(), **kwargs)

    def argv(self, paths, output_dir, extra=()):
        argv = [
            "--universe", str(paths["universe"]),
            "--high-confidence", str(paths["high_confidence"]),
            "--hold", str(paths["hold"]),
            "--filter-summary", str(paths["filter_summary"]),
            "--classification", str(paths["classification"]),
            "--confounders", str(paths["confounders"]),
            "--output-dir", str(output_dir),
        ]
        for name in ("matrix", "motif"):
            if paths[name] is not None:
                argv += ["--%s" % name, str(paths[name])]
        for name, path in paths["pool_external_layers"].items():
            argv += ["--pool-external-layer", "%s=%s" % (name, path)]
        argv.extend(extra)
        return argv

    def run_main(self, argv):
        """Run ``main`` and return ``(exit_code, printed_stdout)``."""
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = self.module.main(argv)
        return code, buffer.getvalue()

    def test_missing_input_path_raises(self):
        paths = self.fixture()
        output_dir = self.root / "out"
        for flag, key in (
            ("--universe", "universe"),
            ("--hold", "hold"),
            ("--filter-summary", "filter_summary"),
            ("--classification", "classification"),
        ):
            argv = self.argv(paths, output_dir)
            index = argv.index(flag)
            argv[index + 1] = str(self.root / "does_not_exist.tsv")
            with self.assertRaises(ValueError, msg=flag):
                self.module.main(argv)

    def test_missing_pool_external_layer_path_raises(self):
        paths = self.fixture()
        argv = self.argv(paths, self.root / "out")
        argv += ["--pool-external-layer", "strong=%s" % (self.root / "nope.tsv")]
        with self.assertRaises(ValueError):
            self.module.main(argv)

    def test_output_dir_refusal_for_a_non_empty_directory(self):
        paths = self.fixture()
        output_dir = self.root / "occupied"
        output_dir.mkdir()
        (output_dir / "keep.txt").write_text("do not overwrite me\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            self.run_main(self.argv(paths, output_dir))
        self.assertEqual((output_dir / "keep.txt").read_text(encoding="utf-8"), "do not overwrite me\n")

    def test_output_dir_is_created_when_absent(self):
        paths = self.fixture()
        output_dir = self.root / "fresh" / "results"
        code, _ = self.run_main(self.argv(paths, output_dir))
        self.assertEqual(code, 0)
        self.assertTrue(output_dir.is_dir())

    def test_main_prints_a_json_summary(self):
        paths = self.fixture()
        output_dir = self.root / "results"
        code, printed = self.run_main(self.argv(paths, output_dir))
        self.assertEqual(code, 0)
        payload = json.loads(printed)
        self.assertEqual(payload["residual_count"], 6)
        self.assertEqual(payload["sum_attribution_counts"], 6)
        self.assertTrue(payload["attribution_counts_sum_to_residual"])
        self.assertEqual(payload["unexplained_count"], 1)
        self.assertFalse(payload["check_hold_residual_silenced"])
        self.assertEqual(payload["check_hold_residual_status"], "mismatch")
        self.assertIsNone(payload["pending_buckets"].get("pool_external_boundary"))

    def test_end_to_end_summary_closes_and_is_written(self):
        paths = self.fixture()
        output_dir = self.root / "results"
        code, _ = self.run_main(self.argv(paths, output_dir))
        self.assertEqual(code, 0)
        residual = residual_of(self.module, paths)
        self.assertEqual(len(residual), 6)
        summary_path = output_dir / ("residual_%d_summary.json" % len(residual))
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        self.assertEqual(payload["residual_count"], len(residual))
        self.assertEqual(payload["sum_attribution_counts"], len(residual))
        self.assertTrue(payload["attribution_counts_sum_to_residual"])
        self.assertEqual(
            payload["attribution_counts"]["profile_unassigned_no_discovery_claim"], 3
        )
        self.assertEqual(payload["attribution_counts"]["discovery_claim_not_unique"], 1)
        self.assertEqual(payload["attribution_counts"]["score_tier_not_upgraded"], 1)
        self.assertEqual(payload["attribution_counts"]["unexplained"], 1)
        self.assertEqual(payload["unexplained_count"], 1)
        self.assertEqual(payload["unexplained_accessions"], ["GCA_1|A_8"])
        self.assertFalse(payload["check_hold_residual_silenced"])
        self.assertEqual(payload["check_hold_residual_status"], "mismatch")
        self.assertEqual(payload["confounder_flagged_count"], 1)
        self.assertEqual(payload["pool_external_layers_supplied"], ["strong"])
        self.assertEqual(payload["pool_external_layer_intersection_with_residual"], {"strong": 0})
        self.assertEqual(
            payload["attribution_counts"]["pool_external_boundary"], 0
        )
        self.assertFalse(payload["pending_buckets"])
        self.assertEqual(
            payload["hold_reason_anchor"]["observed_hold_reason_counts_nonexclusive"],
            {"localization_conflict": 1},
        )
        self.assertTrue(payload["hold_reason_anchor"]["matches_frozen_filter_summary"])
        self.assertEqual(payload["assignment_completeness"]["status"], "ok")
        self.assertEqual(payload["registry_range_check"]["status"], "ok")
        self.assertEqual(payload["sanity_anchor"]["hold_reason_union_count"], 1)

    def test_hold_residual_alarm_is_reported_not_silenced(self):
        paths = self.fixture()
        output_dir = self.root / "results"
        self.run_main(self.argv(paths, output_dir))
        check = json.loads((output_dir / "check_hold_residual.json").read_text(encoding="utf-8"))
        self.assertEqual(check["status"], "mismatch")
        self.assertEqual(check["residual_unexplained_count"], check["residual_count"])
        self.assertEqual(check["residual_bucket_names"], [])

    def test_residual_rows_and_attribution_table_are_written(self):
        paths = self.fixture()
        output_dir = self.root / "results"
        self.run_main(self.argv(paths, output_dir))
        residual = residual_of(self.module, paths)
        residual_tsv = output_dir / ("residual_%d.tsv" % len(residual))
        lines = residual_tsv.read_text(encoding="utf-8").strip().split("\n")
        self.assertEqual(len(lines), len(residual) + 1)
        header = lines[0].split("\t")
        for column in ("accession", "layer", "major_superfamily", "score_tier", "attribution"):
            self.assertIn(column, header, column)
        accessions = [line.split("\t")[header.index("accession")] for line in lines[1:]]
        self.assertEqual(len(accessions), len(set(accessions)))
        self.assertEqual(sorted(accessions), residual)

        attribution_tsv = output_dir / ("residual_%d_attribution.tsv" % len(residual))
        rows = attribution_tsv.read_text(encoding="utf-8").strip().split("\n")
        attribution_header = rows[0].split("\t")
        self.assertEqual(attribution_header[:3], ["attribution", "count", "primary_disposition"])
        counted = 0
        listed = set()
        for line in rows[1:]:
            fields = line.split("\t")
            counted += int(fields[attribution_header.index("count")])
            for accession in fields[attribution_header.index("accessions")].split(";"):
                if accession:
                    listed.add(accession)
        self.assertEqual(counted, len(residual))
        self.assertEqual(listed, set(residual))
        self.assertEqual(len(rows) - 1, len(self.module.ATTRIBUTION_ORDER))

    def test_disposition_proposal_is_written_and_covers_every_attribution(self):
        paths = self.fixture()
        output_dir = self.root / "results"
        self.run_main(self.argv(paths, output_dir))
        residual = residual_of(self.module, paths)
        proposal = (output_dir / ("residual_%d_disposition_proposal.md" % len(residual))).read_text(
            encoding="utf-8"
        )
        for name in self.module.ATTRIBUTION_ORDER:
            self.assertIn(name, proposal, name)
        self.assertIn("primary_disposition", proposal)
        self.assertIn("::", proposal)

    def test_input_manifest_carries_real_hashes_for_every_input(self):
        paths = self.fixture()
        output_dir = self.root / "results"
        self.run_main(self.argv(paths, output_dir))
        manifest = (output_dir / "input_manifest.tsv").read_text(encoding="utf-8").strip().split("\n")
        header = manifest[0].split("\t")
        self.assertEqual(header, ["name", "path", "size", "sha256", "status"])
        rows = [dict(zip(header, line.split("\t"))) for line in manifest[1:]]
        names = {row["name"] for row in rows}
        for name in ("universe", "high_confidence", "hold", "filter_summary", "classification",
                     "matrix", "motif", "confounders"):
            self.assertIn(name, names, name)
        for row in rows:
            if row["path"] and Path(row["path"]).is_file():
                self.assertEqual(row["status"], "verified", row["name"])
                self.assertEqual(len(row["sha256"]), 64, row["name"])
            else:
                self.assertEqual(row["sha256"], "pending", row["name"])

    def test_nothing_is_written_outside_the_output_dir(self):
        paths = self.fixture()
        before = {p for p in self.root.rglob("*")}
        output_dir = self.root / "results"
        self.run_main(self.argv(paths, output_dir))
        after = {p for p in self.root.rglob("*")}
        for path in after - before:
            self.assertTrue(
                output_dir in path.parents or path == output_dir,
                "written outside --output-dir: %s" % path,
            )

    def test_residual_is_a_set_difference_not_aggregate_subtraction(self):
        # a high-confidence row outside the universe makes the aggregate
        # arithmetic disagree with the set difference; the script must report the
        # set difference and must surface the mismatch instead of hiding it.
        entries = self.entries()
        entries.append(
            {
                "accession": "GCA_9|OUTSIDE_1",
                "flow": "hc",
                "family": "ePhaZ",
                "in_universe": False,
            }
        )
        paths = build_fixture(self.root, entries)
        output_dir = self.root / "results"
        self.run_main(self.argv(paths, output_dir))
        residual = residual_of(self.module, paths)
        payload = json.loads(
            (output_dir / ("residual_%d_summary.json" % len(residual))).read_text(encoding="utf-8")
        )
        self.assertEqual(payload["residual_count"], len(residual))
        self.assertEqual(payload["residual_count"], 6)
        self.assertEqual(
            payload["aggregate_subtraction_count"],
            payload["universe_count"]
            - payload["high_confidence_count"]
            - payload["hold_count"],
        )
        self.assertNotEqual(payload["aggregate_subtraction_count"], payload["residual_count"])
        self.assertFalse(payload["aggregate_subtraction_agrees_with_the_set_difference"])
        self.assertEqual(payload["check_hold_residual_status"], "mismatch")
        self.assertIn("outside_universe_count", payload["check_hold_residual_problems_source"])

    def test_residual_is_read_through_the_reconcile_module(self):
        module = self.module
        self.assertIs(module.check_hold_residual, module.reconcile.check_hold_residual)
        self.assertIs(module.read_accession_list, module.reconcile.read_accession_list)
        self.assertIs(module.sha256_file, module.reconcile.sha256_file)

    def test_exit_code_guards(self):
        module = self.module
        self.assertEqual(module.exit_code_for(True, "mismatch", 616), 0)
        self.assertEqual(module.exit_code_for(False, "mismatch", 616), 1)
        self.assertEqual(module.exit_code_for(True, "ok", 616), 2)
        self.assertEqual(module.exit_code_for(True, "ok", 0), 0)


class RegistryRangeTests(unittest.TestCase):
    def test_filter_summary_registry_is_checked_against_the_frozen_filter(self):
        module = load_module()
        frozen = set(module.load_frozen_filter_registry())
        payload = module.check_registry_range(frozen | {"not_a_superfamily"}, frozen)
        self.assertEqual(payload["status"], "mismatch")
        self.assertIn("not_a_superfamily", payload["outside_the_frozen_registry"])
        ok = module.check_registry_range(frozen, frozen)
        self.assertEqual(ok["status"], "ok")
        self.assertEqual(ok["outside_the_frozen_registry"], [])


if __name__ == "__main__":
    unittest.main()
