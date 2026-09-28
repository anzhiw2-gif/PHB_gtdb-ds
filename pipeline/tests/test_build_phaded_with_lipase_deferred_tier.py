"""Tests for the with-lipase deferred structural-validation tier builder.

Written before the implementation (test-first, 2026-09-19).  The pool-external
high-confidence filter excluded the ``intracellular nPHASCL with lipase box``
superfamily (1,206,655 pool-external hits) from the 35,558 high-confidence set.
Those excluded hits are NOT discarded: they are preserved as a deferred tier
for later structural validation (Foldseek vs PDB 8YNV / predicted-structure
screening).  This module builds that tier and its provenance.

Load-bearing invariants:

* the trained-profile priority rule is exactly "trained best E-value strictly
  smaller than the discovery best E-value" (an equal or weaker trained hit does
  NOT override the discovery assignment) — this reproduces the published
  exclusion count 1,238,812 - 32,157 = 1,206,655;
* the deferred tier is labelled recall-only / ``deferred_structural_validation``
  and never feeds the high-confidence count;
* the pool-internal with-lipase rows are extracted by intersecting the
  discovery classification with the frozen candidate universe and keep the
  classification's own columns.
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "build_phaded_with_lipase_deferred_tier.py"
)

WITH_LIPASE = "intracellular nPHASCL with lipase box"

MANIFEST = (
    "profile_id\tprofile_kind\tphaded_superfamily\tphaded_family_id\tmodel_status\n"
    "family_DED_hfam_2_6cf580a7126e\tfamily\t" + WITH_LIPASE + "\tDED_hfam_2\treference_only\n"
    "family_DED_hfam_4_3c28e8cee1c4\tfamily\tintracellular nPHAMCL\tDED_hfam_4\ttrained\n"
    "family_DED_hfam_52_22ef55bacdb2\tfamily\textracellular dPHASCL type 1\tDED_hfam_52\ttrained\n"
    "family_DED_hfam_7_7fa5e8040d92\tfamily\textracellular native-SCL/PhaZ7-like\tDED_hfam_7\treference_only\n"
    "superfamily_extracellular_dPHASCL_type_1_b97713ebe577\tsuperfamily\textracellular dPHASCL type 1\t\ttrained\n"
)

#: F3 fixture: the SAME superfamily assignments as MANIFEST, expressed with the
#: post-redesign ``model_layer`` column.  ``family_DED_hfam_2_6cf580a7126e`` is the
#: trap: it is a DISCOVERY-layer HMM (recall-only, uncalibrated) whose legacy
#: ``model_status`` still reads ``trained``, so a ``model_status``-keyed overlay
#: would let it move a pool-external protein out of the deferred tier.
MANIFEST_LAYERED = (
    "profile_id\tprofile_kind\tphaded_superfamily\tphaded_family_id"
    "\tmodel_layer\tmodel_status\tfunctional_calibration_status\n"
    "family_DED_hfam_2_6cf580a7126e\tfamily\t" + WITH_LIPASE
    + "\tDED_hfam_2\tdiscovery_hmm_uncalibrated\ttrained\tnot_function_calibrated\n"
    "family_DED_hfam_4_3c28e8cee1c4\tfamily\tintracellular nPHAMCL"
    "\tDED_hfam_4\tsequence_family_hmm_validated\ttrained\tcandidate_gate_passed_not_promoted\n"
    "family_DED_hfam_52_22ef55bacdb2\tfamily\textracellular dPHASCL type 1"
    "\tDED_hfam_52\tsequence_family_hmm_validated\ttrained\tnot_function_calibrated\n"
    "family_DED_hfam_7_7fa5e8040d92\tfamily\textracellular native-SCL/PhaZ7-like"
    "\tDED_hfam_7\treference_query_only\treference_only\tnot_function_calibrated\n"
    "superfamily_extracellular_dPHASCL_type_1_b97713ebe577\tsuperfamily"
    "\textracellular dPHASCL type 1\t\tsequence_family_hmm_validated\ttrained"
    "\tnot_function_calibrated\n"
)

DISCOVERY = (
    "protein_id\tfamilies_hit\tbest_family\tbest_evalue\n"
    "GCA_1|CTX_1_1\tDED_hfam_2\tDED_hfam_2\t4.9e-14\n"
    "GCA_1|CTX_1_2\tDED_hfam_2\tDED_hfam_2\t5e-07\n"
    "GCA_1|CTX_1_3\tDED_hfam_2\tDED_hfam_2\t1e-80\n"
    "GCA_2|CTX_2_1\tDED_hfam_2\tDED_hfam_2\t2e-10\n"
    "GCA_3|CTX_3_1\tDED_hfam_7\tDED_hfam_7\t1e-100\n"
)

TRAINED = (
    "protein_id\tbest_model\tbest_evalue\n"
    # trained stronger than discovery -> override, leaves with-lipase
    "GCA_1|CTX_1_1\tfamily_DED_hfam_4_3c28e8cee1c4\t1e-40\n"
    # trained weaker than discovery -> no override, stays with-lipase
    "GCA_1|CTX_1_2\tfamily_DED_hfam_4_3c28e8cee1c4\t0.1\n"
    # trained equal to discovery -> strict less fails, stays with-lipase
    "GCA_1|CTX_1_3\tfamily_DED_hfam_4_3c28e8cee1c4\t1e-80\n"
    # trained-only protein (no discovery row)
    "GCA_4|CTX_4_1\tsuperfamily_extracellular_dPHASCL_type_1_b97713ebe577\t3e-50\n"
)

#: F3 fixture: the scan's best "trained" model for every row is the
#: DISCOVERY-layer profile ``family_DED_hfam_2_6cf580a7126e``.  Under the legacy
#: ``model_status == "trained"`` rule that recall-only model would take over the
#: assignment; under ``model_layer`` it cannot move anything.
TRAINED_DISCOVERY_LAYER = (
    "protein_id\tbest_model\tbest_evalue\n"
    "GCA_1|CTX_1_1\tfamily_DED_hfam_2_6cf580a7126e\t1e-90\n"
    "GCA_5|CTX_5_1\tfamily_DED_hfam_2_6cf580a7126e\t1e-90\n"
)

CLASSIFICATION = (
    "protein_id\tsuperfamily_claim\tsuperfamily_confidence\tbest_evalue\n"
    "GCA_9|CTX_9_1\t" + WITH_LIPASE + "\tunique\t2e-20\n"
    "GCA_9|CTX_9_2\tambiguous:" + WITH_LIPASE + "|intracellular nPHAMCL\tambiguous_2\t2e-20\n"
    "GCA_9|CTX_9_3\tintracellular nPHASCL without lipase box\tunique\t2e-20\n"
)

CANDIDATE_FAA = (
    ">GCA_9|CTX_9_1 description\nMKWV\n"
    ">GCA_9|CTX_9_2 description\nMKWV\n"
    ">GCA_9|CTX_9_3 description\nMKWV\n"
)


MATRIX = (
    "accession\tgenome\tphaded_superfamily_best\tprofile_evidence_status"
    "\tsignalp_class\tarchitecture_consistency\tinterpro_status\n"
    # trained-matrix branch fires with a different superfamily -> not deferred
    "GCA_9|CTX_9_1\tGCA_9\textracellular dPHASCL type 1\tprofile_trained_hit"
    "\tSP\tconsistent\tinterpro_supported\n"
    # matrix branch itself claims with-lipase (ambiguous family) -> deferred
    "GCA_9|CTX_9_2\tGCA_9\tintracellular nPHASCL with lipase box"
    "\tprofile_ambiguous_family\tOTHER\tconsistent\tinterpro_supported\n"
    # no matrix row -> discovery-unique fallback -> deferred
    "GCA_9|CTX_9_4\tGCA_9\tintracellular nPHASCL with lipase box"
    "\tprofile_trained_hit\tOTHER\tconsistent\tinterpro_supported\n"
)

CONFOUNDERS = (
    "accession\tgenome\n"
    "GCA_9|CTX_9_4\tGCA_9\n"
)

CLASSIFICATION2 = (
    "protein_id\tsuperfamily_claim\tsuperfamily_confidence\tbest_evalue\n"
    "GCA_9|CTX_9_1\t" + WITH_LIPASE + "\tunique\t2e-20\n"
    "GCA_9|CTX_9_2\t" + WITH_LIPASE + "\tunique\t2e-20\n"
    "GCA_9|CTX_9_3\t" + WITH_LIPASE + "\tunique\t2e-20\n"
    "GCA_9|CTX_9_4\t" + WITH_LIPASE + "\tunique\t2e-20\n"
)

CANDIDATE_FAA2 = (
    ">GCA_9|CTX_9_1 description\nMKWV\n"
    ">GCA_9|CTX_9_2 description\nMKWV\n"
    ">GCA_9|CTX_9_3 description\nMKWV\n"
    ">GCA_9|CTX_9_4 description\nMKWV\n"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "build_phaded_with_lipase_deferred_tier", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


class ManifestMappingTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_wl_deferred_"))
        self.manifest = self.tmp / "manifest.tsv"
        write_text(self.manifest, MANIFEST_LAYERED)

    def test_family_superfamily_mapping(self):
        fam2sf = self.module.read_family_superfamily(self.manifest)
        self.assertEqual(fam2sf["DED_hfam_2"], WITH_LIPASE)
        self.assertEqual(fam2sf["DED_hfam_7"], "extracellular native-SCL/PhaZ7-like")
        # superfamily-only rows have an empty family id and must not appear
        self.assertNotIn("", fam2sf)

    def test_trained_model_mapping(self):
        model2sf = self.module.read_trained_model_superfamily(self.manifest)
        self.assertEqual(model2sf["family_DED_hfam_4_3c28e8cee1c4"], "intracellular nPHAMCL")
        self.assertEqual(
            model2sf["superfamily_extracellular_dPHASCL_type_1_b97713ebe577"],
            "extracellular dPHASCL type 1",
        )
        # the discovery-layer profile is not part of the validated overlay, even
        # though its legacy model_status column still reads "trained"
        self.assertNotIn("family_DED_hfam_2_6cf580a7126e", model2sf)
        # a reference_query_only profile (no HMM) is not part of the overlay
        self.assertNotIn("family_DED_hfam_7_7fa5e8040d92", model2sf)

    def test_legacy_manifest_yields_no_overlay_models(self):
        """Absent model_layer => non-discriminating fallback, never "trained"."""
        legacy = self.tmp / "legacy.tsv"
        write_text(legacy, MANIFEST)
        self.assertEqual(self.module.read_trained_model_superfamily(legacy), {})


class ClassificationRuleTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_wl_deferred_"))
        self.discovery = self.tmp / "discovery.tsv"
        self.trained = self.tmp / "trained.tsv"
        write_text(self.discovery, DISCOVERY)
        write_text(self.trained, TRAINED)

    def classify(self, manifest_text=MANIFEST_LAYERED):
        manifest = self.tmp / "manifest.tsv"
        write_text(manifest, manifest_text)
        fam2sf = self.module.read_family_superfamily(manifest)
        model2sf = self.module.read_trained_model_superfamily(manifest)
        trained = self.module.load_trained_best(self.trained)
        rows, counts = self.module.classify_pool_external(
            self.discovery, trained, fam2sf, model2sf
        )
        return rows, counts

    def test_override_rule_is_strictly_stronger_trained(self):
        rows, _ = self.classify()
        by_id = {r["protein_id"]: r for r in rows}
        # trained 1e-40 < discovery 4.9e-14 -> trained superfamily wins
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["superfamily"], "intracellular nPHAMCL")
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["override_by_trained"], "1")
        # trained 0.1 > discovery 5e-07 -> discovery wins
        self.assertEqual(by_id["GCA_1|CTX_1_2"]["superfamily"], WITH_LIPASE)
        self.assertEqual(by_id["GCA_1|CTX_1_2"]["override_by_trained"], "0")
        # trained 1e-80 == discovery 1e-80 -> strict less fails, discovery wins
        self.assertEqual(by_id["GCA_1|CTX_1_3"]["superfamily"], WITH_LIPASE)
        self.assertEqual(by_id["GCA_1|CTX_1_3"]["override_by_trained"], "0")

    def test_trained_only_protein_keeps_trained_superfamily(self):
        rows, _ = self.classify()
        by_id = {r["protein_id"]: r for r in rows}
        self.assertEqual(
            by_id["GCA_4|CTX_4_1"]["superfamily"], "extracellular dPHASCL type 1"
        )
        self.assertEqual(by_id["GCA_4|CTX_4_1"]["override_by_trained"], "1")

    def test_deferred_tier_contains_exactly_with_lipase_rows(self):
        rows, counts = self.classify()
        deferred = [r for r in rows if r["superfamily"] == WITH_LIPASE]
        self.assertEqual(
            {r["protein_id"] for r in deferred},
            {"GCA_1|CTX_1_2", "GCA_1|CTX_1_3", "GCA_2|CTX_2_1"},
        )
        self.assertEqual(counts.get(WITH_LIPASE, 0), 3)
        # the overridden row is recorded in the counts under its new superfamily
        self.assertEqual(counts["intracellular nPHAMCL"], 1)

    def test_genome_split_from_protein_id(self):
        self.assertEqual(self.module.genome_of("GCA_1|CTX_1_1"), "GCA_1")


class PoolInternalExtractionTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_wl_deferred_"))
        self.cls = self.tmp / "classification.tsv"
        self.faa = self.tmp / "candidates.faa"
        write_text(self.cls, CLASSIFICATION)
        write_text(self.faa, CANDIDATE_FAA)

    def test_pool_internal_extraction_unique_only(self):
        rows = self.module.extract_pool_internal_with_lipase(self.cls, self.faa)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["protein_id"], "GCA_9|CTX_9_1")
        # ambiguous and other-superfamily rows are excluded
        ids = {r["protein_id"] for r in rows}
        self.assertNotIn("GCA_9|CTX_9_2", ids)
        self.assertNotIn("GCA_9|CTX_9_3", ids)


class PoolInternalFilterConsistentTests(unittest.TestCase):
    """The 13-row filter-consistent pool-internal set: same assignment
    precedence as filter_phaded_high_confidence.py (matrix branch first, then
    discovery-unique fallback), no common/architecture screening."""

    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_wl_deferred_"))
        self.cls = self.tmp / "classification.tsv"
        self.matrix = self.tmp / "matrix.tsv"
        self.conf = self.tmp / "confounders.tsv"
        write_text(self.cls, CLASSIFICATION2)
        write_text(self.matrix, MATRIX)
        write_text(self.conf, CONFOUNDERS)

    def test_matrix_branch_precedence(self):
        rows = self.module.extract_pool_internal_filter_consistent(
            self.cls, self.matrix, self.conf
        )
        ids = {r["protein_id"] for r in rows}
        # matrix trained hit for a different superfamily wins -> excluded
        self.assertNotIn("GCA_9|CTX_9_1", ids)
        # matrix ambiguous-family with-lipase -> included
        self.assertIn("GCA_9|CTX_9_2", ids)
        # no matrix row, discovery-unique with-lipase -> included
        self.assertIn("GCA_9|CTX_9_3", ids)
        # matrix branch claiming with-lipase -> included; confounder membership
        # is a screening flag only and never removes a deferred row
        self.assertIn("GCA_9|CTX_9_4", ids)
        self.assertEqual(len(rows), 3)


class ModelLayerTests(unittest.TestCase):
    """F3: the trained overlay is keyed on ``model_layer``, never ``model_status``.

    ``model_status`` still reads ``trained`` for discovery-layer HMMs, so keying
    the overlay on it would let a recall-only model reassign a pool-external
    protein out of the deferred structural-validation tier.
    """

    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_wl_layer_"))
        self.manifest = self.tmp / "manifest.tsv"
        write_text(self.manifest, MANIFEST_LAYERED)

    def test_only_validated_layers_form_the_trained_overlay(self):
        model2sf = self.module.read_trained_model_superfamily(self.manifest)
        self.assertEqual(set(model2sf),
                         {"family_DED_hfam_4_3c28e8cee1c4",
                          "family_DED_hfam_52_22ef55bacdb2",
                          "superfamily_extracellular_dPHASCL_type_1_b97713ebe577"})
        # family_DED_hfam_2 is the discovery-layer trap: legacy word "trained",
        # layer discovery_hmm_uncalibrated -> never an overlay model.
        self.assertNotIn("family_DED_hfam_2_6cf580a7126e", model2sf)
        # no-HMM reference profile -> never an overlay model
        self.assertNotIn("family_DED_hfam_7_7fa5e8040d92", model2sf)

    def test_unknown_declared_layer_fails_closed(self):
        broken = self.tmp / "broken.tsv"
        write_text(broken, MANIFEST_LAYERED.replace(
            "discovery_hmm_uncalibrated", "discovery_hmm"))
        with self.assertRaisesRegex(ValueError, "model_layer"):
            self.module.read_trained_model_superfamily(broken)

    def classify(self, trained_text, manifest=None):
        """Run the full manifest read + classification path like ``main`` does."""
        manifest = manifest or self.manifest
        self.trained = self.tmp / "trained.tsv"
        self.discovery = self.tmp / "discovery.tsv"
        write_text(self.trained, trained_text)
        write_text(self.discovery, DISCOVERY)
        report = self.module.read_manifest_summary(manifest)
        return self.module.classify_pool_external(
            self.discovery,
            self.module.load_trained_best(self.trained),
            self.module.read_family_superfamily(manifest),
            self.module.read_trained_model_superfamily(manifest),
            known_profile_ids=report["profile_ids"],
            profile2sf=self.module.read_profile_id_superfamily(manifest),
        )

    def test_discovery_layer_hmm_cannot_move_a_protein_out_of_the_deferred_tier(self):
        rows, counts = self.classify(TRAINED_DISCOVERY_LAYER)
        by_id = {r["protein_id"]: r for r in rows}
        # GCA_1|CTX_1_1 is a with-lipase discovery hit whose strongest "trained"
        # hit is the DISCOVERY-layer DED_hfam_2 model (1e-90 << 4.9e-14).  Keyed
        # on model_status="trained" that would override the assignment and
        # silently delete the row from the deferred tier; keyed on model_layer it
        # cannot move anything.
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["superfamily"], WITH_LIPASE)
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["override_by_trained"], "0")
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["trained_best_model"], "")
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["discovery_layer_overlay_suppressed"],
                         "family_DED_hfam_2_6cf580a7126e")
        # a protein seen ONLY by the discovery-layer model keeps that model's
        # family superfamily (recall may add, never delete) and is never
        # recorded as a validated overlay
        self.assertEqual(by_id["GCA_5|CTX_5_1"]["superfamily"], WITH_LIPASE)
        self.assertEqual(by_id["GCA_5|CTX_5_1"]["override_by_trained"], "0")
        self.assertEqual(by_id["GCA_5|CTX_5_1"]["model_layer"], "discovery_hmm_uncalibrated")
        # no with-lipase row is deleted: GCA_1|CTX_1_1 stays in the tier
        self.assertEqual(
            {r["protein_id"] for r in rows if r["superfamily"] == WITH_LIPASE},
            {"GCA_1|CTX_1_1", "GCA_1|CTX_1_2", "GCA_1|CTX_1_3", "GCA_2|CTX_2_1",
             "GCA_5|CTX_5_1"},
        )
        self.assertEqual(counts.get(WITH_LIPASE, 0), 5)

    def test_validated_overlay_still_assigns_where_it_is_validated(self):
        rows, _ = self.classify(TRAINED)
        by_id = {r["protein_id"]: r for r in rows}
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["superfamily"], "intracellular nPHAMCL")
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["override_by_trained"], "1")
        self.assertEqual(by_id["GCA_4|CTX_4_1"]["superfamily"],
                         "extracellular dPHASCL type 1")

    def test_legacy_manifest_falls_back_to_non_discriminating(self):
        """A manifest without ``model_layer`` cannot express a validated layer."""
        legacy = self.tmp / "legacy.tsv"
        write_text(legacy, MANIFEST)
        # fail-closed: nothing silently becomes an overlay model...
        self.assertEqual(self.module.read_trained_model_superfamily(legacy), {})
        # ...and the fallback is reported rather than assumed
        report = self.module.read_model_layer_report(legacy)
        self.assertFalse(report["model_layer_column_present"])
        self.assertEqual(report["layer_source"], "absent_model_layer_column")
        self.assertEqual(report["fallback_layer"], "reference_query_only")
        self.assertIn("never assumed validated", report["fallback_note"])
        # a legacy "trained" row can never override the discovery assignment
        rows, counts = self.classify(TRAINED, manifest=legacy)
        by_id = {r["protein_id"]: r for r in rows}
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["superfamily"], WITH_LIPASE)
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["override_by_trained"], "0")
        self.assertEqual(by_id["GCA_1|CTX_1_1"]["discovery_layer_overlay_suppressed"],
                         "family_DED_hfam_4_3c28e8cee1c4")
        # and the candidate survives: it is labelled, never dropped
        self.assertIn("GCA_1|CTX_1_1", by_id)
        # TRAINED has no GCA_5 row, so the legacy tier holds 4 with-lipase rows
        self.assertEqual(counts.get(WITH_LIPASE, 0), 4)
        # the pool-external protein seen only by a legacy "trained" profile keeps
        # its declared superfamily instead of being dropped
        self.assertEqual(by_id["GCA_4|CTX_4_1"]["superfamily"],
                         "extracellular dPHASCL type 1")
        self.assertEqual(by_id["GCA_4|CTX_4_1"]["override_by_trained"], "0")

    def test_declared_layer_column_is_reported_as_such(self):
        report = self.module.read_model_layer_report(self.manifest)
        self.assertTrue(report["model_layer_column_present"])
        self.assertEqual(report["layer_source"], "model_layer_column")
        self.assertEqual(report["fallback_note"], "")
        self.assertEqual(report["model_layer_counts"],
                         {"discovery_hmm_uncalibrated": 1,
                          "sequence_family_hmm_validated": 3,
                          "reference_query_only": 1})

    def test_overlay_never_promotes_the_deferred_layer(self):
        """No layer may relabel the deferred tier as validated or calibrated."""
        self.assertEqual(self.module.MODEL_LAYER, "discovery_hmm_uncalibrated")
        self.assertNotIn(self.module.MODEL_LAYER,
                         self.module.TRAINED_OVERLAY_MODEL_LAYERS)
        self.assertEqual(self.module.DEFERRED_REASON, "deferred_structural_validation")
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn('"trained_profile"', source)


class ModelLayerCliTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_wl_cli_"))
        self.discovery = self.tmp / "discovery.tsv"
        self.trained = self.tmp / "trained.tsv"
        self.manifest = self.tmp / "manifest.tsv"
        self.cls = self.tmp / "classification.tsv"
        self.faa = self.tmp / "candidates.faa"
        write_text(self.discovery, DISCOVERY)
        write_text(self.trained, TRAINED)
        write_text(self.manifest, MANIFEST_LAYERED)
        write_text(self.cls, CLASSIFICATION)
        write_text(self.faa, CANDIDATE_FAA)

    def test_main_reports_layer_and_keeps_the_deferred_tier_intact(self):
        import json
        out_dir = self.tmp / "out"
        code = self.module.main([
            "--discovery", str(self.discovery), "--trained", str(self.trained),
            "--profile-manifest", str(self.manifest),
            "--pool-classification", str(self.cls), "--candidate-faa", str(self.faa),
            "--out-dir", str(out_dir),
        ])
        self.assertEqual(code, 0)
        summary = json.loads((out_dir / "deferred_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["deferred_reason"], "deferred_structural_validation")
        self.assertEqual(summary["model_layer"], "discovery_hmm_uncalibrated")
        self.assertEqual(summary["trained_overlay_model_layers"],
                         ["calibrated_candidate_model", "sequence_family_hmm_validated"])
        self.assertEqual(summary["decision_field"], "model_layer")
        self.assertFalse(summary["legacy_model_status_is_a_decision_input"])
        self.assertEqual(summary["with_lipase_pool_external_deferred"], 3)
        self.assertNotIn("high_confidence", summary)
        self.assertNotIn("candidate_count", summary)
        # the deferred rows are preserved, marked and never deleted
        deferred_path = out_dir / "pool_external_with_lipase_deferred.tsv"
        self.assertTrue(deferred_path.exists())
        with deferred_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 3)
        self.assertEqual({row["model_layer"] for row in rows},
                         {"discovery_hmm_uncalibrated"})
        self.assertNotIn("GCA_1|CTX_1_1", {row["protein_id"] for row in rows})
        self.assertIn("GCA_1|CTX_1_2", {row["protein_id"] for row in rows})


if __name__ == "__main__":
    unittest.main()
