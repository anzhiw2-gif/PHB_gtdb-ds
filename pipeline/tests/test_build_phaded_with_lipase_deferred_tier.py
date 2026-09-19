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

import importlib.util
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
        write_text(self.manifest, MANIFEST)

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
        # reference-only profiles are not part of the trained overlay
        self.assertNotIn("family_DED_hfam_2_6cf580a7126e", model2sf)


class ClassificationRuleTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_wl_deferred_"))
        self.discovery = self.tmp / "discovery.tsv"
        self.trained = self.tmp / "trained.tsv"
        write_text(self.discovery, DISCOVERY)
        write_text(self.trained, TRAINED)

    def classify(self):
        manifest = self.tmp / "manifest.tsv"
        write_text(manifest, MANIFEST)
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


if __name__ == "__main__":
    unittest.main()
