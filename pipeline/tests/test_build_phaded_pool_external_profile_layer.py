"""Tests for the pool-external profile-score layer builder.

Written before the implementation (test-first).  §11.9.4.2 and §11.9.4.3 of
the project handoff require that the pool-external layer:

* be split from the pool-internal high-confidence set and NEVER carry the
  "high_confidence" label (its file name and its per-row layer marker must not
  claim high confidence);
* carry an explicit E-value tier so weak "kept-entirely" specific-family hits
  (E >= 1e-10) are visibly marked "score-only weak" instead of looking like the
  strong hits they sit next to;
* reproduce the published 35,558 pool-external anchors: specific families
  (Cys / dPHAMCL / PhaZ7-like / periplasmic) kept whole, broad families
  (type 1 / type 2 / nPHAMCL) kept only when effective best E < 1e-30, and
  with-lipase excluded (deferred to structural validation).
"""

import importlib.util
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "build_phaded_pool_external_profile_layer.py"
)

CYS = "intracellular nPHASCL without lipase box"
TYPE1 = "extracellular dPHASCL type 1"
TYPE2 = "extracellular dPHASCL type 2"
NPHA = "intracellular nPHAMCL"
DPHA = "extracellular dPHAMCL"
PHAZ7 = "extracellular native-SCL/PhaZ7-like"
PERI = "periplasmic PHA depolymerases"
WITHLIP = "intracellular nPHASCL with lipase box"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "build_phaded_pool_external_profile_layer", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def row(sf, disc_e, trained_e="", override="0"):
    return {
        "protein_id": "P",
        "genome": "G",
        "superfamily": sf,
        "discovery_best_family": "DED_hfam_x",
        "discovery_families_hit": "DED_hfam_x",
        "discovery_best_evalue": disc_e,
        "trained_best_model": "",
        "trained_best_evalue": trained_e,
        "override_by_trained": override,
        "model_layer": "discovery_hmm_uncalibrated",
    }


class ScoreTierTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_boundaries(self):
        self.assertEqual(self.module.score_tier(1e-40), "strong")
        self.assertEqual(self.module.score_tier(1e-30), "mid")
        self.assertEqual(self.module.score_tier(1e-15), "mid")
        self.assertEqual(self.module.score_tier(1e-10), "weak")
        self.assertEqual(self.module.score_tier(1e-5), "weak")
        self.assertEqual(self.module.score_tier(1.0), "weak")


class SelectAndTierTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def select(self, *rows):
        return self.module.select_and_tier(list(rows))

    def test_specific_family_kept_whole_even_weak(self):
        kept = self.select(row(CYS, "1e-5"))
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]["score_tier"], "weak")
        self.assertEqual(kept[0]["evidence_layer"], "score_only")

    def test_broad_family_weak_excluded(self):
        self.assertEqual(self.select(row(TYPE1, "1e-5")), [])
        self.assertEqual(self.select(row(TYPE2, "1e-9")), [])
        self.assertEqual(self.select(row(NPHA, "1e-10")), [])

    def test_broad_family_strong_kept(self):
        kept = self.select(row(TYPE1, "1e-40"))
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]["score_tier"], "strong")

    def test_trained_override_sets_effective_evalue(self):
        kept = self.select(row(TYPE1, "1e-5", trained_e="1e-40", override="1"))
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]["score_tier"], "strong")
        self.assertEqual(kept[0]["best_evalue"], "1e-40")

    def test_with_lipase_always_excluded(self):
        self.assertEqual(self.select(row(WITHLIP, "1e-60")), [])

    def test_output_rows_carry_pending_evidence_markers(self):
        kept = self.select(row(DPHA, "1e-40"))
        self.assertEqual(kept[0]["architecture_status"], "pending")
        self.assertEqual(kept[0]["localization_status"], "pending")
        self.assertEqual(kept[0]["interpro_status"], "pending")

    def test_filename_never_claims_high_confidence(self):
        self.assertNotIn("high_confidence", self.module.OUTPUT_NAME)
        self.assertNotIn("high_confidence", self.module.EVIDENCE_LAYER)


if __name__ == "__main__":
    unittest.main()
