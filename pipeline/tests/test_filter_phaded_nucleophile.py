"""Tests for nucleophile-type unification and explicit hold marking.

Written before the implementation (test-first).  The pool-in families have
asymmetric nucleophile evidence: Ser families carry a GxSxG lipase box (direct
Ser-nucleophile marker), but the Cys superfamily is typed only indirectly
("no GxSxG" + PF06850).  This adds:

* a ``nucleophile_type`` column (ser / cys / ahsmg / conflict / undetermined);
* a ``hold_candidates.tsv`` table listing every excluded row with an explicit
  ``hold_reason`` so contradiction rows (Cys-with-GxSxG, Ser-without-GxSxG,
  localization mismatch) are never silently dropped.
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1] / "scripts" / "filter_phaded_high_confidence.py"
)

CYS = "intracellular nPHASCL without lipase box"
TYPE1 = "extracellular dPHASCL type 1"
NPHA = "intracellular nPHAMCL"
PHAZ7 = "extracellular native-SCL/PhaZ7-like"

MOTIF = (
    "accession\tlipase_box_state\tlipase_box_x1\tsbd_pf06850_binding_state\tahsmg_state\n"
    "CYS_PASS\tnot_detected_pattern\t\tdetected\tnot_detected_pattern\n"
    "T1_PASS\tsupported\tF\tnot_detected\tnot_detected_pattern\n"
    "CYS_CONFLICT\tsupported\tF\tdetected\tnot_detected_pattern\n"
    "NPHA_CONFLICT\tnot_detected_pattern\t\tnot_detected\tnot_detected_pattern\n"
)

MATRIX = (
    "accession\tgenome\tphaded_superfamily_best\tprofile_evidence_status"
    "\tsignalp_class\tarchitecture_consistency\tinterpro_status\n"
    "CYS_PASS\tG\t" + CYS + "\tprofile_trained_hit\tpending\tconsistent\tinterpro_supported\n"
    "T1_PASS\tG\t" + TYPE1 + "\tprofile_trained_hit\tSP\tconsistent\tinterpro_supported\n"
    "CYS_CONFLICT\tG\t" + CYS + "\tprofile_trained_hit\tpending\tconsistent\tinterpro_supported\n"
    "NPHA_CONFLICT\tG\t" + NPHA + "\tprofile_trained_hit\tpending\tconsistent\tinterpro_supported\n"
)

DOMAIN = "accession\tcatalytic_domain_type\nT1_PASS\ttype1_verified\n"

CLASSIFICATION = "protein_id\tsuperfamily_claim\tsuperfamily_confidence\tbest_evalue\n"

CONFOUNDERS = "accession\tgenome\n"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "filter_phaded_high_confidence", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


class NucleophileTypeHelperTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()

    def test_ser_cys_ahsmg_conflict(self):
        f = self.module.nucleophile_type
        self.assertEqual(f(CYS, "not_detected_pattern", ""), "cys")
        self.assertEqual(f(CYS, "supported", ""), "conflict")
        self.assertEqual(f(TYPE1, "supported", ""), "ser")
        self.assertEqual(f(TYPE1, "not_detected_pattern", ""), "conflict")
        self.assertEqual(f(NPHA, "not_detected_pattern", ""), "conflict")
        self.assertEqual(f(PHAZ7, "", "supported"), "ahsmg")
        self.assertEqual(f(PHAZ7, "", "not_detected_pattern"), "conflict")


class NucleophileFilterTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_hc_nuc_"))
        self.motif = self.tmp / "motif.tsv"
        self.matrix = self.tmp / "matrix.tsv"
        self.domain = self.tmp / "domain.tsv"
        self.cls = self.tmp / "cls.tsv"
        self.conf = self.tmp / "conf.tsv"
        write_text(self.motif, MOTIF)
        write_text(self.matrix, MATRIX)
        write_text(self.domain, DOMAIN)
        write_text(self.cls, CLASSIFICATION)
        write_text(self.conf, CONFOUNDERS)
        self.out = self.tmp / "out"
        self.module.filter_candidates(
            self.cls, self.motif, self.domain, self.matrix, self.conf, self.out
        )

    def read_rows(self, name):
        with (self.out / name).open(encoding="utf-8", newline="") as h:
            return {r["accession"]: r for r in csv.DictReader(h, delimiter="\t")}

    def test_pass_rows_carry_nucleophile_type(self):
        rows = self.read_rows("high_confidence_candidates.tsv")
        self.assertEqual(rows["CYS_PASS"]["nucleophile_type"], "cys")
        self.assertEqual(rows["T1_PASS"]["nucleophile_type"], "ser")

    def test_conflict_rows_are_held_not_passed(self):
        passed = self.read_rows("high_confidence_candidates.tsv")
        self.assertNotIn("CYS_CONFLICT", passed)
        self.assertNotIn("NPHA_CONFLICT", passed)
        held = self.read_rows("hold_candidates.tsv")
        self.assertEqual(held["CYS_CONFLICT"]["hold_reason"], "nucleophile_conflict")
        self.assertEqual(held["NPHA_CONFLICT"]["hold_reason"], "nucleophile_conflict")

    def test_summary_records_hold_breakdown(self):
        data = json.loads((self.out / "filter_summary.json").read_text(encoding="utf-8"))
        self.assertIn("hold_breakdown", data)
        self.assertEqual(data["hold_breakdown"]["nucleophile_conflict"], 2)


if __name__ == "__main__":
    unittest.main()
