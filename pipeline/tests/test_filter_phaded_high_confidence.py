"""Tests for the PhaDED high-confidence filter's honest-typing label.

Written before the implementation (test-first).  §11.9.4.1 of the project
handoff requires that ``intracellular nPHASCL without lipase box`` (Cys)
candidates be labelled as *inferred* Cys type whose catalytic residue is NOT
per-sequence verified: "no GxSxG + PF06850 + strong Cys-family profile" is a
presumption, not a verified catalytic Cys, and an un-detected Ser type without
GxSxG is not excluded.

Load-bearing invariants:

* every output row carries ``catalytic_residue_verification``;
* Cys-superfamily rows are ``inferred_cys_not_verified``;
* every other superfamily row is ``motif_level`` (residue identified at motif
  level, not biochemically verified per sequence);
* the summary JSON carries an ``honesty_notes.cys_typing`` block spelling out
  the caveat.
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

MOTIF = (
    "accession\tlipase_box_state\tlipase_box_x1\tsbd_pf06850_binding_state\n"
    "ACC_CYS\tnot_detected_pattern\t\tdetected\n"
    "ACC_T1\tsupported\tF\tnot_detected\n"
)

MATRIX = (
    "accession\tgenome\tphaded_superfamily_best\tprofile_evidence_status"
    "\tsignalp_class\tarchitecture_consistency\tinterpro_status\n"
    "ACC_CYS\tG\t" + CYS + "\tprofile_trained_hit\tpending\tconsistent\tinterpro_supported\n"
    "ACC_T1\tG\t" + TYPE1 + "\tprofile_trained_hit\tSP\tconsistent\tinterpro_supported\n"
)

DOMAIN = "accession\tcatalytic_domain_type\nACC_T1\ttype1_verified\n"

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


class HonestTypingLabelTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_hc_filter_"))
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
        self.summary = self.module.filter_candidates(
            self.cls, self.motif, self.domain, self.matrix, self.conf, self.out
        )

    def read_rows(self):
        path = self.out / "high_confidence_candidates.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            return {r["accession"]: r for r in csv.DictReader(handle, delimiter="\t")}

    def test_cys_rows_labelled_inferred_not_verified(self):
        rows = self.read_rows()
        self.assertIn("ACC_CYS", rows)
        self.assertEqual(
            rows["ACC_CYS"]["catalytic_residue_verification"],
            "inferred_cys_not_verified",
        )

    def test_ser_motif_rows_labelled_motif_level(self):
        rows = self.read_rows()
        self.assertIn("ACC_T1", rows)
        self.assertEqual(
            rows["ACC_T1"]["catalytic_residue_verification"],
            "motif_level",
        )

    def test_summary_carries_cys_honesty_note(self):
        data = json.loads(
            (self.out / "filter_summary.json").read_text(encoding="utf-8")
        )
        self.assertIn("honesty_notes", data)
        self.assertIn("inferred_cys_not_verified", data["honesty_notes"]["cys_typing"])


if __name__ == "__main__":
    unittest.main()
