"""Tests for the Cys nucleophile motif evidence annotator.

Written before the implementation (test-first, 2026-09-20).

Scientific basis
----------------
The `intracellular nPHASCL without lipase box` (Cys) superfamily has, until now,
carried only an *inferred* catalytic residue (``inferred_cys_not_verified``): the
Cys typing rests on "no GxSxG + PF06850 + strong Cys-family profile", never on a
catalytic-residue observation.

Two literature rules are available:

* Knoll 2009 states that the Cys-type nucleophile has a **hydrophobic residue at
  position -1** (``cysteine-1``);
* the only experimentally verified Cys-type nucleophile is **Ralstonia eutropha
  PhaZ1 Cys183** (site-directed mutagenesis: C183A and D355A/H388Q abolish
  activity while C183S retains it; PMID 16233560), whose local context is
  ``V182-C183-Q184``.

This module reports, per candidate, whether the sequence carries a cysteine whose
``-1/+1`` context matches that single literature anchor.  It is a **sequence
context match**, NOT a catalytic-activity verification: every output row must say
so explicitly, and the calibration is fail-closed (if the pattern does not
reproduce the anchor, nothing is written).
"""

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "annotate_phaded_cys_nucleophile_motif.py"
)

# anchor fixture: an 8 aa sequence with V-C-Q at positions 3-5
ANCHOR_SEQ = "MKVCQAAK"
REF_FASTA = (
    ">DED_hfam_65_0001|CAJ92291.1\n" + ANCHOR_SEQ + "\n"
    ">DED_hfam_4_0002|AAM63408.1\nMKLCDAAK\n"
)
CAND_FASTA = (
    ">cand_exact\nMKVCQAAK\n"
    ">cand_loose\nMKLCQAAK\n"
    ">cand_minus1\nMKVCDAAK\n"
    ">cand_none\nMKDCAAAK\n"
    ">cand_nocys\nMKDAAAAK\n"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "annotate_phaded_cys_nucleophile_motif", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


class CysContextTests(unittest.TestCase):
    def setUp(self):
        self.m = load_module()

    def test_contexts_are_one_based_with_flanks(self):
        ctxs = self.m.cys_contexts("MKVCQAAK")
        self.assertEqual(len(ctxs), 1)
        c = ctxs[0]
        self.assertEqual(c["position"], 4)
        self.assertEqual(c["minus1"], "V")
        self.assertEqual(c["plus1"], "Q")
        self.assertEqual(c["tripeptide"], "VCQ")

    def test_terminal_cysteines_have_empty_flanks(self):
        ctxs = self.m.cys_contexts("CA")
        self.assertEqual(ctxs[0]["minus1"], "")
        self.assertEqual(ctxs[0]["plus1"], "A")
        ctxs = self.m.cys_contexts("AC")
        self.assertEqual(ctxs[0]["minus1"], "A")
        self.assertEqual(ctxs[0]["plus1"], "")


class CysMotifTierTests(unittest.TestCase):
    def setUp(self):
        self.m = load_module()

    def tier(self, seq):
        return self.m.classify_motif(seq)["tier"]

    def test_exact_anchor_context(self):
        self.assertEqual(self.tier("MKVCQAAK"), "anchor_exact_vcq")

    def test_loose_hydrophobic_c_q(self):
        # L is hydrophobic but not the anchor's V
        self.assertEqual(self.tier("MKLCQAAK"), "anchor_context_hydrophobic_c_q")

    def test_minus1_hydrophobic_without_glutamine(self):
        self.assertEqual(self.tier("MKVCDAAK"), "minus1_hydrophobic_only")

    def test_no_hydrophobic_minus1(self):
        self.assertEqual(self.tier("MKDCAAAK"), "not_detected")

    def test_no_cysteine(self):
        self.assertEqual(self.tier("MKDAAAAK"), "no_cysteine")

    def test_exact_beats_loose_when_both_present(self):
        # V-C-Q and L-C-Q in one sequence -> highest tier wins
        out = self.m.classify_motif("MKVCQAAKLCQAAK")
        self.assertEqual(out["tier"], "anchor_exact_vcq")
        self.assertGreaterEqual(len(out["positions"]), 1)


class AnchorCalibrationTests(unittest.TestCase):
    def setUp(self):
        self.m = load_module()

    def test_anchor_is_reproduced(self):
        res = self.m.verify_anchor({"CAJ92291.1": ANCHOR_SEQ}, "CAJ92291.1", 183)
        self.assertTrue(res["reproduced"])
        self.assertEqual(res["observed_positions"], [4])

    def test_missing_anchor_accession_fails_closed(self):
        res = self.m.verify_anchor({"OTHER": ANCHOR_SEQ}, "CAJ92291.1", 183)
        self.assertFalse(res["reproduced"])

    def test_pattern_absent_at_anchor_fails_closed(self):
        res = self.m.verify_anchor({"CAJ92291.1": "MKDCAAAK"}, "CAJ92291.1", 183)
        self.assertFalse(res["reproduced"])


class EndToEndTests(unittest.TestCase):
    def setUp(self):
        self.m = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_cys_motif_"))
        self.ref = self.tmp / "ref.faa"
        self.cand = self.tmp / "cand.faa"
        write_text(self.ref, REF_FASTA)
        write_text(self.cand, CAND_FASTA)
        self.out = self.tmp / "out"

    def run_annotate(self, extra=None):
        return self.m.annotate(
            reference_fasta=self.ref,
            candidate_fastas={"pool_in": self.cand},
            output_dir=self.out,
            anchor_accession="CAJ92291.1",
            anchor_position=4,  # fixture numbering
            **(extra or {}),
        )

    def test_writes_evidence_table_with_tiers(self):
        self.run_annotate()
        with (self.out / "cys_nucleophile_motif_evidence.tsv").open(
            encoding="utf-8", newline=""
        ) as h:
            rows = {r["accession"]: r for r in csv.DictReader(h, delimiter="\t")}
        self.assertEqual(rows["cand_exact"]["cys_nucleophile_motif"], "anchor_exact_vcq")
        self.assertEqual(rows["cand_loose"]["cys_nucleophile_motif"], "anchor_context_hydrophobic_c_q")
        self.assertEqual(rows["cand_minus1"]["cys_nucleophile_motif"], "minus1_hydrophobic_only")
        self.assertEqual(rows["cand_none"]["cys_nucleophile_motif"], "not_detected")
        self.assertEqual(rows["cand_nocys"]["cys_nucleophile_motif"], "no_cysteine")

    def test_no_row_ever_claims_verified_catalysis(self):
        self.run_annotate()
        with (self.out / "cys_nucleophile_motif_evidence.tsv").open(
            encoding="utf-8", newline=""
        ) as h:
            for r in csv.DictReader(h, delimiter="\t"):
                self.assertEqual(r["catalytic_activity_verified"], "false")
                self.assertEqual(r["independent_validation"], "pending")

    def test_fail_closed_when_anchor_not_reproduced(self):
        # anchor accession present but the pattern is absent -> refuse to write
        bad = self.tmp / "bad_ref.faa"
        write_text(bad, ">DED_hfam_65_0001|CAJ92291.1\nMKDCAAAK\n")
        with self.assertRaises(SystemExit):
            self.m.annotate(
                reference_fasta=bad,
                candidate_fastas={"pool_in": self.cand},
                output_dir=self.out / "bad",
                anchor_accession="CAJ92291.1",
                anchor_position=4,
            )
        self.assertFalse((self.out / "bad" / "cys_nucleophile_motif_evidence.tsv").exists())

    def test_calibration_json_records_separation(self):
        self.run_annotate()
        cal = json.loads(
            (self.out / "cys_nucleophile_motif_calibration.json").read_text(encoding="utf-8")
        )
        self.assertTrue(cal["anchor_reproduced"])
        self.assertEqual(cal["anchor_accession"], "CAJ92291.1")
        self.assertIn("reference_panel_separation", cal)


    def test_candidate_accession_keeps_full_header(self):
        """Candidate ids carry a genome prefix ('GCA_1|CTX_1_5') and must be
        emitted verbatim -- truncating at '|' breaks every downstream join."""
        cand = self.tmp / "genome.faa"
        write_text(cand, ">GCA_1|CTX_1_5\nMKVCQAAK\n")
        out = self.tmp / "out_full"
        self.m.annotate(
            reference_fasta=self.ref,
            candidate_fastas={"pool_in": cand},
            output_dir=out,
            anchor_accession="CAJ92291.1",
            anchor_position=4,
        )
        with (out / "cys_nucleophile_motif_evidence.tsv").open(
            encoding="utf-8", newline=""
        ) as h:
            rows = list(csv.DictReader(h, delimiter="\t"))
        self.assertEqual(rows[0]["accession"], "GCA_1|CTX_1_5")

    def test_reference_lookup_still_resolves_piped_accession(self):
        """The reference header is 'id|accession'; the anchor must resolve."""
        res = self.m.verify_anchor(
            {"CAJ92291.1": ANCHOR_SEQ, "DED_hfam_65_0001|CAJ92291.1": ANCHOR_SEQ},
            "CAJ92291.1",
            4,
        )
        self.assertTrue(res["reproduced"])


class CliTests(unittest.TestCase):
    """The CLI must accept repeatable POOL=PATH pairs (found by a real-data run)."""

    def setUp(self):
        self.m = load_module()
        self.tmp = Path(tempfile.mkdtemp(prefix="dsh_cys_cli_"))
        self.ref = self.tmp / "ref.faa"
        self.cand = self.tmp / "cand.faa"
        write_text(self.ref, REF_FASTA)
        write_text(self.cand, CAND_FASTA)
        self.out = self.tmp / "out"

    def test_cli_accepts_pool_equals_path(self):
        rc = self.m.main([
            "--reference-fasta", str(self.ref),
            "--candidate-fasta", f"pool_in={self.cand}",
            "--anchor-accession", "CAJ92291.1",
            "--anchor-position", "4",
            "--output-dir", str(self.out),
        ])
        self.assertEqual(rc, 0)
        self.assertTrue((self.out / "cys_nucleophile_motif_evidence.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
