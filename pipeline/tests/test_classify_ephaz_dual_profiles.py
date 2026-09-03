import unittest

from pipeline.scripts.classify_ephaz_dual_profiles import (
    classify_accession,
    classify_hits,
    parse_domtblout,
    parse_fasta_accession,
)


class DualProfileClassificationTests(unittest.TestCase):
    def test_separates_profile_specific_hits(self):
        phb = {"P": {"evalue": 1e-40, "bitscore": 120.0}, "M": {"evalue": 1e-6, "bitscore": 40.0}}
        mcl = {"P": {"evalue": 1e-3, "bitscore": 20.0}, "M": {"evalue": 1e-40, "bitscore": 130.0}}
        self.assertEqual(classify_accession(phb["P"], mcl["P"]), "PHB_like")
        self.assertEqual(classify_accession(phb["M"], mcl["M"]), "MCL_like")

    def test_marks_close_dual_hits_ambiguous(self):
        self.assertEqual(
            classify_accession({"evalue": 1e-20, "bitscore": 100.0}, {"evalue": 1e-20, "bitscore": 96.0}, margin=10.0),
            "ambiguous",
        )

    def test_no_hit_and_threshold_are_explicit(self):
        self.assertEqual(classify_accession(None, None), "no_hit")
        self.assertEqual(classify_accession({"evalue": 1e-4, "bitscore": 200.0}, {"evalue": 1e-30, "bitscore": 50.0}), "MCL_like")

    def test_classify_hits_is_deterministic(self):
        rows = classify_hits({"B": {"evalue": 1e-20, "bitscore": 100}, "A": {"evalue": 1e-20, "bitscore": 100}}, {}, threshold=1e-5)
        self.assertEqual([row["accession"] for row in rows], ["A", "B"])
        self.assertTrue(all(row["classification"] == "PHB_like" for row in rows))

    def test_legacy_single_mcl_profile_keeps_mcl_like_label(self):
        rows = classify_hits({}, {"M": {"evalue": 1e-30, "bitscore": 100.0}}, accessions=["M"])
        self.assertEqual(rows[0]["classification"], "MCL_like")

    def test_includes_unhit_probe_accessions_as_no_hit(self):
        rows = classify_hits({"HIT": {"evalue": 1e-20, "bitscore": 100}}, {}, accessions=["MISS", "HIT"])
        self.assertEqual([row["accession"] for row in rows], ["HIT", "MISS"])
        self.assertEqual(rows[1]["classification"], "no_hit")

    def test_normalizes_uniprot_fasta_headers_without_losing_accession(self):
        self.assertEqual(parse_fasta_accession("sp|P12345|PHAZ_EXAMPLE description"), "P12345")
        self.assertEqual(parse_fasta_accession("tr|A0A123|PHAZ_EXAMPLE description"), "A0A123")
        self.assertEqual(parse_fasta_accession("AZSS01000334.1:12616-13485(-)|candidate"), "AZSS01000334.1:12616-13485(-)")

    def test_classifies_named_mcl_subfamilies_and_reports_best_profile(self):
        rows = classify_hits(
            {"P": {"evalue": 1e-20, "bitscore": 30.0}, "W": {"evalue": 1e-6, "bitscore": 20.0}},
            {"classical": {"P": {"evalue": 1e-30, "bitscore": 120.0}}, "lipase_associated": {"W": {"evalue": 1e-40, "bitscore": 130.0}}},
            accessions=["P", "W"],
        )
        by_accession = {row["accession"]: row for row in rows}
        self.assertEqual(by_accession["P"]["classification"], "MCL_classical")
        self.assertEqual(by_accession["P"]["best_mcl_subfamily"], "classical")
        self.assertEqual(by_accession["W"]["classification"], "MCL_lipase_associated")
        self.assertEqual(by_accession["W"]["best_mcl_subfamily"], "lipase_associated")
        self.assertIn("mcl_classical_evalue", by_accession["P"])

    def test_parse_domtblout_merges_hmm_intervals_into_coverage(self):
        row = "target|panel - 200 profile - 100 1e-20 80 0.0 1 2 1e-20 1e-20 80 0.0 1 40 1 40 1 40 0.99"
        row2 = "target|panel - 200 profile - 100 1e-20 70 0.0 2 2 1e-10 1e-10 70 0.0 51 100 50 100 50 100 0.99"
        hits = parse_domtblout(row + "\n" + row2)
        self.assertAlmostEqual(hits["target"]["coverage"], 0.9)
        self.assertEqual(hits["target"]["bitscore"], 80.0)

    def test_coverage_filter_is_opt_in_and_rejects_partial_hit(self):
        hit = {"evalue": 1e-40, "bitscore": 200.0, "coverage": 0.637}
        self.assertEqual(classify_accession(None, hit), "MCL_like")
        self.assertEqual(classify_accession(None, hit, min_domain_coverage=0.9), "no_hit")

    def test_named_profile_coverage_filter_preserves_profile_label(self):
        rows = classify_hits(
            {},
            {"Streptomyces": {"J7K890": {"evalue": 1e-115, "bitscore": 375.6, "coverage": 0.637}}},
            accessions=["J7K890"],
            min_domain_coverage=0.9,
        )
        self.assertEqual(rows[0]["classification"], "no_hit")

    def test_profile_label_does_not_duplicate_mcl_prefix(self):
        rows = classify_hits(
            {},
            {"mcl_classical": {"P": {"evalue": 1e-20, "bitscore": 50.0, "coverage": 1.0}}},
            accessions=["P"],
        )
        self.assertEqual(rows[0]["classification"], "MCL_classical")


if __name__ == "__main__":
    unittest.main()
