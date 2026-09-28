import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "stratify_phaded_lipase_box_x1.py"
ALL_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
HYDROPHOBIC = "ACFILMVWY"
NON_HYDROPHOBIC = "DEHKNQRST"  # polar/charged; Gly and Pro are handled separately
SPECIAL = "GP"

EVIDENCE_HEADER = ["accession", "genome", "lipase_box_x1", "motif_panel_status"]


def load_module():
    spec = importlib.util.spec_from_file_location("stratify_phaded_lipase_box_x1", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_evidence(path, rows):
    lines = ["\t".join(EVIDENCE_HEADER)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class X1PolicyTests(unittest.TestCase):
    """Knoll 2009 x1 criterion: lipase/esterase x1 is polar, PHA depolymerase x1 is hydrophobic."""

    def test_twenty_amino_acids_classified_consistent_with_knoll_2009(self):
        module = load_module()
        self.assertEqual(set(ALL_AMINO_ACIDS), set(module.STANDARD_AMINO_ACIDS))
        self.assertEqual(len(module.STANDARD_AMINO_ACIDS), 20)
        self.assertEqual(module.HYDROPHOBIC_RESIDUES, frozenset(HYDROPHOBIC))
        for residue in HYDROPHOBIC:
            self.assertEqual(module.classify_x1(residue), "hydrophobic", residue)
        for residue in NON_HYDROPHOBIC:
            self.assertEqual(module.classify_x1(residue), "non_hydrophobic", residue)

    def test_glycine_and_proline_are_non_hydrophobic_with_distinct_labels(self):
        module = load_module()
        for residue in SPECIAL:
            self.assertEqual(module.classify_x1(residue), "non_hydrophobic", residue)
        self.assertEqual(module.x1_group("G"), "glycine")
        self.assertEqual(module.x1_group("P"), "proline")
        self.assertNotEqual(module.x1_group("G"), module.x1_group("P"))
        for residue in NON_HYDROPHOBIC:
            self.assertIn(module.x1_group(residue), {"acidic", "basic", "polar_uncharged"}, residue)
        for residue in HYDROPHOBIC:
            self.assertEqual(module.x1_group(residue), "hydrophobic", residue)

    def test_classification_is_exhaustive_and_disjoint(self):
        module = load_module()
        self.assertEqual(module.X1_CLASSES, ("hydrophobic", "non_hydrophobic", "not_tested"))
        hydrophobic = set(module.HYDROPHOBIC_RESIDUES)
        non_hydrophobic = set(module.NON_HYDROPHOBIC_RESIDUES)
        separately_annotated = set(module.SEPARATELY_ANNOTATED_RESIDUES)
        self.assertEqual(separately_annotated, set(SPECIAL))
        self.assertEqual(hydrophobic & non_hydrophobic, set())
        self.assertEqual(non_hydrophobic & separately_annotated, set())
        self.assertEqual(hydrophobic & separately_annotated, set())
        self.assertEqual(hydrophobic | non_hydrophobic | separately_annotated, set(ALL_AMINO_ACIDS))
        for residue in ALL_AMINO_ACIDS:
            self.assertIn(module.classify_x1(residue), module.X1_CLASSES)

    def test_unknown_or_ambiguous_residue_fails_closed(self):
        module = load_module()
        for residue in ("X", "B", "Z", "J", "U", "O", "*", "s", "AD"):
            with self.assertRaises(ValueError, msg=residue):
                module.classify_x1(residue)

    def test_missing_residue_is_not_tested_and_never_a_negative(self):
        module = load_module()
        for missing in ("", "   ", "\t"):
            self.assertEqual(module.classify_x1(missing), "not_tested", repr(missing))
        self.assertEqual(module.x1_group(""), "not_tested")
        self.assertEqual(module.NEGATIVE_CLASSES, frozenset())
        for x1_class in module.X1_CLASSES:
            self.assertFalse(module.is_biological_negative(x1_class), x1_class)
        self.assertEqual(module.confounder_flag("not_tested"), "not_assessed")
        self.assertEqual(module.confounder_flag("hydrophobic"), "none")

    def test_non_hydrophobic_cap_never_allows_high_candidate_only(self):
        module = load_module()
        self.assertEqual(module.recommended_confidence_cap("hydrophobic"), "high_candidate_only")
        for residue in NON_HYDROPHOBIC + SPECIAL:
            x1_class = module.classify_x1(residue)
            cap = module.recommended_confidence_cap(x1_class)
            self.assertNotEqual(cap, "high_candidate_only", residue)
        self.assertEqual(
            module.recommended_confidence_cap("non_hydrophobic"), module.CAP_NON_HYDROPHOBIC
        )
        # an untested residue is uninformative, so it must not be capped like a passing hydrophobic x1
        untested_cap = module.recommended_confidence_cap("not_tested")
        self.assertNotIn(untested_cap, {"high_candidate_only", module.CAP_NON_HYDROPHOBIC})
        with self.assertRaises(ValueError):
            module.recommended_confidence_cap("unknown_class")

    def test_confounder_flag_only_for_non_hydrophobic(self):
        module = load_module()
        self.assertEqual(
            module.confounder_flag("non_hydrophobic"), module.CONFOUNDER_FLAG_NON_HYDROPHOBIC
        )
        self.assertEqual(module.confounder_flag("hydrophobic"), "none")
        self.assertEqual(module.confounder_flag("not_tested"), "not_assessed")
        with self.assertRaises(ValueError):
            module.confounder_flag("unknown_class")


class StratifyRunTests(unittest.TestCase):
    def _synthetic(self, root):
        evidence = root / "motif_candidate_evidence.tsv"
        write_evidence(
            evidence,
            [
                ["p6", "G1", "L", "supported"],
                ["p2", "G1", "D", "supported"],
                ["p5", "G1", "", "not_detected_pattern"],
                ["p1", "G2", "A", "supported"],
                ["p4", "G2", "P", "supported"],
                ["p3", "G3", "G", "supported"],
            ],
        )
        return evidence

    def test_end_to_end_outputs_are_sorted_and_self_consistent(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            evidence = self._synthetic(root)
            before = evidence.read_bytes()
            output = root / "results"
            manifest = module.stratify(evidence, output)

            lines = (output / "x1_stratification.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                lines[0],
                "\t".join(
                    [
                        "accession",
                        "lipase_box_x1",
                        "x1_class",
                        "confounder_flag",
                        "recommended_confidence_cap",
                        "x1_group",
                    ]
                ),
            )
            accessions = [line.split("\t")[0] for line in lines[1:]]
            self.assertEqual(accessions, ["p1", "p2", "p3", "p4", "p5", "p6"])
            by_accession = {line.split("\t")[0]: line.split("\t") for line in lines[1:]}
            self.assertEqual(by_accession["p1"][2], "hydrophobic")
            self.assertEqual(by_accession["p2"][2], "non_hydrophobic")
            self.assertEqual(by_accession["p2"][3], module.CONFOUNDER_FLAG_NON_HYDROPHOBIC)
            self.assertEqual(by_accession["p3"][5], "glycine")
            self.assertEqual(by_accession["p4"][5], "proline")
            self.assertEqual(by_accession["p5"][2], "not_tested")
            self.assertEqual(by_accession["p5"][4], module.CAP_NOT_TESTED)
            self.assertEqual(by_accession["p6"][4], "high_candidate_only")
            self.assertEqual(evidence.read_bytes(), before)

            self.assertEqual(manifest["total_candidates"], 6)
            self.assertEqual(manifest["tested_count"], 5)
            self.assertEqual(manifest["hydrophobic_count"], 2)
            self.assertEqual(manifest["non_hydrophobic_count"], 3)
            self.assertEqual(manifest["not_tested_count"], 1)
            self.assertEqual(manifest["negative_call_count"], 0)
            self.assertEqual(
                manifest["hydrophobic_count"]
                + manifest["non_hydrophobic_count"]
                + manifest["not_tested_count"],
                manifest["total_candidates"],
            )
            self.assertEqual(manifest["residue_counts"]["D"], 1)

            summary = (output / "x1_summary.tsv").read_text(encoding="utf-8")
            self.assertIn("total_candidates\t6", summary)
            self.assertIn("non_hydrophobic\t3", summary)
            self.assertIn("not_tested\t1", summary)
            self.assertIn("negative_call_count\t0", summary)

            confounders = (output / "confounder_candidates.tsv").read_text(encoding="utf-8")
            confounder_lines = confounders.splitlines()
            self.assertEqual(len(confounder_lines), 4)
            self.assertEqual(
                [line.split("\t")[0] for line in confounder_lines[1:]], ["p2", "p3", "p4"]
            )
            self.assertEqual(confounder_lines[0].split("\t")[:2], ["accession", "genome"])
            self.assertNotIn("p5", confounders)
            self.assertFalse((output / "x1_manifest.json").read_text(encoding="utf-8").count("pending"))

    def test_unknown_residue_and_duplicate_accession_fail_closed(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            bad = root / "bad.tsv"
            write_evidence(bad, [["p1", "G1", "X", "supported"]])
            with self.assertRaises(ValueError):
                module.stratify(bad, root / "out_bad")
            duplicate = root / "dup.tsv"
            write_evidence(duplicate, [["p1", "G1", "A", "supported"], ["p1", "G1", "L", "supported"]])
            with self.assertRaises(ValueError):
                module.stratify(duplicate, root / "out_dup")

    def test_missing_required_column_fails_closed(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            broken = root / "broken.tsv"
            broken.write_text("accession\tgenome\np1\tG1\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                module.stratify(broken, root / "out")


if __name__ == "__main__":
    unittest.main()
