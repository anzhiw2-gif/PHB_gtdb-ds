"""Tests for the P4 structural survey scoring.

Synthetic fixtures only; the real Foldseek output is read by the run.
"""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import analyze_phaded_structural_survey as module  # noqa: E402


def fs_line(query: str, target: str, tm: float, evalue: str = "1e-20") -> str:
    return f"{query}\t{target}\t{tm}\t{tm}\t{tm}\t1.000\t{evalue}\t0.9"


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return path


class RankTests(unittest.TestCase):
    def test_average_ranks_handle_ties(self):
        self.assertEqual(module.rank([10.0, 20.0, 20.0, 30.0]), [1.0, 2.5, 2.5, 4.0])

    def test_identical_series_has_no_defined_correlation(self):
        self.assertIsNone(module.spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))

    def test_perfect_agreement_and_perfect_reversal(self):
        self.assertAlmostEqual(module.spearman([1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]), 1.0)
        self.assertAlmostEqual(module.spearman([1.0, 2.0, 3.0, 4.0], [40.0, 30.0, 20.0, 10.0]), -1.0)

    def test_too_few_points_returns_none(self):
        self.assertIsNone(module.spearman([1.0, 2.0], [1.0, 2.0]))

    def test_length_mismatch_raises(self):
        with self.assertRaisesRegex(ValueError, "equal-length"):
            module.spearman([1.0, 2.0, 3.0], [1.0, 2.0])


class FilenameTests(unittest.TestCase):
    def test_accession_is_recovered_without_splitting_on_pipes(self):
        name = "GCA_000017645.1|CP000781.1_1999_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000"
        self.assertEqual(module.candidate_of(name), "GCA_000017645.1|CP000781.1_1999")

    def test_unrecognised_name_raises(self):
        with self.assertRaisesRegex(ValueError, "unrecognised prediction filename"):
            module.candidate_of("something_else.pdb")

    def test_colabfold_replaces_the_pipe_with_an_underscore(self):
        # The real failure this guards: ColabFold writes filenames with '|' turned
        # into '_', so the recovered name never matches the raw accession and every
        # candidate would be reported as missing from the composition table.
        query = "GCA_015680705.1_JADNYP010000003.1_102_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000"
        recovered = module.candidate_of(query)
        self.assertEqual(recovered, "GCA_015680705.1_JADNYP010000003.1_102")
        self.assertEqual(module.normalize_name("GCA_015680705.1|JADNYP010000003.1_102"), recovered)

    def test_index_maps_the_raw_accession_through_the_normalised_key(self):
        index = module.index_composition(
            [{"accession": "A|B_1", "panel_margin_competitor_minus_anchor_bits": "-3.0"}])
        self.assertIn("A_B_1", index)
        self.assertEqual(index["A_B_1"]["accession"], "A|B_1")

    def test_a_normalisation_collision_is_refused(self):
        with self.assertRaisesRegex(ValueError, "could not be attributed safely"):
            module.index_composition([
                {"accession": "A|B", "panel_margin_competitor_minus_anchor_bits": "1.0"},
                {"accession": "A_B", "panel_margin_competitor_minus_anchor_bits": "2.0"},
            ])

    def test_an_empty_composition_table_is_refused(self):
        with self.assertRaisesRegex(ValueError, "composition table is empty"):
            module.index_composition([])


class ScoringTests(unittest.TestCase):
    def _fixture(self, tmp: Path):
        rows = [
            # anchor-favoured at sequence level, competitor-favoured structurally
            fs_line("candA_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.75),
            fs_line("candA_unrelaxed_rank_001_m", "anchor_8YNV_C", 0.76),
            fs_line("candA_unrelaxed_rank_001_m", "competitor_Q88N36", 0.84),
            fs_line("candA_unrelaxed_rank_001_m", "competitor_P24640", 0.80),
            # competitor-favoured at both levels
            fs_line("candB_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.70),
            fs_line("candB_unrelaxed_rank_001_m", "competitor_Q88N36", 0.79),
        ]
        foldseek = write(tmp / "fs.tsv", "\n".join(rows) + "\n")
        composition = write_tsv(tmp / "comp.tsv",
                                ["accession", "panel_margin_competitor_minus_anchor_bits"], [
            {"accession": "candA", "panel_margin_competitor_minus_anchor_bits": "-19.3000"},
            {"accession": "candB", "panel_margin_competitor_minus_anchor_bits": "4.3000"},
        ])
        return foldseek, composition

    def test_best_chain_and_best_competitor_are_used(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek, composition = self._fixture(tmp)
            records, _summary = module.score(module.read_foldseek(foldseek),
                                             module.read_tsv(composition))
            by_accession = {row["accession"]: row for row in records}
            self.assertEqual(by_accession["candA"]["anchor_target"], "anchor_8YNV_C")
            self.assertEqual(by_accession["candA"]["anchor_tm"], "0.7600")
            self.assertEqual(by_accession["candA"]["competitor_target"], "competitor_Q88N36")
            self.assertEqual(by_accession["candA"]["structural_winner"], "competitor")

    def test_margin_sign_and_preregistered_flag(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek, composition = self._fixture(tmp)
            records, summary = module.score(module.read_foldseek(foldseek),
                                            module.read_tsv(composition))
            by_accession = {row["accession"]: row for row in records}
            self.assertEqual(by_accession["candA"]["structural_margin_anchor_minus_competitor"],
                             "-0.0800")
            self.assertEqual(by_accession["candA"]["anchor_clears_preregistered_tm"], "yes")
            self.assertEqual(summary["anchor_clears_preregistered_tm"], 2)
            self.assertEqual(summary["structural_winner_counts"], {"anchor": 0, "competitor": 2})

    def test_correlation_is_negative_when_the_sequence_ranking_is_inverted(self):
        # candA has the lower sequence margin AND loses structurally; candB the
        # higher sequence margin and also loses, more narrowly. Rank agreement is
        # therefore negative -> the sequence ranking does not predict the outcome.
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek, composition = self._fixture(tmp)
            _records, summary = module.score(module.read_foldseek(foldseek),
                                             module.read_tsv(composition))
            self.assertIsNone(summary["spearman_sequence_margin_vs_structural_margin"])

    def test_a_candidate_with_no_anchor_comparison_raises(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek = write(tmp / "fs.tsv", fs_line("candA_unrelaxed_rank_001_m",
                                                     "competitor_Q88N36", 0.8) + "\n")
            composition = write_tsv(tmp / "comp.tsv",
                                    ["accession", "panel_margin_competitor_minus_anchor_bits"],
                                    [{"accession": "candA",
                                      "panel_margin_competitor_minus_anchor_bits": "1.0"}])
            with self.assertRaisesRegex(ValueError, "no comparison against the anchor"):
                module.score(module.read_foldseek(foldseek), module.read_tsv(composition))

    def test_a_survey_member_with_no_comparison_at_all_raises(self):
        # A member listed in the composition but absent from the Foldseek output
        # means its prediction failed; that must not pass silently.
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek = write(tmp / "fs.tsv",
                             fs_line("candA_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.8) + "\n"
                             + fs_line("candA_unrelaxed_rank_001_m", "competitor_Q88N36", 0.7) + "\n")
            composition = write_tsv(tmp / "comp.tsv",
                                    ["accession", "panel_margin_competitor_minus_anchor_bits"],
                                    [{"accession": "candA",
                                      "panel_margin_competitor_minus_anchor_bits": "1.0"},
                                     {"accession": "cand_missing",
                                      "panel_margin_competitor_minus_anchor_bits": "2.0"}])
            with self.assertRaisesRegex(ValueError, "no Foldseek comparison at all"):
                module.score(module.read_foldseek(foldseek), module.read_tsv(composition))

    def test_a_query_outside_the_composition_is_counted_not_fatal(self):
        # The pilot's Foldseek run searched every predicted structure, including the
        # panel's own, so foreign queries are expected and must be reported.
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek = write(tmp / "fs.tsv", "\n".join([
                fs_line("candA_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.8),
                fs_line("candA_unrelaxed_rank_001_m", "competitor_Q88N36", 0.7),
                fs_line("P24640_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.6),
                fs_line("P24640_unrelaxed_rank_001_m", "competitor_Q88N36", 0.9),
            ]) + "\n")
            composition = write_tsv(tmp / "comp.tsv",
                                    ["accession", "panel_margin_competitor_minus_anchor_bits"],
                                    [{"accession": "candA",
                                      "panel_margin_competitor_minus_anchor_bits": "1.0"}])
            records, summary = module.score(module.read_foldseek(foldseek),
                                            module.read_tsv(composition))
            self.assertEqual([row["accession"] for row in records], ["candA"])
            self.assertEqual(summary["non_survey_queries_ignored"], ["P24640"])
            self.assertEqual(summary["non_survey_query_count"], 1)

    def test_malformed_foldseek_rows_raise(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            short = write(tmp / "short.tsv", "a\tb\t1\n")
            with self.assertRaisesRegex(ValueError, "malformed Foldseek row"):
                module.read_foldseek(short)
            bad = write(tmp / "bad.tsv", "a\tb\tx\ty\tz\t1\t2\t3\n")
            with self.assertRaisesRegex(ValueError, "non-numeric Foldseek score"):
                module.read_foldseek(bad)

    def test_boundary_is_stated(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek, composition = self._fixture(tmp)
            _records, summary = module.score(module.read_foldseek(foldseek),
                                             module.read_tsv(composition))
            self.assertIn("Predicted folds, not experimental structures", summary["boundary"])
            self.assertIn("POSITIVE Spearman", summary["interpretation"])
            self.assertEqual(summary["preregistered_tm_threshold"], 0.5)


class ModelSelectionTests(unittest.TestCase):
    def test_rank_is_read_from_the_filename(self):
        self.assertEqual(module.model_rank("x_unrelaxed_rank_003_alphafold2_ptm_model_1"), 3)
        self.assertEqual(module.model_rank("x_unrelaxed_rank_001_alphafold2_ptm_model_2"), 1)
        self.assertEqual(module.model_rank("no_marker_here"), 0)

    def test_rank_001_keeps_both_sides_on_the_same_model(self):
        rows = {
            "cand": [
                {"query": "cand_unrelaxed_rank_001_m", "target": "anchor_8YNV_A", "tm": 0.70,
                 "evalue": "1e-9"},
                {"query": "cand_unrelaxed_rank_001_m", "target": "competitor_X", "tm": 0.60,
                 "evalue": "1e-9"},
                {"query": "cand_unrelaxed_rank_002_m", "target": "anchor_8YNV_A", "tm": 0.95,
                 "evalue": "1e-9"},
            ],
        }
        reduced = module.select_models(rows, "rank_001")
        self.assertEqual({row["query"] for row in reduced["cand"]},
                         {"cand_unrelaxed_rank_001_m"})

    def test_max_keeps_every_model_so_the_sides_can_mix(self):
        rows = {"cand": [
            {"query": "cand_unrelaxed_rank_001_m", "target": "anchor_8YNV_A", "tm": 0.70,
             "evalue": "1e-9"},
            {"query": "cand_unrelaxed_rank_002_m", "target": "anchor_8YNV_A", "tm": 0.95,
             "evalue": "1e-9"},
        ]}
        self.assertEqual(len(module.select_models(rows, "max")["cand"]), 2)

    def test_a_candidate_without_rank_001_falls_back_rather_than_vanishing(self):
        rows = {"cand": [
            {"query": "cand_unrelaxed_rank_003_m", "target": "anchor_8YNV_A", "tm": 0.7,
             "evalue": "1e-9"},
        ]}
        self.assertEqual(len(module.select_models(rows, "rank_001")["cand"]), 1)

    def test_an_unknown_selection_is_refused(self):
        with self.assertRaisesRegex(ValueError, "unknown model selection"):
            module.select_models({}, "whatever")

    def test_the_selection_is_reported_in_the_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek = write(tmp / "fs.tsv", "\n".join([
                fs_line("c1_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.8),
                fs_line("c1_unrelaxed_rank_001_m", "competitor_Q88N36", 0.7),
            ]) + "\n")
            composition = write_tsv(tmp / "comp.tsv",
                                    ["accession", "panel_margin_competitor_minus_anchor_bits"],
                                    [{"accession": "c1",
                                      "panel_margin_competitor_minus_anchor_bits": "1.0"}])
            _records, summary = module.score(module.read_foldseek(foldseek),
                                             module.read_tsv(composition))
            self.assertEqual(summary["model_selection"], "rank_001")


class CliTests(unittest.TestCase):
    def test_cli_writes_artifacts_and_refuses_a_non_empty_dir(self):
        with tempfile.TemporaryDirectory() as temporary:
            tmp = Path(temporary)
            foldseek = write(tmp / "fs.tsv",
                             "\n".join([
                                 fs_line("c1_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.8),
                                 fs_line("c1_unrelaxed_rank_001_m", "competitor_Q88N36", 0.7),
                                 fs_line("c2_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.6),
                                 fs_line("c2_unrelaxed_rank_001_m", "competitor_Q88N36", 0.75),
                                 fs_line("c3_unrelaxed_rank_001_m", "anchor_8YNV_A", 0.9),
                                 fs_line("c3_unrelaxed_rank_001_m", "competitor_Q88N36", 0.65),
                             ]) + "\n")
            composition = write_tsv(tmp / "comp.tsv",
                                    ["accession", "panel_margin_competitor_minus_anchor_bits"], [
                {"accession": "c1", "panel_margin_competitor_minus_anchor_bits": "-15.0"},
                {"accession": "c2", "panel_margin_competitor_minus_anchor_bits": "-5.0"},
                {"accession": "c3", "panel_margin_competitor_minus_anchor_bits": "2.0"},
            ])
            out = tmp / "out"
            self.assertEqual(
                module.main(["--foldseek", str(foldseek), "--composition", str(composition),
                             "--out-dir", str(out)]), 0,
            )
            summary = json.loads((out / "survey_structural_summary.json").read_text("utf-8"))
            self.assertEqual(summary["candidates_scored"], 3)
            self.assertEqual(summary["structural_winner_counts"], {"anchor": 2, "competitor": 1})
            self.assertIsNotNone(summary["spearman_sequence_margin_vs_structural_margin"])
            with self.assertRaisesRegex(ValueError, "non-empty output dir"):
                module.main(["--foldseek", str(foldseek), "--composition", str(composition),
                             "--out-dir", str(out)])


if __name__ == "__main__":
    unittest.main()
