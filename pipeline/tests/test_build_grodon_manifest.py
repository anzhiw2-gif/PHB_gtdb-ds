import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "deploy" / "20260920_phaded_grodon_growth_01" / "scripts" / "build_grodon_manifest.py"

V2_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "grodon_reanalysis_v2.py"


def load_module():
    spec = importlib.util.spec_from_file_location("build_grodon_manifest", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_v2_module():
    spec = importlib.util.spec_from_file_location("grodon_reanalysis_v2_manifest", V2_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tax(genus="TestGenus", domain="Bacteria", phylum="Pseudomonadota", cls="Gammaproteobacteria", order="TestOrder", family="TestFamily", species="sp"):
    return {
        "domain": domain,
        "phylum": phylum,
        "class": cls,
        "order": order,
        "family": family,
        "genus": genus,
        "species": species,
    }


def gid(i):
    # 9-digit GTDB accession: GCA_000000001.1
    return f"GCA_000000{i:03d}.1"


class BuildGrodonManifestTests(unittest.TestCase):
    def test_taxonomy_string_parsing(self):
        module = load_module()
        parsed = module.parse_taxonomy_string(
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Test;f__Fam;g__GenA;s__sp"
        )
        self.assertEqual(parsed["domain"], "Bacteria")
        self.assertEqual(parsed["genus"], "GenA")
        self.assertEqual(parsed["phylum"], "Pseudomonadota")

    def test_genome_path_pattern(self):
        module = load_module()
        path = module.genome_path_from_id("GCA_000000001.1")
        normalized = path.replace("\\", "/")
        self.assertTrue(normalized.endswith("GCA_000000001.1_genomic.fna.gz"))
        self.assertIn("/GCA/000/000/001/", normalized)
        with self.assertRaises(ValueError):
            module.genome_path_from_id("not_an_accession")

    def test_manifest_balances_within_genus_and_excludes_pool_genomes(self):
        module = load_module()
        degrader = {
            gid(1): {"genome_id": gid(1), "major_superfamily": "intracellular nPHASCL without lipase box", "n_candidates": "3", "group": "intracellular"},
            gid(2): {"genome_id": gid(2), "major_superfamily": "extracellular dPHASCL type 1", "n_candidates": "1", "group": "extracellular"},
        }
        exclusion = {gid(4)}
        tax_by_gid = {gid(i): tax() for i in range(1, 6)}
        rows, stats = module.build_manifest(
            degrader, exclusion, tax_by_gid, seed=42, max_per_genus=0, file_exists=lambda p: True
        )
        # genus has 2 pos, 2 usable neg (gid3, gid5; gid4 excluded) -> balanced 2/2
        self.assertEqual(stats["manifest_rows"], 4)
        self.assertEqual(stats["manifest_positive"], 2)
        self.assertEqual(stats["manifest_negative"], 2)
        statuses = {r["genome_id"]: r["phaZ_status"] for r in rows}
        self.assertEqual(statuses[gid(1)], "phaZ_positive")
        neg_ids = {r["genome_id"] for r in rows if r["phaZ_status"] == "phaZ_negative"}
        self.assertNotIn(gid(4), neg_ids)
        self.assertIn(gid(3), neg_ids)
        for r in rows:
            if r["phaZ_status"] == "phaZ_negative":
                self.assertEqual(r["group"], "control")
                self.assertEqual(r["major_subtype"], "none")

    def test_manifest_filters_genomes_without_fasta(self):
        module = load_module()
        degrader = {
            gid(1): {"genome_id": gid(1), "major_superfamily": "X", "n_candidates": "2", "group": "intracellular"},
            gid(2): {"genome_id": gid(2), "major_superfamily": "X", "n_candidates": "1", "group": "intracellular"},
        }
        tax_by_gid = {gid(i): tax() for i in range(1, 5)}
        # gid(2) FASTA missing; gid(3) control FASTA missing
        missing = {gid(2), gid(3)}
        rows, stats = module.build_manifest(
            degrader, set(), tax_by_gid, seed=42, max_per_genus=0,
            file_exists=lambda p: Path(p).name.replace("_genomic.fna.gz", "") not in missing,
        )
        # only gid(1) usable pos, gid(4) usable neg -> 1 pair
        self.assertEqual(stats["manifest_rows"], 2)
        self.assertEqual(stats["manifest_positive"], 1)
        self.assertEqual(stats["degrader_genomes_missing_fasta"], 1)
        pos = [r["genome_id"] for r in rows if r["phaZ_status"] == "phaZ_positive"]
        self.assertEqual(pos, [gid(1)])

    def test_manifest_skips_genus_without_controls(self):
        module = load_module()
        degrader = {gid(1): {"genome_id": gid(1), "major_superfamily": "X", "n_candidates": "1", "group": "intracellular"}}
        tax_by_gid = {gid(1): tax("LonelyGenus"), gid(2): tax("OtherGenus")}
        rows, stats = module.build_manifest(
            degrader, set(), tax_by_gid, seed=42, max_per_genus=0, file_exists=lambda p: True
        )
        self.assertEqual(rows, [])
        self.assertEqual(stats["degrader_genomes_excluded_no_control"], 1)
        self.assertEqual(stats["genera_skipped_no_control"], 1)

    def test_manifest_seed_determinism(self):
        module = load_module()
        degrader = {gid(i): {"genome_id": gid(i), "major_superfamily": "X", "n_candidates": "1", "group": "intracellular"} for i in range(1, 6)}
        tax_by_gid = {gid(i): tax() for i in range(1, 12)}
        rows_a, _ = module.build_manifest(
            degrader, set(), tax_by_gid, seed=42, max_per_genus=0, file_exists=lambda p: True
        )
        rows_b, _ = module.build_manifest(
            degrader, set(), tax_by_gid, seed=42, max_per_genus=0, file_exists=lambda p: True
        )
        ids_a = [r["genome_id"] for r in rows_a]
        ids_b = [r["genome_id"] for r in rows_b]
        self.assertEqual(ids_a, ids_b)
        self.assertEqual(len(rows_a), 10)


class BuildGrodonManifestV2Tests(unittest.TestCase):
    """Task 12 additions: the v2 module supersedes the legacy deploy script.

    Everything asserted here is against ``pipeline/scripts/grodon_reanalysis_v2.py``.
    The legacy deploy tests above are untouched on purpose: the frozen
    2026-09-20 deploy is historical evidence and must keep behaving exactly as it
    did, while the dated v2 deploy copies the corrected module.
    """

    def v2_carriers(self, module, count=3):
        return {
            gid(i): {
                "candidate_family": "DED_hfam_70",
                "candidate_count": i,
                "group": "intracellular" if i % 2 else "extracellular",
            }
            for i in range(1, count + 1)
        }

    def test_v2_manifest_emits_candidate_carrier_language(self):
        module = load_v2_module()
        carriers = self.v2_carriers(module, 2)
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 8)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=2,
            file_exists=lambda _: True,
        )
        self.assertEqual(
            {row["candidate_detection_status"] for row in rows},
            {"candidate_gene_carrier", "candidate_not_detected_under_defined_search"},
        )
        self.assertEqual(stats["manifest_carrier"], 2)
        self.assertEqual(stats["manifest_control"], 2)
        self.assertNotIn("phaZ_status", module.MANIFEST_COLUMNS)
        for row in rows:
            self.assertNotIn(
                "degrader", " ".join(str(value) for value in row.values()).lower()
            )

    def test_v2_manifest_asserts_zero_positive_control_overlap(self):
        module = load_v2_module()
        carriers = self.v2_carriers(module, 2)
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 8)}
        rows, _ = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=2,
            file_exists=lambda _: True,
        )
        self.assertEqual(
            module.carrier_genomes(rows) & module.control_genomes(rows), set()
        )
        self.assertTrue(module.sets_are_disjoint(rows))
        with self.assertRaisesRegex(ValueError, "positive/control overlap"):
            module.build_manifest(
                carriers=carriers,
                exclusion={gid(2)},
                tax_by_gid=tax_by_gid,
                seed=42,
                max_per_genus=2,
                file_exists=lambda _: True,
            )

    def test_v2_manifest_keeps_the_proven_balancing_and_fasta_rules(self):
        module = load_v2_module()
        carriers = self.v2_carriers(module, 2)
        tax_by_gid = {gid(i): tax("GenA") for i in range(1, 8)}
        missing = {gid(4)}
        rows, stats = module.build_manifest(
            carriers=carriers,
            exclusion=set(),
            tax_by_gid=tax_by_gid,
            seed=42,
            max_per_genus=1,
            file_exists=lambda path: Path(path).name.replace(
                "_genomic.fna.gz", ""
            ) not in missing,
        )
        controls = module.control_genomes(rows)
        self.assertTrue(controls)
        self.assertEqual(controls & missing, set())
        self.assertEqual(stats["matched_pairs"], 1)
        self.assertGreaterEqual(stats["control_genomes_missing_fasta"], 0)

    def test_v2_manifest_never_reuses_a_legacy_phenotype_column(self):
        module = load_v2_module()
        self.assertIn("candidate_detection_status", module.MANIFEST_COLUMNS)
        self.assertIn("candidate_count", module.MANIFEST_COLUMNS)
        self.assertNotIn("phaZ_status", module.MANIFEST_COLUMNS)
        self.assertNotIn("phaZ_validated_count", module.MANIFEST_COLUMNS)


if __name__ == "__main__":
    unittest.main()
