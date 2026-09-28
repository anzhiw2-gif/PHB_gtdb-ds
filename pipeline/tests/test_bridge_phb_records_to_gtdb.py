import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "bridge_phb_records_to_gtdb.py"


class GtdbBridgeTests(unittest.TestCase):
    def _load_module(self):
        spec = importlib.util.spec_from_file_location("bridge_phb_records_to_gtdb", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def _write_tsv(self, path, fields, rows):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_exact_genome_accession_is_direct_bridge(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature = root / "literature.tsv"
            taxonomy = root / "taxonomy.tsv"
            output = root / "bridge.tsv"
            self._write_tsv(literature, ["record_id", "organism_name", "strain_name", "protein_accession", "gene_accession", "genome_accession", "ncbi_taxid", "gtdb_taxonomy", "source_database", "source_database_version", "retrieval_date", "primary_doi", "pmid", "pmcid", "experimental_substrate", "experimental_result", "evidence_type", "evidence_status", "exact_identifier_scope", "sequence_sha256", "notes"], [{
                "record_id": "lit_1", "organism_name": "Paucimonas lemoignei", "strain_name": "PHA1", "protein_accession": "Q51871", "gene_accession": "", "genome_accession": "GCA_000123.1", "ncbi_taxid": "", "gtdb_taxonomy": "", "source_database": "PubMed", "source_database_version": "2026-09-13", "retrieval_date": "2026-09-13", "primary_doi": "10.1/example", "pmid": "1", "pmcid": "", "experimental_substrate": "PHB", "experimental_result": "degraded", "evidence_type": "primary_experiment", "evidence_status": "direct_literature_supported", "exact_identifier_scope": "genome", "sequence_sha256": "", "notes": "",
            }])
            self._write_tsv(taxonomy, ["genome_accession", "gtdb_taxonomy"], [{"genome_accession": "GCA_000123.1", "gtdb_taxonomy": "d__Bacteria;p__Example"}])
            rows = module.bridge_records(literature, taxonomy, output)
            self.assertEqual(rows[0]["match_class"], "exact_genome_accession")
            self.assertEqual(rows[0]["evidence_tier"], "direct_literature_supported")
            self.assertEqual(rows[0]["gtdb_genome_accession"], "GCA_000123.1")

    def test_experimental_positive_with_ncbi_assembly_alias_is_direct_bridge(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature = root / "literature.tsv"
            taxonomy = root / "taxonomy.tsv"
            output = root / "bridge.tsv"
            fields = ["record_id", "organism_name", "strain_name", "protein_accession", "gene_accession", "genome_accession", "ncbi_taxid", "gtdb_taxonomy", "source_database", "source_database_version", "retrieval_date", "primary_doi", "pmid", "pmcid", "experimental_substrate", "experimental_result", "evidence_type", "evidence_status", "exact_identifier_scope", "sequence_sha256", "notes"]
            self._write_tsv(literature, fields, [{"record_id": "lit_alias", "organism_name": "Cupriavidus necator", "strain_name": "H16", "protein_accession": "Q0K9H3", "gene_accession": "H16_A2251", "genome_accession": "NC_008313.1", "ncbi_taxid": "381666", "gtdb_taxonomy": "", "source_database": "UniProtKB;Europe PMC", "source_database_version": "v", "retrieval_date": "2026-09-13", "primary_doi": "10.1/example", "pmid": "1", "pmcid": "", "experimental_substrate": "PHB", "experimental_result": "degraded", "evidence_type": "direct_experiment", "evidence_status": "experimental_positive", "exact_identifier_scope": "protein_accession+genome_accession+strain", "sequence_sha256": "", "notes": ""}])
            self._write_tsv(taxonomy, ["genome_accession", "gtdb_genome_accession", "gtdb_taxonomy", "organism_name", "strain_name"], [{"genome_accession": "NC_008313.1", "gtdb_genome_accession": "RS_GCF_000009285.1", "gtdb_taxonomy": "d__Bacteria;p__Pseudomonadota", "organism_name": "Cupriavidus necator", "strain_name": "H16"}])
            rows = module.bridge_records(literature, taxonomy, output)
            self.assertEqual(rows[0]["match_class"], "exact_genome_accession")
            self.assertEqual(rows[0]["evidence_tier"], "direct_literature_supported")
            self.assertEqual(rows[0]["gtdb_genome_accession"], "RS_GCF_000009285.1")

    def test_experimental_positive_exact_accession_is_direct_bridge(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature = root / "literature.tsv"
            taxonomy = root / "taxonomy.tsv"
            output = root / "bridge.tsv"
            fields = [
                "record_id", "organism_name", "strain_name", "protein_accession",
                "gene_accession", "genome_accession", "ncbi_taxid", "gtdb_taxonomy",
                "source_database", "source_database_version", "retrieval_date",
                "primary_doi", "pmid", "pmcid", "experimental_substrate",
                "experimental_result", "evidence_type", "evidence_status",
                "exact_identifier_scope", "sequence_sha256", "notes",
            ]
            self._write_tsv(literature, fields, [{
                "record_id": "lit_experimental_positive",
                "organism_name": "Cupriavidus necator",
                "strain_name": "H16",
                "protein_accession": "Q0K9H3",
                "gene_accession": "H16_A2251",
                "genome_accession": "NC_008313.1",
                "ncbi_taxid": "381666",
                "gtdb_taxonomy": "",
                "source_database": "UniProtKB;Europe PMC",
                "source_database_version": "v",
                "retrieval_date": "2026-09-13",
                "primary_doi": "10.1/example",
                "pmid": "1",
                "pmcid": "",
                "experimental_substrate": "PHB",
                "experimental_result": "degraded",
                "evidence_type": "direct_experiment",
                "evidence_status": "experimental_positive",
                "exact_identifier_scope": "protein_accession+genome_accession+strain",
                "sequence_sha256": "",
                "notes": "",
            }])
            self._write_tsv(taxonomy, [
                "genome_accession", "gtdb_taxonomy",
            ], [{
                "genome_accession": "NC_008313.1",
                "gtdb_taxonomy": "d__Bacteria;p__Pseudomonadota",
            }])
            rows = module.bridge_records(literature, taxonomy, output)
            self.assertEqual(rows[0]["match_class"], "exact_genome_accession")
            self.assertEqual(rows[0]["evidence_tier"], "direct_literature_supported")

    def test_taxonomy_only_match_cannot_create_tier_a(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature = root / "literature.tsv"
            taxonomy = root / "taxonomy.tsv"
            output = root / "bridge.tsv"
            fields = ["record_id", "organism_name", "strain_name", "protein_accession", "gene_accession", "genome_accession", "ncbi_taxid", "gtdb_taxonomy", "source_database", "source_database_version", "retrieval_date", "primary_doi", "pmid", "pmcid", "experimental_substrate", "experimental_result", "evidence_type", "evidence_status", "exact_identifier_scope", "sequence_sha256", "notes"]
            self._write_tsv(literature, fields, [{"record_id": "lit_2", "organism_name": "Paucimonas lemoignei", "strain_name": "PHA1", "protein_accession": "", "gene_accession": "", "genome_accession": "", "ncbi_taxid": "", "gtdb_taxonomy": "", "source_database": "PubMed", "source_database_version": "2026-09-13", "retrieval_date": "2026-09-13", "primary_doi": "10.1/example", "pmid": "2", "pmcid": "", "experimental_substrate": "PHB", "experimental_result": "degraded", "evidence_type": "primary_experiment", "evidence_status": "direct_literature_supported", "exact_identifier_scope": "strain", "sequence_sha256": "", "notes": ""}])
            self._write_tsv(taxonomy, ["genome_accession", "gtdb_taxonomy", "organism_name", "strain_name"], [{"genome_accession": "GCA_000999.1", "gtdb_taxonomy": "d__Bacteria;p__Example", "organism_name": "Paucimonas lemoignei", "strain_name": "other"}])
            rows = module.bridge_records(literature, taxonomy, output)
            self.assertEqual(rows[0]["match_class"], "taxonomy_only_match")
            self.assertNotEqual(rows[0]["evidence_tier"], "direct_literature_supported")

    def test_unresolved_record_is_retained(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature = root / "literature.tsv"
            taxonomy = root / "taxonomy.tsv"
            output = root / "bridge.tsv"
            fields = ["record_id", "organism_name", "strain_name", "protein_accession", "gene_accession", "genome_accession", "ncbi_taxid", "gtdb_taxonomy", "source_database", "source_database_version", "retrieval_date", "primary_doi", "pmid", "pmcid", "experimental_substrate", "experimental_result", "evidence_type", "evidence_status", "exact_identifier_scope", "sequence_sha256", "notes"]
            self._write_tsv(literature, fields, [{"record_id": "lit_3", "organism_name": "Unknown", "strain_name": "", "protein_accession": "", "gene_accession": "", "genome_accession": "GCA_MISSING", "ncbi_taxid": "", "gtdb_taxonomy": "", "source_database": "PubMed", "source_database_version": "2026-09-13", "retrieval_date": "2026-09-13", "primary_doi": "10.1/example", "pmid": "3", "pmcid": "", "experimental_substrate": "PHB", "experimental_result": "degraded", "evidence_type": "primary_experiment", "evidence_status": "direct_literature_supported", "exact_identifier_scope": "genome", "sequence_sha256": "", "notes": ""}])
            self._write_tsv(taxonomy, ["genome_accession", "gtdb_taxonomy"], [])
            rows = module.bridge_records(literature, taxonomy, output)
            self.assertEqual(rows[0]["match_class"], "unresolved_no_gtdb_record")
            self.assertEqual(rows[0]["evidence_tier"], "hold_unresolved")

    def test_one_literature_record_retains_all_matching_gtdb_rows(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature = root / "literature.tsv"
            taxonomy = root / "taxonomy.tsv"
            output = root / "bridge.tsv"
            fields = ["record_id", "organism_name", "strain_name", "protein_accession", "gene_accession", "genome_accession", "ncbi_taxid", "gtdb_taxonomy", "source_database", "source_database_version", "retrieval_date", "primary_doi", "pmid", "pmcid", "experimental_substrate", "experimental_result", "evidence_type", "evidence_status", "exact_identifier_scope", "sequence_sha256", "notes"]
            self._write_tsv(literature, fields, [{"record_id": "lit_multi", "organism_name": "Example bacterium", "strain_name": "", "protein_accession": "", "gene_accession": "", "genome_accession": "", "ncbi_taxid": "", "gtdb_taxonomy": "", "source_database": "PubMed", "source_database_version": "v", "retrieval_date": "2026-09-13", "primary_doi": "10.1/x", "pmid": "4", "pmcid": "", "experimental_substrate": "PHB", "experimental_result": "degraded", "evidence_type": "primary_experiment", "evidence_status": "direct_literature_supported", "exact_identifier_scope": "strain", "sequence_sha256": "", "notes": ""}])
            self._write_tsv(taxonomy, ["genome_accession", "gtdb_taxonomy", "organism_name", "strain_name"], [{"genome_accession": "GCA_1", "gtdb_taxonomy": "d__Bacteria", "organism_name": "Example bacterium", "strain_name": "A"}, {"genome_accession": "GCA_2", "gtdb_taxonomy": "d__Bacteria", "organism_name": "Example bacterium", "strain_name": "B"}])
            rows = module.bridge_records(literature, taxonomy, output)
            self.assertEqual([r["gtdb_genome_accession"] for r in rows], ["GCA_1", "GCA_2"])
            self.assertTrue(all(r["evidence_tier"] != "direct_literature_supported" for r in rows))

    def test_manifest_binds_inputs_and_output(self):
        module = self._load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            literature = root / "literature.tsv"
            taxonomy = root / "taxonomy.tsv"
            output = root / "bridge.tsv"
            manifest = root / "manifest.json"
            self._write_tsv(literature, ["record_id", "organism_name"], [{"record_id": "lit", "organism_name": "x"}])
            self._write_tsv(taxonomy, ["genome_accession"], [])
            module.bridge_records(literature, taxonomy, output)
            payload = module.write_manifest(manifest, literature, taxonomy, output, "GTDB R232")
            self.assertEqual(payload["gtdb_release"], "GTDB R232")
            self.assertTrue(payload["files"]["literature"]["sha256"])
            self.assertTrue(json.loads(manifest.read_text(encoding="utf-8"))["files"]["output"]["size"] > 0)


if __name__ == "__main__":
    unittest.main()
