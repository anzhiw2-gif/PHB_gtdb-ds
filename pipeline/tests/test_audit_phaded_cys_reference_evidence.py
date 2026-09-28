"""Tests for the Cys reference evidence audit (Task 2026-09-17, Step 2).

Written before the implementation (test-first).  The module under test is a
read-only audit: it joins the three externally-declared accessions
(``O87189`` / ``Q7WT48`` / ``Q7WT49``) to the frozen 723-row reference ledger,
records verbatim retrieved database evidence with retrieval metadata and
hashes, and never writes into the frozen ledger.
"""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_phaded_cys_reference_evidence.py"

CYS = "intracellular nPHASCL without lipase box"

LEDGER_HEADER = [
    "reference_id", "accession", "sequence_sha256", "sequence_source",
    "source_database", "source_database_version", "retrieval_date", "organism",
    "taxonomy_id", "phaded_superfamily", "phaded_family_id",
    "reported_localization", "substrate_class", "substrate_detail",
    "experimental_assay", "experimental_result", "evidence_status",
    "positive_negative_control", "catalytic_residues", "lipase_box_or_ahsmg",
    "oxyanion_hole_evidence", "domain_architecture", "primary_doi", "pmid",
    "pmcid", "notes",
]

RETRIEVAL_HEADER = [
    "accession", "uniprot_entry_name", "protein_name", "gene_name",
    "embl_protein_accession", "sequence_length", "source_pmid", "source_doi",
    "evidence_class", "experimental_positive", "evidence_statement",
    "task5_panel_label", "retrieval_date", "retrieval_url",
    "retrieved_file", "retrieval_http_status",
]


def load_module():
    spec = importlib.util.spec_from_file_location("audit_phaded_cys_reference_evidence", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_text(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8", newline="\n")


def write_tsv(path, header, rows):
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join(str(row.get(col, "")) for col in header))
    write_text(path, "\n".join(lines) + "\n")


UNIPROT_O87189 = """ID   O87189_CUPNE            Unreviewed;       419 AA.
AC   O87189;
DE   SubName: Full=Intracellular PHB depolymerase {ECO:0000313|EMBL:BAA33394.1};
GN   Name=phaZre {ECO:0000313|EMBL:BAA33394.1};
OS   Cupriavidus necator (Alcaligenes eutrophus) (Ralstonia eutropha).
RX   PubMed=11114905; DOI=10.1128/JB.183.1.94-100.2001;
DR   EMBL; AB017612; BAA33394.1; -; Genomic_DNA.
DR   Pfam; PF06850; PHB_depo_C; 1.
SQ   SEQUENCE   419 AA;  47317 MW;
     MLYQLHEFQR SILHPLTAWA
//
"""


class AuditModuleTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    # ---- reference index -------------------------------------------------
    def test_reference_index_splits_reference_id_and_accession(self):
        faa = self.root / "ref.faa"
        write_text(faa, ">DED_hfam_65_0001|CAJ92291.1\nMKT\n>DED_hfam_61_0029|AAP74580.1\nMKT\n")
        by_id, by_acc = self.module.load_reference_index(faa)
        self.assertEqual(sorted(by_id), ["DED_hfam_61_0029", "DED_hfam_65_0001"])
        self.assertEqual(by_acc["CAJ92291.1"], "MKT")

    # ---- UniProt parsing -------------------------------------------------
    def test_parse_uniprot_record_extracts_cross_references(self):
        parsed = self.module.parse_uniprot_record(UNIPROT_O87189)
        self.assertEqual(parsed["accession"], "O87189")
        self.assertEqual(parsed["entry_name"], "O87189_CUPNE")
        self.assertEqual(parsed["sequence_length"], 419)
        self.assertEqual(parsed["embl_protein_accessions"], ["BAA33394.1"])
        self.assertEqual(parsed["pubmed_ids"], ["11114905"])
        self.assertEqual(parsed["doi"], ["10.1128/JB.183.1.94-100.2001"])
        self.assertIn("Intracellular PHB depolymerase", parsed["protein_name"])

    # ---- sequence identity ----------------------------------------------
    def test_exact_match_targets_finds_byte_identical_reference(self):
        seq = "MLYQLHEFQR"
        by_acc = {"CAJ92291.1": "MLYQLHEFQR", "AAP74580.1": "MHIAARKNQN"}
        self.assertEqual(self.module.exact_match_targets(seq, by_acc), ["CAJ92291.1"])

    def test_exact_match_targets_reports_truncation_relations(self):
        by_acc = {"CAJ93939.1": "LYHAYQIYADMILPACTLAEL"}
        relations = self.module.sequence_relations("MHIAARKNQNSGDRAMLYHAYQIYADMILPACTLAEL", by_acc)
        self.assertEqual(relations["CAJ93939.1"], "target_is_c_terminal_subsequence_of_query")

    # ---- audit rows ------------------------------------------------------
    def _fixture(self):
        ledger = self.root / "ledger.tsv"
        write_tsv(ledger, LEDGER_HEADER, [
            {"reference_id": "DED_hfam_65_0001", "accession": "CAJ92291.1",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_65",
             "evidence_status": "experimental_positive", "reported_localization": "intracellular"},
            {"reference_id": "DED_hfam_61_0029", "accession": "AAP74580.1",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_61",
             "evidence_status": "annotation_only", "reported_localization": "intracellular"},
            {"reference_id": "DED_hfam_66_0001", "accession": "AAP74581.1",
             "phaded_superfamily": CYS, "phaded_family_id": "DED_hfam_66",
             "evidence_status": "annotation_only", "reported_localization": "intracellular"},
        ])
        faa = self.root / "ref.faa"
        write_text(faa, ">DED_hfam_65_0001|CAJ92291.1\nAAAA\n"
                        ">DED_hfam_61_0029|AAP74580.1\nCCCC\n"
                        ">DED_hfam_66_0001|AAP74581.1\nDDDD\n")
        retrieval = self.root / "retrieval.tsv"
        write_tsv(retrieval, RETRIEVAL_HEADER, [
            {"accession": "O87189", "sequence_length": 419, "source_pmid": "11114905",
             "evidence_class": "cloned_gene_in_vitro_activity_and_in_vivo_localization",
             "experimental_positive": "true", "embl_protein_accession": "BAA33394.1",
             "retrieval_date": "2026-09-17", "retrieval_http_status": "200"},
            {"accession": "Q7WT49", "sequence_length": 419, "source_pmid": "12813072",
             "evidence_class": "deletion_strain_in_vivo_phenotype",
             "experimental_positive": "true", "embl_protein_accession": "AAP74580.1",
             "retrieval_date": "2026-09-17", "retrieval_http_status": "200"},
            {"accession": "Q7WT48", "sequence_length": 407, "source_pmid": "12813072",
             "evidence_class": "gene_identified_role_not_established",
             "experimental_positive": "false", "embl_protein_accession": "AAP74581.1",
             "retrieval_date": "2026-09-17", "retrieval_http_status": "200"},
        ])
        # sequences that make O87189==CAJ92291.1 and Q7WT49==AAP74580.1
        by_acc = {"CAJ92291.1": "AAAA", "AAP74580.1": "CCCC", "AAP74581.1": "DDDD"}
        records = {
            "O87189": dict(self.module.parse_uniprot_record(UNIPROT_O87189),
                           sequence="AAAA", retrieved_file="uniprot_O87189.txt",
                           retrieved_file_sha256="deadbeef", retrieval_http_status="200"),
            "Q7WT49": {"accession": "Q7WT49", "entry_name": "Q7WT49_CUPNE",
                       "protein_name": "Intracellular PHB depolymerase PhaZ2",
                       "sequence_length": 419, "embl_protein_accessions": ["AAP74580.1"],
                       "pubmed_ids": ["12813072"], "doi": [], "sequence": "CCCC",
                       "retrieved_file": "uniprot_Q7WT49.txt",
                       "retrieved_file_sha256": "cafebabe", "retrieval_http_status": "200"},
            "Q7WT48": {"accession": "Q7WT48", "entry_name": "Q7WT48_CUPNE",
                       "protein_name": "Intracellular PHB depolymerase PhaZ3",
                       "sequence_length": 407, "embl_protein_accessions": ["AAP74581.1"],
                       "pubmed_ids": ["12813072"], "doi": [], "sequence": "DDDD",
                       "retrieved_file": "uniprot_Q7WT48.txt",
                       "retrieved_file_sha256": "feedface", "retrieval_http_status": "200"},
        }
        return ledger, faa, retrieval, by_acc, records

    def test_build_audit_rows_marks_duplicate_of_frozen_positive(self):
        ledger, faa, retrieval, by_acc, records = self._fixture()
        rows = self.module.build_audit_rows(
            retrieval_rows=self.module.read_tsv(retrieval),
            uniprot_records=records,
            reference_by_accession=by_acc,
            ledger_rows=self.module.read_tsv(ledger),
        )
        by_accession = {r["accession"]: r for r in rows}
        o87189 = by_accession["O87189"]
        self.assertEqual(o87189["sequence_identical_to_ledger_accession"], "CAJ92291.1")
        self.assertEqual(o87189["ledger_reference_id"], "DED_hfam_65_0001")
        self.assertEqual(o87189["ledger_evidence_status"], "experimental_positive")
        self.assertEqual(o87189["is_alias_of_existing_ledger_record"], "true")
        self.assertEqual(o87189["experimental_positive"], "true")
        self.assertEqual(o87189["source_pmid"], "11114905")

        q7wt49 = by_accession["Q7WT49"]
        self.assertEqual(q7wt49["ledger_reference_id"], "DED_hfam_61_0029")
        self.assertEqual(q7wt49["ledger_evidence_status"], "annotation_only")
        self.assertEqual(q7wt49["experimental_positive"], "true")
        self.assertEqual(q7wt49["is_alias_of_existing_ledger_record"], "true")

        q7wt48 = by_accession["Q7WT48"]
        self.assertEqual(q7wt48["ledger_reference_id"], "DED_hfam_66_0001")
        self.assertEqual(q7wt48["experimental_positive"], "false")

    def test_every_row_records_retrieval_provenance(self):
        ledger, faa, retrieval, by_acc, records = self._fixture()
        rows = self.module.build_audit_rows(
            retrieval_rows=self.module.read_tsv(retrieval),
            uniprot_records=records,
            reference_by_accession=by_acc,
            ledger_rows=self.module.read_tsv(ledger),
        )
        for row in rows:
            self.assertTrue(row["retrieval_date"])
            self.assertTrue(row["retrieved_file_sha256"])
            self.assertEqual(row["retrieval_http_status"], "200")
            self.assertIn("does_not_modify_frozen_ledger", row["boundary_note"])

    # ---- per-run evidence statistics ------------------------------------
    def test_summary_counts_distinct_experimental_positive_proteins(self):
        ledger, faa, retrieval, by_acc, records = self._fixture()
        rows = self.module.build_audit_rows(
            retrieval_rows=self.module.read_tsv(retrieval),
            uniprot_records=records,
            reference_by_accession=by_acc,
            ledger_rows=self.module.read_tsv(ledger),
        )
        summary = self.module.summarize_evidence(rows, self.module.read_tsv(ledger))
        self.assertEqual(summary["frozen_ledger_cys_reference_count"], 3)
        self.assertEqual(summary["frozen_ledger_cys_experimental_positive_count"], 1)
        self.assertEqual(summary["audited_accessions"], ["O87189", "Q7WT48", "Q7WT49"])
        # O87189 is a byte-identical alias of the frozen positive CAJ92291.1,
        # while Q7WT49 maps to AAP74580.1 (annotation_only in the frozen ledger
        # but experimentally supported in the literature).  Two distinct
        # ledger records therefore carry experimental support, and the frozen
        # ledger itself is untouched.
        self.assertEqual(summary["distinct_experimental_positive_cys_records"], 2)
        self.assertEqual(summary["experimental_positive_records_added_by_this_audit"],
                         ["DED_hfam_61_0029"])
        self.assertEqual(summary["gate_positive_requirement"], 3)
        self.assertEqual(summary["gate_positive_requirement_met"], False)
        self.assertEqual(summary["frozen_ledger_modified"], False)

    def test_cli_writes_outputs_into_run_dir_only(self):
        ledger, faa, retrieval, by_acc, records = self._fixture()
        run_dir = self.root / "runs" / "20260917_phaded_cys_discovery_hmm_01"
        (run_dir / "inputs" / "external_evidence").mkdir(parents=True)
        rows = self.module.build_audit_rows(
            retrieval_rows=self.module.read_tsv(retrieval),
            uniprot_records=records,
            reference_by_accession=by_acc,
            ledger_rows=self.module.read_tsv(ledger),
        )
        out_tsv = run_dir / "results" / "reference_evidence_audit.tsv"
        out_json = run_dir / "results" / "cys_evidence_statistics_per_run.json"
        summary = self.module.summarize_evidence(rows, self.module.read_tsv(ledger))
        self.module.write_audit_outputs(out_tsv, out_json, rows, summary)
        self.assertTrue(out_tsv.is_file())
        payload = json.loads(out_json.read_text(encoding="utf-8"))
        self.assertEqual(payload["frozen_ledger_modified"], False)
        self.assertEqual(payload["boundary"], "candidate_only")


if __name__ == "__main__":
    unittest.main()
