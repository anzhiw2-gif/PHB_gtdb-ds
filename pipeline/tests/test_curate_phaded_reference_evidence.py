"""Tests for the F13 reference-evidence curation script (synthetic fixtures only).

The curation script turns the frozen 723-row PhaDED reference ledger plus the
30-row primary-source amendment into a *new* ledger that carries the 15 v2
evidence columns.  Every fixture here is synthetic and lives in a temporary
directory: the frozen ``runs/`` trees are neither read nor written by this module.

The reference ids used by the fixtures are real DED reference ids *because* the
curation table of the script is keyed by them; each fixture copies the documented
facts of that row (seed GI, primary paper, evidence type) from the frozen ledger
and the amendment.  The values a test asserts on (grouping, counting, grades,
pending sentinels) are then produced by the real code paths.

The invariants pinned here are the ones the design spec makes non-negotiable:

* append-only — no original column value of the frozen ledger may change;
* no fabricated provenance — a positive without an amendment row is refused;
* grade ``A`` can never be ``functional_calibration_eligible=true``;
* the functional-positive count is the number of unique ``independence_group``
  values, never the row count;
* ``discovery_training_eligible=true`` requires ``sequence_integrity=complete``;
* a missing input or a non-empty output directory is refused, not overwritten.
"""

from __future__ import annotations

import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "pipeline" / "scripts" / "curate_phaded_reference_evidence.py"
VALIDATOR_SCRIPT = REPO / "pipeline" / "scripts" / "validate_phaded_reference_ledger.py"
EVIDENCE_COLUMNS_CONFIG = REPO / "pipeline" / "config" / "phaded_reference_evidence_columns.tsv"

# The frozen 26-column ledger contract, in file order.
FROZEN_FIELDS = [
    "reference_id", "accession", "sequence_sha256", "sequence_source", "source_database",
    "source_database_version", "retrieval_date", "organism", "taxonomy_id",
    "phaded_superfamily", "phaded_family_id", "reported_localization", "substrate_class",
    "substrate_detail", "experimental_assay", "experimental_result", "evidence_status",
    "positive_negative_control", "catalytic_residues", "lipase_box_or_ahsmg",
    "oxyanion_hole_evidence", "domain_architecture", "primary_doi", "pmid", "pmcid", "notes",
]

# The v2 evidence columns, in the order declared by pipeline/config/phaded_reference_evidence_columns.tsv.
EVIDENCE_FIELDS = [
    "accession_version", "gene_id", "study_id", "experiment_unit_id", "independence_group",
    "experimental_evidence_grade", "assay_directness", "experimental_substrate_class",
    "experimental_system", "negative_control_description", "catalytic_mutant_evidence",
    "localization_evidence", "sequence_integrity", "discovery_training_eligible",
    "functional_calibration_eligible",
]

PROVENANCE_FIELDS = [
    "reference_id", "accession", "organism", "phaded_superfamily", "phaded_family_id",
    "reported_localization", "substrate_class", "seed_gi", "ledger_pmid", "ledger_doi",
    "primary_pmid", "primary_doi", "primary_title", "evidence_type", "confidence",
    "source_url", "note",
]

SUPERFAMILY_EXTRACELLULAR_1 = "extracellular dPHASCL type 1"
SUPERFAMILY_INTRACELLULAR = "intracellular nPHAMCL"

FAMILIES = [
    (SUPERFAMILY_EXTRACELLULAR_1, "DED_hfam_52"),
    (SUPERFAMILY_EXTRACELLULAR_1, "DED_hfam_53"),
    (SUPERFAMILY_INTRACELLULAR, "DED_hfam_4"),
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sequence_for(accession: str, length: int = 48, *, leading_m: bool = True, ambiguous: bool = False) -> str:
    """A deterministic synthetic protein; distinct accessions always give distinct sequences."""
    digest = hashlib.sha256(accession.encode("utf-8")).hexdigest()
    body = "".join(
        "ACDEFGHIKLMNPQRSTVWY"[int(digest[index] + digest[index + 1], 16) % 20]
        for index in range(length - 1)
    )
    if ambiguous:
        body = body[:-1] + "X"
    return ("M" if leading_m else "A") + body


def ledger_row(
    reference_id: str,
    accession: str,
    family_id: str,
    status: str,
    sequence: str,
    *,
    superfamily: str = SUPERFAMILY_EXTRACELLULAR_1,
    organism: str = "Escherichia coli",
    seed_gi: str | None = None,
    notes: str | None = None,
    substrate_class: str = "SCL",
    substrate_detail: str = "denatured PHA",
) -> dict:
    """One synthetic frozen-ledger row that satisfies the frozen ledger contract."""
    positive = status in {"experimental_positive", "experimental_negative"}
    if notes is None:
        notes = (
            f"DED family page http://fixture.invalid/hfam; Knoll 2009 Table 1 seed GI {seed_gi}"
            if positive and seed_gi
            else "DED family page http://fixture.invalid/hfam; fixture annotation record"
        )
    return {
        "reference_id": reference_id,
        "accession": accession,
        "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        "sequence_source": "http://fixture.invalid/seq.txt",
        "source_database": "PhaDED/DED",
        "source_database_version": "fixture snapshot",
        "retrieval_date": "2026-09-28",
        "organism": organism,
        "taxonomy_id": "562",
        "phaded_superfamily": superfamily,
        "phaded_family_id": family_id,
        "reported_localization": "extracellular",
        "substrate_class": substrate_class,
        "substrate_detail": substrate_detail,
        "experimental_assay": "fixture assay" if positive else "pending",
        "experimental_result": "fixture depolymerase activity" if positive else "pending",
        "evidence_status": status,
        "positive_negative_control": "positive_seed" if positive else "not_assessed",
        "catalytic_residues": "pending",
        "lipase_box_or_ahsmg": "pending",
        "oxyanion_hole_evidence": "pending",
        "domain_architecture": "pending",
        "primary_doi": "10.1186/1471-2105-10-89",
        "pmid": "19296857",
        "pmcid": "PMC2666664",
        "notes": notes,
    }


def provenance_row(
    reference_id: str,
    accession: str,
    family_id: str,
    seed_gi: str,
    *,
    primary_pmid: str = "",
    primary_doi: str = "",
    evidence_type: str = "cloned_activity",
    confidence: str = "high",
    note: str = "fixture note",
    primary_title: str = "fixture primary characterization",
    superfamily: str = SUPERFAMILY_EXTRACELLULAR_1,
    organism: str = "Escherichia coli",
) -> dict:
    return {
        "reference_id": reference_id,
        "accession": accession,
        "organism": organism,
        "phaded_superfamily": superfamily,
        "phaded_family_id": family_id,
        "reported_localization": "extracellular",
        "substrate_class": "SCL",
        "seed_gi": seed_gi,
        "ledger_pmid": "19296857",
        "ledger_doi": "10.1186/1471-2105-10-89",
        "primary_pmid": primary_pmid,
        "primary_doi": primary_doi,
        "primary_title": primary_title,
        "evidence_type": evidence_type,
        "confidence": confidence,
        "source_url": (
            "https://pubmed.ncbi.nlm.nih.gov/%s/" % primary_pmid
            if primary_pmid
            else "https://fixture.invalid/genbank"
        ),
        "note": note,
    }


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


class Fixture:
    """A synthetic ledger + provenance amendment + FASTA + family definitions."""

    def __init__(
        self,
        root: Path,
        ledger_rows: list[dict],
        provenance_rows: list[dict],
        sequences: dict[str, str],
        families: list[tuple[str, str]] | None = None,
    ) -> None:
        self.root = root
        self.ledger_rows = ledger_rows
        self.provenance_rows = provenance_rows
        self.sequences = sequences
        self.families = families or list(FAMILIES)
        self.ledger = root / "phaded_reference_ledger.tsv"
        self.provenance = root / "positive_primary_source_amendment.tsv"
        self.fasta = root / "phaded_reference.faa"
        self.family_definitions = root / "phaded_family_definitions.tsv"
        self.write_ledger()
        self.write_provenance()
        write_tsv(
            self.family_definitions,
            ["phaded_superfamily", "phaded_family_id"],
            [{"phaded_superfamily": sf, "phaded_family_id": fam} for sf, fam in self.families],
        )
        self.write_fasta()

    def write_ledger(self) -> None:
        write_tsv(self.ledger, FROZEN_FIELDS, self.ledger_rows)

    def write_provenance(self) -> None:
        write_tsv(self.provenance, PROVENANCE_FIELDS, self.provenance_rows)

    def write_fasta(self) -> None:
        lines: list[str] = []
        for row in self.ledger_rows:
            lines.append(">" + row["reference_id"] + "|" + row["accession"])
            lines.append(self.sequences[row["reference_id"]])
        self.fasta.write_text("\n".join(lines) + "\n", encoding="ascii")

    def kwargs(self, output_dir: Path, **overrides) -> dict:
        kwargs = {
            "ledger": self.ledger,
            "provenance": self.provenance,
            "family_definitions": self.family_definitions,
            "reference_fasta": self.fasta,
            "output_dir": output_dir,
            "evidence_columns": EVIDENCE_COLUMNS_CONFIG,
            "expected_family_count": len(self.families),
            "write_report": False,
        }
        kwargs.update(overrides)
        return kwargs


def same_paper_fixture(root: Path) -> Fixture:
    """DED_hfam_52_0004 + DED_hfam_52_0007: two genes, one primary paper (PMID 7836292)."""
    accessions = {"DED_hfam_52_0004": "AAA65703.1", "DED_hfam_52_0007": "AAA65705.1"}
    sequences = {reference: sequence_for(accession) for reference, accession in accessions.items()}
    ledger_rows = [
        ledger_row(reference, accession, "DED_hfam_52", "experimental_positive", sequences[reference],
                   seed_gi=seed_gi, organism="Paucimonas lemoignei")
        for reference, accession, seed_gi in (
            ("DED_hfam_52_0004", "AAA65703.1", "531464"),
            ("DED_hfam_52_0007", "AAA65705.1", "531466"),
        )
    ]
    provenance_rows = [
        provenance_row(reference, accession, "DED_hfam_52", seed_gi, primary_pmid="7836292",
                       organism="Paucimonas lemoignei")
        for reference, accession, seed_gi in (
            ("DED_hfam_52_0004", "AAA65703.1", "531464"),
            ("DED_hfam_52_0007", "AAA65705.1", "531466"),
        )
    ]
    return Fixture(root, ledger_rows, provenance_rows, sequences)


def different_paper_fixture(root: Path) -> Fixture:
    """DED_hfam_52_0001 + DED_hfam_52_0002: two genes, two different primary papers."""
    sequences = {
        "DED_hfam_52_0001": sequence_for("P12625.1"),
        "DED_hfam_52_0002": sequence_for("AAB40611.1"),
    }
    ledger_rows = [
        ledger_row("DED_hfam_52_0001", "P12625.1", "DED_hfam_52", "experimental_positive",
                   sequences["DED_hfam_52_0001"], seed_gi="130019", organism="Ralstonia pickettii"),
        ledger_row("DED_hfam_52_0002", "AAB40611.1", "DED_hfam_52", "experimental_positive",
                   sequences["DED_hfam_52_0002"], seed_gi="1777951", organism="Alcaligenes faecalis"),
    ]
    provenance_rows = [
        provenance_row("DED_hfam_52_0001", "P12625.1", "DED_hfam_52", "130019",
                       primary_pmid="2644188", organism="Ralstonia pickettii"),
        provenance_row("DED_hfam_52_0002", "AAB40611.1", "DED_hfam_52", "1777951",
                       primary_pmid="9177489", organism="Alcaligenes faecalis"),
    ]
    return Fixture(root, ledger_rows, provenance_rows, sequences)


def duplicate_gene_fixture(root: Path) -> Fixture:
    """DED_hfam_53_0001 (characterized) + DED_hfam_53_0008 (same gene duplicate) + one annotation row."""
    sequences = {
        "DED_hfam_53_0001": sequence_for("O82950"),
        "DED_hfam_53_0008": sequence_for("ACG63775.1"),
        "DED_hfam_4_0009": sequence_for("Q99999.1"),
    }
    ledger_rows = [
        ledger_row("DED_hfam_53_0001", "O82950", "DED_hfam_53", "experimental_positive",
                   sequences["DED_hfam_53_0001"], seed_gi="75538924",
                   organism="Pseudomonas stutzeri"),
        ledger_row("DED_hfam_53_0008", "ACG63775.1", "DED_hfam_53", "experimental_positive",
                   sequences["DED_hfam_53_0008"], seed_gi="75538924",
                   organism="Pseudomonas stutzeri"),
        ledger_row("DED_hfam_4_0009", "Q99999.1", "DED_hfam_4", "annotation_only",
                   sequences["DED_hfam_4_0009"], superfamily=SUPERFAMILY_INTRACELLULAR),
    ]
    provenance_rows = [
        provenance_row("DED_hfam_53_0001", "O82950", "DED_hfam_53", "75538924",
                       primary_pmid="9872779", organism="Pseudomonas stutzeri"),
        provenance_row("DED_hfam_53_0008", "ACG63775.1", "DED_hfam_53", "75538924",
                       evidence_type="same_gene_duplicate", confidence="none",
                       primary_title="(Direct submission; P. stutzeri PHB depolymerase)",
                       note="GenBank re-sequencing of the P. stutzeri PHB depolymerase (same gene "
                            "as O82950). No independent paper. NOT an independent positive.",
                       organism="Pseudomonas stutzeri"),
    ]
    return Fixture(root, ledger_rows, provenance_rows, sequences)


class CurationTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.module = load_module(SCRIPT, "curate_phaded_reference_evidence")
        self.validator = load_module(VALIDATOR_SCRIPT, "validate_phaded_reference_ledger")
        self._temporary = tempfile.TemporaryDirectory()
        self.root = Path(self._temporary.name)
        self.addCleanup(self._temporary.cleanup)
        self.output_dir = self.root / "results"

    def read_curated(self, output_dir: Path | None = None):
        path = (output_dir or self.output_dir) / "phaded_reference_ledger_curated_v2.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return list(reader.fieldnames or []), list(reader)

    def summary(self, output_dir: Path | None = None) -> dict:
        path = (output_dir or self.output_dir) / "reference_evidence_independence_summary.json"
        return json.loads(path.read_text(encoding="utf-8"))

    def curated_by_reference(self, payload: dict) -> dict[str, dict]:
        return {row["reference_id"]: row for row in payload["curated_rows"]}


class AppendOnlyTests(CurationTestCase):
    def test_curated_ledger_preserves_every_original_column_and_value(self):
        fixture = same_paper_fixture(self.root)
        self.module.run(**fixture.kwargs(self.output_dir))
        fields, curated = self.read_curated()
        self.assertEqual(fields[: len(FROZEN_FIELDS)], FROZEN_FIELDS)
        self.assertEqual(fields[len(FROZEN_FIELDS):], EVIDENCE_FIELDS)
        self.assertEqual(len(curated), len(fixture.ledger_rows))
        for original, new in zip(fixture.ledger_rows, curated):
            for field in FROZEN_FIELDS:
                self.assertEqual(new[field], original[field], field)
        # The frozen bytes are a literal prefix of every curated line.
        original_lines = fixture.ledger.read_text(encoding="utf-8").splitlines()
        curated_lines = (self.output_dir / "phaded_reference_ledger_curated_v2.tsv").read_text(
            encoding="utf-8"
        ).splitlines()
        self.assertEqual(curated_lines[0], original_lines[0] + "\t" + "\t".join(EVIDENCE_FIELDS))
        for original_line, curated_line in zip(original_lines[1:], curated_lines[1:]):
            self.assertTrue(curated_line.startswith(original_line + "\t"))
            self.assertEqual(len(curated_line.split("\t")), len(FROZEN_FIELDS) + len(EVIDENCE_FIELDS))

    def test_append_only_guard_rejects_a_changed_original_value(self):
        fixture = same_paper_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        curated = [dict(row) for row in payload["curated_rows"]]
        curated[0]["organism"] = "Tampered organism"
        with self.assertRaisesRegex(ValueError, "append-only"):
            self.module.assert_append_only(FROZEN_FIELDS, fixture.ledger_rows, curated)

    def test_append_only_guard_rejects_a_whitespace_only_change(self):
        fixture = same_paper_fixture(self.root)
        curated = [dict(row) for row in fixture.ledger_rows]
        curated[0]["notes"] = curated[0]["notes"] + " "
        with self.assertRaisesRegex(ValueError, "append-only"):
            self.module.assert_append_only(FROZEN_FIELDS, fixture.ledger_rows, curated)


class ProvenanceTests(CurationTestCase):
    def test_positive_without_provenance_row_is_refused(self):
        fixture = same_paper_fixture(self.root)
        write_tsv(fixture.provenance, PROVENANCE_FIELDS, fixture.provenance_rows[:1])
        with self.assertRaisesRegex(ValueError, "DED_hfam_52_0007"):
            self.module.run(**fixture.kwargs(self.output_dir))

    def test_provenance_row_without_ledger_positive_is_refused(self):
        fixture = same_paper_fixture(self.root)
        extra = provenance_row("DED_hfam_52_0999", "ZZZ99999.1", "DED_hfam_52", "1",
                               primary_pmid="1111111")
        write_tsv(fixture.provenance, PROVENANCE_FIELDS, fixture.provenance_rows + [extra])
        with self.assertRaisesRegex(ValueError, "DED_hfam_52_0999"):
            self.module.run(**fixture.kwargs(self.output_dir))

    def test_seed_gi_disagreement_between_ledger_note_and_amendment_is_refused(self):
        fixture = same_paper_fixture(self.root)
        fixture.ledger_rows[0]["notes"] = fixture.ledger_rows[0]["notes"].replace("531464", "999999")
        fixture.write_ledger()
        with self.assertRaisesRegex(ValueError, "seed GI"):
            self.module.run(**fixture.kwargs(self.output_dir))

    def test_study_anchor_disagreeing_with_the_amendment_is_refused(self):
        fixture = same_paper_fixture(self.root)
        fixture.provenance_rows[0]["primary_pmid"] = "1111111"
        fixture.write_provenance()
        with self.assertRaisesRegex(ValueError, "study anchor disagrees"):
            self.module.run(**fixture.kwargs(self.output_dir))

    def test_amendment_family_call_disagreeing_with_the_ledger_is_refused(self):
        fixture = same_paper_fixture(self.root)
        fixture.provenance_rows[0]["phaded_family_id"] = "DED_hfam_99"
        fixture.write_provenance()
        with self.assertRaisesRegex(ValueError, "phaded_family_id disagrees"):
            self.module.run(**fixture.kwargs(self.output_dir))


class GradeTests(CurationTestCase):
    def test_grade_a_rows_can_never_be_functional_calibration_eligible(self):
        fixture = duplicate_gene_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        by_reference = self.curated_by_reference(payload)
        duplicate = by_reference["DED_hfam_53_0008"]
        self.assertEqual(duplicate["experimental_evidence_grade"], "A")
        self.assertEqual(duplicate["assay_directness"], "annotation")
        self.assertEqual(duplicate["functional_calibration_eligible"], "false")
        annotation = by_reference["DED_hfam_4_0009"]
        self.assertEqual(annotation["experimental_evidence_grade"], "A")
        self.assertEqual(annotation["functional_calibration_eligible"], "false")
        # The validator is the second, independent guard for the same invariant.
        row = dict(annotation)
        row["functional_calibration_eligible"] = "true"
        with self.assertRaisesRegex(ValueError, "functional_calibration_eligible"):
            self.validator.validate_reference_row(row)

    def test_documented_evidence_types_map_to_documented_grades(self):
        fixture = duplicate_gene_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        by_reference = self.curated_by_reference(payload)
        self.assertEqual(by_reference["DED_hfam_53_0001"]["experimental_evidence_grade"], "E2")
        self.assertEqual(by_reference["DED_hfam_53_0001"]["assay_directness"], "direct")
        self.assertEqual(by_reference["DED_hfam_53_0001"]["functional_calibration_eligible"], "true")
        self.assertEqual(by_reference["DED_hfam_53_0008"]["experimental_evidence_grade"], "A")
        self.assertEqual(self.module.grade_for_evidence_type("purified_activity"), "E3")
        self.assertEqual(self.module.grade_for_evidence_type("cloned_activity"), "E2")
        self.assertEqual(self.module.grade_for_evidence_type("sequence_only"), "A")
        with self.assertRaisesRegex(ValueError, "undocumented amendment evidence_type"):
            self.module.grade_for_evidence_type("in_silico_hunch")

    def test_annotation_only_rows_are_grade_a_with_explicit_pending_units(self):
        fixture = duplicate_gene_fixture(self.root)
        self.module.run(**fixture.kwargs(self.output_dir))
        _fields, curated = self.read_curated()
        annotation = next(row for row in curated if row["reference_id"] == "DED_hfam_4_0009")
        self.assertEqual(annotation["experimental_evidence_grade"], "A")
        self.assertEqual(annotation["assay_directness"], "annotation")
        self.assertEqual(annotation["functional_calibration_eligible"], "false")
        for field in ("study_id", "experiment_unit_id", "independence_group",
                      "experimental_substrate_class", "experimental_system", "gene_id"):
            self.assertEqual(annotation[field], "pending", field)
        classes = {entry["column"]: entry for entry in self.summary()["pending_value_classes"]}
        self.assertIn("no experimental study", classes["study_id"]["reason"])
        self.assertEqual(classes["study_id"]["count"], 1)
        self.assertEqual(classes["study_id"]["row_class"], "annotation_only")

    def test_positives_carry_no_pending_evidence_value(self):
        fixture = duplicate_gene_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        self.assertEqual(payload["positives_with_pending_values"], [])
        for row in payload["curated_rows"]:
            if row["evidence_status"] != "experimental_positive":
                continue
            for field in EVIDENCE_FIELDS:
                self.assertNotEqual(row[field].strip().lower(), "pending", field)

    def test_substrate_class_is_composed_from_the_frozen_substrate_fields(self):
        fixture = duplicate_gene_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        by_reference = self.curated_by_reference(payload)
        self.assertEqual(
            by_reference["DED_hfam_53_0001"]["experimental_substrate_class"], "SCL_denatured_PHA"
        )
        self.assertIn(
            "inherited DED label",
            next(
                row["experimental_substrate_class_basis"]
                for row in payload["derivation_rows"]
                if row["reference_id"] == "DED_hfam_53_0008"
            ),
        )


class IndependenceTests(CurationTestCase):
    def test_two_rows_sharing_an_independence_group_count_as_one(self):
        fixture = same_paper_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        summary = payload["independence_summary"]
        self.assertEqual(summary["functional_positive_row_count"], 2)
        self.assertEqual(summary["functional_positive_count"], 1)
        self.assertEqual(summary["unique_independence_group_count"], 1)
        groups = {row["independence_group"] for row in payload["curated_rows"]
                  if row["evidence_status"] == "experimental_positive"}
        self.assertEqual(groups, {"IG_PMID_7836292"})
        units = {row["experiment_unit_id"] for row in payload["curated_rows"]
                 if row["evidence_status"] == "experimental_positive"}
        self.assertEqual(len(units), 2, "different genes in one paper need distinct experiment units")

    def test_two_rows_with_different_groups_count_as_two(self):
        fixture = different_paper_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        summary = payload["independence_summary"]
        self.assertEqual(summary["functional_positive_row_count"], 2)
        self.assertEqual(summary["functional_positive_count"], 2)
        self.assertEqual(
            summary["functional_positive_independence_groups"], ["IG_PMID_2644188", "IG_PMID_9177489"]
        )

    def test_same_gene_rows_share_one_group_and_count_once(self):
        fixture = duplicate_gene_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        curated = self.curated_by_reference(payload)
        group = curated["DED_hfam_53_0001"]["independence_group"]
        self.assertEqual(curated["DED_hfam_53_0008"]["independence_group"], group)
        self.assertEqual(group, "IG_GI_75538924")
        self.assertEqual(
            curated["DED_hfam_53_0008"]["experiment_unit_id"],
            curated["DED_hfam_53_0001"]["experiment_unit_id"],
        )
        summary = payload["independence_summary"]
        # Two rows in the group (a duplicate record plus the characterized gene) but only the
        # characterized row can be counted as a functional positive.
        self.assertEqual(summary["counts_by_independence_group"][group], 2)
        self.assertEqual(summary["functional_positive_by_group"][group], 1)
        self.assertEqual(summary["functional_positive_count"], 1)

    def test_group_derivation_is_recorded_per_positive_row(self):
        fixture = same_paper_fixture(self.root)
        derivation = self.root / "positive_curation_derivation.tsv"
        self.module.run(**fixture.kwargs(self.output_dir, derivation_table=derivation))
        with derivation.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        for row in rows:
            self.assertTrue(row["independence_basis"])
            self.assertTrue(row["independence_group"])
            self.assertTrue(row["grade_basis"])
            self.assertTrue(row["study_id_basis"])
            self.assertEqual(row["seed_gi_source"], "ledger_note_and_amendment")
            self.assertEqual(row["independence_binding"], "same_primary_study")
            self.assertEqual(row["independence_members"].count("DED_hfam_52_"), 2)
        self.assertIn(
            "same_study_key=PMID:7836292",
            rows[0]["independence_documented_bindings"],
        )


class SequenceIntegrityTests(CurationTestCase):
    def test_discovery_training_requires_a_complete_sequence(self):
        module = self.module
        self.assertEqual(module.sequence_integrity(sequence_for("A", leading_m=True)), "complete")
        self.assertEqual(module.sequence_integrity(sequence_for("A", leading_m=False)), "partial")
        self.assertEqual(module.sequence_integrity(sequence_for("A", ambiguous=True)), "unresolved")
        self.assertEqual(module.sequence_integrity("MSTAG*"), "complete")
        self.assertEqual(module.sequence_integrity("MST*AG"), "unresolved")
        self.assertEqual(module.sequence_integrity(""), "unresolved")

    def test_discovery_training_eligible_column_follows_sequence_integrity(self):
        fixture = duplicate_gene_fixture(self.root)
        fixture.sequences["DED_hfam_4_0009"] = sequence_for("Q99999.1", leading_m=False)
        fixture.ledger_rows[2]["sequence_sha256"] = hashlib.sha256(
            fixture.sequences["DED_hfam_4_0009"].encode("ascii")
        ).hexdigest()
        fixture.write_ledger()
        fixture.write_fasta()
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        curated = self.curated_by_reference(payload)
        self.assertEqual(curated["DED_hfam_4_0009"]["sequence_integrity"], "partial")
        self.assertEqual(curated["DED_hfam_4_0009"]["discovery_training_eligible"], "false")
        self.assertEqual(curated["DED_hfam_53_0001"]["sequence_integrity"], "complete")
        self.assertEqual(curated["DED_hfam_53_0001"]["discovery_training_eligible"], "true")
        self.assertFalse(self.module.family_call_consistent(
            {"phaded_superfamily": "not a superfamily", "phaded_family_id": "DED_hfam_53"}
        ))
        self.assertFalse(self.module.family_call_consistent(
            {"phaded_superfamily": SUPERFAMILY_EXTRACELLULAR_1, "phaded_family_id": "unresolved"}
        ))
        self.assertTrue(self.module.family_call_consistent(
            {"phaded_superfamily": SUPERFAMILY_EXTRACELLULAR_1, "phaded_family_id": "DED_hfam_53"}
        ))


class ContractTests(CurationTestCase):
    def test_missing_input_is_refused(self):
        fixture = same_paper_fixture(self.root)
        for name in ("ledger", "provenance", "family_definitions", "reference_fasta"):
            with self.subTest(name=name):
                kwargs = fixture.kwargs(self.output_dir)
                kwargs[name] = self.root / f"missing_{name}.tsv"
                with self.assertRaisesRegex(ValueError, "missing_"):
                    self.module.run(**kwargs)

    def test_missing_evidence_columns_config_is_refused(self):
        fixture = same_paper_fixture(self.root)
        kwargs = fixture.kwargs(self.output_dir, evidence_columns=self.root / "missing_columns.tsv")
        with self.assertRaisesRegex(ValueError, "evidence columns"):
            self.module.run(**kwargs)

    def test_non_empty_output_directory_is_refused(self):
        fixture = same_paper_fixture(self.root)
        self.output_dir.mkdir(parents=True)
        (self.output_dir / "existing.txt").write_text("do not overwrite\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "not empty"):
            self.module.run(**fixture.kwargs(self.output_dir))
        self.assertEqual(
            (self.output_dir / "existing.txt").read_text(encoding="utf-8"), "do not overwrite\n"
        )

    def test_main_returns_non_zero_for_a_refused_run(self):
        fixture = same_paper_fixture(self.root)
        argv = [
            "--ledger", str(fixture.ledger),
            "--provenance", str(fixture.provenance),
            "--family-definitions", str(fixture.family_definitions),
            "--reference-fasta", str(fixture.fasta),
            "--output-dir", str(self.output_dir),
            "--evidence-columns", str(EVIDENCE_COLUMNS_CONFIG),
            "--expected-family-count", str(len(fixture.families)),
            "--no-report",
        ]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(self.module.main(argv), 0)
            # Second run into the same directory must refuse instead of overwriting.
            self.assertEqual(self.module.main(argv), 2)

    def test_strict_validator_gate_passes_on_a_curated_fixture(self):
        fixture = duplicate_gene_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        report = self.validator.validate_ledger(
            self.output_dir / "phaded_reference_ledger_curated_v2.tsv",
            fixture.family_definitions,
            fixture.fasta,
            evidence_column_mode="strict",
            expected_family_count=len(fixture.families),
        )
        self.assertEqual(report["status"], "valid")
        self.assertEqual(report["pending_evidence_columns"], [])
        self.assertEqual(report["functional_positive_count"], 1)
        self.assertEqual(payload["strict_gate"]["status"], "valid")
        self.assertEqual(
            payload["family_functional_positive_counts"],
            report["functional_calibration_gate"]["family_functional_positive_counts"],
        )

    def test_summary_reports_per_family_positives_and_distinct_genera(self):
        fixture = duplicate_gene_fixture(self.root)
        payload = self.module.run(**fixture.kwargs(self.output_dir))
        # The validator reports every family of the ledger, including families with zero
        # curated independent positives, so a zero is visible instead of omitted.
        self.assertEqual(
            payload["family_functional_positive_counts"], {"DED_hfam_4": 0, "DED_hfam_53": 1}
        )
        self.assertEqual(
            payload["family_distinct_genus_counts"], {"DED_hfam_4": 0, "DED_hfam_53": 1}
        )
        self.assertEqual(payload["families_meeting_minimum_independent_positives"], [])
        self.assertEqual(payload["minimum_independent_positive_count"], 3)
        summary = self.summary()
        self.assertEqual(
            summary["family_functional_positive_counts"], {"DED_hfam_4": 0, "DED_hfam_53": 1}
        )
        self.assertIn("counts_by_independence_group", summary["independence_summary"])
        self.assertEqual(summary["positive_count_basis"], "unique_independence_group")

    def test_input_contract_records_paths_sizes_and_hashes(self):
        fixture = duplicate_gene_fixture(self.root)
        contract = self.root / "input_contract.json"
        self.module.run(**fixture.kwargs(self.output_dir, input_contract=contract))
        payload = json.loads(contract.read_text(encoding="utf-8"))
        expected_paths = {
            "frozen_reference_ledger": fixture.ledger,
            "positive_primary_source_amendment": fixture.provenance,
            "family_definitions": fixture.family_definitions,
            "reference_fasta": fixture.fasta,
            "evidence_columns_config": EVIDENCE_COLUMNS_CONFIG,
        }
        for key, expected in expected_paths.items():
            entry = payload["inputs"][key]
            self.assertEqual(Path(entry["path"]).resolve(), expected.resolve(), key)
            self.assertEqual(entry["size"], expected.stat().st_size, key)
            self.assertEqual(
                entry["sha256"], hashlib.sha256(expected.read_bytes()).hexdigest(), key
            )
        self.assertEqual(
            payload["outputs"]["curated_ledger"]["row_count"], len(fixture.ledger_rows)
        )
        self.assertEqual(
            payload["frozen_ledger_sha256_before"], payload["frozen_ledger_sha256_after"]
        )
        self.assertTrue(payload["frozen_ledger_unchanged"])
        self.assertEqual(
            payload["outputs"]["curated_ledger"]["sha256"],
            hashlib.sha256(
                (self.output_dir / "phaded_reference_ledger_curated_v2.tsv").read_bytes()
            ).hexdigest(),
        )

    def test_report_is_written_with_the_curation_sections(self):
        fixture = duplicate_gene_fixture(self.root)
        report = self.root / "reference_evidence_curation_report.md"
        self.module.run(**fixture.kwargs(self.output_dir, report=report, write_report=True))
        text = report.read_text(encoding="utf-8")
        for heading in (
            "# PhaDED reference-ledger evidence curation (F13)",
            "## 3. Grade distribution",
            "## 4. Per-positive curation",
            "## 5. Independence-group derivation",
            "## 6. Per-family independent-positive counts",
            "## 7. Values left `pending`",
            "## 8. Acceptance gate (strict evidence-column mode)",
        ):
            self.assertIn(heading, text)
        self.assertIn("DED_hfam_53_0008", text)


if __name__ == "__main__":
    unittest.main()
