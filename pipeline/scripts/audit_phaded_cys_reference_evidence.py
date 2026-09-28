#!/usr/bin/env python3
"""Audit the three externally-declared accessions O87189 / Q7WT48 / Q7WT49
against the frozen 723-row PhaDED reference ledger (2026-09-17, Step 2).

This script is *read-only* with respect to every frozen artefact.  It joins
externally retrieved database evidence (UniProt records + PubMed abstracts,
stored verbatim in the run's ``inputs/external_evidence/``) to the frozen
ledger and emits:

* ``results/reference_evidence_audit.tsv`` - one row per audited accession with
  its evidence status, source PMID, retrieval provenance and its relation to
  the existing ledger records;
* ``results/cys_evidence_statistics_per_run.json`` - the Cys superfamily
  evidence statistics *as measured in this run*, which explicitly does not
  modify the frozen ledger.

Boundary: every field recorded here is label/sequence evidence only.  Nothing in
this file is a validated PHB/PHA degradation phenotype, a family call, or a
registry entry.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
from pathlib import Path

CYS_SUPERFAMILY = "intracellular nPHASCL without lipase box"
AUDIT_ACCESSIONS = ("O87189", "Q7WT48", "Q7WT49")
GATE_POSITIVE_REQUIREMENT = 3
BOUNDARY = "candidate_only"
BOUNDARY_NOTE = (
    "candidate_only evidence audit recorded in this run only; "
    "does_not_modify_frozen_ledger; label/sequence evidence only, "
    "not a validated PHB/PHA degradation phenotype"
)

AUDIT_FIELDS = [
    "accession",
    "uniprot_entry_name",
    "protein_name",
    "gene_name",
    "organism",
    "sequence_length",
    "source_pmid",
    "source_doi",
    "evidence_class",
    "experimental_positive",
    "evidence_statement",
    "primary_literature_finding",
    "embl_protein_accession",
    "sequence_identical_to_ledger_accession",
    "same_protein_ledger_records",
    "ledger_reference_id",
    "ledger_superfamily",
    "ledger_family_id",
    "ledger_evidence_status",
    "in_frozen_ledger",
    "is_alias_of_existing_ledger_record",
    "is_duplicate_of_existing_frozen_positive",
    "task5_panel_label",
    "task5_panel_label_assessment",
    "retrieval_date",
    "retrieval_url",
    "retrieved_file",
    "retrieved_file_sha256",
    "retrieval_http_status",
    "relation_to_existing_ledger",
    "boundary_note",
]


# --------------------------------------------------------------------------- io
def read_tsv(path: os.PathLike[str] | str) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: os.PathLike[str] | str, fieldnames: list[str], rows: list[dict]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: os.PathLike[str] | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_of_protein_sequence(sequence: str) -> str:
    return sha256_text(sequence)


# ------------------------------------------------------------------- reference
def load_reference_index(path: os.PathLike[str] | str) -> tuple[dict[str, str], dict[str, str]]:
    """Return (by_reference_id, by_accession) for ``>reference_id|accession`` FASTAs."""
    by_id: dict[str, str] = {}
    by_accession: dict[str, str] = {}
    header: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8", newline="") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    _store(header, parts, by_id, by_accession)
                header, parts = line[1:], []
            elif line.strip():
                parts.append(line.strip())
        if header is not None:
            _store(header, parts, by_id, by_accession)
    return by_id, by_accession


def _store(header: str, parts: list[str], by_id: dict[str, str], by_accession: dict[str, str]) -> None:
    sequence = "".join(parts)
    reference_id, accession = (header.split("|", 1) + [""])[:2]
    by_id[reference_id] = sequence
    by_accession[accession] = sequence


def exact_match_targets(sequence: str, reference_by_accession: dict[str, str]) -> list[str]:
    if not sequence:
        return []
    return sorted(acc for acc, target in reference_by_accession.items() if target and target == sequence)


def sequence_relations(sequence: str, reference_by_accession: dict[str, str]) -> dict[str, str]:
    """Classify the relation of ``sequence`` to every reference sequence."""
    relations: dict[str, str] = {}
    if not sequence:
        return relations
    for accession, target in reference_by_accession.items():
        if not target:
            continue
        if target == sequence:
            relations[accession] = "identical"
        elif target in sequence:
            relations[accession] = "target_is_c_terminal_subsequence_of_query"
        elif sequence in target:
            relations[accession] = "query_is_subsequence_of_target"
    return relations


# -------------------------------------------------------------------- uni prot
_ECO_RE = re.compile(r"\{ECO:[^}]*\}")


def _clean(value: str) -> str:
    return _ECO_RE.sub("", value).strip().strip("{").strip("}").strip()


def parse_uniprot_record(text: str) -> dict:
    """Parse the subset of a UniProtKB flat-file record used by this audit."""
    record: dict = {
        "accession": "",
        "entry_name": "",
        "protein_name": "",
        "gene_name": "",
        "organism": "",
        "sequence_length": None,
        "embl_protein_accessions": [],
        "pubmed_ids": [],
        "doi": [],
        "sequence": "",
    }
    sequence_parts: list[str] = []
    in_sequence = False
    for line in text.splitlines():
        code, _, rest = line[:2], None, line[5:].strip() if len(line) > 5 else ""
        if code == "ID":
            fields = rest.split()
            if fields:
                record["entry_name"] = fields[0]
            match = re.search(r"(\d+)\s+AA", rest)
            if match:
                record["sequence_length"] = int(match.group(1))
            in_sequence = False
        elif code == "AC":
            if not record["accession"]:
                record["accession"] = rest.split(";")[0].strip()
        elif code == "DE":
            if "SubName: Full=" in rest and not record["protein_name"]:
                record["protein_name"] = _clean(rest.split("SubName: Full=", 1)[1])
            elif "Full=" in rest and not record["protein_name"]:
                record["protein_name"] = _clean(rest.split("Full=", 1)[1])
        elif code == "GN":
            match = re.search(r"Name=([^\s;{]+)", rest)
            if match and not record["gene_name"]:
                record["gene_name"] = match.group(1)
        elif code == "OS":
            record["organism"] = (record["organism"] + " " + rest).strip()
        elif code == "RX":
            match = re.search(r"PubMed=(\d+)", rest)
            if match:
                record["pubmed_ids"].append(match.group(1))
            match = re.search(r"DOI=([^;]+)", rest)
            if match:
                record["doi"].append(match.group(1).strip())
        elif code == "DR":
            fields = [f.strip() for f in rest.split(";")]
            if fields and fields[0] == "EMBL" and len(fields) > 2:
                record["embl_protein_accessions"].append(fields[2])
        elif code == "SQ":
            in_sequence = True
            match = re.search(r"SEQUENCE\s+(\d+)\s+AA", rest)
            if match:
                record["sequence_length"] = int(match.group(1))
        elif in_sequence:
            if code == "//":
                in_sequence = False
            else:
                sequence_parts.append(re.sub(r"[^A-Za-z]", "", line))
    record["sequence"] = "".join(sequence_parts).upper()
    return record


def load_uniprot_records(evidence_dir: os.PathLike[str] | str) -> dict[str, dict]:
    records: dict[str, dict] = {}
    base = Path(evidence_dir)
    for accession in AUDIT_ACCESSIONS:
        path = base / f"uniprot_{accession}.txt"
        if not path.is_file():
            continue
        record = parse_uniprot_record(path.read_text(encoding="utf-8"))
        record["retrieved_file"] = path.name
        record["retrieved_file_sha256"] = sha256_file(path)
        records[accession] = record
    return records


# ------------------------------------------------------------------ audit rows
def _ledger_by_accession(ledger_rows: list[dict]) -> dict[str, dict]:
    return {row.get("accession", ""): row for row in ledger_rows if row.get("accession")}


def build_audit_rows(
    *,
    retrieval_rows: list[dict],
    uniprot_records: dict[str, dict],
    reference_by_accession: dict[str, str],
    ledger_rows: list[dict],
) -> list[dict]:
    ledger_by_accession = _ledger_by_accession(ledger_rows)
    reference_id_by_accession = {
        row["accession"]: row.get("reference_id", "") for row in ledger_rows if row.get("accession")
    }
    rows: list[dict] = []
    for claim in retrieval_rows:
        accession = claim.get("accession", "")
        record = uniprot_records.get(accession, {})
        sequence = record.get("sequence") or ""
        relations = sequence_relations(sequence, reference_by_accession)
        identical = [acc for acc, rel in relations.items() if rel == "identical"]
        embl = claim.get("embl_protein_accession") or (
            record.get("embl_protein_accessions") or [""])[0]

        # the ledger record this UniProt entry corresponds to: prefer a
        # byte-identical reference, else the EMBL cross-reference
        match_accession = ""
        if identical:
            match_accession = sorted(identical)[0]
        elif embl in ledger_by_accession:
            match_accession = embl
        ledger_row = ledger_by_accession.get(match_accession, {})

        same_protein: list[str] = []
        for acc, rel in sorted(relations.items()):
            if rel in {"identical", "target_is_c_terminal_subsequence_of_query"}:
                same_protein.append(f"{reference_id_by_accession.get(acc, '')}|{acc}")

        in_ledger = "true" if ledger_row else "false"
        is_alias = "true" if (identical and ledger_row) or embl in ledger_by_accession else "false"
        ledger_status = ledger_row.get("evidence_status", "")
        is_duplicate_of_positive = (
            "true" if ledger_status == "experimental_positive" else "false"
        )
        if identical:
            relation = (
                f"sequence_identical_to_frozen_ledger_record_{match_accession}"
                if ledger_row else "sequence_identical_to_reference_fasta_but_no_ledger_row"
            )
        elif match_accession:
            relation = f"embl_cross_reference_matches_frozen_ledger_record_{match_accession}"
        else:
            relation = "no_matching_frozen_ledger_record"

        rows.append({
            "accession": accession,
            "uniprot_entry_name": claim.get("uniprot_entry_name") or record.get("entry_name", ""),
            "protein_name": claim.get("protein_name") or record.get("protein_name", ""),
            "gene_name": claim.get("gene_name") or record.get("gene_name", ""),
            "organism": record.get("organism", ""),
            "sequence_length": claim.get("sequence_length") or (
                record.get("sequence_length") or len(sequence) or ""),
            "source_pmid": claim.get("source_pmid", ""),
            "source_doi": claim.get("source_doi", ""),
            "evidence_class": claim.get("evidence_class", ""),
            "experimental_positive": claim.get("experimental_positive", ""),
            "evidence_statement": claim.get("evidence_statement", ""),
            "primary_literature_finding": claim.get("primary_literature_finding", ""),
            "embl_protein_accession": embl,
            "sequence_identical_to_ledger_accession": match_accession if identical else "",
            "same_protein_ledger_records": ";".join(same_protein),
            "ledger_reference_id": ledger_row.get("reference_id", ""),
            "ledger_superfamily": ledger_row.get("phaded_superfamily", ""),
            "ledger_family_id": ledger_row.get("phaded_family_id", ""),
            "ledger_evidence_status": ledger_status,
            "in_frozen_ledger": in_ledger,
            "is_alias_of_existing_ledger_record": is_alias,
            "is_duplicate_of_existing_frozen_positive": is_duplicate_of_positive,
            "task5_panel_label": claim.get("task5_panel_label", ""),
            "task5_panel_label_assessment": claim.get("task5_panel_label_assessment", ""),
            "retrieval_date": claim.get("retrieval_date", ""),
            "retrieval_url": claim.get("retrieval_url", ""),
            "retrieved_file": record.get("retrieved_file", ""),
            "retrieved_file_sha256": record.get("retrieved_file_sha256", ""),
            "retrieval_http_status": (
                claim.get("retrieval_http_status") or claim.get("retrieved_http_status", "")),
            "relation_to_existing_ledger": relation,
            "boundary_note": BOUNDARY_NOTE,
        })
    return rows


def summarize_evidence(rows: list[dict], ledger_rows: list[dict]) -> dict:
    cys = [r for r in ledger_rows if r.get("phaded_superfamily") == CYS_SUPERFAMILY]
    frozen_positive_ids = {
        r["reference_id"] for r in cys if r.get("evidence_status") == "experimental_positive"
    }
    audit_positive_ids = {
        r["ledger_reference_id"] for r in rows
        if r.get("experimental_positive") == "true" and r.get("ledger_reference_id")
    }
    added = sorted(audit_positive_ids - frozen_positive_ids)
    distinct = sorted(frozen_positive_ids | audit_positive_ids)
    aliases = sorted(
        r["accession"] for r in rows if r.get("is_alias_of_existing_ledger_record") == "true"
    )
    return {
        "schema": "phaded-cys-reference-evidence-audit-v1",
        "frozen_ledger_cys_reference_count": len(cys),
        "frozen_ledger_cys_experimental_positive_count": len(frozen_positive_ids),
        "frozen_ledger_cys_experimental_positive_reference_ids": sorted(frozen_positive_ids),
        "audited_accessions": sorted(r["accession"] for r in rows),
        "audited_accessions_that_are_frozen_ledger_aliases": aliases,
        "distinct_experimental_positive_cys_records": len(distinct),
        "distinct_experimental_positive_cys_reference_ids": distinct,
        "experimental_positive_records_added_by_this_audit": added,
        "records_with_unresolved_function": sorted(
            r["ledger_reference_id"] for r in rows
            if r.get("experimental_positive") == "false" and r.get("ledger_reference_id")
        ),
        "gate_positive_requirement": GATE_POSITIVE_REQUIREMENT,
        "gate_positive_requirement_met": len(distinct) >= GATE_POSITIVE_REQUIREMENT,
        "frozen_ledger_modified": False,
        "frozen_ledger_edited_in_place": "no",
        "evidence_layer": "reference_layer_external_database_audit",
        "not_a_phenotype_claim": True,
    }


def write_audit_outputs(out_tsv, out_json, rows: list[dict], summary: dict) -> dict:
    write_tsv(out_tsv, AUDIT_FIELDS, rows)
    payload = dict(summary)
    payload["boundary"] = BOUNDARY
    payload["audit_tsv"] = str(out_tsv)
    payload["audit_row_count"] = len(rows)
    payload["generated_at_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    target = Path(out_json)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--reference-faa", required=True)
    parser.add_argument("--retrieval-table", required=True)
    parser.add_argument("--evidence-dir", default=None)
    args = parser.parse_args(argv)

    run_dir = Path(args.run_dir)
    evidence_dir = args.evidence_dir or (run_dir / "inputs" / "external_evidence")
    _by_id, by_accession = load_reference_index(args.reference_faa)
    rows = build_audit_rows(
        retrieval_rows=read_tsv(args.retrieval_table),
        uniprot_records=load_uniprot_records(evidence_dir),
        reference_by_accession=by_accession,
        ledger_rows=read_tsv(args.ledger),
    )
    summary = summarize_evidence(rows, read_tsv(args.ledger))
    payload = write_audit_outputs(
        run_dir / "results" / "reference_evidence_audit.tsv",
        run_dir / "results" / "cys_evidence_statistics_per_run.json",
        rows,
        summary,
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
