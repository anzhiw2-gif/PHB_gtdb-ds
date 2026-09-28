#!/usr/bin/env python3
"""Fail-closed validation for the accession-level PHB literature ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Sequence


REQUIRED_COLUMNS = {
    "record_id", "organism_name", "strain_name", "protein_accession", "gene_accession",
    "genome_accession", "ncbi_taxid", "gtdb_taxonomy", "source_database",
    "source_database_version", "retrieval_date", "primary_doi", "pmid", "pmcid",
    "experimental_substrate", "experimental_result", "evidence_type", "evidence_status",
    "exact_identifier_scope", "sequence_sha256", "notes",
}
ALLOWED_STATUSES = {
    "experimental_positive", "experimental_negative", "annotation_only",
    "challenge_control", "pending_review",
}
ALLOWED_TYPES = {
    "direct_experiment", "other_pha_experiment", "generic_polyester_experiment",
    "annotation_only", "challenge_control", "structure_only", "pending_review",
}
SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")


def _read(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"ledger is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return fields, list(reader)


def validate_ledger(ledger_path: str | Path, output_path: str | Path | None = None) -> dict:
    ledger = Path(ledger_path)
    fields, rows = _read(ledger)
    missing = sorted(REQUIRED_COLUMNS - set(fields))
    if missing:
        raise ValueError("ledger missing required columns: " + ",".join(missing))
    seen_records: set[str] = set()
    seen_proteins: set[str] = set()
    statuses: Counter[str] = Counter()
    types: Counter[str] = Counter()
    failures: list[str] = []
    for index, row in enumerate(rows, start=2):
        record_id = (row.get("record_id") or "").strip()
        protein = (row.get("protein_accession") or "").strip()
        if not record_id:
            raise ValueError(f"row {index}: missing record_id")
        if record_id in seen_records:
            raise ValueError(f"duplicate record_id: {record_id}")
        seen_records.add(record_id)
        if not protein:
            raise ValueError(f"row {index}: missing protein_accession")
        if protein in seen_proteins:
            raise ValueError(f"duplicate protein_accession: {protein}")
        seen_proteins.add(protein)
        status = (row.get("evidence_status") or "").strip()
        if status not in ALLOWED_STATUSES:
            raise ValueError(f"{protein}: invalid evidence_status: {status}")
        citation = any((row.get(key) or "").strip() for key in ("primary_doi", "pmid", "pmcid"))
        source_database = (row.get("source_database") or "").strip()
        if not citation and not (status in {"annotation_only", "challenge_control", "pending_review"} and source_database):
            raise ValueError(f"{protein}: missing citation")
        evidence_type = (row.get("evidence_type") or "").strip()
        if evidence_type not in ALLOWED_TYPES:
            raise ValueError(f"{protein}: invalid evidence_type: {evidence_type}")
        digest = (row.get("sequence_sha256") or "").strip()
        if digest != "pending" and not SHA256.fullmatch(digest):
            raise ValueError(f"{protein}: invalid sequence_sha256")
        scope = (row.get("exact_identifier_scope") or "").strip()
        if not scope or scope in {"name_only", "genus_only", "species_only"}:
            if status == "experimental_positive":
                raise ValueError(f"{protein}: experimental_positive lacks exact identifier scope")
        if status == "experimental_positive":
            substrate = (row.get("experimental_substrate") or "").lower()
            result = (row.get("experimental_result") or "").strip()
            if "phb" not in substrate and "poly(3-hydroxybutyrate)" not in substrate:
                raise ValueError(f"{protein}: experimental_positive must name PHB substrate")
            if not result:
                raise ValueError(f"{protein}: experimental_positive lacks experimental_result")
            if evidence_type != "direct_experiment":
                raise ValueError(f"{protein}: experimental_positive requires direct_experiment type")
        statuses[status] += 1
        types[evidence_type] += 1
    report = {
        "ledger": str(ledger),
        "row_count": len(rows),
        "evidence_status_counts": dict(sorted(statuses.items())),
        "evidence_type_counts": dict(sorted(types.items())),
        "duplicate_accessions": [],
        "validation_failures": failures,
        "source_file_sha256": hashlib.sha256(ledger.read_bytes()).hexdigest(),
        "status": "valid",
    }
    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger")
    parser.add_argument("--output")
    args = parser.parse_args(argv)
    report = validate_ledger(args.ledger, args.output)
    print(f"validated {report['row_count']} PHB literature ledger rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
