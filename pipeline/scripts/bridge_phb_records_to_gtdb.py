#!/usr/bin/env python3
"""Bridge accession-level PHB literature records to local GTDB records.

The bridge is deliberately conservative: only an exact accession match can
retain direct literature support. Organism/strain name matches are retained as
taxonomy-only context and never promote a record to Tier A.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Iterable, Sequence


OUTPUT_COLUMNS = [
    "record_id", "literature_genome_accession", "literature_protein_accession",
    "literature_gene_accession", "organism_name", "strain_name",
    "gtdb_genome_accession", "gtdb_taxonomy", "gtdb_protein_accession",
    "gtdb_gene_accession", "match_class", "mapping_method", "evidence_tier",
    "literature_evidence_status", "exact_identifier_scope", "source_database",
    "primary_doi", "pmid", "retrieval_date", "notes",
]

TIER_A = "direct_literature_supported"
# The literature ledger uses ``experimental_positive`` while older callers
# may already provide the normalized Tier A label.
TIER_A_EVIDENCE_STATUSES = {"experimental_positive", TIER_A}


def _read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _value(row: dict[str, str], key: str) -> str:
    value = (row.get(key) or "").strip()
    if value:
        return value
    # Existing Tier1 summary uses ``genome``; accept the alias while keeping
    # the normalized output column as gtdb_genome_accession.
    if key == "genome_accession":
        return (row.get("genome") or "").strip()
    return ""


def _accession_aliases(row: dict[str, str]) -> set[str]:
    """Return normalized assembly/accession aliases exposed by GTDB metadata."""
    keys = (
        "genome_accession", "genome", "gtdb_genome_accession",
        "accession", "ncbi_genbank_assembly_accession",
        "ncbi_refseq_assembly_accession", "ncbi_wgs_master",
    )
    aliases = {_value(row, key) for key in keys}
    return {value for value in aliases if value and value.casefold() != "none"}


def _norm(value: str) -> str:
    return " ".join(value.casefold().split())


def _taxonomy_matches(lit: dict[str, str], tax_rows: list[dict[str, str]]) -> list[tuple[str, dict[str, str], str, str]]:
    genome = _value(lit, "genome_accession")
    protein = _value(lit, "protein_accession")
    gene = _value(lit, "gene_accession")
    organism = _norm(_value(lit, "organism_name"))
    strain = _norm(_value(lit, "strain_name"))
    matches: list[tuple[str, dict[str, str], str, str]] = []
    if genome:
        for row in tax_rows:
            if genome in _accession_aliases(row):
                matches.append(("exact_genome_accession", row, "genome_accession", "exact genome accession"))
        if matches:
            return matches
    if protein:
        for row in tax_rows:
            if _value(row, "protein_accession") == protein and _value(row, "genome_accession"):
                matches.append(("exact_protein_accession", row, "protein_accession", "protein accession bound to GTDB genome"))
        if matches:
            return matches
    if gene:
        for row in tax_rows:
            if _value(row, "gene_accession") == gene and _value(row, "genome_accession"):
                matches.append(("exact_gene_accession", row, "gene_accession", "gene accession bound to GTDB genome"))
        if matches:
            return matches
    if organism:
        for row in tax_rows:
            row_org = _norm(_value(row, "organism_name"))
            row_strain = _norm(_value(row, "strain_name"))
            # Names are context only; even a species match cannot establish
            # strain identity, so retain all organism-level matches together.
            if row_org and row_org == organism:
                matches.append(("taxonomy_only_match", row, "organism_strain_name", "name-only taxonomy match"))
        if matches:
            return matches
    return [("unresolved_no_gtdb_record", {}, "none", "no exact accession or taxonomy record")]


def _tier(lit: dict[str, str], match_class: str) -> str:
    status = _value(lit, "evidence_status").casefold()
    if match_class.startswith("exact_") and status in TIER_A_EVIDENCE_STATUSES:
        return TIER_A
    if match_class.startswith("exact_"):
        return "genome_supported_candidate"
    if match_class == "taxonomy_only_match":
        return "sequence_only_candidate"
    return "hold_unresolved"


def bridge_records(literature_path: str | Path, taxonomy_path: str | Path, output_path: str | Path) -> list[dict[str, str]]:
    literature = _read_tsv(literature_path)
    taxonomy = _read_tsv(taxonomy_path)
    rows: list[dict[str, str]] = []
    for lit in literature:
        for match_class, tax, method, note in _taxonomy_matches(lit, taxonomy):
            row = {name: "" for name in OUTPUT_COLUMNS}
            row.update({
                "record_id": _value(lit, "record_id"),
                "literature_genome_accession": _value(lit, "genome_accession"),
                "literature_protein_accession": _value(lit, "protein_accession"),
                "literature_gene_accession": _value(lit, "gene_accession"),
                "organism_name": _value(lit, "organism_name"),
                "strain_name": _value(lit, "strain_name"),
                "gtdb_genome_accession": _value(tax, "gtdb_genome_accession") or _value(tax, "genome_accession") or _value(tax, "genome") or _value(tax, "accession"),
                "gtdb_taxonomy": _value(tax, "gtdb_taxonomy"),
                "gtdb_protein_accession": _value(tax, "protein_accession"),
                "gtdb_gene_accession": _value(tax, "gene_accession"),
                "match_class": match_class,
                "mapping_method": method,
                "evidence_tier": _tier(lit, match_class),
                "literature_evidence_status": _value(lit, "evidence_status"),
                "exact_identifier_scope": _value(lit, "exact_identifier_scope"),
                "source_database": _value(lit, "source_database"),
                "primary_doi": _value(lit, "primary_doi"),
                "pmid": _value(lit, "pmid"),
                "retrieval_date": _value(lit, "retrieval_date"),
                "notes": "; ".join(filter(None, [_value(lit, "notes"), note])),
            })
            rows.append(row)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _file_binding(path: str | Path) -> dict[str, object]:
    p = Path(path)
    digest = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
    return {"path": str(p.resolve()), "size": p.stat().st_size if p.is_file() else None, "sha256": digest, "status": "verified" if digest else "pending"}


def write_manifest(manifest_path: str | Path, literature_path: str | Path, taxonomy_path: str | Path, output_path: str | Path, gtdb_release: str, retrieval_date: str | None = None, metadata_path: str | Path | None = None) -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": "1.0",
        "generated_at": date.today().isoformat(),
        "gtdb_release": gtdb_release,
        "retrieval_date": retrieval_date or date.today().isoformat(),
        "mapping_policy": "exact accession only for direct_literature_supported; name-only matches never Tier A",
        "source_run_ids": ["20260909_phaded_candidate_mapping_01"],
        "files": {
            "literature": _file_binding(literature_path),
            "taxonomy": _file_binding(taxonomy_path),
            "metadata": _file_binding(metadata_path) if metadata_path else {"path": "pending", "size": None, "sha256": None, "status": "pending"},
            "output": _file_binding(output_path),
        },
    }
    path = Path(manifest_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--literature", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--gtdb-release", default="pending")
    parser.add_argument("--retrieval-date")
    parser.add_argument("--metadata")
    args = parser.parse_args(argv)
    rows = bridge_records(args.literature, args.taxonomy, args.output)
    write_manifest(args.manifest, args.literature, args.taxonomy, args.output, args.gtdb_release, args.retrieval_date, args.metadata)
    print(f"wrote {len(rows)} bridge rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
