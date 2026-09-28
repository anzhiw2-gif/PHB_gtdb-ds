#!/usr/bin/env python3
"""Materialize a source-bound PhaDED/DED 38-family reference library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable


KNOLL_DOI = "10.1186/1471-2105-10-89"
KNOLL_PMCID = "PMC2666664"
DED_BASE = "http://www.ded.uni-stuttgart.de"
ALLOWED_SUPERFAMILIES = {
    "intracellular nPHASCL without lipase box",
    "intracellular nPHASCL with lipase box",
    "periplasmic PHA depolymerases",
    "intracellular nPHAMCL",
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
    "extracellular native-SCL/PhaZ7-like",
    "extracellular dPHAMCL",
}
LEDGER_FIELDS = [
    "reference_id", "accession", "sequence_sha256", "sequence_source", "source_database",
    "source_database_version", "retrieval_date", "organism", "taxonomy_id", "phaded_superfamily",
    "phaded_family_id", "reported_localization", "substrate_class", "substrate_detail",
    "experimental_assay", "experimental_result", "evidence_status", "positive_negative_control",
    "catalytic_residues", "lipase_box_or_ahsmg", "oxyanion_hole_evidence", "domain_architecture",
    "primary_doi", "pmid", "pmcid", "notes",
]


def _clean(value: str) -> str:
    value = html.unescape(re.sub(r"<[^>]+>", " ", value))
    return re.sub(r"\s+", " ", value).strip()


def normalize_superfamily(label: str) -> str:
    lower = label.lower()
    if "i-nphascl" in lower and "no lipase" in lower:
        return "intracellular nPHASCL without lipase box"
    if "i-nphascl" in lower:
        return "intracellular nPHASCL with lipase box"
    if "i-nphamcl" in lower:
        return "intracellular nPHAMCL"
    if "p-nphascl" in lower or "periplasm" in lower:
        return "periplasmic PHA depolymerases"
    if "e-dphascl" in lower and "type 1" in lower:
        return "extracellular dPHASCL type 1"
    if "e-dphascl" in lower and "type 2" in lower:
        return "extracellular dPHASCL type 2"
    if "e-nphascl" in lower:
        return "extracellular native-SCL/PhaZ7-like"
    if "e-dphamcl" in lower:
        return "extracellular dPHAMCL"
    raise ValueError(f"unknown DED superfamily label: {label}")


def parse_family_index(path: Path) -> list[dict[str, str]]:
    text = path.read_text(encoding="latin1")
    pattern = re.compile(r'<a href="index\.pl\?page=hfam&id=(\d+)">(.*?)</a>', re.I | re.S)
    definitions = []
    seen: set[str] = set()
    for family_numeric_id, raw_label in pattern.findall(text):
        if family_numeric_id in seen:
            continue
        seen.add(family_numeric_id)
        label = _clean(raw_label)
        definitions.append({
            "source_hfam_id": family_numeric_id,
            "source_family_label": label,
            "phaded_superfamily": normalize_superfamily(label),
            "phaded_family_id": f"DED_hfam_{family_numeric_id}",
            "source_family_url": f"{DED_BASE}/cgi-bin/DEDTAX/index.pl?page=hfam&id={family_numeric_id}",
            "source_sequence_url": f"{DED_BASE}/CONTENT/seq{family_numeric_id}.txt",
            "source_alignment_url": f"{DED_BASE}/CONTENT/aln{family_numeric_id}.aln",
            "source_hmm_url": f"{DED_BASE}/CONTENT/aln{family_numeric_id}.hmm",
        })
    return definitions


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    accession = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if accession is not None:
                if not chunks:
                    raise ValueError(f"empty DED sequence: {accession}")
                records.append((accession, "".join(chunks)))
            accession = line[1:].split()[0]
            if not accession:
                raise ValueError(f"empty DED FASTA header: {path}")
            chunks = []
        else:
            if accession is None:
                raise ValueError(f"sequence before header in {path}")
            chunks.append(line)
    if accession is not None:
        if not chunks:
            raise ValueError(f"empty DED sequence: {accession}")
        records.append((accession, "".join(chunks)))
    return records


def seed_gis_from_knoll_xml(path: Path | None) -> set[str]:
    if path is None:
        return set()
    root = ET.fromstring(path.read_text(encoding="utf-8"))
    table = next((node for node in root.iter() if node.tag.endswith("table-wrap") and node.attrib.get("id") == "T1"), None)
    if table is None:
        raise ValueError("Knoll XML does not contain Table 1 (T1)")
    seed_gis: set[str] = set()
    for row in table.iter():
        if not row.tag.endswith("tr"):
            continue
        cells = ["".join(cell.itertext()).strip() for cell in row if cell.tag.endswith("td")]
        if cells and re.fullmatch(r"\d+(?:\s*\([^)]*\))?", cells[0]):
            seed_gis.add(re.match(r"\d+", cells[0]).group(0))
    return seed_gis


def _row_records(page_text: str) -> list[tuple[str, str, set[str]]]:
    rows = []
    for raw_row in re.findall(r"<TR>(.*?)</TR>", page_text, flags=re.I | re.S):
        cells = re.findall(r"<TD[^>]*>(.*?)</TD>", raw_row, flags=re.I | re.S)
        if len(cells) < 4:
            continue
        organism = _clean(cells[0])
        identifiers = set(re.findall(r"(?:[A-Z]{1,6}\d{3,9}(?:\.\d+)?|\d{5,12})", cells[3]))
        accessions = {item for item in identifiers if not item.isdigit()}
        rows.append((organism, ";".join(sorted(accessions)), identifiers))
    return rows


def _record_metadata(page_text: str, accession: str, seed_gis: set[str]) -> tuple[str, str]:
    for organism, accessions, identifiers in _row_records(page_text):
        if accession in accessions:
            matched = sorted(seed_gis.intersection(identifiers))
            return organism, ";".join(matched)
    return "pending", ""


def _substrate(superfamily: str) -> tuple[str, str, str]:
    localization = "intracellular" if superfamily.startswith("intracellular") else "periplasmic" if superfamily.startswith("periplasmic") else "extracellular"
    substrate = "MCL" if "MCL" in superfamily else "SCL"
    detail = "native PHA" if "native" in superfamily or superfamily.startswith("intracellular") or superfamily.startswith("periplasmic") else "denatured PHA"
    return localization, substrate, detail


def materialize(raw_dir: str | Path, output_dir: str | Path, *, expected_family_count: int = 38, seed_gis: set[str] | None = None, retrieval_date: str | None = None) -> dict[str, object]:
    raw_dir, output_dir = Path(raw_dir), Path(output_dir)
    retrieval_date = retrieval_date or date.today().isoformat()
    seed_gis = seed_gis or set()
    definitions = parse_family_index(raw_dir / "family_index.html")
    if len(definitions) != expected_family_count:
        raise ValueError(f"DED family index contains {len(definitions)} families; expected {expected_family_count}")
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    fasta_records: list[tuple[str, str]] = []
    seen_hashes: dict[str, tuple[str, str]] = {}
    duplicate_rows: list[dict[str, str]] = []
    for definition in definitions:
        family_numeric_id = definition["source_hfam_id"]
        seq_path = raw_dir / f"seq{family_numeric_id}.txt"
        page_path = raw_dir / f"family_{family_numeric_id}.html"
        if not seq_path.is_file() or not seq_path.stat().st_size:
            raise ValueError(f"missing DED family FASTA: {seq_path}")
        page_text = page_path.read_text(encoding="latin1") if page_path.is_file() else ""
        for ordinal, (accession, sequence) in enumerate(parse_fasta(seq_path), start=1):
            digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            if digest in seen_hashes:
                duplicate_rows.append({
                    "family_id": definition["phaded_family_id"], "accession": accession,
                    "sequence_sha256": digest, "kept_reference_id": seen_hashes[digest][0],
                    "kept_family_id": seen_hashes[digest][1],
                })
                continue
            reference_id = f"{definition['phaded_family_id']}_{ordinal:04d}"
            seen_hashes[digest] = (reference_id, definition["phaded_family_id"])
            organism, matched_gi = _record_metadata(page_text, accession, seed_gis)
            status = "experimental_positive" if matched_gi else "annotation_only"
            localization, substrate_class, substrate_detail = _substrate(definition["phaded_superfamily"])
            note = f"DED family page {definition['source_family_url']}; DED sequence file {definition['source_sequence_url']}"
            if matched_gi:
                note += f"; Knoll 2009 Table 1 seed GI {matched_gi}"
            rows.append({
                "reference_id": reference_id, "accession": accession, "sequence_sha256": digest,
                "sequence_source": definition["source_sequence_url"], "source_database": "PhaDED/DED",
                "source_database_version": "DED live snapshot; Knoll et al. 2009",
                "retrieval_date": retrieval_date, "organism": organism, "taxonomy_id": "pending",
                "phaded_superfamily": definition["phaded_superfamily"], "phaded_family_id": definition["phaded_family_id"],
                "reported_localization": localization, "substrate_class": substrate_class,
                "substrate_detail": substrate_detail, "experimental_assay": "Knoll 2009 Table 1 seed designation" if matched_gi else "pending",
                "experimental_result": "validated PHA depolymerase activity reported by Knoll et al. 2009" if matched_gi else "pending",
                "evidence_status": status, "positive_negative_control": "positive_seed" if matched_gi else "not_assessed",
                "catalytic_residues": "pending", "lipase_box_or_ahsmg": "pending",
                "oxyanion_hole_evidence": "pending", "domain_architecture": "pending",
                "primary_doi": KNOLL_DOI, "pmid": "19296857", "pmcid": KNOLL_PMCID, "notes": note,
            })
            fasta_records.append((reference_id, accession, sequence))
    definitions_fields = [
        "phaded_superfamily", "phaded_family_id", "source_hfam_id", "source_family_label",
        "source_family_url", "source_sequence_url", "source_alignment_url", "source_hmm_url",
        "source_family_page_sha256", "source_sequence_file_sha256", "source_sequence_count",
    ]
    with (output_dir / "phaded_family_definitions.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=definitions_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for definition in definitions:
            family_numeric_id = definition["source_hfam_id"]
            page_path = raw_dir / f"family_{family_numeric_id}.html"
            seq_path = raw_dir / f"seq{family_numeric_id}.txt"
            record = dict(definition)
            record.update({
                "source_family_page_sha256": hashlib.sha256(page_path.read_bytes()).hexdigest(),
                "source_sequence_file_sha256": hashlib.sha256(seq_path.read_bytes()).hexdigest(),
                "source_sequence_count": str(sum(row["phaded_family_id"] == definition["phaded_family_id"] for row in rows)),
            })
            writer.writerow(record)
    with (output_dir / "phaded_reference_ledger.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with (output_dir / "phaded_reference.faa").open("w", encoding="ascii", newline="\n") as handle:
        for reference_id, accession, sequence in fasta_records:
            handle.write(f">{reference_id}|{accession}\n{sequence}\n")
    with (output_dir / "ded_duplicate_sequence_records.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["family_id", "accession", "sequence_sha256", "kept_reference_id", "kept_family_id"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(duplicate_rows)
    result = {
        "retrieval_date": retrieval_date, "family_count": len(definitions),
        "raw_sequence_entry_count": sum(len(parse_fasta(raw_dir / f"seq{d['source_hfam_id']}.txt")) for d in definitions),
        "reference_count": len(rows), "duplicate_sequence_count": len(duplicate_rows),
        "experimental_positive_count": sum(row["evidence_status"] == "experimental_positive" for row in rows),
        "annotation_only_count": sum(row["evidence_status"] == "annotation_only" for row in rows),
        "source_database": "PhaDED/DED", "knoll_doi": KNOLL_DOI, "knoll_pmcid": KNOLL_PMCID,
    }
    (output_dir / "materialization.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--knoll-xml", type=Path)
    parser.add_argument("--retrieval-date", default=None)
    parser.add_argument("--expected-family-count", type=int, default=38)
    args = parser.parse_args()
    result = materialize(
        args.raw_dir, args.output_dir,
        expected_family_count=args.expected_family_count,
        seed_gis=seed_gis_from_knoll_xml(args.knoll_xml),
        retrieval_date=args.retrieval_date,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
