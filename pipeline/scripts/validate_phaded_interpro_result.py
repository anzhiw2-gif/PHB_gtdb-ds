#!/usr/bin/env python3
"""Validate an InterProScan TSV against its accession-bound FASTA input."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Iterable


EXPECTED_COLUMNS = 15
HEX_MD5 = re.compile(r"^[0-9a-f]{32}$", re.IGNORECASE)
IPR_ACCESSION = re.compile(r"^IPR\d+$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                if current in records or not chunks:
                    raise ValueError(f"duplicate or empty FASTA record: {current}")
                records[current] = "".join(chunks)
            current = line[1:].split(None, 1)[0]
            if not current:
                raise ValueError("empty FASTA accession")
            chunks = []
        else:
            if current is None:
                raise ValueError("sequence appears before FASTA header")
            chunks.append(line)
    if current is not None:
        if current in records or not chunks:
            raise ValueError(f"duplicate or empty FASTA record: {current}")
        records[current] = "".join(chunks)
    if not records:
        raise ValueError("FASTA is empty")
    return records


def read_accessions(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="ascii").splitlines() if line.strip()]


def read_excluded(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row.get("accession", "").strip(): row.get("reason", "").strip() for row in rows if row.get("accession", "").strip()}


def _top(counter: Counter[str], limit: int = 20) -> dict[str, int]:
    return {key: value for key, value in counter.most_common(limit)}


def validate_paths(result_path: Path, fasta_path: Path, accession_path: Path, excluded_path: Path) -> dict[str, object]:
    sequences = read_fasta(fasta_path)
    candidates = read_accessions(accession_path)
    excluded = read_excluded(excluded_path)
    hard_errors: list[str] = []
    if len(candidates) != len(set(candidates)):
        hard_errors.append("candidate accession list contains duplicates")
    candidate_set = set(candidates)
    fasta_set = set(sequences)
    if not fasta_set <= candidate_set:
        hard_errors.append("FASTA contains accession outside candidate list")
    if not set(excluded) <= candidate_set:
        hard_errors.append("excluded table contains accession outside candidate list")
    if set(excluded) & fasta_set:
        hard_errors.append("excluded accession is present in normalized FASTA")

    field_errors = 0
    md5_mismatches = 0
    length_mismatches = 0
    coordinate_errors = 0
    numeric_errors = 0
    status_errors = 0
    date_errors = 0
    extra_accessions: set[str] = set()
    result_accessions: set[str] = set()
    exact_rows: Counter[tuple[str, ...]] = Counter()
    analyses: Counter[str] = Counter()
    signatures: Counter[str] = Counter()
    interpros: Counter[str] = Counter()
    rows = 0
    with result_path.open(encoding="utf-8-sig", newline="") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            rows += 1
            exact_rows[tuple(fields)] += 1
            if len(fields) != EXPECTED_COLUMNS:
                field_errors += 1
                continue
            accession, protein_md5, length_text, analysis, signature = fields[:5]
            result_accessions.add(accession)
            if accession not in sequences:
                extra_accessions.add(accession)
            else:
                sequence = sequences[accession]
                if not HEX_MD5.fullmatch(protein_md5) or protein_md5.lower() != hashlib.md5(sequence.encode("ascii")).hexdigest():
                    md5_mismatches += 1
                try:
                    if int(length_text) != len(sequence):
                        length_mismatches += 1
                except ValueError:
                    length_mismatches += 1
            try:
                start, end = int(fields[6]), int(fields[7])
                seq_len = int(length_text)
                if start < 1 or end < start or end > seq_len:
                    coordinate_errors += 1
            except ValueError:
                coordinate_errors += 1
            if fields[8] != "-":
                try:
                    float(fields[8])
                except ValueError:
                    numeric_errors += 1
            if fields[9] not in {"T", "F"}:
                status_errors += 1
            try:
                datetime.strptime(fields[10], "%d-%m-%Y")
            except ValueError:
                date_errors += 1
            if fields[11] != "-" and not IPR_ACCESSION.fullmatch(fields[11]):
                hard_errors.append(f"invalid InterPro accession at line {line_number}")
            analyses[analysis] += 1
            signatures[signature] += 1
            if fields[11] != "-":
                interpros[fields[11]] += 1

    duplicate_rows = sum(count - 1 for count in exact_rows.values() if count > 1)
    missing = candidate_set - result_accessions
    missing_excluded = missing & set(excluded)
    missing_unreported = missing - set(excluded)
    hard_count = len(hard_errors) + field_errors + md5_mismatches + length_mismatches + coordinate_errors + numeric_errors + status_errors + date_errors + len(extra_accessions)
    if hard_count:
        status = "validation_failed"
    elif missing_unreported:
        status = "validated_with_unreported_accessions"
    else:
        status = "validated_candidate_only"
    return {
        "status": status,
        "result_path": str(result_path.resolve()),
        "result_size": result_path.stat().st_size,
        "result_sha256": sha256(result_path),
        "fasta_path": str(fasta_path.resolve()),
        "fasta_size": fasta_path.stat().st_size,
        "fasta_sha256": sha256(fasta_path),
        "candidate_count": len(candidate_set),
        "fasta_accession_count": len(fasta_set),
        "result_row_count": rows,
        "result_accession_count": len(result_accessions),
        "excluded_count": len(excluded),
        "missing_count": len(missing),
        "missing_excluded_count": len(missing_excluded),
        "missing_unreported_count": len(missing_unreported),
        "missing_unreported_examples": sorted(missing_unreported)[:20],
        "extra_accession_count": len(extra_accessions),
        "duplicate_exact_row_count": duplicate_rows,
        "field_error_count": field_errors,
        "md5_mismatch_count": md5_mismatches,
        "length_mismatch_count": length_mismatches,
        "coordinate_error_count": coordinate_errors,
        "numeric_error_count": numeric_errors,
        "status_error_count": status_errors,
        "date_error_count": date_errors,
        "analysis_counts": _top(analyses),
        "signature_counts": _top(signatures),
        "interpro_counts": _top(interpros),
        "errors": hard_errors[:20],
        "phenotype_boundary": "InterProScan annotations validate sequence/domain records only; they do not prove PHB/PHA degradation phenotype.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--accessions", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = validate_paths(args.result, args.fasta, args.accessions, args.excluded)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, sort_keys=True))
    return 0 if report["status"] != "validation_failed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
