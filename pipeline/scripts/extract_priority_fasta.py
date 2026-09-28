#!/usr/bin/env python3
"""Extract an ordered, accession-bound FASTA for a priority review table."""

from __future__ import annotations

import argparse
from pathlib import Path

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                accession = header.split(None, 1)[0]
                if accession in records or not chunks:
                    raise ValueError(f"duplicate or empty FASTA record: {accession}")
                records[accession] = (header, "".join(chunks))
            header, chunks = line[1:], []
        else:
            if header is None:
                raise ValueError("sequence appears before FASTA header")
            chunks.append(line)
    if header is not None:
        accession = header.split(None, 1)[0]
        if accession in records or not chunks:
            raise ValueError(f"duplicate or empty FASTA record: {accession}")
        records[accession] = (header, "".join(chunks))
    return records


def sanitize_sequence(sequence: str, accession: str, sanitize_terminal_stops: bool = False) -> tuple[str, bool]:
    sanitized = False
    if sanitize_terminal_stops and sequence.endswith("*"):
        sequence = sequence.rstrip("*")
        sanitized = True
    if any(residue not in AMINO_ACIDS for residue in sequence):
        raise ValueError(f"invalid amino acid in accession: {accession}")
    if not sequence:
        raise ValueError(f"empty sequence after normalization: {accession}")
    return sequence, sanitized


def extract(
    fasta: Path,
    accessions: list[str],
    output: Path,
    sanitize_terminal_stops: bool = False,
    exclude_invalid: Path | None = None,
) -> int:
    records = read_fasta(fasta)
    if len(accessions) != len(set(accessions)):
        raise ValueError("priority accession list contains duplicates")
    missing = [accession for accession in accessions if accession not in records]
    if missing:
        raise ValueError(f"accessions missing from FASTA: {','.join(missing[:5])}")
    output.parent.mkdir(parents=True, exist_ok=True)
    excluded: list[tuple[str, str]] = []
    with output.open("w", encoding="ascii", newline="\n") as handle:
        for accession in accessions:
            header, sequence = records[accession]
            try:
                sequence, _ = sanitize_sequence(sequence, accession, sanitize_terminal_stops)
            except ValueError:
                if exclude_invalid is None:
                    raise
                excluded.append((accession, "invalid_amino_acid"))
                continue
            handle.write(f">{header}\n{sequence}\n")
    if exclude_invalid is not None:
        exclude_invalid.parent.mkdir(parents=True, exist_ok=True)
        with exclude_invalid.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("accession\treason\n")
            for accession, reason in excluded:
                handle.write(f"{accession}\t{reason}\n")
    return len(accessions) - len(excluded)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--accessions", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sanitize-terminal-stops", action="store_true")
    parser.add_argument("--exclude-invalid", type=Path)
    args = parser.parse_args(argv)
    accessions = [line.strip() for line in args.accessions.read_text(encoding="utf-8").splitlines() if line.strip()]
    print(extract(args.fasta, accessions, args.output, args.sanitize_terminal_stops, args.exclude_invalid))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
