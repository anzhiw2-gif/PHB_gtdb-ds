#!/usr/bin/env python3
"""Append independently curated boundary proteins to a run-scoped panel."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def _read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                _add(records, header, chunks)
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError(f"sequence precedes FASTA header: {path}")
        else:
            chunks.append(line.replace("-", "").replace(".", "").upper())
    if header is not None:
        _add(records, header, chunks)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _add(records: dict[str, tuple[str, str]], header: str, chunks: list[str]) -> None:
    fields = header.split("|")
    accession = fields[1].split()[0] if fields and fields[0].lower() in {"sp", "tr", "ref", "gi"} and len(fields) > 1 else fields[0].split()[0]
    if accession in records or not chunks:
        raise ValueError(f"invalid or duplicate accession: {accession}")
    records[accession] = (header, "".join(chunks))


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def extend_panel(panel: Path | str, additions: Path | str, ledger: Path | str, panel_name: str) -> list[str]:
    panel_path, additions_path, ledger_path = Path(panel), Path(additions), Path(ledger)
    existing = _read_fasta(panel_path)
    new = _read_fasta(additions_path)
    overlap = sorted(set(existing) & set(new))
    if overlap:
        raise ValueError(f"accession already present in panel: {', '.join(overlap)}")
    with panel_path.open("a", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(new):
            header, sequence = new[accession]
            handle.write(f">{header}\n{sequence}\n")
    with ledger_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("accession\tpanel\tsequence_length\tsequence_sha256\tsource_fasta\tsource_sha256\n")
        source_sha = hashlib.sha256(additions_path.read_bytes()).hexdigest()
        for accession in sorted(new):
            _, sequence = new[accession]
            handle.write(f"{accession}\t{panel_name}\t{len(sequence)}\t{_sha(sequence)}\t{additions_path.resolve()}\t{source_sha}\n")
    return sorted(new)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--additions", required=True, type=Path)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--panel-name", required=True)
    args = parser.parse_args(argv)
    for accession in extend_panel(args.panel, args.additions, args.ledger, args.panel_name):
        print(accession)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
