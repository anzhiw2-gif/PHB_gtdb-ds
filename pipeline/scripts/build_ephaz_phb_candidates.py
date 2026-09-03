#!/usr/bin/env python3
"""Build small, run-scoped PHB-focused ePhaZ seed candidates."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


DEFAULT_INDEPENDENT = ["AAB40611.1", "O24719", "A0A8W8", "Q9LBN6", "Q5YEW3"]


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                _add_record(records, header, chunks)
            header, chunks = line[1:], []
        elif header is None:
            raise ValueError(f"sequence precedes FASTA header in {path}")
        else:
            chunks.append(line.replace("-", "").replace(".", "").upper())
    if header is not None:
        _add_record(records, header, chunks)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _add_record(records: dict[str, tuple[str, str]], header: str, chunks: list[str]) -> None:
    accession = header.split("|", 1)[0].split()[0]
    if accession in records or not chunks:
        raise ValueError(f"invalid or duplicate FASTA accession: {accession}")
    records[accession] = (header, "".join(chunks))


def _write(path: Path, records: dict[str, tuple[str, str]], accessions: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(accessions):
            header, sequence = records[accession]
            handle.write(f">{header}\n{sequence}\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_candidates(core_fasta: Path | str, independent_fasta: Path | str, outdir: Path | str, independent_accessions: list[str] | None = None) -> dict[str, list[str]]:
    core = read_fasta(Path(core_fasta))
    independent = read_fasta(Path(independent_fasta))
    selected = independent_accessions or DEFAULT_INDEPENDENT
    missing = [acc for acc in selected if acc not in independent]
    if missing:
        raise ValueError(f"missing independent accession(s): {', '.join(missing)}")
    overlap = sorted(set(core) & set(independent))
    if overlap:
        raise ValueError(f"accession appears in both FASTA files: {', '.join(overlap)}")
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=False)
    no_q = sorted(acc for acc in core if acc != "Q51718")
    plus = sorted(no_q + selected)
    _write(out / "PHB-focused_no_Q51718.faa", core, no_q)
    merged = {**core, **independent}
    _write(out / "PHB-focused_plus_independent.faa", merged, plus)
    return {"PHB-focused_no_Q51718": no_q, "PHB-focused_plus_independent": plus}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-fasta", required=True, type=Path)
    parser.add_argument("--independent-fasta", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--independent-accession", action="append", dest="accessions")
    args = parser.parse_args(argv)
    result = build_candidates(args.core_fasta, args.independent_fasta, args.outdir, args.accessions)
    for name, accessions in result.items():
        path = args.outdir / f"{name}.faa"
        print(f"{name}\t{len(accessions)}\t{sha256_file(path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
