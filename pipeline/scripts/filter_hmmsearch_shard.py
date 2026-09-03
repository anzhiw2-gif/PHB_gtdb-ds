#!/usr/bin/env python3
"""Remove HMMER-incompatible protein records and record them as tool exclusions."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


MAX_TARGET_AA = 100000
EXCLUSION_LEDGER_NAME = "overlength_exclusions.tsv"


def filter_shard(source: Path, destination: Path, exclusions: Path) -> int:
    """Copy records at or below the HMMER target limit and ledger longer records."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    exclusions.parent.mkdir(parents=True, exist_ok=True)
    temp = destination.with_name(destination.name + ".tmp")
    count = 0
    excluded = []
    header = None
    sequence: list[str] = []

    def flush(handle) -> None:
        nonlocal count, header, sequence
        if header is None:
            return
        seq = "".join(sequence).replace(" ", "").replace("\t", "")
        if len(seq) > MAX_TARGET_AA:
            excluded.append((header[1:].split()[0], len(seq)))
        else:
            handle.write(header)
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")
            count += 1
        header = None
        sequence = []

    with source.open("r", encoding="utf-8") as src, temp.open("w", encoding="utf-8", newline="\n") as dst:
        for line in src:
            if line.startswith(">"):
                flush(dst)
                header = line
            elif header is not None:
                sequence.append(line.strip())
            elif line.strip():
                raise ValueError(f"sequence data before FASTA header in {source}")
        flush(dst)
    os.replace(temp, destination)
    with exclusions.open("a", encoding="utf-8", newline="") as ledger:
        for accession, length in excluded:
            ledger.write(f"{source.name}\t{accession}\t{length}\t>{MAX_TARGET_AA} aa HMMER tool limit\n")
    return len(excluded)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("exclusions", type=Path)
    args = parser.parse_args()
    filter_shard(args.source, args.destination, args.exclusions)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
