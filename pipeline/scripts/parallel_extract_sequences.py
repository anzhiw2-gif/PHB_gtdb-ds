#!/usr/bin/env python3
"""Extract requested protein records from one FASTA shard."""
from __future__ import annotations

import argparse
from pathlib import Path


def extract_shard(shard: Path, wanted: set[str], out: Path) -> int:
    count = 0
    current = None
    seq: list[str] = []
    with shard.open() as handle, out.open("w") as output:
        def flush() -> None:
            nonlocal count
            if current in wanted:
                output.write(f">{current}\n{''.join(seq)}\n")
                count += 1
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                current = line[1:].split()[0]
                seq = []
            elif current is not None:
                seq.append(line.strip())
        flush()
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=Path, required=True)
    parser.add_argument("--ids", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    wanted = {line.strip() for line in args.ids.open() if line.strip()}
    print(extract_shard(args.shard, wanted, args.out), flush=True)


if __name__ == "__main__":
    main()
