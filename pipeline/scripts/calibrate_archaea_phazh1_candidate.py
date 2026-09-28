#!/usr/bin/env python3
"""Run a bounded candidate-only hmmsearch calibration on fixed panels."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def search_command(hmm: Path, fasta: Path, tblout: Path, cpu: int = 1) -> list[str]:
    if cpu < 1:
        raise ValueError("cpu must be positive")
    return ["hmmsearch", "--noali", "--cpu", str(cpu), "--tblout", str(tblout), str(hmm), str(fasta)]


def parse_tblout(path: Path) -> list[dict[str, object]]:
    hits = []
    if not path.is_file():
        raise FileNotFoundError(path)
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 6:
            continue
        try:
            hits.append({"target": fields[0], "evalue": float(fields[4]), "bitscore": float(fields[5])})
        except ValueError:
            continue
    return hits


def calibrate(hmm: Path, positive: Path, negative: Path, outdir: Path, *, cpu: int = 1, hmmsearch: str = "hmmsearch") -> list[dict[str, object]]:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, fasta in (("positive", positive), ("negative", negative)):
        tblout = outdir / f"{label}.tblout"
        command = search_command(hmm, fasta, tblout, cpu)
        subprocess.run([hmmsearch, "--noali", "--cpu", str(cpu), "--tblout", str(tblout), str(hmm), str(fasta)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for hit in parse_tblout(tblout):
            rows.append({"panel": label, **hit})
    output = outdir / "calibration.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["panel", "target", "evalue", "bitscore"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmm", type=Path, required=True)
    parser.add_argument("--positive", type=Path, required=True)
    parser.add_argument("--negative", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--hmmsearch", default="hmmsearch")
    args = parser.parse_args()
    calibrate(args.hmm, args.positive, args.negative, args.outdir, cpu=args.cpu, hmmsearch=args.hmmsearch)


if __name__ == "__main__":
    main()
