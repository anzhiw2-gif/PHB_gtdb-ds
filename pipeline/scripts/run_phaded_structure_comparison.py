#!/usr/bin/env python3
"""Run bounded Foldseek comparisons for an accession-bound PhaDED pair table."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path


OUTPUT_FIELDS = [
    "candidate_id", "reference_model_id", "query", "target", "fident", "alnlen",
    "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits",
    "prob", "alntmscore", "qtmscore", "ttmscore", "rmsd", "status", "raw_output",
]
ALIGN_FIELDS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend",
    "tstart", "tend", "evalue", "bits", "prob", "alntmscore", "qtmscore", "ttmscore", "rmsd",
]


def normalize_fieldnames(fieldnames: list[str]) -> list[str]:
    if not fieldnames:
        return []
    return [fieldnames[0].lstrip("\ufeff")] + fieldnames[1:]


def validate_threads(threads: int, *, server: bool) -> None:
    if threads < 1:
        raise ValueError("threads must be >= 1")
    if server and threads > 40:
        raise ValueError("server threads must be <= 40")


def validate_pair(row: dict[str, str]) -> None:
    for key in ("candidate_pdb", "reference_pdb"):
        path = Path(row.get(key, ""))
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"{key}: {path}")


def tool_version(executable: str) -> str:
    result = subprocess.run([executable, "version"], check=True, capture_output=True, text=True)
    return result.stdout.strip() or result.stderr.strip()


def run_pairs(pairs: Path, executable: str, output: Path, tmp: Path, threads: int, *, server: bool) -> dict[str, object]:
    validate_threads(threads, server=server)
    output.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)
    version = tool_version(executable)
    rows: list[dict[str, str]] = []
    with pairs.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = normalize_fieldnames(reader.fieldnames or [])
        required = {"candidate_id", "reference_model_id", "candidate_pdb", "reference_pdb"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing pair columns: {sorted(missing)}")
        for index, pair in enumerate(reader, start=1):
            validate_pair(pair)
            raw = output / f"pair_{index:03d}.m8"
            command = [
                executable, "easy-search", pair["candidate_pdb"], pair["reference_pdb"],
                str(raw), str(tmp / f"tmp_{index:03d}"), "--threads", str(threads),
                "--format-output",
                "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore,qtmscore,ttmscore,rmsd",
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"Foldseek failed for {pair['candidate_id']}: {result.stderr.strip()}")
            lines = raw.read_text(encoding="utf-8").splitlines() if raw.exists() else []
            values = lines[0].split("\t") if lines else [""] * len(ALIGN_FIELDS)
            values += [""] * (len(ALIGN_FIELDS) - len(values))
            row = dict(zip(ALIGN_FIELDS, values[: len(ALIGN_FIELDS)]))
            row.update({
                "candidate_id": pair["candidate_id"],
                "reference_model_id": pair["reference_model_id"],
                "status": "completed",
                "raw_output": str(raw),
            })
            rows.append(row)
    result_path = output / "comparison_results.tsv"
    with result_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "pair_count": len(rows), "completed_count": len(rows), "tool": "Foldseek",
        "tool_version": version, "threads": threads, "server_thread_limit": 40 if server else None,
        "results": str(result_path), "candidate_only": True,
    }
    (output / "comparison_execution_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--foldseek", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--tmp", required=True, type=Path)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--server", action="store_true")
    args = parser.parse_args()
    run_pairs(args.pairs, args.foldseek, args.output, args.tmp, args.threads, server=args.server)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
