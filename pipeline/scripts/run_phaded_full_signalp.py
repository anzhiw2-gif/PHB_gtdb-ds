#!/usr/bin/env python3
"""Run a dated, candidate-only SignalP 6 full-library scan."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


MAX_THREADS = 40


def build_signalp_command(
    fasta: Path,
    output_dir: Path,
    *,
    torch_threads: int,
    write_processes: int,
) -> list[str]:
    if torch_threads < 1 or write_processes < 1:
        raise ValueError("thread and writer counts must be positive")
    if torch_threads + write_processes > MAX_THREADS:
        raise ValueError("combined SignalP threads exceed project limit of 40")
    if not fasta.is_file() or fasta.is_symlink():
        raise ValueError(f"FASTA is not a regular file: {fasta}")
    return [
        "signalp6",
        "--fastafile",
        str(fasta),
        "--output_dir",
        str(output_dir),
        "--format",
        "txt",
        "--organism",
        "other",
        "--mode",
        "fast",
        "--torch_num_threads",
        str(torch_threads),
        "--write_procs",
        str(write_processes),
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastafile", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--torch-threads", type=int, default=32)
    parser.add_argument("--write-processes", type=int, default=8)
    args = parser.parse_args(argv)
    command = build_signalp_command(
        args.fastafile,
        args.output_dir,
        torch_threads=args.torch_threads,
        write_processes=args.write_processes,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, check=False)
    summary = {
        "status": "completed" if completed.returncode == 0 else "failed",
        "returncode": completed.returncode,
        "command": command,
        "candidate_only": True,
        "phenotype_boundary": "SignalP is localization evidence only and does not prove PHB/PHA degradation phenotype.",
    }
    (args.output_dir / "signalp_run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
