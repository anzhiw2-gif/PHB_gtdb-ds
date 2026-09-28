#!/usr/bin/env python3
"""Run bounded Pfam hmmscan stages for the PhaDED architecture evidence run."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scan_command(hmmscan: str, database: Path, fasta: Path, domtblout: Path, cpu: int) -> list[str]:
    return [hmmscan, "--noali", "--cpu", str(cpu), "--cut_ga", "--domtblout", str(domtblout), str(database), str(fasta)]


def run_scan(hmmscan: str, database: Path, fasta: Path, domtblout: Path, cpu: int) -> dict[str, object]:
    if cpu < 1 or cpu > 40:
        raise ValueError("cpu must be between 1 and 40")
    domtblout.parent.mkdir(parents=True, exist_ok=True)
    command = scan_command(hmmscan, database, fasta, domtblout, cpu)
    completed = subprocess.run(command, capture_output=True, text=True, check=True)
    return {
        "command": command,
        "stdout_bytes": len((completed.stdout or "").encode("utf-8")),
        "stderr_bytes": len((completed.stderr or "").encode("utf-8")),
        "database": {"path": str(database), "size": database.stat().st_size, "sha256": sha256(database)},
        "fasta": {"path": str(fasta), "size": fasta.stat().st_size, "sha256": sha256(fasta)},
        "domtblout": {"path": str(domtblout), "size": domtblout.stat().st_size, "sha256": sha256(domtblout)},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("reference", "candidate"), required=True)
    parser.add_argument("--hmmscan", default="${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb/bin/hmmscan")
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--domtblout", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=40)
    args = parser.parse_args()
    result = run_scan(args.hmmscan, args.database, args.fasta, args.domtblout, args.cpu)
    payload = {
        "schema_version": 1,
        "stage": args.stage,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "hmmscan": {"path": args.hmmscan, "version": subprocess.run([args.hmmscan, "-h"], capture_output=True, text=True, check=True).stdout.splitlines()[1]},
        "cpu": args.cpu,
        "scan": result,
        "phenotype_boundary": "Pfam domain evidence is candidate architecture support, not phenotype validation.",
    }
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    args.provenance.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage": args.stage, "domtblout": str(args.domtblout), "size": args.domtblout.stat().st_size}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
