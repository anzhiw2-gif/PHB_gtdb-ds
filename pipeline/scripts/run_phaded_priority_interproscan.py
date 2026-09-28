#!/usr/bin/env python3
"""Run a bounded InterProScan confirmation for priority PhaDED candidates."""

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


def verify_input(path: Path, expected_sha256: str) -> int:
    actual = sha256(path)
    if actual.lower() != expected_sha256.lower():
        raise ValueError(f"SHA-256 mismatch for {path}: expected {expected_sha256}, observed {actual}")
    return path.stat().st_size


def interpro_command(interproscan: str, fasta: Path, output: Path, temp_dir: Path, cpu: int) -> list[str]:
    if cpu < 1 or cpu > 40:
        raise ValueError("cpu must be between 1 and 40")
    return [interproscan, "-i", str(fasta.resolve()), "-o", str(output.resolve()), "-f", "TSV", "--iprlookup", "-cpu", str(cpu), "-T", str(temp_dir.resolve())]


def run(interproscan: str, fasta: Path, output: Path, temp_dir: Path, cpu: int, expected_sha256: str) -> dict[str, object]:
    verify_input(fasta, expected_sha256)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    command = interpro_command(interproscan, fasta, output, temp_dir, cpu)
    completed = subprocess.run(command, capture_output=True, text=True, check=True)
    if not output.is_file():
        raise RuntimeError("InterProScan completed without creating output")
    return {
        "status": "completed_candidate_only",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "cpu": cpu,
        "input": {"path": str(fasta.resolve()), "size": fasta.stat().st_size, "sha256": sha256(fasta)},
        "output": {"path": str(output.resolve()), "size": output.stat().st_size, "sha256": sha256(output)},
        "stdout_bytes": len((completed.stdout or "").encode("utf-8")),
        "stderr_bytes": len((completed.stderr or "").encode("utf-8")),
        "phenotype_boundary": "InterProScan domain annotations are secondary candidate architecture evidence and do not prove PHB/PHA degradation.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interproscan", default="${PHB_REMOTE_ROOT}/software/bin/interproscan.sh")
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--temp-dir", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--cpu", type=int, default=40)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args(argv)
    result = run(args.interproscan, args.fasta, args.output, args.temp_dir, args.cpu, args.expected_sha256)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "output": result["output"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
