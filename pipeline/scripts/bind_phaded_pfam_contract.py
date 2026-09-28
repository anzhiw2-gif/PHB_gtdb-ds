#!/usr/bin/env python3
"""Bind Pfam database and scan policy to the architecture evidence contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--database-path", required=True)
    parser.add_argument("--database-size", type=int, required=True)
    parser.add_argument("--database-sha256", required=True)
    parser.add_argument("--hmmscan-path", required=True)
    parser.add_argument("--hmmscan-size", type=int, required=True)
    parser.add_argument("--hmmscan-sha256", required=True)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    contract["status"] = "prepared_candidate_only"
    contract["domain_evidence_policy"] = {
        "database": "Pfam-A",
        "database_path": args.database_path,
        "database_size": args.database_size,
        "database_sha256": args.database_sha256,
        "hmmscan_path": args.hmmscan_path,
        "hmmscan_size": args.hmmscan_size,
        "hmmscan_sha256": args.hmmscan_sha256,
        "hmmscan_version": "HMMER 3.4 (Aug 2023)",
        "cutoff": "--cut_ga",
        "server_cpu_max": 40,
        "reference_first": True,
        "interpretation": "Reference-derived Pfam fingerprints define targeted architecture support; absent or untested domains remain pending, never biological negatives.",
    }
    args.contract.write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": contract["status"], "database": args.database_path}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
