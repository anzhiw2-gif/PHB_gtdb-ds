#!/usr/bin/env python3
"""Audit PhaDED profile manifest status and trained artifact hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit(manifest_path: Path, profiles_dir: Path, output: Path) -> dict[str, object]:
    with manifest_path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len({row["profile_id"] for row in rows}) != len(rows):
        raise ValueError("profile manifest contains duplicate profile_id")
    trained = [row for row in rows if row["model_status"] == "trained"]
    reference_only = [row for row in rows if row["model_status"] == "reference_only"]
    artifacts = []
    for row in trained:
        checks = []
        for suffix, field in (("alignment.faa", "alignment_sha256"), ("hmm", "hmm_sha256")):
            path = profiles_dir / f"{row['profile_id']}.{suffix}"
            if not path.is_file():
                raise ValueError(f"trained profile artifact missing: {path}")
            actual = sha256(path)
            if actual != row[field]:
                raise ValueError(f"trained profile artifact hash mismatch: {path}")
            checks.append({"path": str(path.resolve()), "size": path.stat().st_size, "sha256": actual})
        artifacts.append({"profile_id": row["profile_id"], "checks": checks})
    report = {
        "schema_version": 1,
        "manifest": {"path": str(manifest_path.resolve()), "size": manifest_path.stat().st_size, "sha256": sha256(manifest_path)},
        "profile_count": len(rows),
        "trained_count": len(trained),
        "reference_only_count": len(reference_only),
        "trained_profiles": artifacts,
        "reference_only_profiles": [row["profile_id"] for row in reference_only],
        "status": "verified",
        "interpretation": "reference_only families are retained as reference vocabulary and were not used as calibrated HMMs.",
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--profiles-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.manifest, args.profiles_dir, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
