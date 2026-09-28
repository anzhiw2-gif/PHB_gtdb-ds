#!/usr/bin/env python3
"""Remove redundant hmmscan stdout from a retained provenance JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def compact(source: Path, output: Path) -> None:
    payload = json.loads(source.read_text(encoding="utf-8"))
    scan = payload.get("scan", {})
    if isinstance(scan, dict):
        scan.pop("stdout", None)
        scan.pop("stderr", None)
        payload["scan"] = scan
    payload["stdout_removed"] = True
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    compact(args.source, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
