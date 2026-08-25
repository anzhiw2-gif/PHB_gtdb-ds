#!/usr/bin/env python3
"""Create isolated run directories and auditable input contracts."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Mapping


RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
GTDB_INPUT_NAMES = ("taxonomy", "metadata", "tree")


def validate_run_id(run_id: str) -> str:
    """Validate a run identifier before it is used as a directory name."""
    if not isinstance(run_id, str) or not RUN_ID_RE.fullmatch(run_id):
        raise ValueError(
            "run_id must be 1-64 ASCII characters matching [A-Za-z0-9][A-Za-z0-9._-]*"
        )
    if run_id in {".", ".."} or ".." in run_id:
        raise ValueError("run_id must not contain path traversal")
    return run_id


def create_run_layout(root: os.PathLike[str] | str, run_id: str) -> Path:
    """Create and return ``root/runs/<run_id>`` with stable output folders."""
    run_id = validate_run_id(run_id)
    run_dir = Path(root).expanduser() / "runs" / run_id
    if run_dir.exists() or run_dir.is_symlink():
        raise FileExistsError(f"run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=False)
    for name in ("logs", "inputs", "results"):
        (run_dir / name).mkdir(exist_ok=True)
    return run_dir


def sha256_file(path: os.PathLike[str] | str) -> str:
    """Return the SHA-256 digest of a non-directory file."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _path_string(path: os.PathLike[str] | str | None) -> str | None:
    if path is None:
        return None
    return os.path.abspath(os.path.expanduser(os.fspath(path)))


def _describe_input(
    path: os.PathLike[str] | str | None, *, missing_status: str
) -> dict[str, object | None]:
    path_string = _path_string(path)
    record: dict[str, object | None] = {
        "path": path_string,
        "status": missing_status,
        "sha256": None,
        "size": None,
    }
    if path_string is None:
        return record
    candidate = Path(path_string)
    if candidate.is_file() and candidate.stat().st_size > 0:
        record.update(
            status="verified",
            sha256=sha256_file(candidate),
            size=candidate.stat().st_size,
        )
    return record


def build_input_contract(
    run_dir: os.PathLike[str] | str,
    *,
    run_id: str | None = None,
    gtdb_inputs: Mapping[str, os.PathLike[str] | str | None] | None = None,
    inputs: Mapping[str, os.PathLike[str] | str | None] | None = None,
) -> dict[str, object]:
    """Build a JSON-serializable input contract without inventing missing hashes."""
    run_dir = Path(run_dir).expanduser()
    declared_gtdb = dict(gtdb_inputs or {})
    gtdb = {
        name: _describe_input(declared_gtdb.get(name), missing_status="pending")
        for name in GTDB_INPUT_NAMES
    }
    other = {
        name: _describe_input(path, missing_status="missing")
        for name, path in (inputs or {}).items()
    }
    records = [*gtdb.values(), *other.values()]
    status = "verified" if all(item["status"] == "verified" for item in records) else "pending"
    return {
        "schema_version": "1.0",
        "run_id": validate_run_id(run_id if run_id is not None else run_dir.name),
        "run_dir": str(run_dir.resolve()),
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "status": status,
        "gtdb": gtdb,
        "inputs": other,
    }


def write_input_contract(
    run_dir: os.PathLike[str] | str,
    *,
    run_id: str | None = None,
    gtdb_inputs: Mapping[str, os.PathLike[str] | str | None] | None = None,
    inputs: Mapping[str, os.PathLike[str] | str | None] | None = None,
) -> dict[str, object]:
    """Write ``input_contract.json`` in an existing run directory."""
    run_dir = Path(run_dir).expanduser()
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)
    contract = build_input_contract(
        run_dir, run_id=run_id, gtdb_inputs=gtdb_inputs, inputs=inputs
    )
    output = run_dir / "input_contract.json"
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(contract, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return contract


__all__ = [
    "GTDB_INPUT_NAMES",
    "RUN_ID_RE",
    "build_input_contract",
    "create_run_layout",
    "sha256_file",
    "validate_run_id",
    "write_input_contract",
]
