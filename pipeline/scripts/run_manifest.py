#!/usr/bin/env python3
"""run_manifest.py — 主管线运行清单固化（轻量）

将 run_pipeline.sh 追加的 JSONL 步骤记录（results/run_manifest.jsonl）固化为
results/run_manifest.json，并补充：
  - 元数据（生成时间、主机、git commit）
  - 输入文件 SHA-256（--inputs 列表）
  - 输出文件 SHA-256（--outputs 列表，缺失文件记为 null 并计为 missing）

用法:
  python run_manifest.py finalize \
      --jsonl results/run_manifest.jsonl \
      --out   results/run_manifest.json \
      --inputs  data/hmms/ePhaZ.hmm ... \
      --outputs results/tables/tier1_genome_family.tsv ...

schema_version: 1.0
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys


class ManifestError(RuntimeError):
    """Raised when a finalized manifest cannot prove required files exist."""


def validate_steps(steps):
    """Reject incomplete, failed, or ambiguous step records before finalization."""
    if not steps:
        raise ManifestError("no step records supplied; refusing to finalize an empty run")
    names = set()
    for record in steps:
        name = record.get("step") if isinstance(record, dict) else None
        if not name:
            raise ManifestError(f"invalid step record: {record!r}")
        if name in names:
            raise ManifestError(f"duplicate step record: {name}")
        names.add(name)
        if record.get("exit_code") != 0:
            raise ManifestError(f"step {name} did not succeed: exit_code={record.get('exit_code')!r}")
    return steps


def sha256(path: str):
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_paths(paths, label, *, allow_missing=False):
    """Return SHA-256 values for required paths, failing closed on missing/empty files."""
    hashes = {}
    missing = []
    for path in paths:
        path = os.fspath(path)
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            missing.append(path)
            continue
        hashes[path] = sha256(path)
    if missing and not allow_missing:
        raise ManifestError(f"required {label} missing or empty: {', '.join(missing)}")
    if allow_missing:
        for path in missing:
            hashes[path] = None
    return hashes


def source_bundle_sha256(source_files):
    """Return a deterministic digest over source paths and their file hashes."""
    digest = hashlib.sha256()
    for path, file_hash in sorted(source_files.items()):
        digest.update(os.fspath(path).encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_hash.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def collect_environment():
    """Capture execution identity without relying on optional external tools."""
    return {
        "entrypoint": os.path.abspath(sys.argv[0]),
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "cwd": os.getcwd(),
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
    }


def _step_commands(steps):
    return [record.get("command") for record in steps if record.get("command")]


def _validate_provenance(
    steps, source_files, gtdb_inputs, hmm_inputs, environment, commands, input_contract
):
    if not source_files:
        raise ManifestError("strict provenance requires non-empty source_files")
    if not gtdb_inputs and not input_contract:
        raise ManifestError("strict provenance requires GTDB inputs or an input contract")
    if not hmm_inputs:
        raise ManifestError("strict provenance requires non-empty hmm_inputs")
    if not environment:
        raise ManifestError("strict provenance requires a non-empty environment")
    if not commands or len(commands) != len(steps):
        raise ManifestError("strict provenance requires one command record per step")
    for record in steps:
        command = record.get("command")
        if not command:
            raise ManifestError(f"strict provenance missing command for step {record.get('step')!r}")


def build_manifest(
    steps,
    inputs,
    outputs,
    final_step=None,
    *,
    source_files=None,
    gtdb_inputs=None,
    hmm_inputs=None,
    environment=None,
    commands=None,
    strict=False,
    run_id=None,
    run_root=None,
    input_contract=None,
    allow_pending_gtdb=False,
):
    """Build a manifest, optionally including the successful finalizer step.

    The finalizer itself cannot be appended to JSONL after writing the JSON
    manifest without making the two records disagree.  Callers can therefore
    provide the already-started final step here; it is validated like every
    other step and becomes part of the immutable output.
    """
    all_steps = list(steps)
    if final_step is not None:
        all_steps.append(final_step)
    validate_steps(all_steps)
    source_hashes = validate_paths(source_files or [], "source") if source_files else {}
    gtdb_hashes = (
        validate_paths(gtdb_inputs or [], "GTDB input", allow_missing=allow_pending_gtdb)
        if gtdb_inputs else {}
    )
    hmm_hashes = validate_paths(hmm_inputs or [], "HMM input") if hmm_inputs else {}
    environment = collect_environment() if environment is None else environment
    commands = _step_commands(all_steps) if commands is None else commands
    contract_payload = None
    if input_contract:
        contract_path = os.fspath(input_contract)
        if not os.path.isfile(contract_path) or os.path.getsize(contract_path) == 0:
            raise ManifestError(f"input contract missing or empty: {contract_path}")
        with open(contract_path, encoding="utf-8") as handle:
            contract_payload = json.load(handle)
        if contract_payload.get("schema_version") != "1.0":
            raise ManifestError("unsupported input contract schema")
        if run_id is not None and contract_payload.get("run_id") != run_id:
            raise ManifestError("input contract run_id does not match manifest run_id")
        if run_root is not None:
            contract_root = contract_payload.get("run_dir")
            if not contract_root or os.path.realpath(contract_root) != os.path.realpath(run_root):
                raise ManifestError("input contract run_dir does not match manifest run_root")
    if strict:
        _validate_provenance(
            all_steps, source_hashes, gtdb_hashes, hmm_hashes, environment, commands,
            contract_payload,
        )
    bundle_hash = source_bundle_sha256(source_hashes) if source_hashes else None
    provenance = {
        "source_bundle_sha256": bundle_hash,
        "source_files": source_hashes,
        "environment": environment,
        "gtdb_inputs": gtdb_hashes,
        "hmm_inputs": hmm_hashes,
        "input_contract": contract_payload,
        "run_id": run_id,
        "run_root": os.path.realpath(run_root) if run_root else None,
        "commands": commands,
        "strict": bool(strict),
    }
    return {
        "schema_version": "1.1",
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "git_commit": git_commit(),
        "steps": all_steps,
        "inputs": validate_paths(inputs, "input") if inputs else {},
        "outputs": validate_paths(outputs, "output") if outputs else {},
        # Keep these fields at the top level for simple audit tooling.  The
        # nested copy is convenient for consumers that treat provenance as a
        # single object, and both representations are generated together.
        **provenance,
        "provenance": provenance,
    }


def git_commit():
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip() or None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("finalize")
    ap.add_argument("--jsonl", default="results/run_manifest.jsonl")
    ap.add_argument("--out", default="results/run_manifest.json")
    ap.add_argument("--inputs", nargs="*", default=[])
    ap.add_argument("--outputs", nargs="*", default=[])
    ap.add_argument("--final-step-name")
    ap.add_argument("--final-step-note", default="")
    ap.add_argument("--final-step-started")
    ap.add_argument("--source-files", nargs="*", default=[])
    ap.add_argument("--gtdb-inputs", nargs="*", default=[])
    ap.add_argument("--hmm-inputs", nargs="*", default=[])
    ap.add_argument("--input-contract")
    ap.add_argument("--run-id")
    ap.add_argument("--run-root")
    ap.add_argument("--allow-pending-gtdb", action="store_true")
    ap.add_argument("--strict-provenance", "--strict", dest="strict_provenance", action="store_true")
    ap.add_argument("--final-step-command", nargs="*", default=[])
    args = ap.parse_args()

    steps = []
    if os.path.exists(args.jsonl):
        with open(args.jsonl, encoding="utf-8-sig") as f:
            for line in f:
                line = line.strip()
                if line:
                    steps.append(json.loads(line))

    final_step = None
    if args.final_step_name:
        final_step = {
            "step": args.final_step_name,
            "exit_code": 0,
            "started": args.final_step_started or dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "ended": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "note": args.final_step_note,
        }
        if args.final_step_command:
            final_step["command"] = args.final_step_command
    manifest = build_manifest(
        steps,
        args.inputs,
        args.outputs,
        final_step=final_step,
        source_files=args.source_files,
        gtdb_inputs=args.gtdb_inputs,
        hmm_inputs=args.hmm_inputs,
        strict=args.strict_provenance,
        run_id=args.run_id,
        run_root=args.run_root,
        input_contract=args.input_contract,
        allow_pending_gtdb=args.allow_pending_gtdb,
    )

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(f"run_manifest.json -> {args.out}（{len(steps)} 步）")


if __name__ == "__main__":
    try:
        main()
    except ManifestError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
