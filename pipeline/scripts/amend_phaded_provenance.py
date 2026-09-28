#!/usr/bin/env python3
"""Create a dated provenance correction for completed PhaDED evidence runs."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from datetime import date
from pathlib import Path
from typing import Mapping


PHENOTYPE_BOUNDARY = (
    "All profile, domain, motif, localization, structure, phylogeny, and "
    "Foldseek fields remain candidate-only and do not prove PHB/PHA degradation phenotype."
)

RUN_CONTEXT_SCRIPT = Path(__file__).resolve().parent / "run_context.py"
#: the word this script used to put in the contract's top-level ``status``; it is a
#: *completion* state, not a contract status, so it now travels in
#: ``run_completion_status`` and the top-level ``status`` stays in the shared
#: ``run_context`` domain (``verified`` / ``pending``).
COMPLETION_STATUS_KEY = "run_completion_status"


def load_run_context() -> object:
    """Load the shared contract implementation by path (never a private copy).

    The Phase 1 P1 defect was exactly that this script wrote its own bespoke
    contract dict instead of delegating to ``run_context``; loading the sibling
    module by path keeps the run layout verifiable without any import-path
    assumption about the caller's working directory.
    """
    spec = importlib.util.spec_from_file_location("run_context", RUN_CONTEXT_SCRIPT)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise RuntimeError(f"cannot load {RUN_CONTEXT_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_run_input_contract(
    *,
    run_dir: Path,
    run_id: str,
    manifest: Mapping[str, object],
    execution: Mapping[str, object],
) -> dict[str, object]:
    """Write ``input_contract.json`` through ``run_context`` and add provenance.

    The contract skeleton comes from ``run_context.write_input_contract`` (so the
    core keys, the ``gtdb`` block, the per-input ``status`` and the top-level
    status domain cannot drift from the shared implementation), and the
    provenance-only payload is added as extra top-level keys afterwards.  The
    amendment's completion word never occupies the top-level ``status``.
    """
    run_context = load_run_context()
    contract = run_context.write_input_contract(
        run_dir,
        run_id=run_id,
        gtdb_inputs={name: None for name in run_context.GTDB_INPUT_NAMES},
        inputs={
            name: Path(str(entry["path"]))
            for name, entry in dict(manifest["inputs"]).items()
        },
    )
    inputs_status = sorted({str(item["status"]) for item in contract["inputs"].values()})
    if inputs_status != ["verified"]:
        raise ValueError(
            "every declared input must re-verify against disk before the contract is "
            f"written; observed statuses: {inputs_status}"
        )
    completion_status = str(manifest["status"])
    contract.update(
        {
            "execution": dict(execution),
            "phenotype_boundary": PHENOTYPE_BOUNDARY,
            COMPLETION_STATUS_KEY: completion_status,
            "run_inputs_status": "verified",
            "status_note": (
                "the top-level status follows pipeline/scripts/run_context.py and is 'pending' "
                "only because the three GTDB input slots (taxonomy/metadata/tree) are declared "
                "pending: this amendment consumes no GTDB file. Every declared input is "
                f"verified. The amendment completion state '{completion_status}' is carried in "
                f"{COMPLETION_STATUS_KEY}, not in status."
            ),
        }
    )
    output = Path(run_dir) / "input_contract.json"
    output.write_text(
        json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return contract


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bound_file(path: Path) -> dict[str, object]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"required regular file missing: {path}")
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}


def build_provenance_amendment(
    *,
    run_id: str,
    output_dir: Path,
    signalp_summary: Path,
    foldseek_manifest: Path,
    final_matrix: Path,
    server_observation: Mapping[str, object],
) -> dict[str, object]:
    signalp = json.loads(signalp_summary.read_text(encoding="utf-8"))
    foldseek = json.loads(foldseek_manifest.read_text(encoding="utf-8"))
    if signalp.get("status") != "completed" or signalp.get("returncode") != 0:
        raise ValueError("SignalP summary is not a successful completed run")
    if foldseek.get("status") not in {"completed_candidate_only", "completed"}:
        raise ValueError("Foldseek manifest is not a completed candidate-only run")
    server = str(server_observation.get("server", "")).strip()
    account = str(server_observation.get("account", "")).strip()
    if not server or not account:
        raise ValueError("server and account are required")
    signalp_threads = int(server_observation.get("signalp_threads", 0))
    if signalp_threads < 1 or signalp_threads > 40:
        raise ValueError("SignalP thread count must be between 1 and 40")

    output_dir.mkdir(parents=True, exist_ok=True)
    inputs = {
        "signalp_run_summary": bound_file(signalp_summary),
        "foldseek_merge_manifest": bound_file(foldseek_manifest),
        "final_foldseek_matrix": bound_file(final_matrix),
    }
    execution = {
        "mode": "server_dated_deploy",
        "server": server,
        "account": account,
        "signalp_returncode": int(signalp["returncode"]),
        "signalp_threads": signalp_threads,
        "signalp_command": signalp.get("command", []),
        "gpu_execution": str(server_observation.get("gpu_execution", "not_recorded")),
        "gpu_observation": str(server_observation.get("gpu_observation", "not_recorded")),
        "server_observed_at": str(server_observation.get("observed_at", "not_recorded")),
        "logical_cpu": server_observation.get("logical_cpu", "pending"),
        "available_memory": server_observation.get("available_memory", "pending"),
        "deploy_path": str(server_observation.get("deploy_path", "pending")),
    }
    outputs = {
        "provenance_amendment": {
            "path": str((output_dir / "provenance_amendment.json").resolve()),
        },
        "provenance_readme": {
            "path": str((output_dir / "provenance_readme.md").resolve()),
        },
    }
    manifest = {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at": date.today().isoformat(),
        "status": "completed_candidate_only",
        "inputs": inputs,
        "execution": execution,
        "outputs": outputs,
        "correction": {
            "previous_execution_label": "local_read_only/server_not_used",
            "corrected_execution_label": "SignalP completed on T141 dated deploy; Foldseek merged locally from accession-bound output",
            "no_historical_run_overwritten": True,
        },
        "profile_boundary": "No new family call or profile promotion is authorized by this amendment.",
        "phenotype_boundary": PHENOTYPE_BOUNDARY,
    }
    amendment_path = output_dir / "provenance_amendment.json"
    amendment_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    outputs["provenance_amendment"].update(
        {"size": amendment_path.stat().st_size, "sha256": sha256(amendment_path)}
    )
    readme = output_dir / "provenance_readme.md"
    readme.write_text(
        "\n".join(
            [
                f"# {run_id}",
                "",
                "本 amendment 只修正执行 provenance，不重算、不覆盖历史结果。",
                "",
                f"- SignalP：在 `{execution['server']}` 的 dated deploy 上完成，返回码 `{execution['signalp_returncode']}`，线程 `{execution['signalp_threads']}`。",
                f"- Foldseek：使用已完成的 accession-bound 结果，仅作为新增结构证据列合并；未测试的全库记录保持 `foldseek_not_tested_full_library`。",
                "- 不新增 family 判定，不解除 architecture conflict，不把计算证据写成 PHB/PHA 表型证明。",
                "",
                PHENOTYPE_BOUNDARY,
                "",
            ]
        ),
        encoding="utf-8",
    )
    outputs["provenance_readme"].update({"size": readme.stat().st_size, "sha256": sha256(readme)})
    manifest["outputs"] = outputs
    amendment_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_run_input_contract(
        run_dir=output_dir.parent,
        run_id=run_id,
        manifest=manifest,
        execution=execution,
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--signalp-summary", type=Path, required=True)
    parser.add_argument("--foldseek-manifest", type=Path, required=True)
    parser.add_argument("--final-matrix", type=Path, required=True)
    parser.add_argument("--server", required=True)
    parser.add_argument("--account", required=True)
    parser.add_argument("--signalp-threads", type=int, required=True)
    parser.add_argument("--logical-cpu", type=int, default=0)
    parser.add_argument("--available-memory", default="pending")
    parser.add_argument("--deploy-path", default="pending")
    parser.add_argument("--observed-at", default="not_recorded")
    parser.add_argument("--gpu-execution", default="not_recorded")
    parser.add_argument("--gpu-observation", default="not_recorded")
    args = parser.parse_args(argv)
    result = build_provenance_amendment(
        run_id=args.run_id,
        output_dir=args.output_dir,
        signalp_summary=args.signalp_summary,
        foldseek_manifest=args.foldseek_manifest,
        final_matrix=args.final_matrix,
        server_observation={
            "server": args.server,
            "account": args.account,
            "signalp_threads": args.signalp_threads,
            "logical_cpu": args.logical_cpu,
            "available_memory": args.available_memory,
            "deploy_path": args.deploy_path,
            "observed_at": args.observed_at,
            "gpu_execution": args.gpu_execution,
            "gpu_observation": args.gpu_observation,
        },
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
