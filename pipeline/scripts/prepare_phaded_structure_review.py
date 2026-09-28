#!/usr/bin/env python3
"""Prepare candidate-only structure-review inputs from triaged decisions.

This script prepares a review packet.  It never runs a structure predictor or a
structure comparison, it never deletes, demotes or excludes a candidate, and it
never lets the with-lipase pool-external deferred layer (documented scale
1,206,655 pool-external hits, disposition ``deferred_structure_review``) enter
any count it reports.

Two measured defects were removed here (task F5 of
``docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign-followup.md``):

1. **Hard-coded candidate expectation.**  The module used to abort unless the
   triaged selection held exactly the original pilot's nineteen candidates,
   which made the tool fail closed for every other target set -- including the
   nPHAMCL panel that packet P4 needs.  The expectation is now an explicit
   caller contract:
   ``prepare(..., expected_candidates=N, expected_candidates_source="...")`` or
   ``--expected-candidates N``.  With no declared expectation the observed
   count is measured and recorded in ``candidate_count_contract`` and nothing
   is enforced.  A declared expectation that does not match still raises, and
   the message names the declared source instead of a bare literal.

2. **Asserted tool availability.**  The manifest used to declare
   ``structure_predictor.status = "pending_tool"`` with the reason "no predictor
   or structure-comparison executable is available".  That was never measured
   and directly contradicted the measured server fact that Foldseek is
   available and self-tested (see
   ``docs/T141_20260920_phaded_structure_verification_status.md`` section 1).
   Availability is now produced by :func:`probe_structure_tools`, a pure
   function over injected collaborators, which records per tool: status, the
   located path, the version string the executable prints for its own
   ``--version``/``-h``, the SHA-256 of the executable file (``pending`` when it
   cannot be computed -- never invented), the probe command, the lookup trace
   and a reason.  The manifest carries that measured report verbatim under
   ``structure_predictor.tools``.

The manifest keeps two independent booleans so that a probe outcome can never be
read as a computed result: ``tool_available`` says an executable was located and
answered its own version command, ``prediction_performed`` says whether any
prediction or comparison actually ran.  This script never runs one, so
``prediction_performed`` is always ``False`` here, and an available Foldseek is
never presented as a completed comparison.

Status vocabulary of ``structure_predictor.status``:

``available``
    At least one checked tool was located and answered its version probe.
``unavailable``
    The probe ran to completion and located no tool.  This is a measured
    statement about the runtime that built the manifest only; it is not a claim
    about any other environment (for example the server).
``pending_tool``
    Reserved for a probe that genuinely failed: the probe raised, a located
    executable did not answer its version command, the caller skipped the probe,
    or no probe record was produced.  ``reason`` then names the probe command
    and its failure.

Per-tool records additionally carry ``measurement_complete``.  Availability
requires a located path plus a version string the executable itself produced; a
hash that could not be computed leaves ``sha256`` at ``pending`` and
``measurement_complete`` false without downgrading availability, because the
hash is a provenance field and not the availability measurement.

Every pre-existing flag (``--decisions``, ``--candidates``, ``--reference``,
``--output-dir``), output filename and manifest key is preserved; new keys were
added rather than renamed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
from collections.abc import Callable, Iterable, Mapping, Sequence
from functools import partial
from pathlib import Path


KEEP = {"retain_dPHASCL1_like_candidate", "retain_extracellular_dPHASCL2_like_candidate"}

#: Disposition of the with-lipase pool-external deferred layer.  Rows carrying it
#: stay deferred, are never selected here and never enter any count.
DEFERRED_DISPOSITION = "deferred_structure_review"

#: Columns that may carry a primary disposition in a decisions table.
DISPOSITION_FIELDS = ("primary_disposition", "disposition", "candidate_disposition")

#: Tools whose availability is measured at runtime, with the executables to look
#: up and the flag each one accepts for printing its own version.
STRUCTURE_TOOLS: dict[str, dict[str, tuple[str, ...]]] = {
    "Foldseek": {"executables": ("foldseek",), "version_args": ("--version",)},
    "ColabFold": {"executables": ("colabfold_batch", "colabfold"), "version_args": ("--version",)},
    "AlphaFold": {"executables": ("run_alphafold.py", "run_alphafold.sh", "alphafold"), "version_args": ("--help",)},
    "ESMFold": {"executables": ("esm-fold", "esm_fold", "fold.py"), "version_args": ("--help",)},
    "Boltz": {"executables": ("boltz",), "version_args": ("--version",)},
}

PROBE_TIMEOUT_SECONDS = 20.0
PROBE_OUTPUT_LIMIT = 200

SHA256_HEX = re.compile(r"\A[0-9a-f]{64}\Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header = None
    chunks: list[str] = []
    for line in path.read_text(encoding="ascii").splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records[header.split()[0]] = (header, "".join(chunks))
            header, chunks = line[1:], []
        elif line:
            chunks.append(line)
    if header is not None:
        records[header.split()[0]] = (header, "".join(chunks))
    return records


def write_records(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


# ---------------------------------------------------------------------------
# Measured structure-tool availability probe
# ---------------------------------------------------------------------------


def _default_path_exists(path: Path) -> bool:
    try:
        return Path(path).is_file()
    except OSError:
        return False


def _run_command(command: Sequence[str], timeout: float = PROBE_TIMEOUT_SECONDS) -> tuple[int, str, str]:
    completed = subprocess.run(
        list(command), capture_output=True, text=True, timeout=timeout, check=False
    )
    return completed.returncode, completed.stdout or "", completed.stderr or ""


def _first_line(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped[:PROBE_OUTPUT_LIMIT]
    return ""


def _probe_one_tool(
    tool: str,
    spec: Mapping[str, Sequence[str]],
    *,
    which: Callable[[str], str | None],
    path_exists: Callable[[Path], bool],
    runner: Callable[[Sequence[str]], tuple[int, str, str]],
    sha256_file: Callable[[Path], str],
    declared_path: str | None,
) -> dict[str, object]:
    """Measure one tool.  Never raises and never guesses an unmeasured value."""
    executables = tuple(str(name) for name in spec.get("executables", ()))
    version_args = tuple(str(arg) for arg in spec.get("version_args", ("--version",)))
    record: dict[str, object] = {
        "tool": tool,
        "status": "pending",
        "path": None,
        "executable": None,
        "executables_checked": list(executables),
        "version": None,
        "version_command": None,
        "version_args": list(version_args),
        "sha256": "pending",
        "sha256_status": "not_attempted",
        "probe_command": "",
        "probe_lookup": "",
        "probe_exit_code": None,
        "probe_output": None,
        "measurement_complete": False,
        "reason": "",
    }
    lookups: list[str] = []
    lookup_errors: list[str] = []
    located: str | None = None

    if declared_path:
        lookups.append(f"declared path {declared_path}")
        try:
            found_declared = bool(path_exists(Path(declared_path)))
        except Exception as error:  # a broken collaborator is a probe failure
            found_declared = False
            lookup_errors.append(f"declared path probe raised {type(error).__name__}: {error}")
            lookups.append(lookup_errors[-1])
        if found_declared:
            located = str(declared_path)

    if located is None:
        for name in executables:
            lookups.append(f"which({name})")
            try:
                found_on_path = which(name)
            except Exception as error:
                found_on_path = None
                lookup_errors.append(f"which({name}) raised {type(error).__name__}: {error}")
                lookups.append(lookup_errors[-1])
            if found_on_path:
                located = str(found_on_path)
                record["executable"] = name
                break

    record["probe_lookup"] = "; ".join(lookups)

    if located is None:
        record["probe_command"] = f"which {executables[0]}" if executables else f"which {tool}"
        if lookup_errors:
            record["status"] = "pending"
            record["reason"] = "executable lookup failed: " + "; ".join(lookup_errors)
        else:
            record["status"] = "unavailable"
            record["reason"] = (
                f"no executable found for {tool}: looked up {', '.join(executables) or tool} "
                "on PATH and no declared --tool-path matched"
            )
        return record

    record["path"] = located
    probe_command = " ".join([located, *version_args])
    record["probe_command"] = probe_command
    record["version_command"] = probe_command

    hash_error: str | None = None
    try:
        digest = str(sha256_file(Path(located)))
    except Exception as error:
        hash_error = f"{type(error).__name__}: {error}"
    else:
        if SHA256_HEX.match(digest):
            record["sha256"] = digest
            record["sha256_status"] = "computed"
        else:
            hash_error = f"sha256 helper returned a non-digest value ({digest!r})"
    if hash_error is not None:
        record["sha256"] = "pending"
        record["sha256_status"] = "pending"

    try:
        exit_code, stdout, stderr = runner([located, *version_args])
    except Exception as error:
        record["status"] = "pending"
        record["reason"] = (
            f"located at {located} but probe command '{probe_command}' raised "
            f"{type(error).__name__}: {error}"
        )
        return record

    record["probe_exit_code"] = int(exit_code)
    output = _first_line(stdout) or _first_line(stderr)
    record["probe_output"] = output

    if int(exit_code) != 0:
        record["status"] = "pending"
        record["reason"] = (
            f"located at {located} but probe command '{probe_command}' exited {int(exit_code)}"
            f"; output {output!r}"
        )
        return record

    if not output:
        record["status"] = "pending"
        record["reason"] = (
            f"located at {located} but probe command '{probe_command}' exited 0 without a "
            "version string, so availability is not confirmed"
        )
        return record

    record["status"] = "available"
    record["version"] = output
    record["measurement_complete"] = hash_error is None
    reason = f"located at {located}; '{probe_command}' exited 0 with version {output!r}"
    if hash_error is not None:
        reason += f"; sha256 pending: {hash_error}"
    record["reason"] = reason
    return record


def probe_structure_tools(
    *,
    which: Callable[[str], str | None] | None = None,
    path_exists: Callable[[Path], bool] | None = None,
    runner: Callable[[Sequence[str]], tuple[int, str, str]] | None = None,
    sha256_file: Callable[[Path], str] | None = None,
    tool_paths: Mapping[str, str] | None = None,
    tool_specs: Mapping[str, Mapping[str, Sequence[str]]] | None = None,
    timeout: float = PROBE_TIMEOUT_SECONDS,
) -> dict[str, dict[str, object]]:
    """Measure which structure tools are available in *this* runtime.

    Every collaborator is injected so tests need no real tool: ``which`` looks an
    executable up (default ``shutil.which``), ``path_exists`` checks a declared
    path, ``runner`` executes the version command and returns
    ``(exit_code, stdout, stderr)``, ``sha256_file`` hashes the executable file.
    ``tool_paths`` maps a tool name to an explicitly declared executable path
    (for example the server Foldseek install), which is probed before PATH.

    Returns one record per tool keyed by tool name, with ``status`` equal to
    ``available`` (located and its own version command succeeded),
    ``unavailable`` (not found) or ``pending`` (located but not confirmable, or
    the lookup itself failed).  Nothing is asserted that was not measured: a
    hash that cannot be computed is recorded as ``pending``.
    """
    lookup = which if which is not None else shutil.which
    exists = path_exists if path_exists is not None else _default_path_exists
    execute = runner if runner is not None else partial(_run_command, timeout=timeout)
    digest = sha256_file if sha256_file is not None else sha256
    specs = tool_specs if tool_specs is not None else STRUCTURE_TOOLS
    declared = {str(name): str(value) for name, value in (tool_paths or {}).items()}

    report: dict[str, dict[str, object]] = {}
    for tool, spec in specs.items():
        report[str(tool)] = _probe_one_tool(
            str(tool),
            spec,
            which=lookup,
            path_exists=exists,
            runner=execute,
            sha256_file=digest,
            declared_path=declared.get(str(tool)),
        )
    return report


def summarize_tool_probe(
    report: Mapping[str, Mapping[str, object]] | None,
    *,
    failure: str | None = None,
    skipped_note: str | None = None,
) -> dict[str, object]:
    """Turn measured per-tool records into the ``structure_predictor`` block.

    ``tool_available`` and ``prediction_performed`` are independent: a located
    predictor is never reported as a performed prediction.  ``pending_tool`` is
    emitted only when the probe genuinely failed to measure anything.
    """
    tools = {str(name): dict(record) for name, record in (report or {}).items()}
    available = sorted(name for name, record in tools.items() if record.get("status") == "available")
    pending = sorted(name for name, record in tools.items() if record.get("status") == "pending")
    unavailable = sorted(
        name for name, record in tools.items() if record.get("status") == "unavailable"
    )
    unmeasured = sorted(set(tools) - set(available) - set(pending) - set(unavailable))
    probe_commands = sorted(
        {
            str(record.get("probe_command") or "")
            for record in tools.values()
            if record.get("probe_command")
        }
    )
    tool_available = bool(available)

    if tool_available:
        status = "available"
        reason = "measured at runtime: " + "; ".join(
            f"{name} at {tools[name].get('path')} "
            f"(version {tools[name].get('version')!r}, sha256 {tools[name].get('sha256')})"
            for name in available
        )
    elif failure is not None:
        status = "pending_tool"
        reason = (
            f"tool probe failed before a measured record existed: {failure}; "
            "no probe command completed"
        )
    elif skipped_note is not None:
        status = "pending_tool"
        reason = skipped_note
    elif pending or unmeasured:
        status = "pending_tool"
        reason = "tool probe did not complete: " + "; ".join(
            f"{name}: '{tools[name].get('probe_command')}' -> {tools[name].get('reason')}"
            for name in [*pending, *unmeasured]
        )
    elif not tools:
        status = "pending_tool"
        reason = "tool probe produced no records; no probe command was executed"
    else:
        status = "unavailable"
        reason = "probe ran to completion and measured no available tool: " + "; ".join(
            f"{name}: {tools[name].get('reason')}" for name in unavailable
        )

    probe_performed = failure is None and skipped_note is None and bool(tools)
    return {
        "status": status,
        "candidates_checked": sorted(tools),
        "tools": tools,
        "available_tools": available,
        "pending_tools": pending,
        "unavailable_tools": unavailable,
        "tool_available": tool_available,
        "prediction_performed": False,
        "prediction_evidence": None,
        "prediction_performed_note": (
            "This script prepares review inputs only. No structure prediction and no structure "
            "comparison was executed, so tool availability is not evidence of any computed "
            "structural result and must never be reported as one."
        ),
        "probe_command": "; ".join(probe_commands) if probe_commands else "not executed",
        "probe_environment": {
            "probe_performed": probe_performed,
            "platform": platform.platform(),
            "python": platform.python_version(),
            "executable": sys.executable,
            "path_lookup": "shutil.which",
            "scope": (
                "measurement applies to the runtime that built this manifest; it is not a claim "
                "about any other environment"
            ),
        },
        "reason": reason,
    }


# ---------------------------------------------------------------------------
# Preparation
# ---------------------------------------------------------------------------


def prepare(
    decisions: Path,
    candidates: Path,
    output_dir: Path,
    reference: Path | None = None,
    *,
    expected_candidates: int | None = None,
    expected_candidates_source: str | None = None,
    tool_probe: Callable[[], Mapping[str, Mapping[str, object]]] | None = None,
    probe_tools: bool = True,
) -> dict[str, object]:
    """Freeze structure-review inputs and record what was actually measured.

    ``expected_candidates`` is the caller's candidate-count contract; there is no
    built-in expectation.  When it is declared but does not match the observed
    selection the call raises and the message names
    ``expected_candidates_source``.  ``tool_probe`` overrides the runtime tool
    probe (tests inject it); ``probe_tools=False`` skips probing and records
    ``pending_tool`` plus the reason instead of claiming an availability.
    """
    if expected_candidates is not None:
        if isinstance(expected_candidates, bool) or not isinstance(expected_candidates, int):
            raise ValueError(
                f"expected candidate count must be a non-negative integer, observed {expected_candidates!r}"
            )
        if expected_candidates < 0:
            raise ValueError(
                f"expected candidate count must be a non-negative integer, observed {expected_candidates}"
            )
    if DEFERRED_DISPOSITION in KEEP:
        raise ValueError(f"{DEFERRED_DISPOSITION} must never be part of the selection set")

    all_rows = read_tsv(decisions)
    rows = [row for row in all_rows if row.get("phylo_decision") in KEEP]
    observed = len(rows)
    source = expected_candidates_source or (
        "caller-declared expectation" if expected_candidates is not None else None
    )
    if expected_candidates is not None and observed != expected_candidates:
        raise ValueError(
            "candidate count contract not satisfied: declared source "
            f"{source!r} expects {expected_candidates} candidate(s), observed {observed}"
        )

    def is_deferred(row: Mapping[str, str]) -> bool:
        return any(row.get(field) == DEFERRED_DISPOSITION for field in DISPOSITION_FIELDS)

    deferred_rows = [row for row in all_rows if is_deferred(row)]
    deferred_selected = [row for row in rows if is_deferred(row)]
    if deferred_selected:
        raise ValueError(
            f"{DEFERRED_DISPOSITION} rows must never enter the structure-review selection"
        )

    accessions = [row.get("source_accession", "") for row in rows]
    if not all(accessions) or len(set(accessions)) != len(accessions):
        raise ValueError("candidate accessions must be non-empty and unique")
    candidate_records = read_fasta(candidates)
    missing = [accession for accession in accessions if accession not in candidate_records]
    if missing:
        raise ValueError("candidate FASTA missing: " + ",".join(missing[:3]))
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = output_dir / "structure_candidates.faa"
    write_records(candidate_path, [(accession, candidate_records[accession][1]) for accession in accessions])
    selection_path = output_dir / "structure_candidates.tsv"
    fields = sorted({key for row in rows for key in row})
    with selection_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    refs_written = 0
    reference_path = output_dir / "nearest_reference_controls.faa"
    if reference is not None:
        ref_records = read_fasta(reference)
        by_accession = {header.split("|", 1)[-1]: (header, sequence) for header, sequence in ref_records.values()}
        controls = []
        for row in rows:
            accession = row.get("nearest_reference_accession", "")
            record = by_accession.get(accession)
            if record is not None:
                controls.append((f"reference|{record[0]}", record[1]))
        write_records(reference_path, controls)
        refs_written = len(controls)

    probe_failure: str | None = None
    skipped_note: str | None = None
    report: Mapping[str, Mapping[str, object]] = {}
    if probe_tools:
        target = tool_probe if tool_probe is not None else probe_structure_tools
        try:
            report = target()
        except Exception as error:
            report = {}
            probe_failure = f"{type(error).__name__}: {error}"
    else:
        skipped_note = (
            "tool availability probe skipped by the caller (--skip-tool-probe); no probe command "
            "was executed and no availability claim is made"
        )
    structure_predictor = summarize_tool_probe(report, failure=probe_failure, skipped_note=skipped_note)

    contract = {
        "observed_candidates": observed,
        "expected_candidates": expected_candidates,
        "expected_candidates_source": source,
        "enforced": expected_candidates is not None,
        "status": "satisfied" if expected_candidates is not None else "no_fixed_expectation",
        "note": (
            "No candidate count is hard-coded. Supply --expected-candidates (or the caller's own "
            "input contract) to enforce one; without it the observed count is only measured."
        ),
    }
    candidate_only = {
        "script_role": "prepare_structure_review_inputs_only",
        "deletes_candidates": False,
        "demotes_candidates": False,
        "excludes_candidates": False,
        "deferred_disposition": DEFERRED_DISPOSITION,
        "deferred_rows_observed": len(deferred_rows),
        "deferred_rows_selected": len(deferred_selected),
        "deferred_layer_included_in_counts": False,
        "note": (
            "The with-lipase pool-external deferred layer stays deferred_structure_review: it is "
            "never selected here, never counted, and never deleted, demoted or excluded."
        ),
    }

    manifest = {
        "status": "planned_not_run",
        "candidate_count": observed,
        "decision_counts": {decision: sum(row.get("phylo_decision") == decision for row in rows) for decision in sorted(KEEP)},
        "nearest_reference_controls_written": refs_written,
        "inputs": {"decisions": {"path": str(decisions.resolve()), "size": decisions.stat().st_size, "sha256": sha256(decisions)}, "candidate_fasta": {"path": str(candidates.resolve()), "size": candidates.stat().st_size, "sha256": sha256(candidates)}},
        "outputs": {"candidate_fasta": {"path": str(candidate_path.resolve()), "size": candidate_path.stat().st_size, "sha256": sha256(candidate_path)}, "selection_tsv": {"path": str(selection_path.resolve()), "size": selection_path.stat().st_size, "sha256": sha256(selection_path)}},
        "structure_predictor": structure_predictor,
        "phenotype_boundary": "Predicted or retrieved structures would provide candidate architecture evidence only; they cannot validate PHB/PHA degradation phenotype.",
        "candidate_count_contract": contract,
        "candidate_only_semantics": candidate_only,
        "tool_available": structure_predictor["tool_available"],
        "prediction_performed": structure_predictor["prediction_performed"],
        "structure_evidence_boundary": "Tool availability is measured, not asserted, and is never a computed structural result; no prediction or comparison was performed by this script.",
    }
    manifest_path = output_dir / "structure_review_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def parse_tool_paths(values: Iterable[str], tool_specs: Mapping[str, object] | None = None) -> dict[str, str]:
    """Parse repeated ``TOOL=PATH`` declarations into a mapping."""
    specs = tool_specs if tool_specs is not None else STRUCTURE_TOOLS
    known = {str(name).lower(): str(name) for name in specs}
    parsed: dict[str, str] = {}
    for value in values:
        tool, separator, path = str(value).partition("=")
        tool, path = tool.strip(), path.strip()
        if not separator or not tool or not path:
            raise ValueError(f"--tool-path expects TOOL=PATH, observed {value!r}")
        canonical = known.get(tool.lower())
        if canonical is None:
            raise ValueError(
                f"--tool-path names an unknown tool {tool!r}; known tools: {', '.join(sorted(known.values()))}"
            )
        parsed[canonical] = path
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decisions", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--expected-candidates",
        type=int,
        default=None,
        help=(
            "Declared number of high-information candidates. Optional: with no declared "
            "expectation the observed count is only measured and recorded."
        ),
    )
    parser.add_argument(
        "--expected-candidates-source",
        default=None,
        help="Name of the document/contract the declared expectation comes from (recorded in the manifest).",
    )
    parser.add_argument(
        "--tool-path",
        action="append",
        default=[],
        metavar="TOOL=PATH",
        help="Declared executable path for a structure tool, probed before the PATH lookup. Repeatable.",
    )
    parser.add_argument(
        "--skip-tool-probe",
        action="store_true",
        help=(
            "Do not probe tool availability. The manifest then records pending_tool with the reason "
            "that the probe was skipped; it never claims availability or unavailability."
        ),
    )
    parser.add_argument("--probe-timeout", type=float, default=PROBE_TIMEOUT_SECONDS)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    tool_paths = parse_tool_paths(args.tool_path)
    tool_probe = None
    if not args.skip_tool_probe:
        tool_probe = partial(
            probe_structure_tools, tool_paths=tool_paths, timeout=args.probe_timeout
        )
    expected_source = args.expected_candidates_source
    if expected_source is None and args.expected_candidates is not None:
        expected_source = "--expected-candidates"
    manifest = prepare(
        args.decisions,
        args.candidates,
        args.output_dir,
        args.reference,
        expected_candidates=args.expected_candidates,
        expected_candidates_source=expected_source,
        tool_probe=tool_probe,
        probe_tools=not args.skip_tool_probe,
    )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
