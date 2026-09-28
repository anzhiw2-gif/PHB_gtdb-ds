#!/usr/bin/env python3
"""Score the bound PhaDED candidate profiles, keyed on ``model_layer`` (F3).

Model layer semantics (AGENTS.md, evidence-model redesign Task 5)
-----------------------------------------------------------------
Every profile carries two *independent* fields:

* ``model_layer`` in {:data:`REFERENCE_QUERY_ONLY`, :data:`DISCOVERY_LAYER`,
  :data:`SEQUENCE_FAMILY_LAYER`, :data:`CALIBRATED_LAYER`} — what the model *is*;
* ``functional_calibration_status`` — whether the *function* has been calibrated.
  It is never inferred from ``model_layer`` and vice versa; only an explicitly
  authorized promotion produces ``calibrated_candidate_model``.

Only :data:`DISCRIMINATING_MODEL_LAYERS` (``sequence_family_hmm_validated`` and
``calibrated_candidate_model``) may discriminate between candidates.  A
``discovery_hmm_uncalibrated`` profile is **recall-only**: it may score existing
intermediates, but it must never filter, drop, demote or exclude a candidate and
must never produce a family call (measured sensitivity on the 563 x1 confounder
hits is 78.69% — no discriminating power).

Legacy ``model_status`` fallback (documented, pinned by test)
-------------------------------------------------------------
``model_status`` is still written for backward compatibility, where it reads
``trained`` both for validated family HMMs and for discovery-layer HMMs.  This
script therefore never consults it for a decision; it is only carried through to
the audit table verbatim.  When the manifest **predates** ``model_layer`` (the
column is absent), the deliberately conservative default is applied: an absent
layer is read as :data:`REFERENCE_QUERY_ONLY` — *non-discriminating* — and every
such profile is reported as unscoreable instead of being assumed validated.
That is fail-closed on purpose: a frozen 2026-09 manifest whose rows say
``trained`` will score 0 validated profiles and report why, rather than silently
promoting nine recall-only/discovery HMMs into discriminating family models.
The fallback is stated in ``hmmsearch_provenance.json`` (``layer_source``,
``fallback_note``), never applied silently.

Outputs
-------
* ``results/phaded_profile_scores.tsv`` — discriminating layers only;
* ``results/phaded_profile_scores_recall_only.tsv`` — optional
  (``--include-discovery-layer``) discovery-layer scores, explicitly labelled
  ``recall_only=true`` with an empty ``family_call``;
* ``results/model_layer_audit.tsv`` — every manifest profile with its layer, the
  legacy word, whether it was discriminating and how it was used;
* ``results/hmmsearch_provenance.json`` — layer counts, the fallback note and the
  unscoreable profiles.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

#: The governed model-layer vocabulary (four layers, fixed by AGENTS.md).
REFERENCE_QUERY_ONLY = "reference_query_only"
DISCOVERY_LAYER = "discovery_hmm_uncalibrated"
SEQUENCE_FAMILY_LAYER = "sequence_family_hmm_validated"
CALIBRATED_LAYER = "calibrated_candidate_model"
MODEL_LAYERS = frozenset({
    REFERENCE_QUERY_ONLY,
    DISCOVERY_LAYER,
    SEQUENCE_FAMILY_LAYER,
    CALIBRATED_LAYER,
})

#: The only layers that may assign a sequence family / discriminate a candidate.
DISCRIMINATING_MODEL_LAYERS = frozenset({SEQUENCE_FAMILY_LAYER, CALIBRATED_LAYER})
#: Layers that may only recall: they never filter, drop, demote or call a family.
NON_DISCRIMINATING_MODEL_LAYERS = frozenset({REFERENCE_QUERY_ONLY, DISCOVERY_LAYER})

assert DISCRIMINATING_MODEL_LAYERS | NON_DISCRIMINATING_MODEL_LAYERS == MODEL_LAYERS, (
    "the discriminating/non-discriminating partition must cover the model-layer vocabulary"
)
assert not (DISCRIMINATING_MODEL_LAYERS & NON_DISCRIMINATING_MODEL_LAYERS), (
    "a model layer cannot be both discriminating and recall-only"
)

#: Layer assumed when a manifest has no ``model_layer`` column at all.  Absent
#: layer => non-discriminating: "legacy trained" is never read as "validated".
ABSENT_LAYER_FALLBACK = REFERENCE_QUERY_ONLY
ABSENT_LAYER_NOTE = (
    "profile_manifest.tsv has no model_layer column: every profile is treated as "
    f"'{ABSENT_LAYER_FALLBACK}' (non-discriminating) and reported as unscoreable. "
    "A legacy model_status value of 'trained' is never assumed validated; re-emit the "
    "manifest with model_layer to score validated sequence-family models."
)
LAYER_SOURCE_COLUMN = "model_layer_column"
LAYER_SOURCE_ABSENT = "absent_model_layer_column"

UNSCOREABLE_REASON_REFERENCE_ONLY = "reference_query_only_has_no_hmm"

USAGE_VALIDATED = "validated_sequence_family"
USAGE_RECALL_ONLY = "recall_only"
USAGE_REJECTED = "rejected_unscoreable"

RECALL_ONLY_TABLE = "phaded_profile_scores_recall_only.tsv"
AUDIT_TABLE = "model_layer_audit.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# model_layer resolution
# ---------------------------------------------------------------------------

def resolve_model_layer(row: dict[str, str]) -> str:
    """Return the row's declared ``model_layer``, fail-closed on anything else.

    ``model_status`` is *not* a fallback: it is deliberately ignored, because it
    reads ``trained`` for discovery-layer HMMs as well.  A row without a declared
    layer (legacy manifest) resolves to :data:`ABSENT_LAYER_FALLBACK`.
    """
    declared = (row.get("model_layer") or "").strip()
    if not declared:
        return ABSENT_LAYER_FALLBACK
    if declared not in MODEL_LAYERS:
        profile_id = (row.get("profile_id") or "").strip() or "<unknown profile>"
        raise ValueError(
            f"profile {profile_id!r} declares model_layer={declared!r}, which is not a "
            f"recognized model layer; expected one of {sorted(MODEL_LAYERS)}"
        )
    return declared


def is_discriminating(layer: str) -> bool:
    """True only for layers allowed to assign a sequence family."""
    assert layer in MODEL_LAYERS, f"unresolved model layer {layer!r} reached the decision"
    return layer in DISCRIMINATING_MODEL_LAYERS


def read_profile_manifest(manifest: Path) -> "tuple[list[dict[str, str]], dict[str, object]]":
    """Read the profile manifest, resolving ``model_layer`` for every row.

    Returns ``(rows, report)``.  ``rows`` gain a resolved ``model_layer`` plus the
    bookkeeping fields ``layer_source`` and ``functional_calibration_status`` (the
    latter carried verbatim, ``unknown`` when the column is absent — never
    inferred from the layer).  ``report`` documents where the layer came from so
    the fallback can be read straight out of the run provenance.
    """
    with manifest.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        raw_rows = list(reader)
    column_present = "model_layer" in fieldnames
    layer_source = LAYER_SOURCE_COLUMN if column_present else LAYER_SOURCE_ABSENT

    rows: list[dict[str, str]] = []
    counts: dict[str, int] = {layer: 0 for layer in sorted(MODEL_LAYERS)}
    for raw in raw_rows:
        layer = resolve_model_layer(raw)
        row = dict(raw)
        row["model_layer"] = layer
        row["layer_source"] = layer_source
        row["functional_calibration_status"] = (
            (raw.get("functional_calibration_status") or "").strip() or "unknown"
        )
        counts[layer] += 1
        rows.append(row)

    report: dict[str, object] = {
        "model_layer_column_present": column_present,
        "layer_source": layer_source,
        "fallback_layer": "" if column_present else ABSENT_LAYER_FALLBACK,
        "fallback_note": "" if column_present else ABSENT_LAYER_NOTE,
        "model_layer_counts": {layer: n for layer, n in counts.items() if n},
        "manifest_columns": fieldnames,
    }
    return rows, report


def select_scoreable_profiles(
    rows: "list[dict[str, str]]", *, include_discovery_layer: bool = False
) -> dict[str, list[dict[str, str]]]:
    """Split manifest rows into validated / recall-only / unscoreable.

    ``validated`` is the only group that may feed a discriminating score table;
    ``discovery`` is recall-only and is kept in a separate, labelled table;
    ``rejected`` holds profiles that cannot be scored at all (no HMM).
    """
    selection: dict[str, list[dict[str, str]]] = {
        "validated": [], "discovery": [], "rejected": [],
    }
    for row in rows:
        layer = row["model_layer"]
        if layer in DISCRIMINATING_MODEL_LAYERS:
            selection["validated"].append(row)
        elif layer == DISCOVERY_LAYER:
            # Recall-only: scored (optionally) but never discriminating.
            if include_discovery_layer:
                selection["discovery"].append(row)
        else:
            rejected = dict(row)
            rejected["unscoreable_reason"] = UNSCOREABLE_REASON_REFERENCE_ONLY
            selection["rejected"].append(rejected)
    return selection


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

def _default_runner(command, capture_output: bool = True, text: bool = True, check: bool = True):
    return subprocess.run(command, capture_output=capture_output, text=text, check=check)


def _score_profile(
    row: dict[str, str],
    *,
    inputs: Path,
    results: Path,
    logs: Path,
    fasta: Path,
    cpu: int,
    hmmsearch: str,
    runner,
) -> "tuple[list[dict[str, object]], dict[str, object]]":
    profile_id = row["profile_id"]
    layer = row["model_layer"]
    hmm = inputs / f"{profile_id}.hmm"
    tblout = results / f"{profile_id}.tblout"
    stdout_log = logs / f"{profile_id}.stdout.log"
    stderr_log = logs / f"{profile_id}.stderr.log"
    command = [hmmsearch, "--noali", "--cpu", str(cpu), "--tblout", str(tblout),
               "-E", "1e6", "--incE", "1e6", str(hmm), str(fasta)]
    completed = runner(command, capture_output=True, text=True, check=True)
    stdout_log.write_text(completed.stdout or "", encoding="utf-8")
    stderr_log.write_text(completed.stderr or "", encoding="utf-8")
    rows: list[dict[str, object]] = []
    hits = 0
    with tblout.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(f"malformed tblout row in {tblout}: {line.rstrip()}")
            row_out: dict[str, object] = {
                "accession": fields[0], "profile_id": profile_id, "score": float(fields[5]),
            }
            if layer == DISCOVERY_LAYER:
                # Recall-only rows are labelled and can never carry a family call.
                row_out.update({
                    "model_layer": layer,
                    "recall_only": "true",
                    "family_call": "",
                })
            rows.append(row_out)
            hits += 1
    command_record = {
        "profile_id": profile_id,
        "model_layer": layer,
        "command": command,
        "cpu": cpu,
        "hmm": {"path": str(hmm), "size": hmm.stat().st_size, "sha256": sha256(hmm)},
        "tblout": {"path": str(tblout), "size": tblout.stat().st_size, "sha256": sha256(tblout)},
        "hits": hits,
    }
    return rows, command_record


def _write_tsv(path: Path, fieldnames: "list[str]", rows: "list[dict[str, object]]") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _audit_rows(
    selection: "dict[str, list[dict[str, str]]]",
    *,
    include_discovery_layer: bool,
) -> "list[dict[str, str]]":
    usage: dict[str, str] = {}
    for row in selection["validated"]:
        usage[row["profile_id"]] = USAGE_VALIDATED
    for row in selection["discovery"]:
        usage[row["profile_id"]] = USAGE_RECALL_ONLY
    for row in selection["rejected"]:
        usage[row["profile_id"]] = USAGE_REJECTED
    ordered = selection["validated"] + selection["discovery"] + selection["rejected"]
    out: list[dict[str, str]] = []
    for row in ordered:
        layer = row["model_layer"]
        if layer == DISCOVERY_LAYER and not include_discovery_layer:
            # Not scored this run, but still recall-only by construction.
            usage.setdefault(row["profile_id"], USAGE_RECALL_ONLY)
        out.append({
            "profile_id": row["profile_id"],
            "profile_kind": (row.get("profile_kind") or "").strip(),
            "model_layer": layer,
            "layer_source": row["layer_source"],
            # Legacy column, carried verbatim for audit only — never a decision.
            "model_status": (row.get("model_status") or "").strip(),
            "functional_calibration_status": row["functional_calibration_status"],
            "discriminating": "true" if is_discriminating(layer) else "false",
            "usage": usage[row["profile_id"]],
            "unscoreable_reason": row.get("unscoreable_reason", ""),
        })
    return out


def _tool_version(executor, hmmsearch: str) -> str:
    """Version banner of the bound binary (``-h``); blank when unparseable."""
    banner = executor([hmmsearch, "-h"], capture_output=True, text=True, check=True).stdout or ""
    for line in banner.splitlines():
        if line.startswith("#"):
            return line.strip()
    return ""


def run(
    deploy_dir: Path,
    run_dir: Path,
    *,
    cpu: int,
    hmmsearch: str,
    runner=None,
    include_discovery_layer: bool = False,
) -> dict[str, object]:
    if cpu < 1 or cpu > 40:
        raise ValueError("cpu must be between 1 and 40")
    executor = runner if runner is not None else _default_runner
    inputs = deploy_dir / "inputs"
    results = run_dir / "results" / "hmmsearch"
    logs = run_dir / "logs"
    results.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    manifest_path = inputs / "profile_manifest.tsv"
    fasta = inputs / "candidate_union.faa"

    manifest_rows, layer_report = read_profile_manifest(manifest_path)
    selection = select_scoreable_profiles(
        manifest_rows, include_discovery_layer=include_discovery_layer
    )
    validated = selection["validated"]
    discovery = selection["discovery"]
    rejected = selection["rejected"]

    score_rows: list[dict[str, object]] = []
    recall_rows: list[dict[str, object]] = []
    commands: list[dict[str, object]] = []
    for row in validated + discovery:
        rows_out, command_record = _score_profile(
            row, inputs=inputs, results=results, logs=logs, fasta=fasta,
            cpu=cpu, hmmsearch=hmmsearch, runner=executor,
        )
        commands.append(command_record)
        if row["model_layer"] == DISCOVERY_LAYER:
            recall_rows.extend(rows_out)
        else:
            score_rows.extend(rows_out)

    score_path = run_dir / "results" / "phaded_profile_scores.tsv"
    _write_tsv(score_path, ["accession", "profile_id", "score"], score_rows)
    recall_path = run_dir / "results" / RECALL_ONLY_TABLE
    _write_tsv(recall_path, ["accession", "profile_id", "score", "model_layer",
                             "recall_only", "family_call"], recall_rows)
    audit_path = run_dir / "results" / AUDIT_TABLE
    _write_tsv(audit_path, ["profile_id", "profile_kind", "model_layer", "layer_source",
                            "model_status", "functional_calibration_status",
                            "discriminating", "usage", "unscoreable_reason"],
               _audit_rows(selection, include_discovery_layer=include_discovery_layer))

    unscoreable = [
        {"profile_id": row["profile_id"], "model_layer": row["model_layer"],
         "reason": row["unscoreable_reason"]}
        for row in rejected
    ]
    provenance = {
        "schema_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "hmmsearch": {"path": hmmsearch, "version": _tool_version(executor, hmmsearch)},
        "cpu": cpu,
        "decision_field": "model_layer",
        "legacy_model_status_is_a_decision_input": False,
        "model_layer_column_present": layer_report["model_layer_column_present"],
        "layer_source": layer_report["layer_source"],
        "fallback_layer": layer_report["fallback_layer"],
        "fallback_note": layer_report["fallback_note"],
        "model_layer_counts": layer_report["model_layer_counts"],
        "include_discovery_layer": include_discovery_layer,
        "validated_profile_count": len(validated),
        "discovery_profile_count": len(discovery),
        "unscoreable_profiles": unscoreable,
        "candidate_fasta": {"path": str(fasta), "size": fasta.stat().st_size, "sha256": sha256(fasta)},
        "profile_manifest": {"path": str(manifest_path), "size": manifest_path.stat().st_size,
                             "sha256": sha256(manifest_path)},
        "commands": commands,
        "score_table": {"path": str(score_path), "size": score_path.stat().st_size,
                        "sha256": sha256(score_path), "rows": len(score_rows)},
        "recall_only_score_table": {"path": str(recall_path), "size": recall_path.stat().st_size,
                                    "sha256": sha256(recall_path), "rows": len(recall_rows),
                                    "discriminating": False},
        "model_layer_audit": {"path": str(audit_path), "size": audit_path.stat().st_size,
                              "sha256": sha256(audit_path)},
    }
    (run_dir / "results" / "hmmsearch_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return {
        "validated_profiles": len(validated),
        "discovery_profiles": len(discovery),
        "unscoreable_profiles": len(unscoreable),
        "validated_score_rows": len(score_rows),
        "recall_only_score_rows": len(recall_rows),
        "model_layer_column_present": layer_report["model_layer_column_present"],
        "layer_source": layer_report["layer_source"],
        "fallback_note": layer_report["fallback_note"],
    }


def main(argv: "list[str] | None" = None, *, runner=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--deploy-dir", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=40)
    parser.add_argument("--hmmsearch", default="${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb/bin/hmmsearch")
    parser.add_argument(
        "--include-discovery-layer", action="store_true",
        help="additionally score discovery_hmm_uncalibrated profiles into a separate "
             "recall-only table; they never enter the discriminating score table and "
             "never filter, drop or demote a candidate",
    )
    args = parser.parse_args(argv)
    print(json.dumps(run(args.deploy_dir, args.run_dir, cpu=args.cpu, hmmsearch=args.hmmsearch,
                         runner=runner, include_discovery_layer=args.include_discovery_layer),
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
