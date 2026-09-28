#!/usr/bin/env python3
"""Parse the bounded Cys targeted-recall HMM scores (2026-09-17).

The discovery-layer HMM (``cys_discovery_uncalibrated``) was scored against
the 100 existing scan shards of the frozen scan 13 input space.  This parser
rebuilds the recall tables and quantifies how many discovery hits fall
**outside** the frozen 109,087-protein candidate universe.

Load-bearing invariants (never relax):

* the discovery threshold is frozen at ``E < 1e-5`` and ``classify_hits``
  refuses any other value;
* an output row never carries a family call and never changes a
  ``subtype_call``; every row is labelled ``discovery_hmm_uncalibrated``.

Because the discovery HMM has no discriminative power, every recalled row is a
review-priority-queue entry, never specificity evidence.  Annotations are only
taken from existing in-pool tables; out-of-pool proteins are ``not_available``
and are never fabricated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

DISCOVERY_EVALUE_THRESHOLD = 1e-5
MODEL_LAYER = "discovery_hmm_uncalibrated"
RECALL_POOL_Z = 109087
NOT_AVAILABLE = "not_available"
DEFINITIONAL_ZERO = "definitional_zero_outside_task6_scope"
NO_DISCRIMINATION_NOTE = (
    "discovery layer has no discriminative power (443/563 x1 confounders hit); "
    "review priority only, never specificity evidence"
)

POOL_SCOPES = {"full_pool", "declared_shard_sample"}

RECALL_HIT_FIELDS = [
    "protein_id",
    "best_evalue",
    "best_score",
    "n_reported_rows",
    "in_existing_candidate_pool",
    "candidate_assignment_status",
    "candidate_prior_superfamily",
    "sbd_pf06850_binding_state",
    "lipase_box_state",
    "in_task6_tag_layer",
    "model_layer",
    "family_call_made",
    "subtype_call_impact",
]

NEW_CANDIDATE_FIELDS = [
    "protein_id",
    "best_evalue",
    "best_score",
    "n_reported_rows",
    "in_existing_candidate_pool",
    "task6_tag_layer_claimed",
    "task6_overlap_basis",
    "review_priority_only",
    "discriminative_power_caveat",
    "model_layer",
    "family_call_made",
    "new_family_call_made",
    "subtype_call_impact",
]

CALL_COLUMN_NAMES = {"family_call", "subtype_call", "new_family_call"}
BOUNDARY_DECLARATION_NAMES = {
    "family_call_made",
    "new_family_call_made",
    "subtype_call_impact",
    "subtype_call_rows_changed",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bound_file(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"path": str(path), "status": "pending", "sha256": None, "size": None}
    return {
        "path": str(path),
        "status": "verified",
        "sha256": sha256_file(path),
        "size": path.stat().st_size,
    }


def parse_tblout(path: Path) -> dict[str, dict[str, Any]]:
    """Return {target_name: {best_e, best_score, n_reported_rows}}.

    GTDB protein identifiers contain ``|`` as part of the identifier
    (``>GCA_xxx|contig_N``); the target name must NOT be split on ``|``.
    Comment lines (``#``) and rows with fewer than six whitespace-separated
    fields are skipped.
    """
    result: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            target = fields[0]
            try:
                e_value = float(fields[4])
                score = float(fields[5])
            except ValueError:
                continue
            row = result.get(target)
            if row is None:
                result[target] = {
                    "best_e": e_value,
                    "best_score": score,
                    "n_reported_rows": 1,
                }
            else:
                row["n_reported_rows"] += 1
                if e_value < row["best_e"]:
                    row["best_e"] = e_value
                    row["best_score"] = score
    return result


def classify_hits(
    parsed: Mapping[str, Mapping[str, Any]], threshold: float
) -> dict[str, list[str]]:
    """Split parsed targets into hits (best_e < threshold) and non-hits.

    The threshold must equal the pinned constant; any other value raises.
    """
    if threshold != DISCOVERY_EVALUE_THRESHOLD:
        raise ValueError(
            "discovery threshold is frozen at E < 1e-5; refusing threshold=%r" % threshold
        )
    hits = [t for t, r in parsed.items() if r["best_e"] < threshold]
    non_hits = [t for t, r in parsed.items() if r["best_e"] >= threshold]
    return {"hits": hits, "non_hits": non_hits}


def rescale_evalue(evalue: float, from_z: float, to_z: float) -> float:
    """Rescale an E-value linearly with database size Z (E is proportional to Z)."""
    if from_z <= 0 or to_z <= 0:
        raise ValueError("database sizes must be positive")
    return evalue * (to_z / from_z)


def _annotation_cells(
    protein_id: str,
    motif_rows: Mapping[str, Mapping[str, str]],
    in_pool: bool,
) -> dict[str, str]:
    if in_pool and protein_id in motif_rows:
        row = motif_rows[protein_id]
        return {
            "candidate_assignment_status": row.get(
                "candidate_assignment_status", NOT_AVAILABLE
            ),
            "candidate_prior_superfamily": row.get(
                "candidate_prior_superfamily", NOT_AVAILABLE
            ),
            "sbd_pf06850_binding_state": row.get(
                "sbd_pf06850_binding_state", NOT_AVAILABLE
            ),
            "lipase_box_state": row.get("lipase_box_state", NOT_AVAILABLE),
        }
    return {
        "candidate_assignment_status": NOT_AVAILABLE,
        "candidate_prior_superfamily": NOT_AVAILABLE,
        "sbd_pf06850_binding_state": NOT_AVAILABLE,
        "lipase_box_state": NOT_AVAILABLE,
    }


def build_hit_rows(
    pool_ids: Iterable[str],
    tblout_hits: Mapping[str, Mapping[str, Any]],
    candidate_universe: set[str],
    motif_rows: Mapping[str, Mapping[str, str]],
    tag_layer: set[str],
) -> list[dict[str, str]]:
    """One row per discovery hit (best_e < 1e-5), sorted by protein_id."""
    rows: list[dict[str, str]] = []
    for protein_id in sorted(set(pool_ids) & set(tblout_hits)):
        record = tblout_hits[protein_id]
        if record["best_e"] >= DISCOVERY_EVALUE_THRESHOLD:
            continue
        in_pool = protein_id in candidate_universe
        cells = _annotation_cells(protein_id, motif_rows, in_pool)
        rows.append(
            {
                "protein_id": protein_id,
                "best_evalue": "%.3g" % record["best_e"],
                "best_score": "%.1f" % record["best_score"],
                "n_reported_rows": str(record["n_reported_rows"]),
                "in_existing_candidate_pool": "true" if in_pool else "false",
                "candidate_assignment_status": cells["candidate_assignment_status"],
                "candidate_prior_superfamily": cells["candidate_prior_superfamily"],
                "sbd_pf06850_binding_state": cells["sbd_pf06850_binding_state"],
                "lipase_box_state": cells["lipase_box_state"],
                "in_task6_tag_layer": "true" if protein_id in tag_layer else "false",
                "model_layer": MODEL_LAYER,
                "family_call_made": "false",
                "subtype_call_impact": "none",
            }
        )
    return rows


def build_new_candidate_rows(hit_rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    """New candidates = hits outside the existing pool.

    The Task 6 tag layer only covers in-pool proteins, so the overlap of new
    candidates with the tag layer is definitionally zero.
    """
    rows: list[dict[str, str]] = []
    for hit in hit_rows:
        if hit["in_existing_candidate_pool"] != "false":
            continue
        rows.append(
            {
                "protein_id": hit["protein_id"],
                "best_evalue": hit["best_evalue"],
                "best_score": hit["best_score"],
                "n_reported_rows": hit["n_reported_rows"],
                "in_existing_candidate_pool": "false",
                "task6_tag_layer_claimed": "false",
                "task6_overlap_basis": DEFINITIONAL_ZERO,
                "review_priority_only": "true",
                "discriminative_power_caveat": NO_DISCRIMINATION_NOTE,
                "model_layer": MODEL_LAYER,
                "family_call_made": "false",
                "new_family_call_made": "false",
                "subtype_call_impact": "none",
            }
        )
    return rows


def build_summary(
    *,
    pool_size: int,
    pool_scope: str,
    sampling_statement: str,
    hits: int,
    new_candidates: int,
    in_pool_hits: int,
    tag_layer_size: int,
    new_candidate_tag_overlap: int,
    strict_hits_at_pool_z: int,
    pool_z: int,
    scored_from_annotation_table: int,
) -> dict[str, Any]:
    if pool_scope not in POOL_SCOPES:
        raise ValueError("unknown pool scope: %r" % pool_scope)
    if pool_scope == "declared_shard_sample" and not sampling_statement:
        raise ValueError("sampled scope requires a sampling statement")
    if strict_hits_at_pool_z > hits:
        raise ValueError("strict variant cannot exceed the reported hit count")
    if pool_size <= 0:
        raise ValueError("pool size must be positive")
    return {
        "model_layer": MODEL_LAYER,
        "pool_scope": pool_scope,
        "sampling_statement": sampling_statement,
        "sampling_declared": pool_scope == "declared_shard_sample",
        "is_full_library_rescan": False,
        "hmmsearch_database_size_Z": RECALL_POOL_Z,
        "actual_pool_size_z": pool_z,
        "counts": {
            "recall_pool_proteins_scored": pool_size,
            "recall_pool_hits": hits,
            "new_cys_candidates": new_candidates,
            "in_pool_hits": in_pool_hits,
            "strict_hits_at_actual_pool_z": strict_hits_at_pool_z,
            "tag_layer_size": tag_layer_size,
            "new_candidate_tag_overlap": new_candidate_tag_overlap,
            "scored_from_annotation_table": scored_from_annotation_table,
        },
        "hit_rate": hits / pool_size,
        "new_candidate_share_of_hits": new_candidates / hits if hits else 0.0,
        "family_call_made": False,
        "new_family_call_made": False,
        "subtype_call_rows_changed": 0,
        "calibration_gate_met": False,
        "calibration_gate_unchanged": True,
    }


def _read_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def _read_motif_rows(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        header_line = handle.readline()
        if not header_line:
            return rows
        header = header_line.rstrip("\n").split("\t")
        wanted = [
            "accession",
            "candidate_assignment_status",
            "candidate_prior_superfamily",
            "sbd_pf06850_binding_state",
            "lipase_box_state",
        ]
        missing = [c for c in wanted if c not in header]
        if missing:
            raise ValueError("motif table missing columns: %s" % ",".join(missing))
        indices = [header.index(c) for c in wanted]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            values = [fields[i] if i < len(fields) else "" for i in indices]
            rows[values[0]] = dict(zip(wanted, values))
    return rows


def _task6_tag_layer(motif_rows: Mapping[str, Mapping[str, str]]) -> set[str]:
    tag: set[str] = set()
    for accession, row in motif_rows.items():
        if (
            row.get("sbd_pf06850_binding_state") == "detected"
            and row.get("lipase_box_state") != "supported"
        ):
            tag.add(accession)
    return tag


def build_arg_namespace(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--pool-faa", type=Path, required=True)
    parser.add_argument("--tblout", type=Path, required=True)
    parser.add_argument("--candidate-faa", type=Path, required=True)
    parser.add_argument("--motif-tsv", type=Path, required=True)
    parser.add_argument("--pool-scope", required=True, choices=sorted(POOL_SCOPES))
    parser.add_argument("--actual-pool-size-z", type=int, required=True)
    parser.add_argument("--server-snapshot", type=Path, default=None)
    parser.add_argument("--tool", default="pending")
    parser.add_argument("--command", default="pending")
    parser.add_argument("--input-record", action="append", default=[])
    return parser.parse_args(argv)


def _write_tsv(path: Path, fields: list[str], rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(f, "")) for f in fields) + "\n")


def scan_output_call_columns(
    scores: Iterable[str],
    claims: Iterable[str],
) -> dict[str, Any]:
    fields = list(scores) + list(claims)
    call_columns = [f for f in fields if f in CALL_COLUMN_NAMES]
    boundary_declaration_columns = any(
        f in BOUNDARY_DECLARATION_NAMES for f in fields
    )
    return {
        "call_columns": call_columns,
        "boundary_declaration_columns": boundary_declaration_columns,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    run_dir: Path = args.run_dir
    results_dir = run_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    pool_ids = _read_fasta_ids(args.pool_faa)
    if not pool_ids:
        raise ValueError("recall pool is empty")
    candidate_universe = set(_read_fasta_ids(args.candidate_faa))
    parsed = parse_tblout(args.tblout)
    if not parsed:
        raise ValueError("no pool record is reported in the tblout")
    motif_rows = _read_motif_rows(args.motif_tsv)
    tag_layer = _task6_tag_layer(motif_rows)

    hit_rows = build_hit_rows(
        pool_ids=pool_ids,
        tblout_hits=parsed,
        candidate_universe=candidate_universe,
        motif_rows=motif_rows,
        tag_layer=tag_layer,
    )
    new_rows = build_new_candidate_rows(hit_rows)

    in_pool_hits = sum(1 for r in hit_rows if r["in_existing_candidate_pool"] == "true")
    new_candidates = len(new_rows)
    tag_overlap = sum(
        1 for r in new_rows if r["task6_tag_layer_claimed"] == "true"
    )

    summary = build_summary(
        pool_size=len(pool_ids),
        pool_scope=args.pool_scope,
        sampling_statement=getattr(args, "sampling_statement", ""),
        hits=len(hit_rows),
        new_candidates=new_candidates,
        in_pool_hits=in_pool_hits,
        tag_layer_size=len(tag_layer),
        new_candidate_tag_overlap=tag_overlap,
        strict_hits_at_pool_z=len(hit_rows),
        pool_z=args.actual_pool_size_z,
        scored_from_annotation_table=sum(
            1 for r in hit_rows if r["in_existing_candidate_pool"] == "true"
        ),
    )

    _write_tsv(results_dir / "recall_hits.tsv", RECALL_HIT_FIELDS, hit_rows)
    _write_tsv(results_dir / "new_cys_candidates.tsv", NEW_CANDIDATE_FIELDS, new_rows)

    tool_records: dict[str, str] = {}
    for record in args.input_record:
        if "=" in record:
            key, value = record.split("=", 1)
            tool_records[key.strip()] = value.strip()

    manifest: dict[str, Any] = {
        "schema_version": "1.0",
        "run_id": run_dir.name,
        "model_layer": MODEL_LAYER,
        "registry_modified": False,
        "family_call_made": False,
        "new_family_call_made": False,
        "subtype_call_rows_changed": 0,
        "tool": args.tool,
        "command": args.command,
        "threads": {"max_single_task_threads_allowed": 40},
        "inputs": {
            "recall_pool": _bound_file(args.pool_faa),
            "tblout": _bound_file(args.tblout),
            "candidate_universe": _bound_file(args.candidate_faa),
            "motif_table": _bound_file(args.motif_tsv),
            "server_load_snapshot": (
                _bound_file(args.server_snapshot) if args.server_snapshot else None
            ),
        },
        "gtdb": {
            "taxonomy": {"path": None, "status": "pending", "sha256": None, "size": None},
            "metadata": {"path": None, "status": "pending", "sha256": None, "size": None},
            "tree": {"path": None, "status": "pending", "sha256": None, "size": None},
        },
    }
    for key, value in tool_records.items():
        manifest["inputs"][key] = _bound_file(Path(value))

    (results_dir / "recall_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (results_dir / "recall_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {"counts": summary["counts"], "summary": summary, "manifest": manifest}


def main(argv: list[str] | None = None) -> int:
    args = build_arg_namespace(argv)
    payload = run(args)
    print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
