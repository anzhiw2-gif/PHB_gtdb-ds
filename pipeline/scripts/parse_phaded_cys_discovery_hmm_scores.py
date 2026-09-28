#!/usr/bin/env python3
"""Parse Cys discovery-layer HMM scores and evaluate the pre-registered gates
(2026-09-17, Step 7).

Consumes the archived HMMER ``tblout`` files, rebuilds the 109,087-row score
table, derives the unassigned-claim list, evaluates G1-G4 exactly as
pre-registered, and writes:

* ``results/cys_discovery_scores.tsv``   - one row per candidate (109,087)
* ``results/cys_discovery_claims.tsv``   - hits that are not already assigned
* ``results/gate_evaluation.json``       - G1-G4 measured outcome (fail-closed)
* ``results/discovery_manifest.json``    - hashes, tool versions, commands,
  threads and the server load snapshot

The pre-registered discovery threshold (E-value < 1e-5) is a module constant and
``classify_hits`` refuses any other value: the threshold can never be relaxed
after a failed gate.

Boundary: hits are ``discovery_hmm_uncalibrated`` candidate-only evidence of
candidate homology or functional potential.  No family call is produced and no
``subtype_call`` row is changed.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path

MODEL_LAYER = "discovery_hmm_uncalibrated"
DISCOVERY_EVALUE_THRESHOLD = 1e-5
BOUNDARY = "candidate_only"
BOUNDARY_TEXT = (
    "candidate_only discovery-layer homology evidence from an uncalibrated HMM; "
    "not a family judgement and not a validated PHB/PHA degradation phenotype"
)

CYS_SUPERFAMILY = "intracellular nPHASCL without lipase box"
WITH_LIPASE_BOX_SUPERFAMILY = "intracellular nPHASCL with lipase box"
EXTRACELLULAR_SUPERFAMILIES = (
    "extracellular dPHASCL type 1",
    "extracellular dPHASCL type 2",
    "extracellular dPHAMCL",
    "extracellular native-SCL/PhaZ7-like",
)
MCC_CONTROL_ACCESSION = "Q84C08"

PRIMARY_SCORE_SOURCE = "pre_registered_threshold_run"
SUPPLEMENTARY_SCORE_SOURCE = "supplementary_unfiltered_run"
UNREPORTED_E = "not_reported"

SCORE_FIELDS = [
    "accession",
    "best_E",
    "best_score",
    "discovery_hit",
    "n_reported_rows",
    "score_source",
    "model_layer",
    "family_call_made",
    "subtype_call_impact",
]
CLAIM_FIELDS = [
    "accession",
    "best_E",
    "best_score",
    "candidate_prior_superfamily",
    "candidate_assignment_status",
    "unassigned_or_no_prior",
    "task6_tag_layer_claimed",
    "task6_tag_layer_rule",
    "claim_layer",
    "model_layer",
    "family_call_made",
    "new_family_call_made",
    "subtype_call_impact",
    "boundary_text",
]
FORBIDDEN_OUTPUT_COLUMN_TOKENS = ("family_call", "subtype_call")


# --------------------------------------------------------------------------- io
def sha256_file(path: os.PathLike[str] | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def read_tsv(path: os.PathLike[str] | str) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: os.PathLike[str] | str, fieldnames: list[str], rows: list[dict]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


# -------------------------------------------------------------------- parsing
def parse_tblout(path: os.PathLike[str] | str, *, alias_separator: str = "|") -> dict[str, dict]:
    """Return a hit table keyed by target name, plus ``|``-split aliases.

    HMMER reports the *whole* FASTA header up to the first whitespace as the
    target name.  The frozen reference FASTA uses ``>reference_id|accession``
    headers while every downstream table is keyed by ``accession``, so a lookup
    that only used the raw name would silently find nothing.  Each target name
    containing ``alias_separator`` is therefore additionally indexed under each
    of its fields (aliases share one entry object, so domain counts are not
    duplicated).
    """
    hits: dict[str, dict] = {}
    with open(path, encoding="utf-8", newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            target = fields[0]
            try:
                evalue = float(fields[4])
                score = float(fields[5])
            except ValueError:
                continue
            entry = hits.get(target)
            if entry is None:
                hits[target] = {"best_e": evalue, "best_score": score, "n_reported_rows": 1}
            else:
                entry["n_reported_rows"] += 1
                if evalue < entry["best_e"]:
                    entry["best_e"] = evalue
                    entry["best_score"] = score
    if alias_separator:
        for target in [k for k in hits if alias_separator in k]:
            entry = hits[target]
            for alias in target.split(alias_separator):
                if alias and alias not in hits:
                    hits[alias] = entry
    return hits


def classify_hits(hits: dict[str, dict], *, threshold: float) -> dict[str, dict]:
    """Split a parsed tblout into hits and non-hits at a pinned threshold."""
    if threshold != DISCOVERY_EVALUE_THRESHOLD:
        raise ValueError(
            "the pre-registered discovery threshold is frozen at %r; refusing threshold %r"
            % (DISCOVERY_EVALUE_THRESHOLD, threshold)
        )
    passing: dict[str, dict] = {}
    failing: dict[str, dict] = {}
    for target, entry in hits.items():
        (passing if entry["best_e"] < threshold else failing)[target] = entry
    return {"hits": passing, "non_hits": failing}


def read_candidate_accessions(path: os.PathLike[str] | str) -> list[str]:
    accessions: list[str] = []
    with open(path, encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith(">"):
                accessions.append(line[1:].strip().split()[0])
    return accessions


def build_score_rows(
    candidates: list[str],
    primary: dict[str, dict],
    supplementary: dict[str, dict],
) -> list[dict]:
    if not candidates:
        raise ValueError("candidate accession list is empty")
    classified = classify_hits(primary, threshold=DISCOVERY_EVALUE_THRESHOLD)
    primary_hits = classified["hits"]
    rows: list[dict] = []
    missing: list[str] = []
    for accession in candidates:
        in_primary = accession in primary_hits
        source = None
        entry = None
        if in_primary:
            entry = primary_hits[accession]
            source = PRIMARY_SCORE_SOURCE
        elif accession in supplementary:
            entry = supplementary[accession]
            source = SUPPLEMENTARY_SCORE_SOURCE
        if entry is None:
            missing.append(accession)
            rows.append({
                "accession": accession,
                "best_E": UNREPORTED_E,
                "best_score": UNREPORTED_E,
                "discovery_hit": "false",
                "n_reported_rows": "",
                "score_source": UNREPORTED_E,
                "model_layer": MODEL_LAYER,
                "family_call_made": "false",
                "subtype_call_impact": "none",
            })
            continue
        rows.append({
            "accession": accession,
            "best_E": format_evalue(entry["best_e"]),
            "best_score": format_score(entry["best_score"]),
            "discovery_hit": "true" if in_primary else "false",
            "n_reported_rows": entry.get("n_reported_rows", ""),
            "score_source": source,
            "model_layer": MODEL_LAYER,
            "family_call_made": "false",
            "subtype_call_impact": "none",
        })
    if missing and len(missing) == len(candidates):
        raise ValueError(
            "no score available for any candidate (missing %d of %d)"
            % (len(missing), len(candidates))
        )
    return rows


def format_evalue(value: float) -> str:
    if value == 0.0:
        return "0.0"
    return "%.3g" % value


def format_score(value: float) -> str:
    return "%.1f" % value


# ------------------------------------------------------------- gate evaluation
def _accessions_by_superfamily(ledger_rows: list[dict]) -> dict[str, set[str]]:
    buckets: dict[str, set[str]] = {}
    for row in ledger_rows:
        buckets.setdefault(row.get("phaded_superfamily", ""), set()).add(row.get("accession", ""))
    return buckets


def _hits_in(population: set[str], *sources: dict[str, dict]) -> list[dict]:
    out: list[dict] = []
    for accession in sorted(population):
        for source in sources:
            entry = source.get(accession)
            if entry and entry.get("best_e") is not None and entry["best_e"] < DISCOVERY_EVALUE_THRESHOLD:
                out.append({
                    "accession": accession,
                    "best_E": format_evalue(entry["best_e"]),
                    "best_score": format_score(entry["best_score"]),
                })
                break
    return out


def evaluate_g1(
    *,
    reference_hits: dict[str, dict],
    ledger_rows: list[dict],
    detected_controls: dict[str, dict],
) -> dict:
    buckets = _accessions_by_superfamily(ledger_rows)
    extracellular: set[str] = set()
    for superfamily in EXTRACELLULAR_SUPERFAMILIES:
        extracellular |= buckets.get(superfamily, set())
    lipase_box_superfamily = buckets.get(WITH_LIPASE_BOX_SUPERFAMILY, set())

    control_population = {MCC_CONTROL_ACCESSION}
    population = extracellular | control_population

    unexplained = _hits_in(population, reference_hits, detected_controls)
    cross_talk = _hits_in(lipase_box_superfamily, reference_hits, detected_controls)
    for entry in cross_talk:
        entry["classification"] = "explained_cross_talk"
        entry["explanation"] = (
            "adjacent intracellular superfamily that does carry a lipase box; a hit here is "
            "expected cross-talk and is disclosed, not counted against G1"
        )
    return {
        "gate_id": "G1",
        "gate_name": "specificity",
        "threshold": "unexplained_hits == 0",
        "population_size": len(population),
        "extracellular_reference_count": len(extracellular),
        "mcc_control_count": len(control_population),
        "unexplained_hits": len(unexplained),
        "unexplained_hit_list": unexplained,
        "explained_cross_talk_population_size": len(lipase_box_superfamily),
        "explained_cross_talk_count": len(cross_talk),
        "explained_cross_talk": cross_talk,
        "explained_cross_talk_counted_as_gate_failure": "false",
        "outcome": "passed" if not unexplained else "failed",
    }


def evaluate_g2(
    *, reference_hits: dict[str, dict], training_accessions: list[str], threshold: float
) -> dict:
    denominators = [a for a in training_accessions if a]
    called = [
        a for a in denominators
        if reference_hits.get(a) and reference_hits[a]["best_e"] < DISCOVERY_EVALUE_THRESHOLD
    ]
    recall = (len(called) / len(denominators)) if denominators else 0.0
    missed = sorted(set(denominators) - set(called))
    return {
        "gate_id": "G2",
        "gate_name": "sensitivity",
        "threshold": "recall >= %.2f" % threshold,
        "denominator": len(denominators),
        "recalled": len(called),
        "recall": recall,
        "missed_accessions": missed,
        "outcome": "passed" if recall >= threshold else "failed",
        "measurement_scope": "in_sample_not_held_out",
        "measurement_scope_note": (
            "the denominator is the integrity-filter-passing training alignment membership; the "
            "frozen ledger contains no held-out Cys positive, so this is a self-consistency "
            "measurement and NOT a generalisation estimate"
        ),
    }


def evaluate_g3(
    *,
    before_sha256: str,
    after_sha256: str,
    labelled_outputs: list[str],
    unlabelled_outputs: list[str],
) -> dict:
    unchanged = bool(before_sha256) and before_sha256 == after_sha256
    labelled = not unlabelled_outputs
    return {
        "gate_id": "G3",
        "gate_name": "stratification",
        "threshold": "sha256_before == sha256_after AND unlabelled_outputs == 0",
        "formal_scan_models_tsv_sha256_before": before_sha256,
        "formal_scan_models_tsv_sha256_after": after_sha256,
        "formal_scan_models_tsv_unchanged": unchanged,
        "labelled_output_count": len(labelled_outputs),
        "unlabelled_outputs": list(unlabelled_outputs),
        "model_layer": MODEL_LAYER,
        "outcome": "passed" if (unchanged and labelled) else "failed",
    }


def evaluate_g4(
    *, score_columns: list[str], claim_columns: list[str], subtype_call_rows_changed: int
) -> dict:
    """Evaluate the *pre-registered* G4 metric: ``subtype_call_rows_changed == 0``.

    The pre-registered threshold is the row count.  The column scan is a
    supporting observation and distinguishes three cases explicitly, because a
    naive substring test on "family_call" also matches the boundary columns that
    *declare the absence* of a call (``family_call_made=false``):

    * ``call_columns``                 - columns that would carry a call (fail)
    * ``boundary_declaration_columns`` - columns asserting no call is made
    * ``literal_token_columns``        - every column whose name contains the
      tokens, reported verbatim so the literal reading is never hidden
    """
    all_columns = [*score_columns, *claim_columns]
    literal = sorted(
        column for column in all_columns
        if any(token in column for token in FORBIDDEN_OUTPUT_COLUMN_TOKENS)
    )
    boundary = sorted(
        column for column in literal
        if column.endswith("_made") or column.endswith("_impact")
    )
    calling = sorted(set(literal) - set(boundary))
    return {
        "gate_id": "G4",
        "gate_name": "no_family_call",
        "threshold": "subtype_call_rows_changed == 0",
        "registered_metric": "subtype_call_rows_changed",
        "call_columns": calling,
        "boundary_declaration_columns": boundary,
        "literal_token_columns": literal,
        "literal_token_scan_note": (
            "boundary_declaration_columns contain the tokens 'family_call'/'subtype_call' "
            "because they assert the absence of a call; they are listed verbatim and are "
            "not treated as call columns"
        ),
        "subtype_call_rows_changed": subtype_call_rows_changed,
        "registry_files_written": [],
        "family_call_made": False,
        "new_family_call_made": False,
        "outcome": "passed" if (not calling and subtype_call_rows_changed == 0) else "failed",
    }


# ------------------------------------------------------------------- claims
def build_claim_rows(
    *,
    hit_rows: list[dict],
    motif_rows: dict[str, dict],
    tag_layer_accessions: set[str],
    tag_layer_rule: str = (
        "sbd_pf06850_binding_state == detected AND lipase_box_state == not_detected_pattern"
    ),
) -> list[dict]:
    claims: list[dict] = []
    for hit in hit_rows:
        accession = hit["accession"]
        motif = motif_rows.get(accession, {})
        status = motif.get("candidate_assignment_status", "")
        prior = motif.get("candidate_prior_superfamily", "")
        unassigned = status == "unassigned_PhaDED_like" or prior == "" or not status
        if not unassigned:
            continue
        claims.append({
            "accession": accession,
            "best_E": hit.get("best_E", ""),
            "best_score": hit.get("best_score", ""),
            "candidate_prior_superfamily": prior,
            "candidate_assignment_status": status,
            "unassigned_or_no_prior": "true",
            "task6_tag_layer_claimed": "true" if accession in tag_layer_accessions else "false",
            "task6_tag_layer_rule": tag_layer_rule,
            "claim_layer": MODEL_LAYER,
            "model_layer": MODEL_LAYER,
            "family_call_made": "false",
            "new_family_call_made": "false",
            "subtype_call_impact": "none",
            "boundary_text": BOUNDARY_TEXT,
        })
    return claims


def evaluate_confounder_burden(*, records: list[str], hit_sources: list[dict[str, dict]]) -> dict:
    """Broad-specificity check against an explicitly declared record universe.

    ``records`` must be the FASTA headers of the scored set: counting keys of a
    parsed tblout would over-count, because every ``|``-delimited name is
    additionally indexed under its component aliases.
    """
    records = [name for name in records if name]
    hit_list: list[dict] = []
    for record in records:
        for source in hit_sources:
            entry = source.get(record)
            if entry and entry.get("best_e") is not None and entry["best_e"] < DISCOVERY_EVALUE_THRESHOLD:
                hit_list.append({
                    "accession": record,
                    "best_E": format_evalue(entry["best_e"]),
                    "best_score": format_score(entry["best_score"]),
                })
                break
    return {
        "records": len(records),
        "hits": len(hit_list),
        "hit_rate": (len(hit_list) / len(records)) if records else 0.0,
        "hit_list": hit_list,
        "measurement_role": "broad_specificity_observation_not_a_pre_registered_gate",
    }


def build_tag_layer_accessions(motif_rows: dict[str, dict]) -> set[str]:
    return {
        accession for accession, row in motif_rows.items()
        if row.get("sbd_pf06850_binding_state") == "detected"
        and row.get("lipase_box_state") == "not_detected_pattern"
    }


# ------------------------------------------------------------------ labelling
LABEL_TOKENS = ("discovery_hmm_uncalibrated", "cys_discovery_uncalibrated")
DISCOVERY_LAYER_OUTPUTS = (
    "inputs/training_set.faa",
    "inputs/training_set_manifest.tsv",
    "inputs/training_set_summary.json",
    "inputs/falsification_preregistration.tsv",
    "inputs/falsification_preregistration.json",
    "inputs/falsification_preregistration.md",
    "results/cys_discovery.hmm",
    "results/cys_discovery_training.aln",
    "results/cys_discovery_scores.tsv",
    "results/cys_discovery_claims.tsv",
    "results/gate_evaluation.json",
    "results/discovery_manifest.json",
)


def scan_labelled_outputs(paths) -> tuple[list[str], list[str]]:
    """Split paths into (carrying the discovery-layer label, not carrying it)."""
    labelled: list[str] = []
    unlabelled: list[str] = []
    for path in paths:
        candidate = Path(path)
        if not candidate.is_file():
            unlabelled.append(str(path))
            continue
        text = candidate.read_bytes().decode("utf-8", errors="ignore")
        (labelled if any(token in text for token in LABEL_TOKENS) else unlabelled).append(str(path))
    return labelled, unlabelled


# --------------------------------------------------------------------- driving
def _load_motif_rows(path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["accession"]] = row
    return out


def run(args) -> dict:
    run_dir = Path(args.run_dir)
    results = run_dir / "results"

    primary = parse_tblout(args.candidate_tblout)
    supplementary = parse_tblout(args.candidate_tblout_max) if args.candidate_tblout_max else {}
    reference = parse_tblout(args.reference_tblout)
    reference_max = parse_tblout(args.reference_tblout_max) if args.reference_tblout_max else {}
    confounder = parse_tblout(args.confounder_tblout) if args.confounder_tblout else {}
    confounder_max = parse_tblout(args.confounder_tblout_max) if args.confounder_tblout_max else {}

    candidates = read_candidate_accessions(args.candidate_faa)
    motif_rows = _load_motif_rows(args.motif_tsv)
    ledger_rows = read_tsv(args.ledger)

    score_rows = build_score_rows(candidates, primary, supplementary)
    write_tsv(results / "cys_discovery_scores.tsv", SCORE_FIELDS, score_rows)

    hit_rows = [row for row in score_rows if row["discovery_hit"] == "true"]
    tag_layer = build_tag_layer_accessions(motif_rows)
    claim_rows = build_claim_rows(hit_rows=hit_rows, motif_rows=motif_rows,
                                  tag_layer_accessions=tag_layer)
    write_tsv(results / "cys_discovery_claims.tsv", CLAIM_FIELDS, claim_rows)

    # G1 also uses the MCC control, which is scored inside the reference panel run
    control_hits = {a: e for a, e in reference.items() if a == MCC_CONTROL_ACCESSION}
    control_hits.update({a: e for a, e in reference_max.items() if a == MCC_CONTROL_ACCESSION})

    training_manifest = read_tsv(args.training_manifest)
    training_accessions = [
        row["accession"] for row in training_manifest if row["include_in_training"] == "true"
    ]

    g1 = evaluate_g1(reference_hits=reference, ledger_rows=ledger_rows,
                     detected_controls=control_hits)
    g2 = evaluate_g2(reference_hits={**reference_max, **reference},
                     training_accessions=training_accessions,
                     threshold=float(args.sensitivity_threshold))
    # the outputs that already exist at G3 evaluation time are scanned for the
    # discovery-layer label; gate_evaluation.json and discovery_manifest.json are
    # written after this check and carry "model_layer" by construction
    label_paths = [run_dir / name for name in DISCOVERY_LAYER_OUTPUTS
                   if (run_dir / name).is_file()]
    labelled, unlabelled = scan_labelled_outputs(label_paths)
    g3 = evaluate_g3(
        before_sha256=args.formal_scan_models_before_sha256,
        after_sha256=sha256_file(args.formal_scan_models),
        labelled_outputs=labelled,
        unlabelled_outputs=unlabelled,
    )
    g4 = evaluate_g4(score_columns=SCORE_FIELDS, claim_columns=CLAIM_FIELDS,
                     subtype_call_rows_changed=int(args.subtype_call_rows_changed))

    gates = [g1, g2, g3, g4]
    falsified = [gate["gate_id"] for gate in gates if gate["outcome"] == "failed"]
    gate_payload = {
        "schema": "phaded-cys-discovery-gate-evaluation-v1",
        "run_id": run_dir.name,
        "model_layer": MODEL_LAYER,
        "discovery_evalue_threshold": "1e-05",
        "hmmsearch_database_size_Z": args.database_size_Z,
        "threshold_lowered_after_failure": False,
        "gates": gates,
        "falsified_gate_ids": falsified,
        "status": "falsified" if falsified else "all_gates_passed",
        "calibration_gate_unchanged": True,
        "calibration_gate_met": False,
        "boundary": BOUNDARY,
        "boundary_text": BOUNDARY_TEXT,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }
    with (results / "gate_evaluation.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(gate_payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    confounder_burden = evaluate_confounder_burden(
        records=read_candidate_accessions(args.confounder_faa) if args.confounder_faa else [],
        hit_sources=[confounder, confounder_max],
    )
    reference_records = read_candidate_accessions(args.reference_faa) if args.reference_faa else []
    confounder_hits = confounder_burden["hit_list"]
    hit_accessions = {row["accession"] for row in hit_rows}
    claims_overlap = sum(1 for row in claim_rows if row["task6_tag_layer_claimed"] == "true")

    manifest = {
        "schema": "phaded-cys-discovery-manifest-v1",
        "run_id": run_dir.name,
        "model_layer": MODEL_LAYER,
        "boundary": BOUNDARY,
        "boundary_text": BOUNDARY_TEXT,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "inputs": args.input_records,
        "outputs": {str(p): sha256_file(p) for p in [
            results / "cys_discovery_scores.tsv",
            results / "cys_discovery_claims.tsv",
            results / "gate_evaluation.json",
        ]},
        "tools": args.tool_versions,
        "commands": args.command,
        "threads": {"hmmsearch_cpu": 40, "hmmbuild_cpu": 40, "mafft_thread": 40,
                    "max_single_task_threads_allowed": 40},
        "server": args.server_snapshot,
        "counts": {
            "candidates_scored": len(score_rows),
            "candidate_discovery_hits": len(hit_rows),
            "claims_unassigned_or_no_prior": len(claim_rows),
            "claims_overlapping_task6_tag_layer": claims_overlap,
            "task6_tag_layer_size": len(tag_layer),
            "hits_overlapping_task6_tag_layer": len(hit_accessions & tag_layer),
            "reference_records_scored": len(reference_records) or len(reference),
            "confounder_hits": confounder_burden["hits"],
            "confounder_set_scored": confounder_burden["records"],
            "confounder_hit_rate": confounder_burden["hit_rate"],
            "training_sequences": len(training_accessions),
        },
        "confounder_burden": {
            "records": confounder_burden["records"],
            "hits": confounder_burden["hits"],
            "hit_rate": confounder_burden["hit_rate"],
            "measurement_role": confounder_burden["measurement_role"],
            "note": (
                "the x1-hydrophobicity confounder set is a lipase/esterase-suspicion layer, not a "
                "validated negative set; a high hit rate here is a breadth caveat, not a gate failure"
            ),
        },
        "confounder_hit_list": confounder_hits[:200],
        "confounder_hit_list_truncated": len(confounder_hits) > 200,
        "family_call_made": False,
        "new_family_call_made": False,
        "subtype_call_rows_changed": 0,
        "registry_modified": False,
        "training_equals_calibration": False,
    }
    with (results / "discovery_manifest.json").open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return gate_payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--candidate-tblout", required=True)
    parser.add_argument("--candidate-tblout-max", default=None)
    parser.add_argument("--reference-tblout", required=True)
    parser.add_argument("--reference-tblout-max", default=None)
    parser.add_argument("--confounder-tblout", default=None)
    parser.add_argument("--confounder-tblout-max", default=None)
    parser.add_argument("--confounder-faa", default=None,
                        help="the scored confounder FASTA; defines the record universe for the "
                             "broad-specificity count (avoids alias over-counting)")
    parser.add_argument("--reference-faa", default=None,
                        help="the scored reference FASTA; defines the reference record universe")
    parser.add_argument("--candidate-faa", required=True)
    parser.add_argument("--motif-tsv", required=True)
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--training-manifest", required=True)
    parser.add_argument("--formal-scan-models", required=True)
    parser.add_argument("--formal-scan-models-before-sha256", required=True)
    parser.add_argument("--sensitivity-threshold", default="0.90")
    parser.add_argument("--subtype-call-rows-changed", default="0")
    parser.add_argument("--database-size-Z", default="109087")
    parser.add_argument("--input-record", action="append", default=[],
                        help="name=path, hashed into the manifest")
    parser.add_argument("--tool", action="append", default=[], help="name=version string")
    parser.add_argument("--command", action="append", default=[])
    parser.add_argument("--server-snapshot", default=None,
                        help="path to a JSON file with the server load snapshot")
    args = parser.parse_args(argv)

    args.input_records = {}
    for item in args.__dict__.get("input_record", []):
        name, _, path = item.partition("=")
        args.input_records[name] = (
            {"path": path, "sha256": sha256_file(path)} if Path(path).is_file()
            else {"path": path, "sha256": "pending", "status": "missing"}
        )
    args.tool_versions = {}
    for item in args.__dict__.get("tool", []):
        name, _, version = item.partition("=")
        args.tool_versions[name] = version
    args.server_snapshot = (
        json.loads(Path(args.server_snapshot).read_text(encoding="utf-8"))
        if args.server_snapshot and Path(args.server_snapshot).is_file() else "pending"
    )
    payload = run(args)
    print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    return 3 if payload["falsified_gate_ids"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
