#!/usr/bin/env python3
"""Build the Cys discovery-layer training set (2026-09-17, Step 3).

The input is the *frozen* 723-row PhaDED reference ledger restricted to the
``intracellular nPHASCL without lipase box`` superfamily (276 records).  An
explicit integrity filter is applied and every record is written to a manifest
with its keep/drop decision and the reason.  Nothing frozen is modified.

Integrity criteria (all constants are recorded in the summary):

``localization``
    ``reported_localization`` must be ``intracellular``.  Records that are
    annotated as secreted/extracellular are dropped.
``signalp``
    The frozen reference layer contains **no** SignalP run (0 reference
    accessions appear in the candidate-layer SignalP table), so this evidence
    is recorded as ``pending`` and the sequence is retained with an annotation,
    as required by the run contract.
``lipase_box``
    The superfamily is *defined* as having no lipase box.  Any member whose
    frozen mapping-manifest row reports ``lipase_box_state=supported`` for a
    superfamily that does not expect one is an annotation conflict with the
    family definition (the lipase box is the hallmark of the extracellular
    dPHASCL architecture).  The flag is reproduced independently here by
    scanning the sequence for the canonical Knoll 2009 pattern ``Gx1Sx2G``;
    both signals must agree for a record to be dropped.
``length``
    Project precedent (``pipeline/scripts/split_ephaz_seeds.py``) sets
    ``minimum_length_for_core = 200`` with the policy
    ``retain_in_broad_and_review``.  Records below the minimum are retained in
    the run inputs and in the manifest but are excluded from the core training
    alignment.

Boundary: this is a *discovery layer* (``discovery_hmm_uncalibrated``) input.
Training on DED family members - including ``annotation_only`` records -
reproduces what the PHA Depolymerase Engineering Database itself does (Knoll
2009 provides an HMMER profile for every family for in silico identification).
It is not calibration, it grants no registry eligibility, and it makes no
family call.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
from pathlib import Path

CYS_SUPERFAMILY = "intracellular nPHASCL without lipase box"
MODEL_LAYER = "discovery_hmm_uncalibrated"
MIN_CORE_LENGTH = 200
SIGNALP_PENDING = "pending_reference_layer_has_no_signalp_run"
LIPASE_BOX_PATTERN = re.compile(r"G.S.G")
BOUNDARY = "candidate_only"

EXCLUSION_LOCALIZATION = "reported_localization_not_intracellular"
EXCLUSION_LIPASE_BOX = "lipase_box_annotation_conflict_with_superfamily"
EXCLUSION_SHORT = "short_sequence_below_core_minimum"

MANIFEST_FIELDS = [
    "reference_id",
    "accession",
    "phaded_family_id",
    "evidence_status",
    "reported_localization",
    "signalp_evidence",
    "length",
    "lipase_box_state",
    "lipase_box_expectation_from_superfamily",
    "canonical_lipase_box_positions",
    "decision",
    "exclusion_reason",
    "include_in_training",
    "sequence_sha256",
    "sequence_sha256_is_of_protein_sequence",
    "model_layer",
    "training_role",
]


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


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_reference_sequences(path: os.PathLike[str] | str) -> dict[str, str]:
    """Return {reference_id: sequence} for ``>reference_id|accession`` FASTAs."""
    sequences: dict[str, str] = {}
    header: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8", newline="") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    sequences[header.split("|", 1)[0]] = "".join(parts)
                header, parts = line[1:], []
            elif line.strip():
                parts.append(line.strip())
        if header is not None:
            sequences[header.split("|", 1)[0]] = "".join(parts)
    return sequences


def load_reference_accessions(path: os.PathLike[str] | str) -> dict[str, str]:
    accessions: dict[str, str] = {}
    with open(path, encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                reference_id, accession = (header.split("|", 1) + [""])[:2]
                accessions[reference_id] = accession
    return accessions


def detect_canonical_lipase_box(sequence: str) -> list[int]:
    """1-based positions of the glycine that starts a canonical ``Gx1Sx2G``."""
    return [match.start() + 1 for match in LIPASE_BOX_PATTERN.finditer(sequence)]


def build_training_set(
    *,
    ledger_rows: list[dict],
    mapping_rows: list[dict],
    reference_by_id: dict[str, str],
    min_core_length: int = MIN_CORE_LENGTH,
) -> tuple[list[dict], list[dict]]:
    mapping_by_id = {row.get("reference_id", ""): row for row in mapping_rows}
    rows: list[dict] = []
    fasta_rows: list[dict] = []

    for ledger_row in ledger_rows:
        if ledger_row.get("phaded_superfamily") != CYS_SUPERFAMILY:
            continue
        reference_id = ledger_row.get("reference_id", "")
        accession = ledger_row.get("accession", "")
        sequence = reference_by_id.get(reference_id, "")
        mapping = mapping_by_id.get(reference_id, {})
        positions = detect_canonical_lipase_box(sequence)
        lipase_state = mapping.get("lipase_box_state", "")
        expectation = mapping.get("lipase_box_expectation_from_superfamily", "")

        exclusion = ""
        if ledger_row.get("reported_localization") != "intracellular":
            exclusion = EXCLUSION_LOCALIZATION
        elif lipase_state == "supported" and expectation != "expected_for_superfamily" and positions:
            exclusion = EXCLUSION_LIPASE_BOX
        elif len(sequence) < min_core_length:
            exclusion = EXCLUSION_SHORT

        include = "false" if exclusion else "true"
        rows.append({
            "reference_id": reference_id,
            "accession": accession,
            "phaded_family_id": ledger_row.get("phaded_family_id", ""),
            "evidence_status": ledger_row.get("evidence_status", ""),
            "reported_localization": ledger_row.get("reported_localization", ""),
            "signalp_evidence": SIGNALP_PENDING,
            "length": len(sequence),
            "lipase_box_state": lipase_state,
            "lipase_box_expectation_from_superfamily": expectation,
            "canonical_lipase_box_positions": ";".join(str(p) for p in positions),
            "decision": "training_core" if include == "true" else "excluded",
            "exclusion_reason": exclusion,
            "include_in_training": include,
            "sequence_sha256": sha256_text(sequence),
            "sequence_sha256_is_of_protein_sequence": "true",
            "model_layer": MODEL_LAYER,
            "training_role": "unsupervised_family_alignment_member",
        })
        if include == "true":
            fasta_rows.append({
                "reference_id": reference_id,
                "accession": accession,
                "header": f"{reference_id}|{accession}|{MODEL_LAYER}",
                "sequence": sequence,
            })
    return rows, fasta_rows


def write_training_set(
    out_faa: os.PathLike[str] | str,
    out_manifest: os.PathLike[str] | str,
    out_summary: os.PathLike[str] | str,
    *,
    ledger_rows: list[dict],
    mapping_rows: list[dict],
    reference_by_id: dict[str, str],
    min_core_length: int = MIN_CORE_LENGTH,
    source_paths: dict[str, str] | None = None,
) -> dict:
    rows, fasta_rows = build_training_set(
        ledger_rows=ledger_rows,
        mapping_rows=mapping_rows,
        reference_by_id=reference_by_id,
        min_core_length=min_core_length,
    )
    faa_path = Path(out_faa)
    faa_path.parent.mkdir(parents=True, exist_ok=True)
    with faa_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in fasta_rows:
            handle.write(f">{row['header']}\n{row['sequence']}\n")
    write_tsv(out_manifest, MANIFEST_FIELDS, rows)

    reasons: dict[str, int] = {}
    for row in rows:
        key = row["exclusion_reason"] or "retained"
        reasons[key] = reasons.get(key, 0) + 1
    summary = {
        "schema": "phaded-cys-discovery-training-set-v1",
        "model_layer": MODEL_LAYER,
        "superfamily": CYS_SUPERFAMILY,
        "candidate_superfamily_records": len(rows),
        "training_sequence_count": len(fasta_rows),
        "excluded_count": len(rows) - len(fasta_rows),
        "exclusion_reason_counts": dict(sorted(reasons.items())),
        "criteria": {
            "localization_required": "intracellular",
            "signalp_evidence": SIGNALP_PENDING,
            "lipase_box_conflict_rule": (
                "frozen lipase_box_state=supported in a superfamily that does not expect a "
                "lipase box AND independently reproduced canonical Gx1Sx2G"
            ),
            "canonical_lipase_box_pattern": "Gx1Sx2G",
            "min_core_length": min_core_length,
            "min_core_length_precedent": (
                "pipeline/scripts/split_ephaz_seeds.py criteria.minimum_length_for_core"
            ),
            "short_sequence_policy": "retain_in_manifest_and_broad_not_in_core_training",
            "training_source": "frozen DED family members (Knoll 2009 style profile construction)",
        },
        "frozen_ledger_modified": False,
        "family_call_made": False,
        "gate_values_changed": False,
        "training_equals_calibration": False,
        "registry_eligible": False,
        "boundary": BOUNDARY,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "outputs": {str(out_faa): sha256_text(faa_path.read_text(encoding="utf-8"))},
    }
    if source_paths:
        summary["source_inputs"] = source_paths
    target = Path(out_summary)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--reference-faa", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--min-core-length", type=int, default=MIN_CORE_LENGTH)
    args = parser.parse_args(argv)

    run_dir = Path(args.run_dir)
    summary = write_training_set(
        run_dir / "inputs" / "training_set.faa",
        run_dir / "inputs" / "training_set_manifest.tsv",
        run_dir / "inputs" / "training_set_summary.json",
        ledger_rows=read_tsv(args.ledger),
        mapping_rows=read_tsv(args.mapping),
        reference_by_id=load_reference_sequences(args.reference_faa),
        min_core_length=args.min_core_length,
        source_paths={
            "ledger": args.ledger,
            "reference_faa": args.reference_faa,
            "mapping": args.mapping,
        },
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
