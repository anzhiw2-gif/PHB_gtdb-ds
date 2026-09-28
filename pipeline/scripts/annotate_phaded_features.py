#!/usr/bin/env python3
"""Annotate candidate sequence/architecture evidence without phenotype calls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Mapping


OUTPUT_FIELDS = [
    "accession", "sequence_integrity", "signalp_evidence", "signalp_class",
    "lipase_box_state", "lipase_box_coordinates", "ahsmg_state", "ahsmg_coordinates",
    "catalytic_ser_cys_state", "catalytic_ser_cys_coordinates", "his_state", "his_coordinates",
    "asp_state", "asp_coordinates", "oxyanion_hole_state", "oxyanion_hole_coordinates",
    "sbd_state", "sbd_coordinates", "linker_state", "linker_coordinates", "lid_state", "lid_coordinates",
    "pfam_interpro_state", "architecture_consistency", "feature_evidence_status",
]
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _coordinates(match: re.Match[str] | None) -> str:
    if match is None:
        return ""
    return f"{match.start() + 1}-{match.end()}"


def _domain_value(domain_records: Mapping[str, object], key: str) -> tuple[str, str]:
    value = domain_records.get(key)
    if value is None:
        return "no_record", ""
    if isinstance(value, Mapping):
        state = str(value.get("state", "present"))
        coordinates = str(value.get("coordinates", ""))
        return state, coordinates
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return str(value[0]), str(value[1])
    return "present", str(value)


def annotate_sequence(
    accession: str,
    sequence: str,
    phaded_superfamily: str,
    *,
    signalp_record: Mapping[str, object] | None = None,
    domain_records: Mapping[str, object] | None = None,
) -> dict[str, str]:
    if not accession.strip():
        raise ValueError("accession is required")
    sequence = sequence.strip().upper()
    if not sequence:
        raise ValueError("sequence is required")
    sequence_integrity = "valid"
    if any(residue not in AMINO_ACIDS for residue in sequence):
        if sequence.endswith("*") and sequence.count("*") == 1 and all(residue in AMINO_ACIDS for residue in sequence[:-1]):
            sequence_integrity = "terminal_stop_only" if sequence.startswith("M") else "possible_N_truncation"
        else:
            sequence_integrity = "invalid_internal_character"
    elif not sequence.startswith("M"):
        sequence_integrity = "possible_N_truncation"
    signalp_evidence = "no_record"
    signalp_class = "pending"
    if signalp_record is not None:
        signalp_evidence = str(signalp_record.get("status", "provided"))
        signalp_class = str(signalp_record.get("type", "OTHER"))
    lipase_match = re.search(r"G.S.G", sequence)
    ahsmg_match = re.search(r"AHSMG", sequence)
    phaz7 = "PhaZ7" in phaded_superfamily or "native-SCL" in phaded_superfamily
    if phaz7:
        lipase_state, lipase_coordinates = "not_expected_for_superfamily", ""
    elif lipase_match:
        lipase_state, lipase_coordinates = "present", _coordinates(lipase_match)
    else:
        lipase_state, lipase_coordinates = "not_detected", ""
    ahsmg_state = "present" if ahsmg_match else "not_detected"
    ahsmg_coordinates = _coordinates(ahsmg_match)
    catalytic_state = "supported" if lipase_match and not phaz7 else "pending"
    catalytic_coordinates = _coordinates(lipase_match) if catalytic_state == "supported" else ""
    domain = domain_records
    if domain is None:
        sbd_state = linker_state = lid_state = pfam_state = "pending"
        sbd_coordinates = linker_coordinates = lid_coordinates = ""
    else:
        sbd_state, sbd_coordinates = _domain_value(domain, "sbd")
        linker_state, linker_coordinates = _domain_value(domain, "linker")
        lid_state, lid_coordinates = _domain_value(domain, "lid")
        pfam_state, _ = _domain_value(domain, "pfam_interpro")
    his_state = asp_state = oxyanion_state = "pending"
    his_coordinates = asp_coordinates = oxyanion_coordinates = ""
    if domain is not None:
        oxyanion_state, oxyanion_coordinates = _domain_value(domain, "oxyanion_hole")
        his_state, his_coordinates = _domain_value(domain, "his")
        asp_state, asp_coordinates = _domain_value(domain, "asp")
    if phaz7 and ahsmg_match:
        architecture = "consistent"
    elif phaz7:
        architecture = "pending"
    elif lipase_match and domain is not None:
        architecture = "partial"
    else:
        architecture = "pending"
    pending_fields = ("pending" in {pfam_state, lid_state, his_state, asp_state, oxyanion_state} or signalp_evidence == "no_record")
    evidence_status = "pending" if pending_fields else ("supported" if architecture == "consistent" else "partial")
    return {
        "accession": accession, "sequence_integrity": sequence_integrity,
        "signalp_evidence": signalp_evidence, "signalp_class": signalp_class,
        "lipase_box_state": lipase_state, "lipase_box_coordinates": lipase_coordinates,
        "ahsmg_state": ahsmg_state, "ahsmg_coordinates": ahsmg_coordinates,
        "catalytic_ser_cys_state": catalytic_state, "catalytic_ser_cys_coordinates": catalytic_coordinates,
        "his_state": his_state, "his_coordinates": his_coordinates, "asp_state": asp_state,
        "asp_coordinates": asp_coordinates, "oxyanion_hole_state": oxyanion_state,
        "oxyanion_hole_coordinates": oxyanion_coordinates, "sbd_state": sbd_state,
        "sbd_coordinates": sbd_coordinates, "linker_state": linker_state,
        "linker_coordinates": linker_coordinates, "lid_state": lid_state,
        "lid_coordinates": lid_coordinates, "pfam_interpro_state": pfam_state,
        "architecture_consistency": architecture, "feature_evidence_status": evidence_status,
    }


def annotate_rows(rows: list[dict[str, str]], sequences: Mapping[str, str]) -> list[dict[str, str]]:
    candidate_accessions = [row.get("accession", "").strip() for row in rows]
    if any(not accession for accession in candidate_accessions) or len(candidate_accessions) != len(set(candidate_accessions)):
        raise ValueError("duplicate or missing candidate accession")
    if set(candidate_accessions) != set(sequences):
        raise ValueError("FASTA accession set mismatch")
    output = []
    for row in rows:
        accession = row.get("accession", "").strip()
        if not accession or accession not in sequences:
            raise ValueError(f"candidate sequence missing: {accession}")
        declared_hash = row.get("sequence_sha256", "").strip().lower()
        if declared_hash and hashlib.sha256(sequences[accession].strip().upper().encode("ascii")).hexdigest() != declared_hash:
            raise ValueError(f"sequence SHA-256 mismatch: {accession}")
        features = annotate_sequence(
            accession, sequences[accession], row.get("phaded_superfamily", ""),
            signalp_record={"type": row["signalp_class"], "status": row.get("signalp_evidence", "provided")}
            if row.get("signalp_class") else None,
        )
        output.append(features)
    return output


def read_fasta(path: Path) -> dict[str, str]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"candidate FASTA is not a regular file: {path}")
    records: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                if current in records or not chunks:
                    raise ValueError(f"duplicate or empty FASTA record: {current}")
                records[current] = "".join(chunks)
            current = line[1:].split()[0]
            if not current:
                raise ValueError("empty FASTA header")
            chunks = []
        else:
            if current is None:
                raise ValueError("sequence appears before FASTA header")
            chunks.append(line)
    if current is not None:
        if current in records or not chunks:
            raise ValueError(f"duplicate or empty FASTA record: {current}")
        records[current] = "".join(chunks)
    if not records:
        raise ValueError("candidate FASTA is empty")
    return records


def merge_assignment_evidence(
    assignment_rows: list[Mapping[str, object]],
    feature_rows: list[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Join feature evidence without rewriting the original mapping call."""
    by_accession: dict[str, Mapping[str, object]] = {}
    for row in feature_rows:
        accession = str(row.get("accession", "")).strip()
        if not accession or accession in by_accession:
            raise ValueError(f"duplicate or missing feature accession: {accession}")
        by_accession[accession] = row
    output: list[dict[str, object]] = []
    for row in assignment_rows:
        accession = str(row.get("accession", "")).strip()
        if not accession or accession not in by_accession:
            raise ValueError(f"assignment has no matching feature row: {accession}")
        feature = by_accession[accession]
        merged = dict(row)
        merged["architecture_consistency"] = feature.get("architecture_consistency", "pending")
        merged["feature_evidence_status"] = feature.get("feature_evidence_status", "pending")
        merged["assignment_review"] = (
            "review_required" if merged["architecture_consistency"] == "conflicting" else "no_conflict_recorded"
        )
        output.append(merged)
    if len(output) != len(by_accession):
        raise ValueError("feature/assignment accession sets differ")
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--sequences", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    with args.candidates.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sequences = read_fasta(args.sequences)
    features = annotate_rows(rows, sequences)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(features)
    manifest = {
        "status": "candidate-only",
        "candidate_count": len(features),
        "inputs": {
            "candidate_table": {"path": str(args.candidates.resolve()), "sha256": sha256(args.candidates)},
            "candidate_fasta": {"path": str(args.sequences.resolve()), "sha256": sha256(args.sequences)},
        },
        "output": {"path": str(args.output.resolve()), "sha256": sha256(args.output)},
    }
    (args.output.parent / "feature_annotation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"candidate_count": len(features), "status": "candidate-only"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
