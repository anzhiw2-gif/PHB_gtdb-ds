#!/usr/bin/env python3
"""Audit PhaDED catalytic-pattern and architecture evidence conservatively.

The output deliberately distinguishes sequence-pattern support from a validated
catalytic mechanism.  Missing reference residue mappings remain pending and
are never converted into biological negatives.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = set("LIVMFWAY")
LIPASE_BOX = re.compile(r"G([A-Z])S([A-Z])G")
ASP_PATTERN = re.compile(r"G[A-Z]{2}DYTV")
HIS_PATTERN = re.compile(r"GM[A-Z]H[A-Z]{2}P[A-Z]{2}G")
OXY_PATTERN = re.compile(r"H[A-Z]?G[A-Z]?C[A-Z]?Q")
AHSMG_PATTERN = re.compile(r"AHSMG")
PROFILE_FAMILY_RE = re.compile(r"DED_hfam_(\d+)")
ALIGNMENT_ROW_RE = re.compile(r"^\s*([^\s]+)\s+([A-Za-z*-]+)(?:\s+\d+)?\s*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
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
                records[current] = "".join(chunks).upper()
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
        records[current] = "".join(chunks).upper()
    if not records:
        raise ValueError(f"FASTA is empty: {path}")
    return records


def _coord(match: re.Match[str] | None) -> str:
    return f"{match.start() + 1}-{match.end()}" if match else ""


def _first(pattern: re.Pattern[str], sequence: str) -> re.Match[str] | None:
    return pattern.search(sequence.rstrip("*"))


def _pfam_coordinates(pfam_coordinates: str, accession: str) -> str:
    for token in (pfam_coordinates or "").split(";"):
        if token.startswith(accession + ":"):
            return token.split(":", 1)[1]
    return ""


def inspect_alignment(path: Path) -> dict[str, str]:
    """Classify an existing CLUSTAL-like file without treating a header as an alignment."""
    text = path.read_text(encoding="ascii", errors="replace")
    rows = []
    for line in text.splitlines():
        match = ALIGNMENT_ROW_RE.match(line)
        if match and match.group(1).upper() not in {"CLUSTAL", "MUSCLE", "PROBCONS"}:
            rows.append(match.group(1))
    unique_rows = sorted(set(rows))
    return {
        "alignment_status": "usable_reference_alignment" if unique_rows else "empty_or_header_only",
        "alignment_sequence_rows": str(len(unique_rows)),
    }


def _expected_residue(subtype: str) -> str:
    if "nPHASCL without lipase box" in subtype:
        return "Cys"
    return "Ser"


def _architecture_state(subtype: str, pfam_accessions: str, pfam_coordinates: str) -> dict[str, str]:
    accessions = {item.strip() for item in (pfam_accessions or "").split(";") if item.strip()}
    tested = bool(accessions)
    has_sbd = "PF06850" in accessions
    if "extracellular dPHAMCL" in subtype:
        sbd = "not_expected_for_subtype"
        linker = "not_expected_for_subtype"
        lid = "pending_reference_annotation"
    else:
        sbd = "supported" if has_sbd else ("not_detected_in_tested_pfam" if tested else "pending_tool_input")
        linker = "pending_reference_annotation"
        lid = "pending_reference_annotation" if "nPHAMCL" in subtype else "not_expected_for_subtype"
    return {
        "sbd_state": sbd,
        "sbd_coordinates": _pfam_coordinates(pfam_coordinates, "PF06850") if has_sbd else "",
        "linker_state": linker,
        "linker_coordinates": "",
        "lid_state": lid,
        "lid_coordinates": "",
    }


def audit_sequence(
    accession: str,
    sequence: str,
    subtype: str,
    *,
    pfam_accessions: str = "",
    pfam_coordinates: str = "",
) -> dict[str, str]:
    """Return accession-bound motif evidence without phenotype promotion."""
    if not accession.strip() or not sequence.strip():
        raise ValueError("accession and sequence are required")
    sequence = sequence.strip().upper()
    clean = sequence.rstrip("*")
    lipase = _first(LIPASE_BOX, clean)
    asp = _first(ASP_PATTERN, clean)
    his = _first(HIS_PATTERN, clean)
    oxy = _first(OXY_PATTERN, clean)
    ahsmg = _first(AHSMG_PATTERN, clean)
    expected_residue = _expected_residue(subtype)
    no_lipase_family = "nPHASCL without lipase box" in subtype
    phaz7 = "native-SCL/PhaZ7" in subtype
    if no_lipase_family:
        lipase_state = "not_expected_for_subtype"
        lipase_coord = ""
        catalytic_state = "pending_reference_annotation"
        catalytic_coord = ""
    elif phaz7:
        lipase_state = "not_expected_for_subtype"
        lipase_coord = ""
        catalytic_state = "supported" if ahsmg else "pending_reference_annotation"
        catalytic_coord = _coord(ahsmg)
    else:
        lipase_state = "supported" if lipase else "not_detected_pattern"
        lipase_coord = _coord(lipase)
        catalytic_state = "supported" if lipase else "not_detected_pattern"
        catalytic_coord = f"{lipase.start() + 3}" if lipase else ""

    def catalytic_pattern_state(match: re.Match[str] | None) -> str:
        return "supported" if match else "not_detected_pattern"

    asp_state = catalytic_pattern_state(asp)
    his_state = catalytic_pattern_state(his)
    oxy_state = catalytic_pattern_state(oxy)
    if lipase and oxy and "extracellular dPHASCL type 1" in subtype and oxy.start() > lipase.start():
        oxy_state = "conflict_relative_position"
    if lipase and oxy and "extracellular dPHASCL type 2" in subtype and oxy.start() < lipase.start():
        oxy_state = "conflict_relative_position"

    architecture = _architecture_state(subtype, pfam_accessions, pfam_coordinates)
    pending: list[str] = []
    if no_lipase_family:
        pending.append("cys_his_asp_reference_mapping_pending")
    if architecture["linker_state"] == "pending_reference_annotation":
        pending.append("linker_profile_not_bound")
    if architecture["lid_state"] == "pending_reference_annotation":
        pending.append("lid_profile_not_bound")
    if "extracellular dPHASCL" in subtype and architecture["sbd_state"] != "supported":
        pending.append("sbd_reference_or_pfam_binding_pending")
    catalytic_states = [lipase_state, catalytic_state, his_state, asp_state, oxy_state]
    supported_count = sum(value == "supported" for value in catalytic_states)
    if no_lipase_family:
        panel_status = "pending_reference_annotation"
    elif supported_count == len(catalytic_states) and not pending:
        panel_status = "complete_pattern_panel"
    elif supported_count >= 2:
        panel_status = "partial_catalytic_pattern_panel"
    else:
        panel_status = "partial_lipase_box_only" if lipase else "not_tested_full_library"
    if oxy_state == "conflict_relative_position":
        panel_status = "conflict_relative_position"
        pending.append("oxyanion_relative_position_conflict")
    return {
        "accession": accession,
        "motif_reference_subtype": subtype,
        "motif_expected_catalytic_residue": expected_residue,
        "lipase_box_state": lipase_state,
        "lipase_box_coordinates": lipase_coord,
        "lipase_box_x1": lipase.group(1) if lipase else "",
        "catalytic_ser_cys_state": catalytic_state,
        "catalytic_ser_cys_coordinates": catalytic_coord,
        "his_state": his_state,
        "his_coordinates": _coord(his),
        "asp_state": asp_state,
        "asp_coordinates": _coord(asp),
        "oxyanion_hole_state": oxy_state,
        "oxyanion_hole_coordinates": _coord(oxy),
        "ahsmg_state": "supported" if ahsmg else "not_detected_pattern",
        "ahsmg_coordinates": _coord(ahsmg),
        "sbd_state": architecture["sbd_state"],
        "sbd_coordinates": architecture["sbd_coordinates"],
        "linker_state": architecture["linker_state"],
        "linker_coordinates": architecture["linker_coordinates"],
        "lid_state": architecture["lid_state"],
        "lid_coordinates": architecture["lid_coordinates"],
        "motif_panel_status": panel_status,
        "motif_completeness": "pattern_complete" if panel_status == "complete_pattern_panel" else "pattern_partial",
        "motif_evidence_level": "motif_pattern_only",
        "motif_pending_reason": ";".join(sorted(set(pending))) or "reference_pattern_and_domain_binding_available",
        "motif_phenotype_boundary": "sequence-pattern/domain evidence only; not catalytic validation and not PHB/PHA phenotype proof",
    }


def _index(rows: Iterable[Mapping[str, str]], key: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        value = str(row.get(key, "")).strip()
        if value:
            if value in indexed:
                raise ValueError(f"duplicate {key}: {value}")
            indexed[value] = dict(row)
    return indexed


def audit_reference_panel(
    ledger_rows: list[Mapping[str, str]],
    reference_sequences: Mapping[str, str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in ledger_rows:
        accession = str(source.get("accession", "")).strip()
        if not accession or accession not in reference_sequences:
            continue
        result = audit_sequence(
            accession,
            reference_sequences[accession],
            str(source.get("phaded_superfamily", "")),
        )
        result.update({
            "reference_id": source.get("reference_id", ""),
            "phaded_family_id": source.get("phaded_family_id", ""),
            "evidence_status": source.get("evidence_status", ""),
            "primary_doi": source.get("primary_doi", ""),
        })
        rows.append(result)
    return rows


def build_alignment_manifest(
    ledger_rows: list[Mapping[str, str]], alignment_dir: Path | None,
) -> list[dict[str, str]]:
    families: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in ledger_rows:
        subtype = str(row.get("phaded_superfamily", "")).strip()
        family = str(row.get("phaded_family_id", "")).strip()
        if subtype and family:
            families[(subtype, family)].add(str(row.get("reference_id", "")))
    output: list[dict[str, str]] = []
    for (subtype, family), refs in sorted(families.items()):
        match = PROFILE_FAMILY_RE.search(family)
        alignment = alignment_dir / f"aln{match.group(1)}.aln" if alignment_dir and match else None
        if alignment and alignment.exists():
            quality = inspect_alignment(alignment)
            status = quality["alignment_status"]
            path = str(alignment.resolve())
            size = str(alignment.stat().st_size)
            digest = sha256(alignment)
            sequence_rows = quality["alignment_sequence_rows"]
        else:
            status = "pending_alignment_file"
            path = str(alignment.resolve()) if alignment else "pending"
            size = "pending"
            digest = "pending"
            sequence_rows = "pending"
        output.append({
            "phaded_superfamily": subtype,
            "phaded_family_id": family,
            "reference_count": str(len(refs)),
            "reference_alignment": path,
            "alignment_status": status,
            "alignment_sequence_rows": sequence_rows,
            "alignment_size": size,
            "alignment_sha256": digest,
        })
    return output


def _write_tsv(path: Path, rows: list[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    if not fields:
        fields = ["status"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_audit(
    assignments_path: Path,
    features_path: Path,
    pfam_path: Path,
    candidate_fasta: Path,
    reference_ledger: Path,
    reference_fasta: Path,
    output_dir: Path,
    run_id: str,
    alignment_dir: Path | None = None,
) -> dict[str, object]:
    assignments = read_tsv(assignments_path)
    features = _index(read_tsv(features_path), "accession")
    pfam = _index(read_tsv(pfam_path), "accession")
    sequences = read_fasta(candidate_fasta)
    if set(sequences) != {str(row.get("accession", "")).strip() for row in assignments}:
        raise ValueError("candidate FASTA and assignment accession sets differ")
    candidate_rows: list[dict[str, str]] = []
    for assignment in assignments:
        accession = str(assignment.get("accession", "")).strip()
        subtype = str(assignment.get("phaded_superfamily_best", "")).strip()
        pfam_row = pfam.get(accession, {})
        result = audit_sequence(
            accession,
            sequences[accession],
            subtype,
            pfam_accessions=pfam_row.get("pfam_accessions", ""),
            pfam_coordinates=pfam_row.get("pfam_coordinates", ""),
        )
        feature = features.get(accession, {})
        result.update({
            "genome": assignment.get("genome", ""),
            "assignment_status": assignment.get("assignment_status", ""),
            "signalp_evidence": feature.get("signalp_evidence", "pending"),
            "signalp_class": feature.get("signalp_class", "pending"),
            "pfam_accessions": pfam_row.get("pfam_accessions", ""),
            "pfam_coordinates": pfam_row.get("pfam_coordinates", ""),
            "architecture_consistency": pfam_row.get("architecture_consistency", "pending"),
        })
        candidate_rows.append(result)
    reference_sequences_raw = read_fasta(reference_fasta)
    reference_sequences = {
        key.split("|", 1)[1] if "|" in key else key: sequence
        for key, sequence in reference_sequences_raw.items()
    }
    reference_rows = audit_reference_panel(read_tsv(reference_ledger), reference_sequences)
    alignments = build_alignment_manifest(read_tsv(reference_ledger), alignment_dir)
    counts = {
        field: dict(sorted(Counter(row.get(field, "") for row in candidate_rows).items()))
        for field in ("motif_panel_status", "motif_evidence_level", "sbd_state", "linker_state", "lid_state", "oxyanion_hole_state")
    }
    subtype_summary: list[dict[str, object]] = []
    by_subtype: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in reference_rows:
        by_subtype[row.get("motif_reference_subtype", "")].append(row)
    for subtype, rows in sorted(by_subtype.items()):
        subtype_summary.append({
            "phaded_superfamily": subtype,
            "reference_count": len(rows),
            "lipase_box_supported": sum(row.get("lipase_box_state") == "supported" for row in rows),
            "his_pattern_supported": sum(row.get("his_state") == "supported" for row in rows),
            "asp_pattern_supported": sum(row.get("asp_state") == "supported" for row in rows),
            "oxyanion_pattern_supported": sum(row.get("oxyanion_hole_state") == "supported" for row in rows),
            "sbd_domain_supported": sum(row.get("sbd_state") == "supported" for row in rows),
            "reference_annotation_status": "sequence_pattern_only; catalytic residues and domain roles require accession-bound curation",
        })
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "motif_candidate_evidence.tsv", candidate_rows)
    _write_tsv(output_dir / "reference_motif_panel.tsv", reference_rows)
    _write_tsv(output_dir / "reference_motif_summary.tsv", subtype_summary)
    _write_tsv(output_dir / "reference_alignment_manifest.tsv", alignments)
    manifest = {
        "schema_version": "1.0",
        "run_id": run_id,
        "status": "completed_candidate_only",
        "candidate_count": len(candidate_rows),
        "reference_count": len(reference_rows),
        "counts": counts,
        "reference_summary": subtype_summary,
        "inputs": {
            name: {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}
            for name, path in {
                "assignments": assignments_path, "features": features_path, "pfam": pfam_path,
                "candidate_fasta": candidate_fasta, "reference_ledger": reference_ledger,
                "reference_fasta": reference_fasta,
            }.items()
        },
        "outputs": {
            name: {"path": str((output_dir / name).resolve()), "size": (output_dir / name).stat().st_size, "sha256": sha256(output_dir / name)}
            for name in ("motif_candidate_evidence.tsv", "reference_motif_panel.tsv", "reference_motif_summary.tsv", "reference_alignment_manifest.tsv")
        },
        "phenotype_boundary": "Motif patterns, Pfam architecture, and reference alignments are candidate evidence only; they do not prove PHB/PHA degradation.",
    }
    (output_dir / "motif_audit_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--pfam", type=Path, required=True)
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--reference-ledger", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--alignment-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args(argv)
    print(json.dumps(run_audit(
        args.assignments, args.features, args.pfam, args.candidate_fasta,
        args.reference_ledger, args.reference_fasta, args.output_dir, args.run_id,
        args.alignment_dir,
    ), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
