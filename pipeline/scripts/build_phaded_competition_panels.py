#!/usr/bin/env python3
"""Build four accession-bound PhaDED competition phylogeny input panels.

The output is intentionally limited to candidate architecture/homology evidence.
Each panel has separate candidate/reference FASTA files plus a combined FASTA for
tree tools; candidate roles are retained in headers and a membership ledger.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable


PANEL_LABELS = {
    "extracellular_explicit": "A_retain_high_priority_architecture",
    "phaC_vs_phaZ": "E_hold_synthase_like_structure",
    "structure_anomaly": "E_hold_structure_assignment_review",
}
EXPECTED_COUNTS = {
    "extracellular_explicit": 4,
    "phaC_vs_phaZ": 12,
    "structure_anomaly": 8,
    "ephaz_competition": 18,
}
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
REQUIRED_CANDIDATE_COLUMNS = {"accession", "phaded_family_best", "phaded_superfamily_best"}
REQUIRED_REFERENCE_COLUMNS = {"reference_id", "accession", "phaded_family_id", "phaded_superfamily"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"input is not a regular file: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_candidate_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Read candidate FASTA keyed by the first header token."""
    records: dict[str, tuple[str, str]] = {}
    current: str | None = None
    header = ""
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                _store_sequence(records, current, header, chunks)
            header = line[1:].strip()
            current = header.split(None, 1)[0]
            chunks = []
        else:
            if current is None:
                raise ValueError(f"sequence before FASTA header: {path}")
            chunks.append(line)
    if current is not None:
        _store_sequence(records, current, header, chunks)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _store_sequence(records: dict[str, tuple[str, str]], key: str, header: str, chunks: list[str]) -> None:
    if key in records or not chunks:
        raise ValueError(f"duplicate or empty FASTA record: {key}")
    sequence = "".join(chunks).upper()
    if sequence.endswith("*"):
        sequence = sequence.rstrip("*")
    records[key] = (header, sequence)


def read_reference_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Read PhaDED FASTA keyed by reference_id, validating accession binding."""
    records: dict[str, tuple[str, str]] = {}
    current: str | None = None
    accession = ""
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                _store_reference(records, current, accession, chunks)
            fields = line[1:].strip().split("|", 1)
            if len(fields) != 2 or not fields[0] or not fields[1]:
                raise ValueError("reference FASTA header must be >reference_id|accession")
            current, accession, chunks = fields[0], fields[1], []
        else:
            if current is None:
                raise ValueError(f"sequence before FASTA header: {path}")
            chunks.append(line)
    if current is not None:
        _store_reference(records, current, accession, chunks)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _store_reference(records: dict[str, tuple[str, str]], key: str, accession: str, chunks: list[str]) -> None:
    if key in records or not chunks:
        raise ValueError(f"duplicate or empty reference FASTA record: {key}")
    sequence = "".join(chunks).upper()
    if sequence.endswith("*"):
        sequence = sequence.rstrip("*")
    records[key] = (accession, sequence)


def _validated_candidate_sequence(accession: str, sequence: str) -> str:
    if not sequence or any(residue not in AMINO_ACIDS for residue in sequence):
        raise ValueError(f"invalid amino acid in selected candidate: {accession}")
    return sequence


def _normal_family(value: str) -> str | None:
    match = re.search(r"(?:family_)?(DED_hfam_\d+)(?:_|$)", value or "")
    return match.group(1) if match else None


def _load_reference_rows(ledger_path: Path, fasta_path: Path) -> tuple[list[dict[str, str]], dict[str, tuple[str, str]]]:
    rows = read_tsv(ledger_path)
    if not rows:
        raise ValueError("reference ledger is empty")
    missing = sorted(REQUIRED_REFERENCE_COLUMNS - set(rows[0]))
    if missing:
        raise ValueError("reference ledger missing required columns: " + ",".join(missing))
    fasta = read_reference_fasta(fasta_path)
    seen: set[str] = set()
    for row in rows:
        reference_id = row.get("reference_id", "").strip()
        if not reference_id or reference_id in seen:
            raise ValueError(f"duplicate or missing reference_id: {reference_id!r}")
        if reference_id not in fasta:
            raise ValueError(f"reference missing from FASTA: {reference_id}")
        if fasta[reference_id][0] != row.get("accession", "").strip():
            raise ValueError(f"reference accession mismatch: {reference_id}")
        seen.add(reference_id)
    extra = set(fasta) - seen
    if extra:
        raise ValueError(f"FASTA contains references absent from ledger: {len(extra)}")
    return rows, fasta


def _pick_reference_ids(row: dict[str, str], reference_rows: Iterable[dict[str, str]]) -> tuple[list[str], list[str]]:
    family = _normal_family(row.get("phaded_family_best", ""))
    superfamily = row.get("phaded_superfamily_best", "").strip()
    matching = [r for r in reference_rows if family and r.get("phaded_family_id", "").strip() == family]
    if not matching and superfamily:
        matching = [r for r in reference_rows if r.get("phaded_superfamily", "").strip() == superfamily]
    if not matching:
        raise ValueError(f"no PhaDED references match candidate {row.get('accession', '')}")
    # Keep source-ledger order so the resulting input is deterministic and auditable.
    return [r["reference_id"].strip() for r in matching], []


def _candidate_row_maps(priority_path: Path, dossier_path: Path, interpro_path: Path) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    priority = read_tsv(priority_path)
    if not priority:
        raise ValueError("priority candidate set is empty")
    missing = sorted(REQUIRED_CANDIDATE_COLUMNS - set(priority[0]))
    if missing:
        raise ValueError("priority candidate set missing required columns: " + ",".join(missing))
    priority_map: dict[str, dict[str, str]] = {}
    for row in priority:
        accession = row.get("accession", "").strip()
        if not accession or accession in priority_map:
            raise ValueError(f"duplicate or missing priority accession: {accession!r}")
        priority_map[accession] = row
    dossier = read_tsv(dossier_path)
    if not dossier or "provisional_label" not in dossier[0] or "review_order" not in dossier[0]:
        raise ValueError("conflicting dossier missing provisional_label/review_order")
    interpro = read_tsv(interpro_path)
    if not interpro or "review_stage" not in interpro[0]:
        raise ValueError("InterPro priority table missing review_stage")
    selected: dict[str, dict[str, str]] = {}
    for row in dossier:
        accession = row.get("accession", "").strip()
        if accession in priority_map:
            merged = dict(priority_map[accession]); merged.update(row)
            selected[accession] = merged
    for row in interpro:
        accession = row.get("accession", "").strip()
        if row.get("review_stage", "").strip() == "ePhaZ_competition_review" and accession in priority_map:
            merged = dict(priority_map[accession]); merged.update(row)
            selected.setdefault(accession, merged)
    return priority_map, selected


def _write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


def build_panels(
    priority_path: Path,
    dossier_path: Path,
    interpro_path: Path,
    candidate_fasta_path: Path,
    reference_fasta_path: Path,
    reference_ledger_path: Path,
    output_dir: Path,
    *,
    expected_counts: dict[str, int] | None = None,
) -> dict[str, object]:
    priority_map, selected_map = _candidate_row_maps(priority_path, dossier_path, interpro_path)
    candidate_fasta = read_candidate_fasta(candidate_fasta_path)
    reference_rows, reference_fasta = _load_reference_rows(reference_ledger_path, reference_fasta_path)

    panels: dict[str, list[dict[str, str]]] = {name: [] for name in (*PANEL_LABELS, "ephaz_competition")}
    for accession, row in selected_map.items():
        label = row.get("provisional_label", "").strip()
        if label in PANEL_LABELS.values():
            panel = next(name for name, value in PANEL_LABELS.items() if value == label)
            panels[panel].append(row)
        if row.get("review_stage", "").strip() == "ePhaZ_competition_review":
            panels["ephaz_competition"].append(row)
    for panel in panels:
        panels[panel].sort(key=lambda row: (int(row.get("review_order", "0") or 0), row.get("accession", "")))
    seen: dict[str, str] = {}
    for panel, rows in panels.items():
        for row in rows:
            accession = row["accession"].strip()
            if accession in seen:
                raise ValueError(f"candidate accession appears in multiple panels: {accession} ({seen[accession]}, {panel})")
            seen[accession] = panel
            if accession not in candidate_fasta:
                raise ValueError(f"candidate missing from FASTA: {accession}")
            candidate_fasta[accession] = (candidate_fasta[accession][0], _validated_candidate_sequence(accession, candidate_fasta[accession][1]))
    if expected_counts:
        for panel, expected in expected_counts.items():
            observed = len(panels.get(panel, []))
            if observed != expected:
                raise ValueError(f"{panel} candidate count {observed} != expected {expected}")

    output_dir.mkdir(parents=True, exist_ok=True)
    panel_manifest: dict[str, object] = {}
    for panel, rows in panels.items():
        candidate_records: list[tuple[str, str]] = []
        reference_records: list[tuple[str, str]] = []
        members: list[dict[str, str]] = []
        reference_ids: list[str] = []
        role = "candidate_positive" if panel == "extracellular_explicit" else "candidate_challenge"
        for row in rows:
            accession = row["accession"].strip()
            candidate_records.append((f"candidate|{panel}|{accession}", candidate_fasta[accession][1]))
            members.append({"panel": panel, "role": role, "member_id": accession, "source_accession": accession, "family": row.get("phaded_family_best", ""), "superfamily": row.get("phaded_superfamily_best", "")})
            targets, _ = _pick_reference_ids(row, reference_rows)
            for reference_id in targets:
                if reference_id in reference_ids:
                    continue
                reference_ids.append(reference_id)
                accession_ref, sequence = reference_fasta[reference_id]
                if not sequence or any(residue not in AMINO_ACIDS for residue in sequence):
                    raise ValueError(f"invalid amino acid in selected reference: {reference_id}")
                reference_records.append((f"reference_target|{panel}|{reference_id}|{accession_ref}", sequence))
                refrow = next(item for item in reference_rows if item["reference_id"] == reference_id)
                members.append({"panel": panel, "role": "reference_target", "member_id": reference_id, "source_accession": accession_ref, "family": refrow.get("phaded_family_id", ""), "superfamily": refrow.get("phaded_superfamily", "")})
        combined = candidate_records + reference_records
        candidate_path = output_dir / f"panel_{panel}.candidate.faa"
        reference_path = output_dir / f"panel_{panel}.references.faa"
        combined_path = output_dir / f"panel_{panel}.phylo.faa"
        members_path = output_dir / f"panel_{panel}.members.tsv"
        _write_fasta(candidate_path, candidate_records)
        _write_fasta(reference_path, reference_records)
        _write_fasta(combined_path, combined)
        with members_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["panel", "role", "member_id", "source_accession", "family", "superfamily"], delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(members)
        panel_manifest[panel] = {
            "candidate_count": len(candidate_records),
            "reference_target_count": len(reference_records),
            "leaf_count": len(combined),
            "candidate_fasta": file_record(candidate_path),
            "reference_fasta": file_record(reference_path),
            "phylo_fasta": file_record(combined_path),
            "membership_tsv": file_record(members_path),
            "challenge_control_policy": "candidate_positive only for explicit extracellular bucket; all other candidates are candidate_challenge; reference_target records are kept separate in their own FASTA.",
        }
    manifest = {
        "schema_version": 1,
        "status": "candidate-only",
        "inputs": {name: file_record(path) for name, path in {
            "priority_candidate_set_v2": priority_path,
            "conflicting_review_dossier": dossier_path,
            "interpro_priority_candidates": interpro_path,
            "candidate_union_fasta": candidate_fasta_path,
            "phaded_reference_fasta": reference_fasta_path,
            "phaded_reference_ledger": reference_ledger_path,
        }.items()},
        "panels": panel_manifest,
        "candidate_counts": {panel: len(rows) for panel, rows in panels.items()},
        "phenotype_boundary": "Phylogenetic placement, domains, motifs, and localization are candidate homology/architecture evidence and do not validate PHB/PHA degradation phenotype.",
    }
    (output_dir / "competition_panel_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--priority-set", type=Path, required=True)
    parser.add_argument("--dossier", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-ledger", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = build_panels(
        args.priority_set, args.dossier, args.interpro, args.candidate_fasta,
        args.reference_fasta, args.reference_ledger, args.output_dir,
        expected_counts=EXPECTED_COUNTS,
    )
    print(json.dumps({"status": manifest["status"], "candidate_counts": manifest["candidate_counts"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
