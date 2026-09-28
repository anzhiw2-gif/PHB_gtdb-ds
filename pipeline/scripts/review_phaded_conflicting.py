#!/usr/bin/env python3
"""Produce a deterministic, candidate-only review dossier for Pfam conflicts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
TYPE1_MOTIF_OUTSIDE = {
    "GCA_041168165.1|JBEELK010000002.1_2",
    "GCA_049285915.1|JBAQAW010000024.1_32",
    "GCA_964643805.1|CAZXUD010000048.1_9",
    "GCF_039710355.1|NZ_JBDKWZ010000014.1_90",
}
TYPE1_ASSIGNMENT_UNCERTAIN = "GCA_964540355.1|CAZIMU010000007.1_4"
TYPE2_LOW_GAP = {
    "GCA_005787565.1|SYAP01000036.1_14",
    "GCF_021172045.1|NZ_CP087164.1_2925",
    "GCA_963731355.1|CAVSUS010000013.1_102",
}


def sequence_integrity(sequence: str) -> str:
    sequence = sequence.strip().upper()
    if sequence.endswith("*") and all(residue in AMINO_ACIDS for residue in sequence[:-1]):
        return "terminal_stop_only"
    if all(residue in AMINO_ACIDS for residue in sequence):
        return "valid"
    return "invalid_internal_character"


def _numeric(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _coordinates(value: str) -> list[tuple[int, int, str]]:
    result = []
    for token in filter(None, value.split(";")):
        match = re.match(r"([^:]+):(.*)$", token)
        if not match:
            continue
        pfam = match.group(1)
        for interval in match.group(2).split(","):
            bounds = interval.split("-")
            if len(bounds) == 2:
                try:
                    result.append((int(bounds[0]), int(bounds[1]), pfam))
                except ValueError:
                    pass
    return result


def motif_inside_pfam(feature: dict[str, str]) -> str:
    motif = feature.get("lipase_box_coordinates", "")
    if not motif or "-" not in motif:
        return "not_applicable"
    start, end = (int(value) for value in motif.split("-", 1))
    return "yes" if any(domain_start <= start and end <= domain_end for domain_start, domain_end, _ in _coordinates(feature.get("pfam_coordinates", ""))) else "no"


def provisional_label(row: dict[str, str]) -> tuple[str, str]:
    accession = row.get("accession", "")
    superfamily = row.get("phaded_superfamily_best", "")
    interpro = row.get("interpro_secondary_state", "")
    if interpro == "pha_synthase_like_alternative":
        return "E_hold_synthase_like_structure", "InterPro/PANTHER indicates PHA synthase-like evidence; do not treat as depolymerase."
    if interpro == "explicit_pha_related_domain" and superfamily == "extracellular dPHASCL type 2":
        return "A_retain_high_priority_architecture", "Explicit PHA-related domain evidence supports an extracellular PhaZ-like architecture, but type boundary remains unresolved."
    if accession in TYPE1_MOTIF_OUTSIDE or accession == TYPE1_ASSIGNMENT_UNCERTAIN or accession in TYPE2_LOW_GAP:
        return "E_hold_structure_assignment_review", "Catalytic motif/domain or PhaDED assignment confidence is unresolved; require structure or phylogenetic review."
    if superfamily == "extracellular dPHASCL type 1":
        return "B_cross_fold_architecture_review", "PhaDED type 1 assignment conflicts with a non-target alpha/beta-hydrolase fingerprint; retain as homolog candidate only."
    if superfamily == "extracellular dPHASCL type 2":
        return "B_family_unresolved_architecture_review", "Extracellular hydrolase/PhaZ-like architecture is plausible, but family/type assignment is unresolved."
    return "C_generic_or_unresolved", "Evidence does not support a resolved PhaDED family interpretation."


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = None
    chunks: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if accession is not None:
                records[accession] = "".join(chunks)
            accession, chunks = line[1:].split()[0], []
        elif line:
            chunks.append(line)
    if accession is not None:
        records[accession] = "".join(chunks)
    return records


def review(features_path: Path, secondary_path: Path, fasta_path: Path) -> list[dict[str, str]]:
    features = {row["accession"]: row for row in _read(features_path)}
    secondary = [row for row in _read(secondary_path) if row.get("review_stage") == "conflicting"]
    sequences = _read_fasta(fasta_path)
    if len(secondary) != 83:
        raise ValueError(f"expected 83 conflicting records, observed {len(secondary)}")
    output = []
    for row in secondary:
        feature = features.get(row["accession"])
        if feature is None:
            raise ValueError(f"missing feature row: {row['accession']}")
        merged = dict(row)
        merged["sequence_integrity_review"] = sequence_integrity(sequences.get(row["accession"], ""))
        merged["motif_inside_pfam"] = motif_inside_pfam(feature)
        merged["pfam_coordinates_review"] = feature.get("pfam_coordinates", "")
        merged["lipase_box_review"] = feature.get("lipase_box_coordinates", "")
        label, reason = provisional_label(merged)
        merged["provisional_label"] = label
        merged["provisional_reason"] = reason
        merged["next_action"] = {
            "A_retain_high_priority_architecture": "structure_or_experimental_shortlist",
            "B_cross_fold_architecture_review": "structure_or_phylogeny_before_experiment",
            "B_family_unresolved_architecture_review": "family_boundary_review_then_structure",
            "E_hold_synthase_like_structure": "separate_synthase_vs_depolymerase_review",
            "E_hold_structure_assignment_review": "structure_or_phylogeny_review",
        }.get(label, "retain_candidate_only")
        merged["phenotype_boundary"] = "candidate_only_not_phenotype_validation"
        output.append(merged)
    return output


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--secondary", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    rows = review(args.features, args.secondary, args.fasta)
    fields = sorted({key for row in rows for key in row})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    manifest = {"status": "completed_candidate_only", "conflicting_count": len(rows), "inputs": {"features": {"path": str(args.features.resolve()), "sha256": sha256(args.features)}, "secondary": {"path": str(args.secondary.resolve()), "sha256": sha256(args.secondary)}, "fasta": {"path": str(args.fasta.resolve()), "sha256": sha256(args.fasta)}}, "output": {"path": str(args.output.resolve()), "sha256": sha256(args.output)}, "phenotype_boundary": "This is a computational review dossier; no label establishes PHB/PHA degradation activity."}
    (args.output.parent / "conflicting_review_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "conflicting_count": len(rows)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
