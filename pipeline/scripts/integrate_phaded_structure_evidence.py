#!/usr/bin/env python3
"""Integrate structure, corrected feature, Pfam and phylogeny evidence.

The output is a candidate-only review table.  No field in this script is a
phenotype call or a training-set authorization.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Iterable


PHENOTYPE_BOUNDARY = "candidate_only"
COORDINATE_RE = re.compile(r":(\d+)-(\d+)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def index_unique(rows: Iterable[dict[str, str]], key: str, label: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        value = row.get(key, "")
        if not value:
            raise ValueError(f"{label}: empty {key}")
        if value in indexed:
            raise ValueError(f"{label}: duplicate {key}={value}")
        indexed[value] = row
    return indexed


def structure_quality(band: str) -> str:
    return {
        "high_structure_confidence": "high",
        "intermediate_structure_confidence": "intermediate",
        "low_structure_confidence": "low",
    }.get(band, "unknown")


def summarize_pdb(path: Path) -> dict[str, str]:
    """Summarize per-residue pLDDT stored in the PDB B-factor column."""
    values: list[float] = []
    for line in path.read_text(encoding="ascii", errors="strict").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            values.append(float(line[60:66]))
    if not values:
        raise ValueError(f"no CA residues found in {path}")
    low = [value < 70.0 for value in values]
    longest = current = 0
    for is_low in low:
        current = current + 1 if is_low else 0
        longest = max(longest, current)
    edge = min(20, len(values))
    return {
        "residue_count": str(len(values)),
        "mean_plddt": f"{sum(values) / len(values):.2f}",
        "low_fraction": f"{sum(low) / len(values):.3f}",
        "longest_low_run": str(longest),
        "n_terminal_20_mean": f"{sum(values[:edge]) / edge:.2f}",
        "c_terminal_20_mean": f"{sum(values[-edge:]) / edge:.2f}",
    }


def pfam_span(coordinates: str, length: str) -> dict[str, str]:
    intervals = [(int(start), int(end)) for start, end in COORDINATE_RE.findall(coordinates or "")]
    if not intervals or not length:
        return {"pfam_span_start": "", "pfam_span_end": "", "pfam_span_coverage": ""}
    start = min(item[0] for item in intervals)
    end = max(item[1] for item in intervals)
    return {"pfam_span_start": str(start), "pfam_span_end": str(end), "pfam_span_coverage": f"{(end - start + 1) / int(length):.3f}"}


def build_review_row(
    structure: dict[str, str],
    feature: dict[str, str],
    pfam: dict[str, str],
    phylo: dict[str, str],
    interpro: dict[str, str] | None = None,
) -> dict[str, str]:
    interpro = interpro or {}
    quality = structure_quality(structure.get("structure_confidence_band", ""))
    architecture = pfam.get("architecture_consistency", "pending") or "pending"
    integrity = feature.get("sequence_integrity", "pending") or "pending"
    if quality == "low" or integrity == "possible_N_truncation":
        tier = "low_confidence_or_truncation_hold"
        action = "hold_for_sequence_and_structure_review"
    elif architecture == "conflicting":
        tier = "architecture_conflict_manual_review"
        action = "inspect_pfam_coordinates_and_structure_domain_boundary"
    elif quality == "intermediate":
        tier = "intermediate_structure_manual_review"
        action = "inspect_fold_completeness_and_domain_boundary"
    else:
        tier = "high_confidence_candidate_only"
        action = "retain_candidate_only_for_downstream_comparison"

    support = float(phylo.get("candidate_reference_mrca_support") or 0)
    family_status = pfam.get("assignment_status", "")
    family_gap = float(pfam.get("family_score_gap") or 0)
    conflict_codes = []
    if architecture == "conflicting" or pfam.get("pfam_interpro_state") == "tested_other_superfamily_hit":
        conflict_codes.append("PFAM_OTHER_SUPERFAMILY")
    if family_status == "ambiguous" or family_gap < 1.0:
        conflict_codes.append("FAMILY_ASSIGNMENT_AMBIGUOUS")
    if integrity == "possible_N_truncation":
        conflict_codes.append("POSSIBLE_N_TRUNCATION")
    if quality == "low":
        conflict_codes.append("LOW_STRUCTURE")
    if feature.get("signalp_class", "OTHER") == "OTHER":
        conflict_codes.append("SIGNALP_NOT_IN_LAYER_CONTRACT")
    if integrity == "terminal_stop_only":
        conflict_codes.append("TERMINAL_STOP_CLEANED")
    hard_conflict = "PFAM_OTHER_SUPERFAMILY" in conflict_codes or "FAMILY_ASSIGNMENT_AMBIGUOUS" in conflict_codes
    if hard_conflict:
        merged_priority = "P0_architecture_or_integrity_hold"
    elif "LOW_STRUCTURE" in conflict_codes:
        merged_priority = "P3_low_structure_hold"
    elif "POSSIBLE_N_TRUNCATION" in conflict_codes:
        merged_priority = "P2_truncation_hold"
    elif quality == "high" and support >= 90:
        merged_priority = "P1_high_architecture_phylo_candidate"
    else:
        merged_priority = "P2_intermediate_candidate"

    return {
        "source_accession": structure["source_accession"],
        "panel": structure.get("panel", ""),
        "superfamily": structure.get("superfamily", ""),
        "phylo_decision": structure.get("phylo_decision", ""),
        "nearest_reference_family": phylo.get("nearest_reference_family", ""),
        "nearest_reference_accession": phylo.get("nearest_reference_accession", ""),
        "nearest_reference_distance": phylo.get("nearest_reference_distance", ""),
        "mrca_support": phylo.get("candidate_reference_mrca_support", ""),
        "sequence_length_cleaned": structure.get("sequence_length_cleaned", ""),
        "best_model": structure.get("best_model", ""),
        "best_model_mean_plddt": structure.get("best_model_mean_plddt", ""),
        "best_model_ptm": structure.get("best_model_ptm", ""),
        "best_model_max_pae": structure.get("best_model_max_pae", ""),
        "structure_quality": quality,
        "sequence_integrity": integrity,
        "signalp_class": feature.get("signalp_class", "pending") or "pending",
        "lipase_box_state": feature.get("lipase_box_state", "pending") or "pending",
        "catalytic_ser_cys_state": feature.get("catalytic_ser_cys_state", "pending") or "pending",
        "his_state": feature.get("his_state", "pending") or "pending",
        "asp_state": feature.get("asp_state", "pending") or "pending",
        "oxyanion_hole_state": feature.get("oxyanion_hole_state", "pending") or "pending",
        "pfam_accessions": pfam.get("pfam_accessions", "") or "none_recorded",
        "pfam_coordinates": pfam.get("pfam_coordinates", "") or "",
        **pfam_span(pfam.get("pfam_coordinates", ""), structure.get("sequence_length_cleaned", "")),
        "phaded_superfamily_best": pfam.get("phaded_superfamily_best", ""),
        "phaded_superfamily_second": pfam.get("phaded_superfamily_second", ""),
        "superfamily_score_gap": pfam.get("superfamily_score_gap", ""),
        "phaded_family_best": pfam.get("phaded_family_best", ""),
        "phaded_family_second": pfam.get("phaded_family_second", ""),
        "architecture_evidence": architecture,
        "assignment_review": pfam.get("assignment_review", "pending") or "pending",
        "assignment_status": pfam.get("assignment_status", "pending") or "pending",
        "family_score_gap": pfam.get("family_score_gap", ""),
        "interpro_secondary_state": interpro.get("interpro_secondary_state", "pending") or "pending",
        "interpro_explicit_ids": interpro.get("interpro_explicit_ids", "") or "",
        "interpro_signatures": interpro.get("interpro_signatures", "") or "",
        "interpro_hit_count": interpro.get("interpro_hit_count", "") or "",
        "review_tier": tier,
        "merged_priority": merged_priority,
        "conflict_codes": ";".join(conflict_codes) or "none_recorded",
        "next_action": action,
        "phenotype_claim": PHENOTYPE_BOUNDARY,
    }


def integrate(structure_path: Path, feature_path: Path, pfam_path: Path, phylo_path: Path, output_dir: Path, pdb_dir: Path | None = None, interpro_path: Path | None = None) -> dict[str, object]:
    structures = read_tsv(structure_path)
    features = index_unique(read_tsv(feature_path), "accession", "features")
    pfams = index_unique(read_tsv(pfam_path), "accession", "pfam")
    interpros = index_unique(read_tsv(interpro_path), "accession", "interpro") if interpro_path is not None else {}
    phylo_rows = read_tsv(phylo_path)
    # One candidate can occur in multiple panels.  Prefer the explicit
    # extracellular panel when present, otherwise use the ePhaZ competition panel.
    phylo_by_accession: dict[str, dict[str, str]] = {}
    panel_rank = {"extracellular_explicit": 0, "ephaz_competition": 1}
    for row in phylo_rows:
        accession = row.get("source_accession", "")
        if not accession:
            continue
        current = phylo_by_accession.get(accession)
        if current is None or panel_rank.get(row.get("panel", ""), 99) < panel_rank.get(current.get("panel", ""), 99):
            phylo_by_accession[accession] = row

    rows = []
    for structure in structures:
        accession = structure.get("source_accession", "")
        if not accession:
            raise ValueError("structure summary contains an empty source_accession")
        if accession not in features:
            raise ValueError(f"missing feature row for {accession}")
        if accession not in pfams:
            raise ValueError(f"missing Pfam row for {accession}")
        if accession not in phylo_by_accession:
            raise ValueError(f"missing phylogeny row for {accession}")
        rows.append(build_review_row(structure, features[accession], pfams[accession], phylo_by_accession[accession], interpros.get(accession)))

    if len({row["source_accession"] for row in rows}) != len(rows):
        raise ValueError("structure summary must contain unique candidate accessions")
    output_dir.mkdir(parents=True, exist_ok=True)
    if pdb_dir is not None:
        for row in rows:
            pdb_name = row["best_model"].replace("_scores_rank_", "_unrelaxed_rank_").replace(".json", ".pdb")
            metrics = summarize_pdb(pdb_dir / pdb_name)
            row.update({f"pdb_{key}": value for key, value in metrics.items()})

    table_path = output_dir / "phaded_structure_evidence_review.tsv"
    fields = list(rows[0]) if rows else []
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary_path = output_dir / "structure_evidence_review_summary.tsv"
    counts = Counter(row["review_tier"] for row in rows)
    priority_counts = Counter(row["merged_priority"] for row in rows)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["review_tier", "candidate_count"])
        writer.writerows(sorted(counts.items()))

    pdb_models = []
    if pdb_dir is not None:
        for row in rows:
            pdb_name = row["best_model"].replace("_scores_rank_", "_unrelaxed_rank_").replace(".json", ".pdb")
            pdb_path = pdb_dir / pdb_name
            pdb_models.append({"path": str(pdb_path.resolve()), "size": pdb_path.stat().st_size, "sha256": sha256(pdb_path)})

    manifest = {
        "schema_version": "1.0",
        "status": "completed_candidate_only",
        "candidate_count": len(rows),
        "inputs": {
            name: {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}
            for name, path in {
                "structure_summary": structure_path,
                "corrected_feature_annotation": feature_path,
                "pfam_evidence": pfam_path,
                "phylogeny_decisions": phylo_path,
                **({"interpro_secondary_evidence": interpro_path} if interpro_path is not None else {}),
            }.items()
        },
        "outputs": {
            "review_table": {"path": str(table_path.resolve()), "size": table_path.stat().st_size, "sha256": sha256(table_path)},
            "review_summary": {"path": str(summary_path.resolve()), "size": summary_path.stat().st_size, "sha256": sha256(summary_path)},
        },
        "pdb_confidence_source": ({"path": str(pdb_dir.resolve()), "status": "bound_directory"} if pdb_dir is not None else {"status": "not_requested"}),
        "pdb_models": pdb_models,
        "review_tier_counts": dict(sorted(counts.items())),
        "merged_priority_counts": dict(sorted(priority_counts.items())),
        "phenotype_boundary": "Predicted structures and existing sequence/domain/phylogeny evidence remain candidate-only and do not validate PHB/PHA degradation.",
        "training_boundary": "No candidate is added to ePhaZ_curated_core or any positive training set.",
    }
    manifest_path = output_dir / "structure_evidence_review_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    run_dir = output_dir.parent
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "inputs").mkdir(parents=True, exist_ok=True)
    contract = {
        "schema_version": "1.0",
        "run_id": run_dir.name,
        "status": "completed_candidate_only",
        "inputs": manifest["inputs"],
        "outputs": manifest["outputs"],
        "phenotype_boundary": manifest["phenotype_boundary"],
        "training_boundary": manifest["training_boundary"],
    }
    (run_dir / "input_contract.json").write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure-summary", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--pfam", type=Path, required=True)
    parser.add_argument("--phylo", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pdb-dir", type=Path)
    parser.add_argument("--interpro", type=Path)
    args = parser.parse_args()
    print(json.dumps(integrate(args.structure_summary, args.features, args.pfam, args.phylo, args.output_dir, args.pdb_dir, args.interpro), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
