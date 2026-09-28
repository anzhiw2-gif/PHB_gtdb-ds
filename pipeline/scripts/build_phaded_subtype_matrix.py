#!/usr/bin/env python3
"""Build a conservative accession-level PhaDED subtype evidence matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping


PROFILE_FAMILY_RE = re.compile(r"(DED_hfam_\d+)")
EXPORT_CLASSES = {"SP", "LIPO", "TAT", "TATLIPO"}
NEW_FIELDS = [
    "profile_model_status", "profile_evidence_status",
    "domain_evidence_status", "domain_pfam_state", "domain_interpro_state",
    "motif_lipase_box", "motif_catalytic_ser_cys", "motif_his", "motif_asp",
    "motif_oxyanion_hole", "motif_sbd", "motif_linker", "motif_lid",
    "motif_evidence_status", "motif_completeness",
    "localization_evidence_status", "structure_evidence_status_v2",
    "phylogeny_evidence_status_v2", "manual_review_status",
    "subtype_call", "subtype_confidence", "subtype_evidence_boundary",
]


def _index(rows: Iterable[Mapping[str, str]], key: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        value = str(row.get(key, "")).strip()
        if not value:
            continue
        if value in result:
            raise ValueError(f"duplicate accession in evidence table: {value}")
        result[value] = dict(row)
    return result


def _family_id(value: str) -> str:
    match = PROFILE_FAMILY_RE.search(value or "")
    return match.group(1) if match else ""


def _profile_status(assignment: Mapping[str, str], profiles: list[Mapping[str, str]]) -> tuple[str, str]:
    status = assignment.get("assignment_status", "")
    if status == "unassigned_PhaDED_like":
        return "unassigned", "profile_unassigned"
    if status == "ambiguous_superfamily":
        return "ambiguous", "profile_ambiguous_superfamily"
    if status == "ambiguous_family":
        return "ambiguous", "profile_ambiguous_family"
    family = _family_id(assignment.get("phaded_family_best", ""))
    superfamily = assignment.get("phaded_superfamily_best", "")
    matches = [
        row for row in profiles
        if (family and row.get("phaded_family_id", "") == family)
        or (superfamily and row.get("phaded_superfamily", "") == superfamily)
    ]
    if any(row.get("model_status") == "trained" for row in matches):
        return "trained", "profile_trained_hit"
    if matches:
        return "reference_only", "profile_reference_only_hit"
    return "unknown", "profile_model_binding_pending"


def _domain_status(pfam: Mapping[str, str], interpro: Mapping[str, str]) -> tuple[str, str, str]:
    pfam_state = pfam.get("architecture_consistency", "pending") or "pending"
    interpro_state = interpro.get("interpro_secondary_state", "no_interpro_hit") or "no_interpro_hit"
    if pfam_state == "conflicting":
        return "domain_conflict", pfam_state, interpro_state
    if pfam_state == "partial":
        return "domain_supported_partial", pfam_state, interpro_state
    if interpro_state == "explicit_pha_related_domain":
        return "domain_supported_explicit_interpro", pfam_state, interpro_state
    if interpro.get("interpro_status") == "interpro_supported":
        return "domain_supported_generic_interpro", pfam_state, interpro_state
    if interpro.get("interpro_status") == "interpro_tool_input_excluded":
        return "domain_not_tested_tool_input_excluded", pfam_state, interpro_state
    return "domain_pending", pfam_state, interpro_state


def _motif_status(feature: Mapping[str, str]) -> tuple[dict[str, str], str, str]:
    if feature.get("motif_panel_status"):
        fields = {
            "motif_lipase_box": feature.get("lipase_box_state", "pending") or "pending",
            "motif_catalytic_ser_cys": feature.get("catalytic_ser_cys_state", "pending") or "pending",
            "motif_his": feature.get("his_state", "pending") or "pending",
            "motif_asp": feature.get("asp_state", "pending") or "pending",
            "motif_oxyanion_hole": feature.get("oxyanion_hole_state", "pending") or "pending",
            "motif_sbd": feature.get("sbd_state", "pending") or "pending",
            "motif_linker": feature.get("linker_state", "pending") or "pending",
            "motif_lid": feature.get("lid_state", "pending") or "pending",
        }
        return fields, str(feature["motif_panel_status"]), str(feature.get("motif_completeness", "pattern_partial"))
    fields = {
        "motif_lipase_box": feature.get("lipase_box_state", "pending") or "pending",
        "motif_catalytic_ser_cys": feature.get("catalytic_ser_cys_state", "pending") or "pending",
        "motif_his": feature.get("his_state", "pending") or "pending",
        "motif_asp": feature.get("asp_state", "pending") or "pending",
        "motif_oxyanion_hole": feature.get("oxyanion_hole_state", "pending") or "pending",
        "motif_sbd": feature.get("sbd_state", "pending") or "pending",
        "motif_linker": feature.get("linker_state", "pending") or "pending",
        "motif_lid": feature.get("lid_state", "pending") or "pending",
    }
    catalytic = {fields[key] for key in ("motif_lipase_box", "motif_catalytic_ser_cys", "motif_his", "motif_asp", "motif_oxyanion_hole")}
    if catalytic == {"supported"}:
        return fields, "complete_catalytic_motif", "complete"
    if fields["motif_lipase_box"] == "present" and fields["motif_catalytic_ser_cys"] == "supported":
        return fields, "partial_lipase_box_only", "partial"
    if all(value in {"not_expected_for_superfamily", "not_detected"} for value in catalytic):
        return fields, "motif_not_detected", "not_supported"
    return fields, "not_tested_full_library", "not_tested"


def _localization_status(feature: Mapping[str, str]) -> str:
    evidence = feature.get("signalp_evidence", "")
    signal_class = feature.get("signalp_class", "")
    if evidence in {"accepted_layer_signalp", "signalp_supported"} and signal_class in EXPORT_CLASSES:
        return "export_signal_supported"
    if evidence in {"accepted_layer_signalp", "signalp_supported"}:
        return "localization_unassigned"
    if evidence == "signalp_tool_input_excluded":
        return "localization_tool_input_excluded"
    if evidence == "signalp_not_reported":
        return "localization_not_reported"
    return "not_tested_full_library"


def _structure_status(row: Mapping[str, str] | None) -> str:
    if not row:
        return "not_tested_full_library"
    if row.get("architecture_evidence") == "conflicting" or "architecture_conflict" in row.get("review_tier", ""):
        return "structure_conflict_candidate_only"
    if row.get("review_tier") == "high_confidence_candidate_only" or row.get("structural_evidence_band") in {"moderate", "moderate_high"}:
        return "structure_support_candidate_only"
    return "structure_review_candidate_only"


def _phylogeny_status(row: Mapping[str, str] | None) -> str:
    if not row:
        return "not_tested_full_library"
    decision = row.get("phylo_decision", "") or row.get("phylogeny_decision", "")
    if decision.startswith("retain_"):
        return "phylogeny_support_candidate_only"
    if decision.startswith("unresolved"):
        return "phylogeny_unresolved_candidate_only"
    return "phylogeny_review_candidate_only"


def _base_subtype(assignment: Mapping[str, str]) -> str:
    superfamily = assignment.get("phaded_superfamily_best", "")
    if "extracellular dPHASCL type 1" in superfamily:
        return "dPHASCL1_like_candidate"
    if "extracellular dPHASCL type 2" in superfamily:
        return "dPHASCL2_like_candidate"
    if "extracellular dPHAMCL" in superfamily:
        return "dPHAMCL_like_candidate"
    if "intracellular nPHAMCL" in superfamily:
        return "nPHAMCL_like_candidate"
    if "intracellular nPHASCL" in superfamily:
        return "nPHASCL_like_candidate"
    return "unresolved_PhaDED_like"


def build_matrix(
    assignments: list[Mapping[str, str]],
    features: list[Mapping[str, str]],
    pfam: list[Mapping[str, str]],
    interpro: list[Mapping[str, str]],
    profiles: list[Mapping[str, str]],
    structure: list[Mapping[str, str]] | None = None,
    phylogeny: list[Mapping[str, str]] | None = None,
    manual: list[Mapping[str, str]] | None = None,
    base_ledger: list[Mapping[str, str]] | None = None,
) -> list[dict[str, str]]:
    assignments_by = _index(assignments, "accession")
    features_by = _index(features, "accession")
    pfam_by = _index(pfam, "accession")
    interpro_by = _index(interpro, "accession")
    structure_by = _index(structure or [], "source_accession")
    phylo_by = _index(phylogeny or [], "source_accession")
    manual_by = _index(manual or [], "candidate_id")
    base_by = _index(base_ledger or [], "accession")
    if set(assignments_by) != set(features_by):
        raise ValueError("assignment and feature accession sets differ")
    if not set(pfam_by).issubset(assignments_by) or not set(interpro_by).issubset(assignments_by):
        raise ValueError("evidence table contains accession outside assignments")

    output: list[dict[str, str]] = []
    for accession, assignment in assignments_by.items():
        feature = features_by[accession]
        pfam_row = pfam_by.get(accession, {})
        interpro_row = interpro_by.get(accession, {"interpro_status": "interpro_not_reported", "interpro_secondary_state": "no_interpro_hit"})
        profile_model_status, profile_evidence = _profile_status(assignment, profiles)
        domain_evidence, pfam_state, interpro_state = _domain_status(pfam_row, interpro_row)
        motifs, motif_status, motif_completeness = _motif_status(feature)
        structure_status = _structure_status(structure_by.get(accession))
        phylogeny_status = _phylogeny_status(phylo_by.get(accession))
        manual_row = manual_by.get(accession, {})
        manual_decision = manual_row.get("manual_review_decision", "")
        manual_status = "not_in_19_candidate_review"
        if manual_decision.startswith("hold_"):
            manual_status = manual_decision
        elif manual_decision:
            manual_status = "retain_candidate_only"

        if manual_decision == "hold_architecture_conflict" or domain_evidence == "domain_conflict":
            subtype_call, confidence = "hold_architecture_conflict", "hold"
        elif manual_decision == "hold_gene_model_or_structure":
            subtype_call, confidence = "hold_gene_model_or_structure", "hold"
        elif assignment.get("assignment_status") == "unassigned_PhaDED_like":
            subtype_call, confidence = "unassigned_PhaDED_like", "unresolved"
        elif assignment.get("assignment_status") == "ambiguous_superfamily":
            subtype_call, confidence = "ambiguous_superfamily", "unresolved"
        elif assignment.get("assignment_status") == "ambiguous_family":
            subtype_call, confidence = "ambiguous_family_within_superfamily", "unresolved"
        else:
            subtype_call = _base_subtype(assignment)
            if profile_model_status == "trained" and domain_evidence.startswith("domain_supported") and motif_completeness == "complete":
                confidence = "high_candidate_only"
            elif profile_model_status == "trained" and domain_evidence.startswith("domain_supported"):
                confidence = "moderate_candidate_only"
            elif profile_model_status in {"trained", "reference_only"}:
                confidence = "low_candidate_only"
            else:
                confidence = "unresolved"

        row = dict(base_by.get(accession, {}))
        row.update(assignment)
        row.update(motifs)
        row.update({
            "profile_model_status": profile_model_status,
            "profile_evidence_status": profile_evidence,
            "domain_evidence_status": domain_evidence,
            "domain_pfam_state": pfam_state,
            "domain_interpro_state": interpro_state,
            "architecture_consistency": pfam_state,
            "motif_evidence_status": motif_status,
            "motif_completeness": motif_completeness,
            "localization_evidence_status": _localization_status(feature),
            "structure_evidence_status_v2": structure_status,
            "phylogeny_evidence_status_v2": phylogeny_status,
            "manual_review_status": manual_status,
            "subtype_call": subtype_call,
            "subtype_confidence": confidence,
            "subtype_evidence_boundary": "candidate_only; homology/domain/motif/localization/structure/phylogeny evidence is not phenotype validation",
            "interpro_status": interpro_row.get("interpro_status", "interpro_not_reported"),
            "interpro_evidence_status": interpro_row.get("interpro_status", "interpro_not_reported"),
            "interpro_secondary_state": interpro_row.get("interpro_secondary_state", "no_interpro_hit"),
            "interpro_hit_count": interpro_row.get("interpro_hit_count", "0"),
            "interpro_signatures": interpro_row.get("interpro_signatures", ""),
            "interpro_ids": interpro_row.get("interpro_ids", ""),
        })
        output.append(row)
    return output


def read_tsv(path: Path, fieldnames: list[str] | None = None) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, fieldnames=fieldnames, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_matrix(args: argparse.Namespace) -> dict[str, object]:
    assignment_rows = read_tsv(args.assignments)
    feature_rows = read_tsv(args.features)
    pfam_rows = read_tsv(args.pfam)
    interpro_rows = read_tsv(args.interpro)
    profile_rows = read_tsv(args.profiles)
    base_rows = read_tsv(args.base_ledger) if args.base_ledger else []
    optional = lambda path, key: read_tsv(path) if path else []
    rows = build_matrix(assignment_rows, feature_rows, pfam_rows, interpro_rows, profile_rows, optional(args.structure, "structure"), optional(args.phylogeny, "phylogeny"), optional(args.manual, "manual"), base_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    counts = {field: dict(Counter(row.get(field, "") for row in rows)) for field in ("subtype_call", "subtype_confidence", "profile_evidence_status", "domain_evidence_status", "motif_evidence_status", "localization_evidence_status", "structure_evidence_status_v2", "phylogeny_evidence_status_v2", "interpro_status")}
    manifest = {
        "status": "completed_candidate_only",
        "candidate_count": len(rows),
        "counts": counts,
        "inputs": {},
        "output": {"path": str(args.output.resolve()), "size": args.output.stat().st_size, "sha256": sha256(args.output)},
        "phenotype_boundary": "Subtype calls are computational candidate classifications and do not validate PHB/PHA degradation phenotype.",
    }
    for name, path in (("assignments", args.assignments), ("features", args.features), ("pfam", args.pfam), ("interpro", args.interpro), ("profiles", args.profiles), ("base_ledger", args.base_ledger), ("structure", args.structure), ("phylogeny", args.phylogeny), ("manual", args.manual)):
        if path:
            manifest["inputs"][name] = {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}
    manifest_path = args.output.parent / "subtype_matrix_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("assignments", "features", "pfam", "interpro", "profiles", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    for name in ("structure", "phylogeny", "manual", "base_ledger"):
        parser.add_argument(f"--{name}", type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(write_matrix(args), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
