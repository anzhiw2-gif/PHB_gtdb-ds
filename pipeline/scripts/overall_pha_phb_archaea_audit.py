#!/usr/bin/env python3
"""Combine existing PHA evidence layers with a conservative PHB focus flag.

This is a candidate-only reporting step. It does not infer phenotype from
homology and keeps archaeal PhaZh1-like hits in an independent audit table.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def _bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def classify_pha_row(row: dict[str, str], validated_phb: set[str] | None = None,
                     validated_mcl: set[str] | None = None) -> dict[str, str]:
    """Return the original overall layer plus a deliberately conservative PHB flag."""
    validated_phb = validated_phb or set()
    validated_mcl = validated_mcl or set()
    accession = row.get("accession", "")
    assignment = row.get("assignment", "").strip().lower()
    if accession in validated_mcl or "mcl" in assignment:
        focus, reason = "MCL_or_nonPHB_review", "MCL_lipase_associated_or_MCL_assignment"
        subtype = "iPhaZ_MCL_competition_evidence" if row.get("family") == "iPhaZ" else "unresolved_competition_evidence"
    elif accession in validated_phb:
        focus, reason = "PHB_candidate_supported", "independent_PHB_panel_match"
        subtype = "PHB_reference_supported_evidence"
    elif row.get("family") in {"ePhaZ", "iPhaZ"}:
        focus, reason = "PHB_candidate_review", "no_subtype_evidence"
        if row.get("family") == "ePhaZ" and "nonsecreted" in row.get("layer", ""):
            subtype = "ePhaZ_core_nonsecreted_evidence"
        elif row.get("family") == "ePhaZ" and "secreted" in row.get("layer", ""):
            subtype = "ePhaZ_SCL_secreted_evidence"
        elif row.get("family") == "ePhaZ" and "tier2_review" in row.get("layer", ""):
            subtype = "ePhaZ_tier2_review_evidence"
        elif row.get("family") == "ePhaZ" and "competition_review" in row.get("layer", ""):
            subtype = "ePhaZ_competition_review_evidence"
        elif row.get("family") == "iPhaZ":
            subtype = ("iPhaZ_tier2_review_evidence" if "tier2_review" in row.get("layer", "")
                       else "iPhaZ_intracellular_evidence")
        else:
            subtype = "unresolved_pha_evidence"
    else:
        focus, reason = "not_PHB_depolymerase", "background_family"
        subtype = "background_family_evidence"
    return {
        "overall_pha_layer": row.get("layer", ""),
        "pha_subtype_evidence": subtype,
        "phb_focus": focus,
        "phb_exclusion_reason": reason,
    }


def audit_archaeal_row(row: dict[str, str]) -> dict[str, str]:
    """Add explicit audit status without converting candidates into phenotype calls."""
    layer = row.get("layer", "")
    domain = row.get("domain", "").strip().lower()
    conflict = _bool(row.get("patatin_conflict", False))
    try:
        coverage = float(row.get("coverage", 0.0))
        evalue = float(row.get("evalue", "inf"))
    except (TypeError, ValueError):
        coverage, evalue = 0.0, float("inf")
    markers = {x.strip() for x in row.get("nearby_markers", "").replace(",", ";").split(";") if x.strip()}
    if not layer:
        if domain not in {"archaea", "d__archaea"} or conflict:
            layer = "non_PhaZh1_patatin"
        elif evalue > 1e-5 or coverage < 0.6:
            layer = "archaeal_patatin_exploratory"
        elif coverage >= 0.8 and _bool(row.get("complete", False)) and markers & {"BdhA", "PhaJ"}:
            layer = "PhaZh1_like_high"
        elif _bool(row.get("complete", False)) and markers & {"BdhA", "PhaJ", "PhaC", "PhaE", "phasin"}:
            layer = "PhaZh1_like_review"
        else:
            layer = "archaeal_patatin_exploratory"
    if domain not in {"archaea", "d__archaea"} or conflict:
        status = "excluded"
    elif layer == "PhaZh1_like_high":
        status = "review_required"
    elif layer == "PhaZh1_like_review":
        status = "review_required"
    else:
        status = "exploratory"
    return {
        "layer": layer,
        "audit_status": status,
        "phenotype_claim": "candidate_only",
    }


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta_accessions(path: Path) -> set[str]:
    """Read FASTA identifiers without loading sequences into the audit table."""
    return {line[1:].split(None, 1)[0] for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith(">")}


def n_terminal_gxsxg(path: Path, limit: int = 80) -> dict[str, str]:
    """Record an N-terminal patatin-associated motif as audit evidence only."""
    sequences: dict[str, list[str]] = {}
    accession = ""
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            accession = line[1:].split(None, 1)[0]
            sequences[accession] = []
        elif accession:
            sequences[accession].append(line.strip())
    return {key: ("present" if re.search(r"G.S.G", "".join(value)[:limit], re.I) else "absent")
            for key, value in sequences.items()}


def _accessions(path: Path, classification: str) -> set[str]:
    if not path or not path.exists():
        return set()
    accepted = {classification}
    if classification == "MCL_lipase_associated":
        # The frozen candidate calibration predates the reporting label and
        # records the same competing signal as `MCL_like`.
        accepted.add("MCL_like")
    return {r["accession"] for r in _read(path) if r.get("classification") in accepted}


def run(layer_tsv: Path, arch_tsv: Path, taxonomy_tsv: Path, context_tsv: Path,
        output_dir: Path, dual_profile_tsv: Path | None = None, arch_fasta: Path | None = None) -> dict[str, object]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to write into non-empty output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    phb = _accessions(dual_profile_tsv, "PHB_like") if dual_profile_tsv else set()
    mcl = _accessions(dual_profile_tsv, "MCL_lipase_associated") if dual_profile_tsv else set()
    pha_rows = _read(layer_tsv)
    out_rows = []
    for row in pha_rows:
        result = dict(row)
        result.update(classify_pha_row(row, phb, mcl))
        out_rows.append(result)
    fields = list(out_rows[0]) if out_rows else ["accession", "family", "layer", "overall_pha_layer", "phb_focus", "phb_exclusion_reason"]
    with (output_dir / "overall_pha_phb.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(out_rows)
    counts = Counter(r["phb_focus"] for r in out_rows)

    taxonomy = {r.get("genome", ""): r for r in _read(taxonomy_tsv)}
    context = {r.get("hit_locus", ""): r for r in _read(context_tsv)}
    sequence_ids = fasta_accessions(arch_fasta) if arch_fasta else set()
    nterm_motifs = n_terminal_gxsxg(arch_fasta) if arch_fasta else {}
    arch_rows = []
    for row in _read(arch_tsv):
        merged = dict(row)
        tax = taxonomy.get(row.get("genome", ""), {})
        ctx = context.get(row.get("locus", ""), {})
        target = f"{row.get('genome', '')}|{row.get('locus', '')}"
        merged.update({"phylum": tax.get("phylum", "pending"), "class": tax.get("class", "pending"),
                       "taxonomy": tax.get("gtdb_taxonomy", "pending"),
                       "context_markers": ctx.get("nearby_markers", row.get("nearby_markers", "")),
                       "sequence_present": "true" if target in sequence_ids else "false" if arch_fasta else "pending",
                       "nterm_gxsxg_motif": nterm_motifs.get(target, "pending"),
                       "domain_annotation": "candidate_PhaZh1_like_HMM;independent_domain_pending"})
        merged.update(audit_archaeal_row(merged))
        arch_rows.append(merged)
    afields = list(arch_rows[0]) if arch_rows else ["genome", "locus", "layer", "audit_status", "phenotype_claim"]
    with (output_dir / "archaeal_phazh1_full_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=afields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(arch_rows)
    audit_counts = Counter(r["audit_status"] for r in arch_rows)
    manifest = {"schema_version": 1, "status": "completed_candidate_only",
                "purpose": "overall PHA layering followed by conservative PHB focus and full archaeal PhaZh1-like audit",
                "counts": {"pha_records": len(out_rows), "archaeal_records": len(arch_rows),
                           "phb_focus": dict(sorted(counts.items())), "archaeal_audit": dict(sorted(audit_counts.items()))},
                "inputs": {str(p.name): {"path": str(p.resolve()), "bytes": p.stat().st_size, "sha256": _sha(p)}
                           for p in (layer_tsv, arch_tsv, taxonomy_tsv, context_tsv, arch_fasta) if p and p.exists()},
                "dual_profile_input": None if dual_profile_tsv is None else {"path": str(dual_profile_tsv.resolve()), "bytes": dual_profile_tsv.stat().st_size, "sha256": _sha(dual_profile_tsv)},
                "formal_registry_modified": False, "formal_scan_started": False,
                "phenotype_boundary": "All homology, domain, signal, neighborhood, and tree evidence remains candidate-only; archaeal hits are not included in the bacterial PHB main conclusion."}
    (output_dir / "analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--layer-tsv", type=Path, required=True)
    parser.add_argument("--arch-tsv", type=Path, required=True)
    parser.add_argument("--taxonomy-tsv", type=Path, required=True)
    parser.add_argument("--context-tsv", type=Path, required=True)
    parser.add_argument("--dual-profile-tsv", type=Path)
    parser.add_argument("--arch-fasta", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run(args.layer_tsv, args.arch_tsv, args.taxonomy_tsv, args.context_tsv, args.output_dir, args.dual_profile_tsv, args.arch_fasta)


if __name__ == "__main__":
    main()
