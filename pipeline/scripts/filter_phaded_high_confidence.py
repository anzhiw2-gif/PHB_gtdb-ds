#!/usr/bin/env python3
"""Filter the PhaDED candidate pool into a per-superfamily high-confidence set.

Applies the preregistered multi-evidence criteria (filter_criteria.md): common
criteria (InterPro supported, no architecture conflict, strong profile, unique
superfamily, hydrophobic x1) plus superfamily-specific architecture criteria
(localization, lipase box, catalytic-domain geometry, lid, PF06850, AHSMG).

Fail-closed: missing/pending evidence never passes a criterion.

DEPRECATED / FROZEN v1 SEMANTICS -- COMPARISON INPUT ONLY
---------------------------------------------------------
This module is the v1 filter.  Its behaviour is frozen: the 36,611 / 36,559 /
4,131 / 40,690 counts and every column it emits are historical evidence and
must not be changed, re-run or rewritten.  Keep using this module only to
reproduce the v1 comparison set.

Its semantics are superseded.  In particular ``NT_AHSMG = "ahsmg"`` below makes
AHSMG a *third* nucleophile type alongside Ser and Cys, which is exactly the
error v2 corrects: AHSMG is a Ser-containing motif class, not a nucleophile
identity.  v2 (``build_phaded_candidate_catalog_v2.py``) therefore classifies
with ``nucleophile_identity`` in {ser, cys, unresolved} crossed with
``motif_class`` (GxSxG / AHSMG / Cys-associated / other / unresolved), with
AHSMG recorded as a Ser motif.  v2 also keeps SignalP transport predictions
separate from localization evidence and keeps SBD out of the type 1 / type 2
definition, neither of which this frozen v1 module does.

Do not add new semantics here.  Extend ``build_phaded_candidate_catalog_v2.py``
instead.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

EXPORT = {"SP", "LIPO", "TAT", "TATLIPO"}
NO_EXPORT = {"OTHER", "pending"}
HYDROPHOBIC = {"A", "C", "F", "I", "L", "M", "V", "W", "Y"}

# Honest-typing labels (§11.9.4.1): the Cys superfamily's catalytic residue is
# inferred (absence of GxSxG + PF06850 + strong Cys-family profile), NOT
# per-sequence verified; every other family is motif-level (residue identified
# at motif level, not biochemically verified per sequence).
CYS_SUPERFAMILY = "intracellular nPHASCL without lipase box"
CATALYTIC_RESIDUE_INFERRED_CYS = "inferred_cys_not_verified"
CATALYTIC_RESIDUE_MOTIF = "motif_level"

# Nucleophile typing (direct marker, not alignment transfer): Ser families carry
# a GxSxG lipase box (direct Ser-nucleophile marker); the Cys superfamily has no
# GxSxG and is typed by its family assignment; the PhaZ7-like superfamily uses
# the AHSMG motif instead.  A Ser family with a missing GxSxG, or the Cys family
# carrying a GxSxG, is a conflict (evidence contradicts the expected nucleophile).
#
# FROZEN v1 SEMANTICS -- DO NOT EXTEND.  ``NT_AHSMG`` is kept only so the v1
# comparison set reproduces byte-for-byte; treating AHSMG as a third nucleophile
# type is the defect v2 fixes.  In v2, AHSMG is a Ser motif class emitted through
# ``nucleophile_identity`` x ``motif_class`` (see the module docstring).
PHAZ7_SUPERFAMILY = "extracellular native-SCL/PhaZ7-like"
NT_SER = "ser"
NT_CYS = "cys"
NT_AHSMG = "ahsmg"
NT_CONFLICT = "conflict"
NT_UNDETERMINED = "undetermined"


def nucleophile_type(superfamily: str, lipase_box_state: str, ahsmg_state: str) -> str:
    """Family-level nucleophile typing for one candidate."""
    if superfamily == CYS_SUPERFAMILY:
        return NT_CONFLICT if lipase_box_state == "supported" else NT_CYS
    if superfamily == PHAZ7_SUPERFAMILY:
        return NT_AHSMG if ahsmg_state == "supported" else NT_CONFLICT
    if lipase_box_state == "supported":
        return NT_SER
    if lipase_box_state == "not_detected_pattern":
        return NT_CONFLICT
    return NT_UNDETERMINED


def nucleophile_family(superfamily: str, ntype: str, cl, mc) -> str:
    """Direct nucleophile marker: family HMM for Cys, GxSxG for Ser, AHSMG for PhaZ7."""
    if ntype == NT_CYS:
        fam = ""
        if cl is not None:
            fam = cl.get("families_in_best_superfamily", "") or ""
        if not fam:
            fam = mc.get("phaded_family_best", "") or ""
        return fam or "cys_family"
    if ntype == NT_SER:
        return "GxSxG"
    if ntype == NT_AHSMG:
        return "AHSMG"
    return ""

# superfamily -> (export_requirement, architecture_requirements)
# export_requirement: "export" | "no_export" | "any"
# architecture_requirements: list of field=value pairs evaluated against the
# merged per-candidate evidence dict.
SUPERFAMILY_CRITERIA = {
    "extracellular dPHASCL type 1": (
        "export",
        ["lipase_box_state=supported", "x1_hydrophobic=1", "catalytic_domain_type=type1_verified"],
    ),
    "extracellular dPHASCL type 2": (
        "export",
        ["lipase_box_state=supported", "x1_hydrophobic=1", "catalytic_domain_type=type2_verified"],
    ),
    "intracellular nPHASCL without lipase box": (
        "no_export",
        ["lipase_box_state=not_detected_pattern", "pf06850=detected"],
    ),
    "intracellular nPHASCL with lipase box": (
        "no_export",
        ["lipase_box_state=supported", "x1_hydrophobic=1"],
    ),
    "intracellular nPHAMCL": (
        "no_export",
        ["lipase_box_state=supported", "x1_hydrophobic=1", "lid_state=supported"],
    ),
    "extracellular dPHAMCL": (
        "export",
        ["lipase_box_state=supported", "x1_hydrophobic=1"],
    ),
    "periplasmic PHA depolymerases": (
        "export",
        ["catalytic_domain_type=type2_verified"],
    ),
    "extracellular native-SCL/PhaZ7-like": (
        "export",
        ["ahsmg_state=supported"],
    ),
}


def load_tsv(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("accession") or row.get("protein_id") or ""
            if key:
                out[key] = row
    return out


def load_confounders(path: Path) -> set[str]:
    out: set[str] = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for field in ("accession", "protein_id"):
                if row.get(field):
                    out.add(row[field])
                    break
    return out


def evaluate(criteria: str, evidence: dict[str, str]) -> bool:
    field, expected = criteria.split("=", 1)
    if field == "x1_hydrophobic":
        value = evidence.get("lipase_box_x1", "")
        return value in HYDROPHOBIC
    if field == "pf06850":
        return evidence.get("sbd_pf06850_binding_state", "") == "detected"
    return evidence.get(field, "") == expected


def filter_candidates(
    classification: Path,
    motif: Path,
    domain_type: Path,
    matrix: Path,
    confounders: Path,
    out_dir: Path,
) -> dict[str, object]:
    cls = load_tsv(classification)
    motif_rows = load_tsv(motif)
    type_rows = load_tsv(domain_type)
    conf = load_confounders(confounders)

    # matrix: signalp_class, architecture_consistency, interpro_status,
    # profile_evidence_status, phaded_superfamily_best
    matrix_cols = {}
    with matrix.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            acc = row.get("accession", "")
            if acc:
                matrix_cols[acc] = row

    out_dir.mkdir(parents=True, exist_ok=True)

    stats: dict[str, dict[str, int]] = defaultdict(
        lambda: {"total": 0, "pass": 0, "fail_localization": 0, "fail_architecture": 0, "fail_common": 0}
    )
    rows_out: list[dict[str, str]] = []
    hold_rows: list[dict[str, str]] = []

    for acc, mrow in motif_rows.items():
        mc = matrix_cols.get(acc, {})
        cl = cls.get(acc)

        # Superfamily assignment precedence: the trained-profile-derived matrix
        # superfamily FIRST (a hit at family OR superfamily level counts as a
        # trained-profile hit), then the discovery-layer classification for the
        # previously-unassigned candidates (unique superfamily only).
        sf = None
        strong_profile = False
        best_e_str = ""
        profile_evidence = ""
        matrix_sf = mc.get("phaded_superfamily_best", "")
        matrix_pev = mc.get("profile_evidence_status", "")
        if matrix_sf in SUPERFAMILY_CRITERIA and matrix_pev in ("profile_trained_hit", "profile_ambiguous_family"):
            sf = matrix_sf
            strong_profile = True
            profile_evidence = matrix_pev
        elif cl is not None and cl.get("superfamily_confidence") == "unique" and cl.get("superfamily_claim", "").strip() in SUPERFAMILY_CRITERIA:
            sf = cl.get("superfamily_claim", "").strip()
            best_e_str = cl.get("best_evalue", "")
            try:
                best_e = float(best_e_str or "1")
            except ValueError:
                best_e = 1.0
            strong_profile = best_e < 1e-10
            profile_evidence = "discovery_hit"

        if sf is None:
            continue

        stats[sf]["total"] += 1

        sig = mc.get("signalp_class", "")
        trow = type_rows.get(acc)
        ctype = trow.get("catalytic_domain_type", "") if trow else ""
        lb = mrow.get("lipase_box_state", "")
        ahsmg = mrow.get("ahsmg_state", "")
        ntype = nucleophile_type(sf, lb, ahsmg)
        nfam = nucleophile_family(sf, ntype, cl, mc)

        def hold_row(reason: str) -> dict[str, str]:
            return {
                "accession": acc,
                "genome": mc.get("genome", ""),
                "superfamily": sf,
                "nucleophile_type": ntype,
                "hold_reason": reason,
                "lipase_box_state": lb,
                "sbd_pf06850_binding_state": mrow.get("sbd_pf06850_binding_state", ""),
                "signalp_class": sig,
                "catalytic_domain_type": ctype,
                "interpro_status": mc.get("interpro_status", ""),
                "profile_evidence_status": profile_evidence,
            }

        # common criteria
        ok_common = (
            mc.get("interpro_status", "") == "interpro_supported"
            and mc.get("architecture_consistency", "") != "conflicting"
            and strong_profile
            and acc not in conf
        )
        if not ok_common:
            stats[sf]["fail_common"] += 1
            hold_rows.append(hold_row("common_criteria"))
            continue

        export_req, arch_criteria = SUPERFAMILY_CRITERIA[sf]
        if export_req == "export" and sig not in EXPORT:
            stats[sf]["fail_localization"] += 1
            hold_rows.append(hold_row("localization_conflict"))
            continue
        if export_req == "no_export" and sig not in NO_EXPORT:
            stats[sf]["fail_localization"] += 1
            hold_rows.append(hold_row("localization_conflict"))
            continue

        evidence = dict(mrow)
        if trow:
            evidence["catalytic_domain_type"] = ctype

        all_ok = True
        for crit in arch_criteria:
            if not evaluate(crit, evidence):
                all_ok = False
                break
        if not all_ok:
            stats[sf]["fail_architecture"] += 1
            hold_rows.append(
                hold_row("nucleophile_conflict" if ntype == NT_CONFLICT else "architecture_criteria")
            )
            continue

        stats[sf]["pass"] += 1
        rows_out.append(
            {
                "accession": acc,
                "genome": mc.get("genome", ""),
                "superfamily": sf,
                "signalp_class": sig,
                "profile_best_evalue": best_e_str,
                "profile_evidence_status": profile_evidence,
                "lipase_box_state": lb,
                "lipase_box_x1": mrow.get("lipase_box_x1", ""),
                "catalytic_domain_type": ctype,
                "catalytic_residue_verification": (
                    CATALYTIC_RESIDUE_INFERRED_CYS if sf == CYS_SUPERFAMILY else CATALYTIC_RESIDUE_MOTIF
                ),
                "nucleophile_type": ntype,
                "nucleophile_family": nfam,
                "lid_state": mrow.get("lid_state", ""),
                "sbd_pf06850_binding_state": mrow.get("sbd_pf06850_binding_state", ""),
                "ahsmg_state": ahsmg,
                "interpro_status": mc.get("interpro_status", ""),
                "high_confidence": "1",
            }
        )

    fieldnames = [
        "accession", "genome", "superfamily", "signalp_class", "profile_best_evalue",
        "profile_evidence_status", "lipase_box_state", "lipase_box_x1",
        "catalytic_domain_type", "catalytic_residue_verification",
        "nucleophile_type", "nucleophile_family",
        "lid_state", "sbd_pf06850_binding_state",
        "ahsmg_state", "interpro_status", "high_confidence",
    ]
    out_path = out_dir / "high_confidence_candidates.tsv"
    with out_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows_out)

    hold_fieldnames = [
        "accession", "genome", "superfamily", "nucleophile_type", "hold_reason",
        "lipase_box_state", "sbd_pf06850_binding_state", "signalp_class",
        "catalytic_domain_type", "interpro_status", "profile_evidence_status",
    ]
    hold_path = out_dir / "hold_candidates.tsv"
    with hold_path.open("w", encoding="utf-8", newline="\n") as handle:
        w = csv.DictWriter(handle, fieldnames=hold_fieldnames, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(hold_rows, key=lambda r: (r["superfamily"], r["accession"])))

    hold_breakdown: dict[str, int] = defaultdict(int)
    for h in hold_rows:
        hold_breakdown[h["hold_reason"]] += 1

    summary = {
        "criteria": "see inputs/filter_criteria.md",
        "total_high_confidence": len(rows_out),
        "total_held": len(hold_rows),
        "hold_breakdown": dict(sorted(hold_breakdown.items())),
        "per_superfamily": {sf: dict(s) for sf, s in sorted(stats.items())},
        "honesty_notes": {
            "cys_typing": (
                "intracellular nPHASCL without lipase box rows are labelled "
                "catalytic_residue_verification=inferred_cys_not_verified: the "
                "Cys typing is inferred from absence of a GxSxG lipase box plus "
                "PF06850 plus a strong Cys-family profile; the catalytic residue "
                "is NOT per-sequence verified and an undetected Ser type without "
                "GxSxG is not excluded."
            ),
            "general": (
                "no catalytic residue is biochemically verified per sequence; "
                "non-Cys rows are catalytic_residue_verification=motif_level "
                "(residue identified at motif level)."
            ),
        },
    }
    (out_dir / "filter_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--motif", type=Path, required=True)
    parser.add_argument("--domain-type", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--confounders", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(
        filter_candidates(args.classification, args.motif, args.domain_type,
                          args.matrix, args.confounders, args.out_dir),
        indent=2, ensure_ascii=False,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
