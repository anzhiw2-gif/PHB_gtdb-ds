#!/usr/bin/env python3
"""Build a dated, fail-closed amendment for the PhaDED reference panel.

This stage reconciles accession records already present in the immutable
reference ledger with NCBI Protein snapshots.  It never creates a family call
and never treats annotation-only records as formal negatives.

Rule #3 (axis-anchored credit): ``by_profile`` below only sees ledger records whose
family is one of the reference-only families, so an operator may additionally pass the
discrimination-axis panel (``--axis-panel`` + ``--axis-definitions``).  Only then are
the extra ``axis_anchored_*`` columns appended; the legacy counters are untouched, and
without those two inputs the readiness output is byte-for-byte what it was before.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import importlib.util
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


AXIS_PANEL_REQUIRED = {
    "axis_id", "axis_end", "role", "accession", "evidence_status", "anchor_evidence_kind",
    "family_binding_status", "bound_family", "bound_reference_only_profile_id",
    "anchor_justification", "source_doi_or_pmid",
}
DEFAULT_AXIS_CREDIT_SCRIPT = "credit_phaded_axis_anchored_controls.py"


def _load_credit_module(script_name: str = DEFAULT_AXIS_CREDIT_SCRIPT):
    """Load the axis-credit module from this script's own directory."""
    path = Path(__file__).resolve().parent / script_name
    spec = importlib.util.spec_from_file_location(Path(script_name).stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _base(accession: str) -> str:
    return accession.strip().split(".", 1)[0]


def infer_run_id(output_dir: str | Path) -> str:
    """Return the dated run directory containing ``results``."""
    path = Path(output_dir).resolve()
    if path.name == "results":
        return path.parent.name
    if path.parent.name == "results":
        return path.parent.parent.name
    raise ValueError(f"output directory must be a run results directory: {output_dir}")


def parse_genpept(lines: list[str]) -> list[dict[str, str]]:
    """Parse the small metadata subset needed from GenPept flat files."""
    records: list[dict[str, str]] = []
    block: list[str] = []
    for line in lines + ["//"]:
        if line.strip() == "//":
            if block:
                joined = "\n".join(block)
                acc = re.search(r"^ACCESSION\s+([^\s]+)", joined, re.MULTILINE)
                version = re.search(r"^VERSION\s+([^\s]+)", joined, re.MULTILINE)
                definition = re.search(
                    r"^DEFINITION\s+(.+?)(?=\n(?:ACCESSION|VERSION|DBLINK|KEYWORDS|SOURCE)\b)",
                    joined,
                    re.MULTILINE | re.DOTALL,
                )
                if acc:
                    pmids = re.findall(r"^\s*PUBMED\s+(\d+)", joined, re.MULTILINE)
                    records.append(
                        {
                            "accession": (version.group(1) if version else acc.group(1)).strip(),
                            "definition": re.sub(r"\s+", " ", (definition.group(1) if definition else "")).strip(),
                            "pmids": ",".join(dict.fromkeys(pmids)),
                        }
                    )
            block = []
        else:
            block.append(line.rstrip("\r\n"))
    return records


def load_ledger(path: str | Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    exact: dict[str, str] = {}
    bases: dict[str, list[str]] = defaultdict(list)
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row.get("accession", "").strip()
            family = row.get("phaded_family_id", "").strip()
            if not accession or not family:
                continue
            exact.setdefault(accession, family)
            if family not in bases[_base(accession)]:
                bases[_base(accession)].append(family)
    return exact, dict(bases)


def bind_accession(accession: str, exact: dict[str, str], bases: dict[str, list[str]]) -> tuple[str, str]:
    accession = accession.strip()
    if accession in exact:
        return "resolved_existing_family", exact[accession]
    families = bases.get(_base(accession), [])
    if len(families) == 1:
        return "resolved_existing_family_versionless", families[0]
    if len(families) > 1:
        return "ambiguous_existing_family", ""
    return "unresolved", ""


def role_for(evidence_status: str, binding_scope: str) -> str:
    if evidence_status == "experimental_positive":
        return "independent_positive"
    if evidence_status == "experimental_negative":
        return "family_resolved_negative"
    return "challenge_control"


def eligibility_for(evidence_status: str, binding_scope: str) -> str:
    if evidence_status == "annotation_only":
        return "not_eligible_annotation_only"
    if evidence_status == "experimental_negative" and binding_scope.startswith("resolved"):
        return "eligible_negative"
    if evidence_status == "experimental_positive" and binding_scope.startswith("resolved"):
        return "eligible_positive_pending_independence"
    return "not_eligible_unresolved_family"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_amendment(profile_path: str | Path, ledger_path: str | Path, gp_paths: list[str | Path], output_dir: str | Path, *, axis_panel: str | Path | None = None, axis_definitions: str | Path | None = None) -> dict[str, object]:
    profile_path, ledger_path, output_dir = Path(profile_path), Path(ledger_path), Path(output_dir)
    if (axis_panel is None) != (axis_definitions is None):
        raise ValueError("axis_anchored credit requires both axis_panel and axis_definitions")
    profiles = _read_tsv(profile_path)
    reference_only = [row for row in profiles if row.get("model_status", "").strip() == "reference_only"]
    if len(reference_only) != 37:
        raise ValueError(f"expected 37 reference-only profiles, observed {len(reference_only)}")
    target_by_family = {row.get("phaded_family_id", "").strip(): row for row in reference_only if row.get("profile_kind") == "family" and row.get("phaded_family_id", "").strip()}
    exact, bases = load_ledger(ledger_path)
    ledger_rows = _read_tsv(ledger_path)
    gp: dict[str, dict[str, str]] = {}
    for path_text in gp_paths:
        path = Path(path_text)
        gp_records = parse_genpept(path.read_text(encoding="utf-8", errors="replace").splitlines())
        for record in gp_records:
            gp.setdefault(record["accession"], record)
            gp.setdefault(_base(record["accession"]), record)

    training = {row.get("phaded_family_id", "").strip(): {x.strip().split(".", 1)[0] for x in row.get("training_accessions", "").split(";") if x.strip()} for row in profiles}
    records: list[dict[str, str]] = []
    for source in ledger_rows:
        family = source.get("phaded_family_id", "").strip()
        if family not in target_by_family:
            continue
        accession = source.get("accession", "").strip()
        binding, bound_family = bind_accession(accession, exact, bases)
        metadata = gp.get(accession) or gp.get(_base(accession), {})
        status = source.get("evidence_status", "").strip() or "pending_review"
        role = role_for(status, binding)
        eligibility = eligibility_for(status, binding)
        if status == "experimental_positive" and _base(accession) in training.get(family, set()):
            eligibility = "not_eligible_training_accession"
        if status == "experimental_positive" and source.get("pmid", "").strip() == "19296857":
            eligibility = "not_eligible_same_publication_lineage"
        records.append({
            "accession": accession,
            "profile_id": target_by_family[family].get("profile_id", "").strip(),
            "phaded_family_id": family,
            "phaded_superfamily": source.get("phaded_superfamily", "").strip(),
            "family_binding_status": binding,
            "role": role,
            "evidence_status": status,
            "eligibility": eligibility,
            "definition": metadata.get("definition", ""),
            "pmid": source.get("pmid", "").strip(),
            "doi": source.get("primary_doi", "").strip(),
            "sequence_sha256": source.get("sequence_sha256", "").strip(),
            "source_database": source.get("source_database", "").strip(),
            "new_family_call_created": "false",
        })

    output_dir.mkdir(parents=True, exist_ok=True)
    panel_path = output_dir / "panel_evidence_amendment.tsv"
    readiness_path = output_dir / "profile_calibration_readiness.tsv"
    report_path = output_dir / "panel_acquisition_amendment_report.json"
    fields = list(records[0]) if records else ["accession"]
    with panel_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)

    readiness: list[dict[str, str]] = []
    by_profile: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in records:
        by_profile[row["profile_id"]].append(row)
    for profile in sorted(reference_only, key=lambda row: row.get("profile_id", "")):
        pid = profile.get("profile_id", "").strip()
        rows = by_profile.get(pid, [])
        positive = sum(row["eligibility"] == "eligible_positive_pending_independence" for row in rows)
        negative = sum(row["eligibility"] == "eligible_negative" for row in rows)
        challenge = sum(row["role"] == "challenge_control" and row["family_binding_status"].startswith("resolved") for row in rows)
        existing_positive = sum(row.get("evidence_status") == "experimental_positive" for row in rows)
        readiness.append({
            "profile_id": pid,
            "profile_kind": profile.get("profile_kind", "").strip(),
            "phaded_superfamily": profile.get("phaded_superfamily", "").strip(),
            "phaded_family_id": profile.get("phaded_family_id", "").strip(),
            "existing_experimental_positive": str(existing_positive),
            "external_bound_positive": str(positive),
            "external_bound_negative": str(negative),
            "external_bound_challenge": str(challenge),
            "heldout_positive": "0",
            "family_resolved_negative_status": "sufficient" if negative else "missing",
            "challenge_status": "sufficient" if challenge else "missing",
            "calibration_decision": "reference_only_insufficient_panel" if profile.get("profile_kind") == "family" else "planned_not_run_superfamily_requires_family_resolution",
            "new_family_call_created": "false",
        })
    axis_credit = None
    if axis_panel is not None:
        credit_module = _load_credit_module()
        axis_rows = _read_tsv(Path(axis_panel))
        if not axis_rows:
            raise ValueError(f"axis panel is empty: {axis_panel}")
        missing_axis_columns = sorted(AXIS_PANEL_REQUIRED - set(axis_rows[0]))
        if missing_axis_columns:
            raise ValueError("axis panel missing columns: " + ",".join(missing_axis_columns))
        definitions = credit_module.load_axis_definitions(axis_definitions)
        family_superfamily: dict[str, str] = {}
        for row in profiles:
            family = row.get("phaded_family_id", "").strip()
            superfamily = row.get("phaded_superfamily", "").strip()
            if not family or not superfamily:
                continue
            previous = family_superfamily.get(family)
            if previous is not None and previous != superfamily:
                raise ValueError(f"{family}: conflicting superfamily {previous!r} vs {superfamily!r}")
            family_superfamily[family] = superfamily
        axis_credit = credit_module.attribute_axis_credit(axis_rows, definitions, reference_only, family_superfamily)
        for row in readiness:
            payload = axis_credit["profiles"].get(row["profile_id"], {})
            row["axis_anchored_negative_count"] = str(payload.get("axis_anchored_negative_count", 0))
            row["axis_anchored_challenge_count"] = str(payload.get("axis_anchored_challenge_count", 0))
            row["axis_anchor_basis_status"] = payload.get("axis_anchor_basis_status", "not_audited")
            row["axis_anchored_negative_credit_ids"] = ";".join(payload.get("axis_anchored_negative_credit_ids", []))
            row["axis_anchored_challenge_credit_ids"] = ";".join(payload.get("axis_anchored_challenge_credit_ids", []))
    readiness_fields = list(readiness[0])
    with readiness_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=readiness_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(readiness)
    report = {
        "schema_version": 1,
        "run_id": infer_run_id(output_dir),
        "status": "completed_candidate_only",
        "inputs": {"profiles": {"path": str(profile_path.resolve()), "sha256": _sha256(profile_path)}, "ledger": {"path": str(ledger_path.resolve()), "sha256": _sha256(ledger_path)}, "genpept_files": [{"path": str(Path(p).resolve()), "sha256": _sha256(Path(p))} for p in gp_paths]},
        "summary": {
            "reference_only_profiles": len(reference_only),
            "panel_records": len(records),
            "resolved_records": sum(row["family_binding_status"].startswith("resolved") for row in records),
            "formal_negative_records": sum(row["eligibility"] == "eligible_negative" for row in records),
            "independent_positive_records": sum(row["eligibility"] == "eligible_positive_pending_independence" for row in records),
            "new_family_call_created": False,
        },
        "outputs": {"panel": {"path": str(panel_path.resolve()), "sha256": _sha256(panel_path)}, "readiness": {"path": str(readiness_path.resolve()), "sha256": _sha256(readiness_path)}},
        "boundaries": ["Only accession-bound existing families are retained.", "Annotation-only records are never formal negatives.", "Same-publication-lineage positives are not independent held-out positives.", "No new family call, profile training, or phenotype claim is created."],
    }
    if axis_credit is not None:
        report["axis_anchored_credit"] = {
            "convention": "axis_anchored_credit_v1",
            "legacy_counters_unchanged": True,
            "source_files": {
                "axis_panel": {"path": str(Path(axis_panel).resolve()), "sha256": _sha256(Path(axis_panel))},
                "axis_definitions": {"path": str(Path(axis_definitions).resolve()), "sha256": _sha256(Path(axis_definitions))},
            },
            "counts": axis_credit["counts"],
            "boundaries": [
                "Axis-anchored credit is an additional, separately named counting channel; it never feeds the family-bound counters.",
                "A record only counts when axis, anchor_justification and source DOI/PMID are all present and its evidence kind is pre-registered for that axis end.",
                "A record only credits a profile that holds a family-resolved anchor on the opposite end of the same axis.",
            ],
        }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--genpept", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--axis-panel", type=Path, default=None, help="optional discrimination-axis panel for axis_anchored credit")
    parser.add_argument("--axis-definitions", type=Path, default=None, help="optional axis_reachability.json carrying the frozen axis definitions")
    args = parser.parse_args()
    print(json.dumps(build_amendment(args.profiles, args.ledger, args.genpept, args.output_dir, axis_panel=args.axis_panel, axis_definitions=args.axis_definitions), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
