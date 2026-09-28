#!/usr/bin/env python3
"""Summarize InterProScan secondary evidence without changing PhaDED assignments."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EXPLICIT_SIGNATURES = {"PF10503", "PF06850", "NF050057", "NF050058", "NF050059", "TIGR01840", "TIGR01849"}
EXPLICIT_INTERPRO = {"IPR010126"}
GENERIC_SIGNATURES = {"PF00756", "PF00561", "PF12697", "PF02230"}
SYNTHASE_SIGNATURES = {"PTHR36837", "IPR051321"}


def classify_hits(hits: list[dict[str, str]]) -> str:
    if not hits:
        return "no_interpro_hit"
    if any(hit.get("signature", "") in SYNTHASE_SIGNATURES or hit.get("interpro_id", "") in SYNTHASE_SIGNATURES or "synthase" in (hit.get("description", "") + " " + hit.get("interpro_desc", "")).lower() for hit in hits):
        return "pha_synthase_like_alternative"
    if any(hit.get("signature", "") in EXPLICIT_SIGNATURES or hit.get("interpro_id", "") in EXPLICIT_INTERPRO for hit in hits):
        return "explicit_pha_related_domain"
    if any(hit.get("signature", "") in GENERIC_SIGNATURES or "esterase" in hit.get("description", "").lower() or "hydrolase" in hit.get("description", "").lower() for hit in hits):
        return "generic_hydrolase_or_esterase"
    return "other_interpro_domain"


def summarize(candidates: list[dict[str, str]], hits: list[dict[str, str]], excluded: dict[str, str] | None = None) -> list[dict[str, str]]:
    by_accession: dict[str, list[dict[str, str]]] = {}
    for hit in hits:
        by_accession.setdefault(hit.get("accession", ""), []).append(hit)
    excluded = excluded or {}
    output = []
    for candidate in candidates:
        accession = candidate.get("accession", "")
        state = "not_submitted_tool_input_limit" if accession in excluded else classify_hits(by_accession.get(accession, []))
        row = dict(candidate)
        row["interpro_secondary_state"] = state
        row["interpro_hit_count"] = str(len(by_accession.get(accession, [])))
        row["interpro_explicit_ids"] = ";".join(sorted({hit.get("interpro_id", "") for hit in by_accession.get(accession, []) if hit.get("interpro_id", "")}))
        row["interpro_signatures"] = ";".join(sorted({hit.get("signature", "") for hit in by_accession.get(accession, []) if hit.get("signature", "")}))
        output.append(row)
    return output


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    with args.candidates.open(encoding="utf-8-sig", newline="") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    with args.interpro.open(encoding="utf-8-sig", newline="") as handle:
        hits = list(csv.DictReader(handle, fieldnames=["accession", "md5", "length", "analysis", "signature", "description", "start", "end", "evalue", "status", "date", "interpro_id", "interpro_desc", "go", "pathways"], delimiter="\t"))
    with args.excluded.open(encoding="utf-8-sig", newline="") as handle:
        excluded = {row["accession"]: row["reason"] for row in csv.DictReader(handle, delimiter="\t")}
    rows = summarize(candidates, hits, excluded)
    fields = sorted({key for row in rows for key in row})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    manifest = {"status": "completed_candidate_only", "candidate_count": len(rows), "inputs": {"candidates": {"path": str(args.candidates.resolve()), "sha256": sha256(args.candidates)}, "interpro": {"path": str(args.interpro.resolve()), "sha256": sha256(args.interpro)}, "excluded": {"path": str(args.excluded.resolve()), "sha256": sha256(args.excluded)}}, "output": {"path": str(args.output.resolve()), "sha256": sha256(args.output)}, "phenotype_boundary": "InterProScan states are secondary candidate architecture evidence and do not prove PHB/PHA degradation."}
    (args.output.parent / "interpro_summary_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "candidate_count": len(rows)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
