#!/usr/bin/env python3
"""Aggregate validated InterProScan rows into one accession-level table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping


EXPLICIT_SIGNATURES = {"PF10503", "PF06850", "NF050057", "NF050058", "NF050059", "TIGR01840", "TIGR01849"}
EXPLICIT_INTERPRO = {"IPR010126"}
GENERIC_SIGNATURES = {"PF00756", "PF00561", "PF12697", "PF02230"}
SYNTHASE_SIGNATURES = {"PTHR36837", "IPR051321"}
HITS_FIELDS = (
    "accession", "md5", "length", "analysis", "signature", "description",
    "start", "end", "evalue", "status", "date", "interpro_id",
    "interpro_desc", "go", "pathways",
)
OUTPUT_FIELDS = (
    "accession", "interpro_status", "interpro_secondary_state", "interpro_hit_count",
    "interpro_analyses", "interpro_signatures", "interpro_ids", "interpro_dates",
)


def _classify_hits(hits: list[Mapping[str, str]]) -> str:
    if not hits:
        return "no_interpro_hit"
    if any(
        hit.get("signature", "") in SYNTHASE_SIGNATURES
        or hit.get("interpro_id", "") in SYNTHASE_SIGNATURES
        or "synthase" in (hit.get("description", "") + " " + hit.get("interpro_desc", "")).lower()
        for hit in hits
    ):
        return "pha_synthase_like_alternative"
    if any(hit.get("signature", "") in EXPLICIT_SIGNATURES or hit.get("interpro_id", "") in EXPLICIT_INTERPRO for hit in hits):
        return "explicit_pha_related_domain"
    if any(
        hit.get("signature", "") in GENERIC_SIGNATURES
        or "esterase" in hit.get("description", "").lower()
        or "hydrolase" in hit.get("description", "").lower()
        for hit in hits
    ):
        return "generic_hydrolase_or_esterase"
    return "other_interpro_domain"


def _unique_values(hits: Iterable[Mapping[str, str]], key: str) -> str:
    return ";".join(sorted({str(hit.get(key, "")).strip() for hit in hits if str(hit.get(key, "")).strip() and hit.get(key, "") != "-"}))


def summarize_interpro(
    candidates: list[Mapping[str, str]],
    hits: Iterable[Mapping[str, str]],
    excluded: Mapping[str, str],
) -> list[dict[str, str]]:
    accessions = [str(row.get("accession", "")).strip() for row in candidates]
    if any(not accession for accession in accessions) or len(accessions) != len(set(accessions)):
        raise ValueError("candidate accession list contains duplicate or empty accession")
    candidate_set = set(accessions)
    by_accession: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for hit in hits:
        accession = str(hit.get("accession", "")).strip()
        if not accession:
            raise ValueError("InterPro hit has empty accession")
        if accession not in candidate_set:
            raise ValueError(f"InterPro hit is outside candidate set: {accession}")
        by_accession[accession].append(hit)
    if not set(excluded).issubset(candidate_set):
        raise ValueError("tool exclusion contains accession outside candidate set")

    output: list[dict[str, str]] = []
    for accession in accessions:
        accession_hits = by_accession.get(accession, [])
        if accession in excluded:
            status = "interpro_tool_input_excluded"
        elif accession in by_accession:
            status = "interpro_supported"
        else:
            status = "interpro_not_reported"
        output.append({
            "accession": accession,
            "interpro_status": status,
            "interpro_secondary_state": _classify_hits(accession_hits),
            "interpro_hit_count": str(len(accession_hits)),
            "interpro_analyses": _unique_values(accession_hits, "analysis"),
            "interpro_signatures": _unique_values(accession_hits, "signature"),
            "interpro_ids": _unique_values(accession_hits, "interpro_id"),
            "interpro_dates": _unique_values(accession_hits, "date"),
        })
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


def write_summary(candidates_path: Path, interpro_path: Path, excluded_path: Path, output_path: Path) -> dict[str, object]:
    candidates = read_tsv(candidates_path)
    hits = read_tsv(interpro_path, list(HITS_FIELDS))
    excluded_rows = read_tsv(excluded_path)
    excluded = {row.get("accession", "").strip(): row.get("reason", "").strip() for row in excluded_rows if row.get("accession", "").strip()}
    rows = summarize_interpro(candidates, hits, excluded)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_FIELDS), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "status": "completed_candidate_only",
        "candidate_count": len(rows),
        "status_counts": {},
        "inputs": {
            "candidates": {"path": str(candidates_path.resolve()), "size": candidates_path.stat().st_size, "sha256": sha256(candidates_path)},
            "interpro": {"path": str(interpro_path.resolve()), "size": interpro_path.stat().st_size, "sha256": sha256(interpro_path)},
            "excluded": {"path": str(excluded_path.resolve()), "size": excluded_path.stat().st_size, "sha256": sha256(excluded_path)},
        },
        "output": {"path": str(output_path.resolve()), "size": output_path.stat().st_size, "sha256": sha256(output_path)},
        "phenotype_boundary": "InterProScan states are sequence/domain evidence only and do not prove PHB/PHA degradation phenotype.",
    }
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["interpro_status"]] = counts.get(row["interpro_status"], 0) + 1
    manifest["status_counts"] = counts
    (output_path.parent / "interpro_summary_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--interpro", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(write_summary(args.candidates, args.interpro, args.excluded, args.output), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
