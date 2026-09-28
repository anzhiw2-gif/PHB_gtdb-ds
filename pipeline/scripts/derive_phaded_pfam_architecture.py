#!/usr/bin/env python3
"""Derive and apply a targeted Pfam fingerprint for PhaDED architecture review."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_domtblout(lines: Iterable[str]) -> list[dict[str, object]]:
    """Parse HMMER domtblout, preserving query IDs and 1-based alignment coordinates."""
    hits: list[dict[str, object]] = []
    for raw in lines:
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 22:
            continue
        try:
            hits.append({
                "pfam_name": fields[0],
                "pfam_accession": fields[1].split(".", 1)[0],
                "accession": fields[3],
                "i_evalue": float(fields[12]),
                "bitscore": float(fields[13]),
                "coordinates": f"{int(fields[17])}-{int(fields[18])}",
                "hmm_coordinates": f"{int(fields[15])}-{int(fields[16])}",
                "description": " ".join(fields[22:]),
            })
        except (ValueError, IndexError):
            continue
    return hits


def classify_candidate_domains(
    best_superfamily: str,
    pfam_accessions: set[str],
    expected: Mapping[str, set[str]],
) -> dict[str, str]:
    if not best_superfamily or best_superfamily not in expected:
        return {"pfam_interpro_state": "pending_no_superfamily", "architecture_consistency": "pending", "assignment_review": "review_required"}
    target = expected[best_superfamily]
    if not target:
        return {"pfam_interpro_state": "pending_no_reference_fingerprint", "architecture_consistency": "pending", "assignment_review": "no_conflict_recorded"}
    if pfam_accessions & target:
        return {"pfam_interpro_state": "tested_targeted_pfam_hit", "architecture_consistency": "partial", "assignment_review": "no_conflict_recorded"}
    if pfam_accessions and any(pfam_accessions & other for sf, other in expected.items() if sf != best_superfamily):
        return {"pfam_interpro_state": "tested_other_superfamily_hit", "architecture_consistency": "conflicting", "assignment_review": "review_required"}
    return {"pfam_interpro_state": "tested_no_targeted_pfam_hit", "architecture_consistency": "pending", "assignment_review": "no_conflict_recorded"}


def reference_fingerprint(reference_hits: Iterable[Mapping[str, object]], reference_superfamilies: Mapping[str, str], *, min_reference_count: int = 2) -> dict[str, set[str]]:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for hit in reference_hits:
        accession = str(hit["accession"]).split("|", 1)[0]
        sf = reference_superfamilies.get(accession)
        if sf:
            counts[sf][str(hit["pfam_accession"])] += 1
    return {sf: {pfam for pfam, count in counter.items() if count >= min_reference_count} for sf, counter in counts.items()}


def write_fingerprint(path: Path, fingerprint: Mapping[str, set[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["phaded_superfamily", "pfam_accession"])
        for sf in sorted(fingerprint):
            for pfam in sorted(fingerprint[sf]):
                writer.writerow([sf, pfam])


def selected_name_rows(hits: Iterable[Mapping[str, object]], selected_accessions: set[str]) -> list[tuple[str, str]]:
    """Return one deterministic Pfam accession-to-NAME mapping for hmmfetch."""
    pairs = {
        (str(hit["pfam_accession"]), str(hit["pfam_name"]))
        for hit in hits
        if str(hit["pfam_accession"]) in selected_accessions
    }
    return sorted(pairs)


def apply(candidate_hits: Iterable[Mapping[str, object]], assignments_path: Path, fingerprint: Mapping[str, set[str]], output: Path) -> dict[str, object]:
    hits_by_accession: dict[str, set[str]] = defaultdict(set)
    coords: dict[tuple[str, str], list[str]] = defaultdict(list)
    for hit in candidate_hits:
        accession = str(hit["accession"])
        pfam = str(hit["pfam_accession"])
        hits_by_accession[accession].add(pfam)
        coords[(accession, pfam)].append(str(hit["coordinates"]))
    with assignments_path.open(encoding="utf-8-sig", newline="") as handle:
        assignments = list(csv.DictReader(handle, delimiter="\t"))
    fields = list(assignments[0]) + ["pfam_accessions", "pfam_coordinates", "pfam_interpro_state", "architecture_consistency", "assignment_review"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in assignments:
            accession = row["accession"]
            pfams = hits_by_accession.get(accession, set())
            evidence = classify_candidate_domains(row.get("phaded_superfamily_best", ""), pfams, fingerprint)
            result = dict(row)
            result["pfam_accessions"] = ";".join(sorted(pfams))
            result["pfam_coordinates"] = ";".join(f"{pfam}:{','.join(coords[(accession, pfam)])}" for pfam in sorted(pfams))
            result.update(evidence)
            writer.writerow(result)
    return {"candidate_count": len(assignments), "candidate_domain_hit_accessions": sum(bool(v) for v in hits_by_accession.values()), "fingerprint_superfamilies": len(fingerprint), "output": {"path": str(output.resolve()), "size": output.stat().st_size, "sha256": sha256(output)}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    fp = sub.add_parser("fingerprint")
    fp.add_argument("--reference-domtblout", type=Path, required=True)
    fp.add_argument("--reference-ledger", type=Path, required=True)
    fp.add_argument("--output", type=Path, required=True)
    names = sub.add_parser("names")
    names.add_argument("--domtblout", type=Path, required=True)
    names.add_argument("--fingerprint", type=Path, required=True)
    names.add_argument("--output", type=Path, required=True)
    ap = sub.add_parser("apply")
    ap.add_argument("--candidate-domtblout", type=Path, required=True)
    ap.add_argument("--assignments", type=Path, required=True)
    ap.add_argument("--fingerprint", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "fingerprint":
        with args.reference_domtblout.open(encoding="utf-8", errors="replace") as handle:
            hits = parse_domtblout(handle)
        with args.reference_ledger.open(encoding="utf-8-sig", newline="") as handle:
            ledger = {row["reference_id"]: row["phaded_superfamily"] for row in csv.DictReader(handle, delimiter="\t")}
        fp_data = reference_fingerprint(hits, ledger)
        write_fingerprint(args.output, fp_data)
        print(json.dumps({"reference_hits": len(hits), "fingerprint": {key: sorted(value) for key, value in fp_data.items()}}, sort_keys=True))
    elif args.command == "names":
        with args.domtblout.open(encoding="utf-8", errors="replace") as handle:
            hits = parse_domtblout(handle)
        with args.fingerprint.open(encoding="utf-8-sig", newline="") as handle:
            selected = {row["pfam_accession"] for row in csv.DictReader(handle, delimiter="\t")}
        rows = selected_name_rows(hits, selected)
        if {accession for accession, _ in rows} != selected:
            raise ValueError("fingerprint Pfam accessions cannot all be resolved to NAME")
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["pfam_accession", "pfam_name"])
            writer.writerows(rows)
        print(json.dumps({"pfam_count": len(rows), "output": str(args.output)}, sort_keys=True))
    else:
        with args.candidate_domtblout.open(encoding="utf-8", errors="replace") as handle:
            hits = parse_domtblout(handle)
        with args.fingerprint.open(encoding="utf-8-sig", newline="") as handle:
            fp_data: dict[str, set[str]] = defaultdict(set)
            for row in csv.DictReader(handle, delimiter="\t"):
                fp_data[row["phaded_superfamily"]].add(row["pfam_accession"])
        print(json.dumps(apply(hits, args.assignments, fp_data, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
