"""Audit hold-candidate proteins against authoritative GenBank CDS features."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable

from Bio import SeqIO


def _clean_protein(value: str) -> str:
    return "".join(str(value).split()).rstrip("*").upper()


def _qual(feature, name: str) -> str:
    values = feature.qualifiers.get(name, [])
    return str(values[0]) if values else ""


def _location_text(feature) -> str:
    return str(feature.location)


def _partial_location(feature) -> bool:
    text = _location_text(feature)
    return "<" in text or ">" in text


def _feature_row(record, feature, translation: str) -> dict[str, object]:
    location = feature.location
    return {
        "contig_accession": record.id,
        "contig_length": len(record.seq),
        "location": _location_text(feature),
        "start_1based": int(location.start) + 1,
        "end_1based": int(location.end),
        "strand": int(location.strand or 0),
        "partial_location": _partial_location(feature),
        "locus_tag": _qual(feature, "locus_tag"),
        "protein_id": _qual(feature, "protein_id"),
        "gene": _qual(feature, "gene"),
        "product": _qual(feature, "product"),
        "translation_length": len(translation),
        "translation_sha256": hashlib.sha256(translation.encode("ascii")).hexdigest(),
        "translation": translation,
    }


def iter_cds(gbff: Path) -> Iterable[dict[str, object]]:
    with gbff.open(encoding="utf-8", errors="replace") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                raw = _qual(feature, "translation")
                translation = _clean_protein(raw)
                if not translation:
                    continue
                yield _feature_row(record, feature, translation)


def summarize_gbff(gbff: Path, candidates: dict[str, str]) -> list[dict[str, object]]:
    wanted = {key: _clean_protein(value) for key, value in candidates.items()}
    matches: list[dict[str, object]] = []
    for row in iter_cds(gbff):
        for candidate_id, sequence in wanted.items():
            if row["translation"] == sequence:
                result = dict(row)
                result["candidate_id"] = candidate_id
                result.pop("translation", None)
                matches.append(result)
    return matches


def summarize_context(gbff: Path, candidates: dict[str, str]) -> dict[str, object]:
    """Return exact matches plus nearby CDS features for ORF-context review."""
    wanted = {key: _clean_protein(value) for key, value in candidates.items()}
    all_rows = list(iter_cds(gbff))
    matches: list[dict[str, object]] = []
    for row in all_rows:
        for candidate_id, sequence in wanted.items():
            if row["translation"] == sequence:
                result = dict(row)
                result["candidate_id"] = candidate_id
                result.pop("translation", None)
                matches.append(result)

    nearby: list[dict[str, object]] = []
    for match in matches:
        for row in all_rows:
            if row["contig_accession"] != match["contig_accession"]:
                continue
            if row["start_1based"] == match["start_1based"] and row["end_1based"] == match["end_1based"]:
                continue
            gap = max(
                int(match["start_1based"]) - int(row["end_1based"]) - 1,
                int(row["start_1based"]) - int(match["end_1based"]) - 1,
                0,
            )
            overlap = not (
                int(row["end_1based"]) < int(match["start_1based"])
                or int(match["end_1based"]) < int(row["start_1based"])
            )
            if overlap or gap <= 5000:
                item = dict(row)
                item.pop("translation", None)
                item["candidate_id"] = match["candidate_id"]
                item["gap_bp"] = gap
                item["overlap"] = overlap
                nearby.append(item)

    pairwise: list[dict[str, object]] = []
    for i, left in enumerate(matches):
        for right in matches[i + 1 :]:
            if left["contig_accession"] != right["contig_accession"]:
                continue
            gap = max(
                int(left["start_1based"]) - int(right["end_1based"]) - 1,
                int(right["start_1based"]) - int(left["end_1based"]) - 1,
                0,
            )
            pairwise.append(
                {
                    "left_candidate_id": left["candidate_id"],
                    "right_candidate_id": right["candidate_id"],
                    "gap_bp": gap,
                    "overlap": not (
                        int(left["end_1based"]) < int(right["start_1based"])
                        or int(right["end_1based"]) < int(left["start_1based"])
                    ),
                }
            )
    return {"matches": matches, "nearby_cds": nearby, "pairwise_matches": pairwise}


def _read_candidates(fasta: Path, ids: list[str] | None) -> dict[str, str]:
    selected = set(ids or [])
    result: dict[str, str] = {}
    with fasta.open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if selected and record.id not in selected:
                continue
            if not selected and not ("JALYXY010000192.1_6" in record.id or "NZ_CP083952.1_5059" in record.id or "NZ_CP083952.1_5067" in record.id):
                continue
            result[record.id] = str(record.seq)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-fasta", type=Path, required=True)
    parser.add_argument("--gbff", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--candidate-id", action="append", dest="candidate_ids")
    args = parser.parse_args()

    candidates = _read_candidates(args.candidate_fasta, args.candidate_ids)
    context = summarize_context(args.gbff, candidates)
    rows = context["matches"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "candidate_id", "contig_accession", "contig_length", "location",
        "start_1based", "end_1based", "strand", "partial_location",
        "locus_tag", "protein_id", "gene", "product", "translation_length",
        "translation_sha256",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    args.json.write_text(
        json.dumps(
            {"candidate_count": len(candidates), "match_count": len(rows), "candidates": sorted(candidates), "rows": rows, "nearby_cds": context["nearby_cds"], "pairwise_matches": context["pairwise_matches"]},
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
