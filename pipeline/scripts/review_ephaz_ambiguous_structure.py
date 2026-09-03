#!/usr/bin/env python3
"""Create a conservative, auditable structural review of sampled candidates.

The review is sequence/domain evidence only.  It does not assert a PHB
degradation phenotype and does not mutate any control or candidate table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


STRONG_IEVALUE = 1e-3
REQUIRED_SAMPLE = {"accession", "length", "ephaz_bitscore", "iphaz_bitscore", "assignment"}


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if start > end:
            continue
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def coverage_from_intervals(intervals: list[tuple[int, int]], model_length: int) -> float:
    if model_length <= 0:
        return 0.0
    return min(1.0, sum(end - start + 1 for start, end in merge_intervals(intervals)) / model_length)


def classify_architecture(
    ephaz_intervals: list[tuple[int, int]],
    iphaz_intervals: list[tuple[int, int]],
    ephaz_coverage: float,
    iphaz_coverage: float,
) -> str:
    """Assign a conservative architecture label from strong HMM segments."""
    e_n = any(start <= 150 and end >= 40 for start, end in ephaz_intervals)
    e_c = any(start >= 150 and end >= 180 for start, end in ephaz_intervals)
    e_two = e_n and e_c
    i_central = any(start <= 236 and end >= 300 for start, end in iphaz_intervals)
    if e_two and i_central:
        return "mixed_cross_family"
    if e_two and not i_central:
        return "ePhaZ_consistent"
    if i_central and not e_two:
        return "iPhaZ_consistent"
    if e_n and not i_central:
        return "partial_ePhaZ_signal"
    if ephaz_coverage >= 0.35 and iphaz_coverage < 0.30:
        return "partial_ePhaZ_signal"
    if iphaz_coverage >= 0.30:
        return "partial_iPhaZ_signal"
    return "insufficient_structural_support"


def review_decision(integrity_status: str, architecture_evidence: str) -> str:
    if integrity_status != "complete":
        return "pending_manual"
    if architecture_evidence == "iPhaZ_consistent":
        return "provisional_iPhaZ_challenge"
    if architecture_evidence == "ePhaZ_consistent":
        return "provisional_ePhaZ_review"
    return "pending_manual"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": _sha256(path)}


def _read_sample(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_SAMPLE - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"sample TSV missing columns: {sorted(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows or len({row["accession"] for row in rows}) != len(rows):
        raise ValueError("sample TSV must contain unique non-empty accessions")
    if any(row["assignment"] != "ambiguous" for row in rows):
        raise ValueError("sample TSV must contain only ambiguous rows")
    return rows


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    sequence: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                key = header.split(None, 1)[0]
                if key in records:
                    raise ValueError(f"duplicate FASTA accession: {key}")
                records[key] = "".join(sequence)
            header, sequence = line[1:], []
        elif line:
            if header is None:
                raise ValueError("FASTA sequence precedes header")
            sequence.append(line)
    if header is not None:
        key = header.split(None, 1)[0]
        if key in records:
            raise ValueError(f"duplicate FASTA accession: {key}")
        records[key] = "".join(sequence)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _read_signalp(path: Path) -> dict[str, str]:
    predictions = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) >= 2:
                predictions[fields[0].strip()] = fields[1].strip()
    return predictions


def _read_domtblout(paths: list[Path], accessions: set[str]) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = defaultdict(lambda: {"all": [], "strong": [], "qlen": 0.0, "best_score": 0.0})
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip() or raw.lstrip().startswith("#"):
                    continue
                fields = raw.split()
                if len(fields) < 22 or fields[0] not in accessions:
                    continue
                try:
                    accession = fields[0]
                    qlen = float(fields[5])
                    i_evalue = float(fields[12])
                    score = float(fields[13])
                    start, end = int(fields[15]), int(fields[16])
                except (IndexError, ValueError):
                    continue
                if qlen <= 0 or start < 1 or end < start:
                    continue
                item = result[accession]
                item["qlen"] = qlen
                item["best_score"] = max(float(item["best_score"]), score)
                item["all"].append((start, end))
                if i_evalue <= STRONG_IEVALUE:
                    item["strong"].append((start, end))
    return result


def _read_neighborhood(path: Path, accessions: set[str]) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = f"{(row.get('genome') or '').strip()}|{(row.get('hit_locus') or '').strip()}"
            if key in accessions:
                result[key].append(row)
    return result


def review(sample_tsv: Path | str, sample_faa: Path | str, ephaz_domtblout: Path | str, iphaz_domtblout: Path | str, signalp: Path | str, neighborhood: Path | str, outdir: Path | str) -> dict[str, object]:
    paths = [Path(item) for item in (sample_tsv, sample_faa, ephaz_domtblout, iphaz_domtblout, signalp, neighborhood)]
    rows = _read_sample(paths[0])
    accessions = {row["accession"] for row in rows}
    sequences = _read_fasta(paths[1])
    missing = accessions - set(sequences)
    if missing:
        raise ValueError(f"sample accessions missing from FASTA: {sorted(missing)[:3]}")
    ephaz = _read_domtblout([paths[2]], accessions)
    iphaz = _read_domtblout([paths[3]], accessions)
    signalp_predictions = _read_signalp(paths[4])
    neighborhood_rows = _read_neighborhood(paths[5], accessions)
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=False)
    reviewed = []
    for row in rows:
        accession = row["accession"]
        sequence = sequences[accession]
        integrity = "complete"
        if not sequence.startswith("M"):
            integrity = "possible_N_truncation"
        elif "*" in sequence[:-1]:
            integrity = "internal_stop"
        e_item, i_item = ephaz.get(accession, {}), iphaz.get(accession, {})
        e_strong = merge_intervals(list(e_item.get("strong", [])))
        i_strong = merge_intervals(list(i_item.get("strong", [])))
        e_all = merge_intervals(list(e_item.get("all", [])))
        i_all = merge_intervals(list(i_item.get("all", [])))
        e_cov = coverage_from_intervals(e_strong, int(e_item.get("qlen", 317) or 317))
        i_cov = coverage_from_intervals(i_strong, int(i_item.get("qlen", 455) or 455))
        architecture = classify_architecture(e_strong, i_strong, e_cov, i_cov)
        decision = review_decision(integrity, architecture)
        nrows = neighborhood_rows.get(accession, [])
        marker_rows = [item for item in nrows if (item.get("marker_family") or "").strip()]
        marker_families = sorted({item.get("marker_family", "").strip() for item in marker_rows if item.get("marker_family", "").strip()})
        distances = []
        for item in marker_rows:
            try:
                distances.append(int(float(item.get("distance_bp", ""))))
            except (TypeError, ValueError):
                pass
        reviewed.append({
            "accession": accession,
            "length": row["length"],
            "delta_abs": row.get("delta_abs", ""),
            "ephaz_bitscore": row["ephaz_bitscore"],
            "iphaz_bitscore": row["iphaz_bitscore"],
            "sequence_start": sequence[:1],
            "terminal_stop": "yes" if sequence.endswith("*") else "no",
            "internal_stop": "yes" if "*" in sequence[:-1] else "no",
            "integrity_status": integrity,
            "signalp_prediction": signalp_predictions.get(accession, "pending"),
            "ephaz_reported_segments": len(e_all),
            "iphaz_reported_segments": len(i_all),
            "ephaz_strong_segments": len(e_strong),
            "iphaz_strong_segments": len(i_strong),
            "ephaz_strong_hmm_coverage": f"{e_cov:.6f}",
            "iphaz_strong_hmm_coverage": f"{i_cov:.6f}",
            "ephaz_strong_hmm_blocks": ";".join(f"{s}-{e}" for s, e in e_strong),
            "iphaz_strong_hmm_blocks": ";".join(f"{s}-{e}" for s, e in i_strong),
            "architecture_evidence": architecture,
            "neighborhood_status": "marker_present" if marker_rows else ("record_no_marker" if nrows else "no_record"),
            "neighborhood_marker_families": ";".join(marker_families),
            "neighborhood_best_distance_bp": "" if not distances else str(min(abs(value) for value in distances)),
            "manual_decision": decision,
            "review_note": "HMM/domain evidence only; no phenotype claim",
        })
    fields = list(reviewed[0])
    tsv = output / "ambiguous_structural_review.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(reviewed)
    counts = defaultdict(int)
    for item in reviewed:
        counts[item["manual_decision"]] += 1
    arch_counts = defaultdict(int)
    for item in reviewed:
        arch_counts[item["architecture_evidence"]] += 1
    metadata = {
        "schema_version": 1,
        "tool": "review_ephaz_ambiguous_structure.py",
        "parameters": {"strong_domain_i_evalue": STRONG_IEVALUE, "ephaz_model_length": 317, "iphaz_model_length": 455},
        "counts": {"sampled": len(reviewed), "manual_decision": dict(sorted(counts.items())), "architecture_evidence": dict(sorted(arch_counts.items()))},
        "inputs": {name: _file_record(path) for name, path in zip(("sample_tsv", "sample_faa", "ephaz_domtblout", "iphaz_domtblout", "signalp", "neighborhood"), paths)},
        "outputs": {"ambiguous_structural_review.tsv": _file_record(tsv)},
        "interpretation": "Provisional architecture labels are sequence/domain evidence and are not phenotype validation.",
    }
    meta = output / "structural_review_metadata.json"
    meta.write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    meta.write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-tsv", required=True)
    parser.add_argument("--sample-faa", required=True)
    parser.add_argument("--ephaz-domtblout", required=True)
    parser.add_argument("--iphaz-domtblout", required=True)
    parser.add_argument("--signalp", required=True)
    parser.add_argument("--neighborhood", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args(argv)
    try:
        metadata = review(**vars(args))
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(metadata["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
