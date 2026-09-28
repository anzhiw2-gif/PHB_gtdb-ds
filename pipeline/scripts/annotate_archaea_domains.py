#!/usr/bin/env python3
"""Run a bounded, candidate-only cross-HMM annotation for archaeal PhaZh1 candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path


ALLOWED_LAYERS = {"PhaZh1_like_high", "PhaZh1_like_review"}
ALLOWED_ROLES = {"positive_domain", "conflict_domain", "context_marker"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_candidate_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome", "locus", "layer", "phylum", "nearby_markers"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"candidate table missing columns: {', '.join(sorted(missing))}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    seen = set()
    for row in rows:
        key = f"{row['genome']}|{row['locus']}"
        if not row["genome"] or not row["locus"]:
            raise ValueError("candidate table contains empty genome or locus")
        if key in seen:
            raise ValueError(f"duplicate locus: {key}")
        if row["layer"] not in ALLOWED_LAYERS:
            raise ValueError(f"unsupported candidate layer: {row['layer']}")
        seen.add(key)
    if not rows:
        raise ValueError("candidate table is empty")
    return rows


def read_model_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"family", "path", "version", "sha256", "role"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"model manifest missing columns: {', '.join(sorted(missing))}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError("model manifest is empty")
    families = set()
    positive = 0
    for row in rows:
        model = Path(row["path"])
        if not row["family"] or row["family"] in families:
            raise ValueError(f"missing or duplicate model family: {row['family']}")
        if row["role"] not in ALLOWED_ROLES:
            raise ValueError(f"unsupported model role: {row['role']}")
        if not row["version"] or row["version"].lower() == "pending":
            raise ValueError(f"model version missing for {row['family']}")
        if not model.is_file() or model.stat().st_size == 0:
            raise ValueError(f"model file missing or empty: {model}")
        observed = sha256_file(model)
        if row["sha256"].lower() != observed.lower():
            raise ValueError(f"model sha256 mismatch for {row['family']}")
        families.add(row["family"])
        positive += row["role"] == "positive_domain"
    if not positive:
        raise ValueError("model manifest has no positive_domain model")
    return rows


def parse_domtblout(
    lines: list[str], *, family_override: str | None = None, max_i_evalue: float | None = None
) -> list[dict[str, object]]:
    hits = []
    for raw in lines:
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.split()
        if len(fields) < 22:
            continue
        try:
            i_evalue = float(fields[12])
            if max_i_evalue is not None and i_evalue > max_i_evalue:
                continue
            hits.append({
                "accession": fields[0],
                "family": family_override or fields[3],
                "i_evalue": i_evalue,
                "bitscore": float(fields[13]),
            })
        except ValueError:
            continue
    return hits


def classify_domain_status(
    hit_counts: dict[str, int],
    conflict_families: set[str],
    supported_families: set[str] | None = None,
) -> str:
    if any(hit_counts.get(family, 0) for family in conflict_families):
        return "domain_conflict"
    supported_families = supported_families or {"ArchPhaZ_patatin"}
    if any(hit_counts.get(family, 0) for family in supported_families):
        return "domain_supported"
    return "no_registered_domain_hit"


def search_command(hmmsearch: str, model: Path, fasta: Path, domtblout: Path) -> list[str]:
    return [
        hmmsearch, "--noali", "-E", "1e-5", "--domE", "1e-5", "--domtblout",
        str(domtblout), str(model), str(fasta),
    ]


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    sequence: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(sequence)
            header, sequence = line[1:].split()[0], []
        elif line:
            if header is None:
                raise ValueError(f"FASTA sequence precedes header: {path}")
            sequence.append(line)
    if header is not None:
        records[header] = "".join(sequence)
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _write_fasta(path: Path, records: dict[str, str], accessions: set[str]) -> None:
    missing = sorted(accessions - set(records))
    if missing:
        raise ValueError(f"candidate accessions missing from FASTA: {', '.join(missing[:3])}")
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(accessions):
            handle.write(f">{accession}\n{records[accession]}\n")


def annotate(
    candidate_tsv: Path,
    source_fasta: Path,
    positive_controls: Path,
    negative_controls: Path,
    model_manifest: Path,
    outdir: Path,
    *,
    hmmsearch: str = "hmmsearch",
) -> dict[str, object]:
    rows = read_candidate_table(candidate_tsv)
    models = read_model_manifest(model_manifest)
    source_records = _read_fasta(source_fasta)
    accessions = {f"{row['genome']}|{row['locus']}" for row in rows}
    outdir.mkdir(parents=True, exist_ok=False)
    inputs = outdir / "inputs"
    raw = outdir / "raw"
    inputs.mkdir(); raw.mkdir()
    candidate_faa = inputs / "candidate_high_review.faa"
    _write_fasta(candidate_faa, source_records, accessions)
    targets = {
        "candidate": candidate_faa,
        "positive_control": positive_controls,
        "negative_control": negative_controls,
    }
    for label, fasta in targets.items():
        if not fasta.is_file() or fasta.stat().st_size == 0:
            raise ValueError(f"{label} FASTA missing or empty: {fasta}")
    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    calibration = []
    commands = []
    for model in models:
        model_path = Path(model["path"])
        for label, fasta in targets.items():
            output = raw / f"{model['family']}.{label}.domtblout"
            command = search_command(hmmsearch, model_path, fasta, output)
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            commands.append(command)
            parsed = parse_domtblout(
                output.read_text(encoding="utf-8", errors="replace").splitlines(),
                family_override=model["family"],
                max_i_evalue=1e-5,
            )
            for hit in parsed:
                if label == "candidate":
                    hits[str(hit["accession"])].append(hit)
                else:
                    calibration.append({"panel": label, **hit})
    conflict_families = {row["family"] for row in models if row["role"] == "conflict_domain"}
    supported_families = {row["family"] for row in models if row["role"] == "positive_domain"}
    output_rows = []
    for row in rows:
        accession = f"{row['genome']}|{row['locus']}"
        per_family: dict[str, int] = defaultdict(int)
        best: dict[str, float] = {}
        for hit in hits.get(accession, []):
            family = str(hit["family"])
            per_family[family] += 1
            best[family] = max(best.get(family, float("-inf")), float(hit["bitscore"]))
        output_rows.append({
            **row,
            "accession": accession,
            "registered_model_annotation": classify_domain_status(
                per_family, conflict_families, supported_families
            ),
            "registered_model_hits": ";".join(sorted(per_family)),
            "registered_model_best_bitscore": ";".join(f"{name}:{best[name]:.3f}" for name in sorted(best)),
            "independent_domain_annotation": "pending_no_registered_pfam_or_interpro_database",
            "phenotype_claim": "candidate_only",
        })
    result_tsv = outdir / "domain_annotation.tsv"
    with result_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)
    calibration_tsv = outdir / "positive_negative_calibration.tsv"
    with calibration_tsv.open("w", encoding="utf-8", newline="") as handle:
        fields = ["panel", "accession", "family", "i_evalue", "bitscore"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(calibration)
    metadata = {
        "status": "completed_candidate_only",
        "candidate_counts": {layer: sum(row["layer"] == layer for row in rows) for layer in sorted(ALLOWED_LAYERS)},
        "models": models,
        "commands": commands,
        "independent_domain_annotation": "pending_no_registered_pfam_or_interpro_database",
        "phenotype_boundary": "Registered HMM cross-annotation and neighborhood evidence remain candidate-only and do not establish PHB degradation phenotype.",
    }
    (outdir / "annotation_manifest.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-tsv", type=Path, required=True)
    parser.add_argument("--source-fasta", type=Path, required=True)
    parser.add_argument("--positive-controls", type=Path, required=True)
    parser.add_argument("--negative-controls", type=Path, required=True)
    parser.add_argument("--model-manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--hmmsearch", default="hmmsearch")
    args = parser.parse_args()
    print(json.dumps(annotate(**vars(args)), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
