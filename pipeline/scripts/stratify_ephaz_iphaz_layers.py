#!/usr/bin/env python3
"""Create auditable ePhaZ/iPhaZ candidate layers from fixed tier inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def layer_for(row: dict[str, str]) -> str:
    family = row.get("family", "").strip()
    source = row.get("source_layer", "").strip().lower()
    signal = row.get("signal_type", "").strip().upper()
    assignment = row.get("assignment", "").strip().lower()
    if family == "ePhaZ":
        if assignment == "ambiguous":
            return "ePhaZ_competition_review"
        if source == "ephaz_broad_discovery":
            return "ePhaZ_broad_discovery"
        if source == "tier2":
            return "ePhaZ_tier2_review"
        if source != "ephaz_curated_core":
            return "ePhaZ_review_pending"
        if signal in {"SP", "LIPO", "TAT", "TATLIPO"}:
            return "ePhaZ_curated_core_secreted"
        return "ePhaZ_curated_core_nonsecreted"
    if family == "iPhaZ":
        if assignment == "ambiguous":
            return "iPhaZ_competition_review"
        if source == "tier1":
            return "iPhaZ_tier1"
        if source == "tier2":
            return "iPhaZ_tier2_review"
        return "iPhaZ_review_pending"
    raise ValueError(f"unsupported family: {family}")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records = {}
    header, sequence = None, []
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if header is not None:
                key = header.split(None, 1)[0]
                if key in records:
                    raise ValueError(f"duplicate FASTA accession: {key}")
                records[key] = (header, "".join(sequence))
            header, sequence = line[1:], []
        elif line:
            if header is None:
                raise ValueError(f"sequence precedes header: {path}")
            sequence.append(line)
    if header is not None:
        key = header.split(None, 1)[0]
        if key in records:
            raise ValueError(f"duplicate FASTA accession: {key}")
        records[key] = (header, "".join(sequence))
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    return records


def _signal_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    result = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if not {"accession", "type"}.issubset(fields):
            raise ValueError("SignalP table requires accession and type columns")
        for row in reader:
            accession = row["accession"].strip()
            if accession in result:
                raise ValueError(f"duplicate SignalP accession: {accession}")
            result[accession] = row["type"].strip()
    return result


def _competition_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    result = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if not {"accession", "assignment"}.issubset(fields):
            raise ValueError("competition table requires accession and assignment columns")
        for row in reader:
            accession = row["accession"].strip()
            if accession in result:
                raise ValueError(f"duplicate competition accession: {accession}")
            result[accession] = row["assignment"].strip()
    return result


def _layer_priority(layer: str) -> int:
    """Return deterministic evidence priority for duplicate protein accessions."""
    return {
        "ePhaZ_competition_review": 0,
        "iPhaZ_competition_review": 0,
        "ePhaZ_curated_core_secreted": 1,
        "ePhaZ_curated_core_nonsecreted": 2,
        "iPhaZ_tier1": 1,
        "ePhaZ_tier2_review": 3,
        "iPhaZ_tier2_review": 3,
        "ePhaZ_broad_discovery": 4,
        "ePhaZ_review_pending": 5,
        "iPhaZ_review_pending": 5,
    }.get(layer, 99)


def stratify(ephaz_faa: Path, iphaz_faa: Path, output_dir: Path, signalp: Path | None = None, ephaz_tier2: Path | None = None, iphaz_tier2: Path | None = None, competition: Path | None = None) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=False)
    signal = _signal_map(signalp)
    competition_assignments = _competition_map(competition)
    rows = []
    inputs = [("ePhaZ", ephaz_faa, "ePhaZ_curated_core"), ("iPhaZ", iphaz_faa, "tier1")]
    if ephaz_tier2 is not None:
        inputs.append(("ePhaZ", ephaz_tier2, "tier2"))
    if iphaz_tier2 is not None:
        inputs.append(("iPhaZ", iphaz_tier2, "tier2"))
    for family, path, source in inputs:
        for accession, (header, sequence) in sorted(read_fasta(path).items()):
            signal_type = signal.get(accession, "OTHER") if family == "ePhaZ" else ""
            assignment = competition_assignments.get(accession, "")
            row = {"accession": accession, "genome": accession.split("|", 1)[0], "family": family, "source_layer": source, "signal_type": signal_type, "assignment": assignment, "length": str(len(sequence))}
            row["layer"] = layer_for(row)
            rows.append(row)
    raw_row_count = len(rows)
    selected = {}
    for row in rows:
        key = (row["family"], row["accession"])
        incumbent = selected.get(key)
        if incumbent is None or (_layer_priority(row["layer"]), row["source_layer"]) < (_layer_priority(incumbent["layer"]), incumbent["source_layer"]):
            selected[key] = row
    rows = sorted(selected.values(), key=lambda row: (row["family"], row["accession"]))
    columns = ["accession", "genome", "family", "source_layer", "signal_type", "assignment", "length", "layer"]
    with (output_dir / "protein_layers.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    genome_layers = defaultdict(set)
    genome_counts = Counter()
    for row in rows:
        genome_layers[(row["genome"], row["family"])].add(row["layer"])
        genome_counts[(row["genome"], row["family"])] += 1
    genome_rows = []
    for (genome, family), layers in sorted(genome_layers.items()):
        genome_rows.append({"genome": genome, "family": family, "layers": ";".join(sorted(layers)), "protein_count": str(genome_counts[(genome, family)])})
    with (output_dir / "genome_layers.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "family", "layers", "protein_count"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(genome_rows)
    counts = Counter(row["layer"] for row in rows)
    metadata_inputs = {"ephaz_faa": {"path": str(ephaz_faa.resolve()), "bytes": ephaz_faa.stat().st_size, "sha256": _sha256(ephaz_faa)}, "iphaz_faa": {"path": str(iphaz_faa.resolve()), "bytes": iphaz_faa.stat().st_size, "sha256": _sha256(iphaz_faa)}}
    for name, path in (("ephaz_tier2", ephaz_tier2), ("iphaz_tier2", iphaz_tier2), ("signalp", signalp), ("competition", competition)):
        metadata_inputs[name] = None if path is None else {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": _sha256(path)}
    metadata = {"schema_version": 1, "status": "candidate-only", "family_scope": ["ePhaZ", "iPhaZ"], "counts": dict(sorted(counts.items())), "inputs": metadata_inputs, "outputs": ["protein_layers.tsv", "genome_layers.tsv"], "deduplication": {"key": ["family", "accession"], "priority": "competition_review > curated/tier1 > tier2 > broad/pending", "raw_rows": raw_row_count, "retained_rows": len(rows), "duplicate_rows_removed": raw_row_count - len(rows)}, "formal_registry_modified": False, "formal_scan_started": False}
    (output_dir / "stratification_metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ephaz-faa", type=Path, required=True)
    parser.add_argument("--iphaz-faa", type=Path, required=True)
    parser.add_argument("--signalp", type=Path)
    parser.add_argument("--ephaz-tier2", type=Path)
    parser.add_argument("--iphaz-tier2", type=Path)
    parser.add_argument("--competition", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    stratify(args.ephaz_faa, args.iphaz_faa, args.output_dir, args.signalp, args.ephaz_tier2, args.iphaz_tier2, args.competition)


if __name__ == "__main__":
    main()
