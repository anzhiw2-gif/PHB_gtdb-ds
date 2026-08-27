#!/usr/bin/env python3
"""Export the archaeal subset of ArchPhaZ_hydrolase tier1 evidence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


FAMILY = "ArchPhaZ_hydrolase"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path):
    record = None
    sequence = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if record is not None:
                    yield record, "".join(sequence)
                record, sequence = line[1:].strip(), []
            elif record is not None:
                sequence.append(line.strip())
        if record is not None:
            yield record, "".join(sequence)


def export_archaea(table_path, fasta_path, out_faa, out_tsv, provenance_path):
    selected = {}
    with open(table_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome", "family", "copies", "gtdb_taxonomy", "phylum", "class"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing tier table columns: {', '.join(sorted(missing))}")
        for row in reader:
            taxonomy = row["gtdb_taxonomy"]
            if row["family"] == FAMILY and taxonomy.startswith("d__Archaea;"):
                key = row["genome"]
                selected[key] = row

    records = {}
    for header, sequence in read_fasta(fasta_path):
        genome, sep, locus = header.partition("|")
        if not sep or not genome or not locus:
            raise ValueError(f"invalid FASTA header (expected genome|locus): {header}")
        if genome in selected:
            records[(genome, locus)] = (header, sequence)

    expected = sum(int(row["copies"]) for row in selected.values())
    if len(records) != expected:
        raise ValueError(f"missing FASTA records: expected {expected}, found {len(records)}")

    out_faa = Path(out_faa)
    out_tsv = Path(out_tsv)
    provenance_path = Path(provenance_path)
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(records.items())
    with open(out_faa, "w", encoding="utf-8") as handle:
        for (_, _), (header, sequence) in ordered:
            handle.write(f">{header}\n{sequence}\n")
    columns = ["genome", "family", "locus", "copies", "gtdb_taxonomy", "phylum", "class"]
    with open(out_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (genome, locus), _ in ordered:
            row = dict(selected[genome])
            row["locus"] = locus
            writer.writerow({column: row.get(column, "") for column in columns})

    payload = {
        "family": FAMILY,
        "filter": {"domain": "d__Archaea", "taxonomy_field": "gtdb_taxonomy", "family_field": "family"},
        "records": len(ordered),
        "genomes": len({genome for genome, _ in records}),
        "inputs": {"tier1_table": str(table_path), "tier1_table_sha256": sha256(table_path), "tier1_faa": str(fasta_path), "tier1_faa_sha256": sha256(fasta_path)},
        "outputs": {"faa_sha256": sha256(out_faa), "table_sha256": sha256(out_tsv)},
    }
    provenance_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", default="results/tables/tier1_genome_family.tsv")
    parser.add_argument("--fasta", default="data/screen/tiers/ArchPhaZ_hydrolase_tier1.faa")
    parser.add_argument("--outdir", default="results/archphaz_hydrolase_archaea")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    payload = export_archaea(Path(args.table), Path(args.fasta), outdir / "ArchPhaZ_hydrolase_tier1_archaea.faa", outdir / "ArchPhaZ_hydrolase_tier1_archaea.tsv", outdir / "provenance.json")
    print(f"exported {payload['records']} sequences from {payload['genomes']} archaeal genomes")


if __name__ == "__main__":
    main()
