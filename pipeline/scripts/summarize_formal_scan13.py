#!/usr/bin/env python3
"""Summarize completed formal scan 13 without changing raw HMMER outputs."""
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


CORE = {"ePhaZ_curated_core": "ePhaZ", "ePhaZ_broad_discovery": "ePhaZ",
        "iPhaZ": "iPhaZ", "OH": "OH", "ArchPhaZ_hydrolase": "ArchPhaZ_hydrolase"}


def load_registry(path: Path) -> dict[str, tuple[float, float]]:
    out = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["model"]] = (float(row["threshold"].replace("e-", "1e-")), float(row["min_cov"]))
    return out


def taxonomy(path: Path) -> dict[str, str]:
    out = {}
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                out[p[0]] = p[1]
    return out


def lookup_tax(acc: str, tax: dict[str, str]) -> str:
    return tax.get(acc, tax.get(("GB_" if acc.startswith("GCA_") else "RS_") + acc, "unclassified"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", type=Path, required=True)
    ap.add_argument("--registry", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--taxonomy", type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    reg = load_registry(args.registry)
    rows = []
    with args.hits.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            model = row["family"]
            if model not in reg:
                continue
            ev, cov = float(row["E-value"]), float(row.get("cov", "0"))
            threshold, min_cov = reg[model]
            if ev <= threshold and cov >= min_cov:
                protein = row["protein"]
                genome = protein.split("|", 1)[0]
                rows.append((model, genome, protein, ev, cov, float(row["score"])))

    model_counts = Counter(r[0] for r in rows)
    model_proteins = defaultdict(set)
    model_genomes = defaultdict(set)
    genome_families = defaultdict(set)
    for model, genome, protein, *_ in rows:
        model_proteins[model].add(protein); model_genomes[model].add(genome)
        if model in CORE:
            genome_families[genome].add(CORE[model])

    with (args.outdir / "model_summary.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["model", "accepted_hits", "unique_proteins", "unique_genomes"])
        for model in reg:
            w.writerow([model, model_counts[model], len(model_proteins[model]), len(model_genomes[model])])

    with (args.outdir / "genome_family.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["genome", "family", "models"])
        for genome in sorted(genome_families):
            for fam in sorted(genome_families[genome]):
                models = sorted(m for m in CORE if CORE[m] == fam and genome in model_genomes[m])
                w.writerow([genome, fam, ",".join(models)])

    union = Counter()
    for genome, fams in genome_families.items():
        union["core_union"] += 1
        for fam in fams: union[fam] += 1
    with (args.outdir / "genome_union_summary.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["set", "genomes"])
        for key in ["ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase", "core_union"]:
            w.writerow([key, union[key]])

    # Genome-level co-occurrence among the four core families.
    combos = Counter("+".join(sorted(fams)) for fams in genome_families.values())
    with (args.outdir / "core_cooccurrence.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["family_set", "genomes"])
        for combo, n in sorted(combos.items(), key=lambda x: (-x[1], x[0])):
            w.writerow([combo, n])

    if args.taxonomy:
        tax = taxonomy(args.taxonomy)
        with (args.outdir / "phylum_family.tsv").open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t"); w.writerow(["phylum", "family", "genomes"])
            counts = Counter()
            for genome, fams in genome_families.items():
                t = lookup_tax(genome, tax); parts = t.split(";")
                phylum = parts[1].removeprefix("p__") if len(parts) > 1 else "unknown"
                for fam in fams: counts[(phylum, fam)] += 1
            for (phylum, fam), n in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
                w.writerow([phylum, fam, n])

    (args.outdir / "summary.txt").write_text(
        f"accepted_rows={len(rows)}\ncore_union_genomes={union['core_union']}\n"
        "scope=raw HMM hits accepted by registry thresholds; not tier1 validation\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
