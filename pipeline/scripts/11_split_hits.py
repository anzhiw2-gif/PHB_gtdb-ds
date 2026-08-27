#!/usr/bin/env python3
"""Split locus-level hits into deterministic genome batches for parallel clustering."""
import argparse
import csv
import os
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--batches", type=int, default=80)
    args = ap.parse_args()
    if args.batches < 1:
        raise SystemExit("--batches must be positive")
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.hits, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome", "locus", "family"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"missing columns: {sorted(required - set(reader.fieldnames or []))}")
        rows = defaultdict(list)
        for row in reader:
            rows[row["genome"]].append(row)
        fields = reader.fieldnames
    genomes = sorted(rows)
    buckets = [[] for _ in range(args.batches)]
    for index, genome in enumerate(genomes):
        buckets[index % args.batches].append(genome)
    manifest = []
    for index, bucket in enumerate(buckets):
        batch = os.path.join(args.outdir, f"batch_{index:03d}")
        os.makedirs(batch, exist_ok=True)
        path = os.path.join(batch, "hits.tsv")
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for genome in bucket:
                writer.writerows(rows[genome])
        with open(os.path.join(batch, "genomes.txt"), "w") as handle:
            handle.write("\n".join(bucket) + ("\n" if bucket else ""))
        manifest.append({"batch": f"batch_{index:03d}", "genomes": len(bucket), "loci": sum(len(rows[g]) for g in bucket)})
    with open(os.path.join(args.outdir, "batches.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["batch", "genomes", "loci"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    print(f"split {len(genomes)} genomes / {sum(len(v) for v in rows.values())} loci into {args.batches} batches")


if __name__ == "__main__":
    main()
