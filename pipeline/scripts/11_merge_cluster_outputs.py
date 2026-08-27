#!/usr/bin/env python3
"""Merge parallel cluster batch outputs and recompute unique support counts."""
import argparse
import csv
import os
from collections import Counter, defaultdict


def read_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-root", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--marker-families", required=True)
    args = ap.parse_args()
    with open(os.path.join(args.batch_root, "batches.tsv"), newline="") as handle:
        batches = [row["batch"] for row in csv.DictReader(handle, delimiter="\t")]
    contexts = []
    audits = []
    for batch in batches:
        root = os.path.join(args.batch_root, batch, "results", "tables")
        ctx = os.path.join(root, "cluster_context.tsv")
        audit = os.path.join(root, "cluster_locus_audit.tsv")
        if not os.path.isfile(ctx) or not os.path.isfile(audit):
            raise SystemExit(f"missing batch outputs: {batch}")
        contexts.extend(read_rows(ctx))
        audits.extend(read_rows(audit))
    bad = [row for row in audits if row.get("status") != "analyzed"]
    if bad:
        counts = Counter(row.get("status") for row in bad)
        raise SystemExit("incomplete cluster batches: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    out_tables = os.path.join(args.outdir, "tables")
    os.makedirs(out_tables, exist_ok=True)
    context_cols = ["genome", "contig", "hit_locus", "hit_family", "hit_start", "hit_end", "hit_strand", "marker_locus", "marker_family", "marker_bitscore", "distance_bp", "direction", "arbitration"]
    with open(os.path.join(out_tables, "cluster_context.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=context_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(contexts, key=lambda r: tuple(r.get(k, "") for k in ("genome", "hit_locus", "marker_locus"))))
    audit_cols = ["genome", "locus", "family", "status"]
    with open(os.path.join(out_tables, "cluster_locus_audit.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(audits, key=lambda r: tuple(r.get(k, "") for k in ("genome", "locus", "family"))))
    by_key = defaultdict(list)
    for row in contexts:
        by_key[(row["hit_family"], row["marker_family"])].append(row)
    with open(os.path.join(out_tables, "cluster_summary.tsv"), "w", newline="") as handle:
        fields = ["hit_family", "marker_family", "marker_hits", "supporting_loci", "supporting_genomes"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (hit_family, marker_family), rows in sorted(by_key.items(), key=lambda item: (-len(item[1]), item[0])):
            writer.writerow({"hit_family": hit_family, "marker_family": marker_family, "marker_hits": len(rows), "supporting_loci": len({r["hit_locus"] for r in rows}), "supporting_genomes": len({r["genome"] for r in rows})})
    with open(os.path.join(out_tables, "cluster_genome_audit.tsv"), "w", newline="") as handle:
        fields = ["genome", "requested_loci", "analyzed_loci", "not_analyzed_statuses"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        grouped = defaultdict(list)
        for row in audits:
            grouped[row["genome"]].append(row)
        for genome in sorted(grouped):
            rows = grouped[genome]
            writer.writerow({"genome": genome, "requested_loci": len(rows), "analyzed_loci": sum(r["status"] == "analyzed" for r in rows), "not_analyzed_statuses": "none"})
    with open(os.path.join(args.outdir, "cluster_parallel_metadata.tsv"), "w") as handle:
        handle.write("field\tvalue\n")
        handle.write(f"batches\t{len(batches)}\nmarker_families\t{args.marker_families}\naudited_loci\t{len(audits)}\ncontext_rows\t{len(contexts)}\n")
    print(f"merged {len(batches)} batches, {len(audits)} audited loci, {len(contexts)} context rows")


if __name__ == "__main__":
    main()
