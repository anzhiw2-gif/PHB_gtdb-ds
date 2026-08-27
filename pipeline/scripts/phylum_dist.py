#!/usr/bin/env python3
"""phylum_dist.py — 对任意 FASTA 子集统计门（phylum）水平分布

输入 FASTA 的 header 第一段（| 前）为 GTDB genome accession（GCA_/GCF_）。
读 GTDB taxonomy（服务器 ~/GTDB/taxonomy/）映射到门，输出分布 TSV。

用法（服务器 T141）:
  ~/miniconda3/envs/phb_gtdb/bin/python scripts/phylum_dist.py \
      --faa data/screen/tiers/ePhaZ_tier1_signalpeptide.faa \
      --out results/tables/ePhaZ_signalpeptide_phylum.tsv
"""
import argparse
import os
from collections import Counter

TAX_BAC = os.path.expanduser("~/GTDB/taxonomy/bac120_taxonomy_r232.tsv")
TAX_AR = os.path.expanduser("~/GTDB/taxonomy/ar53_taxonomy_r232.tsv")


def load_taxonomy():
    tax = {}
    for tf in (TAX_BAC, TAX_AR):
        if os.path.exists(tf):
            for line in open(tf):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    tax[p[0]] = p[1]
    return tax


def lookup_tax(g, tax):
    if g in tax:
        return tax[g]
    if g.startswith("GCA_"):
        return tax.get("GB_" + g, "unclassified")
    if g.startswith("GCF_"):
        return tax.get("RS_" + g, "unclassified")
    return "unclassified"


def phylum_of(tax_str):
    parts = tax_str.split(";")
    if len(parts) > 1:
        return parts[1].replace("p__", "")
    return "unclassified"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--faa", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    tax = load_taxonomy()

    # 读 FASTA，提取 genome
    genomes = []
    for line in open(args.faa):
        if line.startswith(">"):
            genomes.append(line[1:].strip().split("|")[0])

    # 映射到门
    phyla = Counter()
    for g in genomes:
        phyla[phylum_of(lookup_tax(g, tax))] += 1

    # 写 TSV
    with open(args.out, "w") as f:
        f.write("phylum\tgenomes\n")
        for ph, n in phyla.most_common():
            f.write(f"{ph}\t{n}\n")

    print(f"序列数: {len(genomes)}")
    print(f"门分布（top 15，共 {len(phyla)} 门）:")
    for ph, n in phyla.most_common(15):
        print(f"  {ph}: {n} ({n/len(genomes)*100:.1f}%)")
    print(f"结果 -> {args.out}")


if __name__ == "__main__":
    main()
