#!/usr/bin/env python3
"""09e_patatin_phylum.py — 过滤后 patatin 的门分布"""
import os
from collections import Counter

TAX_BAC = os.path.expanduser("~/GTDB/taxonomy/bac120_taxonomy_r232.tsv")
TAX_AR = os.path.expanduser("~/GTDB/taxonomy/ar53_taxonomy_r232.tsv")
tax = {}
for tf in [TAX_BAC, TAX_AR]:
    for line in open(tf):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            tax[p[0]] = p[1]

def lookup(g):
    if g in tax: return tax[g]
    if g.startswith("GCA_"): return tax.get("GB_"+g, "unclassified")
    if g.startswith("GCF_"): return tax.get("RS_"+g, "unclassified")
    return "unclassified"

# 过滤后 patatin 基因组
genomes = set()
for line in open("results/tables/patatin_filtered.tsv"):
    if line.startswith("genome"): continue
    genomes.add(line.split("\t")[0])

phylum = Counter()
arch = Counter()
for g in genomes:
    t = lookup(g)
    parts = t.split(";")
    ph = parts[1].replace("p__", "") if len(parts) > 1 else "unknown"
    dom = parts[0].replace("d__", "") if len(parts) > 0 else "unknown"
    phylum[ph] += 1
    arch[dom] += 1

print(f"过滤后 patatin 基因组: {len(genomes)}")
print(f"\n域分布: {dict(arch)}")
print("\n门分布（top 15）:")
for ph, n in phylum.most_common(15):
    print(f"  {ph}: {n}")
