#!/usr/bin/env python3
"""09f_archaeal_patatin.py — 提取古菌 patatin 解聚酶子集（真正的 PhaZh1 型）"""
import os
from collections import Counter

TAX_AR = os.path.expanduser("~/GTDB/taxonomy/ar53_taxonomy_r232.tsv")
tax = {}
for line in open(TAX_AR):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2:
        tax[p[0]] = p[1]

def is_archaea(g):
    t = tax.get(g) or tax.get("GB_"+g) or tax.get("RS_"+g)
    return t is not None and t.startswith("d__Archaea")

def arch_phylum(g):
    t = tax.get(g) or tax.get("GB_"+g) or tax.get("RS_"+g)
    if not t: return "unknown"
    parts = t.split(";")
    return parts[1].replace("p__", "") if len(parts) > 1 else "unknown"

# 从 cluster 过滤后的 patatin 中提取古菌
faa = "data/screen/tiers/ArchPhaZ_patatin_tier1_cluster.faa"
out = "data/screen/tiers/ArchPhaZ_patatin_archaea.faa"
arch_genomes = Counter()
kept = 0
cur_hdr = None
cur_seq = []
with open(faa) as fin, open(out, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            if cur_hdr and is_archaea(cur_hdr.split("|")[0]):
                fout.write(">" + cur_hdr + "\n" + "".join(cur_seq) + "\n")
                kept += 1
                arch_genomes[cur_hdr.split("|")[0]] += 1
            cur_hdr = line[1:].strip()
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_hdr and is_archaea(cur_hdr.split("|")[0]):
        fout.write(">" + cur_hdr + "\n" + "".join(cur_seq) + "\n")
        kept += 1
        arch_genomes[cur_hdr.split("|")[0]] += 1

print(f"古菌 patatin 解聚酶（PhaZh1 型）: {kept} 序列 / {len(arch_genomes)} 基因组")

# 门分布
phylum = Counter()
for g in arch_genomes:
    phylum[arch_phylum(g)] += 1
print("\n古菌门分布:")
for ph, n in phylum.most_common():
    print(f"  {ph}: {n} 基因组")

# 纲分布（更细）
cls = Counter()
for g in arch_genomes:
    t = tax.get(g) or tax.get("GB_"+g) or tax.get("RS_"+g)
    if t:
        parts = t.split(";")
        if len(parts) > 2:
            cls[parts[2].replace("c__", "")] += 1
print("\n古菌纲分布:")
for c, n in cls.most_common(10):
    print(f"  {c}: {n}")
