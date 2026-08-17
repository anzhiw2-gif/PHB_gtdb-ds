#!/usr/bin/env python3
"""09d_patatin_filter.py — patatin 家族基因簇过滤
原理：真正的古菌 PhaZh1 型解聚酶是颗粒结合蛋白，与 PHA 合酶（PhaC/PhaE）
及颗粒蛋白（PhaP）共定位于 PHA 代谢基因簇。因此：
  - patatin 命中基因组的"真解聚酶"候选 = 该基因组同时含 PhaC（PHA 合酶）
  - 无 PhaC 的 patatin 命中 = 广谱磷脂酶（背景），过滤掉
输出: results/tables/patatin_filtered.tsv + 过滤后的序列
"""
import os
from collections import defaultdict

SCREEN = "data/screen"
TIER = "data/screen/tiers"
OUT = "results/tables"
os.makedirs(OUT, exist_ok=True)

# 1. 找出有 PhaC（PHA 合酶）命中的基因组（来自 hits_all 或 genome_hits）
phac_genomes = set()
# 用 hits_all.tsv（含全部家族命中，PhaC 家族）
hits_path = os.path.join(SCREEN, "hits_all.tsv")
print("读取 hits_all.tsv 中 PhaC 命中基因组...")
with open(hits_path) as f:
    header = f.readline().rstrip("\n").split("\t")
    fi = header.index("family")
    pi = header.index("protein")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(fi, pi):
            continue
        if p[fi] == "PhaC":
            g = p[pi].split("|")[0]
            phac_genomes.add(g)
print(f"含 PhaC 的基因组: {len(phac_genomes)}")

# 2. 过滤 patatin tier1
pat_faa = os.path.join(TIER, "ArchPhaZ_patatin_tier1.faa")
kept = []
dropped = []
cur_hdr = None
cur_seq = []
with open(pat_faa) as f:
    for line in f:
        if line.startswith(">"):
            if cur_hdr:
                g = cur_hdr.split("|")[0]
                (kept if g in phac_genomes else dropped).append((cur_hdr, "".join(cur_seq)))
            cur_hdr = line[1:].strip()
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_hdr:
        g = cur_hdr.split("|")[0]
        (kept if g in phac_genomes else dropped).append((cur_hdr, "".join(cur_seq)))

print(f"\npatatin tier1 总命中: {len(kept)+len(dropped)}")
print(f"含 PhaC（真解聚酶候选）: {len(kept)}")
print(f"无 PhaC（广谱磷脂酶，过滤）: {len(dropped)}")

# 写过滤后序列
out_faa = os.path.join(TIER, "ArchPhaZ_patatin_tier1_cluster.faa")
with open(out_faa, "w") as f:
    for hdr, seq in kept:
        f.write(">" + hdr + "\n" + seq + "\n")

# 统计过滤后的基因组与门分布
kept_genomes = set(h.split("|")[0] for h, _ in kept)
print(f"过滤后基因组数: {len(kept_genomes)}")

# 写 TSV
with open(f"{OUT}/patatin_filtered.tsv", "w") as f:
    f.write("genome\tpatatin_copies\n")
    cnt = defaultdict(int)
    for h, _ in kept:
        cnt[h.split("|")[0]] += 1
    for g, c in sorted(cnt.items()):
        f.write(f"{g}\t{c}\n")

print(f"\n结果: {out_faa}")
print(f"表: {OUT}/patatin_filtered.tsv")
