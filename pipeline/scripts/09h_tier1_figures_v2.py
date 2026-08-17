#!/usr/bin/env python3
"""09h_tier1_figures_v2.py — 重新绘制门水平热图
修复：① 排除细菌 patatin 背景（广谱磷脂酶）② 对数归一化让核心家族可见
      ③ 单独列示古菌 patatin（PhaZh1 型，620 基因组）
输入: results/tables/tier1_phylum_distribution.tsv（核心家族）
      results/tables/patatin_filtered.tsv（过滤后 patatin，含古菌子集）
"""
import os
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

OUT = "results/figures"
os.makedirs(OUT, exist_ok=True)

# 1. 核心家族门分布（排除 ArchPhaZ_patatin 细菌背景）
df = pd.read_csv("results/tables/tier1_phylum_distribution.tsv", sep="\t")
core = df[df["family"] != "ArchPhaZ_patatin"]
pivot = core.pivot(index="phylum", columns="family", values="genomes").fillna(0)

# 2. 古菌 patatin 门分布（从过滤后 patatin 提取，限定古菌）
TAX_AR = os.path.expanduser("~/GTDB/taxonomy/ar53_taxonomy_r232.tsv")
tax = {}
for line in open(TAX_AR):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2:
        tax[p[0]] = p[1]

def arch_phylum(g):
    t = tax.get(g) or tax.get("GB_"+g) or tax.get("RS_"+g)
    if not t:
        return None
    parts = t.split(";")
    return parts[1].replace("p__", "") if len(parts) > 1 else None

arch_pat = Counter()
for line in open("results/tables/patatin_filtered.tsv"):
    if line.startswith("genome"):
        continue
    g = line.split("\t")[0]
    ph = arch_phylum(g)
    if ph:
        arch_pat[ph] += 1

# 将古菌 patatin 作为单独列
pivot["Archaea_patatin"] = 0
for ph, n in arch_pat.items():
    pivot.loc[ph, "Archaea_patatin"] = n

# 排序：按核心家族总基因组数
pivot["_total"] = pivot.sum(axis=1)
pivot = pivot.sort_values("_total", ascending=False).drop(columns="_total")
pivot = pivot.head(25)

# 重新排列列：核心在前，古菌 patatin 最后
cols = [c for c in pivot.columns if c != "Archaea_patatin"] + ["Archaea_patatin"]
pivot = pivot[cols]

# 对数归一化热图
fig, ax = plt.subplots(figsize=(10, 11))
norm = mcolors.LogNorm(vmin=1, vmax=max(pivot.values.max(), 10))
im = ax.imshow(pivot.values, cmap="viridis", aspect="auto", norm=norm)
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels(pivot.index, fontsize=8)
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=9)
# 标注数值
for i in range(pivot.shape[0]):
    for j in range(pivot.shape[1]):
        v = int(pivot.values[i, j])
        if v > 0:
            ax.text(j, i, f"{v:,}", ha="center", va="center", fontsize=6.5,
                    color="white" if v > 1000 else "black")
cbar = plt.colorbar(im, ax=ax, label="genomes (log scale)")
ax.set_title("PHB degradation genes (tier1) across GTDB phyla\n"
             "(Archaea_patatin = PhaZh1-type, cluster-filtered)")
plt.tight_layout()
plt.savefig(f"{OUT}/tier1_phylum_heatmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("figure: tier1_phylum_heatmap.png (rebuilt)")

# 3. 家族基因组数条形图（核心家族 + 古菌 patatin）
fam_total = core.groupby("family")["genomes"].sum()
fam_total["Archaea_patatin"] = sum(arch_pat.values())
fam_total = fam_total.sort_values()
fig, ax = plt.subplots(figsize=(7, 4.5))
bars = ax.barh(fam_total.index, fam_total.values, color="#2a9d8f")
for i, v in enumerate(fam_total.values):
    ax.text(v, i, f" {int(v):,}", va="center", fontsize=9)
ax.set_xlabel("genomes (tier1)")
ax.set_xscale("log")
ax.set_title("PHB degradation gene families (tier1)")
plt.tight_layout()
plt.savefig(f"{OUT}/tier1_family_counts.png", dpi=150, bbox_inches="tight")
plt.close()
print("figure: tier1_family_counts.png (rebuilt)")

print("\n=== 数据摘要 ===")
print(fam_total.to_string())
print(f"\n古菌 patatin 门分布: {dict(arch_pat)}")
