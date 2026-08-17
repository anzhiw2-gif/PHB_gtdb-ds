#!/usr/bin/env python3
"""09c_tier1_figures.py — tier1 分布图（门水平热图 + 家族条形图）"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = "results/figures"
os.makedirs(OUT, exist_ok=True)
df = pd.read_csv("results/tables/tier1_phylum_distribution.tsv", sep="\t")

# 门水平热图（top 20 门 × 家族）
pivot = df.pivot(index="phylum", columns="family", values="genomes").fillna(0)
pivot["total"] = pivot.sum(axis=1)
pivot = pivot.sort_values("total", ascending=False).drop(columns="total").head(20)

fig, ax = plt.subplots(figsize=(9, 10))
im = ax.imshow(pivot.values, cmap="YlGnBu", aspect="auto")
ax.set_yticks(range(len(pivot.index)), pivot.index, fontsize=9)
ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=45, ha="right")
for i in range(pivot.shape[0]):
    for j in range(pivot.shape[1]):
        v = int(pivot.values[i, j])
        if v > 0:
            ax.text(j, i, str(v), ha="center", va="center", fontsize=7)
plt.colorbar(im, label="genomes")
ax.set_title("PHB degradation genes (tier1) across GTDB phyla")
plt.tight_layout()
plt.savefig(f"{OUT}/tier1_phylum_heatmap.png", dpi=150)
plt.close()
print("figure: tier1_phylum_heatmap.png")

# 家族基因组数条形图（不含 patatin，核心解聚酶）
core = df[~df["family"].str.contains("patatin")]
fam_count = core.groupby("family")["genomes"].sum().sort_values()
fig, ax = plt.subplots(figsize=(7, 4))
fam_count.plot(kind="barh", ax=ax, color="#2a9d8f")
for i, v in enumerate(fam_count):
    ax.text(v, i, f" {int(v):,}", va="center", fontsize=9)
ax.set_xlabel("genomes (tier1)")
ax.set_title("Core PHB degradation gene families (tier1)")
plt.tight_layout()
plt.savefig(f"{OUT}/tier1_family_counts.png", dpi=150)
plt.close()
print("figure: tier1_family_counts.png")

# patatin 单独标注（背景家族）
pat = df[df["family"].str.contains("patatin")]
print(f"\npatatin 背景家族: {pat['genomes'].sum():,} 基因组（需基因簇过滤）")
