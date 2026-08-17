#!/usr/bin/env python3
"""09a_tier1_summary.py — 从 tier1 序列生成基因组×家族统计与门水平分布
输入: data/screen/tiers/{fam}_tier1.faa
      ~/GTDB/taxonomy/bac120_taxonomy_r232.tsv (只读)
      ~/GTDB/taxonomy/ar53_taxonomy_r232.tsv (只读)
输出: results/tables/tier1_genome_family.tsv
      results/tables/tier1_phylum_distribution.tsv
"""
import os
from collections import Counter, defaultdict

TIER_DIR = "data/screen/tiers"
TAX_BAC = os.path.expanduser("~/GTDB/taxonomy/bac120_taxonomy_r232.tsv")
TAX_AR = os.path.expanduser("~/GTDB/taxonomy/ar53_taxonomy_r232.tsv")
OUT = "results/tables"
os.makedirs(OUT, exist_ok=True)

# 读 GTDB 分类（键带前缀 GB_GCA_/RS_GCF_，我们的 ID 无前缀）
tax = {}
for tf in [TAX_BAC, TAX_AR]:
    if os.path.exists(tf):
        for line in open(tf):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                tax[p[0]] = p[1]


def lookup_tax(g: str) -> str:
    """我们的 ID（GCA_xxx/GCF_xxx）→ GTDB 分类"""
    if g in tax:
        return tax[g]
    if g.startswith("GCA_"):
        return tax.get("GB_" + g, "unclassified")
    if g.startswith("GCF_"):
        return tax.get("RS_" + g, "unclassified")
    return "unclassified"

FAMILIES = ["ePhaZ", "iPhaZ", "OH", "ArchPhaZ_hydrolase", "ArchPhaZ_patatin"]

genome_fam = Counter()   # (genome, family) -> 拷贝数
genome_set = defaultdict(set)  # genome -> set(families)
for fam in FAMILIES:
    faa = os.path.join(TIER_DIR, f"{fam}_tier1.faa")
    if not os.path.exists(faa):
        continue
    cnt = Counter()
    for line in open(faa):
        if line.startswith(">"):
            g = line[1:].split("|")[0]
            cnt[g] += 1
            genome_set[g].add(fam)
    for g, c in cnt.items():
        genome_fam[(g, fam)] = c
    print(f"{fam}: {len(cnt)} 基因组, {sum(cnt.values())} 序列")

# 写 genome-family 表
with open(f"{OUT}/tier1_genome_family.tsv", "w") as f:
    f.write("genome\tfamily\tcopies\tgtdb_taxonomy\tphylum\tclass\n")
    for (g, fam), c in sorted(genome_fam.items()):
        t = lookup_tax(g)
        parts = t.split(";")
        phylum = parts[1].replace("p__", "") if len(parts) > 1 else "unknown"
        cls = parts[2].replace("c__", "") if len(parts) > 2 else "unknown"
        f.write(f"{g}\t{fam}\t{c}\t{t}\t{phylum}\t{cls}\n")

# 门水平分布
phylum_fam = Counter()
for (g, fam), c in genome_fam.items():
    t = lookup_tax(g)
    parts = t.split(";")
    phylum = parts[1].replace("p__", "") if len(parts) > 1 else "unknown"
    phylum_fam[(phylum, fam)] += 1

with open(f"{OUT}/tier1_phylum_distribution.tsv", "w") as f:
    f.write("phylum\tfamily\tgenomes\n")
    for (ph, fam), n in sorted(phylum_fam.items(), key=lambda x: -x[1]):
        f.write(f"{ph}\t{fam}\t{n}\n")

# 打印摘要
print("\n=== 门水平分布（tier1，基因组数）===")
phylum_total = Counter()
for (ph, fam), n in phylum_fam.items():
    phylum_total[ph] += n
for ph, n in phylum_total.most_common(20):
    fams = {fam: c for (p, fam), c in phylum_fam.items() if p == ph}
    fam_str = ", ".join(f"{f}:{c}" for f, c in sorted(fams.items(), key=lambda x: -x[1])[:4])
    print(f"  {ph}: {n} 基因组 [{fam_str}]")

print(f"\n总基因组数(任一核心家族 tier1): {len(genome_set)}")
multi = sum(1 for g, fs in genome_set.items() if len(fs) >= 2)
print(f"多家族共存基因组: {multi}")
