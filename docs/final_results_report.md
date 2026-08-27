# GTDB 全库 PHB 降解基因系统分析 — 结果报告

> 日期：2026-08-16
> 服务器：T141（10.16.1.141），工作区 /home/data/haoyu/PHB_gtdb-ds
> 数据：GTDB R232 代表基因组（199,923 个，含细菌+古菌，只读）
> 方法依据：文献调研（PubMed 728 篇 + Europe PMC 全文 + OpenAlex/Web）
>
> **结论边界（重要）**：本文所有“解聚酶/降解基因”指**候选同源蛋白（功能潜力）**，
> 系 HMM 同源性 + 基序 + 邻域证据。已完成 SignalP6 序列层预测，但仍未经酶活/遗传/表型实验验证；不表述为
> “功能基因”或“实证发现”。项目状态、数据流契约与待办见 [docs/STATUS.md](STATUS.md)。

## 1. 方法摘要

```
GTDB R232 基因组（199,923 个）
  → Pyrodigal 蛋白预测（~4.7 亿蛋白）
  → 过滤超长伪影序列（>100K aa，剔除 3 条）
  → 9 家族 HMM 全库筛选（HMMER 3.4, E<1e-5）
  → 命中处理（多家族仲裁、每基因组去重）
  → 功能验证（催化三联体/lipase box 疏水x1/NAD 基序/patatin 二元组）
  → 三级重评分（curated 金标准 HMM：E<1e-20 严格 / E<1e-10 中等）
  → 门水平分布 + 系统发育
```

**家族分类体系**（PhaDED, Knoll 2009, BMC Bioinformatics 10:89）：
8 超家族 → 38 同源家族，按定位×底物×催化特征划分。

## 2. 核心结果（tier1 严格集）

### 2.1 家族计数

| 家族 | 功能 | tier1 序列 | tier1 基因组 |
|------|------|:---:|:---:|
| **ePhaZ** | 胞外 PHB/PHA 解聚酶（EC 3.1.1.75/76） | 38,692 | 27,839 |
| **iPhaZ** | 胞内 PHA 解聚酶（含 Cys 型/Ser 型/周质） | 32,926 | 25,920 |
| **OH** | 3HB 寡聚体水解酶（EC 3.1.1.22，min-cov 0.6 排除尼龙水解酶） | 1,429 | 1,410 |
| **ArchPhaZ_hydrolase** | 古菌经典 PHB 解聚酶家族酯酶 | 1,292 | 1,236 |
| **ArchPhaZ_patatin** | 古菌 patatin 样（PhaZh1 型，基因组级 PhaC 共现过滤后） | **1,372** | **620** |

**核心解聚酶 tier1 合计 ~75,650 条序列；去重后 44,814 个基因组（22.4%，不含广谱 patatin）**
**其中古菌 PHB 降解候选同源蛋白：~2,664 序列 / ~1,856 基因组（hydrolase + patatin 古菌子集，不受 OH 变化影响）**
> 注：上表"tier1 基因组"为各家族检出基因组数之和（含多家族共存重复计数）；
> 去重后 ePhaZ/iPhaZ/OH/ArchPhaZ_hydrolase 四类真解聚酶覆盖 44,814 个基因组。
> 广谱 patatin 折叠蛋白（77,424 基因组，多为磷脂酶）不计入核心解聚酶，见 §2.2。
> OH 家族按校准结论加 `--min-cov 0.6`（排除尼龙水解酶假阳性），tier1 由 1,465→1,429
> 序列、1,444→1,410 基因组；四家族去重 44,821→44,814（见 STATUS.md §9）。

**ePhaZ SignalP 细分**（SignalP6 fast，38,692 条）：**有信号肽 21,856 条（56.5%）**，
其中 Sec/SPI 分泌型 16,545（42.8%）、脂蛋白 Lipo/SPII 5,170（13.4%）、Tat 122（0.3%）、
TatLipo 19；**无信号肽（OTHER）16,836 条（43.5%）**。N 端截断分析：无信号肽组仅 2.8%
不以 Met 开头（vs 有信号肽组 0.3%）、无内部终止符、长度中位数与有信号肽组相近
（364 vs 344 aa）→ **截断不是主因**，约 43.5% 无信号肽主要是"催化域同源但非分泌"的
变体。故"胞外解聚酶"标签需审慎表述（详见 STATUS.md）。结果：
`results/signalp/ePhaZ/prediction_results.txt`、`results/tables/signalp_ePhaZ_analysis.tsv`。
**高置信胞外（有信号肽）子集已导出**：`data/screen/tiers/ePhaZ_tier1_signalpeptide.faa`
（21,856 条，header 带 SignalP 类型）+ `results/tables/ePhaZ_signalp_subset.tsv`（明细），
并按类型拆分为 `ePhaZ_tier1_signalpeptide_{SP,LIPO,TAT,TATLIPO}.faa`。
**门分布**（有信号肽 21,856 条）：Actinomycetota 30.3%、Pseudomonadota 25.5%、
Bacteroidota 15.3%；SP-only（16,545 条）与 LIPO-only（5,170 条）分别见
`results/tables/ePhaZ_SP_phylum.tsv`、`ePhaZ_LIPO_phylum.tsv`（LIPO 在 Myxococcota_A
11.3%、Chloroflexota 10.2% 相对富集）。

### 2.2 古菌 patatin 型解聚酶分布（PhaZh1 型候选，经基因簇复筛）

| 古菌门/纲 | 基因组 | 意义 |
|-----------|:---:|------|
| **Thermoproteota / Nitrososphaeria** | 323 | 氨氧化古菌（AOA）——**超出文献已知范围的新发现** |
| Halobacteriota / Halobacteria | 266 | 嗜盐古菌（Haloferax 等，文献已表征）|
| Halobacteriota / Methanosarcinia | 23 | 产甲烷古菌 |
| Thermoplasmatota | 8 | 嗜热酸古菌 |

> 文献仅实验表征了 *Haloferax mediterranei* 的 PhaZh1（PMID 25710370）；
> 本研究首次在 GTDB 全库尺度**检出古菌 patatin 候选同源蛋白**扩展至
> Nitrososphaeria（AOA）、Methanosarcinia、Thermoplasmata（**功能潜力，非实证**）。

#### 基因簇共定位复筛（±10kb，11_clusters）

上述 1,372 条 patatin 命中的历史 ±10kb 基因邻域记录使用 PhaC / PhaE / phasin /
BdhA / PhaJ / PHA_gran_rgn 作为标记；该记录的审计状态如下：

- 历史报告记载 340 条（24.8%）邻近 PHB 代谢标记，但现存 `cluster_summary.tsv` 为仅含
  `cooccurring_loci` 的旧 schema，且与历史叙述的计数不一致。该数字仅保留为历史记录；须用最终
  `hits_filtered.tsv` 重跑位点级分析并输出新 schema 后，才能给出唯一位点或唯一基因组结论。

> 判据依据：PhaZh1 与 **bdhA** 成簇（Liu 2015, PMID 25710370），故降解支路
> （BdhA/PhaJ）邻近比合成簇（PhaC/phasin）邻近更能指示真 PhaZh1 型。

#### 生物学 caveat

> patatin 是广谱"非经典脂质水解酶"结构域（兼具脂酶/酯酶/磷脂酶活性），
> 多数命中并非 PHB 解聚酶。且 PhaZh1 型解聚酶**体内角色有限**（敲除不影响动员，
> 存在替代通路），古菌 PHB 动员的**主通路是 PhaJ（烯酰-CoA 水合酶）**。
> 因此古菌 patatin 型"解聚酶"的实际体内贡献应持审慎表述——它是颗粒结合、
> 体外高效的蛋白，但并非古菌 PHB 降解的主要执行者。

### 2.3 门水平分布（tier1 核心解聚酶，去重后，不含广谱 patatin）

> 下表为 ePhaZ/iPhaZ/OH/ArchPhaZ_hydrolase 四类真解聚酶的去重门分布（共 44,814 个
> 基因组）。广谱 patatin 折叠蛋白（77,424 基因组）单独处理，见 §2.2。

| 门 | 基因组 | 主要家族 |
|----|:---:|------|
| Pseudomonadota | 26,850 | iPhaZ(23,559) + ePhaZ(11,995)；以 Figure 3 source data 为准 |
| Actinomycetota | 6,810 | ePhaZ(5,657) |
| Bacteroidota | 3,247 | ePhaZ(3,117) |
| Chloroflexota | 1,200 | ePhaZ(1,155) |
| Acidobacteriota | 867 | ePhaZ(798) |
| Planctomycetota | 817 | ePhaZ(792) |
| Myxococcota_A | 701 | ePhaZ(623) |
| Desulfobacterota | 585 | ePhaZ(437) |
| Gemmatimonadota | 565 | ePhaZ(544) |
| Bacillota | 434 | ePhaZ(383) |
| 其他细菌门 | ~2,440 | 广泛 |
| **Thermoplasmatota（古菌）** | 154 | hydrolase(127) + ePhaZ(35) |
| **Halobacteriota（古菌）** | 112 | ePhaZ(64) + hydrolase(46) |
| **Thermoproteota（古菌）** | 32 | ePhaZ(19) + hydrolase(13) |

## 3. 与文献对照

| 指标 | 本次 | 文献基准 |
|------|------|---------|
| 解聚酶总数 | ~75,650 (tier1 序列)；44,814 基因组(去重) | Viljakainen & Hug 2021：13,869（3078 宏基因组）|
| 主导门 | Pseudomonadota(26,850 基因组) | 一致（V&H 2021：Proteobacteria+Bacteroidota 主导）|
| 古菌检出 | 核心解聚酶 Halobacteriota 112 + Thermoproteota 32 + Thermoplasmatota 154；patatin 古菌子集 620 | 仅 Hfx. mediterranei 2 酶实验表征；Sulfolobus 旁证 |
| 旧项目（14 refs） | — | 6,532（严重低估，印证本次扩展）|

## 4. 关键发现

1. **Pseudomonadota 主导**：胞内型 iPhaZ 在变形菌中特别富集（23,559 基因组），与 PHA 累积菌（Cupriavidus、Pseudomonas、Ralstonia 等均为变形菌）一致。
2. **古菌谱系检出 PHB 降解候选同源蛋白（功能潜力）**：核心解聚酶（ePhaZ/hydrolase）候选在 Halobacteriota 112 基因组、Thermoproteota 32 基因组、Thermoplasmatota 154 基因组检出；另有 patatin 折叠蛋白古菌子集 620 基因组，其局部邻域支持计数仍待最终方案 A 输入与新 schema 重跑。该分布仅代表候选同源蛋白，未经酶活、遗传或表型验证。
3. **patatin 家族需基因簇过滤**：77,424 基因组检出 patatin 折叠蛋白，但它是广谱磷脂酶结构域；
  古菌子集的邻域支持计数仍待按最终输入及新 schema 重跑。PhaZh1 体内角色有限、PhaJ 才是动员
  主路，因此 patatin 型"解聚酶"只能审慎表述为候选。
4. **多家族共存**：30,062 基因组同时含 ≥2 个核心家族（如 ePhaZ+iPhaZ 或 iPhaZ+OH），反映完整的降解-动员通路共现。

## 5. 数据文件

- 命中汇总：data/screen/hits_all.tsv（676 万行）
- tier1 序列：data/screen/tiers/{fam}_tier1.faa
- 基因组×家族：results/tables/tier1_genome_family.tsv
- 门分布：results/tables/tier1_phylum_distribution.tsv
- 图：results/figures/tier1_phylum_heatmap.png, tier1_family_counts.png
- 系统发育树：results/trees_tier1/*.treefile（生成中）

## 6. 局限与后续

- patatin 的历史 PhaC 共现和位点级邻域结果保留供追溯；必须以方案 A 最终输入重跑
  `11_clusters.py` 后，才可作为当前结论。
- 生态元数据表同样属于历史下游输出，需用方案 A 最终输入重新执行 `10_distribution.py` 后再发布。
- 系统发育树：ArchPhaZ_hydrolase 的 IQ-TREE2 树为当前输入；OH 的 1,465 叶树输入已过期，
  必须以 1,429 条 OH tier1 重建。ePhaZ/iPhaZ 两个大族（38,692/32,926 序列）的全量树仍暂缓。
- BdhA/PhaJ 等广谱代谢酶已标注为背景家族，不计入核心解聚酶计数
- 建议后续：tier1 集的 SignalP 胞外/胞内细分
- **PhaZh1 专属 HMM 重建（已尝试，结论：不可行）**：以实验表征的 PhaZh1（I3RBH0）为探针 DIAMOND
  找到 361 条同源序列建 HMM，但重筛古菌 patatin 仍匹配 93.4%（1,282/1,372）——因 patatin 催化结构域
  序列高度保守，同源法无法区分颗粒结合型解聚酶与广谱磷脂酶。故 PhaZh1 型特异性应以**基因簇共定位**
  （以最终 `hits_filtered.tsv` 重跑的 ±10kb 邻域结果）为判据，而非序列 HMM。
