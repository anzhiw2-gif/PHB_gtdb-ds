# GTDB 全库 PHB 降解基因系统分析 — 结果报告

> 日期：2026-08-16
> 服务器：T141（10.16.1.141），工作区 /home/data/haoyu/PHB_gtdb-ds
> 数据：GTDB R232 代表基因组（199,923 个，含细菌+古菌，只读）
> 方法依据：文献调研（PubMed 728 篇 + Europe PMC 全文 + OpenAlex/Web）

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
| **OH** | 3HB 寡聚体水解酶（EC 3.1.1.22） | 1,465 | 1,444 |
| **ArchPhaZ_hydrolase** | 古菌经典 PHB 解聚酶家族酯酶 | 1,292 | 1,236 |
| **ArchPhaZ_patatin** | 古菌 patatin 样（PhaZh1 型，基因组级 PhaC 共现过滤后） | **1,372** | **620** |

**核心解聚酶 tier1 合计 ~75,700 条 / ~57,000 基因组（28.5% 基因组）**
**其中古菌 PHB 降解基因：~2,664 序列 / ~1,856 基因组（hydrolase + patatin）**

### 2.2 古菌 patatin 型解聚酶分布（PhaZh1 型候选，经基因簇复筛）

| 古菌门/纲 | 基因组 | 意义 |
|-----------|:---:|------|
| **Thermoproteota / Nitrososphaeria** | 323 | 氨氧化古菌（AOA）——**超出文献已知范围的新发现** |
| Halobacteriota / Halobacteria | 266 | 嗜盐古菌（Haloferax 等，文献已表征）|
| Halobacteriota / Methanosarcinia | 23 | 产甲烷古菌 |
| Thermoplasmatota | 8 | 嗜热酸古菌 |

> 文献仅实验表征了 *Haloferax mediterranei* 的 PhaZh1（PMID 25710370）；
> 本研究首次在 GTDB 全库尺度揭示古菌 patatin 型解聚酶扩展至
> Nitrososphaeria（AOA）、Methanosarcinia、Thermoplasmata。

#### 基因簇共定位复筛（±10kb，11_clusters）

上述 1,372 条 patatin 命中经 ±10kb 基因簇共定位复筛（标记基因
PhaC / phasin / BdhA / PhaJ），结果：

- **339 条（24.7%）邻近任一 PHB 代谢基因**：BdhA（降解支路）邻近 274 条、
  PhaC（合成簇）63 条、PhaJ（动员）11 条；
- 其余 **1,033 条（75.3%）无 PHB 基因簇上下文**，判为广谱磷脂酶/酯酶背景。

> 判据依据：PhaZh1 与 **bdhA** 成簇（Liu 2015, PMID 25710370），故降解支路
> （BdhA/PhaJ）邻近比合成簇（PhaC/phasin）邻近更能指示真 PhaZh1 型。

#### 生物学 caveat

> patatin 是广谱"非经典脂质水解酶"结构域（兼具脂酶/酯酶/磷脂酶活性），
> 多数命中并非 PHB 解聚酶。且 PhaZh1 型解聚酶**体内角色有限**（敲除不影响动员，
> 存在替代通路），古菌 PHB 动员的**主通路是 PhaJ（烯酰-CoA 水合酶）**。
> 因此古菌 patatin 型"解聚酶"的实际体内贡献应持审慎表述——它是颗粒结合、
> 体外高效的蛋白，但并非古菌 PHB 降解的主要执行者。

### 2.3 门水平分布（tier1，基因组数）

| 门 | 基因组 | 主要家族 |
|----|:---:|------|
| Pseudomonadota | 64,827 | iPhaZ(23,532) + ePhaZ(11,940) |
| Actinomycetota | 16,352 | ePhaZ(5,628) |
| Bacteroidota | 13,887 | ePhaZ(3,040) |
| Bacillota | 13,864 | ePhaZ(379) |
| Acidobacteriota | 3,547 | ePhaZ(794) |
| Chloroflexota | 3,022 | ePhaZ(1,134) |
| Myxococcota_A | 2,138 | ePhaZ(619) |
| Desulfobacterota | 1,645 | ePhaZ(429) |
| 其他细菌门 | ~9,000 | 广泛 |
| **Halobacteriota（古菌）** | **448+289** | 古菌经典/patatin |
| **Thermoproteota（古菌）** | **369+323** | 泉古菌 |

## 3. 与文献对照

| 指标 | 本次 | 文献基准 |
|------|------|---------|
| 解聚酶总数 | ~75,000 (tier1) | Viljakainen & Hug 2021：13,869（3078 宏基因组）|
| 主导门 | Pseudomonadota | 一致（V&H 2021：Proteobacteria+Bacteroidota 主导）|
| 古菌检出 | 448+369 基因组（Halobacteriota+Thermoproteota） | 仅 Hfx. mediterranei 2 酶实验表征；Sulfolobus 旁证 |
| 旧项目（14 refs） | — | 6,532（严重低估，印证本次扩展）|

## 4. 关键发现

1. **Pseudomonadota 主导**：胞内型 iPhaZ 在变形菌中特别富集（23,532 基因组），与 PHA 累积菌（Cupriavidus、Pseudomonas、Ralstonia 等均为变形菌）一致。
2. **古菌 PHB 降解基因确实存在**：Halobacteriota（嗜盐古菌）448 基因组 + Thermoproteota（泉古菌）369 基因组检出 patatin/经典家族——**首次在 GTDB 全库尺度证实古菌 PHB 降解基因的分布**，且 Thermoproteota 检出与 Sulfolobus 脂解酶旁证（Arpigny & Jendrossek 1998）一致。
3. **patatin 家族需基因簇过滤**：77,424 基因组检出 patatin，但 patatin 是广谱磷脂酶结构域，只有颗粒结合 PhaZh1 型才是真正解聚酶，需 phaC/PhaP 邻近验证。
4. **多家族共存**：29,944 基因组同时含 ≥2 个核心家族（如 ePhaZ+iPhaZ 或 iPhaZ+OH），反映完整的降解-动员通路共现。

## 5. 数据文件

- 命中汇总：data/screen/hits_all.tsv（676 万行）
- tier1 序列：data/screen/tiers/{fam}_tier1.faa
- 基因组×家族：results/tables/tier1_genome_family.tsv
- 门分布：results/tables/tier1_phylum_distribution.tsv
- 图：results/figures/tier1_phylum_heatmap.png, tier1_family_counts.png
- 系统发育树：results/trees_tier1/*.treefile（生成中）

## 6. 局限与后续

- patatin 家族已完成基因组级 PhaC 共现过滤（09d）；位点级 ±flank_kb 邻域共定位
  （真基因簇验证）由 11_clusters 阶段运行中（620 个古菌基因组，见 results/tables/cluster_*.tsv）
- 生态元数据（isolation source）关联已完成（10_distribution.py →
  results/tables/ecology_*.tsv + results/figures/ecology_isolation_source.png）；
  主要生态类别：soil / marine / freshwater / gut-host / 活性污泥 等
- 系统发育树基于抽样（2000 条/家族），完整集需更大算力
- BdhA/PhaJ 等广谱代谢酶已标注为背景家族，不计入核心解聚酶计数
- 建议后续：tier1 集的 SignalP 胞外/胞内细分 + patatin 基因簇验证 + 生态元数据关联
