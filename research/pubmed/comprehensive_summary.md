# GTDB 全库 PHB 降解基因系统筛选 —— PubMed 穷尽式文献调研报告

> 检索工具：`pubmed-database` skill（NCBI E-utilities，经 `scripts/pubmed_api.py` 包装，遵守限速规则）
> 检索日期：2026-xx（本会话）
> 检索范围：PubMed 全库（无年代限制，含 1953–2026 年文献）
> 数据文件：`comp_*.json`（80 个查询原始 PMID）、`comp_abstracts_full.json`（748 篇元数据+摘要）、`comp_slim_all.json`（slim 版）
> 总计：**80 个查询**，去重后 **748 篇唯一 PMID**，其中 **729 篇有摘要**

---

## 1. 检索概况表

### A. 家族分类体系（8+ 查询）

| 查询 | 命中 | 文件 |
|---|---|---|
| polyhydroxybutyrate depolymerase classification | 15 | comp_a01 |
| PHA depolymerase family | 12 | comp_a02 |
| polyhydroxyalkanoate depolymerase catalytic domain classification | 1 | comp_a03 |
| Jendrossek polyhydroxyalkanoate depolymerase types | 1 | comp_a04 |
| Jendrossek PHA depolymerase（变体） | 14 | comp_a09 |
| PHA depolymerase substrate-binding domain | 9 | comp_a05 |
| intracellular polyhydroxyalkanoate depolymerase | 18 | comp_a06 |
| polyhydroxybutyrate oligomer hydrolase | 23 | comp_a07 |
| PHA depolymerase database OR bioinformatics | 25 | comp_a08 |
| depolymerase catalytic domain lipase box | 9 | comp_a10 |
| polyhydroxyalkanoate depolymerase classification（变体） | 6 | comp_a11 |
| extracellular polyhydroxybutyrate depolymerase | 30 | comp_a12 |
| polyhydroxybutyrate depolymerase type I type II classification | 0 | comp_t01 |
| PHA depolymerase types classification（变体） | 2 | comp_t06 |
| polyhydroxybutyrate depolymerase domain architecture | 1 | comp_t08 |
| PHB depolymerase catalytic domain type linker substrate-binding domain classification | 2 | comp_u04 |

### B. 细菌已表征酶（按属，10+ 查询）

| 查询 | 命中 | 文件 |
|---|---|---|
| Pseudomonas lemoignei polyhydroxybutyrate depolymerase | 12 | comp_b01 |
| Pseudomonas lemoignei PhaZ | 1 | comp_s02 |
| Alcaligenes faecalis PHB depolymerase | 26 | comp_b02 |
| Comamonas polyhydroxybutyrate depolymerase | 8 | comp_b03 |
| Acidovorax polyhydroxybutyrate depolymerase | 30 | comp_b04 |
| Bacillus polyhydroxybutyrate depolymerase | 18 | comp_b05 |
| Streptomyces PHA depolymerase | 30 | comp_b06 |
| Thermus polyhydroxybutyrate depolymerase | 30 | comp_b07 |
| Paucimonas polyhydroxybutyrate depolymerase | 30 | comp_b08 |
| PhaZ7 Paucimonas | 10 | comp_s01 |
| Cupriavidus necator PHB degradation | 30 | comp_s03 |
| Ralstonia eutropha PHB depolymerase | 26 | comp_s04 |
| Pseudomonas putida PHA depolymerase | 29 | comp_s05 |
| Rhodococcus PHA depolymerase | 2 | comp_s06 |
| Burkholderia PHA depolymerase | 4 | comp_s09 |
| Aeromonas PHA depolymerase | 30 | comp_s10 |
| thermophilic PHB depolymerase | 5 | comp_s07 |
| marine PHA depolymerase | 30 | comp_s08 |
| PhaZ1-7 编号 | 30 | comp_s13 |
| PhaZ family polyhydroxybutyrate | 6 | comp_s14 |
| PHB depolymerase serine hydrolase catalytic triad | 9 | comp_s15 |
| short-chain-length PHA depolymerase | 7 | comp_s22 |
| medium-chain-length PHA depolymerase | 23 | comp_s23 |
| PhaDED / PhaZ extracellular PHB | 13 | comp_s24 |
| PhaZCn / PhaZ1 Cupriavidus | 8 | comp_s26 |

### C. 古菌（10+ 查询）

| 查询 | 命中 | 文件 |
|---|---|---|
| haloarchaea polyhydroxyalkanoate degradation | 30 | comp_c01 |
| Haloferax polyhydroxybutyrate depolymerase OR degradation | 30 | comp_c02 |
| archaea polyhydroxyalkanoate depolymerase | 1 | comp_c03 |
| Natronomonas OR Halobacterium OR Haloarcula PHA degradation | 30 | comp_c04 |
| archaeal polyhydroxyalkanoate degradation | 30 | comp_c05 |
| archaea PHB degradation OR PHB depolymerase | 30 | comp_c06 |
| halophilic archaea polyhydroxyalkanoate | 30 | comp_c07 |
| Haloferax mediterranei PHB | 30 | comp_s19 |
| archaea polyhydroxybutyrate granules degradation | 10 | comp_s20 |
| archaea polyhydroxyalkanoate mobilization OR degradation genes | 30 | comp_s28 |
| archaeal PhaZ / haloarchaeal PHB depolymerase | 0 | comp_s16 |
| Haloferax mediterranei PHA degradation enzymes genes | 4 | comp_t02 |
| halophilic archaea PHA depolymerase gene phaZ | 0 | comp_t03 |
| archaea PHB granules degradation enzymes | 4 | comp_t07 |
| patatin-like depolymerase PHA | 1 | comp_t09 |
| PhaZh1 Haloferax | 1 | comp_u01 |
| Haloferax mediterranei depolymerase OR mobilization | 10 | comp_u03 |

### D. 3HB 代谢与基因簇（5+ 查询）

| 查询 | 命中 | 文件 |
|---|---|---|
| 3-hydroxybutyrate dehydrogenase bdhA | 14 | comp_d01 |
| PHA mobilization enoyl-CoA hydratase | 1 | comp_d02 |
| enoyl-CoA hydratase PHB degradation | 12 | comp_d04 |
| 3-hydroxybutyrate dehydrogenase polyhydroxybutyrate | 30 | comp_d05 |
| polyhydroxybutyrate degradation gene cluster | 30 | comp_d03 |
| oligomer hydrolase intracellular PHB | 15 | comp_s27 |
| 3-hydroxybutyrate oligomer hydrolase | 25 | comp_t05 |
| PHB depolymerase substrate binding domain linker catalytic domain | 11 | comp_t04 |
| PhaZ depolymerase regulation transcription | 5 | comp_u05 |

### E. 综述与方法（6+ 查询）

| 查询 | 命中 | 文件 |
|---|---|---|
| polyhydroxyalkanoate degradation review | 30 | comp_e01 |
| polyhydroxyalkanoate degradation | 30 | comp_e04 |
| PHA depolymerase metagenome OR genome screening | 30（含噪声） | comp_e02 |
| (PHA depolymerase OR PHB depolymerase) AND (metagenome OR genome mining...) | 13 | comp_e02b |
| polyhydroxybutyrate depolymerase structure | 30 | comp_e03 |
| PHB depolymerase crystal structure | 21 | comp_e05 |
| phaZ polyhydroxyalkanoate depolymerase gene | 16 | comp_e06 |
| Penicillium / fungal PHB depolymerase | 30 | comp_s25 |
| PHB depolymerase biodegradation plastic | 23 | comp_s17 |
| Jendrossek[Author] polyhydroxyalkanoate depolymerase | 11 | comp_u02 |
| Jendrossek microbial degradation polyhydroxyalkanoates | 2 | comp_s18 |
| PhaZ depolymerase review | 0 | comp_s12 |

**说明**：个别宽泛查询（comp_c02、comp_c04、comp_e02、comp_s28）被 PROTAC/癌症/蛋白降解类文献污染（"degradation"一词），已通过变体查询 + 人工甄别剔除；表中均标注。

---

## 2. PHA 解聚酶家族分类体系综述

### 2.1 总体框架（胞外 vs 胞内）

PHA 解聚酶（EC 3.1.1.75 PHB depolymerase / EC 3.1.1.76 PHA depolymerase）是一类丝氨酸酯酶，共享 **α/β-hydrolase 折叠** 与 **催化三联体（Ser-His-Asp）**，但序列与底物特异性极为多样 [[19296857]](https://pubmed.ncbi.nlm.nih.gov/19296857/)。

- **胞外 PHA 解聚酶（e-PHA depolymerase，EC 3.1.1.75/76）**：由 PHA 降解菌分泌，将胞外聚合物水解为可吸收的（R）-3-羟基脂肪酸/寡聚体 [[12213937]](https://pubmed.ncbi.nlm.nih.gov/12213937/)、[[9008883]](https://pubmed.ncbi.nlm.nih.gov/9008883/)。
- **胞内 PHA 解聚酶（i-PHA depolymerase / 动员酶）**：PHA 积累菌自身颗粒上的降解系统，参与储存聚合物动员 [[12213937]](https://pubmed.ncbi.nlm.nih.gov/12213937/)、[[11004196]](https://pubmed.ncbi.nlm.nih.gov/11004196/)。

### 2.2 PhaDED 数据库（PHA Depolymerase Engineering Database）—— 权威家族分类

- 2009 年建立，收录 **587 条 PHA 解聚酶序列**，基于序列相似性划分为 **8 个超家族（superfamilies）和 38 个同源家族（homologous families）**；提供多序列比对、profile HMM 与功能位点注释 [[19296857]](https://pubmed.ncbi.nlm.nih.gov/19296857/)。
- 该数据库是"从基因组中 in silico 鉴定新 PHA 解聚酶、分类、预测生化性质、设计酶变体"的标准工具（R. rubrum PhaZ3 即由该库预测后被实验验证为胞内酶 [[21274528]](https://pubmed.ncbi.nlm.nih.gov/21274528/)）。
- 应用实例：在 silico metagenome 挖掘酯酶（含 PHA 解聚酶同源物）[[25502823]](https://pubmed.ncbi.nlm.nih.gov/25502823/)；Isabel 岛苏打湖宏基因组中鉴定 16 个假定的 PHB 解聚酶 [[41702408]](https://pubmed.ncbi.nlm.nih.gov/41702408/)。

### 2.3 胞外解聚酶的域结构与分类

典型 e-PHA 解聚酶为多域结构：**信号肽 — 催化域（CD）— 连接域（linker）— 底物结合域（SBD）** [[7836292]](https://pubmed.ncbi.nlm.nih.gov/7836292/)、[[9371441]](https://pubmed.ncbi.nlm.nih.gov/9371441/)、[[7606660]](https://pubmed.ncbi.nlm.nih.gov/7606660/)。

**(a) 催化域（CD）分类 —— lipase box 位置（type I/II）**
- 依据 **催化域内脂酶盒（lipase box，G-X-S-X-G / PHB depolymerase box G-L-S-X-G）的位置** 划分 type I 与 type II：Ralstonia pickettii T1 酶为 type I，Acidovorax sp. TP4 酶为 type II [[12099829]](https://pubmed.ncbi.nlm.nih.gov/12099829/)。
- 依据底物结合域与催化域的同源性，A. faecalis（PhaZ Afa）、P. stutzeri（PhaZ Pst）、C. acidovorans（PhaZ Cac）曾被归为 **type A/B**（按 lipase box 在催化域中的位置）[[10408639]](https://pubmed.ncbi.nlm.nih.gov/10408639/)。
- 催化域另有 type 1 / type 2 之分的描述：Schlegelella sp. KB1a 酶 = catalytic domain **type 2** – linker Fn3 – SBD **type 1** [[15340791]](https://pubmed.ncbi.nlm.nih.gov/15340791/)；Bacillus sp. NRRL B-14911 酶催化域为 type 1 [[21948827]](https://pubmed.ncbi.nlm.nih.gov/21948827/)。
- 催化三联体实例：A. faecalis T1 **Ser139-Asp214-His273** [[9371441]](https://pubmed.ncbi.nlm.nih.gov/9371441/)、[[8764515]](https://pubmed.ncbi.nlm.nih.gov/8764515/)；P. lemoignei PhaZ5 **Ser138**（G-L-S-S-G）[[8764515]](https://pubmed.ncbi.nlm.nih.gov/8764515/)；Thermus thermophilus HB8 **Ser183-Glu310-His405**（G-X-S-X-G）[[19214501]](https://pubmed.ncbi.nlm.nih.gov/19214501/)；Schlegelella KB1a **Ser20-Asp104-His138**（G-L-S-A-G）[[15340791]](https://pubmed.ncbi.nlm.nih.gov/15340791/)。

**(b) 连接域（linker）类型**
- 常见类型：**fibronectin type III（Fn3）模块** [[9406404]](https://pubmed.ncbi.nlm.nih.gov/9406404/)、[[7606660]](https://pubmed.ncbi.nlm.nih.gov/7606660/)、[[15340791]](https://pubmed.ncbi.nlm.nih.gov/15340791/)；**cadherin-like 重复**（P. stutzeri [[9872779]](https://pubmed.ncbi.nlm.nih.gov/9872779/)、Marinobacter sp. NK-1 [[14607367]](https://pubmed.ncbi.nlm.nih.gov/14607367/)）；海洋菌常见 **Fn3 linker**（Alteromonas [[40500476]](https://pubmed.ncbi.nlm.nih.gov/40500476/)）。
- **新型 linker**：Bacillus sp. NRRL B-14911 酶具 **LD1、LD2 两种新型 linker 域**，与已知三类 linker 均无同源性 [[21948827]](https://pubmed.ncbi.nlm.nih.gov/21948827/)。

**(c) 底物结合域（SBD）类型**
- 经典研究确立 SBD 位于 C 端、对 PHB 颗粒特异吸附（P. lemoignei PhaZ4 C 端 55 aa 缺失即丧失结合 [[8837471]](https://pubmed.ncbi.nlm.nih.gov/8837471/)）。
- **已知 2 类 dPHB 结合域**；Bacillus sp. NRRL B-14911 发现 **SBD1/SBD2 两种新型结合域**（与已知类型无同源性，且彼此无同源性）[[21948827]](https://pubmed.ncbi.nlm.nih.gov/21948827/)。
- SBD 突变功能研究（R. pickettii T1）：Ser/Tyr/Val/Ala/Leu 残基对吸附与降解至关重要，SBD 通过氢键（Ser-OH…羰基）与疏水作用结合 PHB 表面 [[16963553]](https://pubmed.ncbi.nlm.nih.gov/16963553/)、[[20058938]](https://pubmed.ncbi.nlm.nih.gov/20058938/)；吸附符合 Langmuir 等温线，单分子占用面积约 12±4 nm² [[10408639]](https://pubmed.ncbi.nlm.nih.gov/10408639/)。
- SBD 与 linker 删除突变体失去不溶性 PHB 水解能力，但保留水溶性寡聚体水解活性 → 催化域功能独立于 SBD/linker（P. stutzeri [[11710119]](https://pubmed.ncbi.nlm.nih.gov/11710119/)）。

**(d) 无 SBD 的特殊胞外酶（PhaZ7 亚群）**
- Paucimonas lemoignei PhaZ7：**只水解天然（无定形）PHB（nPHB），不能水解结晶态 dPHB**，无 SBD，单域 α/β-hydrolase 折叠；最初被误认为"PHB depolymerase inhibitor（PDI）" [[12099828]](https://pubmed.ncbi.nlm.nih.gov/12099828/)、[[12855176]](https://pubmed.ncbi.nlm.nih.gov/12855176/)。催化三联体 **His47-Ser136-Asp242** + 氧阴离子洞 His306 [[12855176]](https://pubmed.ncbi.nlm.nih.gov/12855176/)；结构上活性位点埋藏，需构象开关（loop 281–295 类似 lid）打开通道 [[18706425]](https://pubmed.ncbi.nlm.nih.gov/18706425/)、[[28370478]](https://pubmed.ncbi.nlm.nih.gov/28370478/)；1.2 Å 原子分辨率结构（迄今最高分辨率解聚酶结构）[[20516591]](https://pubmed.ncbi.nlm.nih.gov/20516591/)；底物结合位点 Y105/Y176/Y189/Y190 由 mutein 结构确定 [[24007310]](https://pubmed.ncbi.nlm.nih.gov/24007310/)。
- 真核：Penicillium funiculosum PHB depolymerase 为**环状排列（circularly permuted）的 α/β-hydrolase**，无独立 SBD，靠表面 13 个疏水残基吸附 [[16405909]](https://pubmed.ncbi.nlm.nih.gov/16405909/)。

### 2.4 胞内解聚酶（i-PHA depolymerase）分类

**模式菌 Ralstonia (Cupriavidus) eutropha H16**：基因组含 **9 个胞内动员相关基因，分 4 类**（PhaZa1–a5、PhaZb、PhaZc、PhaZd1/PhaZd2）[[18156336]](https://pubmed.ncbi.nlm.nih.gov/18156336/)、[[12813072]](https://pubmed.ncbi.nlm.nih.gov/12813072/)、[[27059479]](https://pubmed.ncbi.nlm.nih.gov/27059479/)：
- **PhaZ1（PhaZa1）**：419 aa / 47.3 kDa，无经典脂酶盒；催化三联体 **Cys183-Asp355-His388**（含半胱氨酸，与 PHB 合酶催化残基位置相似）[[11114905]](https://pubmed.ncbi.nlm.nih.gov/11114905/)、[[16233560]](https://pubmed.ncbi.nlm.nih.gov/16233560/)；主要动员酶，受磷酸化调控 [[29678915]](https://pubmed.ncbi.nlm.nih.gov/29678915/)、[[18156336]](https://pubmed.ncbi.nlm.nih.gov/18156336/)。
- **PhaZ2（= 3HB 寡聚体水解酶）**：水解 3HB 寡聚体（endo/exo 混合型）与人工无定形 PHB [[12775684]](https://pubmed.ncbi.nlm.nih.gov/12775684/)。
- **PhaZc**：新型胞内 3HB 寡聚体水解酶，主要在胞质 [[16030206]](https://pubmed.ncbi.nlm.nih.gov/16030206/)。
- **PhaZd**：与 R. pickettii T1 胞外酶催化域同源但**无信号肽/linker/SBD**，比活高于其他胞内酶 [[16199568]](https://pubmed.ncbi.nlm.nih.gov/16199568/)。
- **PhaZd1/PhaZd2（又称 PhaZ6/PhaZ7）**：体外高活性但体内对 PHB 动员无显著作用，生理功能不明 [[24907326]](https://pubmed.ncbi.nlm.nih.gov/24907326/)。
- 除酶外，R. eutropha 动员还经 **crotonyl-CoA 途径（立体选择性生成 (S)-3-羟基丁酰-CoA）** [[23667237]](https://pubmed.ncbi.nlm.nih.gov/23667237/)。

**其他属胞内酶**：
- Rhodospirillum rubrum：**PhaZ1**（周质定位！催化域与胞外酶 type II 同源，Ser42-Asp138-His178）[[15489436]](https://pubmed.ncbi.nlm.nih.gov/15489436/)；**PhaZ3**（PhaDED 预测胞外、实为胞内，新类型，Ca/Mg 抑制，nPHB 特异）[[21274528]](https://pubmed.ncbi.nlm.nih.gov/21274528/)；**ApdA 激活因子**（= phasin，Mms16 同源，17.4 kDa，耐热至 121 ℃ 以上，改造颗粒表面而非蛋白水解）[[14757249]](https://pubmed.ncbi.nlm.nih.gov/14757249/)、[[15060050]](https://pubmed.ncbi.nlm.nih.gov/15060050/)。
- Bacillus megaterium：**PhaZ1**（nPHB 特异、产物为单体 3HB，可水解变性 PHB）[[19561190]](https://pubmed.ncbi.nlm.nih.gov/19561190/)。
- Bacillus thuringiensis：**PhaZ**（原先注释为 3-oxoadipate enol-lactonase PcaD，G-W-S102-M-G 脂酶盒样序列，S102A 失活，无信号肽）[[16936025]](https://pubmed.ncbi.nlm.nih.gov/16936025/)；首个胞内 PhaZ 晶体结构（1.42 Å，α/β-hydrolase 催化域 + 独特 α-螺旋 cap 域，P-1 亚位点只能容纳 1 个 3HB 单体 → 产物以单体为主）[[39592048]](https://pubmed.ncbi.nlm.nih.gov/39592048/)。
- Azotobacter vinelandii：**PhbZ1**（颗粒结合，含 CoA 依赖性硫解式降解，突变株可得均匀高分子量 PHB）[[29435618]](https://pubmed.ncbi.nlm.nih.gov/29435618/)；phasins PhbP2/PhbP3 参与降解（新 phasin 家族）[[41199521]](https://pubmed.ncbi.nlm.nih.gov/41199521/)。
- Paracoccus denitrificans：**i-PhaZ**（靠近 phaC 基因）[[11267773]](https://pubmed.ncbi.nlm.nih.gov/11267773/)；**PhaZc（3HB 寡聚体水解酶，31 kDa）+ Hbd（3HB/3HV 脱氢酶，四聚体）** 参与 PHB/PHV 胞内降解 [[24271169]](https://pubmed.ncbi.nlm.nih.gov/24271169/)、[[11814660]](https://pubmed.ncbi.nlm.nih.gov/11814660/)。
- **mcl-PHA 胞内解聚酶（Pseudomonas 型）**：P. putida KT2442 PhaZ（α/β-hydrolase + **lid 结构**，催化三联体埋于域间，endo/exo）[[17170116]](https://pubmed.ncbi.nlm.nih.gov/17170116/)；lid 为必需（缺失/突变致失活；S184F 变体改变底物特异性、G286R 提高解聚活性）[[41055782]](https://pubmed.ncbi.nlm.nih.gov/41055782/)；P. oleovorans PHA depolymerase 基因在 pha 位点内 [[1989978]](https://pubmed.ncbi.nlm.nih.gov/1989978/)。PHA 降解与合成**紧密偶联**（futile cycle 调控）[[19267463]](https://pubmed.ncbi.nlm.nih.gov/19267463/)、[[23445364]](https://pubmed.ncbi.nlm.nih.gov/23445364/)。

### 2.5 3HB 寡聚体水解酶（oligomer hydrolases）

- **胞外寡聚体水解酶**：A. faecalis（EC 3.1.1.22，68–74 kDa，从羧基端逐个释放单体）[[6626560]](https://pubmed.ncbi.nlm.nih.gov/6626560/)；Pseudomonas sp. A1（72.9 kDa，无 G-X-S-X-G 盒）[[8981982]](https://pubmed.ncbi.nlm.nih.gov/8981982/)。
- **胞内寡聚体水解酶**：Zoogloea ramigera I-16-M（3HB 二聚体水解酶，28–30 kDa，立体特异）[[7285912]](https://pubmed.ncbi.nlm.nih.gov/7285912/)；R. eutropha H16（78.5 kDa，活性位点丝氨酸 S-V-S*-N-G，突变确定）[[16233278]](https://pubmed.ncbi.nlm.nih.gov/16233278/)；Acidovorax sp. SA1（首例克隆的胞内 3HB 寡聚体水解酶基因，292 aa，G-X-S-X-G 盒）[[12070691]](https://pubmed.ncbi.nlm.nih.gov/12070691/)；P. lemoignei（PhaZc 同源物）[[17964488]](https://pubmed.ncbi.nlm.nih.gov/17964488/)。
- 协同模式：**胞外 PHBDP 将 PHB 水解为寡聚体 → 胞外/胞内 3HBOH 将寡聚体水解为单体**（R. pickettii T1 胞外、Acidovorax SA1 胞内、R. eutropha 胞内三种模式）[[15170237]](https://pubmed.ncbi.nlm.nih.gov/15170237/)。

### 2.6 古菌类型

- **唯一已表征的古菌 PHA 解聚酶：PhaZh1（Haloferax mediterranei）**，**patatin 样蛋白**，颗粒结合，水解 nPHB/nPHBV 产 3HB 单体；关键残基 **Gly16、Ser47（经典脂酶盒 G-X-S47-X-G）、Asp195**；与 bdhA 构成基因簇 HFX_6463–6464 [[25710370]](https://pubmed.ncbi.nlm.nih.gov/25710370/)。
- **H. mediterranei PhaJ1**：颗粒结合型 (R)-特异烯脂酰-CoA 水合酶（R-ECH），将 PHA 降解产物 (R)-3-羟基酰基-CoA 脱水为烯脂酰-CoA 接入 β-氧化，介导 PHA 动员 [[27052994]](https://pubmed.ncbi.nlm.nih.gov/27052994/)。
- 古菌 phasin：H. mediterranei PhaP（phaP 位于 phaEC 上游，ΔphaP 颗粒数显著减少）[[22247127]](https://pubmed.ncbi.nlm.nih.gov/22247127/)。
- 注意：**目前 PubMed 尚无其他古菌（Natronomonas、Halobacterium、Haloarcula 等）的 PHA 解聚酶表征论文**；Sulfolobus acidocaldarius 的热稳定脂解酶仅"与 PHA 解聚酶序列相似"（43% 相似 P. oleovorans 假定解聚酶），并非真正的 PHA 解聚酶 [[9785454]](https://pubmed.ncbi.nlm.nih.gov/9785454/)。搜索 comp_c03（1 篇）、comp_t03（0 篇）、comp_s16（0 篇）证实了该结论。

### 2.7 调控与基因簇

- R. pickettii T1：**epdR（TetR 型转录因子）抑制 phaZ 表达**（3HB 诱导解除抑制）[[24146107]](https://pubmed.ncbi.nlm.nih.gov/24146107/)。
- Pseudomonas putida KT2442：**PhaD 调控因子**同时控制 pha 基因表达与 PHA 周转 [[20406286]](https://pubmed.ncbi.nlm.nih.gov/20406286/)。
- 基因簇结构：Rhodobacter capsulatus SB1003 的 pha2-pha1-phaC-orfX 簇（含解聚酶 orfX）[[12194204]](https://pubmed.ncbi.nlm.nih.gov/12194204/)；P. oleovorans 合成/降解基因成簇 [[1989978]](https://pubmed.ncbi.nlm.nih.gov/1989978/)；Sinorhizobium meliloti PHB 降解途径分布在染色体与巨质粒多座位（bhbA-D、bdhA）[[9258668]](https://pubmed.ncbi.nlm.nih.gov/9258668/)。

---

## 3. 已表征酶清单（细菌 + 真菌 + 古菌）

> 说明：仅列出具实验表征的酶；"催化域分类"按原文表述（type I/II、type A/B、type 1/2、PhaDED 家族等）。Accession 仅在摘要明确给出时列出（不编造）。

### 3.1 胞外 scl-PHA（PHB）解聚酶

| 物种 | 酶名/类型 | 催化域分类 | 关键位点 | PMID | Accession |
|---|---|---|---|---|---|
| Pseudomonas lemoignei | PhaZ1（= depolymerase C） | scl 特异 | lipase box G-X-S-X-G；C 端苏氨酸富集（22/36） | [8269961](https://pubmed.ncbi.nlm.nih.gov/8269961/), [7836292](https://pubmed.ncbi.nlm.nih.gov/7836292/) | — |
| P. lemoignei | PhaZ2（= depolymerase B） | scl | — | [7836292](https://pubmed.ncbi.nlm.nih.gov/7836292/), [12099828](https://pubmed.ncbi.nlm.nih.gov/12099828/) | — |
| P. lemoignei | PhaZ3（= depolymerase D） | scl | — | [7836292](https://pubmed.ncbi.nlm.nih.gov/7836292/) | — |
| P. lemoignei | PhaZ4 | scl（原误判为 PHV 酶） | C 端 SBD 关键 | [7836292](https://pubmed.ncbi.nlm.nih.gov/7836292/), [10742216](https://pubmed.ncbi.nlm.nih.gov/10742216/), [8837471](https://pubmed.ncbi.nlm.nih.gov/8837471/) | — |
| P. lemoignei | PhaZ5（= depolymerase A） | scl | Ser138（G-L-S-S-G） | [8764515](https://pubmed.ncbi.nlm.nih.gov/8764515/), [12007811](https://pubmed.ncbi.nlm.nih.gov/12007811/), [19107395](https://pubmed.ncbi.nlm.nih.gov/19107395/) | — |
| P. lemoignei | PhaZ6（真 PHV depolymerase） | scl/PHV 特异 | 43.6 kDa | [10742216](https://pubmed.ncbi.nlm.nih.gov/10742216/) | — |
| P. lemoignei | PhaZ7（PDI/nPHB 特异） | 新亚群 EC 3.1.1.75 | His47-Ser136-Asp242 + His306；Y105/Y176/Y189/Y190 结合位点 | [12099828](https://pubmed.ncbi.nlm.nih.gov/12099828/), [12855176](https://pubmed.ncbi.nlm.nih.gov/12855176/), [18706425](https://pubmed.ncbi.nlm.nih.gov/18706425/), [20516591](https://pubmed.ncbi.nlm.nih.gov/20516591/), [24007310](https://pubmed.ncbi.nlm.nih.gov/24007310/), [28370478](https://pubmed.ncbi.nlm.nih.gov/28370478/) | — |
| Alcaligenes faecalis T1 | PhaZ | type A（lipase box 位置） | Ser139-Asp214-His273；三域 C/F/S | [2644188](https://pubmed.ncbi.nlm.nih.gov/2644188/), [9371441](https://pubmed.ncbi.nlm.nih.gov/9371441/), [8764515](https://pubmed.ncbi.nlm.nih.gov/8764515/), [3942778](https://pubmed.ncbi.nlm.nih.gov/3942778/), [10408639](https://pubmed.ncbi.nlm.nih.gov/10408639/) | — |
| A. faecalis AE122（海洋） | PhaZ | scl | — | [9177489](https://pubmed.ncbi.nlm.nih.gov/9177489/), [7646009](https://pubmed.ncbi.nlm.nih.gov/7646009/) | — |
| Pseudomonas pickettii / Ralstonia pickettii T1 | PhaZ（type I） | type I（lipase box 位置） | SBD 残基 Leu441/Tyr443/Ser445 等 | [8373740](https://pubmed.ncbi.nlm.nih.gov/8373740/), [16963553](https://pubmed.ncbi.nlm.nih.gov/16963553/), [20058938](https://pubmed.ncbi.nlm.nih.gov/20058938/), [18340545](https://pubmed.ncbi.nlm.nih.gov/18340545/), [24146107](https://pubmed.ncbi.nlm.nih.gov/24146107/) | — |
| Comamonas sp. | PhaZCsp | scl；Fn3 linker | 催化三联体 Ser/His/Asp | [7606660](https://pubmed.ncbi.nlm.nih.gov/7606660/) | — |
| Comamonas acidovorans YM1609 | PhaZ Cac | type B | 45 kDa；SBD 结合 | [9406404](https://pubmed.ncbi.nlm.nih.gov/9406404/), [10408639](https://pubmed.ncbi.nlm.nih.gov/10408639/) | — |
| Comamonas testosteroni YM1004 | PhaZ Cte | scl；Fn3 | SBD 吸附 | [9297825](https://pubmed.ncbi.nlm.nih.gov/9297825/) | — |
| Pseudomonas stutzeri | PhaZPst | type A | 催化域识别 ≥2 个单体；SBDI/II + cadherin-like linker | [9872779](https://pubmed.ncbi.nlm.nih.gov/9872779/), [11710119](https://pubmed.ncbi.nlm.nih.gov/11710119/) | — |
| Marinobacter sp. NK-1 | PhaZ | scl | 578 aa；CD+LD+SBD×2 | [14607367](https://pubmed.ncbi.nlm.nih.gov/14607367/) | — |
| Acidovorax sp. TP4 | PhaZ | type II（lipase box 位置） | — | [12099829](https://pubmed.ncbi.nlm.nih.gov/12099829/), [12217032](https://pubmed.ncbi.nlm.nih.gov/12217032/) | — |
| Acidovorax sp. DP5 | 胞外 PHA depolymerase | scl | 碱性 pH 9 / 40 ℃ | [26664741](https://pubmed.ncbi.nlm.nih.gov/26664741/) | — |
| Schlegelella sp. KB1a | PhaZ | CD type 2 – Fn3 – SBD type 1 | Ser20-Asp104-His138（G-L-S-A-G） | [15340791](https://pubmed.ncbi.nlm.nih.gov/15340791/) | — |
| Thermus thermophilus HB8 | TTHA0199 | scl，exo | Ser183-Glu310-His405 | [19214501](https://pubmed.ncbi.nlm.nih.gov/19214501/) | — |
| Lihuaxuella thermophila | LtPHBase | 广谱（PHB/PLA/PCL） | Ser-His-Asp 三联体；70 ℃ 最适；1.2 Å 结构 | [36222314](https://pubmed.ncbi.nlm.nih.gov/36222314/) | — |
| Streptomyces sp. MG（嗜热） | PhaZ | scl，广谱 | 43 kDa；60 ℃/pH 8.5 | [16614903](https://pubmed.ncbi.nlm.nih.gov/16614903/), [15289671](https://pubmed.ncbi.nlm.nih.gov/15289671/) | — |
| Streptomyces exfoliatus K10 | PhaZ（scl） | scl | — | [8810505](https://pubmed.ncbi.nlm.nih.gov/8810505/) | — |
| Streptomyces ascomycinicus | PhaZS（fkbU） | scl，酸性最适 | Ser131-Asp209-His269；48.4 kDa | [23951224](https://pubmed.ncbi.nlm.nih.gov/23951224/), [21845385](https://pubmed.ncbi.nlm.nih.gov/21845385/) | — |
| Bacillus sp. NRRL B-14911 | PhaZ | CD type 1 + LD1/LD2 + SBD1/SBD2 新型 | 新类别代表 | [21948827](https://pubmed.ncbi.nlm.nih.gov/21948827/) | — |
| Bacillus megaterium N-18-25-9 | PhaZ（胞外） | scl | — | [17064368](https://pubmed.ncbi.nlm.nih.gov/17064368/) | — |
| Cupriavidus malaysiensis | CmaPHBd | 胞外 type I | CD+SBD；酸性 pH 5–6；Tm 49.6 ℃ | [41173112](https://pubmed.ncbi.nlm.nih.gov/41173112/) | — |
| Pseudomonas guguanensis | PguPHBd | 胞外 type I | CD+SBD；碱性 pH 9.5 | [41173112](https://pubmed.ncbi.nlm.nih.gov/41173112/) | — |
| Nocardiopsis dassonvillei NCIM 5124 | PHBD | type I（G-L-S-A-G；N 端氧阴离子 His） | 42.46 kDa；pH 7.5/30 ℃ | [39310033](https://pubmed.ncbi.nlm.nih.gov/39310033/), [41151231](https://pubmed.ncbi.nlm.nih.gov/41151231/) | MCK9871921.1 |
| Alteromonas sp. D210916BOD_24 | PhaZ | 海洋型（lipase box 居中 + 双区 SBD + Fn3 linker） | — | [40500476](https://pubmed.ncbi.nlm.nih.gov/40500476/) | — |
| Streptomyces sp. SNG9（海洋） | PhaZ | scl | — | [11770850](https://pubmed.ncbi.nlm.nih.gov/11770850/) | — |
| Nocardiopsis aegyptia sp. nov.（海洋） | PhaZ | scl | — | [16107752](https://pubmed.ncbi.nlm.nih.gov/16107752/) | — |
| Penicillium funiculosum | PHB depolymerase | 环状排列 α/β-hydrolase，无 SBD | Ser39-Asp121-His155；Trp307 识别 | [16405909](https://pubmed.ncbi.nlm.nih.gov/16405909/), [17547455](https://pubmed.ncbi.nlm.nih.gov/17547455/) | — |
| Penicillium citrinum S2 | PhaZPen | scl，真菌 | 20 kDa 糖蛋白 | [21369777](https://pubmed.ncbi.nlm.nih.gov/21369777/) | — |
| Penicillium expansum | PHAZ Pen | scl，exo | 20 kDa | [28324398](https://pubmed.ncbi.nlm.nih.gov/28324398/) | — |
| Penicillium pinophilum | e-PHB depolymerase | scl，真菌 | — | [25328684](https://pubmed.ncbi.nlm.nih.gov/25328684/) | — |
| Microbacterium paraoxydans RZS6 | PHB depolymerase | scl | — | [31211775](https://pubmed.ncbi.nlm.nih.gov/31211775/) | — |
| Stenotrophomonas sp. RZS7 | PHB depolymerase | scl | — | [31910206](https://pubmed.ncbi.nlm.nih.gov/31910206/), [28330251](https://pubmed.ncbi.nlm.nih.gov/28330251/) | — |
| Aeromonas caviae Kuk1-(34) | e-PHB depolymerase | scl | — | [35421107](https://pubmed.ncbi.nlm.nih.gov/35421107/), [37639167](https://pubmed.ncbi.nlm.nih.gov/37639167/) | — |

### 3.2 胞外 mcl-PHA 解聚酶

| 物种 | 酶名 | 特征 | PMID |
|---|---|---|---|
| Pseudomonas alcaligenes LB19 | mcl-PHA depolymerase | 27.6 kDa，pH 9/45 ℃，产单体 | [11888314](https://pubmed.ncbi.nlm.nih.gov/11888314/) |
| P. alcaligenes M4-7 | PhaZPalM4-7 | 28 kDa，pH 9/35 ℃ | [15995648](https://pubmed.ncbi.nlm.nih.gov/15995648/) |
| Pseudomonas sp. RY-1 | mcl-PHA depolymerase | 四聚体 115 kDa，pH 8.5/35 ℃ | [16232726](https://pubmed.ncbi.nlm.nih.gov/16232726/) |
| Xanthomonas sp. JS02 | PHPV depolymerase | 41.7 kDa，芳香族 mcl | [10772473](https://pubmed.ncbi.nlm.nih.gov/10772473/) |
| Streptomyces sp. KJ-72 | mcl-PHA depolymerase | 27.1 kDa，pH 8.7/50 ℃，产二聚体 | [12785312](https://pubmed.ncbi.nlm.nih.gov/12785312/) |
| Streptomyces exfoliatus K10 DSMZ 41693 | PhaZSex2 | 27.6 kDa，endo-exo，产 (R)-3-羟基辛酸单体；可降解功能化 PHACOS | [26156240](https://pubmed.ncbi.nlm.nih.gov/26156240/) |
| Streptomyces venezuelae SO1 | mcl-PHA depolymerase | 27 kDa，pH 碱性/50 ℃，可逆酯合成活性 | [22695803](https://pubmed.ncbi.nlm.nih.gov/22695803/) |
| Streptomyces roseolus SL3 等放线菌 | mcl-PHA depolymerase 新亚群 | 28 kDa，pH 9.5，水解 PCL 与 pNP 酯 | [22865072](https://pubmed.ncbi.nlm.nih.gov/22865072/) |
| Bdellovibrio bacteriovorus HD100 | PhaZBd（Bd3709） | α/β-hydrolase，endo-exo，mcl 特异 | [22706067](https://pubmed.ncbi.nlm.nih.gov/22706067/) |
| Pseudomonas putida KT2442/2440 | PhaZ（胞内 mcl） | lid 结构（见 2.4） | [17170116](https://pubmed.ncbi.nlm.nih.gov/17170116/), [41055782](https://pubmed.ncbi.nlm.nih.gov/41055782/) |
| Pseudomonas resinovorans | i-PhaZ | mcl 解聚酶基因敲除影响 PHA | [12759786](https://pubmed.ncbi.nlm.nih.gov/12759786/) |

### 3.3 胞内酶与寡聚体水解酶

| 物种 | 酶名 | 类型 | 关键位点 | PMID |
|---|---|---|---|---|
| Ralstonia (Cupriavidus) eutropha H16 | PhaZ1/PhaZa1 | 胞内 | Cys183-Asp355-His388 | [11114905](https://pubmed.ncbi.nlm.nih.gov/11114905/), [16233560](https://pubmed.ncbi.nlm.nih.gov/16233560/), [18156336](https://pubmed.ncbi.nlm.nih.gov/18156336/), [29678915](https://pubmed.ncbi.nlm.nih.gov/29678915/) |
| R. eutropha H16 | PhaZ2 | 胞内 3HB 寡聚体水解酶 | — | [12775684](https://pubmed.ncbi.nlm.nih.gov/12775684/), [16233278](https://pubmed.ncbi.nlm.nih.gov/16233278/) |
| R. eutropha H16 | PhaZc | 胞内 3HB 寡聚体水解酶（新） | — | [16030206](https://pubmed.ncbi.nlm.nih.gov/16030206/) |
| R. eutropha H16 | PhaZd | 胞内，与胞外酶 CD 同源 | Ser190-Asp266-His330 | [16199568](https://pubmed.ncbi.nlm.nih.gov/16199568/) |
| R. eutropha H16 | PhaZd1/PhaZd2 | 胞内（体内功能不明） | S190/S193 | [24907326](https://pubmed.ncbi.nlm.nih.gov/24907326/) |
| Rhodospirillum rubrum | PhaZ1（周质） | 胞外酶 type II CD 同源 | Ser42-Asp138-His178 | [15489436](https://pubmed.ncbi.nlm.nih.gov/15489436/) |
| R. rubrum | PhaZ3Rru | 胞内新类型 | — | [21274528](https://pubmed.ncbi.nlm.nih.gov/21274528/) |
| R. rubrum | ApdA（激活因子/phasin） | 颗粒表面蛋白 | Mms16 同源，55% 相同 | [14757249](https://pubmed.ncbi.nlm.nih.gov/14757249/), [15060050](https://pubmed.ncbi.nlm.nih.gov/15060050/) |
| Magnetospirillum gryphiswaldense | Mms16/ApdA | phasin（激活剂） | — | [15774885](https://pubmed.ncbi.nlm.nih.gov/15774885/) |
| Bacillus megaterium | PhaZ1（胞内） | nPHB 特异，产单体 | — | [19561190](https://pubmed.ncbi.nlm.nih.gov/19561190/) |
| Bacillus thuringiensis israelensis | PhaZ（胞内） | 新类型；1.42 Å 结构 | G-W-S102-M-G；S102A 失活 | [16936025](https://pubmed.ncbi.nlm.nih.gov/16936025/), [39592048](https://pubmed.ncbi.nlm.nih.gov/39592048/), [25286954](https://pubmed.ncbi.nlm.nih.gov/25286954/) |
| Azotobacter vinelandii | PhbZ1（胞内） | 含硫解式降解 | — | [29435618](https://pubmed.ncbi.nlm.nih.gov/29435618/) |
| A. vinelandii | PhbP2/PhbP3 | phasin（参与降解） | 新 phasin 家族（PhbP3） | [41199521](https://pubmed.ncbi.nlm.nih.gov/41199521/) |
| Paracoccus denitrificans | i-PhaZ + PhaZc + Hbd | 胞内 + 3HB 寡聚体水解酶 + 脱氢酶 | — | [11267773](https://pubmed.ncbi.nlm.nih.gov/11267773/), [24271169](https://pubmed.ncbi.nlm.nih.gov/24271169/), [11814660](https://pubmed.ncbi.nlm.nih.gov/11814660/) |
| Pseudomonas oleovorans | PHA depolymerase | 胞内 mcl，基因在 pha 位点内 | — | [1989978](https://pubmed.ncbi.nlm.nih.gov/1989978/) |
| Pseudomonas putida KT2442 | PhaZ | 胞内 mcl（paradigmatic） | lid 结构 | [17170116](https://pubmed.ncbi.nlm.nih.gov/17170116/), [19788655](https://pubmed.ncbi.nlm.nih.gov/19788655/), [23445364](https://pubmed.ncbi.nlm.nih.gov/23445364/), [19103481](https://pubmed.ncbi.nlm.nih.gov/19103481/), [25563970](https://pubmed.ncbi.nlm.nih.gov/25563970/), [28952537](https://pubmed.ncbi.nlm.nih.gov/28952537/) |
| Zoogloea ramigera I-16-M | 3HB 二聚体水解酶 | 胞内寡聚体水解酶 | 28–30 kDa | [7285912](https://pubmed.ncbi.nlm.nih.gov/7285912/), [1476778](https://pubmed.ncbi.nlm.nih.gov/1476778/) |
| Pseudomonas sp. A1 | 3HB 寡聚体水解酶（胞外） | EC 3.1.1.22 | 72.9 kDa，无 G-X-S-X-G | [8981982](https://pubmed.ncbi.nlm.nih.gov/8981982/) |
| Alcaligenes faecalis | 3HB 寡聚体水解酶（胞外） | EC 3.1.1.22 | 68–74 kDa | [6626560](https://pubmed.ncbi.nlm.nih.gov/6626560/) |
| Acidovorax sp. SA1 | i3HBOH（胞内） | 首例克隆 | 292 aa；G-X-S-X-G 盒 | [12070691](https://pubmed.ncbi.nlm.nih.gov/12070691/) |
| Paucimonas lemoignei | 3HB 寡聚体水解酶（胞内） | PhaZc 同源 | — | [17964488](https://pubmed.ncbi.nlm.nih.gov/17964488/) |
| Wautersia/R. eutropha | 3HB-oligomer hydrolase（胞内） | — | S-V-S*-N-G | [16233278](https://pubmed.ncbi.nlm.nih.gov/16233278/) |

### 3.4 3HB 脱氢酶（bdhA / BDH）

| 物种 | 酶 | 特征 | PMID |
|---|---|---|---|
| Rhizobium (Sinorhizobium) meliloti | BdhA | 258 aa，短链醇脱氢酶超家族；bhbA-D 多位点 | [9922248](https://pubmed.ncbi.nlm.nih.gov/9922248/), [9258668](https://pubmed.ncbi.nlm.nih.gov/9258668/), [12632261](https://pubmed.ncbi.nlm.nih.gov/12632261/), [16175209](https://pubmed.ncbi.nlm.nih.gov/16175209/) |
| Sinorhizobium sp. NGR234 | BdhA | 91% 相同于 S. meliloti | [15621424](https://pubmed.ncbi.nlm.nih.gov/15621424/) |
| Ralstonia pickettii T1 | BDH1/BDH2/BDH3 | 三种 BDH 生理角色不同 | [16935252](https://pubmed.ncbi.nlm.nih.gov/16935252/), [19219638](https://pubmed.ncbi.nlm.nih.gov/19219638/) |
| Paracoccus denitrificans | Hbd | 四聚体 29 kDa×4；3HB/3HV 均可 | [24271169](https://pubmed.ncbi.nlm.nih.gov/24271169/) |
| Legionella pneumophila | BdhA | 与假定脂酶基因相邻 | [25556866](https://pubmed.ncbi.nlm.nih.gov/25556866/) |
| Haloferax mediterranei | BdhA（推定的） | 与 phaZh1 成簇 HFX_6463–6464 | [25710370](https://pubmed.ncbi.nlm.nih.gov/25710370/) |
| Pseudomonas aeruginosa PAO1 | (R)-3HB 代谢 | PA2005/RpoN 调控 | [26311173](https://pubmed.ncbi.nlm.nih.gov/26311173/) |

### 3.5 古菌酶（含推定）

| 物种 | 酶 | 类型 | 关键位点 | PMID |
|---|---|---|---|---|
| Haloferax mediterranei | **PhaZh1**（patatin-like） | **唯一已表征古菌 PHA 解聚酶**，nPHB/nPHBV → 3HB | Gly16、Ser47（G-X-S47-X-G）、Asp195 | [25710370](https://pubmed.ncbi.nlm.nih.gov/25710370/) |
| H. mediterranei | PhaJ1（R-ECH） | PHA 动员：脱水 (R)-3-羟基酰基-CoA 接入 β-氧化 | PhaJ1–PhaJ5 中仅 PhaJ1 颗粒结合 | [27052994](https://pubmed.ncbi.nlm.nih.gov/27052994/) |
| H. mediterranei | PhaP | phasin（颗粒结构） | phaP 位于 phaEC 上游 | [22247127](https://pubmed.ncbi.nlm.nih.gov/22247127/) |
| Sulfolobus acidocaldarius DSM 639 | 热稳定脂解酶 | 类 PHA depolymerase（43% 相似 P. oleovorans 假定酶） | G-X-S-X-G；314 aa | [9785454](https://pubmed.ncbi.nlm.nih.gov/9785454/) |
| Haloarcula marismortui | PhaE/PhaC（合成酶，非降解） | class III PHA 合酶 | — | [17675423](https://pubmed.ncbi.nlm.nih.gov/17675423/) |
| Haloarcula hispanica | FabG1/PhaB（合成） | acetoacetyl-CoA 还原酶 | — | [19648370](https://pubmed.ncbi.nlm.nih.gov/19648370/) |
| 极端嗜盐古菌 strain 56 | PHB 合酶（合成） | 80 kDa 颗粒结合 | — | [12139978](https://pubmed.ncbi.nlm.nih.gov/12139978/) |

---

## 4. 古菌 PHA 降解酶文献全集

经过 15+ 个古菌定向查询（comp_c01–c07、comp_s19/s20/s28、comp_t02/t03/t07/t09、comp_u01/u03），PubMed 中与**古菌 PHA 降解（而非合成）**直接相关的文献全集如下：

**A. 已表征降解酶（仅 2 个酶，均来自 Haloferax mediterranei）**
1. PhaZh1 — patatin-like PHA 解聚酶，唯一已表征的古菌 PHA 解聚酶 [[25710370]](https://pubmed.ncbi.nlm.nih.gov/25710370/)（Appl Environ Microbiol, 2015）
2. PhaJ1 — (R)-特异烯脂酰-CoA 水合酶，PHA 动员 [[27052994]](https://pubmed.ncbi.nlm.nih.gov/27052994/)（Sci Rep, 2016）

**B. 降解相关蛋白/颗粒生物学**
3. PhaP phasin（H. mediterranei）[[22247127]](https://pubmed.ncbi.nlm.nih.gov/22247127/)
4. H. mediterranei PHA 颗粒染色分析（方法学）[[34069083]](https://pubmed.ncbi.nlm.nih.gov/34069083/)
5. 古菌能量储备代谢通路生物信息学分析 [[30705313]](https://pubmed.ncbi.nlm.nih.gov/30705313/)

**C. 环境降解（haloarchaea 产 PHA 的生物降解）**
6. haloarchaea 产 PHBV 在活性污泥中的环境生物降解 [[27098259]](https://pubmed.ncbi.nlm.nih.gov/27098259/)
7. haloarchaea 产 PHBV 共聚物降解与生物相容性 [[28618347]](https://pubmed.ncbi.nlm.nih.gov/28618347/)

**D. 类 PHA 解聚酶序列（古菌）**
8. Sulfolobus acidocaldarius 热稳定脂解酶（序列相似）[[9785454]](https://pubmed.ncbi.nlm.nih.gov/9785454/)

**E. 古菌 PHA 合成（对照，非降解）**
9. Haloarcula marismortui phaECHm [[17675423]](https://pubmed.ncbi.nlm.nih.gov/17675423/)；H. hispanica FabG1/PhaB [[19648370]](https://pubmed.ncbi.nlm.nih.gov/19648370/)；strain 56 PHB 合酶 [[12139978]](https://pubmed.ncbi.nlm.nih.gov/12139978/)；Halococcus 等 PHA 检测 [[20437233]](https://pubmed.ncbi.nlm.nih.gov/20437233/)；PHA 作为生态货币综述 [[42203387]](https://pubmed.ncbi.nlm.nih.gov/42203387/)

**结论**：古菌 PHA 降解酶的表征研究处于极早期——**PubMed 中仅 2 个 H. mediterranei 酶被实验表征（PhaZh1、PhaJ1）**，而 **Natronomonas、Halobacterium、Haloarcula（除合成外）、Haloquadratum、Natrialba 等属尚无任何 PHA 解聚酶表征报道**（comp_t03 与 comp_s16 查询 0 命中）。这为 GTDB 全库筛查"古菌 PHB 降解基因"提供了明确的知识空白定位：GTDB 古菌门中 PHB 降解基因（PhaZh1 同源、patatin 域、PhaJ 类 R-ECH、bdhA 同源）的分布是一个未被文献覆盖的空白。

---

## 5. 关键论文 URL 列表

全部 748 篇 PMID 的完整 URL 清单见同目录 `comp_pubmed_urls.txt`（每行一个 `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`）。核心分类/综述/方法论文：

- Jendrossek 综述 2002（Annu Rev Microbiol）: https://pubmed.ncbi.nlm.nih.gov/12213937/
- Jendrossek 综述 1996: https://pubmed.ncbi.nlm.nih.gov/9008883/
- 综述 2004（Degradation of microbial polyesters）: https://pubmed.ncbi.nlm.nih.gov/15289671/
- PhaDED 数据库 2009（BMC Bioinformatics）: https://pubmed.ncbi.nlm.nih.gov/19296857/
- P. lemoignei 系统表征: https://pubmed.ncbi.nlm.nih.gov/7836292/
- PhaZ7 结构（JMB）: https://pubmed.ncbi.nlm.nih.gov/18706425/
- PhaZ7 1.2 Å 结构: https://pubmed.ncbi.nlm.nih.gov/20516591/
- P. funiculosum 结构: https://pubmed.ncbi.nlm.nih.gov/16405909/
- B. thuringiensis 胞内 PhaZ 结构（2025）: https://pubmed.ncbi.nlm.nih.gov/39592048/
- R. eutropha 9 基因分类: https://pubmed.ncbi.nlm.nih.gov/18156336/
- PhaZh1（古菌解聚酶）: https://pubmed.ncbi.nlm.nih.gov/25710370/
- PhaJ1（古菌动员）: https://pubmed.ncbi.nlm.nih.gov/27052994/
- 海洋宏基因组 PHB 降解（2025）: https://pubmed.ncbi.nlm.nih.gov/39827799/
- 苏打湖宏基因组 PHB 解聚酶（2026）: https://pubmed.ncbi.nlm.nih.gov/41702408/
- CmaPHBd/PguPHBd 结构域改组（2026）: https://pubmed.ncbi.nlm.nih.gov/41173112/

---

*报告由 pubmed-database skill 生成；所有 PMID/结论均来自检索结果，未编造 accession（仅 MCK9871921.1 为摘要明确给出）。*
