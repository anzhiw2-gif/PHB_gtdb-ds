# PHB 降解基因家族分类方案（v0.3 — 最终版）

> 状态：v0.3 定稿。基于三路穷尽式调研：
> ① Europe PMC 全文（comp_classification.md）
> ② PubMed 728 篇（comprehensive_summary.md）
> ③ OpenAlex/Web 数据库体系（comp_databases_report.md）

## 1. 分类逻辑（两维）

PHA 解聚酶（EC 3.1.1.75 SCL / EC 3.1.1.76 MCL）均属 α/β-水解酶折叠，
催化 Ser-His-Asp/Glu 三联体（个别 Cys-His-Asp/Glu；酸性残基 Asp 与 Glu 可互换，
如 *T. thermophilus* HB8 的 Ser183-Glu310-His405），催化 Ser 在 Gx₁Sx₂G lipase box：

1. **定位×底物状态**：胞内 i-天然颗粒 nPHA / 胞外 e-变性颗粒 dPHA / 周质 / 胞外天然 nPHA
2. **链长**：SCL（3-5 C，PHB/PHV）vs MCL（6-15 C，P(3HO)）

## 2. PhaDED 八大超家族（Knoll 2009；含 38 同源家族、587 蛋白）

| # | 超家族 | 定位/底物 | 催化 | lipase box | 代表（gi/UniProt） |
|---|--------|----------|------|-----------|-------------------|
| 1 | **i-nPHASCL（无 lipase box）** | 胞内/天然 SCL | **Cys-His-Asp** | 无 GxSxG；Cys-1 位疏水 Val | *R. eutropha* PhaZ1/2/3/5、PhaZd（PhaDED 最大族 224 条） |
| 2 | **i-nPHASCL（有 lipase box）** | 胞内/天然 SCL | Ser-His-Asp | x₁=Trp | *B. thuringiensis* |
| 3 | **周质 PHA 解聚酶** | 周质/天然 PHB | Ser-His-Asp | type 2 催化域；x₁=Ile | *R. rubrum* PhaZ1 |
| 4 | **i-nPHAMCL** | 胞内/天然 MCL | Ser-His-Asp | x₁=Val；**lid 结构域**（类脂酶） | *P. putida* KT2440 PhaZKT（S102-D221-H248） |
| 5 | **e-dPHASCL（催化域 type 1）** | 胞外/变性 SCL | Ser-His-Asp | x₁=Leu/Ile（81% 疏水）；oxyanion hole 在 box N 端 | *A. faecalis*、*R. pickettii* T1、*B. megaterium*、*P. lemoignei* 多酶、*P. stutzeri*（PhaDED 最大族 234 条） |
| 6 | **e-dPHASCL（催化域 type 2）** | 胞外/变性 SCL | Ser-His-Asp | x₁ 全部疏水；oxyanion hole 在三元组 C 端 | *Acidovorax* TP4、*Caldimonas*、*Comamonas*、*Delftia*、*Schlegelella*、*S. exfoliatus*、*P. funiculosum*（PDB 2D80） |
| 7 | **e-nPHASCL（PhaZ7 型）** | 胞外/仅天然 SCL | Ser-His-Asp | **AHSMG**（非 GxSxG）；14aa 移动 lid | *P. lemoignei* PhaZ7（PDB 2VTV） |
| 8 | **e-dPHAMCL** | 胞外/变性 MCL | Ser-His-Asp | x₁=Ile；无 SBD/linker | *P. fluorescens* GK13、*P. alcaligenes*、*T. thermophilus* HB8 |

结构域架构（胞外经典）：信号肽 → 催化域 → 连接域（Fn3/Thr-rich/Cad）→ SBD（SBD1/SBD2）
- dPHAMCL 无 SBD/linker（N 端结合位点）
- Pfam 关键标记（**2026-09-20 更正后**）：**PF10503 / IPR010126**（Esterase_PHB，**胞外 type-1/type-2 催化域**；GO=extracellular；头号引用 = A. faecalis 胞外解聚酶 PMID 2644188）；**PF06850 / IPR009656**（PHB_depo_C，**胞内 Cys 型酶的 C 端**，描述为"degrades PHB granules"）；胞外 SBD/连接配件域 = FN3(PF00041)/Ig-like(PF16403,PF17957)/TSP3(PF02412)/CHB(PF13290)/Thr-rich/Cad。

> **更正痕迹（2026-09-20）：** 原句「+ **PF06850**（PHB 解聚酶 C 端 SBD）」把 PF06850 列为**胞外经典架构的 SBD**——**归属有误**。PF06850 是**胞内（Cys 型）** PHB 解聚酶的 C 端标记，不是胞外 type-1 的 SBD；type-1 的催化域是 PF10503。依据：InterPro `IPR009656`（"degrades PHB granules"）与 `IPR010126`（GO extracellular；头号引用 PMID 2644188 = A. faecalis 胞外解聚酶）；诊断 run `20260920_phaded_three_gaps_01` 任务③（verdict=coverage_artifact）。此更正亦使 type-1 高可信度候选「缺 PF06850」不再被读作「缺 SBD」。

## 3. 胞内酶细目（*R. eutropha* 全基因组清单 [B12]）

| 基因 | 功能 | PhaDED 归属 | 蛋白 ID | 位点 |
|------|------|------------|---------|------|
| phaZ1 | 胞内 PHB 解聚酶 | 家族1 | CAJ92291.1 | H16_A1150 |
| phaZ2 | 胞内（效率高） | 家族1 | CAJ93939.1 | H16_A2862 |
| phaZ3 | 推定 | 家族1 | CAJ95139.1 | H16_B0339 |
| phaZ4 | 推定 | 家族1 | AAP85930.1 | PHG178 |
| phaZ5 | 胞内 | 家族1 | CAJ95805.1 | H16_B1014 |
| phaZ6 | PHB 解聚酶 | 家族5（胞外型，胞内基因！） | CAJ96855.1 | H16_B2073 |
| phaZ7 | PHB 解聚酶 | 家族5（胞外型，胞内基因！） | CAJ97183.1 | H16_B2401 |
| phaY1/phaY2 | **寡聚体水解酶**（EC 3.1.1.22） | 不在 PhaDED | CAJ93348.1/CAJ92475.1 | H16_A2251/A1335 |

> ⚠️ 注释陷阱：PhaZ6/Z7 为胞外型序列但属胞内基因；GenBank 注释不可靠，
> **必须以序列聚类为准**。

## 4. 古菌 PHA 降解（独立家族）

| 家族 | 特征 | 代表 |
|------|------|------|
| **ArchPhaZ_patatin**（PhaZh1 型）**⚠️ 仅作 patatin 折叠召回探针** | patatin/磷脂酶 A₂ 折叠；**有自己的 GxSxG 序列 motif，但结构上不是 α/β-水解酶折叠的 lipase box**（见 §4.1 更正痕迹）；颗粒结合 PGAP；体外高效解聚 nPHB(V)，体内角色有限（敲除不影响动员→替代通路） | *Hfx. mediterranei* PhaZh1（I3RBH0/HFX_6463）、*Natronomonas* M1XPT2 |
| **ArchPhaJ**（MaoC R-ECH） | MaoC 域（**PF01575**）；Asp-His 二元组；R 立体专一；(R)-3HA-CoA→烯酰-CoA→β-氧化；**古菌主要动员途径**（96% 含 phaJ 卤古菌具完整 β-氧化） | *Hfx. mediterranei* PhaJ1（HFX_5217，219aa） |
| ArchPhaZ_hydrolase | 经典 PHB 解聚酶家族酯酶（12 条 Halobacteria；**无实验锚点，见 §4.2**） | Haladaptatus 等 |

古菌 PHA 合成：Class III（PhaC+PhaE），仅 scl-PHA。

**非嗜盐古菌旁证（v0.3）**：*Sulfolobus acidocaldarius*（泉古菌）热稳定
脂解酶与 PHA 解聚酶序列显著相似（Arpigny & Jendrossek 1998, FEMS 167:69）
——唯一非嗜盐古菌的序列同源证据（PHA 水解活性未验证）。
Methanosarcina/Thermococcus/Pyrococcus 无实验性 PHA 降解报道（OpenAlex
确认）；宏基因组研究聚焦 phaC 合成基因，未见古菌 contig phaZ 专门报道。

### 4.1 更正痕迹（2026-09-20）—— patatin「非 GxSxG」是事实错误

> 本节由文档审查（交接文档「接下来该做什么」第 4 项）新增。**原表述保留在下方逐字引用中，本表上方已改为更正后措辞。**

- **被更正的原文措辞（逐字引用）：** 「| **ArchPhaZ_patatin**（PhaZh1 型） | patatin/磷脂酶 A₂ 折叠（**非 GxSxG**）；颗粒结合 PGAP；…」
- **为何错误（两条独立证据，均指向具体产物）：**
  1. **一手文献**：[Liu 2015, PMID 25710370](https://pubmed.ncbi.nlm.nih.gov/25710370/)（*Haloferax mediterranei* PhaZh1）原文明写其催化残基为「**Ser47 (in a classical lipase box, G-X-S47-X-G)**」——即 PhaZh1 **有** GxSxG 序列 motif。
  2. **本工作区实测**：`runs/20260906_overall_pha_phb_archaea_audit_04/results/archaeal_phazh1_full_audit.tsv` 的 `nterm_gxsxg_motif` 列，2,307 行中 **`present` 2,118 / `absent` 189**（本次独立逐行重算确认；该 run 的 `results/analysis_manifest.json` 记 `counts.archaeal_records = 2307`）。
- **准确表述（更正后）：** patatin 折叠**有自己的 GxSxG 序列 motif**，但它与 α/β-水解酶折叠的 lipase box **结构语境不同**（patatin 是独立折叠，催化 Ser 位于 patatin 折叠自身的 GxSxG 中）。因此：
  - **不能**用「无 GxSxG」区分 patatin 型与 α/β-水解酶型；
  - **不能**把 §6.2 的 x₁ 疏水性 / lipase box 判据套用到 patatin 候选上打分。

### 4.2 正式降级声明（2026-09-20）—— `ArchPhaZ_patatin` = patatin 折叠**召回探针**

`ArchPhaZ_patatin` 与 `ArchPhaZ_hydrolase` 两个古菌 HMM **不得**被用于声称「古菌 PHB 解聚酶」。降级依据（全部为实测，各指向具体产物）：

1. **校准无法建立零假阳性阈值**：`runs/20260904_archaea_phazh1_candidate_01/results/calibration/calibration.tsv` —— 支持面板 8/8 命中（E = 7.4e-145 … 1.5e-93），但**广谱/冲突对照（`negative_control`）4/4 全部命中**（9.2e-51、6.1e-28、3.7e-23、**5.4e-06**，四条均 ≤ 1e-5）。故该 HMM **独立特异性不足**；未注册进 `pipeline/config/formal_scan_models.tsv`，未启动 GTDB 正式扫描（见同 run `docs/T141_20260904_archaea_phazh1_candidate_01_status.md`）。
2. **实验锚点只有 1 个**：古菌 PHA 解聚酶中唯一被实验鉴定的只有 PhaZh1（Liu 2015，**体外**天然 PHB/PHBV 颗粒水解），且**体内角色有限**——动员主路为 PhaJ β-氧化（Liu 2016, [PMID 27052994](https://pubmed.ncbi.nlm.nih.gov/27052994/)）。patatin 折叠本身是**广谱磷脂酶 A₂/酯酶折叠**，绝大多数成员不是 PHB 解聚酶。
3. **`ArchPhaZ_hydrolase` 同样无实验锚点**：其 HMM 仅 **12 条** Halobacteria 种子（`data/hmms/v2/ArchPhaZ_hydrolase.hmm` 头部 `NSEQ  12`），只有序列相似性支持；唯一非嗜盐古菌旁证是 [Arpigny & Jendrossek 1998, PMID 9785454](https://pubmed.ncbi.nlm.nih.gov/9785454/) 的 *Sulfolobus acidocaldarius* 热稳定脂解酶（与 PHA 解聚酶 43% 相似、含 GxSxG），但**PHA 水解活性从未验证**。
4. **召回结果绝大多数是 patatin 背景**：`runs/20260906_overall_pha_phb_archaea_audit_04/results/analysis_manifest.json` 的 `counts.archaeal_audit` = `review_required` 204（PhaZh1_like_high 197 + review 7）/ `exploratory` 2,027 / `excluded` 76，合计 2,307。

**允许的表述**：命名与角色一律写作「patatin 折叠**召回探针**（recall-only，只召回不筛选）」。
**禁止的表述**：「古菌 PHB 解聚酶」「patatin 型解聚酶」等结论性措辞。

## 5. 寡聚体水解酶（PhaY）

- EC 3.1.1.22；不在 PhaDED 八大家族
- *R. eutropha* PhaY1/PhaY2；产物从 DP1 到 DP2/混合不等
- 下游：(R)-3HB → BdhA → 乙酰乙酸 → AACoA 合成酶 → 乙酰乙酰-CoA → 2×乙酰-CoA

## 6. HMM 建模与验证规则（文献证据）

### 6.0 现成 HMM 资源（可直接下载，v0.3 新增）
| 来源 | 家族 HMM | 说明 |
|------|---------|------|
| **TIGRFAM** | **TIGR01840**（esterase_phb）、**TIGR01849**（PHB_depoly_PhaZ） | 最直接可用的现成 HMM；jhttp://tigrfams.jcvi.org/ |
| **Pfam** | **PF10503**（Esterase PHB depolymerase 催化域，CDD=Esterase_PHB，**胞外 type-1/type-2**）；**PF06850**（PHB_depo_C，**胞内 Cys 型 C 端**，非胞外 SBD——见第 32 行更正痕迹） | 已从服务器 Pfam-A.hmm 提取 |
| **Pfam** | PF09361（phasin_2）、PF01575（MaoC/PhaJ） | 辅助家族 |
| ESTHER | Esterase_phb_PHAZ、PHAZ7_phb_depolymerase 家族页 | 家族成员与 HMM |
| PhaDED | 8 超家族/38 家族各带 profile HMM | 权威分类参照 |
| SCOP | PHB depolymerase-like family（b.69） | 折叠层验证 |

### 6.1 验证四基序（LtPHBase / PDB 8DAJ [T22]，e-PHASCL 适用）
| 基序 | 模式 | LtPHBase 位置 |
|------|------|--------------|
| Ser | `IDXXXXYVXGLSXGG` | ~aa109-124 |
| Asp | `GXXDYTV` | ~aa194-200 |
| His | `GMXHXXPXXG` | ~aa267-274 |
| oxyanion hole | `HGCXQ` | ~aa38-42 |

催化三元组：Ser121-His270-Asp197（LtPHBase 编号）。

### 6.2 lipase box 疏水 x₁ 判别（[K09]）
- PHA 解聚酶 x₁ 几乎全为疏水（Leu/Ile/Trp/Val）；脂酶/酯酶 x₁ 多为极性
- 例外：PhaZ7 型 AHSMG

### 6.3 建模策略（最终，v0.3）
1. **按三层划分建 HMM**：折叠层（α/β-水解酶 vs patatin）→ 超家族层
   （PhaDED 8 超家族扩展为 8 类 + 古菌 3 类）→ 家族层（系统发育二次切分）
2. **现成 HMM 可直接用**：Pfam PF10503（Esterase_PHB）、TIGRFAM
   TIGR01840/TIGR01849、PhaDED/ESTHER 家族 profile；古菌用
   patatin（Pfam Patatin）+ MaoC（PF01575）独立建模
3. 种子带实验证据标签（DOI + UniProt/PDB）
4. 二级验证：四基序 + x₁ 疏水性 + 信号肽（胞外）+ lid（胞内 mcl）
   + 颗粒结合注释（胞内）
5. 跨类别同源处理：PhaZ6/Z7 型按序列聚类归属（注释不可靠）
6. **坑**：patatin 与 α/β-水解酶 phaZ 折叠不同不可共用 HMM；胞外/胞内
   序列相似度低必须用 HMM 而非 BLAST 单阈值；PhaZ6/7 体外活性≠体内功能

### 6.4 最终筛选家族清单（06_screen.sh 已配置）
| 家族 | HMM（v2） | 分类依据 | 用途 |
|------|-----------|---------|------|
| ePhaZ | ePhaZ.hmm（4,458 种子→3,002） | PhaDED 家族5+6+8 | 胞外解聚酶（scl+mcl） |
| iPhaZ | iPhaZ.hmm（152→112） | PhaDED 家族1+2+3+4 | 胞内解聚酶（含周质） |
| OH | OH.hmm（713→~500） | EC 3.1.1.22 / PhaY | 寡聚体水解酶 |
| BdhA | BdhA.hmm（5,903→~4000） | EC 1.1.1.30 | 3HB 脱氢酶 |
| ArchPhaZ_patatin | ArchPhaZ_patatin.hmm（113→103） | patatin 折叠（**有自己的 GxSxG**，见 §4.1） | **patatin 折叠召回探针（只召回不筛选；见 §4.2）** |
| ArchPhaZ_hydrolase | ArchPhaZ_hydrolase.hmm（12） | 经典家族（12 条 Halobacteria，**无实验锚点**） | **召回探针（见 §4.2）** |
| PhaJ（辅助） | PhaJ.hmm（1,071） | PF01575 MaoC | 动员通路（簇背景） |
| phasin（辅助） | phasin.hmm（Pfam PF09361） | Phasin_2 | 颗粒蛋白 |
| PhaC（辅助） | PhaC.hmm（47） | K03821 | 合酶（簇共定位） |

## 7. 与 v2 种子库的对应

| 分类方案家族 | v2 种子（UniProt） | 需补充 |
|-------------|-------------------|--------|
| 家族5+6 e-dPHASCL | ePhaZ.faa (4,458) | 按 type 1/2 分型（待聚类） |
| 家族8 e-dPHAMCL | ePhaZ.faa 中 MCL 部分 | 需 MCL 参考（GK13 等） |
| 家族1 i-nPHASCL Cys 型 | iPhaZ.faa (152) | 需 Cys-His-Asp 校验 |
| 家族2+4 i-nPHASCL Ser 型 / i-nPHAMCL | iPhaZ.faa | lid 域校验 |
| 家族3 周质 | iPhaZ.faa（R. rubrum） | 信号肽校验 |
| 家族7 PhaZ7 型 | ePhaZ.faa（P. lemoignei PhaZ7） | AHSMG 校验 |
| OH/PhaY | OH.faa (713) | EC 3.1.1.22 |
| ArchPhaZ_patatin | ArchPhaZ_patatin.faa (113) | ✅ |
| ArchPhaZ_hydrolase | ArchPhaZ_hydrolase.faa (12) | ✅ |
| ArchPhaJ (MaoC) | PhaJ.faa (1,071) | PF01575 过滤 |
| BdhA | BdhA.faa (5,903) | ✅ |
| phasin | phasin.faa / Pfam PF09361 | ✅ |
| PhaC（簇背景） | PhaC.faa (47) | Class I-III 区分 |

## 8. 参考文献

- [K09] Knoll 2009 PhaDED, BMC Bioinformatics 10:89 — https://europepmc.org/article/PMC/PMC2666664
- [F25] Biodegradation of PHA: current state and future prospects, Front Microbiol 2025 — https://europepmc.org/article/PMC/PMC11893044
- [B25] Biodegradability of PHA biopolyesters in nature, Biodegradation 2025 — https://europepmc.org/article/PMC/PMC12339601
- [L16] Liu 2016 PhaJ1, Sci Rep — https://europepmc.org/article/PMC/PMC4823750
- [Liu15] Liu 2015 PhaZh1, AEM — https://pubmed.ncbi.nlm.nih.gov/25710370/
- [T22] Lihuaxuella thermophila PHB depolymerase 8DAJ, Protein Sci — https://europepmc.org/article/PMC/PMC9601781
- [L25] mclPHA lid, Appl Microbiol Biotechnol — https://europepmc.org/article/PMC/PMC12504323
- [B12] Brigham 2012 R. eutropha PhaZ, PLoS ONE — https://europepmc.org/article/PMC/PMC3430594
- [R26] R. rubrum PHB depolymerization, Microb Cell Fact — https://europepmc.org/article/PMC/PMC12983893
