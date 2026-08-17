# PHA/PHB 降解酶家族分类体系提炼（用于 GTDB 全库 HMM 筛选）

> 本文档基于 Europe PMC 开放获取全文提炼，覆盖：胞外 PHA 解聚酶、胞内 PHA 解聚酶、寡聚体水解酶、古菌 PHA 降解酶，以及各家族保守特征（催化三联体/二元组、lipase box、长度、结构域组成），供构建 GTDB 全库筛选的家族 HMM 体系使用。
>
> **核心文献编号（后续正文以 [编号] 引用）：**
> - **[K09]** Knoll M, Hamm TM, Wagner F, Martinez V, Pleiss J. The PHA Depolymerase Engineering Database (PhaDED). BMC Bioinformatics. 2009;10:89. PMCID: PMC2666664; PMID: 19296857; DOI: 10.1186/1471-2105-10-89
> - **[JH02]** Jendrossek D, Handrick R. Microbial degradation of polyhydroxyalkanoates. Annu Rev Microbiol. 2002;56:403-432. PMID: 12213937; DOI: 10.1146/annurev.micro.56.012302.160838（无 OA 全文，其分类框架经 [F25][K09] 转引确认）
> - **[F25]** Biodegradation of polyhydroxyalkanoates: current state and future prospects. Front Microbiol. 2025;16:1542468. PMCID: PMC11893044; PMID: 40066265; DOI: 10.3389/fmicb.2025.1542468
> - **[B25]** Biodegradability of polyhydroxyalkanoate (PHA) biopolyesters in nature: a review. Biodegradation. 2025. PMCID: PMC12339601; PMID: 40788578; DOI: 10.1007/s10532-025-10164-y（本地全文 ft_biodegradability_pha_review.txt）
> - **[R26]** Gene expression analysis reveals distinct PHB depolymerization mechanisms…in Rhodospirillum rubrum. Microb Cell Fact. 2026. PMCID: PMC12983893; PMID: 41664065; DOI: 10.1186/s12934-026-02946-7
> - **[L25]** Revealing the essential role of the lid in mclPHA intracellular depolymerase from Pseudomonas putida KT2440. Appl Microbiol Biotechnol. 2025. PMCID: PMC12504323; PMID: 41055782; DOI: 10.1007/s00253-025-13605-z
> - **[T22]** Bioplastic degradation by a polyhydroxybutyrate depolymerase from a thermophilic soil bacterium (Lihuaxuella thermophila). Protein Sci. 2022. PMCID: PMC9601781; PMID: 36222314; DOI: 10.1002/pro.4470（PDB: 8DAJ）
> - **[K19]** Polyhydroxyalkanoate Biosynthesis at the Edge of Water Activity—Haloarchaea as Biopolyester Factories. Bioengineering. 2019;6(2):34. PMCID: PMC6631277; PMID: 30995811; DOI: 10.3390/bioengineering6020034
> - **[W19]** Bioinformatics Analysis of Metabolism Pathways of Archaeal Energy Reserves. 2019. PMCID: PMC6355812; PMID: 30705313
> - **[L16]** Liu G, et al. Enoyl-CoA hydratase mediates polyhydroxyalkanoate mobilization in Haloferax mediterranei. Sci Rep. 2016;6:24015. PMCID: PMC4823750; PMID: 27052994; DOI: 10.1038/srep24015
> - **[Liu15]** Liu G, Hou J, Cai S, Zhao D, Cai L, Han J, Zhou J, Xiang H. A patatin-like protein associated with the polyhydroxyalkanoate (PHA) granules of Haloferax mediterranei acts as an efficient depolymerase in the degradation of native PHA. Appl Environ Microbiol. 2015;81(9):3029-3038. PMID: 25710370; DOI: 10.1128/AEM.04269-14（PhaZh1 原始论文；无 OA 全文，经 [L16][K19] 转引）
> - **[B12]** Brigham CJ, et al. Examination of PHB Depolymerases in Ralstonia eutropha: Further Elucidation of the Roles of Enzymes in PHB Homeostasis. PLoS ONE. 2012. PMCID: PMC3430594; PMID: 22537946
> - **[A05]** Abe T, Kobayashi T, Saito T. Properties of a novel intracellular poly(3-hydroxybutyrate) depolymerase with high specific activity (PhaZd) in Wautersia eutropha H16. J Bacteriol. 2005;187:6982-6990. PMCID: PMC1251622; PMID: 16199568
> - **[S01]** Saegusa H, et al. Cloning of an intracellular Poly[D(-)-3-Hydroxybutyrate] depolymerase gene from Ralstonia eutropha H16. J Bacteriol. 2001;183:94-100. PMCID: PMC94854; PMID: 11114905
> - **[H04]** Handrick R, et al. The "intracellular" poly(3-hydroxybutyrate) (PHB) depolymerase of Rhodospirillum rubrum is a periplasm-located protein… J Bacteriol. 2004;186:7243-7253. PMCID: PMC523223; PMID: 15489436

---

## 0. 总体框架：PHA 解聚酶的两维分类逻辑

PHA 解聚酶（EC 3.1.1.75 / EC 3.1.1.76）均属 **α/β-水解酶折叠（α/β-hydrolase fold）** 家族，催化丝氨酸-组氨酸-天冬氨酸（Ser-His-Asp）三元组（个别家族为 Cys-His-Asp），催化 Ser 位于 **Gx₁Sx₂G lipase box** 基序中 [K09][F25]。

分类按两个维度展开（[K09][F25]）：

1. **按底物/定位**（Knoll/PhaDED 八大超家族的基础）：
   - 胞内（i-）降解**天然(native)颗粒**（nPHASCL / nPHAMCL）
   - 胞外（e-）降解**变性(denatured)胞外颗粒**（dPHASCL / dPHAMCL）
   - 周质（periplasmatic）——R. rubrum 特例
   - 例外：胞外但只作用于天然颗粒的 nPHASCL（Paucimonas lemoignei PhaZ7）
2. **按底物链长**：SCL（3-5 C，如 PHB/PHV）vs MCL（6-15 C，如 P(3HO)）

Pfam 对应条目：[K09] 指出 PHA 解聚酶对应 **PF10503**（Esterase PHB depolymerase，即酯酶型 PHB 解聚酶家族）与 **PF06850**（bacterial PHB depolymerase C 端，即底物结合域）。

---

## 1. 胞外 PHA 解聚酶（e-PHA depolymerase）分类

### 1.1 八大超家族（PhaDED，8 superfamilies / 38 homologous families，587 个蛋白）[K09]

| # | 超家族（PhaDED） | 底物/定位 | 催化特征 | 代表（种子序列，gi） |
|---|---|---|---|---|
| 1 | **胞内 nPHASCL（无 lipase box）** | 天然 SCL 颗粒 | Cys-His-Asp 三元组（催化 Cys 替代 Ser），无 GxSxG | R. eutropha H16 胞内 PHB 解聚酶 PhaZ1（gi 3641686 [S01]） |
| 2 | **胞内 nPHASCL（有 lipase box）** | 天然 SCL 颗粒 | Ser-His-Asp；Gx₁Sx₂G 中 x₁=Trp | Bacillus thuringiensis serovar israelensis（gi 75763431） |
| 3 | **周质 PHA 解聚酶** | 周质、天然 PHB | Ser-His-Asp；催化域 type 2；x₁=Ile | Rhodospirillum rubrum（gi 22035160 [H04]） |
| 4 | **胞内 nPHAMCL** | 天然 MCL 颗粒 | Ser-His-Asp；x₁=Val | Pseudomonas oleovorans（gi 130002）、P. putida（gi 21689574） |
| 5 | **胞外 dPHASCL（催化域 type 1）** | 变性 SCL 颗粒 | Ser-His-Asp；oxyanion hole 在 lipase box N 端（类脂酶）；x₁=Leu/Ile（81% 疏水） | Alcaligenes faecalis（gi 1777951）、Ralstonia pickettii T1（gi 130019）、Bacillus megaterium（gi 116744367）、Paucimonas lemoignei 多酶（gi 1730532/7385117/1657610/1621355/531464/531466）、Pseudomonas stutzeri（gi 75538924） |
| 6 | **胞外 dPHASCL（催化域 type 2）** | 变性 SCL 颗粒 | Ser-His-Asp；oxyanion hole 在催化三元组 C 端；x₁ 全部疏水 | Acidovorax sp. TP4（gi 4033618）、Caldimonas manganoxidans（gi 7209864）、Comamonas sp.（gi 565666）、Delftia acidovorans（gi 75340123）、Schlegelella sp. KB1a（gi 47078657）、Streptomyces exfoliatus（gi 1389770）、Penicillium funiculosum（gi 88192747，PDB 2D80） |
| 7 | **胞外 nPHASCL（特殊）** | 仅天然 SCL 颗粒 | 无典型 Gx₁Sx₂G，改为 **AHSMG** 基序；x₁ 无疏水残基；无底物结合域 | Paucimonas lemoignei PhaZ7（gi 15788987，PDB 2VTV） |
| 8 | **胞外 dPHAMCL** | 变性 MCL 颗粒 | Ser-His-Asp；x₁=Ile；无底物结合域（N 端充当结合位点） | Pseudomonas alcaligenes（gi 34452163/29470160）、P. fluorescens GK13（gi 21542177） |

> 注释：PhaDED 中最大两个家族为 #1（224 条，38%）与 #5（234 条，39%）[K09]。家族内细分见 5.1 节保守基序。

### 1.2 催化域类型（type 1 / type 2）[K09][F25]

- **催化域 type 1**：oxyanion hole 位于 **lipase box（GxSxG）N 端**，与经典脂酶一致。
- **催化域 type 2**：oxyanion hole 位于 **催化三元组 C 端**。
- [F25] 另以"lipase box 在序列中部（Type I）vs 近 N 端（Type II）"描述 SCL 解聚酶催化域位置差异（引用 JH02 体系）；两者口径不同，**建模时以 PhaDED 的 type 1/2（oxyanion hole 相对位置）为准**，并将 [F25] 的表述标注为补充。Bacillus 来源 e-PHASCL 为 Type I（lipase box 居中，如 B. megaterium N-18-25-9）；Comamonas acidovorans YM1609 与 Leptothrix sp. HS 为 Type II [F25]。

### 1.3 结构域架构（胞外 dPHASCL 经典四域模型）[F25][K09]

信号肽（N 端）→ 催化域（含 lipase box + oxyanion hole）→ 连接域（linker）→ C 端底物结合域（SBD）

- **信号肽**：约 20-30 aa，Sec 分泌信号（例：LtPHBase 322 aa 全长的 N 端 22 aa 为信号肽 [T22]）。
- **连接域三类** [F25]：**Fn3**（Fibronectin type III；Streptomyces sp. SFB5A、Comamonas acidovorans YM1609、Leptothrix sp. HS）、**Thr-rich**（P. lemoignei）、**Cad**（Cadherin-like；Marinobacter sp. NK-1）。Pseudomonas stutzeri 无连接域。
- **底物结合域两类：SBD1 与 SBD2** [F25]（引用 JH02）。SBD 吸附 PHA 表面，独立于催化域发挥作用（例：A. faecalis 解聚酶可吸附 5 种底物但仅水解 PHB/PHP/P(4HB) 3 种 [F25]）。
- **dPHAMCL 无 SBD 与连接域**：N 端区域兼作底物结合位点 [K09][F25]。
- **e-PHAMCL 结构**：信号肽 + N 端底物结合域 + C 端催化域，277-282 aa [F25]。

### 1.4 胞外家族代表酶的生化特征 [F25]

- **e-PHASCL（EC 3.1.1.75）**：分子量 35-63 kDa；pH 最适宽；温度最适 30-80 °C（Schlegelella thermodepolymerans 至 90 °C；Streptomyces sp. IN1 80 °C/15 min 稳定且耐 pH 12）；受 DTT/β-ME 抑制（依赖二硫键）；受 PMSF/DFP（丝氨酸水解酶抑制剂）与 EDTA 抑制；Tween/Triton/SDS 抑制；Mn²⁺/Fe²⁺/Ni²⁺ 抑制，Ca²⁺/Mg²⁺/Na⁺/K⁺ 激活。
- **e-PHAMCL（EC 3.1.1.76）**：分子量 25-30 kDa（Xanthomonas sp. JS02 41.7 kDa 例外）；pH 8.0-10.0；温度 30-70 °C（T. thermophilus HB8 最耐热）；Pseudomonas 型不受 DTT 影响，Streptomyces/Bdellovibrio 型受 DTT 抑制；大多受 PMSF 抑制（Streptomyces 例外）；主要为疏水氨基酸（芳香族+非极性脂肪族），序列相似度 69-98%。

---

## 2. 胞内 PHA 解聚酶（i-PHA depolymerase）分类

### 2.1 两大催化类型 [K09]

1. **i-nPHASCL（无 lipase box，Cys 型）**——PhaDED 中最大的胞内家族（224 条）。催化三元组 **Cys-His-Asp**，催化 Cys 前的残基几乎全为疏水 Val。代表：
   - **R. eutropha（Cupriavidus necator）PhaZ1**（H16_A1150；gi CAJ92291.1）——模型胞内解聚酶，位于颗粒表面，将 PHB 链上的 3HB-CoA 硫解为单体 [B12][S01]
   - **R. eutropha PhaZ2**（H16_A2862）——胞内降解 PHB，效率高于 PhaZ1，还影响颗粒密度/三维结构 [B12]
   - **R. eutropha PhaZ3 / PhaZ5**（H16_B0339 / H16_B1014）——推定的胞内解聚酶，参与颗粒重塑而非主要动员 [B12]
   - **PhaZd**（Abe 2005，即 PhaZ5 早期命名或独立条目；高比活胞内 PHB 解聚酶）[A05][B12]
2. **i-nPHASCL（有 lipase box，Ser 型）**——Ser-His-Asp；Gx₁Sx₂G 中 x₁=Trp。代表：**Bacillus thuringiensis** 胞内 PHB 解聚酶（gi 75763431）[K09]

### 2.2 胞内 MCL 解聚酶（i-nPHAMCL）

- PhaDED 超家族 #4；Gx₁Sx₂G 中 x₁=Val（几乎全部）[K09]。
- **Pseudomonas 专属特征**：胞内 mcl-PHA 降解为 Pseudomonas 属特异性状 [L25]。
- **模型酶 PhaZKT（P. putida KT2440）**：α/β-水解酶折叠 + **lid 结构域**（类脂酶），催化三元组 **Ser102-Asp221-His248**；lid 是活性必需（删除 FNGIG 环 aa34-38 或 YYWQLF 环 aa190-195 完全失活）；S184F（lid 铰链）改变底物偏好，G286R 提升 mcl-PHA 解聚活性 [L25]。**与胞外 mcl 解聚酶不同（胞外无 lid）**——可作为区分胞内/胞外 mcl 酶的结构标记。
- 胞外 mcl 对照：P. solani（原 P. fluorescens）GK13 PhaZGK13 催化三元组 Ser172-Asp228-His260，无 lid [L25]。

### 2.3 周质 PHB 解聚酶（特殊）[K09][H04][R26]

- **Rhodospirillum rubrum PhaZ1**（gi 22035160）：定位于周质但降解天然 PHB；催化域类胞外酶（type 2）[H04]。
- R. rubrum 共 3 个 PhaZ：**PhaZ1（周质）、PhaZ2（胞质，主要胞内解聚酶）、PhaZ3（胞质，体内可能失活）**，另需 **ApdA**（phasin 样颗粒动员激活因子）激活颗粒供 PhaZ2 作用 [R26]。
- R. rubrum 的 phaZ1 在乙酸培养的降解期高表达（推测作用于膜相关 PHB），phaZ2+apdA 在果糖培养中上调 [R26]。

### 2.4 R. eutropha 全基因组 PhaZ 清单（含 PhaDED 归类）[B12] Table 1

| 基因 | 功能 | PhaDED 超家族 | 蛋白 ID | 位点标签 |
|---|---|---|---|---|
| phaZ1 | 胞内 PHB 解聚酶 | i-nPHAscl（无 lipase box） | CAJ92291.1 | H16_A1150 |
| phaZ2 | 胞内 PHB 解聚酶 | i-nPHAscl（无 lipase box） | CAJ93939.1 | H16_A2862 |
| phaZ3 | 推定胞内 PHB 解聚酶 | i-nPHAscl（无 lipase box） | CAJ95139.1 | H16_B0339 |
| phaZ4 | 推定 PHB 解聚酶 | i-nPHAscl（无 lipase box） | AAP85930.1 | PHG178 |
| phaZ5 | 胞内 PHB 解聚酶 | i-nPHAscl（无 lipase box） | CAJ95805.1 | H16_B1014 |
| phaZ6 | PHB 解聚酶 | e-dPHAscl（催化域 type 1） | CAJ96855.1 | H16_B2073 |
| phaZ7 | PHB 解聚酶 | e-dPHAscl（催化域 type 1） | CAJ97183.1 | H16_B2401 |
| phaY1 | D-(-)-3HB 寡聚体水解酶 | 不在 PhaDED | CAJ93348.1 | H16_A2251 |
| phaY2 | D-(-)-3HB 寡聚体水解酶 | 不在 PhaDED | CAJ92475.1 | H16_A1335 |

> 启示：**R. eutropha 一个基因组内同时含 5 个胞内 nPHASCL（no-lipase-box）型**，说明该家族在基因组中高度冗余；GTDB 全库筛选时同一基因组的多个拷贝均需捕获。

---

## 3. 寡聚体水解酶（oligomer hydrolase）

- **定位**：PHA 降解途径中，解聚酶将聚合物水解为**寡聚体、二聚体与单体**（寡聚体包括 dimer/trimer 等）[F25]。寡聚体水解酶负责将低聚物进一步水解为单体（如 D-(-)-3-羟基丁酸）。
- **EC 号**：3HB 寡聚体水解酶对应 **EC 3.1.1.22**（3-hydroxybutyrate-oligomer hydrolase 类，如 [B12] 中 PhaY1/PhaY2 的注释功能）；胞外主解聚酶为 EC 3.1.1.75（SCL）/3.1.1.76（MCL）[F25][B25]。
- **代表酶**：
  - **R. eutropha PhaY1（H16_A2251）与 PhaY2（H16_A1335）**：注释为 D-(-)-3-hydroxybutyrate oligomer hydrolase；**不在 PhaDED 的 8 大解聚酶超家族内**（N/A）[B12]。
  - 部分胞外解聚酶直接产生单体（DP1，如 S. exfoliatus K10、多数 e-PHASCL），另一些产生二聚体（DP2，如 P. fluorescens GK13、Streptomyces sp. KJ-72、Bdellovibrio bacteriovorus HD100）或单体+二聚体+三聚体混合物（Guo 等；Sadocco 等；Comamonas；Uefuji 等；Blevins 等）[F25]。
- **胞内途径衔接**：胞内降解产物 (R)-3HB-CoA/3HB 进入 **3HB 脱氢酶 + 乙酰乙酸-CoA 连接酶** 途径（R. rubrum 中为 3HB dehydrogenase KUL73_05510 与 AACoA synthase）[R26]。
- **检索备注**：Europe PMC OA 中未检索到专门针对 PHA 寡聚体水解酶分类的综述全文；该类别以 [B12] 的 PhaY 注释与 [F25] 的产物描述为准。EC 3.1.1.22 的分子细节建议后续从 BRENDA/UniProt 与 Saito 组原始文献补充。

---

## 4. 古菌 PHA 降解酶类型

### 4.1 古菌 PHA 代谢总体格局 [W19][L16][K19]

- 427 个古菌参考蛋白组的 HMM 通路筛选证实：古菌含 polyP、PHA、糖原储能途径，**无 TAG/蜡酯**；PHA 代谢与**嗜盐古菌**紧密相关（完整合成途径 51/427 种，含 PhaA/B/C/E/P；降解酶 PhaZ）[W19]。
- 古菌 PHA 合成为 **Class III 型合成酶（PhaC+PhaE 复合体）**，仅产 scl-PHA（PHB/PHBV）[K19][W19]。
- 卤古菌（class Halobacteria，103 个测序物种）：50% 含 phaJ，51% 含 phaC；96% 含 phaJ 的物种同时具备完整 β-氧化四酶组（酰基-CoA 脱氢酶、烯酰-CoA 水合酶、3-羟酰-CoA 脱氢酶、3-酮酰-CoA 硫解酶）[L16]。

### 4.2 patatin 型解聚酶 PhaZh1（古菌特有的颗粒结合解聚酶）[Liu15][L16][K19]

- **PhaZh1**：Haloferax mediterranei 中与 PHA 颗粒结合（PGAP，granule-associated protein）的 **patatin 样蛋白**，是最早鉴定的参与古菌天然颗粒动员（解聚）的酶；体外高效将 PHB(V) 降解为 3-羟基丁酸，是 nPHA 颗粒水解的关键酶 [Liu15][K19][L16]。
- **体内角色有限**：phaZh1 敲除对 Hfx. mediterranei 胞内 PHA 动员无显著影响 → 存在更有效的替代动员途径（即 PhaJ1/β-氧化途径）[L16]。
- 功能注释要点：patatin 样结构域（磷脂酶 A₂ 折叠，α/β-水解酶超家族），与细菌"经典家族酯酶/解聚酶"（GxSxG Ser 型）序列相似性低——**HMM 建模时必须作为独立家族**，不能并入细菌 e-PHASCL。

### 4.3 PhaJ（R 特异性烯酰-CoA 水合酶，MaoC 家族）[L16]

- Hfx. mediterranei 有 5 个推定 R-ECH（PhaJ1-PhaJ5；HFX_1483/2901/5217/6361/6433），均含 **MaoC 样结构域（Pfam PF01575）**，催化二元组 Asp-His（类 PhaJAc 的 Asp31-His36）。
- **仅 PhaJ1（HFX_5217，219 aa）与颗粒结合**（基因簇 HFX_5217-phaR-phaP-phaE-phaC），并主导 PHA 动员：催化 (R)-3-羟基酰基-CoA **脱水** → 烯酰-CoA（可逆反应），将 PHA 动员与 **β-氧化** 衔接；phaJ1 缺失使动员显著下降，回补恢复（62.3% vs 10.6% 降解率）[L16]。
- 细菌中 PhaJ（Aeromonas caviae PhaJAc、Pseudomonas aeruginosa PhaJ1Pa-PhaJ4Pa、R. eutropha 16 个 R-ECH 直系同源）主要用于 PHA **合成**供单体 [L16]；古菌中则主要参与**降解**——建模时注意同一 MaoC 家族跨合成/降解双功能。

### 4.4 古菌"经典家族酯酶"与颗粒蛋白体系 [K19][L16]

- Hfx. mediterranei 的 6 个 PGAPs：PhaJ1（HFX_5217，R-ECH）、PhaP（phasin）、PhaR（调控蛋白，兼转录因子+颗粒结合）、PhaE/PhaC（Class III 合成酶亚基）、**PhaZh1（patatin 解聚酶）** [L16]。
- 目前文献未在古菌中系统报道细菌型胞外 GxSxG 丝氨酸解聚酶；古菌 PHA 降解主要由 **PhaZh1（patatin）+ PhaJ1（MaoC R-ECH → β-氧化）** 双途径完成 [L16]。
- 另注意古菌存在**类 WS/DGAT（蜡酯/甘油三酯合成酶）同源物**但多数缺失 [W19]，与 PHA 降解无直接关系。

---

## 5. 各家族保守特征（HMM 建模与验证规则）

### 5.1 催化特征总表

| 家族 | 催化残基 | lipase box / 关键基序 | 保守特征氨基酸 | 长度/结构域 | 代表 |
|---|---|---|---|---|---|
| e-dPHASCL type 1 | Ser-His-Asp | Gx₁Sx₂G，x₁=Leu/Ile（81% 疏水），x₂=Ala/Ser | oxyanion hole N 端于 lipase box；SBD1/SBD2 | 信号肽+催化域+linker(Fn3/Thr-rich/Cad)+SBD；35-63 kDa | A. faecalis；R. pickettii T1；B. megaterium |
| e-dPHASCL type 2 | Ser-His-Asp | Gx₁Sx₂G，x₁ 全部疏水 | oxyanion hole C 端于催化三元组 | 同上 | Acidovorax TP4；P. funiculosum (2D80) |
| e-dPHAMCL | Ser-His-Asp | Gx₁Sx₂G，x₁=Ile，x₂=Ser（变异：T. thermophilus Gly/Tyr；S. exfoliatus His/Gln） | oxyanion hole 残基多变（His/Ser/Asn/Gln）；无 SBD、无 linker | 信号肽+N 端结合域+C 端催化域；25-30 kDa（277-282 aa） | P. fluorescens GK13；P. alcaligenes；T. thermophilus |
| 胞外 nPHASCL（PhaZ7 型） | Ser-His-Asp | **AHSMG**（非 GxSxG） | 无 SBD；活性位点被 14 aa 移动 lid 遮蔽（5MIX/5MIY 中 aa281-295） | 信号肽+催化域；~17% 同源于 LtPHBase | P. lemoignei PhaZ7（2VTV） |
| i-nPHASCL（无 lipase box） | **Cys-His-Asp** | 无 GxSxG；催化 Cys-1 位为疏水 Val | 与 α/β-水解酶同折叠但 Cys 催化 | 全胞内，颗粒结合 | R. eutropha PhaZ1/Z2/Z3/Z5、PhaZd |
| i-nPHASCL（有 lipase box） | Ser-His-Asp | Gx₁Sx₂G，x₁=Trp | — | 胞内 | B. thuringiensis |
| i-nPHAMCL | Ser-His-Asp | Gx₁Sx₂G，x₁=Val | **lid 结构域**（类脂酶，PhaZKT aa~34-38 与 aa~190-195 区） | 胞内；PhaZKT 三元组 Ser102-Asp221-His248 | P. putida KT2440 PhaZKT |
| 周质 | Ser-His-Asp | Gx₁Sx₂G，x₁=Ile；催化域 type 2 | 类胞外酶 | 周质 | R. rubrum PhaZ1 |
| 古菌 patatin 型（PhaZh1） | 推测 Ser-Asp（patatin 折叠催化二元组/三元组，见文献 [Liu15]） | 无细菌型 GxSxG | patatin/磷脂酶 A₂ 折叠 | 颗粒结合（PGAP） | Hfx. mediterranei PhaZh1 |
| 古菌 MaoC R-ECH（PhaJ1） | Asp-His 催化二元组 | MaoC 结构域（PF01575） | R 立体专一性；(R)-3HA-CoA ↔ 烯酰-CoA | ~219 aa（PhaJ1） | Hfx. mediterranei PhaJ1 |
| 寡聚体水解酶（PhaY） | 见 EC 3.1.1.22 注释 | 不在 PhaDED | D-(-)-3HB 寡聚体 → 单体 | — | R. eutropha PhaY1/PhaY2 |

### 5.2 e-PHASCL 解聚酶四个保守基序（HMM 验证核心规则）[T22]

以 Lihuaxuella thermophila LtPHBase（PDB 8DAJ，1.2 Å；322 aa，含 22 aa 信号肽）为基准：

1. **Ser 基序**：`IDXXXXYVXGLSXGG`（LtPHBase aa109-124：SDDSRRVYAAGLSAGG）——含 lipase box 变体 GLSXGG
2. **Asp 基序**：`GXXDYTV`（aa194-200：GTSDYTV）
3. **His 基序**：`GMXHXXPXXG`（aa267-274：GMGHAWSSG）
4. **oxyanion hole 基序**：`HGCXQ`（aa38-42：HGCTQ）——Cys 型 oxyanion hole

催化三元组：Ser121-His270-Asp197；oxyanion hole Cys40；中央 10 股 β-折叠，两侧 α-螺旋。DALI 结构比对指向酯酶/蛋白酶家族（而非脂酶）[T22]。

> **验证规则建议**：对任一候选 e-PHASCL，应同时检出上述 4 基序且相对位置一致（Ser 基序 ~aa109、Asp ~aa197、His ~aa270、oxyanion ~aa40，LtPHBase 编号），并核对疏水 x₁ 特征（区分于脂酶/酯酶）[K09][T22]。

### 5.3 lipase box 的疏水 x₁ 是 PHA 解聚酶区分于脂酶/酯酶的关键 [K09]

- 脂酶/酯酶 Gx₁Sx₂G 的 x₁ 通常为极性残基；PHA 解聚酶几乎全部为**疏水残基**（SCL 胞外：Leu/Ile；i-nPHASCL-lipase-box：Trp；周质：Ile；i-nPHAMCL：Val；dPHAMCL：Ile）[K09]。
- 例外家族（胞外 nPHASCL，PhaZ7 型）为 AHSMG，与 LED 中 Bacillus 脂酶家族 abH18.01 共享该基序 [K09]。

### 5.4 HMM 建模建议（依据以上文献）

1. **按 PhaDED 八大超家族分别建 profile HMM**（PhaDED 本身即为每个 family/superfamily 提供 HMMER profile，可作初始种子；[K09] 明确说明其用途即"从完整基因组中 in silico 鉴定与分类"）。
2. **必须单独建模的古菌家族**：patatin（PhaZh1）、MaoC R-ECH（PhaJ1，PF01575）——两者与细菌经典家族序列相似性极低 [L16]。
3. **注意跨类别同源**：(a) R. eutropha 的 PhaZ6/PhaZ7 属 e-dPHASCL type 1 但为胞内基因（GenBank 注释与序列归属不一致）；(b) 某些注释为"胞内"的序列（C. taiwanensis、R. eutropha gi:74267419）按序列聚到胞外 dPHASCL type 1；(c) 某 Pseudomonas 序列注释为"胞外"但聚到胞内 nPHAMCL [K09]——**注释不可靠，必须以序列聚类为准**。
4. **Pfam 交叉验证**：PF10503（酯酶型 PHB 解聚酶）+ PF06850（PHB 解聚酶 C 端 SBD）[K09]；MaoC = PF01575 [L16]；α/β-水解酶折叠（如 PF12697/CL0028 等）为全局背景。
5. **阈值与验证**：用 [T22] 的四基序 + [K09] 的 x₁ 疏水性作为 HMM 命中的二级验证；对胞外酶验证信号肽/分泌信号；对胞内酶验证无信号肽且颗粒结合（或用 lid 结构域区分胞内 mcl 酶 [L25]）。

---

## 6. 用到的论文 URL 清单

1. Knoll et al. 2009 (PhaDED): https://europepmc.org/article/PMC/PMC2666664 （DOI: 10.1186/1471-2105-10-89）
2. Jendrossek & Handrick 2002: https://europepmc.org/article/MED/12213937 （DOI: 10.1146/annurev.micro.56.012302.160838；无 OA 全文）
3. Biodegradability of PHA biopolyesters in nature (2025): https://europepmc.org/article/PMC/PMC12339601 （DOI: 10.1007/s10532-025-10164-y）
4. Biodegradation of PHA: current state and future prospects (2025): https://europepmc.org/article/PMC/PMC11893044 （DOI: 10.3389/fmicb.2025.1542468）
5. R. rubrum PHB depolymerization mechanisms (2026): https://europepmc.org/article/PMC/PMC12983893 （DOI: 10.1186/s12934-026-02946-7）
6. Lid in mclPHA intracellular depolymerase P. putida KT2440 (2025): https://europepmc.org/article/PMC/PMC12504323 （DOI: 10.1007/s00253-025-13605-z）
7. Lihuaxuella thermophila PHB depolymerase (2022): https://europepmc.org/article/PMC/PMC9601781 （DOI: 10.1002/pro.4470；PDB 8DAJ）
8. Haloarchaea as Biopolyester Factories (2019): https://europepmc.org/article/PMC/PMC6631277 （DOI: 10.3390/bioengineering6020034）
9. Archaeal Energy Reserves bioinformatics (2019): https://europepmc.org/article/PMC/PMC6355812
10. Enoyl-CoA hydratase mediates PHA mobilization in H. mediterranei (2016): https://europepmc.org/article/PMC/PMC4823750 （DOI: 10.1038/srep24015）
11. Liu et al. 2015 (PhaZh1, patatin-like): https://pubmed.ncbi.nlm.nih.gov/25710370/（AEM 81:3029-3038；DOI: 10.1128/AEM.04269-14；无 OA 全文，经 [L16][K19] 转引）
12. Brigham et al. 2012 (R. eutropha PhaZ): https://europepmc.org/article/PMC/PMC3430594
13. Abe et al. 2005 (PhaZd): https://europepmc.org/article/PMC/PMC1251622
14. Saegusa et al. 2001 (R. eutropha 胞内 PHB 解聚酶): https://europepmc.org/article/PMC/PMC94854
15. Handrick et al. 2004 (R. rubrum 周质 PHB 解聚酶): https://europepmc.org/article/PMC/PMC523223

---

*输出文件：D:\PHB_gtdb-ds\research\europepmc\comp_classification.md*
*配套原始全文：comp_ft_knoll_phaded.txt/xml、comp_ft_biodeg_current.txt、comp_ft_phb_mech_2026.txt、comp_ft_mclpha_lid.txt、comp_ft_thermophile_phb.txt、comp_ft_haloarchaea_pha.txt、comp_ft_archaeal_energy.txt、comp_ft_enoyl_hfx.txt、comp_ft_reutropha_phaz.txt/xml、ft_biodegradability_pha_review.txt（本地既有）*
