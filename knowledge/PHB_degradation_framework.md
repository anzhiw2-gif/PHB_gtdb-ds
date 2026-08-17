# PHB 降解基因知识框架（v0.4，文献调研定稿）

> 状态：v0.4。2026-06-01 四路调研（PubMed + OpenAlex + Europe PMC +
> bioRxiv + Web）全部完成，待核实项已基本清零（遗留低优先项见报告第 9 节）。

## 1. 背景

- PHB（poly(3-hydroxybutyrate)，聚-3-羟基丁酸酯）是 PHA（polyhydroxyalkanoate，
  聚羟基链烷酸酯）家族中最常见的一员，是许多细菌在碳源过剩、营养受限时合成的
  胞内碳/能量储存聚合物（carbonosome 颗粒）。
- 经典生产菌：*Cupriavidus necator* H16（旧名 *Ralstonia eutropha*）、
  *Pseudomonas*、*Bacillus*、*Halomonas*、*Aeromonas*、*Azotobacter* 等。
- 合成途径（3 步）：PhaA（β-酮硫解酶）→ PhaB（乙酰乙酰-CoA 还原酶）→
  PhaC（PHA 合酶）。相关基因常成簇（phaCAB 操纵子）。
- 降解是"塑料生物降解/生物循环"领域的核心环节：PHB 是生物可降解塑料，
  其解聚酶（depolymerase）是酶法回收/降解的关键酶。

## 2. 降解基因目录（候选 gene set）

### 2.1 胞内 PHB 解聚酶（intracellular PHA depolymerase, i-PhaZ）
- 功能：内源颗粒降解/动员（mobilization），碳源耗竭时启动。
- 代表：*Cupriavidus necator* H16 的 PhaZ1（phaZ1, H16_A1435 附近）、PhaZ2、PhaZ3。
- 酶学：EC 3.1.1.75（poly(3-hydroxybutyrate) depolymerase）或相关；
  Ser-Asp-His 催化三联体，α/β-水解酶超家族；部分为 N-端锚定颗粒。

### 2.2 胞外 PHB 解聚酶（extracellular PHA depolymerase, e-PhaZ）
- 功能：水解胞外（环境中的）PHB 聚合物 → 寡聚体 → 单体（3HB/3HA）。
- 代表菌：*Pseudomonas lemoignei*（最早研究）、*Alcaligenes faecalis*、
  *Comamonas*、*Acidovorax*、*Bacillus*、*Streptomyces*、*Paucimonas lemoignei* 等。
- 结构域组织分类（Jendrossek 分类，【待核实】类型 I–IV）：
  - 催化结构域（α/β 水解酶折叠，含 G-X-S-X-G 脂肪酶盒五肽）
  - 连接肽（linker）
  - 底物结合结构域（SBD type 1 / type 2）
  - 信号肽（Sec 通路分泌，Gram 阴性菌）
- 胞外与胞内解聚酶序列相似性低，需分别建立 HMM/参考集。

### 2.3 寡聚体水解酶（oligomer hydrolase, OH）
- 功能：将 (3HB)n 寡聚体水解为单体；EC 3.1.1.22 或相关条目【待核实】。
- 代表：*Pseudomonas lemoignei* 的寡聚体水解酶、*Cupriavidus* 相关酶。

### 2.4 单体代谢（3HB → 乙酰-CoA）
- BdhA / D-3-羟基丁酸脱氢酶（EC 1.1.1.30）：3HB → 乙酰乙酸
- 乙酰乙酸-CoA 转移酶（ScoA/ScoB，EC 2.8.3.x）或乙酰乙酸-CoA 连接酶
  → 乙酰乙酰-CoA
- β-酮硫解酶（EC 2.3.1.9）→ 2 × 乙酰-CoA → TCA 循环
- 注意：与合成途径共用 β-酮硫解酶，筛选时需区分（方向性由底物/调控决定）。

### 2.5 其他相关
- PHA 合酶 PhaC（EC 2.3.1.-）：合成方向，但"降解基因簇"常与合成簇邻近，
  用于基因簇（biosynthetic gene cluster）背景分析。
- PhaP（phasin，颗粒结合蛋白）、PhaR（调控因子）：颗粒结构与调控。
- 细胞内颗粒降解还涉及 PCL/部分解聚酶的同源物（如 PHB depolymerase 家族内的
  结构域组合变异）。
- 部分菌（如 *Pseudomonas*）通过 β-氧化途径降解中长链 PHA 单体（3HA），
  涉及 fadD、fadE 等【待核实是否纳入】。

## 3. 已知分类/注释资源（【待核实】由子代理验证）

- Pfam/InterPro：PHB depolymerase 家族条目（【待核实】具体家族 ID；
  已知 α/β 水解酶相关 Pfam 有 Abhydrolase_1/3/6 等，cutinase 相关条目）
- KEGG：PHA 合成/降解通路、模块号（【待核实】M 模块号）
- eggNOG：PHA 降解相关 OG（【待核实】）
- 专用数据库：PAZy（Plastic Active Enzymes）、PlasticDB、
  PHA depolymerase database（【待核实】是否存在/可访问性）
- CAZy：cutinase/酯酶类可能部分覆盖

## 4. GTDB 侧要点（【待核实】由子代理验证）

- GTDB（Genome Taxonomy Database）：当前 release（R220? 【待核实】规模与日期）
- 可下载文件（【待核实】精确 URL 与版本）：
  - gtdb_genome_reps（代表性基因组，含细菌+古菌）
  - gtdb_proteins_aa_reps（每基因组 .faa 蛋白文件）
  - gtdb_metadata（基因组元数据：GTDB 分类、MIMAG 质量、来源生态等）
  - taxonomy 文件
- GTDB-Tk：本地分类/去冗余工具
- 分析策略候选：
  - 策略 A：用 GTDB 蛋白库（gtdb_proteins_aa_reps）直接做 HMM/DIAMOND 筛选
  - 策略 B：下载代表性基因组，Prodigal 预测后再筛（可保留基因上下文）
  - 生态分布：用 metadata 中的 GTDB 分类与来源信息做 phylum/ecosystem 分布统计

## 5. 待文献核实清单

- [x] Jendrossek 胞外 PHB 解聚酶分类（Annu Rev Microbiol 2002 已确认；
      type I–IV 细节见综述全文）【部分确认】
- [x] phasin_2 结构域 = Pfam **PF09361**（Azospirillum 综述确认）
- [ ] 胞内/胞外解聚酶 Pfam/InterPro 家族 ID（进行中）
- [ ] KEGG 通路/模块号（PHA 合成 M 模块、降解通路）（进行中）
- [x] PHA Depolymerase Engineering Database（PhaDED, 2009）存在（已确认；
      当前可访问性待 Web 核实）
- [ ] GTDB 最新 release 号、规模、蛋白质文件获取方式（进行中）
- [x] 已发表"大规模基因组/宏基因组筛选 PHB/PHA 解聚酶"研究：
      确认存在单菌/宏基因组规模工作，但无 GTDB 全库先例（新颖性确认）
- [ ] 是否有基于 GTDB 的 PHA 相关基因分布研究（进行中）

## 6. v0.2 新增核实内容（2026-06-01 文献调研）

### 6.1 PHB 动员通路（已核实，Azospirillum 综述）
PhaZ 解聚酶 → 3-羟基丁酸 → NAD(P) 依赖 BdhA 脱氢酶（四聚体）→ 乙酰乙酸
→ 乙酰乙酸-CoA 合成酶 → 乙酰乙酰-CoA → β-酮硫解酶 → 2 × 乙酰-CoA →
TCA/乙醛酸/β-氧化。注意：*A. brasilense* 中为"乙酰乙酸-CoA 合成酶"，
部分菌（*C. necator*）用乙酰乙酸-CoA 转移酶（ScoAB）——机制多样。

### 6.2 基因簇组织（已核实）
- *R. eutropha* / *A. brasilense*：phbCAB；*Azotobacter vinelandii*：
  phbBAC；*P. putida*：phaC1ZC2D + phaIF
- 基因可分布于染色体与质粒（多拷贝、跨复制子）
- 大多数 PHA 产生菌编码多个解聚酶同工酶

### 6.3 数据库（已核实存在）
- **PhaDED**（PHA Depolymerase Engineering Database, BMC Bioinformatics
  2009）：PHA 解聚酶系统分类（催化域分类 + 结构域组成），种子序列首选
- **PlasticEnz**（PLoS Comput Biol 2026）：同源 + ML 的塑料降解酶数据库
  与筛选工具（含 PHA 解聚酶），可作交叉验证
- 证据分级范式：demonstrated vs predicted（Martínez-Tobón 2018）

### 6.4 酶分类与结构域（v0.3，Europe PMC 全文核实）
- **e-PhaZ**：SP（信号肽）+ CD（催化域）+ LD（连接域，纤连蛋白 III 型）+
  SBD（C 端底物结合域）；催化三联体 Ser-Asp-His，Ser 在 G-X-S-X-G 脂肪酶盒；
  SBD 保守基序 sxxxHxxAGRa；两步机制（SBD 吸附→CD 水解）；endo/exo 型；
  P(3HB) 解聚酶 EC 3.1.1.75
- **i-PhaZ**：作用于天然颗粒，需颗粒表面蛋白激活；*R. rubrum* PhaZ1
  （周质）、PhaZ2（胞质主酶）、PhaZ3（胞质疑似无活性）；ApdA 激活蛋白、
  PhaR 调控子
- 非特异脂肪酶类也能水解 PHA（旁系同源干扰需注意）

### 6.5 分布规律（v0.3）
- 约 10% 微生物具 PHB 降解能力；e-PhaZ 见于 Acidovorax、Undibacterium、
  Janthinobacterium、Massilia、Duganella、Herbaspirillum、Burkholderia、
  Pseudomonas（mcl）、Ralstonia pickettii、Cupriavidus、Azotobacter、
  Azospirillum、Rhodococcus、Streptomyces 等
- scl/mcl-PHA 解聚酶驱动不同降解菌群落 → 分亚家族建模

### 6.6 最接近的前人研究（v0.3）
- **Viljakainen & Hug 2021**（Environ Microbiol, DOI 10.1111/1462-2920.15409）：
  3078 宏基因组中鉴定 13,869 个推定 PHA 解聚酶；废水系统频率最高、
  海洋与热泉最低；关键类群 Proteobacteria/Bacteroidota——方法蓝本；
  本项目差异化在 GTDB 统一分类框架 + PHB 降解基因目录系统化

### 6.7 v0.4 新增（Web 调研定稿）
- **种子/筛选资源**：
  - PhaDED 托管于 ESTHER（INRAE）：家族 Esterase_phb_PHAZ、
    PHAZ7_phb_depolymerase；HMM 可复用
  - PAZy（Buchholz 2022）：PHB 解聚酶种子 + API + DaRUS 镜像
  - PlasticDB（Gambarini 2022）：塑料降解酶库含 HMM
  - PlasticEnz（2026）：同源 + ML 集成筛选工具
  - NCBI CDD poly(3-hydroxyalkanoate) depolymerase 域模型；
    PDB 4BTV（PhaZ7-3HB）、8YNV（Bt PhaZ）
- **KEGG**：map00640（PHA 代谢通路）、M00843 模块、K03821=phaC；
  ⚠️ M00012 = 乙醛酸循环（非 PHA）
- **GTDB 数据**：最新 R11-RS232（R10 = 715,230 细菌 + 17,245 古菌）；
  下载 data.gtdb.ecogenomic.org/releases/latest/ 的
  gtdb_proteins_aa_reps.tar.gz + metadata；faa→分类学映射见 Forum 599
- **方法蓝本**：Zrimec 2021（mBio，HMMER v3.3 hmmsearch）；
  Microbial Genomics 2025（mgen.0.001814）；plastedma/PDETool 现成流程
- **CAZy**：无专设 PHA 解聚酶家族（以 ESTHER/PAZy/CDD 为准）
- **Pfam**：未检索到明确 PHB 解聚酶家族号（需官网复核）；phasin_2 =
  PF09361 已核实

### 6.8 古菌 PHB 降解基因（v0.5，2026-06-01 核实）

**结论：古菌确实存在 PHB 降解基因，但分两类，机制与细菌不同：**

| 类型 | 代表序列 | 细菌 HMM 覆盖 | 备注 |
|------|---------|-------------|------|
| **经典 PHB 解聚酶家族酯酶**（α/β 水解酶型） | Haladaptatus A0A1M6SXK9、Halomarina A0A6B0GML5、Haloarcula A0A8J7YBF6 等 13 条 | ✅ 12/13 被 ePhaZ/iPhaZ HMM 检出 | 属于经典家族，无需独立建模 |
| **patatin 样解聚酶（PhaZh1 型）** | **I3RBH0**（*Hfx mediterranei* = PhaZh1, 321aa）、**M1XPT2**（*Natronomonas moolapensis*, reviewed, 323aa） | ❌ 0 检出 | 类磷脂酶折叠（Pfam Patatin），须独立建模 |

**文献依据（PubMed）**：
- PMID [25710370](https://pubmed.ncbi.nlm.nih.gov/25710370/)（AEM 2015）：
  *Hfx mediterranei* **PhaZh1**（patatin 样 PHA 解聚酶），定位于颗粒，
  水解 nPHB/nPHBV 产 3HB 单体；催化关键 Ser47（G-X-S47-X-G 脂肪酶盒）
  + Gly16 + Asp195；**与 bdhA 成簇（HFX_6463-6464）**——与细菌
  phaZ+bdhA 结构一致
- PMID [27052994](https://pubmed.ncbi.nlm.nih.gov/27052994/)（Sci Rep 2016）：
  *Hfx mediterranei* **PhaJ1**（颗粒相关烯酰-CoA 水合酶，R-ECH）参与
  PHA 动员（(R)-3-羟基酰基-CoA 脱水）
- PMID [22247127](https://pubmed.ncbi.nlm.nih.gov/22247127/)：古菌
  phasin（PhaP）——已在种子中（I3R9Z2）
- PMID [27098259](https://pubmed.ncbi.nlm.nih.gov/27098259/)：古菌产
  PHBV 的环境生物降解（活性污泥）

**对筛选的意义（已实施）**：
1. patatin 型种子（2 条验证序列）+ Pfam "Patatin" HMM（hmmfetch 自
   Pfam-A.hmm 提取）→ 第 6 个筛选家族 **patatin**（06_screen.sh 已加入）
2. 验证规则：patatin 催化 Ser-Asp 二元组（非三联体），长度 200-500aa；
   细菌中 patatin 命中多为磷脂酶，需 PHA 基因簇上下文二次确认
3. 古菌经典家族酯酶（12/13）会被 ePhaZ/iPhaZ HMM 自动捕获
4. 全库筛选覆盖全部 199,923 基因组（含 ~10,000+ 古菌 reps），
   最终报告对古菌单独统计
