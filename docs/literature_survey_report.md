# GTDB PHB 降解基因系统生信分析 — 文献与项目调研报告

> 版本：v1.0（4 路调研全部完成并整合，2026-06-01）
> 目标：为"基于 GTDB 数据库的 PHB 降解基因系统生信分析"提供文献与方法学基础。
> 数据来源：PubMed（154 篇去重）、OpenAlex（24 篇关键）、Europe PMC
> （13 篇高相关 + 4 篇全文）、bioRxiv（7 篇预印本）、Web 调研（数据库/
> 工具/GTDB 数据核实）。原始数据见 research/ 各子目录。

---

## 目录

1. 项目背景与目标
2. 文献调研方法与检索策略
3. PHB 降解生物学与基因目录（文献综述）
   3.1 PHB 代谢总览（合成-储存-降解）
   3.2 胞内 PHA 解聚酶（i-PhaZ）
   3.3 胞外 PHA 解聚酶（e-PhaZ）与结构域分类
   3.4 寡聚体水解酶与 3HB 单体代谢
   3.5 基因簇组织
   3.6 生态分布
   3.7 基因组/宏基因组层面的筛选研究
   3.8 关键综述
   3.9 关键发现（新颖性判断）
4. 现有数据库与注释资源
   4.1 专用数据库（PhaDED 等）
   4.2 通用注释资源（GTDB-Tk/eggNOG/UniProt/KEGG/iTOL）
   4.3 GTDB 数据资源与获取方式
5. 已有项目与分析方法调研
   5.1 已发表的大规模筛选研究（基因组/宏基因组）
   5.2 方法学范例（HMM/DIAMOND 挖掘流程）
   5.3 系统发育与生态分布分析方法
   5.4 bioRxiv 预印本
6. 关键文献清单
7. 方法学建议与后续分析流程
8. 参考来源汇总（URL）
9. 待核实/待办事项

## 1. 项目背景与目标

### 1.1 背景

- PHB（poly(3-hydroxybutyrate)）是最常见、研究最充分的 PHA 家族成员，
  细菌胞内碳/能量储存聚合物（carbonosome），也是最重要的生物可降解
  塑料之一（生物降解是"生物塑料循环经济"的核心环节）。
- 环境中的 PHB 被微生物分泌的胞外解聚酶（e-PhaZ）水解为单体/寡聚体并
  吸收代谢；胞内颗粒则在碳源耗竭时被胞内解聚酶（i-PhaZ）动员降解。
- GTDB（Genome Taxonomy Database）提供基于基因组系统发育的统一物种
  分类，及其约 47–50 万代表性基因组的蛋白质文件与元数据——为"在全域
  微生物基因组中系统筛查 PHB 降解基因并做生态分布分析"提供了数据基础。

### 1.2 目标

1. 建立 PHB 降解基因（胞内/胞外 PhaZ、寡聚体水解酶、3HB 代谢酶）的
   参考序列目录与 HMM profile（基于文献与 PhaDED 分类体系）。
2. 基于 GTDB 代表性基因组蛋白库进行全库筛选与功能注释（含结构域与
   催化位点双重验证）。
3. 系统发育分析（基因树 + GTDB 物种映射）与生态/分类学分布统计。
4. 交付可复现流程与完整报告。

### 1.3 本报告目的

广泛调研现有文献与已有项目/数据库/工具，明确：已有知识体系、可复用
资源、方法学范例、研究缺口（新颖性），为后续分析流程设计提供依据。

## 2. 文献调研方法与检索策略

| 渠道 | 工具 | 检索内容 | 结果 |
|------|------|----------|------|
| PubMed | NCBI E-utilities | 8 组查询（解聚酶/phaZ/降解/基因组/综述） | 154 篇去重论文，26 篇重点 |
| OpenAlex | OpenAlex API | 12+ 组查询（标题/摘要精准检索 + 高被引排序 + DOI 核实） | 24 篇关键文献 |
| Europe PMC | Europe PMC API | 8 组开放获取检索 + 4 篇综述/论文全文下载 + bioRxiv 预印本调研 | 13 篇高相关 + 7 篇预印本 |
| bioRxiv | Web 搜索 | 预印本检索（重点 Viljakainen & Hug） | 7 篇预印本（详见 5.4） |
| Web | Web 搜索 | 专用数据库（PAZy/PlasticDB/ESTHER）、GTDB 数据资源、已发表项目/标杆路线 | 已核实：无现成 GTDB-PHB 成品流程；GTDB R11-RS232；详见 4.1/4.3/5.2 |

方法学说明：
- PubMed 采用多组同义查询合并去重，批量抓取摘要后按主题筛选；
- OpenAlex 从全文检索（噪声大）转向标题/摘要精准检索（成本低 10 倍），
  对 28 个关键 DOI 做批量核实（28/28 命中）；
- 检索记录与原始数据分别存于 `research/pubmed/`、`research/openalex/`、
  `research/europepmc/`、`research/web/`。

---
<!-- 以下各节由调研结果填充 -->

## 3. PHB 降解生物学与基因目录（文献综述）

> 本节约稿基于 PubMed 检索确认的重点文献（PMID 均来自 2026-06-01 检索结果）。
> PubMed 调研已完成（154 篇去重论文、26 篇重点文献）；Europe PMC 综述全文
> 与其余子代理结果回来后补充细节。

### 3.1 PHB 代谢总览（合成-储存-降解）

- PHB 是 PHA 家族最普遍成员，细菌的碳/能量储存聚合物；合成与降解均与
  营养条件（碳过剩/限氮诱导合成，碳源耗竭诱导动员降解）密切相关。
- 经典综述：
  - "Microbial degradation of polyhydroxyalkanoates"（Annu Rev Microbiol,
    PMID [12213937](https://pubmed.ncbi.nlm.nih.gov/12213937/)）
  - "Biodegradation of polyhydroxyalkanoic acids"（Appl Microbiol Biotechnol,
    PMID [9008883](https://pubmed.ncbi.nlm.nih.gov/9008883/)）
  - "Polyhydroxybutyrate: plastic made and degraded by microorganisms"
    （Rev Environ Contam Toxicol, PMID [9921137](https://pubmed.ncbi.nlm.nih.gov/9921137/)）
- 系统层面："The 'PHAome'"（Trends Biotechnol, PMID
  [26409775](https://pubmed.ncbi.nlm.nih.gov/26409775/)）提出 PHA 代谢"组学"
  视角，强调合成、颗粒形成、降解、调控的基因网络——为基因目录设计提供框架。
- 调控：*Ralstonia eutropha*（*Cupriavidus necator*）PHA 代谢调控研究
  （J Biol Chem, PMID [38969063](https://pubmed.ncbi.nlm.nih.gov/38969063/)）。

### 3.2 胞内 PHA 解聚酶（i-PhaZ）

- 功能：胞内颗粒动员降解；代表系统 *Cupriavidus necator*（*Ralstonia
  eutropha*）PhaZ 家族。
- 关键文献：
  - PhaC1 合酶与 PhaZa1 解聚酶共同作用研究（AEM, PMID
    [29678915](https://pubmed.ncbi.nlm.nih.gov/29678915/)）
  - *Bacillus thuringiensis* phaZ（新型胞内解聚酶，J Bacteriol, PMID
    [16936025](https://pubmed.ncbi.nlm.nih.gov/16936025/)）
  - *Azospirillum brasilense* PhaZ（Arch Microbiol, PMID
    [12898135](https://pubmed.ncbi.nlm.nih.gov/12898135/)）——phoZ 敲除
    丧失降解能力，证实其必需性
  - *Sinorhizobium meliloti* 胞内 PhaZ 鉴定（BMC Microbiol, PMID
    [20346169](https://pubmed.ncbi.nlm.nih.gov/20346169/)）——固氮共生菌
    PHB 循环降解支路
  - *Pseudomonas putida* KT2442 中链长 PHA（mcl-PHA）周转（Environ
    Microbiol, PMID [19788655](https://pubmed.ncbi.nlm.nih.gov/19788655/)；
    J Biol Chem, PMID [17170116](https://pubmed.ncbi.nlm.nih.gov/17170116/)）
  - 注意：*Rhodospirillum rubrum* 的"胞内"PHB 解聚酶实为周质酶（J
    Bacteriol, PMID [15489436](https://pubmed.ncbi.nlm.nih.gov/15489436/)）
    ——挑战胞内/胞外二分法，提示定位注释需谨慎（信号肽预测）。
  - *Cupriavidus necator* 不同营养条件下的颗粒形成与降解（J Basic
    Microbiol, PMID [34342882](https://pubmed.ncbi.nlm.nih.gov/34342882/)）

### 3.3 胞外 PHA 解聚酶（e-PhaZ）与结构域分类

**分类体系（PhaDED, Knoll 2009 — 权威框架）**：587 条序列划为
**8 超家族 / 38 同源家族**（完整方案见 `knowledge/family_classification.md`）：

| 超家族 | 定位/底物 | 催化特征 | 代表 |
|--------|----------|---------|------|
| i-nPHASCL（无 lipase box） | 胞内/天然 SCL | **Cys-His-Asp**（催化 Cys-1 位疏水 Val） | *R. eutropha* PhaZ1/2/3/5、PhaZd |
| i-nPHASCL（有 lipase box） | 胞内/天然 SCL | Ser-His-Asp，x₁=Trp | *B. thuringiensis* |
| 周质 PHA 解聚酶 | 周质/天然 PHB | type 2 催化域，x₁=Ile | *R. rubrum* PhaZ1 |
| i-nPHAMCL | 胞内/天然 MCL | x₁=Val，**lid 结构域** | *P. putida* KT2440 PhaZKT |
| e-dPHASCL（催化域 type 1） | 胞外/变性 SCL | x₁=Leu/Ile；oxyanion hole 在 box N 端 | *A. faecalis*、*R. pickettii* T1、*B. megaterium*、*P. lemoignei* |
| e-dPHASCL（催化域 type 2） | 胞外/变性 SCL | oxyanion hole 在三元组 C 端 | *Acidovorax* TP4、*Comamonas*、*S. exfoliatus*、*P. funiculosum* |
| e-nPHASCL（PhaZ7 型） | 胞外/仅天然 SCL | **AHSMG** 基序（非 GxSxG）；14aa lid | *P. lemoignei* PhaZ7 |
| e-dPHAMCL | 胞外/变性 MCL | x₁=Ile；无 SBD/linker | *P. fluorescens* GK13、*T. thermophilus* HB8 |

- 结构域架构（胞外经典）：信号肽 → 催化域 → 连接域（Fn3/Thr-rich/Cad）→
  SBD（SBD1/SBD2）；dPHAMCL 无 SBD/linker
- **验证四基序**（LtPHBase/PDB 8DAJ）：Ser `IDXXXXYVXGLSXGG`、Asp
  `GXXDYTV`、His `GMXHXXPXXG`、oxyanion hole `HGCXQ`；三联体
  Ser121-His270-Asp197
- **lipase box 疏水 x₁ 是 PHA 解聚酶区别于脂酶/酯酶的关键**
- Pfam 对应：**PF10503**（酯酶型 PHB 解聚酶）+ **PF06850**（SBD C 端）
- ⚠️ 注释陷阱：*R. eutropha* PhaZ6/Z7 属胞外型但为胞内基因；GenBank
  注释不可靠，**必须以序列聚类为准**（来源：PMC2666664, PMC11893044）

- 功能：水解胞外 PHB 聚合物；分布于土壤、海洋、热泉等环境微生物。
- 关键文献：
  - *Comamonas acidovorans* PHB 解聚酶的生化与分子表征（AEM, PMID
    [9406404](https://pubmed.ncbi.nlm.nih.gov/9406404/)）——经典酶学模型
  - *Thermus thermophilus* HB8 胞外 PHB 解聚酶（Appl Microbiol Biotechnol,
    PMID [19214501](https://pubmed.ncbi.nlm.nih.gov/19214501/)）——嗜热酶
  - *Streptomyces ascomycinicus* 新型胞外 PHB 解聚酶与共聚物降解（PLoS ONE,
    PMID [23951224](https://pubmed.ncbi.nlm.nih.gov/23951224/)）
  - *Penicillium funiculosum* PHB 解聚酶晶体结构（J Mol Biol, PMID
    [16405909](https://pubmed.ncbi.nlm.nih.gov/16405909/)）——真菌酶结构
  - 重要概念澄清："To be or not to be a PHB depolymerase: PhaZd1 (PhaZ6)
    and PhaZd2 (PhaZ7)"（AEM, PMID
    [24907326](https://pubmed.ncbi.nlm.nih.gov/24907326/)）——部分 PhaZ
    同源物可能并非真正解聚酶，提示功能验证/结构域筛选的必要性
  - 热嗜土壤菌 PHB 解聚酶降解生物塑料（Protein Sci, PMID
    [36222314](https://pubmed.ncbi.nlm.nih.gov/36222314/)）
  - 海洋：*Alteromonas* 属胞外 PHB 解聚酶鉴定（Mar Biotechnol, PMID
    [40500476](https://pubmed.ncbi.nlm.nih.gov/40500476/)）；海洋细菌 PHA
    降解菌分离（Sci Rep, PMID [40320445](https://pubmed.ncbi.nlm.nih.gov/40320445/)）
  - *Nocardiopsis dassonvillei* PHB 解聚酶过表达与结构建模（Enzyme Microb
    Technol, PMID [41151231](https://pubmed.ncbi.nlm.nih.gov/41151231/)）
  - *Streptomyces exfoliatus* K10 mcl-PHA 解聚酶（脂酶框 Ser-His-Asp, Appl
    Microbiol Biotechnol, PMID [26156240](https://pubmed.ncbi.nlm.nih.gov/26156240/)）
  - 嗜热土壤菌 *Lihuaxuella thermophila* PHB 解聚酶（碱性酶、广谱降解，
    Protein Sci, PMID [36222314](https://pubmed.ncbi.nlm.nih.gov/36222314/)）
- 结构域特征（Europe PMC 综述全文核实）：
  - **e-PhaZ 结构域架构**：N 端信号肽 SP + 催化域 CD + 连接域 LD
    （纤连蛋白 III 型）+ C 端底物结合域 SBD；催化三联体 Ser-Asp-His，
    Ser 位于 **G-X-S-X-G 脂肪酶盒**；SBD 保守基序 **sxxxHxxAGRa**；
    两步机制（SBD 吸附聚合物表面 → CD 水解酯键）；endo 型（产物
    寡聚体 ≤ 五聚体）或 exo 型。
  - **i-PhaZ**：作用于天然颗粒（carbonosome），常需颗粒表面蛋白激活；
    *R. rubrum* 三拷贝系统（Europe PMC 全文核实）：PhaZ1（周质）、
    PhaZ2（胞质、主要胞内解聚酶）、PhaZ3（胞质、体内疑似无活性）；
    ApdA（解聚激活蛋白）、PhaR（颗粒结合转录调控子）。
  - 非特异酯酶（triacylglycerol lipase-like）也能水解 PHA——筛选时
    需注意旁系同源干扰。
- 分类注意：PhaZ 命名跨菌混乱（PhaZ1–7 序号、PcaD→phaZ 等），HMM 筛选
  必须叠加催化位点（lipase box / Ser-Asp-His）验证，不能仅凭基因名。

### 3.7 基因组/宏基因组层面的筛选研究（直接方法学参考）

- **demonstrated vs predicted 分层**：Martínez-Tobón 2018（Appl Microbiol
  Biotechnol, PMID [29951858](https://pubmed.ncbi.nlm.nih.gov/29951858/)）
  ——9 株菌（5 实测 + 4 基因组预测）对比 PHB 降解能力，给出证据分级范式。
- **全基因组→基因定位→功能验证**：Morohoshi 2020（PLoS ONE, PMID
  [32369496](https://pubmed.ncbi.nlm.nih.gov/32369496/)）基于
  *Undibacterium* 完整基因组鉴定新型胞外解聚酶 phaZUD。
- **基因组找同源+跨物种系统发育分布**：Iseki 2025（Mar Biotechnol, PMID
  [40500476](https://pubmed.ncbi.nlm.nih.gov/40500476/)）*Alteromonas*
  解聚酶——与"GTDB 树+基因树"结合思路高度契合。
- **in silico 数据库挖掘**：Leadbeater 2022（Microb Genom, PMID
  [36125959](https://pubmed.ncbi.nlm.nih.gov/36125959/)）鉴定海藻降解型
  生物塑料生产者。
- **环境宏基因组+降解实验**：Omura 2024（Nat Commun, PMID
  [38278791](https://pubmed.ncbi.nlm.nih.gov/38278791/)）深海海底
  757–5552 m 部署实验；Jeon 2023（J Microbiol Biotechnol, PMID
  [37311705](https://pubmed.ncbi.nlm.nih.gov/37311705/)）双层平板+通用引物
  PCR 筛 *Bacillus infantis* phaZ/bdhA；Hachisuka 2025（J Appl Microbiol,
  PMID [40392676](https://pubmed.ncbi.nlm.nih.gov/40392676/)）土壤菌降解
  P(2HB-co-3HB) 共聚物。

### 3.8 关键综述

- **Jendrossek 2002**（Annu Rev Microbiol, PMID
  [12213937](https://pubmed.ncbi.nlm.nih.gov/12213937/)）——奠基综述：
  e-PHA 解聚酶（EC 3.1.1.75/76）与 i-PHA 解聚酶全貌，PhaZ 分类起点
- **Jendrossek 1996**（Appl Microbiol Biotechnol, PMID
  [9008883](https://pubmed.ncbi.nlm.nih.gov/9008883/)）——PHA 生物降解：
  降解菌分离、酶生化、水解机制
- Hankermeyer 1999（Rev Environ Contam Toxicol, PMID
  [9921137](https://pubmed.ncbi.nlm.nih.gov/9921137/)）；Park 2024
  （Biotechnol Adv, PMID [38272380](https://pubmed.ncbi.nlm.nih.gov/38272380/)）
  ——PHA 应用全景；Kliem 2020（Materials, PMID
  [33076314](https://pubmed.ncbi.nlm.nih.gov/33076314/)）——不同环境聚合物
  降解；Bátori 2018（Waste Manag, PMID
  [30455023](https://pubmed.ncbi.nlm.nih.gov/30455023/)）——生物塑料厌氧
  降解；Meng 2014（Curr Opin Biotechnol, PMID
  [24632193](https://pubmed.ncbi.nlm.nih.gov/24632193/)）——聚酯合成-降解
  平衡工程

### 3.9 关键发现（新颖性判断）

- **最接近的前人研究**：Viljakainen & Hug 2021（Environ Microbiol）对
  PHA 降解基因做过系统发育 + 全球分布分析（基于 NCBI/IMG 等基因组集合）。
- 未发现**基于 GTDB 统一分类框架**的 PHB 降解基因系统筛查直接先例：
  GTDB + 酶/基因组挖掘组合检索仅 ~20 条且几乎全是低引用新作。
- 本项目差异化：①GTDB 骨架树 + reps 去冗余 + 元数据生态标签（统一
  分类框架下的全域分布）；②聚焦 PHB 降解基因目录的系统化（胞内/胞外
  PhaZ + 3HB 代谢 + 基因簇共定位），不限于解聚酶单基因。
- 现有工作多为单菌酶学、特定生境分离、或小规模基因组挖掘；降解方向
  文献显著少于合成/生产方向。
- 可借鉴范式：Viljakainen & Hug 2021、PET 酶挖掘（Yoshida 2016
  Science）、宏基因组酶筛选（AEM 2018）。

### 3.4 寡聚体水解酶与 3HB 单体代谢

- 寡聚体水解酶（3HB-oligomer hydrolase）将 (3HB)n 寡聚体水解为单体，
  与解聚酶协同完成 PHB 完全单体化；分布于胞内/胞外/周质。
- 关键文献（PubMed 2026-06-01 补充检索）：
  - 功能整合："Roles of poly(3-hydroxybutyrate) depolymerase and
    3HB-oligomer hydrolase in bacterial PHB metabolism"（Int J Biol
    Macromol?, PMID [15170237](https://pubmed.ncbi.nlm.nih.gov/15170237/)）；
    "oligomer hydrolase 与脱氢酶参与 PHB 动员"（PMID
    [24271169](https://pubmed.ncbi.nlm.nih.gov/24271169/)）；解聚酶与
    寡聚体水解酶关系（PMID [27059479](https://pubmed.ncbi.nlm.nih.gov/27059479/)）
  - 胞内型：*Cupriavidus/Wautersia eutropha* H16（PMID
    [16030206](https://pubmed.ncbi.nlm.nih.gov/16030206/)、
    [16233278](https://pubmed.ncbi.nlm.nih.gov/16233278/)、
    [12070691](https://pubmed.ncbi.nlm.nih.gov/12070691/)）、*Zoogloea
    ramigera* I-16-M（PMID [1476778](https://pubmed.ncbi.nlm.nih.gov/1476778/)、
    [7285912](https://pubmed.ncbi.nlm.nih.gov/7285912/)）
  - 胞外型：*Alcaligenes faecalis*（PMID
    [6626560](https://pubmed.ncbi.nlm.nih.gov/6626560/)）、*Pseudomonas*
    （PMID [8981982](https://pubmed.ncbi.nlm.nih.gov/8981982/)）、
    *Paracoccus denitrificans*（PMID
    [11814660](https://pubmed.ncbi.nlm.nih.gov/11814660/)）
- 胞内动员通路（*Azospirillum* 综述核实，Europe PMC 全文
  `ft_phb_metabolism_azospirillum_review.txt`）：
  PhaZ 解聚酶 → 3-羟基丁酸单体 → NAD(P) 依赖的 3-羟基丁酸脱氢酶
  （BdhA，四聚体）→ 乙酰乙酸 → 乙酰乙酸-CoA 合成酶 → 乙酰乙酰-CoA →
  β-酮硫解酶 → 2 × 乙酰-CoA → TCA/乙醛酸/β-氧化
- 通路酶（β-酮硫解酶、乙酰乙酰-CoA 还原酶、PHB 合酶、BdhA、乙酰乙酸-CoA
  合成酶）在 *A. brasilense* 中组成型表达；NAD(P)H 升高促进 PhaZ 启动降解。
- 分泌相关：*Ralstonia pickettii* T1 解聚酶分泌通路（PMID
  [18340545](https://pubmed.ncbi.nlm.nih.gov/18340545/)）与转录抑制（PMID
  [24146107](https://pubmed.ncbi.nlm.nih.gov/24146107/)）。

### 3.5 基因簇组织

- *R. eutropha* / *A. brasilense*：phbCAB 操纵子（phbC 与 phbAB 可反向）；
  *Azotobacter vinelandii*：phbBAC；*P. putida*：phaC1ZC2D + phaIF 双操纵子。
- 基因可分布在染色体与质粒（*A. baldaniorum* Sp245 的 phbC 在染色体、
  phbCAB 在质粒 4、phbB 拷贝在质粒 1/2）——提示 GTDB 筛选时同源物
  可多拷贝、跨复制子分布。
- 颗粒相关蛋白（GAPs）：合酶、解聚酶、调控子、phasin（phasin_2 结构域
  Pfam **PF09361**）；PhaP 涂层稳定颗粒并控制粒径。
- 大多数 PHA 产生菌编码多个解聚酶同工酶（isoenzymes）。

### 3.6 生态分布

- 深海：微生物分解可生物降解塑料（Nat Commun, PMID
  [38278791](https://pubmed.ncbi.nlm.nih.gov/38278791/)）
- 海冰细菌（AEM, PMID [34160268](https://pubmed.ncbi.nlm.nih.gov/34160268/)）
- 苏打湖 PHB 降解酶（Environ Microbiol Rep, PMID
  [41702408](https://pubmed.ncbi.nlm.nih.gov/41702408/)）
- 分布规律（Europe PMC 综述核实）：
  - 估计约 **10% 微生物**具有 PHB 降解能力；P(3HB) 解聚酶在土壤、淡水、
    海洋（含深海 5500 m）宏基因组中广泛存在
  - 已知携带 e-PhaZ 的类群：*Acidovorax*、*Undibacterium*、
    *Janthinobacterium*、*Massilia*、*Duganella*、*Herbaspirillum*
    （Oxalobacteraceae 内部分布）、*Burkholderia*、*Pseudomonas*
    （mcl-PHA）、*Ralstonia pickettii*、*Cupriavidus*/*R. eutropha*、
    *Azotobacter*、*Azospirillum*、*Rhodococcus*、*Streptomyces* 等
  - scl-PHA 与 mcl-PHA 解聚酶驱动不同降解菌群落 → 挖掘需分亚家族
- 提示：PHB 解聚酶生态分布研究活跃但多集中于特定生境；Viljakainen &
  Hug 2021 已做全球分布分析，基于 GTDB 全域的系统分布仍有差异化空间。

## 4. 现有数据库与注释资源

### 4.1 专用数据库

- **PHA Depolymerase Engineering Database（PhaDED, 2009）** — 最直接相关！
  "The PHA Depolymerase Engineering Database: A systematic analysis tool
  for the diverse family of polyhydroxyalkanoate (PHA) depolymerases"
  （BMC Bioinformatics, 2009, DOI
  [10.1186/1471-2105-10-89](https://doi.org/10.1186/1471-2105-10-89)，
  162 被引；PMID [19296857](https://pubmed.ncbi.nlm.nih.gov/19296857/)）：
  系统收录并分类 PHA 解聚酶家族（催化域分类 + 结构域组成），按底物
  特异性与结构域组织划分多个家族，支持 HMM 家族分配。
  - **托管于 ESTHER**（α/β-水解酶折叠蛋白数据库，INRAE）：
    [Esterase_phb_PHAZ 家族页](https://bioweb.supagro.inrae.fr/ESTHER/family/Esterase_phb_PHAZ)、
    [PHAZ7_phb_depolymerase 家族页](https://bioweb.supagro.inrae.fr/ESTHER/family/PHAZ7_phb_depolymerase)、
    [论文记录页](https://bioweb.supagro.inrae.fr/ESTHER/paper/Knoll_2009_BMC.Bioinformatics_10_89)；
    [Database Commons 记录](https://ngdc.cncb.ac.cn/databasecommons/database/id/3446)
  - 是建立参考序列目录与分家族 HMM 的首要种子来源（序列、家族分类、
    HMM 均可直接复用）。
- **PAZy（Plastics-Active Enzymes Database）**（Buchholz et al. 2022,
  Proteins, DOI [10.1002/prot.26325](https://doi.org/10.1002/prot.26325)）：
  按底物（PET、PHB、PCL 等）分类的塑料活性酶库，含 PHB 解聚酶种子序列、
  结构与活性数据；提供 [API](https://www.pazy.eu/api-docs) 与
  [DaRUS 数据镜像](https://darus.uni-stuttgart.de/dataverse/ibc_tbc_PAZy)。
- **PlasticDB**（Gambarini et al. 2022, Database, DOI
  [10.1093/database/baac008](https://doi.org/10.1093/database/baac008)，PMID
  [35266524](https://pubmed.ncbi.nlm.nih.gov/35266524/)）：塑料生物降解
  微生物与蛋白质库，**含 HMM 模型可直接用于筛选**。
- **PlasticEnz**（PLoS Comput Biol, 2026, DOI
  [10.1371/journal.pcbi.1013892](https://doi.org/10.1371/journal.pcbi.1013892)，
  PMID [41587207](https://pubmed.ncbi.nlm.nih.gov/41587207/)）：同源 +
  机器学习的塑料降解酶综合数据库与筛选工具（含 PHA 解聚酶），可作
  本项目注释交叉验证工具。
- **序列/结构锚点**：
  - BRENDA [EC 3.1.1.75](https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.75)
    poly(3-hydroxybutyrate) depolymerase 条目
  - UniProt [EC 3.1.1.75 reviewed 条目（3 条）](https://www.uniprot.org/uniprotkb?query=(ec:3.1.1.75)%20AND%20(reviewed:true))
  - PDB 结构：[4BTV（PhaZ7-3HB 复合物，*Paucimonas lemoignei*）](https://www.rcsb.org/structure/4BTV)、
    [8YNV（*Bacillus thuringiensis* PhaZ）](https://www.rcsb.org/structure/8YNV)
  - NCBI CDD [poly(3-hydroxyalkanoate) depolymerase 保守域模型](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=XP_039073312.1)
- 基于已表征/预测解聚酶活性的菌株研究（Appl Microbiol Biotechnol 2018,
  PMID [29951858](https://pubmed.ncbi.nlm.nih.gov/29951858/)）提供
  "demonstrated vs predicted" 分层证据思路——可用于给筛选结果分级。

### 4.2 通用注释资源（OpenAlex 核实）

- **GTDB 原始论文系列**：
  - "A standardized bacterial taxonomy based on genome phylogeny"
    （Nat Biotechnol, 2018, DOI
    [10.1038/nbt.4229](https://doi.org/10.1038/nbt.4229)，3964 被引）
  - "GTDB: an ongoing census of bacterial and archaeal diversity"
    （NAR, 2021, DOI
    [10.1093/nar/gkab776](https://doi.org/10.1093/nar/gkab776)，2571 被引）
  - **GTDB-Tk**（Bioinformatics, 2019, DOI
    [10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)，
    5648 被引）与 **GTDB-Tk v2**（Bioinformatics, 2022, DOI
    [10.1093/bioinformatics/btac672](https://doi.org/10.1093/bioinformatics/btac672)，
    2018 被引）
- **eggNOG 5.0**（NAR, 2018, DOI
  [10.1093/nar/gky1085](https://doi.org/10.1093/nar/gky1085)）：直系同源
  功能注释（可查 PHA 相关 OG）
- **UniProt**（NAR, 2016, DOI
  [10.1093/nar/gkw1099](https://doi.org/10.1093/nar/gkw1099)）：种子序列
  与酶学注释（EC 3.1.1.75/76 等）
- **KEGG taxonomy-based analysis**（NAR, 2022, DOI
  [10.1093/nar/gkac963](https://doi.org/10.1093/nar/gkac963)）：通路级
  比较分析（PHA 合成/降解通路）
- **KEGG PHA 通路**（Web 调研核实）：
  - [map00640 Polyhydroxyalkanoate metabolism](https://www.kegg.jp/entry/pathway+map00640)
    ——PHA 合成与降解通路总图
  - 模块 [M00843](https://www.kegg.jp/entry/M00843)；KO **K03821 = phaC/phaC
    合成酶亚基**；⚠️ 注意 **M00012 是乙醛酸循环，不是 PHA 模块，勿误用**
  - 备注：大量 GTDB 基因组自带 KEGG 注释，可交叉查询
- **CAZy**：未发现专设 PHA 解聚酶家族；相关 α/β-水解酶超家族/cutinase
  类可参考 ESTHER + PAZy + NCBI CDD 组织（Web 调研结论）
- **iTOL v6**（NAR, 2024, DOI
  [10.1093/nar/gkae268](https://doi.org/10.1093/nar/gkae268)）：系统发育
  树可视化
- 背景综述：Bioplastics for a circular economy（Nat Rev Mater, 2022,
  DOI [10.1038/s41578-021-00407-8](https://doi.org/10.1038/s41578-021-00407-8)）；
  Microbial and Enzymatic Degradation of Synthetic Plastics（Front
  Microbiol, 2020, DOI
  [10.3389/fmicb.2020.580709](https://doi.org/10.3389/fmicb.2020.580709)）
- 方法学范式参考：PET 降解酶挖掘——Yoshida 2016（Science, DOI
  [10.1126/science.aad6359](https://doi.org/10.1126/science.aad6359)，3342
  被引）；宏基因组酶筛选（AEM 2018, DOI
  [10.1128/aem.02773-17](https://doi.org/10.1128/aem.02773-17)，445 被引）；
  深海海绵放线菌 cutinase 样聚酯酶 BgP 基因组挖掘（Front Microbiol 2022,
  DOI [10.3389/fmicb.2022.888343](https://doi.org/10.3389/fmicb.2022.888343)）

### 4.3 GTDB 数据资源与获取方式（Web 调研核实）

- **最新版本**：R11-RS232（GTDB 论坛公告
  [Announcing GTDB R11-RS232](https://forum.gtdb.ecogenomic.org/t/announcing-gtdb-r11-rs232/826)）；
  R10 时含 **715,230 细菌 + 17,245 古菌** 基因组（约 73.2 万，
  [生物通报道](https://www.ebiotrade.com/newsf/2025-10/20251023083213816.htm)）。
  此前版本：R09-RS220、R226、R232 等；统计页
  https://gtdb.ecogenomic.org/stats/r232 【下载前请核实最新 release 号与规模】
- **官方站点**：https://gtdb.ecogenomic.org/（数据下载：https://data.gtdb.ecogenomic.org/）
- **关键文件**（release 对应目录下）：
  - `gtdb_proteins_aa_reps.tar.gz` — 代表性基因组蛋白文件（每基因组 .faa，
    本项目主数据源）
  - `gtdb_genome_reps.tar.gz` — 代表性基因组（FNA，用于基因簇上下文）
  - `gtdb_metadata.tsv` — 基因组元数据（GTDB 分类、质量、来源信息）
  - `gtdb_taxonomy_metadata.tsv` — 分类元数据
  - 骨架树文件（bac120/ar53 树）【待核实具体文件名】
- **下载路径**（Web 调研核实）：
  - 总清单：https://data.gtdb.ecogenomic.org/releases/latest/FILE_DESCRIPTIONS.txt
  - 示例目录：release232/232.0/genomic_files_reps/
    （https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/）
- **蛋白文件使用实战**：
  - [GTDB Forum：faa header→GTDB 分类学映射](https://forum.gtdb.ecogenomic.org/t/taxonomy-lookup-for-fasta-headers-from-gtdb-proteins-aa-reps-tar-gz/599/3)
  - [James Lingford：GTDB faa → 带分类学的 DIAMOND 数据库](https://www.jameslingford.com/blog/gtdb-to-diamond-taxonomy-database/)
    （可直接套用，顺带获得"酶→物种/谱系"注释能力）
  - [GTDB Forum：获取参考基因组序列](https://forum.gtdb.ecogenomic.org/t/how-to-obtain-reference-genome-sequences-in-gtdb-database/534/3)
- **规模注意**：R10 起 ~73 万基因组、蛋白库达数百 GB 级 → 需评估下载
  带宽/磁盘；本机磁盘空间受沙箱限制未能实测，建议先下载 metadata 与
  reps 蛋白文件，按需扩展。
- **GTDB-Tk**（v2，[GitHub](https://github.com/Ecogenomics/GTDBTk)）：
  本地分类与去冗余；其数据库包含 HMM 与参考基因组（较大，注意版本对应）。
  若直接用 GTDB 现有代表基因组，可跳过 GTDB-Tk，直接利用其分类元数据。

## 5. 已有项目与分析方法调研

### 5.1 已发表的大规模筛选研究（基因组/宏基因组）

**⭐ 最直接相关的前人研究：Viljakainen & Hug 2021**
"The phylogenetic and global distribution of bacterial polyhydroxyalkanoate
bioplastic-degrading genes"（Environmental Microbiology, 2021, DOI
[10.1111/1462-2920.15409](https://doi.org/10.1111/1462-2920.15409)，PMID
[33496062](https://pubmed.ncbi.nlm.nih.gov/33496062/)；预印本 bioRxiv
[10.1101/2020.05.08.085522](https://www.biorxiv.org/content/10.1101/2020.05.08.085522v3.full)）：
- **方法**：筛选 3078 个宏基因组（1914 Gb），鉴定 13,869 个推定 PHA
  解聚酶（分布于 1295 个宏基因组）；另筛 5290 个 MAG 描述系统发育分布。
- **结果**：解聚酶分布不均——废水系统频率最高、海洋与热泉最低；
  系统发育广度远超培养代表；关键类群 Proteobacteria 与 Bacteroidota，
  另见于 Bdellovibrionota、Methylomirabilota、Actinobacteriota、
  Firmicutes、Spirochaetota、Desulfobacterota、Myxococcota、
  Planctomycetota。
- **局限（= 本项目差异化空间）**：①基于宏基因组与 MAG（非 GTDB 统一
  分类框架）；②聚焦 PHA 解聚酶单基因家族（e-PhaZ），未系统覆盖胞内
  PhaZ、寡聚体水解酶与 3HB 代谢通路；③无基因簇共定位与 GTDB 骨架树
  映射。本项目以 GTDB reps 全库 + 完整降解基因目录 + 生态元数据
  系统化为目标。

其余可借鉴范式：
1. **全基因组→基因定位→功能验证**（Morohoshi 2020, PMID
   [32369496](https://pubmed.ncbi.nlm.nih.gov/32369496/)）
2. **同源扫描+跨物种系统发育分布**（Iseki 2025, PMID
   [40500476](https://pubmed.ncbi.nlm.nih.gov/40500476/)）
3. **demonstrated vs predicted 证据分级**（Martínez-Tobón 2018, PMID
   [29951858](https://pubmed.ncbi.nlm.nih.gov/29951858/)）
4. **in silico 数据库挖掘 pipeline**（Leadbeater 2022, PMID
   [36125959](https://pubmed.ncbi.nlm.nih.gov/36125959/)）
5. **环境宏基因组+降解实验验证**（Omura 2024, PMID
   [38278791](https://pubmed.ncbi.nlm.nih.gov/38278791/)）

### 5.2 方法学范例（HMM/DIAMOND 挖掘流程）

- **⭐ Zrimec et al. 2021（mBio, PMID
  [34700384](https://pubmed.ncbi.nlm.nih.gov/34700384/)，DOI
  [10.1128/mbio.02155-21](https://doi.org/10.1128/mbio.02155-21)）**：
  "Plastic-Degrading Potential across the Global Microbiome"——明确用
  **HMMER v3.3 hmmsearch** 在数千万基因级宏基因组数据中筛选塑料降解酶
  同源序列，并与污染水平关联；是本项目全库筛选的直接方法蓝本。
- **Microbial Genomics 2025**（DOI
  [10.1099/mgen.0.001814](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001814)，
  [Zenodo 配套数据](https://zenodo.org/records/15480170)）：在
  "谱系 × 全球生态系统"尺度预测塑料降解潜力，与本项目思路高度同源。
- **PhaDED 分类体系**（PhaD 2009, DOI
  [10.1186/1471-2105-10-89](https://doi.org/10.1186/1471-2105-10-89)）：
  PHA 解聚酶已有系统分类（催化域分类 + 结构域组成），应作为 phaZ
  参考序列分家族建模的骨架。
- **PET 酶基因组/宏基因组挖掘流程**（范式）：Yoshida 2016（Science）、
  AEM 2018 海洋/陆地宏基因组 PET 酶筛选（DOI
  [10.1128/aem.02773-17](https://doi.org/10.1128/aem.02773-17)）——
  BLAST/HMM → 聚类 → 异源表达验证的标准路径。
- **可用流程/工具**（Web 调研核实）：[plastedma](https://gitlink.org.cn/NSCCN/plastedma)
  （宏基因组塑料降解酶检测工作流，bioconda 可装）、
  [PDETool](https://github.com/ozefreitas/PDETool)（塑料降解酶工具）；
  已发表挖掘实例：*S. microflavus* DG19（2025）、*Burkholderia
  vietnamiensis* mcl-PHA 降解基因（Tn-seq + CRISPR-Cas）、*Priestia*
  USM5 双 PHA 代谢、*Photobacterium ganghwense* 全基因组。
- **基因预测 vs 实验验证分层**：Martínez-Tobón 2018（PMID
  [29951858](https://pubmed.ncbi.nlm.nih.gov/29951858/)）。
- **筛选共识范式**（Web 调研）：HMM（HMMER）+ 种子库
  （PAZy/ESTHER/PlasticDB）+ 催化残基验证；速度优先时
  MMseqs2/DIAMOND 预筛 + hmmsearch 精筛。

### 5.3 系统发育与生态分布分析方法

- *Alteromonas* 解聚酶跨物种系统发育分布研究（PMID
  [40500476](https://pubmed.ncbi.nlm.nih.gov/40500476/)）表明：解聚酶基因
  树 + 物种分类映射是标准做法，与 GTDB 骨架树结合可行
- 可视化：iTOL v6（DOI
  [10.1093/nar/gkae268](https://doi.org/10.1093/nar/gkae268)）、ETE3
- 基因家族分析管线：[PEGP](https://github.com/stovc/pegp/)（同源搜索 +
  系统发育 + 域架构 + 共线性 + phyletic pattern）；anvio（比较/泛基因组
  与基因上下文）；BiG-SCAPE 仅适合 BGC 聚类，不作为本任务主工具
  （Web 调研核实）

### 5.4 bioRxiv 预印本（Europe PMC 子代理调研）

| 预印本 | 链接 | 时间 | 相关点 |
|--------|------|------|--------|
| ① 细菌 PHA 生物塑料降解基因的系统发育与全球分布（Viljakainen & Hug） | [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.05.08.085522v3.full) | 2020 | ⭐ 与本项目最相关；已发表于 Environ Microbiol 2021 |
| ② Burkholderia 胞外生物塑料降解基因（Tn-seq + CRISPR-Cas） | [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.20.655153v1.full) | 2025 | 对应 AEM 2026 正式版 |
| ③ PlasticEnz 数据库与筛选工具 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.10.28.685028v1.full) | 2025 | 已发表于 PLoS Comput Biol |
| ④ 重组胞外 PHB 解聚酶表达纯化流程对比 | [bioRxiv](https://www.biorxiv.org/content/10.1101/700252v1.article-info) | 2019 | 异源表达参考 |
| ⑤ *Microbacterium paraoxydans* PHB 解聚酶动力学 | [bioRxiv](https://www.biorxiv.org/content/10.1101/540609v1.article-metrics) | 2019 | 酶学表征参考 |
| ⑥ *Caldimonas/Schlegelella* PHA 菌基因组比较 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.27.559687v1.supplementary-material) | 2023 | 嗜热 PHA 菌注释参考 |
| ⑦ *Caldimonas thermodepolymerans* 基因组编辑工具集 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.22.614348v1.article-metrics) | 2024 | 下游功能验证参考 |

## 6. 关键文献清单

> 合并 PubMed + OpenAlex 检索结果（去除重复项），按类别列出。完整清单见
> `research/pubmed/summary.md`（26 篇重点，含摘要要点）与
> `research/openalex/PHB_PHA_degradation_literature_report.md`（24 篇，含引用数）。

### 奠基综述（PHA 降解分类框架）
- Jendrossek & Handrick 2002, Annu Rev Microbiol, "Microbial degradation of
  polyhydroxyalkanoates" — [PubMed](https://pubmed.ncbi.nlm.nih.gov/12213937/) ·
  [DOI 10.1146/annurev.micro.56.012302.160838](https://doi.org/10.1146/annurev.micro.56.012302.160838)
- Jendrossek 1996, Appl Microbiol Biotechnol — [PubMed](https://pubmed.ncbi.nlm.nih.gov/9008883/)
- Hankermeyer 1999, Rev Environ Contam Toxicol — [PubMed](https://pubmed.ncbi.nlm.nih.gov/9921137/)
- "Review Degradation of microbial polyesters" 2004, Biotechnol Lett —
  [DOI 10.1023/b:bile.0000036599.15302.e5](https://doi.org/10.1023/b:bile.0000036599.15302.e5)

### 酶学与机制
- *Penicillium funiculosum* PHB 解聚酶晶体结构（1.71 Å，α/β 水解酶折叠 +
  Ser-Asp-His 三联体）— [PubMed 16405909](https://pubmed.ncbi.nlm.nih.gov/16405909/)
  · [DOI 10.1016/j.jmb.2005.12.028](https://doi.org/10.1016/j.jmb.2005.12.028)
- *Comamonas acidovorans* 胞外解聚酶 — [PubMed 9406404](https://pubmed.ncbi.nlm.nih.gov/9406404/)
  · [DOI 10.1128/aem.63.12.4844-4852.1997](https://doi.org/10.1128/aem.63.12.4844-4852.1997)
- PhaZd1/PhaZd2 功能辨析（非真解聚酶）— [PubMed 24907326](https://pubmed.ncbi.nlm.nih.gov/24907326/)
- *R. eutropha* PhaC1/PhaZa1 磷酸化调控 — [PubMed 29678915](https://pubmed.ncbi.nlm.nih.gov/29678915/)
- *P. putida* KT2442 mcl-PHA 解聚酶 — [PubMed 17170116](https://pubmed.ncbi.nlm.nih.gov/17170116/) ·
  [19788655](https://pubmed.ncbi.nlm.nih.gov/19788655/) ·
  "A holistic view of PHA metabolism" [DOI 10.1111/1462-2920.12760](https://doi.org/10.1111/1462-2920.12760)
- *Rhodospirillum rubrum* 周质解聚酶 — [PubMed 15489436](https://pubmed.ncbi.nlm.nih.gov/15489436/) ·
  [DOI 10.1007/s00253-011-3096-7](https://doi.org/10.1007/s00253-011-3096-7)
- *Bacillus thuringiensis* 胞内 phaZ — [PubMed 16936025](https://pubmed.ncbi.nlm.nih.gov/16936025/)
- *Azospirillum brasilense* PhaZ — [PubMed 12898135](https://pubmed.ncbi.nlm.nih.gov/12898135/)
- *Sinorhizobium meliloti* 胞内 PhaZ — [PubMed 20346169](https://pubmed.ncbi.nlm.nih.gov/20346169/)
- 嗜热/耐热酶：*Thermus thermophilus* — [PubMed 19214501](https://pubmed.ncbi.nlm.nih.gov/19214501/)；
  *Lihuaxuella thermophila* — [PubMed 36222314](https://pubmed.ncbi.nlm.nih.gov/36222314/)
- *Streptomyces* 系：*S. ascomycinicus* — [PubMed 23951224](https://pubmed.ncbi.nlm.nih.gov/23951224/)；
  *S. exfoliatus* K10 — [PubMed 26156240](https://pubmed.ncbi.nlm.nih.gov/26156240/)
- *Nocardiopsis dassonvillei* 结构建模 — [PubMed 41151231](https://pubmed.ncbi.nlm.nih.gov/41151231/)
- 颗粒蛋白与调控：phasin 综述
  [DOI 10.1021/bm049401n](https://doi.org/10.1021/bm049401n)；*R. eutropha*
  转录调控 — [PubMed 38969063](https://pubmed.ncbi.nlm.nih.gov/38969063/)；
  "PHAome" 概念 — [PubMed 26409775](https://pubmed.ncbi.nlm.nih.gov/26409775/)

### 基因组/宏基因组筛选与挖掘
- demonstrated vs predicted 解聚酶 — [PubMed 29951858](https://pubmed.ncbi.nlm.nih.gov/29951858/)
  · [DOI 10.1007/s00253-018-9153-8](https://doi.org/10.1007/s00253-018-9153-8)
- *Undibacterium* 全基因组→phaZUD — [PubMed 32369496](https://pubmed.ncbi.nlm.nih.gov/32369496/)
- *Alteromonas* 解聚酶系统发育分布 — [PubMed 40500476](https://pubmed.ncbi.nlm.nih.gov/40500476/)
- in silico 数据库挖掘（海藻降解菌）— [PubMed 36125959](https://pubmed.ncbi.nlm.nih.gov/36125959/)
- 深海降解实验 — [PubMed 38278791](https://pubmed.ncbi.nlm.nih.gov/38278791/)
- 海洋宏基因组 bioplastic 降解 — [DOI 10.3389/fmicb.2019.01252](https://doi.org/10.3389/fmicb.2019.01252)
- PET 酶挖掘范式：Yoshida 2016 Science
  [DOI 10.1126/science.aad6359](https://doi.org/10.1126/science.aad6359)；宏基因组
  PET 酶筛选 [DOI 10.1128/aem.02773-17](https://doi.org/10.1128/aem.02773-17)；
  BgP cutinase 样聚酯酶 [DOI 10.3389/fmicb.2022.888343](https://doi.org/10.3389/fmicb.2022.888343)
- *Ralstonia eutropha* H16 基因组 — [DOI 10.1038/nbt1244](https://doi.org/10.1038/nbt1244)

### 数据库与工具
- **PHA Depolymerase Engineering Database (PhaDED)** — [DOI 10.1186/1471-2105-10-89](https://doi.org/10.1186/1471-2105-10-89)
- GTDB 系列：standardized taxonomy [DOI 10.1038/nbt.4229](https://doi.org/10.1038/nbt.4229)；
  census [DOI 10.1093/nar/gkab776](https://doi.org/10.1093/nar/gkab776)；
  GTDB-Tk [DOI 10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)；
  GTDB-Tk v2 [DOI 10.1093/bioinformatics/btac672](https://doi.org/10.1093/bioinformatics/btac672)
- eggNOG 5.0 [DOI 10.1093/nar/gky1085](https://doi.org/10.1093/nar/gky1085)；UniProt
  [DOI 10.1093/nar/gkw1099](https://doi.org/10.1093/nar/gkw1099)；KEGG taxonomy
  [DOI 10.1093/nar/gkac963](https://doi.org/10.1093/nar/gkac963)；iTOL v6
  [DOI 10.1093/nar/gkae268](https://doi.org/10.1093/nar/gkae268)

### 生态与环境降解背景
- 深海、海冰、苏打湖、海洋菌（见 3.6）；PHA 生物塑料环境降解综述：
  Bioplastics for a circular economy [DOI 10.1038/s41578-021-00407-8](https://doi.org/10.1038/s41578-021-00407-8)；
  PHA 环境降解综述（Europe PMC 全文已下载：`research/europepmc/ft_biodegradability_pha_review.txt`）

## 7. 方法学建议与后续分析流程

> 基于文献调研得出的核心设计决策；详细流程草案见 `docs/analysis_plan_draft.md`。

### 7.1 基因目录设计（基于文献证据）

| 基因类 | 代表/证据 | 筛选要点 |
|--------|-----------|----------|
| 胞内 PhaZ（i-PhaZ） | *C. necator* PhaZ1-3、*B. thuringiensis*、*S. meliloti*、*A. brasilense* | 注意周质型（*R. rubrum*）；与胞外酶序列差异大，单独建模 |
| 胞外 PhaZ（e-PhaZ） | *Comamonas acidovorans*、*P. lemoignei*、*Streptomyces*、*Thermus*、*Alteromonas* | 催化域（α/β 水解酶、G-X-S-X-G）+ SBD + 信号肽；按 PhaDED 分类分家族建模 |
| 寡聚体水解酶 | *C. necator* H16、*A. faecalis*、*Paracoccus denitrificans*（PMID 16030206/6626560/11814660） | 与解聚酶协同单体化；分胞内/胞外型 |
| 3HB 代谢 | BdhA、乙酰乙酸-CoA 合成酶/转移酶、β-酮硫解酶 | 与合成共用 β-酮硫解酶，需注意功能区分（合成 vs 降解方向） |
| 颗粒/调控（可选） | PhaP（phasin_2, PF09361）、PhaR | 用于基因簇/代谢岛背景 |

### 7.2 筛选策略（文献支撑的关键决策）

1. **分家族建模**：胞内/胞外解聚酶序列差异大，且 PhaZ 命名混乱
   （PhaZ1–7、PcaD→phaZ）→ 按 PhaDED 分类 + 已表征种子序列分别构建 HMM。
2. **双重验证**：HMM 命中后检查催化三联体（Ser-Asp-His / Ser-His-Asp）与
   lipase box（G-X-S-X-G）保守性；PhaZd1/PhaZd2 案例（PMID 24907326）
   提示同源物可能非真解聚酶 → 需功能特征过滤。
3. **定位注释**：信号肽预测（SignalP）区分胞外/周质/胞内，参考 *R. rubrum*
   案例（PMID 15489436）。
4. **证据分级**：采用 demonstrated vs predicted 分层（PMID 29951858），
   对高置信命中可抽查验证。
5. **先蛋白库后基因组**：GTDB 蛋白质文件（reps）先做 DIAMOND/HMMER 快速
   筛选，对命中基因组再下载做基因簇上下文分析（方案 A→B 两阶段）。

### 7.3 后续流程（待用户确认环境后启动）

1. 数据获取：GTDB release 选择、下载 reps 蛋白库与 metadata【待 Web 核实 URL】
2. 种子序列收集：UniProt/Swiss-Prot + PhaDED + 上述文献酶 → 分家族 HMM
3. 全库筛选：DIAMOND（粗筛）→ HMMER hmmsearch（精筛）→ 催化位点验证
4. 系统发育：MAFFT → trimAl → IQ-TREE2 → ETE3/iTOL 可视化，映射 GTDB 分类
5. 生态分布：按 GTDB 门/纲 + 基因组来源元数据统计检出率
6. 交付：脚本 + 数据表 + 图表 + 报告

### 7.4 计算环境（已实测，见 analysis_plan_draft.md 第 9 节）

本机缺 HMMER/MAFFT/DIAMOND/Prodigal/IQ-TREE 且无 WSL；建议启用 WSL2 +
Miniconda，或使用远程 Linux/HPC 执行全库规模运算。

## 8. 参考来源汇总（URL）

### 数据库与工具
- PhaDED/ESTHER: https://bioweb.supagro.inrae.fr/ESTHER/family/Esterase_phb_PHAZ
  （论文 DOI 10.1186/1471-2105-10-89）
- PAZy: https://www.pazy.eu/ （论文 DOI 10.1002/prot.26325）
- PlasticDB: https://plasticdb.org/ （论文 DOI 10.1093/database/baac008）
- PlasticEnz: DOI 10.1371/journal.pcbi.1013892
- GTDB: https://gtdb.ecogenomic.org/ ；数据 https://data.gtdb.ecogenomic.org/releases/latest/
- GTDB-Tk: https://github.com/Ecogenomics/GTDBTk
- UniProt: https://www.uniprot.org/ ；eggNOG: http://eggnogdb.embl.de/ ；
  KEGG: https://www.kegg.jp/（map00640）；iTOL: https://itol.embl.de/
- BRENDA EC 3.1.1.75: https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.75
- PDB: 4BTV / 8YNV；NCBI CDD poly(3-hydroxyalkanoate) depolymerase
- 流程工具：plastedma、PDETool、PEGP（见 5.2/5.3）

### 文献检索源
- PubMed: https://pubmed.ncbi.nlm.nih.gov/ （检索记录见 research/pubmed/）
- OpenAlex: https://openalex.org/ （检索记录见 research/openalex/）
- Europe PMC: https://europepmc.org/ （检索记录见 research/europepmc/）

### 项目数据文件（本工作区）
- 报告正文：docs/literature_survey_report.md
- PubMed 汇总：research/pubmed/summary.md、methods_log.md
- OpenAlex 汇总：research/openalex/PHB_PHA_degradation_literature_report.md
- Europe PMC 全文：research/europepmc/ft_*.txt
- 知识框架：knowledge/PHB_degradation_framework.md
- 分析流程草案：docs/analysis_plan_draft.md

## 9. 待核实/待办事项

### 已完成（Europe PMC/bioRxiv 调研）
- [x] bioRxiv 相关预印本（7 篇，见 5.4；重点 Viljakainen & Hug 已发表于
      Environ Microbiol 2021, DOI 10.1111/1462-2920.15409）
- [x] e-PhaZ 结构域架构（SP-CD-LD-SBD）、GXSXG 脂肪酶盒、SBD 保守基序、
      EC 3.1.1.75（Europe PMC 全文核实）
- [x] PlasticEnz 数据库工具（PLoS Comput Biol 2026）
- [x] *R. rubrum* PhaZ1/2/3 与 ApdA 分类细节（Microbial Cell Factories 2026 全文）
- [x] 分布规律：约 10% 微生物具 PHB 降解能力；e-PhaZ 已知类群清单
- [x] Viljakainen & Hug 2021 正式发表信息核实（Environ Microbiol）
- [x] 寡聚体水解酶文献（14 篇，见 3.4）

### 已解决（Web 调研 2026-06-01）
- [x] GTDB 最新 release：R11-RS232（R10 = 715,230 细菌 + 17,245 古菌，
      Parks et al. 2025；统计页 gtdb.ecogenomic.org/stats/r232）
- [x] GTDB reps 蛋白库下载：data.gtdb.ecogenomic.org/releases/latest/
      （gtdb_proteins_aa_reps.tar.gz）；faa→分类学映射与 DIAMOND 化实战
- [x] PhaDED 可访问性：托管于 ESTHER（INRAE），家族页可访问、HMM 可复用
- [x] PAZy（Buchholz 2022）/ PlasticDB（Gambarini 2022）存在性与内容、许可
- [x] KEGG：map00640（PHA 代谢通路）、M00843、K03821=phaC；M00012 为
      乙醛酸循环（警告勿误用）
- [x] 大库筛选工具链：HMMER v3.3 hmmsearch（Zrimec 2021）+ DIAMOND/MMseqs2
      预筛；plastedma、PDETool 现成流程
- [x] Pfam：未检索到明确 PHB 解聚酶家族号（需官网复核）；phasin_2 =
      PF09361（已核实）；建议以 ESTHER 家族划分为准

### 遗留低优先复核项（后续在官网确认，不影响流程启动）
- [ ] Pfam PHB 解聚酶家族编号（InterPro 网站复核）
- [ ] KEGG M00843 模块标题、K03513/K17745 具体功能
- [ ] GTDB 骨架树文件在 release 目录中的确切文件名
- [ ] 本机磁盘空间实测（沙箱限制未能测量）

### 待后续阶段执行
- [ ] 种子序列收集（UniProt + PhaDED + 文献酶）→ 分家族 HMM
- [ ] 计算环境搭建（WSL2/miniconda 或远程）
- [ ] GTDB 数据下载与校验
- [ ] 全库筛选 → 双重验证 → 系统发育 → 生态分布
- [ ] 完整分析流程与结果报告交付
