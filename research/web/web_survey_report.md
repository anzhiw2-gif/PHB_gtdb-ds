# 「基于 GTDB 数据库的 PHB 降解基因系统生信分析」Web 调研报告

> 调研方式：系统化 `web_search` 检索（A/B/C/D 四组主题，共 20 余次查询）
> 说明：本报告所有数据库/工具/文献名称、URL 与内容描述均来自 web_search 实际返回结果；未能核实的内容已明确标注"需在官网复核"，未凭记忆补充。
> 生成日期：2026 年（检索到的最新文献为 2025-2026 年发表/预印本）

---

## 目录

1. [专用数据库（A 部分）](#一专用数据库a-部分)
2. [通用功能注释资源（B 部分）](#二通用功能注释资源b-部分)
3. [基因组挖掘方法与已发表项目（C 部分）](#三基因组挖掘方法与已发表项目c-部分)
4. [分析与可视化工具（D 部分）](#四分析与可视化工具d-部分)
5. [方法学建议](#五方法学建议)
6. [参考来源 URL 汇总](#六参考来源-url-汇总)
7. [给主代理的关键结论](#七给主代理的关键结论)

---

## 一、专用数据库（A 部分）

### A1. PHA Depolymerase Engineering Database（PHA 解聚酶工程数据库）★ 最相关

- **论文**：[Knoll et al. 2009, BMC Bioinformatics 10:89](https://pubmed.ncbi.nlm.nih.gov/19296857/)（PMID 19296857），doi: 10.1186/1471-2105-10-89（开放获取）
  - [Springer 全文](https://link.springer.com/article/10.1186/1471-2105-10-89?optIn=true)
  - [Europe PMC (PMC2666664)](http://staging.europepmc.org/backend/articlerender.fcgi?accid=PMC2666664)
- **托管位置**：ESTHER（α/β-水解酶折叠蛋白数据库）系统内
  - [ESTHER 论文记录页](https://bioweb.supagro.inrae.fr/ESTHER/paper/Knoll_2009_BMC.Bioinformatics_10_89)
  - 家族页：[Esterase_phb_PHAZ](https://bioweb.supagro.inrae.fr/ESTHER/family/Esterase_phb_PHAZ)、[PHAZ7_phb_depolymerase](https://bioweb.supagro.inrae.fr/ESTHER/family/PHAZ7_phb_depolymerase)
  - [Database Commons 记录（ID 3446）](https://ngdc.cncb.ac.cn/databasecommons/database/id/3446)
- **内容**：系统整理 PHA 解聚酶多样性，按底物特异性与结构域组织划分为多个家族（含胞外各型与胞内 PhaZ），支持基于 HMM 的家族分配——本项目筛选所需的核心种子资源。
- **复用性**：序列、家族分类、HMM 均可直接用于构建筛选流程；学术免费。

### A2. PAZy（Plastics-Active Enzymes Database，塑料活性酶数据库）

- **入口**：
  - [cbl.uni-stuttgart.de 主页](https://www.cbl.uni-stuttgart.de/doku.php?id=start)
  - [pazy.eu（按底物浏览示例）](https://www.pazy.eu/plastics/pet?substrate=1&limit=50&offset=100)
  - [API 文档](https://www.pazy.eu/api-docs)
  - [数据镜像（DaRUS dataverse: ibc_tbc_PAZy）](https://darus.uni-stuttgart.de/dataverse/ibc_tbc_PAZy)
  - [另一 dataverse 入口（ibtb_pazy）](https://darus.uni-stuttgart.de/dataverse/ibtb_pazy)
- **论文**：[Buchholz et al. 2022, Proteins 90:1443](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.26325)，doi: 10.1002/prot.26325
- **内容**：按底物（PET、PHB、PCL 等）分类的塑料活性酶蛋白家族库，含序列、结构与活性数据；提供 API 和 DaRUS 数据下载。
- **复用性**：可直接下载 PHB 解聚酶种子序列；有[应用实例：从 PAZy 取 105 个 PETase 作种子做深海宏基因组挖掘](https://pmc.ncbi.nlm.nih.gov/articles/PMC12599313/)。

### A3. PlasticDB（塑料生物降解数据库）

- **入口**：[plasticdb.org](https://plasticdb.org/)、[About 页](http://www.plasticdb.org/about)、[Zenodo 版本记录](https://zenodo.org/records/7217453)
- **论文**：[Gambarini et al. 2022, Database (Oxford) 2022:baac008](https://pubmed.ncbi.nlm.nih.gov/35266524/)（PMID 35266524），doi: 10.1093/database/baac008（开放获取，[PMC9216477](https://pmc.ncbi.nlm.nih.gov/articles/PMC9216477/)，[OUP 论文页](https://academic.oup.com/database/article/doi/10.1093/database/baac008/6546196)）
- **内容**：收录与塑料生物降解相关的微生物和蛋白质（含 HMM 模型供筛选）；论文同时给出宏基因组挖掘方法学（基于隐马尔可夫模型的数据库构建与筛选）。
- **复用性**：可下载，HMM 可直接用于筛选流程。
- 补充文献：塑料生物降解综述提及 PlasticDB 的创建背景（[PMC11856541](https://pmc.ncbi.nlm.nih.gov/articles/PMC11856541/)）。

### A4. PHB 降解基因专用数据库

未检索到独立的"PHB 降解基因数据库"，实际可用的替代资源：

- **BRENDA**：[EC 3.1.1.75 poly(3-hydroxybutyrate) depolymerase](https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.75)
  - [序列条目](https://brenda-enzymes.de/enzyme.php?ecno=3.1.1.75&showtm=0&onlyTable=Sequence)
- **UniProt**：
  - [EC 3.1.1.75 且 reviewed 的条目检索（3 条）](https://www.uniprot.org/uniprotkb?query=(ec:3.1.1.75)%20AND%20(reviewed:true))
  - 按底物检索：[((3R)-hydroxybutanoate)(n) 相关条目（122 条）](https://www.uniprot.org/uniprotkb/?query=%22((3R)-hydroxybutanoate)(n)%22)
- **NCBI CDD**：含 [poly(3-hydroxyalkanoate) depolymerase 保守域模型](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=XP_039073312.1)
- **PDB 结构**：
  - [4BTV（PhaZ7 与 3HB 三聚体复合物结构，Paucimonas lemoignei）](https://www.rcsb.org/structure/4BTV)
  - [8YNV（Bacillus thuringiensis PhaZ 结构）](https://www.rcsb.org/structure/8YNV)

### A5. Pfam / InterPro 家族

- 直接检索**未返回**明确的 Pfam PHB 解聚酶家族页面（未核实到家族号，故不列编号以免误导）。
- 可用替代：NCBI CDD 域模型（见 A4）、ESTHER α/β-水解酶折叠家族体系（[综述 Chatonnet et al. 2023](https://www.sciencedirect.com/science/article/abs/pii/S0009279723003381)）。
- **结论**：建议项目组自行在 InterPro/Pfam 网站按 "polyhydroxybutyrate depolymerase" 复核，并以 ESTHER 家族划分为准。

---

## 二、通用功能注释资源（B 部分）

### B6. KEGG

- **通路**：[map00640 Polyhydroxyalkanoate metabolism](https://www.kegg.jp/entry/pathway+map00640)（PHA 合成与降解通路总图）
- **模块**：
  - [M00843 模块条目](https://www.kegg.jp/entry/M00843)（另见 [module 视图](https://www.kegg.jp/module/smal_M00843)）——检索片段未直接给出该模块标题，**具体内容请以 KEGG 页面为准**
  - ⚠️ **注意**：M00012 实际是乙醛酸循环（glyoxylate cycle），不是 PHA 模块，勿误用（[bio2rdf 的 M00012 记录](https://bio2rdf.org/describe/?uri=http://bio2rdf.org/kegg:M00012)）
- **KO**：
  - [K03821 = phbC/phaC，PHA 合成酶亚基 PhaC](https://www.vandepoelelab.be/plaza/versions/plaza_diatoms_01/Kegg/view/K03821)
  - K03513 编号存在（[dnaconda 的 KO K03513 克隆页](https://dnaconda.riken.jp/search/KO/K03/K03513.html)）；K17745 编号存在（[K17745+R09990](https://www.kegg.jp/entry/K17745+R09990)）——两者具体功能请以 KEGG 页面核实
- 大量 GTDB 基因组已带 KEGG 注释，可在 [MGnify 基因组 KEGG 模块接口](https://www.ebi.ac.uk/metagenomics/api/v1/genomes/MGYG000436303/kegg-module) 等处交叉查询。

### B7. eggNOG / UniProt

- UniProt 有 EC 3.1.1.75 reviewed 条目（见 A4）；eggNOG-mapper 类工具适用于对 GTDB 蛋白批量功能注释。
- 本项目可直接用 eggNOG/UniProt 注释结果做交叉验证。

### B8. CAZy

- 检索未发现 CAZy 中专设的 PHA 解聚酶家族；相关酶（cutinase、α/β-水解酶超家族）在文献中常与 PET 水解酶并列讨论：
  - [PET 降解酶均属 α/β-水解酶超家族，类似 cutinase（Environmental Microbiology 综述）](https://enviromicro-journals.onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.15774)
  - [EC 3.1.1.101 PET 水解酶 BRENDA](https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.101)
- **结论**：CAZy 对本项目帮助有限，酶家族组织以 ESTHER + PAZy + NCBI CDD 为主。

---

## 三、基因组挖掘方法与已发表项目（C 部分）

### C9. PHA 解聚酶基因组挖掘（已发表实例）

- [Streptomyces microflavus DG19 的降解活性与基因组挖掘（Biotechnology for the Environment, 2025）](https://link.springer.com/article/10.1186/s44314-025-00024-7)——基于 Pfam/TIGR/COG 注释挖掘 PHA 解聚酶
- [Burkholderia vietnamiensis 胞外中链 PHA（mcl-PHA）降解基因研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13101487/)（[DOI 记录/DOAJ](https://doaj.org/article/43e1e37fa97743dbbc68aa737db26a30)）
- [Priestia sp. USM5 的 PHA 合成与胞外降解双代谢基因组分析](https://www.sciencedirect.com/science/article/abs/pii/S0141391026000960)
- [Undibacterium sp. KW1/YM2 全基因组中新型胞外 PHA 解聚酶鉴定（PMID 32369496）](https://pesquisa.bvsalud.org/portal/resource/es/mdl-32369496)
- 相关背景：[Photobacterium ganghwense C2.2 全基因组（PHA 生产菌，含 phaR/Q/F 调控与降解基因讨论）](https://onlinelibrary.wiley.com/doi/full/10.1002/mbo3.1182)

### C10. GTDB 可下载资源 ★ 项目数据基础

- **数据下载**：
  - [data.gtdb.ecogenomic.org/releases/latest/（FILE_DESCRIPTIONS.txt 文件清单说明）](https://data.gtdb.ecogenomic.org/releases/latest/FILE_DESCRIPTIONS.txt)
  - 当前目录示例：[release232/232.0/genomic_files_reps/](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/)、[release226/226.0/genomic_files_reps/](https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_reps/)——代表基因组蛋白文件（gtdb_proteins_aa_reps.tar.gz）与 genome bundles 均在此
- **规模/统计**：
  - [GTDB R220 统计页](https://gtdb.ecogenomic.org/stats/r220)、[R214 统计页](https://gtdb.ecogenomic.org/stats/r214)
  - 最新论文：[GTDB release 10：715,230 细菌 + 17,245 古菌基因组（Parks et al. 2025）](https://www.semanticscholar.org/paper/GTDB-release-10%3A-a-complete-and-systematic-taxonomy-Parks-Chaumeil/7a716bac654f1366a47f47186d344152ca8d5e88)（[AAU 记录](https://vbn.aau.dk/en/publications/gtdb-release-10-a-complete-and-systematic-taxonomy-for-715230-bac/)、[中文介绍](https://www.ebiotrade.com/newsf/2025-10/20251023083213816.htm)）
- **蛋白文件使用实战**：
  - [GTDB Forum：gtdb_proteins_aa_reps.tar.gz 的 fasta header→分类学映射](https://forum.gtdb.ecogenomic.org/t/taxonomy-lookup-for-fasta-headers-from-gtdb-proteins-aa-reps-tar-gz/599/3)
  - [GTDB Forum：如何获取 GTDB 参考基因组序列](https://forum.gtdb.ecogenomic.org/t/how-to-obtain-reference-genome-sequences-in-gtdb-database/534/3)
  - [James Lingford：把 GTDB faa 蛋白文件转成带分类学的 DIAMOND 数据库](https://www.jameslingford.com/blog/gtdb-to-diamond-taxonomy-database/)（可直接套用）
- 其他生态位资源：
  - [sourmash GTDB RS220 集合](https://sourmash.readthedocs.io/en/stable/databases-md/gtdb220.html)
  - [chem16S R 包内置 GTDB_220 氨基酸组成参考](https://rdrr.io/cran/chem16S/src/inst/RefDB/GTDB_220/genome_AA.R)
  - 已有研究以 GTDB r214 代表基因组为数据基础做代谢分析（[Cell Host & Microbe 2025 示例](https://www.sciencedirect.com/science/article/pii/S193131282500280X)）

### C11–C12. 大规模宏基因组酶挖掘流程与文章 ★ 方法学标杆

- **[Zrimec et al. 2021, mBio 12:e02155-21（PMID 34700384）](https://pubmed.ncbi.nlm.nih.gov/34700384/)**："Plastic-Degrading Potential across the Global Microbiome Correlates with Recent Pollution Trends"
  - [期刊版（mBio）](https://journals.asm.org/doi/10.1128/mbio.02155-21)——明确写出用 **HMMER v3.3 的 hmmsearch** 在宏基因组中筛同源序列
  - [预印本（bioRxiv 2020.12.13.422558）](https://www.biorxiv.org/content/10.1101/2020.12.13.422558v2.full.pdf)
  - 新闻稿称其分析规模达数千万基因级：[Discovery 报道"Thousands of Bacterial Enzymes"](https://www.discovery.com/science/thousands-of-bacterial-enzymes-hold-the-key-to-plastic-pollution)、[Innovation News Network](https://www.innovationnewsnetwork.com/enzymes-degrade-plastic-increasing-association-pollution/16599/)、[Smithsonian](https://www.smithsonianmag.com/smart-news/scientists-say-plastic-degrading-enzymes-are-increasing-in-response-to-pollution-180979250/)
  - 核心思路：构建塑料降解酶 HMM → 全宏基因组 hmmsearch → 与污染水平关联
- **[Predicting the plastic biodegradation potential within microbial lineages and across global ecosystems（Microbial Genomics 2025, mgen.0.001814）](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001814)**，配套数据集 [Zenodo 15480170](https://zenodo.org/records/15480170)——在"谱系 × 全球生态系统"尺度预测塑料降解潜力，与本项目"GTDB 全库筛选"思路高度同源
- **[PlasticEnz：同源 + 机器学习集成鉴定宏基因组塑料降解酶](https://www.biorxiv.org/content/10.1101/2025.10.28.685028v1)**（预印本；[PLOS Comput Biol 图表页 pcbi.1013892](https://journals.plos.org/ploscompbiol/article/figure?id=10.1371/journal.pcbi.1013892.g004)；[中文介绍（ebiotrade）](https://news.ebiotrade.com/2026-1/20260130085459418.htm)；[Semantic Scholar 记录](https://www.semanticscholar.org/paper/PlasticEnz%3A-An-integrated-database-and-screening-to-Krzyn%C3%B3wek-Snoeks/31ade9bc380fd46c4dc253462b1647bd30f244c4)）
- **plastedma**：[宏基因组塑料降解酶检测/注释工作流（NSCCN）](https://gitlink.org.cn/NSCCN/plastedma)，[bioconda 包](https://anaconda.org/bioconda/plastedma/manage)
- **PDETool**：[塑料降解酶工具（GitHub, ozefreitas/PDETool）](https://github.com/ozefreitas/PDETool)
- 其他已发表挖掘：
  - [海洋宏基因组 PHA 合成基因挖掘（Discover Oceans 2026）](https://link.springer.com/article/10.1007/s44289-026-00132-6)
  - [用 PAZy 种子做深海未培养微生物塑料降解挖掘（PMC12599313）](https://pmc.ncbi.nlm.nih.gov/articles/PMC12599313/)
  - [宏基因组挖掘新型酯酶（J Microbiol Biotechnol 2015, PMID 25502823）](https://pubmed.ncbi.nlm.nih.gov/25502823/)

### C13. HMM 筛选方法

- Zrimec 2021 的 HMMER v3.3 hmmsearch 流程（见 C12）
- PlasticDB 提供塑料降解酶 HMM；PHA Depolymerase Engineering Database 提供家族 HMM 分配
- 实例：[Family VIII 酯酶论文用"已报道塑料降解酶的 HMM"鉴定候选酶](https://www.biorxiv.org/content/10.1101/2025.09.01.670660v1.full)
- **结论**：**HMM（HMMER）+ 种子库（PAZy/ESTHER/PlasticDB）+ 催化残基验证**是当前公认范式

---

## 四、分析与可视化工具（D 部分）

### D14. GTDB-Tk

- **仓库**：[github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/gtdbtk)
- **文档**：[ecogenomics.github.io/GTDBTk](https://ecogenomics.github.io/GTDBTk/)、[README](https://raw.githubusercontent.com/Ecogenomics/GTDBTk/120d743952f2d5c6e0a4bbe90410219d4448dce2/README.md)
- **配套**：
  - [可复现的 GTDB-Tk 安装+分类+汇总流水线（gmboowa/gtdbtk-setup-and-classification）](https://github.com/gmboowa/gtdbtk-setup-and-classification)
  - [KBase 的 GTDB-Tk classify_wf 应用](https://kbase.us/applist/apps/kb_gtdbtk/run_kb_gtdbtk_classify_wf/release)
  - [Bactopia 的 gtdb subworkflow](https://bactopia.io/developers/subworkflows/gtdb)
- **用途**：本项目若涉及自组装/新基因组，可用 GTDB-Tk 分类；若只用 GTDB 现有代表基因组，则直接利用其自带分类元数据即可。

### D15. 蛋白质系统发育分析流程

- **标准流程**：MAFFT 多序列比对 → IQ-TREE 2 最大似然建树 → ETE3 树操作/可视化（配套 ETE 工具集：[etetoolkit](https://github.com/etetoolkit/ete)）
- [**PEGP（protein evolutionary genomics pipeline）**](https://github.com/stovc/pegp/)：同源搜索 + 系统发育树 + 域架构注释 + 基因组共线性 + 系统谱系模式/旁系同源鉴定——与本项目"基因家族尺度分析"高度契合
- **anvio**（[anvio.org](https://anvio.org/)，[工具介绍](https://www.bioinformaticshome.com/db/tool/anvio/)）：适合泛基因组/比较基因组与基因上下文分析，可作补充
- **BiG-SCAPE**（[BiG-SCAPE 2.0 论文（Nat Commun 2026）](https://www.nature.com/articles/s41467-026-68733-5)、[主页](https://bigscape-corason.secondarymetabolites.org/about/)）：用于**次级代谢基因簇（BGC）聚类**，**不适合**作为 PHB 降解酶家族系统发育的主工具，仅当关注基因簇结构时可参考

---

## 五、方法学建议

1. **数据获取（GTDB 侧）**：从 [data.gtdb.ecogenomic.org/releases/latest/](https://data.gtdb.ecogenomic.org/releases/latest/FILE_DESCRIPTIONS.txt) 下载 `gtdb_proteins_aa_reps.tar.gz`（代表基因组蛋白）及对应 `bac120_*`/`ar53_*` 元数据（分类、GC%、基因组质量）；按 [Forum 线程 599](https://forum.gtdb.ecogenomic.org/t/taxonomy-lookup-for-fasta-headers-from-gtdb-proteins-aa-reps-tar-gz/599/3) 做 header→GTDB 分类学映射，或直接套用 [James Lingford 的 DIAMOND 化流程](https://www.jameslingford.com/blog/gtdb-to-diamond-taxonomy-database/)（可顺带得到"酶→物种/谱系"的注释能力）。规模参照 GTDB release 10：**约 71.5 万细菌 + 1.7 万古菌代表基因组**。
2. **种子与 HMM 构建**：以 **PHA Depolymerase Engineering Database（ESTHER 家族页）** 的家族序列 + **PAZy** 的 PHB 解聚酶条目 + **PlasticDB** 的 HMM 为种子；用 HMMER 3.x 的 `hmmbuild`/`hmmsearch` 全库扫描（复刻 [Zrimec 2021](https://journals.asm.org/doi/10.1128/mbio.02155-21) 的 HMMER v3.3 做法）。速度优先时用 MMseqs2/DIAMOND 预筛 + hmmsearch 精筛。
3. **候选验证与分类**：按 e-value/覆盖度过滤后，检查 α/β-水解酶催化三联体（Ser-Asp-His）与 lipase box（GxSxG），用 NCBI CDD、KEGG（[map00640](https://www.kegg.jp/entry/pathway+map00640)、EC 3.1.1.75、K03821）、UniProt 交叉验证；按结构域组织区分**胞外型（I–VIII 型，常带底物结合域）与胞内 PhaZ**。
4. **系统发育**：MAFFT（或 Clustal Omega 粗对齐）→ IQ-TREE 2（`-m LG+G4` 类模型，必要时按家族分建）→ ETE3 可视化，用 GTDB 分类学给叶节点标注/着色；可用 PEGP 补域架构与 phyletic pattern。
5. **基因组背景分析**：对命中基因组检查 PHA 合成（phaC/phaA/phaB）、phasin（phaP）、调控（phaR/Q/F）共现，判断降解基因所在代谢背景（参考 [Priestia USM5](https://www.sciencedirect.com/science/article/abs/pii/S0141391026000960) 与 [Photobacterium ganghwense](https://onlinelibrary.wiley.com/doi/full/10.1002/mbo3.1182) 的做法；phasin/颗粒相关蛋白背景见 [PMC4029623](https://pmc.ncbi.nlm.nih.gov/articles/PMC4029623/)）。
6. **文献对照**：与 [Zrimec 2021](https://pubmed.ncbi.nlm.nih.gov/34700384/)、[Microbial Genomics 2025](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001814)、[PlasticEnz](https://www.biorxiv.org/content/10.1101/2025.10.28.685028v1) 的结果交叉比较，突出 GTDB 全覆盖 + 分类学尺度系统发育分析的增量价值。
7. **许可**：GTDB、PAZy（含 API/镜像）、PlasticDB、ESTHER、KEGG（学术用途）均可免费使用，引用相应原始文献。

---

## 六、参考来源 URL 汇总

**专用数据库**
- [Knoll 2009 (BMC Bioinformatics, PubMed)](https://pubmed.ncbi.nlm.nih.gov/19296857/)
- [Knoll 2009 全文 (Springer)](https://link.springer.com/article/10.1186/1471-2105-10-89?optIn=true)
- [ESTHER 主页/论文记录](https://bioweb.supagro.inrae.fr/ESTHER/paper/Knoll_2009_BMC.Bioinformatics_10_89)
- [ESTHER 家族 Esterase_phb_PHAZ](https://bioweb.supagro.inrae.fr/ESTHER/family/Esterase_phb_PHAZ)
- [ESTHER 家族 PHAZ7_phb_depolymerase](https://bioweb.supagro.inrae.fr/ESTHER/family/PHAZ7_phb_depolymerase)
- [Database Commons 3446](https://ngdc.cncb.ac.cn/databasecommons/database/id/3446)
- [PAZy 主页 (cbl.uni-stuttgart.de)](https://www.cbl.uni-stuttgart.de/doku.php?id=start)
- [PAZy 论文 (Wiley, Proteins)](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.26325)
- [PAZy API 文档](https://www.pazy.eu/api-docs)
- [PAZy 数据镜像 (DaRUS ibc_tbc_PAZy)](https://darus.uni-stuttgart.de/dataverse/ibc_tbc_PAZy)
- [PlasticDB](https://plasticdb.org/)
- [PlasticDB About](http://www.plasticdb.org/about)
- [PlasticDB 论文 (PubMed)](https://pubmed.ncbi.nlm.nih.gov/35266524/)
- [PlasticDB 论文 (PMC9216477)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9216477/)
- [PlasticDB Zenodo 记录](https://zenodo.org/records/7217453)
- [BRENDA EC 3.1.1.75](https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.75)
- [UniProt EC 3.1.1.75 reviewed 检索](https://www.uniprot.org/uniprotkb?query=(ec:3.1.1.75)%20AND%20(reviewed:true))
- [NCBI CDD poly(3-hydroxyalkanoate) depolymerase](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=XP_039073312.1)
- [PDB 4BTV (PhaZ7-3HB 复合物)](https://www.rcsb.org/structure/4BTV)
- [PDB 8YNV (Bt PhaZ)](https://www.rcsb.org/structure/8YNV)
- [ESTHER α/β-水解酶综述 2023](https://www.sciencedirect.com/science/article/abs/pii/S0009279723003381)

**功能注释**
- [KEGG map00640](https://www.kegg.jp/entry/pathway+map00640)
- [KEGG M00843](https://www.kegg.jp/entry/M00843)
- [KEGG M00843 (module 视图)](https://www.kegg.jp/module/smal_M00843)
- [K03821 phaC (PLAZA diatoms)](https://www.vandepoelelab.be/plaza/versions/plaza_diatoms_01/Kegg/view/K03821)
- [K17745 (KEGG)](https://www.kegg.jp/entry/K17745+R09990)
- [α/β-水解酶与 PETase 关系综述](https://enviromicro-journals.onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.15774)
- [BRENDA EC 3.1.1.101 (PET 水解酶)](https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.101)

**基因组挖掘 / 已发表项目**
- [Zrimec 2021 (PubMed, mBio)](https://pubmed.ncbi.nlm.nih.gov/34700384/)
- [Zrimec 2021 (mBio 期刊版)](https://journals.asm.org/doi/10.1128/mbio.02155-21)
- [Zrimec 预印本 (bioRxiv)](https://www.biorxiv.org/content/10.1101/2020.12.13.422558v2.full.pdf)
- [Microbial Genomics 2025 (mgen.0.001814)](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001814)
- [Zenodo 15480170 (配套数据集)](https://zenodo.org/records/15480170)
- [PlasticEnz 预印本 (bioRxiv)](https://www.biorxiv.org/content/10.1101/2025.10.28.685028v1)
- [PlasticEnz PLOS Comput Biol 图](https://journals.plos.org/ploscompbiol/article/figure?id=10.1371/journal.pcbi.1013892.g004)
- [plastedma (NSCCN)](https://gitlink.org.cn/NSCCN/plastedma)
- [plastedma (bioconda)](https://anaconda.org/bioconda/plastedma/manage)
- [PDETool (GitHub)](https://github.com/ozefreitas/PDETool)
- [S. microflavus DG19 (2025)](https://link.springer.com/article/10.1186/s44314-025-00024-7)
- [B. vietnamiensis mcl-PHA 降解 (PMC13101487)](https://pmc.ncbi.nlm.nih.gov/articles/PMC13101487/)
- [Priestia USM5 双 PHA 代谢](https://www.sciencedirect.com/science/article/abs/pii/S0141391026000960)
- [Undibacterium KW1/YM2 新型解聚酶 (PMID 32369496)](https://pesquisa.bvsalud.org/portal/resource/es/mdl-32369496)
- [Photobacterium ganghwense C2.2 全基因组](https://onlinelibrary.wiley.com/doi/full/10.1002/mbo3.1182)
- [海洋 PHA 合成基因宏基因组挖掘 (Discover Oceans)](https://link.springer.com/article/10.1007/s44289-026-00132-6)
- [深海未培养微生物塑料降解挖掘 (PMC12599313)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12599313/)
- [Family VIII 酯酶 HMM 筛选 (bioRxiv)](https://www.biorxiv.org/content/10.1101/2025.09.01.670660v1.full)
- [宏基因组挖掘新型酯酶 (PMID 25502823)](https://pubmed.ncbi.nlm.nih.gov/25502823/)

**GTDB**
- [GTDB 数据目录 latest (FILE_DESCRIPTIONS.txt)](https://data.gtdb.ecogenomic.org/releases/latest/FILE_DESCRIPTIONS.txt)
- [release226 genomic_files_reps](https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_reps/)
- [release232 genomic_files_reps](https://data.gtdb.ecogenomic.org/releases/release232/232.0/genomic_files_reps/)
- [GTDB R220 统计](https://gtdb.ecogenomic.org/stats/r220)
- [GTDB R214 统计](https://gtdb.ecogenomic.org/stats/r214)
- [GTDB release 10 论文 (Semantic Scholar)](https://www.semanticscholar.org/paper/GTDB-release-10%3A-a-complete-and-systematic-taxonomy-Parks-Chaumeil/7a716bac654f1366a47f47186d344152ca8d5e88)
- [GTDB Forum: 蛋白文件 header 映射](https://forum.gtdb.ecogenomic.org/t/taxonomy-lookup-for-fasta-headers-from-gtdb-proteins-aa-reps-tar-gz/599/3)
- [GTDB Forum: 获取参考基因组序列](https://forum.gtdb.ecogenomic.org/t/how-to-obtain-reference-genome-sequences-in-gtdb-database/534/3)
- [GTDB faa→DIAMOND 实战 (James Lingford)](https://www.jameslingford.com/blog/gtdb-to-diamond-taxonomy-database/)
- [sourmash GTDB RS220](https://sourmash.readthedocs.io/en/stable/databases-md/gtdb220.html)
- [chem16S GTDB_220 参考](https://rdrr.io/cran/chem16S/src/inst/RefDB/GTDB_220/genome_AA.R)

**工具**
- [GTDB-Tk 仓库 (Ecogenomics)](https://github.com/Ecogenomics/gtdbtk)
- [GTDB-Tk 文档](https://ecogenomics.github.io/GTDBTk/)
- [gmboowa GTDB-Tk 设置与分类流水线](https://github.com/gmboowa/gtdbtk-setup-and-classification)
- [KBase GTDB-Tk classify_wf](https://kbase.us/applist/apps/kb_gtdbtk/run_kb_gtdbtk_classify_wf/release)
- [Bactopia gtdb subworkflow](https://bactopia.io/developers/subworkflows/gtdb)
- [PEGP (protein evolutionary genomics pipeline)](https://github.com/stovc/pegp/)
- [anvio](https://anvio.org/)
- [BiG-SCAPE 2.0 (Nat Commun)](https://www.nature.com/articles/s41467-026-68733-5)
- [BiG-SCAPE CORASON 主页](https://bigscape-corason.secondarymetabolites.org/about/)

**综述 / 背景**
- [Microbial degradation of polyhydroxyalkanoates（Jendrossek & Handrick, 2002）](https://pubmed.ncbi.nlm.nih.gov/12213937/)
- [Make it or break it: PHA synthase and depolymerase 综述（J Polym Environ 2024）](https://link.springer.com/article/10.1007/s10924-024-03474-4)
- [PHB 解聚酶在废物管理中的作用（J Environ Manage 2025）](https://www.sciencedirect.com/science/article/abs/pii/S0301479725009016)
- [PhaP phasins 在 PHB 积累中的作用（BMC Microbiology 2013）](https://pmc.ncbi.nlm.nih.gov/articles/PMC4029623/)

---

## 七、给主代理的关键结论

1. 没有现成的"GTDB PHB 降解基因筛选"成品流程，但**所有积木都已存在且可免费复用**：GTDB 蛋白文件 + ESTHER/PAZy/PlasticDB 种子与 HMM + HMMER + MAFFT/IQ-TREE/ETE3。
2. 最值得模仿的两个先例：**Zrimec 2021（HMM 全库扫描，mBio）** 和 **Microbial Genomics 2025（谱系 × 生态系统尺度）**。
3. GTDB release 10 已有约 73 万代表基因组（71.5 万细菌 + 1.7 万古菌），全部蛋白可在 `gtdb_proteins_aa_reps.tar.gz` 一次性获取，规模可控（MMseqs2/DIAMOND 即可胜任）。
4. 注释与分类锚点：KEGG map00640、EC 3.1.1.75、ESTHER 家族分类体系、NCBI CDD 域模型。
5. 遗留待核实项（后续可在官网确认）：Pfam PHB 解聚酶家族编号；KEGG M00843 模块标题；K03513/K17745 的具体功能定义。
