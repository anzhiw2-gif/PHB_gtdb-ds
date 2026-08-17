# PHA/PHB 降解酶数据库家族体系与古菌 PHA 降解文献调研报告

> 生成日期：2025 年（本报告由 OpenAlex 检索 + web 检索综合而成）
> 检索工具：OpenAlex CLI（`literature-search-openalex` skill，脚本 `C:\Users\HUAWEI\.agents\skills\literature-search-openalex\scripts\openalex_cli.py`）+ web_search
> OpenAlex 原始数据文件：`D:\PHB_gtdb-ds\research\openalex\comp_search*.json`、`comp_work_*.json`

---

## 一、数据库家族体系总表

### 1.1 PhaDED（PHA Depolymerase Engineering Database）— 核心分类体系

**关键文献**：Knoll M, Hamm TM, Wagner F, Martinez V, Pleiss J. *The PHA Depolymerase Engineering Database: A systematic analysis tool for the diverse family of polyhydroxyalkanoate (PHA) depolymerases.* **BMC Bioinformatics** 2009, 10:89. doi:10.1186/1471-2105-10-89（引用 162 次）— [PubMed](https://pubmed.ncbi.nlm.nih.gov/19296857/)、[PMC2666664](https://pmc.ncbi.nlm.nih.gov/articles/PMC2666664/)、[ESTHER 收录页](https://bioweb.supagro.inrae.fr/ESTHER/paper/Knoll_2009_BMC.Bioinformatics_10_89)

| 项目 | 内容 |
|---|---|
| 数据库地址 | http://www.ded.uni-stuttgart.de（现已并入 ESTHER / 斯图加特大学 α/β-水解酶系） |
| 收录规模 | **587 条 PHA 解聚酶序列** |
| 分类层级 | **8 个超家族（superfamilies）→ 38 个同源家族（homologous families）** |
| 判定依据 | 序列相似性（结合保守结构域），非仅 EC 号 |
| 提供的资源 | 每家族多序列比对（MSA）+ **profile hidden Markov models（HMM）** + 功能相关残基注释 |
| 用途 | 基因组 in silico 鉴定新 PHA 解聚酶、家族分类、生化特性预测、酶工程 |

**PhaDED 家族命名/划分逻辑**（依据摘要与后续文献）：
- 所有 PHA 解聚酶共享 **α/β-水解酶折叠（alpha/beta-hydrolase fold）与催化三联体（Ser-His-Asp）**，与脂肪酶/酯酶同源。
- 超家族按"胞外 scl-PHA / 胞外 mcl-PHA / 胞内（颗粒结合）等"的生化类别 + 序列系统发育切分；家族内序列高度同源。
- PhaDED 的分类是构建本任务 HMM 种子库的**权威参照**。

### 1.2 ESTHER 数据库（ESTerases and alpha/beta-Hydrolase Enzymes and Relatives）

- 站点：https://bioweb.supagro.inrae.fr/ESTHER/（INRAE / 蒙彼利埃；与斯图加特数据库合并维护 α/β-水解酶系）
- 综述：Lenfant N, et al. *The ESTHER database on alpha/beta hydrolase fold proteins – An overview of recent developments.*（[ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0009279723003381?fr=RR-2&ref=pdf_download&rr=9cdfec8fac079cf3)、[R Discovery](https://discovery.researcher.life/topic/hydrolase-fold/3927753?page=1&topic_name=Hydrolase+Fold)）
- PHA 解聚酶相关 ESTHER 家族（已核实 URL 存在的家族页面）：
  - **`Esterase_phb_PHAZ`** — https://bioweb.supagro.inrae.fr/ESTHER/family/Esterase_phb_PHAZ
  - **`PHAZ7_phb_depolymerase`** — https://bioweb.supagro.inrae.fr/ESTHER/family/PHAZ7_phb_depolymerase（PhaZ7 来自 *Paucimonas lemoignei*，耐热耐碱 PHB 解聚酶）
  - 另有 `Esterase_phb`（广义 PHB 酯酶组）；Knoll 2009 论文页由 ESTHER 收录。
- ESTHER 家族命名规则：`<家族名>_<底物/功能>`（如 `Esterase_phb_PHAZ`、`PHAZ7_phb_depolymerase`），家族下挂 gene_locus 与序列。

### 1.3 PAZy（Plastics-Active Enzymes Database，塑料活性酶数据库）

- 站点：https://pazy.eu / https://www.cbl.uni-stuttgart.de/doku.php?id=start（斯图加特大学 Pleiss 组）
- 关键文献：Buchholz PCF, et al. *Plastics degradation by hydrolytic enzymes: The plastics-active enzymes database—PAZy.* **Proteins** 2022, 90:1443–1456. doi:10.1002/prot.26325（[Wiley](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.26325)、[ESTHER 收录](https://bioweb.supagro.inrae.fr/ESTHER/paper/Buchholz_2022_Proteins_90_1443)）
- 收录约 **110 个经验证的塑料降解酶**（含 PHB 解聚酶、PETase、MHETase、角质酶等；据 [Microbial Biotechnology 2023 综述](https://www.ovid.com/journals/micbt/fulltext/10.1111/1751-7915.14135~microbial-enzymes-will-offer-limited-solutions-to-the-global)）；可按底物（PHB、PHA、PET 等）浏览条目：[PAZy 按底物列表](https://www.pazy.eu/plastics/pa?substrate=4&limit=50&offset=50)
- 意义：PHA 解聚酶（EC 3.1.1.75）与 PHB 降解相关酶在 PAZy 中作为"塑料活性酶"家族之一收录，可交叉验证 UniProt 条目与实验证据。

### 1.4 BRENDA / 其他酶学数据库 — EC 编号

| EC 编号 | 名称 | 备注 |
|---|---|---|
| **EC 3.1.1.75** | poly(3-hydroxybutyrate) depolymerase（胞外 PHB 解聚酶） | 最主要收录条目；BRENDA [EC 3.1.1.75](https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.75&UniProtAcc=P12625&OrganismID=5160)；[Wikipedia](https://en.wikipedia.org/?curid=14457442) |
| **EC 3.1.1.76** | poly(3-hydroxyoctanoate) depolymerase（胞内 mcl-PHA 解聚酶） | BRENDA [EC 3.1.1.76](https://brenda-enzymes.org/enzyme.php?ecno=3.1.1.76)；MetaCyc [3.1.1.76-RXN](https://guest:guest@microcyc.genoscope.cns.fr/META/NEW-IMAGE?type=REACTION&object=3.1.1.76-RXN)；solcyc [EC-3.1.1.76](https://solcyc.solgenomics.net/META/NEW-IMAGE?type=EC-NUMBER&object=EC-3.1.1.76&&orgids=(LYCO)) |
| EC 3.1.1.22 | 3-hydroxybutyrate oligomer hydrolase（3-羟基丁酸寡聚体水解酶） | BRENDA [EC 3.1.1.22](https://www.brenda-enzymes.de/search_result.php?T%5b1%5d=1&T%5b9%5d=2&W%5b1%5d=3.1.1.22&a=30&l=10,10&orderDesc=1&orderByHTMLField=Organism)；胞内 PHB 降解关键酶之一 |

**EC 3.1.1.75 收录物种（BRENDA/UniProt 核实）**：*Ralstonia pickettii*（UniProt P12625，经典胞外 PHB 解聚酶）、*Pseudomonas lemoignei*（PhaZ1–PhaZ5 系统）、*Comamonas* sp.、*Alcaligenes faecalis*、*Cupriavidus/Ralstonia eutropha*（胞内 i-1/i-2）、*Bacillus* spp.、*Streptomyces* spp.、真菌 *Penicillium funiculosum*（PDB 2D80）等。UniProt 检索式 `(ec:3.1.1.75) AND (reviewed:true)` 返回 3 条 reviewed 主条目（[UniProt 链接](https://www.uniprot.org/uniprotkb?query=(ec:3.1.1.75)%20AND%20(reviewed:true))）。

### 1.5 Pfam / InterPro / CDD / TIGRFAM 家族（含编号核实）

| 数据库 | 家族/条目 | 编号 | 说明 |
|---|---|---|---|
| CDD（NCBI） | **Esterase_PHB** | pfam10503（CDD 链接：[pfam10503](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam10503)） | Pfam PF10503 对应的 CDD 家族，PHB 酯酶/解聚酶催化域 |
| CDD（NCBI） | esterase_phb | **TIGR01840**（[CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/TIGR01840)、[TIGRFAM HMM](http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=TIGR01840)） | TIGRFAM 的 PHB 酯酶家族 HMM |
| CDD（NCBI） | **PHB_depoly_PhaZ** | **TIGR01849**（[TIGRFAM HMM](http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=TIGR01849)、[CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?ascbin=8&maxaln=10&seltype=2&uid=TIGR01849)） | TIGRFAM 的 PHA 解聚酶 PhaZ 家族 HMM |
| SCOP/SUPERFAMILY | PHB depolymerase-like family | SCOP family b.69（[supfam 159747](https://supfam.org/scop/159747)、[SCOP 页面](https://140.114.98.75/scop/data/scop.b.d.baa.b.dh.html)） | α/β-水解酶折叠下的 PHB 解聚酶样家族 |
| SCOPe | PhaZ7 结构域 | PDB 2VTV（[PhaZ7 depolymerase from *Paucimonas lemoignei*](https://pdbj.org/mine/structural_details/2vtv)） | 首个胞外 scl PHB 解聚酶晶体结构 |
| PDB | PHB 解聚酶结构 | 2D80（*Penicillium funiculosum*，[PDBsum](http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=2D80)）、8DAJ（*Lihuaxuella thermophila* 耐热混杂 PHB 解聚酶，[RCSB](https://www.rcsb.org/structure/8DAJ)）、8YNW（*Bacillus thuringiensis* PhaZ，[RCSB](https://www.rcsb.org/structure/8YNW)） | 结构佐证家族划分 |

> ⚠️ 说明：InterPro 中对"polyhydroxybutyrate depolymerase"尚无单一专用 IPR 条目（相关条目散落在 α/β-水解酶折叠与酯酶家族中）；**Pfam PF10503（CDD 收录为 Esterase_PHB）与 TIGRFAM TIGR01840 / TIGR01849 是最可直接下载 HMM 的家族模型**，建议作为 HMM 种子库的主要来源之一。公开检索中未核实到独立"IPR-PHB depolymerase"专属条目，报告以 PFAM/CDD/TIGRFAM/SCOP 编号为准。

---

## 二、古菌 PHA 降解文献全集（OpenAlex 检索 + web 核实）

### 2.1 核心古菌降解文献（按重要性排序）

| # | 文献 | 链接/DOI | 方法/结论要点 |
|---|---|---|---|
| 1 | **Liu G, et al. A Patatin-Like Protein Associated with the PHA Granules of *Haloferax mediterranei* Acts as an Efficient Depolymerase in the Degradation of Native PHA.** AEM 2015, 81(9):3029–38. | doi:10.1128/aem.04269-14；[PubMed 25710370](https://pubmed.ncbi.nlm.nih.gov/25710370/)；[ASM](https://journals.asm.org/doi/full/10.1128/aem.04269-14?af=R&mi=grqqmx&ConceptID=512878&target=topic)；[PMC4393451](https://pmc.ncbi.nlm.nih.gov/articles/PMC4393451/) | 首次鉴定**古菌 PHA 降解酶**：颗粒结合 patatin 样解聚酶 **PhaZh1**（HFX_6463）。体外水解天然 PHA 颗粒（nPHB/nPHBV），主产物 3HB。定点突变证实 **Gly16、Ser47（经典 lipase box G-X-S47-X-G）、Asp195** 必需。`phaZh1` 与 `bdhA`（3HB 脱氢酶）构成基因簇 HFX_6463–6464；**phaZh1 敲除不影响胞内 PHA mobilization → 存在替代降解通路**。 |
| 2 | **Han J, et al. Enoyl-CoA hydratase mediates polyhydroxyalkanoate mobilization in *Haloferax mediterranei*.** Sci Rep 2016, 6:24015. | doi:10.1038/srep24015；[Nature](https://www.nature.com/articles/srep24015)、[Springer 页](https://link.springer.com/article/10.1038/srep24015) | 5 个 (R)-特异性烯脂酰-CoA 水合酶（PhaJ1–PhaJ5）中仅**颗粒结合的 PhaJ1** 参与 PHA mobilization：催化 (R)-3-羟基脂酰-CoA 脱水为烯脂酰-CoA，**衔接 β-氧化**。抑制 β-氧化即抑制 PHA 降解。**96% 含 phaJ 的古菌同时拥有 phaC（PHA 合酶）与全套 β-氧化基因 → PHA 经 β-氧化 mobilization 在嗜盐古菌中普遍**。 |
| 3 | **Lu Q, et al. Identification of Polyhydroxyalkanoates in *Halococcus* and Other Haloarchaeal Species.** Appl Microbiol Biotechnol 2010, 87:1301–1312. | doi:10.1007/s00253-010-2611-6；[Springer](https://link.springer.com/article/10.1007/s00253-010-2611-6) | 20 株嗜盐古菌 PHA 检测（Sudan Black/Nile Blue/Nile Red + TEM + ¹H-NMR）：首次报告 *Halococcus* 各型、*Halorubrum*、嗜碱古菌 *Natronobacterium/Natronococcus*、*Halobacterium noricense* 产 PHB/PHBV；而 *Hbt. salinarum* NRC-1/R1 与 *Hfx. volcanii* 不产。**提供了古菌 PHA 生产（含降解背景）的谱系清单**。 |
| 4 | **Lu Q, et al. Wide Distribution among Halophilic Archaea of a Novel Polyhydroxyalkanoate Synthase Subtype with Homology to Bacterial Type III Synthases.** AEM 2010, 76(22):7811–19. | doi:10.1128/aem.01117-10；[ASM](https://journals.asm.org/doi/10.1128/aem.01117-10) | 28 株 15 属嗜盐古菌中 18 株产 PHB/PHBV；古菌 PHA 合酶为**新亚型 IIIA（type III）**，与细菌 IIIB 明显分开。证明 *Halobacteriaceae* 内 PHA 代谢（合成+降解基因座）广泛分布。 |
| 5 | **Roh H, et al. Bioinformatics Analysis of Metabolism Pathways of Archaeal Energy Reserves.** Sci Rep 2018, 8:1034. | doi:10.1038/s41598-018-37768-0；[Nature](https://www.nature.com/articles/s41598-018-37768-0)、[PMC6355812](https://pmc.ncbi.nlm.nih.gov/articles/PMC6355812/) | 系统盘点古菌储能物质（糖原/PHA 等）代谢通路基因分布——**可用于确认哪些古菌门含 PHA 合成/降解基因**（报告指出 PHA 代谢主要限于 *Euryarchaeota* 的部分类群，尤其 *Halobacteria*）。 |
| 6 | **Arpigny JL, Jendrossek D. A novel heat-stable lipolytic enzyme from *Sulfolobus acidocaldarius* DSM 639 displaying similarity to polyhydroxyalkanoate depolymerases.** FEMS Microbiol Lett 1998, 167(1):69–73. | doi:10.1016/S0378-1097(98)00375-9；[ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0378109798003759)、[Semantic Scholar](https://www.semanticscholar.org/paper/A-novel-heat-stable-lipolytic-enzyme-from-DSM-639-Arpigny-Jendrossek/a432113b295044955fcafa30d6006e37d4b86aac) | **非嗜盐古菌中的直接证据**：*Sulfolobus acidocaldarius*（泉古菌门/Crenarchaeota，嗜热嗜酸）的热稳定脂解酶序列与 PHA 解聚酶（Esterase_phb 家族）显著相似——提示泉古菌中存在 PHA 解聚酶同源序列（是否具 PHA 水解活性待验证）。 |
| 7 | **Jendrossek D, Handrick R. Microbial Degradation of Polyhydroxyalkanoates.** Annu Rev Microbiol 2002, 56:403–432. | doi:10.1146/annurev.micro.56.012302.160838；[PubMed 12213937](https://pubmed.ncbi.nlm.nih.gov/12213937/)；[Semantic Scholar](https://www.semanticscholar.org/paper/Microbial-degradation-of-polyhydroxyalkanoates.-Jendrossek-Handrick/a67a4ff2b41e3328a23e08cfad4045af782df6d5) | 权威综述（引用 692 次）：胞外/胞内解聚酶分类、催化机理、颗粒结合蛋白；指出**古菌 PHA 降解长期是空白**，直到 2015 年 PhaZh1 填补。 |
| 8 | **Poltronieri P, Kumar P. Polyhydroxyalkanoate Biosynthesis at the Edge of Water Activity – Haloarchaea as Biopolyester Factories.** Bioengineering 2019, 6(2):34. | doi:10.3390/bioengineering6020034；[MDPI](https://www.mdpi.com/2076-2607/12/6/1038)（相关综述群） | 综述：嗜盐古菌 PHA 合成机制；含 PHA 降解/mobilization 的讨论（PhaZh1/PhaJ 通路）。 |
| 9 | **Mitra R, et al. Haloarchaea as emerging big players in future polyhydroxyalkanoate bioproduction: Review of trends and perspectives.** Curr Res Biotechnol 2022. | doi:10.1016/j.crbiot.2022.09.002；[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0960308524001603)（相关） | 2022 年综述：嗜盐古菌 PHA 生物合成与 mobilization 现状。 |
| 10 | **Bioplastic Production from Agri-Food Waste through the Use of *Haloferax mediterranei*: A Comprehensive Initial Overview.** Microorganisms 2024, 12(6):1038. | doi:10.3390/microorganisms12061038；[MDPI](https://www.mdpi.com/2076-2607/12/6/1038) | 最新综述：*H. mediterranei* PHA 合成/调控（PhaR 调控子、颗粒相关蛋白）与降解通路总览。 |

### 2.2 重点问题回答：除 haloarchaea 外其他古菌门是否有 PHA 降解基因证据？

**结论：直接实验证据仍仅见于嗜盐古菌（*Haloferax mediterranei* 的 PhaZh1 + PhaJ1）**，但存在以下非嗜盐古菌的序列/旁证：

1. **泉古菌门（Crenarchaeota）**：*Sulfolobus acidocaldarius* DSM 639 热稳定脂解酶与 PHA 解聚酶高度相似（Arpigny & Jendrossek 1998，FEMS Microbiol Lett）——**最直接的序列同源证据**（尚未验证 PHA 水解活性，文献标题即说明"displaying similarity to polyhydroxyalkanoate depolymerases"）。
2. **基因组/生信证据**：Roh 2018（Sci Rep）古菌储能代谢通路分析显示 PHA 通路集中于 *Euryarchaeota*（尤其 *Halobacteria*）；**未见 Methanosarcina/Thermococcus/Pyrococcus 等产甲烷或超嗜热古菌的实验性 PHA 降解报道**——本任务 OpenAlex 检索 "archaeal polyhydroxyalkanoate"（46 条）与 "haloarchaea … degradation OR depolymerase"（50 条）中，无 Methanosarcina/Thermococcus 相关降解论文。
3. **宏基因组证据**：目前公开检索未发现专门报道"宏基因组古菌 contig 上的 phaZ 基因"的论文；现有宏基因组 PHA 研究聚焦**合成基因 phaC**（如西太平洋海洋样品宏基因组挖掘 [Discover Oceans 2026](https://link.springer.com/article/10.1007/s44289-026-00132-6?fromPaywallRec=false)）。Shark Bay 微生物席宏基因组功能复杂性研究（ISME J 2018, doi:10.1038/s41396-018-0208-8）覆盖 PHA 相关功能但非专门针对古菌降解基因。
4. **嗜热微生物宏分析**：*Tracking polyhydroxyalkanoate biosynthesis in thermophilic microorganisms*（Int J Biol Macromol 2025, doi:10.1016/j.ijbiomac.2025.148573；[bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.06.652502v1.full)、[PubMed 41176000](https://pubmed.ncbi.nlm.nih.gov/41176000/)）指出**嗜热古菌属中 PHA 合成基因检出率存在差异**——为后续在 GTDB 古菌基因组中扫 phaZ 提供背景参考。

> **对种子库构建的意义**：古菌 PHA 降解基因（phaZh1 类 patatin + phaJ 类烯脂酰-CoA 水合酶 + bdhA）目前只在 *Haloferax mediterranei* 有功能验证；**建议以 PhaZh1（HFX_6463，patatin/磷脂酶 A2 折叠）为古菌种子，而非胞外 α/β-水解酶型 phaZ**——两者折叠完全不同。

---

## 三、细菌各类型解聚酶代表酶清单（含 UniProt accession）

### 3.1 胞外 scl-PHA（PHB）解聚酶 — EC 3.1.1.75，α/β-水解酶，lipase box G-X-S-X-G

| 酶/基因 | 物种 | UniProt/PDB | 备注 |
|---|---|---|---|
| PhaZ（PhaZ1） | *Ralstonia pickettii* T1（= *Alcaligenes faecalis* T1） | **P12625**（BRENDA 收录；[NCBI protein P12625.1](https://www.ncbi.nlm.nih.gov/protein/P12625.1)） | 经典胞外 PHB 解聚酶，催化域+连接域+底物结合域（SBD） |
| PhaZ1 | *Pseudomonas lemoignei* | [AEM 1995, 177:596–607](https://doi.org/10.1128/jb.177.3.596-607.1995) | 酶 C（PHB 解聚酶 C） |
| PhaZ2 | *P. lemoignei* | 同上 | 酶 B |
| PhaZ3 | *P. lemoignei* | 同上 | 酶 D |
| PhaZ4 | *P. lemoignei* | [AEM 2000, 66(4):1385–1392](https://doi.org/10.1128/aem.66.4.1385-1392.2000) | PHV（3-羟基戊酸）解聚酶 |
| PhaZ5 | *P. lemoignei* | 同上 | 酶 A（PHB） |
| **PhaZ7** | *Paucimonas lemoignei* | PDB **2VTV**；ESTHER 家族 `PHAZ7_phb_depolymerase` | 耐热耐碱、水解无定形 PHB；[JMB 2008 结构](https://www.sciencedirect.com/science/article/abs/pii/S0022283608009479?via%3Dihub) |
| PhaZ | *Comamonas* sp.（JM 与 T1 株） | [JB/JPEP 1993](https://doi.org/10.1007/bf01457653) | 纯化的新型 PHB 解聚酶 |
| PhaZ | *Pseudomonas stutzeri* YM1414 | Ohura 1999 AEM（[ESTHER 页](https://bioweb.supagro.inrae.fr/ESTHER/paper/Ohura_1999_Appl.Environ.Microbiol_65_189)） | SBD 功能分析 |
| PhaZ | *Streptomyces exfoliatus* / *Streptomyces* spp. | [AMB 2012, 93:1975](https://bioweb.supagro.inrae.fr/ESTHER/paper/Garcia-Hidalgo_2012_Appl.Microbiol.Biotechnol_93_1975) | 异源表达于 *Rhodococcus* T104 |
| PhaZ | *Lihuaxuella thermophila* | PDB **8DAJ** | 耐热混杂 PHB 解聚酶（[RCSB](https://www.rcsb.org/structure/8DAJ)） |
| PhaZ | *Bacillus thuringiensis*（胞内型） | PDB **8YNW**（[RCSB](https://www.rcsb.org/structure/8YNW)）；[JB 2006, 188:8883](https://journals.asm.org/doi/full/10.1128/jb.00729-06) | 胞内 PHB 解聚酶新类型 |
| PhaZ | *Penicillium funiculosum*（真菌） | PDB **2D80** | 真菌胞外 PHB 解聚酶 |

### 3.2 胞内/颗粒结合解聚酶（scl 与 mcl）

| 酶/基因 | 物种 | 备注 |
|---|---|---|
| PhaZ（**i-1** / PhaZ1） | *Cupriavidus necator*（*Ralstonia eutropha*）H16 | 胞内 PHB 解聚酶；[AMB Express 2012, 2:26](https://link.springer.com/article/10.1186/2191-0855-2-26)；[PMC3430594](https://pmc.ncbi.nlm.nih.gov/articles/PMC3430594/) |
| PhaZ2（**i-2**） | *C. necator* H16 | 同上；H16 编码 2–3 个胞内解聚酶基因（[Ralstonia eutropha H16 … Depolymerase Genes](https://nufind.nu.edu.sa/EdsRecord/edsair,edsair.doi.dedup.....b811ecf11c60626c11ae20ed11ed5a81)） |
| PhaZd1（**PhaZ6**）/ PhaZd2（**PhaZ7**） | *C. necator* H16 | 高活性 PHB 解聚酶但对 PHB mobilization 无可见作用（[BRENDA 3.1.1.75](https://www.brenda-enzymes.org/result_download.php?a=30&RN=&RNV=&os=&pt=&FNV=&tt=&SYN=&Textmining=&W[1]=3.1.1.75&T[1]=1&T[9]=2&T[12]=1&orderDesc=1&orderByHTMLField=Title&nolimit=1)；Sznajder & Jendrossek） |
| PhaZ | *Rhodospirillum rubrum* | **周质定位**，特异水解天然 PHB，催化域为 type II（[JB 2004, 186:7243](https://doi.org/10.1128/jb.186.21.7243-7253.2004)；[AMB 2011, 91:971](https://doi.org/10.1007/s00253-011-3096-7)） |
| **PhaZ** | *Pseudomonas putida* KT2442/KT2440 | **胞内 mcl-PHA 解聚酶**（α/β-水解酶 + lid 结构）；[JBC 2006](https://doi.org/10.1074/jbc.m608119200)；[AMB 2025 lid 功能](https://link.springer.com/article/10.1007/s00253-025-13605-z?fromPaywallRec=true) |
| PhaZc（3-HB 寡聚体水解酶）+ Hbd（3-HB 脱氢酶） | *Paracoccus denitrificans* | 胞内 PHB/PHV 降解两步酶（[AEM 2013](https://doi.org/10.1128/aem.03396-13)；FEMS 2001 [P. denitrificans phaZ](https://doi.org/10.1111/j.1574-6968.2001.tb10558.x)） |
| mcl 解聚酶 Bd3285 等 | *Bdellovibrio bacteriovorus* | 捕食型水解武器库中的 mcl-PHA 解聚酶（[AEM 2012](https://doi.org/10.1128/aem.01099-12)；[Sci Rep 2016](https://doi.org/10.1038/srep24381)） |
| PhaZ | *Bacillus thuringiensis* | 新类型胞内 PHB 解聚酶（[JB 2006](https://journals.asm.org/doi/full/10.1128/jb.00729-06)） |
| PhaZ | *Paracoccus denitrificans*（胞内 phaZ，靠近 phaC） | [FEMS 2001](https://doi.org/10.1111/j.1574-6968.2001.tb10558.x) |

### 3.3 古菌 patatin 类解聚酶

| 酶/基因 | 物种 | 备注 |
|---|---|---|
| **PhaZh1**（HFX_6463） | *Haloferax mediterranei* | **patatin 样**（磷脂酶 A2/patatin 折叠），颗粒结合，lipase box G-X-S47-X-G；与 bdhA 成簇（[AEM 2015](https://doi.org/10.1128/aem.04269-14)） |
| **PhaJ1**（HFX_？） | *H. mediterranei* | (R)-烯脂酰-CoA 水合酶，颗粒结合，连接 β-氧化（[Sci Rep 2016](https://doi.org/10.1038/srep24015)） |

### 3.4 胞外解聚酶分类（类型 I/II 与 SBD 类型）— 用于 HMM 家族设计

- **催化域（CD）两型**（[Frontiers 2025 综述](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1542468/full)、[europepmc PMC11893044](https://europepmc.org/article/PMC/11893044)）：
  - **Type I**：lipase box 位于序列中部
  - **Type II**：lipase box 位于序列 N 端附近
- **结构域架构**：信号肽 + 催化域（CD）+ 连接域（linker）+ 底物结合域（SBD）。SBD 类型多样：疏水型、fibronectin type III 样、cadherin 样、苏氨酸富集型等（影响底物结合与结晶态 PHB 水解，[Ohura 1999](https://bioweb.supagro.inrae.fr/ESTHER/paper/Ohura_1999_Appl.Environ.Microbiol_65_189)、[Jendrossek 2002 综述](https://pubmed.ncbi.nlm.nih.gov/12213937/)）。
- 催化三联体 **Ser-His-Asp** 为所有家族共有（[Make It or Break It 综述](https://www.preprints.org/manuscript/202411.1773)）。

---

## 四、方法建议：基于数据库家族体系设计 HMM 家族划分方案

### 4.1 家族划分框架（三层）

```
Level 0: 折叠层   α/β-水解酶折叠（全部胞外+多数胞内 phaZ）
                   vs. patatin/磷脂酶折叠（古菌 PhaZh1 类）
Level 1: 超家族层（参照 PhaDED 8 超家族）
                   S1 胞外 scl-PHA 解聚酶（Esterase_phb_PHAZ 类）
                   S2 胞外 mcl-PHA 解聚酶
                   S3 胞内 scl（颗粒结合，如 C. necator i-1/i-2、B. thuringiensis PhaZ）
                   S4 胞内 mcl（P. putida PhaZ 类）
                   S5 3HB 寡聚体水解酶（PhaZc 类，EC 3.1.1.22 相关）
                   S6 真菌 PHB 解聚酶（P. funiculosum 类）
                   S7 古菌 patatin 解聚酶（PhaZh1 类）
                   S8 古菌烯脂酰-CoA 水合酶/其他辅助（PhaJ 类，非水解酶但参与 mobilization）
Level 2: 家族层（参照 PhaDED 38 家族 + TIGRFAM/ESTHER）
```

### 4.2 具体实施建议

1. **种子收集**：
   - 从 **Pfam PF10503（Esterase_PHB）** 与 **TIGRFAM TIGR01840（esterase_phb）、TIGR01849（PHB_depoly_PhaZ）** 直接下载 HMM/种子比对，作为细菌胞外/胞内 phaZ 主干。
   - 用 **PhaDED 587 条序列（8 超家族 38 家族）** 的家族 MSA 作为内部参考（若 DED 站点不可达，用 ESTHER `Esterase_phb_PHAZ`、`PHAZ7_phb_depolymerase` 家族替代）。
   - 古菌种子：PhaZh1（AEM 2015，UniProt 条目可经 HFX_6463 定位）+ PhaJ1（Sci Rep 2016）+ bdhA。
   - 每条种子带实验证据标签（文献 DOI + UniProt/PDB 编号），构建种子-文献映射表。
2. **HMM 构建与划分**：
   - 每超家族单独建 HMM（避免跨折叠污染）；用 `hmmbuild`（HMMER3）+ 序列权重。
   - 先按 PhaDED 超家族划分，再在超家族内用系统发育（IQ-TREE/FastTree）二次切分家族；以"序列一致性 + 催化残基保守 + 结构域架构（CD 类型 I/II、SBD 类型）"为家族判据。
   - 对古菌基因组（GTDB 分类）单独用 patatin HMM 扫描，避免与细菌 α/β-水解酶 HMM 混淆。
3. **阈值与验证**：
   - 用 curated 正集（UniProt reviewed EC 3.1.1.75/76 + PDB 结构序列）与负集（脂肪酶、酯酶、角质酶）定 gathering threshold。
   - 交叉验证：PAZy 已验证酶列表、BRENDA EC 3.1.1.75/76 条目、PDB 结构。
4. **命名规则**（建议沿用社区惯例）：
   - 细菌胞外：`PhaZ_ext_scl` / `PhaZ_ext_mcl`（或直接引用 PhaDED 家族名）
   - 细菌胞内：`PhaZ_int_scl` / `PhaZ_int_mcl`
   - 寡聚体水解酶：`3HB_oligomer_hydrolase`（EC 3.1.1.22）
   - 古菌：`PhaZh1_patatin` / `PhaJ_enoyl_CoA_hydratase` / `BdhA_3HB_dehydrogenase`
5. **已知坑**：
   - 胞外与胞内解聚酶序列相似度低、系统发育混杂 → 必须用 HMM 而非 BLAST 单阈值。
   - `PhaZ6/PhaZ7`（*C. necator*）虽高活性但对体内 mobilization 无作用 → 功能标注需区分"体外活性"与"体内功能"。
   - patatin 类（古菌 PhaZh1）与细菌 α/β-水解酶 phaZ **折叠不同、无序列相似性** → 不可共用 HMM。

---

## 五、所有来源 URL 列表

### OpenAlex 文献（核心）
- https://pubmed.ncbi.nlm.nih.gov/19296857/ （Knoll 2009, PhaDED；doi:10.1186/1471-2105-10-89）
- https://pubmed.ncbi.nlm.nih.gov/25710370/ （Liu 2015, PhaZh1 patatin；doi:10.1128/aem.04269-14）
- https://www.nature.com/articles/srep24015 （Han 2016, PhaJ1；doi:10.1038/srep24015）
- https://journals.asm.org/doi/10.1128/aem.03396-13 （Paracoccus PhaZc/Hbd；doi:10.1128/aem.03396-13）
- https://doi.org/10.1128/jb.186.21.7243-7253.2004 （R. rubrum 周质 PhaZ1）
- https://doi.org/10.1074/jbc.m608119200 （P. putida KT2442 胞内 mcl PhaZ）
- https://doi.org/10.6026/97320630015036 （PHA genes Database 2019）
- https://doi.org/10.1007/s00253-010-2611-6 （Halococcus 等古菌 PHA 鉴定）
- https://doi.org/10.1128/jb.177.3.596-607.1995 （P. lemoignei 5 基因解聚酶系统）
- https://doi.org/10.17516/1997-1389-0024 （PHA 降解综述 2017）
- https://doi.org/10.1128/aem.01117-10 （古菌 IIIA 型 PHA 合酶）
- https://doi.org/10.1111/j.1574-6968.2001.tb10558.x （P. denitrificans 胞内 phaZ）
- https://doi.org/10.1016/j.jbiotec.2016.04.004 （C. necator H16 解聚酶与分子量）
- https://doi.org/10.1007/s00253-011-3096-7 （R. rubrum 新型胞内 PHB 解聚酶）
- https://doi.org/10.1146/annurev.micro.56.012302.160838 （Jendrossek 2002 综述）
- https://doi.org/10.1128/aem.01099-12 （Bdellovibrio mcl 解聚酶）
- https://doi.org/10.1128/aem.66.4.1385-1392.2000 （P. lemoignei PHV 解聚酶）
- https://doi.org/10.1038/s41598-018-37768-0 （古菌储能代谢生信分析）
- https://doi.org/10.1128/aem.03791-14 （R. eutropha 颗粒蛋白组，PhaZ6/7 等）
- https://doi.org/10.1128/jb.00729-06 （B. thuringiensis 胞内 phaZ）
- https://doi.org/10.1038/s41396-018-0208-8 （Shark Bay 宏基因组功能复杂性）
- https://doi.org/10.1007/s44289-026-00132-6 （西太平洋宏基因组 PHA 基因挖掘）
- https://doi.org/10.1016/j.ijbiomac.2025.148573 （嗜热微生物 PHA 生物合成追踪）

### 数据库/家族
- https://bioweb.supagro.inrae.fr/ESTHER/ （ESTHER 主页）
- https://bioweb.supagro.inrae.fr/ESTHER/family/Esterase_phb_PHAZ （家族页）
- https://bioweb.supagro.inrae.fr/ESTHER/family/PHAZ7_phb_depolymerase （家族页）
- https://bioweb.supagro.inrae.fr/ESTHER/paper/Knoll_2009_BMC.Bioinformatics_10_89
- https://bioweb.supagro.inrae.fr/ESTHER/paper/Buchholz_2022_Proteins_90_1443
- https://pazy.eu/ 、https://www.cbl.uni-stuttgart.de/doku.php?id=start （PAZy）
- https://www.pazy.eu/plastics/pa?substrate=4&limit=50&offset=50 （PAZy 按底物）
- https://onlinelibrary.wiley.com/doi/full/10.1002/prot.26325 （PAZy 论文）
- https://www.brenda-enzymes.org/enzyme.php?ecno=3.1.1.75 （EC 3.1.1.75）
- https://brenda-enzymes.org/enzyme.php?ecno=3.1.1.76 （EC 3.1.1.76）
- https://www.brenda-enzymes.de/search_result.php?T%5b1%5d=1&T%5b9%5d=2&W%5b1%5d=3.1.1.22 （EC 3.1.1.22）
- https://www.uniprot.org/uniprotkb?query=(ec:3.1.1.75)%20AND%20(reviewed:true) （UniProt reviewed）
- https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam10503 （CDD Esterase_PHB）
- https://www.ncbi.nlm.nih.gov/Structure/cdd/TIGR01840 （TIGR01840 esterase_phb）
- http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=TIGR01840
- https://www.ncbi.nlm.nih.gov/Structure/cdd/TIGR01849 （TIGR01849 PHB_depoly_PhaZ）
- http://tigrfams.jcvi.org/cgi-bin/HmmReportPage.cgi?acc=TIGR01849
- https://supfam.org/scop/159747 （SCOP PHB depolymerase-like）
- http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=2D80 （PDB 2D80）
- https://www.rcsb.org/structure/8DAJ （PDB 8DAJ）
- https://www.rcsb.org/structure/8YNW （PDB 8YNW）
- https://pdbj.org/mine/structural_details/2vtv （PDB 2VTV PhaZ7）
- https://en.wikipedia.org/?curid=14457442 （PHB depolymerase 维基）

### web 检索佐证（综述/新闻）
- https://europepmc.org/article/PMC/11893044 （PHA 生物降解现状与展望）
- https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1542468/full （Frontiers 2025 综述）
- https://www.preprints.org/manuscript/202411.1773 （Pha 合酶与解聚酶综述）
- https://cyberleninka.ru/article/n/degradation-of-polyhydroxyalkanoate-pha-a-review （PHA 降解综述）
- https://www.sciencedirect.com/science/article/abs/pii/S0301479725009016 （PHB 解聚酶与废物管理 2025）
- https://www.mdpi.com/2076-2607/12/6/1038 （H. mediterranei 综述 2024）
- https://www.biorxiv.org/content/10.1101/2025.05.06.652502v1.full （嗜热微生物 PHA bioRxiv）
- https://www.science.gov/topicpages/a/archaeon+haloferax+mediterranei.html （Science.gov 汇总）
- https://link.springer.com/article/10.1186/s12934-020-01342-z （嗜盐菌 PHA 合成综述）
- https://pmc.ncbi.nlm.nih.gov/articles/PMC3227634/ （Halomonas TD01 PHA 基因比较基因组）

---

*OpenAlex 检索结果 JSON 存于 `D:\PHB_gtdb-ds\research\openalex\comp_search1..7_*.json` 与 `comp_work_*.json`（14 篇关键论文全文元数据+摘要）。所有检索均遵循 OpenAlex skill 规则（resolve→filter、--select 精简、不伪造 ID/DOI）。*
