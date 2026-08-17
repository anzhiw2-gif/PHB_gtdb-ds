# PHB/PHA 降解基因文献调研报告（OpenAlex 检索）

项目背景：基于 GTDB 数据库的 PHB 降解基因系统生信分析
检索日期：2026-06（OpenAlex API，免费 polite pool）
中间数据：本目录下 q*.json 及 key_papers_details.json（均来自 OpenAlex 原始返回，未做任何虚构）

## 一、查询命中数汇总

| # | 检索式 | 方式 | 命中数 |
|---|--------|------|--------|
| Q1 | polyhydroxybutyrate depolymerase | fulltext.search | 1279 |
| Q2 | polyhydroxyalkanoate depolymerase | fulltext.search | 1948 |
| Q3a | PHA degradation genes genome | fulltext.search | 13033（噪声大，弃用） |
| Q3b | depolymerase AND polyhydroxyalkanoate（标题/摘要） | title_and_abstract.search | 262 |
| Q4 | polyhydroxybutyrate degradation（标题/摘要） | title_and_abstract.search | 876（含大量 PLA/组织工程噪声） |
| Q5 | polyhydroxyalkanoate AND depolymerase AND database | title_and_abstract.search | 6 |
| Q6a | "Genome Taxonomy Database"（标题/摘要） | title_and_abstract.search | 6666（偏宽） |
| Q6b | GTDB（标题/摘要） | title_and_abstract.search | 3506 |
| Q6c | GTDB AND mining | title_and_abstract.search | 20 |
| Q7 | plastic degrading enzymes AND genome mining | title_and_abstract.search | 35 |
| Q8 | PHB depolymerase（标题/摘要） | title_and_abstract.search | 499 |
| Q9 | standardized bacterial taxonomy | title_and_abstract.search | 494 |
| Q10 | crystal structure AND polyhydroxybutyrate depolymerase | title_and_abstract.search | 6 |
| Q11 | Ideonella sakaiensis（标题/摘要） | title_and_abstract.search | 409 |
| Q12 | PET hydrolase AND genome mining | title_and_abstract.search | 21 |
| DOI 批量补全 | 28 个 DOI | filter doi | 28/28 |

## 二、关键论文清单（24 篇）

### ① 酶学与机制经典文献（9 篇）
1. Dawes & Senior, 1973. The Role and Regulation of Energy Reserve Polymers in Micro-organisms. *Advances in Microbial Physiology*. 864 次引用. doi:10.1016/s0065-2911(08)60088-0 — PHB 作为微生物储能聚合物的奠基性综述，代谢（合成/降解）调控的经典框架。
2. Jendrossek & Handrick, 2002. Microbial Degradation of Polyhydroxyalkanoates. *Annual Review of Microbiology*. 692 次引用. doi:10.1146/annurev.micro.56.012302.160838 — PHA 降解领域最经典综述：胞外/胞内解聚酶分类、phaZ 基因家族、降解机制与生态意义。
3. (2004) Review Degradation of microbial polyesters. *Biotechnology Letters*. 339 次引用. doi:10.1023/b:bile.0000036599.15302.e5 — 微生物聚酯（PHB/PHA）降解机理综述。
4. (1992) Enzymatic degradation of microbial poly(3-hydroxybutyrate) films. *Die Makromolekulare Chemie*. 232 次引用. doi:10.1002/macp.1992.021930105 — 早期 PHB 薄膜酶解实验经典。
5. Grage et al., 2005. Poly(3-hydroxybutyrate) Granule-Associated Proteins: Impacts on PHB Synthesis and Degradation. *Biomacromolecules*. 226 次引用. doi:10.1021/bm049401n — 颗粒结合蛋白（phaP/phaZ 等）对 PHB 合成与降解的调控。
6. (2005) The Crystal Structure of Polyhydroxybutyrate Depolymerase from Penicillium funiculosum. *Journal of Molecular Biology*. 122 次引用. doi:10.1016/j.jmb.2005.12.028 — 首个 PHB 解聚酶晶体结构，揭示底物识别/催化残基。
7. (1997) Biochemical and molecular characterization of the polyhydroxybutyrate depolymerase of Comamonas acidovorans YM1609. *Applied and Environmental Microbiology*. 83 次引用. doi:10.1128/aem.63.12.4844-4852.1997 — 胞外 PHB 解聚酶基因克隆与生化表征经典。
8. (2011) Biochemical characterization of a new type of intracellular PHB depolymerase from Rhodospirillum rubrum. *Applied Microbiology and Biotechnology*. 27 次引用. doi:10.1007/s00253-011-3096-7 — 新型胞内（颗粒结合型）PHB 解聚酶表征。
9. (2022) Bioplastic degradation by a polyhydroxybutyrate depolymerase from a thermophilic soil bacterium. *Protein Science*. 26 次引用. doi:10.1002/pro.4470 — 嗜热菌来源 PHB 解聚酶新酶发现。

### ② 基因组/宏基因组筛选（8 篇）
10. (2006) Genome sequence of the bioplastic-producing "Knallgas" bacterium Ralstonia eutropha H16. *Nature Biotechnology*. 669 次引用. doi:10.1038/nbt1244 — 模式产 PHB 菌株全基因组，含完整 PHA 代谢基因集（参考基因组）。
11. (2018) New Insights into the Function and Global Distribution of Polyethylene Terephthalate (PET)-Degrading Bacteria and Enzymes in Marine and Terrestrial Metagenomes. *Applied and Environmental Microbiology*. 445 次引用. doi:10.1128/aem.02773-17 — 宏基因组大规模筛选降解酶（PETase 同源物全球分布）的方法学范例。
12. Poblete-Castro et al., 2014. A holistic view of polyhydroxyalkanoate metabolism in Pseudomonas putida. *Environmental Microbiology*. 232 次引用. doi:10.1111/1462-2920.12760 — P. putida 全基因组水平的 PHA 合成-降解代谢网络综述。
13. (2019) Shotgun Metagenomics Reveals the Benthic Microbial Community Response to Plastic and Bioplastic in a Coastal Marine Environment. *Frontiers in Microbiology*. 215 次引用. doi:10.3389/fmicb.2019.01252 — 宏基因组解析生物塑料降解微生物群落。
14. Matsusaki et al., 1998. Cloning and Molecular Analysis of the Poly(3-hydroxybutyrate) and Poly(3-hydroxybutyrate-co-3-hydroxyalkanoate) Biosynthesis Genes in Pseudomonas sp. Strain 61-3. *Journal of Bacteriology*. 214 次引用. doi:10.1128/jb.180.24.6459-6467.1998 — pha 基因簇（含 phaZ 降解基因）克隆分析经典。
15. Yoshida et al., 2016. A bacterium that degrades and assimilates poly(ethylene terephthalate). *Science*. 3342 次引用. doi:10.1126/science.aad6359 — 塑料降解酶发现里程碑（Ideonella sakaiensis/PETase）。
16. (2018) Polyhydroxybutyrate (PHB) biodegradation using bacterial strains with demonstrated and predicted PHB depolymerase activity. *Applied Microbiology and Biotechnology*. 78 次引用. doi:10.1007/s00253-018-9153-8 — 基于"预测解聚酶活性"筛选降解菌株，与本项目"基因预测→功能验证"思路最接近。
17. (2022) Identification of BgP, a Cutinase-Like Polyesterase From a Deep-Sea Sponge-Derived Actinobacterium. *Frontiers in Microbiology*. 23 次引用. doi:10.3389/fmicb.2022.888343 — 基因组挖掘（放线菌）发现新型聚酯酶。

### ③ 数据库与工具（3 篇）
18. Knoll et al., 2009. The PHA Depolymerase Engineering Database: A systematic analysis tool for the diverse family of polyhydroxyalkanoate (PHA) depolymerases. *BMC Bioinformatics*. 162 次引用. doi:10.1186/1471-2105-10-89 — **PHA 解聚酶专用数据库**（PhaD 分类系统），本项目最直接相关的前人数据库。
19. (2022) Bioplastics for a circular economy. *Nature Reviews Materials*. 2047 次引用. doi:10.1038/s41578-021-00407-8 — 生物塑料循环经济综述（背景）。
20. (2020) Microbial and Enzymatic Degradation of Synthetic Plastics. *Frontiers in Microbiology*. 1101 次引用. doi:10.3389/fmicb.2020.580709 — 塑料降解酶全景综述（酶资源盘点）。

### ④ GTDB 相关方法学范例（4 篇）
21. Parks et al., 2018. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nature Biotechnology*. 3964 次引用. doi:10.1038/nbt.4229 — **GTDB 原始论文**：基于基因组系统发育的标准化细菌分类框架。
22. Chaumeil et al., 2019. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics*. 5648 次引用. doi:10.1093/bioinformatics/btz848 — GTDB 分类工具（基因组→GTDB 分类）。
23. Parks et al., 2021. GTDB: an ongoing census of bacterial and archaeal diversity. *Nucleic Acids Research*. 2571 次引用. doi:10.1093/nar/gkab776 — GTDB 数据库持续更新（分类学普查）。
24. Chaumeil et al., 2022. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics*. 2018 次引用. doi:10.1093/bioinformatics/btac672 — GTDB-Tk v2（新版分类流程）。

## 三、要点与启示
- **GTDB 直接用于酶/基因挖掘的已发表范例极少**（Q6c 仅 20 条且多为低引用），说明"基于 GTDB 做 PHB 降解基因系统挖掘"有方法学创新空间；常见做法是 GTDB 分类 + 蛋白家族（Pfam/CAZy 类）注释 + 系统发育（如 Q12/Q14 的 PET 酶挖掘）。
- 与 PHB 降解最直接的核心资源：Jendrossek 综述（②）、PHA Depolymerase Engineering Database（⑬）与"预测解聚酶活性"筛选论文（⑯）。
- 注意：Q4 检索噪声大（"polyhydroxybutyrate degradation" 命中多为 PLA/生物医学材料），PHB 降解相关高相关文献需从 Q3b/Q8/Q10 结果中选取。

## 四、论文 URL 列表（OpenAlex）
- https://openalex.org/W1815987069 (Dawes 1973)
- https://openalex.org/W2164479864 (Jendrossek 2002)
- https://openalex.org/W2023631098 (2004 review)
- https://openalex.org/W1899396383 (1992 enzymatic degradation)
- https://openalex.org/W2087453257 (2005 granule proteins)
- https://openalex.org/W2056256968 (2005 crystal structure)
- https://openalex.org/W2120707072 (1997 Comamonas)
- https://openalex.org/W2026010604 (2011 R. rubrum)
- https://openalex.org/W4304689113 (2022 thermophilic)
- https://openalex.org/W2125645956 (2006 Ralstonia genome)
- https://openalex.org/W2791560780 (2018 PET metagenomes)
- https://openalex.org/W1549936434 (2014 P. putida)
- https://openalex.org/W2955133522 (2019 shotgun metagenomics)
- https://openalex.org/W1592377749 (1998 Pseudomonas 61-3)
- https://openalex.org/W2294707565 (2016 Yoshida Science)
- https://openalex.org/W2809748194 (2018 predicted depolymerase)
- https://openalex.org/W4223437726 (2022 BgP)
- https://openalex.org/W2131164625 (2009 PHA Depolymerase DB)
- https://openalex.org/W4205722716 (2022 bioplastics circular economy)
- https://openalex.org/W3107539693 (2020 synthetic plastics degradation)
- https://openalex.org/W2889019390 (2018 Parks GTDB)
- https://openalex.org/W2986925300 (2019 GTDB-Tk)
- https://openalex.org/W3200103613 (2021 GTDB census)
- https://openalex.org/W4304481015 (2022 GTDB-Tk v2)
