# PHB（聚羟基丁酸酯）降解基因 PubMed 文献调研汇总

> 检索工具：NCBI E-utilities（经 `pubmed-database` 技能包装脚本 `scripts/pubmed_api.py`）
> 检索日期：2026-08-15 ｜ 输出目录：`research/pubmed/`
> 中间文件：`search_0*.json`（原始 PMID 列表）、`*_slim.json`（精简元数据）、`abstracts_all.json`（154 篇摘要）、`selected_top.json`（重点论文摘要）、`count_*.json`（各库命中数）

## 一、检索概况与命中数

| # | 查询式 | PubMed 命中数 | 返回条数 |
|---|--------|--------------|---------|
| 1 | `polyhydroxybutyrate depolymerase` | 191 | 25 |
| 2 | `polyhydroxyalkanoate depolymerase` | 83 | 25 |
| 3 | `PHB depolymerase` | 194 | 20 |
| 4 | `phaZ AND depolymerase` | 87 | 20 |
| 5 | `poly(3-hydroxybutyrate) AND (degradation OR depolymerization)` | 2,415 | 25 |
| 6 | `polyhydroxyalkanoate AND (genome OR metagenome) AND (degradation OR depolymerase)` | 855 | 25 |
| 7 | `PHA degradation AND bacteria` | 3,304 | 25 |
| 8 | `polyhydroxyalkanoate degradation review` | 535 | 15 |

- 8 个查询共返回 180 条 PMID，去重后 **154 篇唯一论文**，全部批量抓取元数据+摘要。
- 命中数通过 `global_database_discovery` 实时统计；其余数据库（pmc/gene/protein 等）命中见 `count_*.json`。
- 注：查询 5/7 命中数高但含较多非降解酶主题（材料/生产/应用），按相关度排序取回 top 结果并人工筛选。

## 二、重点论文清单（26 篇，含摘要分析）

### ① 酶学与机制（胞内/胞外 PHB 解聚酶、PhaZ 分类）

| PMID | 年份 | 期刊 | 第一作者 | 中文标题 | 相关点 |
|------|------|------|---------|---------|--------|
| 16405909 | 2006 | J Mol Biol | Hisano T | 真菌 *Penicillium funiculosum* PHB 解聚酶晶体结构 | 1.71Å 结构：α/β水解酶折叠（环状置换变体），催化三联体 Ser-Asp-His + 氧阴离子洞，阐明底物识别与降解机制 |
| 9406404 | 1997 | Appl Environ Microbiol | Kasuya K | 淡水菌 *Comamonas acidovorans* YM1609 PHB 解聚酶的生化与分子表征 | 胞外 PHB 解聚酶纯化（45 kDa），丝氨酸酯酶（DFP 抑制），降解 P(3HB) 及共聚物 |
| 24907326 | 2014 | Appl Environ Microbiol | Sznajder A | *Ralstonia eutropha* PhaZd1/PhaZd2（PhaZ6/7）到底是不是 PHB 解聚酶 | 胞内解聚酶 PhaZ6/7 体外高活性，但对积累 PHB 动员无可见作用——胞内解聚酶功能冗余问题的经典辨析 |
| 29678915 | 2018 | Appl Environ Microbiol | Juengert JR | *R. eutropha* PHB 聚合酶 PhaC1 与解聚酶 PhaZa1 的磷酸化 | 揭示 PhaC1 Thr373 与 PhaZa1 Ser35 的磷酸化修饰，PHB 积累/降解期翻译后调控 |
| 16936025 | 2006 | J Bacteriol | Tseng CL | *Bacillus thuringiensis* 新型胞内 PHB 解聚酶 phaZ 基因 | 注释为 PcaD 的基因实为新型胞内 PhaZ，无已知解聚酶同源，降解天然/无定形 PHB 颗粒 |
| 20346169 | 2010 | BMC Microbiol | Trainer MA | *Sinorhizobium meliloti* 胞内 PHB 解聚酶 PhaZ 的鉴定 | 共生固氮菌 PHB 循环降解支路唯一未表征酶的克隆与功能验证（原位/异源） |
| 12898135 | 2003 | Arch Microbiol | Kadouri D | *Azospirillum brasilense* PhaZ 及 phaZ 突变体 | phaZ 敲除突变体失去 PHB 降解能力，证明该酶在 PHB 动员中的必需性 |
| 17170116 | 2007 | J Biol Chem | de Eugenio LI | *P. putida* KT2442 phaZ 编码特异性胞内 mcl-PHA 解聚酶的生化证据 | mcl-PHA 胞内解聚酶的"范式性"表征（scl 解聚酶研究多、mcl 少） |
| 19788655 | 2010 | Environ Microbiol | de Eugenio LI | *P. putida* KT2442 mcl-PHA 周转与 PhaZ 的核心作用 | PhaZ 在饥饿期 mcl-PHA 降解与代谢平衡中的关键角色 |
| 19214501 | 2009 | Appl Microbiol Biotechnol | Papaneophytou CP | *Thermus thermophilus* HB8 胞外 PHB 解聚酶 | 嗜热菌胞外解聚酶纯化（42 kDa），基因 TTHA0199，与羧酸酯酶同源 |
| 15489436 | 2004 | J Bacteriol | Handrick R | *Rhodospirillum rubrum* "胞内"PHB 解聚酶实为周质定位蛋白 | 新类型：周质定位、特异降解天然无定形 PHB、结构类似胞外解聚酶——挑战"胞内/胞外"二分法 |
| 23951224 | 2013 | PLoS One | García-Hidalgo J | *Streptomyces ascomycinicus* 新型胞外 PHB 解聚酶 | fkbU 基因（PhaZ_Sa），48.4 kDa，酸性条件下降解 PHB 共聚物 |
| 41151231 | 2026 | Enzyme Microb Technol | Anjulal H | *Nocardiopsis dassonvillei* PHB 解聚酶过表达、生化表征与结构建模 | 密码子优化+重组表达（~50 kDa），分子对接与结构模型 |
| 36222314 | 2022 | Protein Sci | Thomas GM | 嗜热土壤菌 *Lihuaxuella thermophila* PHB 解聚酶降解生物塑料 | 碱性 PHB 解聚酶（LtPHBase）结构与生化：广谱活性（均/共聚物），工业降解潜力 |
| 26156240 | 2015 | Appl Microbiol Biotechnol | Martínez V | *Streptomyces exfoliatus* K10 新型胞外 mcl-PHA 解聚酶 | PhaZSex2：脂酶框（Ser-His-Asp）、α/β水解酶折叠，降解天然及功能化 mcl-PHA |
| 24751505 | 2014 | Int J Biol Macromol | Martínez V | 细胞系统工程生产胞外 PHA 解聚酶 | 重组 *P. fluorescens* GK13 PhaZGK13 生产平台，用于 PHA 降解与 3-羟基脂肪酸回收 |
| 38969063 | 2024 | J Biol Chem | Santolin L | *R. eutropha* PHA 代谢调控：phasin 与解聚酶基因的转录调控子 | 启动子 pull-down+质谱鉴定调控 PHA 稳态的转录因子网络 |
| 34342882 | 2021 | J Basic Microbiol | Nygaard D | *Cupriavidus necator* 不同营养条件下的 PHA 颗粒形成与降解 | 显微+分析手段整合表征颗粒形态/含量/产率（模式菌） |

### ② 基因组/宏基因组筛选研究

| PMID | 年份 | 期刊 | 第一作者 | 中文标题 | 相关点 |
|------|------|------|---------|---------|--------|
| 29951858 | 2018 | Appl Microbiol Biotechnol | Martínez-Tobón DI | 用已验证与**预测** PHB 解聚酶活性菌株进行 PHB 生物降解 | 9 株菌（5 实测+4 预测 PhaZ）对比降解 PHB 膜；预测活性来源（含基因组推断）——本项目方法学直接参考 |
| 32369496 | 2020 | PLoS One | Morohoshi T | 基于 *Undibacterium* sp. KW1/YM2 **完整基因组**鉴定新型胞外 PHA 解聚酶 | 全基因组测序定位 phaZUD 基因，功能验证+系统发育定位 |
| 40500476 | 2025 | Mar Biotechnol | Iseki K | *Alteromonas* 属胞外 P(3HB) 解聚酶鉴定及其**系统发育分布** | 基因组分析找 PhaZ 同源→催化域（lipase box）+底物结合域，跨种系统发育分布分析——与 GTDB 系统发育思路高度契合 |
| 36125959 | 2022 | Microb Genom | Leadbeater DR | 细菌海藻降解生物塑料生产者的 in silico 鉴定 | 挖掘基因/蛋白数据库（数据库挖掘方法学），鉴定利用海藻糖生产 PHA 的菌 |
| 38278791 | 2024 | Nat Commun | Omura T | 深海海底生物可降解塑料的微生物分解 | 757–5552 m 深海部署实验+微生物分析，PHA 类在深海可降解（宏基因组/环境微生物组关联） |
| 34160268 | 2021 | Appl Environ Microbiol | Eronen-Rasimus E | 海冰细菌 *Halomonas* sp. 363 与 *Paracoccus* sp. 392 低温产多种 PHA | 基因组+转录组解析低温 PHA 合成（生产方向，作菌株/基因组背景） |
| 40320445 | 2025 | Sci Rep | Boondaeng A | 海洋 PHA 降解菌分离与酶生产优化 | 海洋垃圾中筛得 6 株产胞外 PHA 解聚酶菌株（clear zone 法） |
| 37311705 | 2023 | J Microbiol Biotechnol | Jeon Y | *Bacillus infantis* 的 P(3HB) 降解 | 双层平板筛菌 + phaZ/bdhA 通用引物 PCR 检测——基因检测筛选流程参考 |
| 40392676 | 2025 | J Appl Microbiol | Hachisuka SI | 土壤细菌降解人工 PHA 共聚物 P(2HB-co-3HB) | 单体序列（random/block）对可降解性的影响，人工共聚物降解菌筛选 |

### ③ 数据库与生物信息学工具/概念

| PMID | 年份 | 期刊 | 第一作者 | 中文标题 | 相关点 |
|------|------|------|---------|---------|--------|
| 26409775 | 2015 | Trends Biotechnol | Chen GQ | "PHAome" | 概念框架：类比 genome/transcriptome，提出 PHA 谱系（单体/均聚/共聚/分子量）的动态全貌——为 PHA 相关组学/数据库构建提供概念基础 |
| 36125959 | 2022 | Microb Genom | Leadbeater DR | in silico 数据库挖掘（见②） | 基因/蛋白数据库挖掘方法，可迁移到 GTDB 基因组扫描流程 |
| 40500476 | 2025 | Mar Biotechnol | Iseki K | 系统发育分布分析（见②） | PhaZ 跨物种系统发育分布+基因组定位方法参考 |
| 38278791 | 2024 | Nat Commun | Omura T | 深海宏基因组（见②） | 深海环境宏基因组+降解实验结合 |

### ④ 综述

| PMID | 年份 | 期刊 | 第一作者 | 中文标题 | 相关点 |
|------|------|------|---------|---------|--------|
| 12213937 | 2002 | Annu Rev Microbiol | Jendrossek D | 聚羟基烷酸酯的微生物降解 | **奠基性综述**：胞外 e-PHA 解聚酶（EC 3.1.1.75/76）与胞内 i-PHA 解聚酶的生化/分子全貌，PhaZ 分类起点 |
| 9008883 | 1996 | Appl Microbiol Biotechnol | Jendrossek D | 聚羟基烷酸（PHA）的生物降解 | 早期系统综述：降解菌分离鉴定方法、解聚酶生化性质、水解机制与调控 |
| 9921137 | 1999 | Rev Environ Contam Toxicol | Hankermeyer CR | PHB：微生物制造并降解的塑料 | PHB/PHBV（Biopol）微生物合成与降解机制综述（含环境归趋） |
| 38272380 | 2024 | Biotechnol Adv | Park H | PHA 不只是生物塑料！ | PHA 应用全景综述（包装/纺织/废水反硝化碳源等），背景材料 |
| 33076314 | 2020 | Materials | Kliem S | 不同环境中聚合物的生物降解综述 | 环境条件×聚合物降解性（含 PHB/PHBV），环境因子背景 |
| 30455023 | 2018 | Waste Manag | Bátori V | 生物塑料的厌氧降解综述 | 厌氧消化条件下 PHB 等生物塑料降解性 |
| 24632193 | 2014 | Curr Opin Biotechnol | Meng DC | 聚酯多样性的工程化 | PHA 合成途径工程（β-氧化减弱等），合成-降解平衡背景 |

## 三、论文 URL（PubMed）

- https://pubmed.ncbi.nlm.nih.gov/12213937/ （Jendrossek 2002 综述）
- https://pubmed.ncbi.nlm.nih.gov/9008883/ （Jendrossek 1996 综述）
- https://pubmed.ncbi.nlm.nih.gov/9921137/ （Hankermeyer 1999 综述）
- https://pubmed.ncbi.nlm.nih.gov/38272380/ （Park 2024 综述）
- https://pubmed.ncbi.nlm.nih.gov/33076314/ （Kliem 2020 综述）
- https://pubmed.ncbi.nlm.nih.gov/30455023/ （Bátori 2018 综述）
- https://pubmed.ncbi.nlm.nih.gov/24632193/ （Meng 2014）
- https://pubmed.ncbi.nlm.nih.gov/16405909/ （Hisano 2006 晶体结构）
- https://pubmed.ncbi.nlm.nih.gov/9406404/ （Kasuya 1997）
- https://pubmed.ncbi.nlm.nih.gov/24907326/ （Sznajder 2014 PhaZ6/7）
- https://pubmed.ncbi.nlm.nih.gov/29678915/ （Juengert 2018 磷酸化）
- https://pubmed.ncbi.nlm.nih.gov/16936025/ （Tseng 2006 Bt phaZ）
- https://pubmed.ncbi.nlm.nih.gov/20346169/ （Trainer 2010 S. meliloti）
- https://pubmed.ncbi.nlm.nih.gov/12898135/ （Kadouri 2003 Azospirillum）
- https://pubmed.ncbi.nlm.nih.gov/17170116/ （de Eugenio 2007 mcl-PHA）
- https://pubmed.ncbi.nlm.nih.gov/19788655/ （de Eugenio 2010 mcl-PHA 周转）
- https://pubmed.ncbi.nlm.nih.gov/19214501/ （Papaneophytou 2009 Thermus）
- https://pubmed.ncbi.nlm.nih.gov/15489436/ （Handrick 2004 R. rubrum）
- https://pubmed.ncbi.nlm.nih.gov/23951224/ （García-Hidalgo 2013 Streptomyces）
- https://pubmed.ncbi.nlm.nih.gov/41151231/ （Anjulal 2026 Nocardiopsis）
- https://pubmed.ncbi.nlm.nih.gov/36222314/ （Thomas 2022 Lihuaxuella）
- https://pubmed.ncbi.nlm.nih.gov/26156240/ （Martínez 2015 mcl-PHA）
- https://pubmed.ncbi.nlm.nih.gov/24751505/ （Martínez 2014 细胞工程）
- https://pubmed.ncbi.nlm.nih.gov/38969063/ （Santolin 2024 转录调控）
- https://pubmed.ncbi.nlm.nih.gov/34342882/ （Nygaard 2021 C. necator）
- https://pubmed.ncbi.nlm.nih.gov/29951858/ （Martínez-Tobón 2018 基因组筛选）
- https://pubmed.ncbi.nlm.nih.gov/32369496/ （Morohoshi 2020 Undibacterium 全基因组）
- https://pubmed.ncbi.nlm.nih.gov/40500476/ （Iseki 2025 Alteromonas 系统发育分布）
- https://pubmed.ncbi.nlm.nih.gov/36125959/ （Leadbeater 2022 in silico 挖掘）
- https://pubmed.ncbi.nlm.nih.gov/38278791/ （Omura 2024 深海降解）
- https://pubmed.ncbi.nlm.nih.gov/34160268/ （Eronen-Rasimus 2021 海冰细菌）
- https://pubmed.ncbi.nlm.nih.gov/40320445/ （Boondaeng 2025 海洋降解菌）
- https://pubmed.ncbi.nlm.nih.gov/37311705/ （Jeon 2023 B. infantis）
- https://pubmed.ncbi.nlm.nih.gov/40392676/ （Hachisuka 2025 人工 PHA 共聚物）
- https://pubmed.ncbi.nlm.nih.gov/26409775/ （Chen 2015 PHAome）

## 四、对 GTDB 项目的方法学启示

1. **PhaZ 基因家族存在"命名混乱"**：Sznajder 2014 与 Tseng 2006 显示 PhaZ 序号（PhaZ1–PhaZ7 等）与功能注释（如 PcaD→phaZ）在不同菌中不一致——GTDB 同源扫描需结合隐马尔可夫模型（HMM）+催化位点（Ser-Asp-His/lipase box）验证，而非仅依赖基因名。
2. **胞内/胞外解聚酶需区分**：Handrick 2004（周质定位）、Jendrossek 2002 综述给出分类框架（e-PHA 解聚酶 EC 3.1.1.75/76 vs i-PHA 解聚酶）；结构域特征（信号肽、底物结合域、连接区）可作为功能分类标记。
3. **基因组/宏基因组筛选已有先例**：Martínez-Tobón 2018（预测 vs 实测活性）、Morohoshi 2020（全基因组定位）、Iseki 2025（系统发育分布）可直接借鉴工作流；phylogenetic distribution 分析模式与 GTDB 树结合是可行思路。
4. **数据库构建参考**：Chen 2015 "PHAome" 提供概念框架；Leadbeater 2022 提供数据库挖掘 pipeline 参考。
5. **文献缺口提示**：未检索到"基于 GTDB 的 PHB 降解基因系统筛查"这一特定组合的直接先例文献——即本项目具有新颖性；6/7 查询中基因组+降解组合命中的多为生产/合成方向论文，降解酶基因组筛查文献相对少。
