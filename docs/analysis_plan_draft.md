# GTDB PHB 降解基因系统生信分析 — 流程设计（v0.4，执行中）

> 状态：v0.4。2026-06-01 已进入实际执行阶段（服务器 T141, <SERVER_HOST>）。
> 完整文献依据见 docs/literature_survey_report.md；服务器部署说明见
> pipeline/README_HPC.md；参数见 pipeline/config/params.yaml。

## 0. 目标

在 GTDB 全基因组集合（代表性基因组）中系统鉴定 PHB 降解相关基因：
- 建立高质量的 PHB 降解基因参考序列集 + HMM profile（胞内/胞外解聚酶、
  寡聚体水解酶、3HB 代谢酶）
- 全库筛选 → 功能注释（含结构域验证、催化三联体检查）
- 系统发育分析（基因树 + 物种树映射）
- 生态/分类学分布统计（基于 GTDB 分类与基因组元数据）
- 交付：可复现流程（脚本）+ 结果报告 + 数据表

## 1. 数据获取（GTDB，已核实）

- 最新 Release：R11-RS232（R10 = 715,230 细菌 + 17,245 古菌，Parks et
  al. 2025）；统计页 https://gtdb.ecogenomic.org/stats/r232
- 下载：https://data.gtdb.ecogenomic.org/releases/latest/ 下
  - gtdb_proteins_aa_reps.tar.gz（主数据源，每基因组 .faa）
  - gtdb_metadata.tsv（分类/质量/来源）、taxonomy 文件
  - （可选）gtdb_genome_reps.tar.gz 用于基因簇上下文
- faa header→GTDB 分类学映射：GTDB Forum 线程 599；
  直接做 DIAMOND 化：James Lingford 博客流程
- 规模：~73 万基因组、蛋白库数百 GB → 评估带宽/磁盘后分块处理

## 2. 参考序列集与 HMM 构建（种子来源已核实）

1. 种子序列来源：
   - **PhaDED/ESTHER**（Esterase_phb_PHAZ、PHAZ7_phb_depolymerase 家族页）
   - **PAZy** PHB 解聚酶条目（API/DaRUS 可下载）
   - **PlasticDB** HMM 模型
   - UniProt reviewed（EC 3.1.1.75，3 条）+ 文献已表征酶
     （详见 knowledge/seed_sequences_checklist.md）
2. 分家族建模（胞内 i-PhaZ / 胞外 e-PhaZ / 寡聚体水解酶 / 3HB 代谢酶）：
   - 去冗余（CD-HIT / MMseqs2 easy-cluster）
   - 多序列比对（MAFFT）
   - HMM（HMMER3 hmmbuild）
   - 阈值校准（对 GTDB 子集做 roc/分布图）
3. 双重验证：结构域（NCBI CDD、ESTHER 家族）+ 催化三联体
   （Ser-Asp-His）+ lipase box（G-X-S-X-G）+ 信号肽（SignalP）

## 3. 全库筛选（范式：Zrimec 2021 mBio 已核实）

- 方案 A（蛋白级，主）：DIAMOND/MMseqs2 预筛（快速）→ HMMER
  v3.3 hmmsearch 精筛（复刻 Zrimec 2021 流程）→ 过滤
  （EVALUE、bias、覆盖度）→ 催化位点验证
- 方案 B（基因级，可选，保留上下文）：代表性基因组 → Prodigal →
  同前筛选 → 提取基因簇（±10 kb 邻域）分析
- 输出：per-genome 检出表（基因 ID、类型、GTDB 物种/门、坐标）
- 可用现成流程参考：plastedma、PDETool；PlasticEnz 交叉验证

## 4. 系统发育分析

- 每类酶：去冗余序列（如 90% ID）→ MAFFT 比对 → trimAl 修剪 →
  IQ-TREE2（LG+G4 或自动选择，1000 UFBoot）→ ETE3/iTOL 可视化
- 标注：催化三联体完整性、结构域组合、GTDB 门/纲、生态来源
- 与 GTDB 骨架树做共系统发育比较（HGT vs 垂直遗传）
- 可选：PEGP 补域架构与 phyletic pattern；anvio 做基因上下文

## 5. 生态分布与统计

- 按 GTDB 门/纲、基因组来源生态（metadata: isolation source）
- 丰度：每类基因在门水平的检出率（% genomes）
- 多拷贝/缺失模式、基因簇共现（与 phaCAB 合成簇邻近关系）
- 与 Viljakainen & Hug 2021（宏基因组 13869 个 PHA 解聚酶）结果对照，
  突出 GTDB 全覆盖 + 分类学尺度增量价值
- 输出图表：堆叠条形图、热图、树标注图

## 6. 交付物

- 流程脚本（Snakemake 或 Python 脚本 + 说明）
- 数据表：筛选结果 tsv、基因簇表、统计表
- 报告：方法、结果、讨论（中文）
- 图：分布图、系统发育树图

## 7. 计算环境

- 当前 Windows 无 HMMER/MAFFT/DIAMOND/Prodigal/IQ-TREE，无 WSL
  （实测）→ 建议 WSL2 + Miniconda 或远程 Linux/HPC
- 工具依赖：HMMER 3.3、MAFFT、DIAMOND/MMseqs2、Prodigal、IQ-TREE2、
  ETE3、CD-HIT、Python（biopython、pandas）、Snakemake
- 规模预估与分块策略（大库需按基因组分块并行）

## 8. 风险与对策

- 假阳性：HMM 阈值过松 → 结构域+催化位点双重验证
- 胞内/胞外解聚酶同源度低 → 分家族建模而非单一模型
- 旁系同源干扰（非特异脂肪酶类）→ 底物/结构域特征过滤
- GTDB 蛋白质文件版本更新 → 固定 release、记录校验和
- 计算量：~73 万基因组的 HMMER/DIAMOND → 分块 + 去冗余后抽样验证

## 9. 执行进度（服务器 T141）

| 步骤 | 状态 | 说明 |
|------|------|------|
| 工作区搭建 | ✅ | ~/PHB_gtdb-ds（边界：外只读） |
| 种子收集 | ✅ | curated 78（manifest）+ v2 扩充种子，9 家族 |
| HMM 构建 | ✅ | 9 家族（ePhaZ/iPhaZ/OH/ArchPhaZ_patatin/ArchPhaZ_hydrolase + BdhA + PhaJ/phasin/PhaC） |
| 蛋白预测 | 🔄 运行中 | Pyrodigal 199,923 基因组, 80 并行, ~50h |
| HMM 全库筛选 | ⏳ | 06_screen.sh（脚本已测通） |
| 命中处理 | ⏳ | 07_process_hits.py（脚本已测通） |
| 功能验证 | ⏳ | 08_validate.py（家族定制规则） |
| 系统发育 | ⏳ | 09_phylogeny.sh |
| 生态分布 | ⏳ | 10_distribution.py |
| 基因簇分析 | ⏳ | 11_clusters.py |
| 最终报告 | ⏳ | |

## 10. 环境评估（2026-06-01 实测）

当前 Windows 工作机（<LOCAL_WORKSPACE>）：
- ✅ 可用：Python 3.12（含 ete3）、uv、PowerShell
- ❌ 未安装：HMMER(hmmsearch/hmmbuild)、MAFFT、DIAMOND、Prodigal、IQ-TREE2、
  CD-HIT、MMseqs2、trimal、Snakemake、SignalP、ClustalO
- ❌ WSL 未安装（wsl.exe --install 可启用）

对策（用户已选定）：**远程 Linux/HPC（T141: <SERVER_HOST>）**
- 服务器：Ubuntu 24.04, 80 核, 1TB 内存, 82TB 可用磁盘
- conda 环境 phb_gtdb：HMMER 3.4 / DIAMOND 2.2.1 / Pyrodigal 3.7.1 /
  MAFFT / trimAl / IQ-TREE / CD-HIT / GNU parallel
- GTDB R232 数据（只读）：199,923 个代表基因组 + bac120 元数据/分类/树 + Pfam
- 数据下载：GTDB（数 TB 级）已由用户完成，不重复下载
