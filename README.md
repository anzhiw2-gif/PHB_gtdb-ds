# PHB_gtdb-ds — GTDB 全库 PHB 降解基因系统生信分析

基于 GTDB（Genome Taxonomy Database）R232 全基因组集合（199,923 个代表基因组，
含细菌+古菌）对 PHB（聚羟基丁酸酯）降解基因的系统筛查与生态/分类学分析。

## 当前状态

**主流程已跑通并出结果**（服务器 T141，2026-08-16）。核心结论：

- **9 家族 HMM 全库筛选**（HMMER 3.4，E<1e-5），覆盖胞内/胞外解聚酶、寡聚体
  水解酶、古菌 patatin/经典酯酶 + 背景代谢（BdhA）+ 辅助家族（PhaJ/phasin/PhaC）。
- **tier1 严格集 ~75,000 条解聚酶序列 / ~57,000 基因组（28%）**，其中
  ePhaZ 38,275、iPhaZ 32,846、OH 1,470、ArchPhaZ_hydrolase 1,292、
  ArchPhaZ_patatin 1,372。
- **古菌 PHB 降解基因首次 GTDB 全库证实**：Halobacteriota + Thermoproteota 检出
  patatin/经典酯酶，且 patatin 型扩展至 Nitrososphaeria（氨氧化古菌/AOA）——
  超出文献已知范围（文献仅实验表征 Haloferax mediterranei PhaZh1）。
- 对比 Viljakainen & Hug 2021（宏基因组 13,869 条）为量级提升。

详见 [docs/final_results_report.md](docs/final_results_report.md)。

## 目录结构

```
├── docs/                          # 报告与计划
│   ├── literature_survey_report.md   # 文献调研报告（阶段一）
│   ├── analysis_plan_draft.md        # 分析流程设计（含执行进度）
│   └── final_results_report.md       # 最终结果报告（主交付物）
├── knowledge/                     # 领域知识框架（家族分类/种子清单）
│   └── family_classification.md      # 家族分类方案（唯一口径，§6.4）
├── pipeline/
│   ├── config/params.yaml            # 参数配置
│   ├── scripts/                      # 01-11 主流程脚本
│   ├── dev/                          # 探索性脚本归档
│   └── seeds/                        # 种子序列 + manifest（78 curated）
├── research/                      # 调研原始数据（PubMed/OpenAlex/EuropePMC/Web）
├── data/                          # GTDB 数据与中间产物（服务器侧，不入 git）
└── results/                       # 结果表/图/树（表与图入 git，大文件走同步脚本）
```

## 家族分类（9 家族）

| 类型 | 家族 |
|------|------|
| 核心解聚酶（进 tier） | ePhaZ、iPhaZ、OH、ArchPhaZ_patatin、ArchPhaZ_hydrolase |
| 背景代谢 | BdhA（3HB 脱氢酶，不计入解聚酶计数） |
| 辅助（簇背景） | PhaJ、phasin、PhaC |

完整分类依据与催化位点规则见
[knowledge/family_classification.md](knowledge/family_classification.md)。

## 复现

1. 环境：`conda env create -f environment.yml`（服务器 T141 已就绪 `phb_gtdb`）。
2. 参数：`pipeline/config/params.yaml`；数据溯源见 `docs/reproducibility.md`。
3. 流程：按 `pipeline/README_HPC.md` 的脚本顺序执行。
4. 结果同步：`bash pipeline/sync_from_server.sh`（需服务器 SSH 访问）。

## 待完成

- [x] 生态元数据（isolation source）关联（`10_distribution.py`，已完成）
- [ ] patatin 位点级 ±flank_kb 基因簇共定位验证（`pipeline/scripts/11_clusters.py`，运行中，620 古菌基因组）
- [ ] ePhaZ SignalP 胞外/胞内细分
- [ ] 完整（非抽样）系统发育树 + 基因树×物种树共进化（HGT 检测）
- [ ] 整理成论文 / 发布 HMM profiles + 命中表（GitHub + Zenodo DOI）
