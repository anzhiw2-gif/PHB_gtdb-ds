# PHB_gtdb-ds — GTDB 全库 PHB 降解基因系统生信分析

> Current snapshot (reviewed 2026-09-03): formal frozen scan 13 and its downstream tier processing are complete. Run-13 results are candidate homology evidence, not phenotype validation. See `docs/STATUS.md` for the authority map.

> Historical Scheme A values remain in the reports for comparability. They must not be mixed with the run-13 frozen split registry or strict tier results.

> Public repository guide: see [docs/PUBLIC_REPOSITORY_GUIDE.md](docs/PUBLIC_REPOSITORY_GUIDE.md) for environment-variable configuration and the sensitive-data boundary.

基于 GTDB（Genome Taxonomy Database）R232 全基因组集合（199,923 个代表基因组，
含细菌+古菌）对 PHB（聚羟基丁酸酯）降解基因的系统筛查与生态/分类学分析。

## 当前状态

> **项目单一状态页见 [docs/STATUS.md](docs/STATUS.md)**：结论边界措辞、数据流契约、
> 线程上限、树状态与待办均以该页为准；本 README 仅作摘要。

**当前发布快照**为 GitHub `main`（审阅日期：2026-09-03）。正式 frozen scan 13 已完成，运行目录和原始
HMMER 输出保留在服务器侧；公开仓库只保留轻量结果、模型、脚本和可复现契约。当前运行证据与
严格 tier 结果见 `docs/CURRENT_STATUS_20260902.md` 和 `docs/T141_20260902_formal_scan13_tier_processing_02_status.md`。

- **9 家族 HMM 全库筛选**（HMMER 3.4，E<1e-5），覆盖胞内/胞外解聚酶、寡聚体
  水解酶、古菌 patatin/经典酯酶 + 背景代谢（BdhA）+ 辅助家族（PhaJ/phasin/PhaC）。
- **run-13 原始 registry 命中**：6,740,900 条 accepted hits，核心四家族（ePhaZ 合并 curated/broad、
  iPhaZ、OH、ArchPhaZ_hydrolase）并集为 **147,690 个基因组**。这是阈值命中层，不是 tier1。
- **run-13 严格 tier1**：38,741 个基因组；6,578 个基因组含至少两个核心家族。各家族为：
  ePhaZ curated 5,080、iPhaZ 25,564、OH 3,446、ArchPhaZ_hydrolase 12,469 个基因组。
- 严格 tier1 序列数分别为 ePhaZ 5,646、iPhaZ 32,226、OH 3,570、ArchPhaZ_hydrolase 14,571；
  ePhaZ broad discovery 520,217 条蛋白记录保持独立，不晋升为严格 tier1。
- **古菌谱系检出 PHB 降解相关候选同源蛋白**（功能潜力，非实证）：核心解聚酶
  （ePhaZ/经典酯酶）候选在 Halobacteriota（112 基因组）、Thermoproteota（32）、
  Thermoplasmatota（154）检出。**注意**：patatin 折叠蛋白（古菌子集 620 基因组）
  是广谱磷脂酶结构域；Figure 5 的邻域 support rate 使用 all-hit candidate-loci denominator，
  不能直接解释为该 620 基因组子集的比例。PhaZh1 体内角色有限、PhaJ 才是动员主路，故"patatin
  型解聚酶"应审慎表述
  （详见 final_results_report.md §2.2）。
- 对比 Viljakainen & Hug 2021（宏基因组 13,869 条）为量级提升。

图源口径：Figure 3 source data 的 Pseudomonadota 为 **26,850**；历史输出中的 26,855 不再作为当前
值。`ArchPhaZ_patatin` 的广谱候选层为 112,926 条蛋白记录，古菌讨论子集为 1,372 loci / 620 genomes，
不能混用。

树状态：OH 旧输入树为 `stale_input`，ePhaZ CD-HIT 树为 `input_not_registered`；ePhaZ/iPhaZ
全量树和 HGT 检测仍暂停，不能把现有抽样/历史树写成完整系统发育证据。

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
| 核心解聚酶（进 tier） | ePhaZ、iPhaZ、OH、ArchPhaZ_hydrolase |
| 广谱折叠（需基因簇过滤） | ArchPhaZ_patatin（patatin 折叠，多为磷脂酶；Figure 5 的 all-hit 邻域率不能当作该子集比例） |
| 背景代谢 | BdhA（3HB 脱氢酶，不计入解聚酶计数） |
| 辅助（簇背景） | PhaJ、phasin、PhaC |

完整分类依据与催化位点规则见
[knowledge/family_classification.md](knowledge/family_classification.md)。

## 复现

### 新版运行目录约定（2026-08-24）

本地新版主编排默认创建一次性目录 `runs/<run_id>/`，不会覆盖仓库根目录的
`data/` 或 `results/`。运行目录包含 `logs/`、`inputs/`、`results/`、
`run_context.env` 和 `input_contract.json`；HMM 目录以只读符号链接接入。

```bash
bash pipeline/scripts/run_pipeline.sh --run-id 20260824T120000Z_demo
bash pipeline/scripts/run_pipeline.sh --legacy-root-results
```

`input_contract.json` 记录 GTDB taxonomy、metadata、tree、全部 HMM、
`environment.yml` 和 `pipeline/config/params.yaml` 的路径、大小和 SHA-256。
缺失 GTDB 文件标记为 `pending`，不生成占位哈希。run-13 的输入契约已记录 GTDB R232
taxonomy/metadata/tree 与模型哈希；原始 GTDB 数据不进入 Git。

1. 环境：`conda env create -f environment.yml`。
2. 参数：`pipeline/config/params.yaml`；数据溯源见 `docs/reproducibility.md`。
3. 流程：按 `pipeline/README_HPC.md` 的脚本顺序执行。
4. 结果同步：`bash pipeline/sync_from_server.sh`（需服务器 SSH 访问）。

## 当前待办

- [x] 生态元数据（isolation source）关联（`10_distribution.py`，已完成）
- [x] 用最终方案 A 输入完成 dated patatin 位点级 ±flank_kb 基因簇共定位（80 批）；Figure 5 仍是 all-hit 邻域层，古菌 patatin 子集需另行按 1,372 loci / 620 genomes 汇总，旧 `cluster_summary.tsv` 不作为最终计数依据
- [x] ePhaZ SignalP 胞外/胞内细分（详见 `docs/STATUS.md`）
- [x] ~~patatin 用 PhaZh1 专属种子重建 HMM~~（已判定不可行：patatin 催化域序列高度保守，
      同源法无法区分；以 ±10kb 基因簇共定位为判据，见 STATUS.md §5）
- [x] 完成 run-13 全库 HMM 扫描及严格 tier1 下游处理（见上述状态文档）
- [ ] 复核 run-13 strict/broad/contextual 分层结果并决定论文统计口径
- [ ] 重建并登记使用当前输入的 OH 树；全量 ePhaZ/iPhaZ 树和 HGT 仍暂停
- [ ] 整理成论文；发布 HMM profiles + 轻量命中表（GitHub Release + Zenodo DOI）

## License

Repository-authored source code and documentation are released under the [MIT License](LICENSE).
GTDB files, external software, datasets, and other third-party materials remain under their
respective upstream licenses and terms.
