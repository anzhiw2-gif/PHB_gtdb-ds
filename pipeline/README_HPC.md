# PHB_gtdb-ds — 服务器分析工作区说明

> 工作区：通过 `PHB_REPO_ROOT` 配置（远程主机信息不写入仓库）
> 边界：**工作区外（~/GTDB、~/PHB_gtdb 等）只读，不修改**
> 数据：GTDB R232 代表基因组（只读，通过 `PHB_GTDB_ROOT` 配置）
> 环境：`conda activate phb_gtdb`（HMMER 3.4 / DIAMOND 2.2.1 / Pyrodigal 3.7.1 /
>        MAFFT / trimAl / IQ-TREE / CD-HIT / GNU parallel）
> 环境定义：`environment.yml`（本仓库根目录）

## 分析目标

基于文献调研（docs/literature_survey_report.md）建立的完整 PHB 降解基因
目录，对 GTDB R232 全库（199,923 个代表基因组）进行系统筛选。

## 家族分类体系（唯一口径，见 knowledge/family_classification.md §6.4）

共 **9 个 HMM 家族** = 5 个核心解聚酶 + 1 个背景代谢 + 3 个辅助：

| 家族 | 类型 | 功能 | 建 HMM 种子（v2，去冗余前→后） |
|------|------|------|-------------------------------|
| ePhaZ | 核心解聚酶 | 胞外 PHB/PHA 解聚酶（EC 3.1.1.75/76） | 4,458 → 3,002 |
| iPhaZ | 核心解聚酶 | 胞内 PHA 解聚酶（Cys/Ser/周质型） | 152 → 112 |
| OH | 核心解聚酶 | 3HB 寡聚体水解酶（EC 3.1.1.22） | 713 → ~500 |
| ArchPhaZ_patatin | 核心解聚酶（古菌） | patatin 样 PhaZh1 型 | 113 → 103 |
| ArchPhaZ_hydrolase | 核心解聚酶（古菌） | 经典 PHB 解聚酶家族酯酶 | 12 |
| BdhA | 背景代谢 | 3HB 脱氢酶（EC 1.1.1.30，不计入解聚酶） | 5,903 → ~4,000 |
| PhaJ | 辅助 | 烯酰-CoA 水合酶（MaoC/PF01575） | 1,071 |
| phasin | 辅助 | 颗粒蛋白（PF09361） | Pfam |
| PhaC | 辅助 | PHA 合酶（簇共定位背景） | 47 |

- **核心解聚酶 = 前 5 个家族**，进入 tier1/tier2/tier3 三级重评分。
- **实验表征的 curated 种子 = 78 条**（`pipeline/seeds/seeds_manifest.tsv`）；
  v2 建 HMM 用的是各家族扩充种子（上表）。
- patatin 家族细菌 HMM 检不出，须用 Pfam Patatin 独立筛选 + 基因簇验证。

## 流程（按序执行）

```bash
conda activate phb_gtdb
cd ~/PHB_gtdb-ds

# 1. 种子收集（02 + 02b/02c/02d/02e 为多轮扩充）
python scripts/02_collect_seeds.py --outdir data/seeds
# ...（其余轮次见各脚本 docstring）

# 2. 分家族整理 + HMM 构建（v2）
python scripts/03_prep_families.py --seeds data/seeds/seeds_family.faa --outdir data/seeds/families
bash scripts/04b_build_hmms_v2.sh --threads 40

# 3. 全库蛋白预测（Pyrodigal, ~40-48h）
nohup bash scripts/05_predict_proteins.sh --threads 80 > results/logs/predict_full.log 2>&1 &

# 3b. 过滤超长伪影序列（>100K aa）
python scripts/filter_long_seqs.py

# 4. 全库 HMM 筛选（9 家族）
bash scripts/06_screen.sh --threads 80 --eval 1e-5

# 5. 命中处理与序列提取
python scripts/07_process_hits.py --hits data/screen/hits_all.tsv
python scripts/07b_extract_seqs.py

# 6. 功能验证 + 三级重评分
python scripts/08_validate.py --signalp 0
bash scripts/08c_tier_rescore.sh

# 7. tier1 汇总 / 系统发育 / 图（09a-09h）
python scripts/09a_tier1_summary.py
bash scripts/09b_tier1_phylogeny.sh --threads 40
python scripts/09c_tier1_figures.py
python scripts/09d_patatin_filter.py   # patatin 基因组级 PhaC 过滤

# 8. 生态/分类学分布
python scripts/10_distribution.py --hits data/screen/genome_hits.tsv

# 9. 基因簇共定位（±flank_kb，含 patatin 二次过滤）
python scripts/11_clusters.py --hits data/screen/hits_filtered.tsv --max-genomes 0
```

> 脚本编号说明：`02b-02e`、`04b`、`06a/06b`、`07b`、`08b/08c`、`09a-09h`
> 是对应主步骤的**变体/细分**（如 02b-d 为种子扩充轮次、04b 为 v2 HMM、
> 09a-h 为 tier1 后处理）。探索性脚本已归档至 `pipeline/dev/`。

## 目录结构

```
PHB_gtdb-ds/
├── config/            # params.yaml（参数）
├── scripts/           # 01-11 主流程脚本
├── dev/               # 探索性/一次性脚本归档
├── seeds/             # 种子序列 + manifest + stats
├── data/              # 蛋白预测/筛选/tier（服务器侧，不入 git）
├── results/
│   ├── tables/        # 结果表
│   ├── trees/         # 系统发育树
│   ├── figures/       # 图表
│   └── logs/          # 运行日志
└── docs/              # 报告与说明
```

## 监控

```bash
# 预测进度
ls data/proteins/per_genome/*.faa.gz | wc -l        # 目标 199,923
tail results/logs/predict_full.log
```
