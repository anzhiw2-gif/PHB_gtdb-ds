# 复现性与数据溯源

> 本文档固化数据版本、环境、参数与运行记录，保证结果可复现。

## 1. 数据源（只读，GTDB R232）

| 项 | 值 |
|----|----|
| Release | **R232**（R11-RS232） |
| 下载页 | https://data.gtdb.ecogenomic.org/releases/latest/ |
| 统计页 | https://gtdb.ecogenomic.org/stats/r232 |
| 代表基因组 | 199,923（细菌 + 古菌） |
| 基因组目录（服务器） | `~/GTDB/gtdb_genomes_reps_r232/database/{GCA|GCF}/...` |
| 元数据 | `~/GTDB/metadata/bac120_metadata_r232.tsv.gz` |
| 分类学 | `~/GTDB/taxonomy/bac120_taxonomy_r232.tsv` |
| 骨架树 | `~/GTDB/GTDB_tree/bac120_r232.tree` |

> **待补**：下载完成后记录各文件的 SHA-256 校验和（`sha256sum`），
> 固定数据版本。示例：
> ```bash
> sha256sum gtdb_proteins_aa_reps.tar.gz gtdb_metadata.tsv > data/checksums.sha256
> ```

## 2. 环境

- 定义：`environment.yml`（本仓库根目录）
- 服务器实例：`conda env create -f environment.yml`，激活 `conda activate phb_gtdb`
- 精确锁定：`conda env export -n phb_gtdb > environment.lock.yml`（建议在服务器固化）

## 3. 参数快照

- `pipeline/config/params.yaml`（筛选 E=1e-5、Pyrodigal meta 模式、tier 阈值 E<1e-20/1e-10、
  每 shard 2000 基因组、系统发育 LG+G4/1000 UFBoot 等）
- tier 定义（`scripts/08c_tier_rescore.sh`）：curated HMM 重评分，
  tier1 E<1e-20、tier2 E<1e-10、tier3 宽模型+通用验证。

## 4. 运行记录

- 预测日志：`results/logs/predict_full.log`
- 筛选日志：`results/logs/screen_{family}.log`
- 各步关键命令见 `pipeline/README_HPC.md`。

## 5. 数据沿革（脚本顺序）

```
GTDB R232 基因组（199,923）
  → 05_predict_proteins.sh（Pyrodigal meta，~4.7 亿蛋白）
  → filter_long_seqs.py（剔除 >100K aa 伪影）
  → 06_screen.sh（9 家族 HMMER E<1e-5）
  → 07_process_hits.py / 07b_extract_seqs.py（过滤/仲裁/去重/提取）
  → 08_validate.py / 08c_tier_rescore.sh（催化位点验证 + 三级重评分）
  → 09a-09h（tier1 汇总/系统发育/图/patatin 过滤）
  → 10_distribution.py（分类/生态分布）
  → 11_clusters.py（±flank_kb 基因簇共定位）
```

## 6. 结果同步（服务器 → 本地）

```bash
bash pipeline/sync_from_server.sh          # 同步 tables/figures/trees + HMM + tier1 序列
bash pipeline/sync_from_server.sh --dry-run
```

- 交付物（HMM profiles、tier1 FASTA、结果表）体积小，建议以 GitHub Release +
  Zenodo DOI 发布；`data/` 大中间产物不入 git。
