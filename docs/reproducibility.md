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

> **待补（高优先）**：下载完成后记录各文件的 SHA-256 校验和，固定数据版本。
> 目标文件：`gtdb_proteins_aa_reps.tar.gz`、`gtdb_metadata.tsv`、taxonomy、骨架树。
> 示例（在服务器执行并提交）：
> ```bash
> sha256sum ~/GTDB/gtdb_genomes_reps_r232/database/checksums.txt \
>   > data/gtdb_checksums.sha256
> ```
> 运行时的输入/输出哈希由 `run_manifest.json` 自动记录（见 §7）。

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
- **运行清单**：当前可验证清单位于 T141 dated 运行目录
  `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260821_schemeA_03/results/run_manifest.json`，由
  `run_manifest.py` 生成，记录每步 exit code、起止时间和输入/输出 SHA-256；其当前 SHA-256 为
  `e991c3bd10a48f8faf9c450f0c17a5a3fb1f0315c256018f0772c5f64f71b2a3`。本地仓库没有同一份可验证
  manifest，不应倒填或伪造。
- **系统发育树清单**：`results/trees_tier1/tree_manifest.tsv` 记录工具、叶数、历史输入和当前输入哈希、抽样名单及状态；`stale_input` 和 `input_not_registered` 的树不得作为当前候选集的树证据。

### 4.1 代码、部署与 manifest 边界

Publication check (2026-08-27): GitHub `main` is `94caa49`, while this local
working tree still contains unpublished research changes. Older commit references
below describe historical audit snapshots and are not the current GitHub revision.

本地 `<LOCAL_WORKSPACE>` 与 GitHub `main` 当前均为 `94caa49`，但本地工作树有未提交修改和未跟踪
审计/结果文件。T141 项目根目录没有 `.git`，主 manifest 的 `git_commit` 为 `null`，因此该清单
并未绑定到 Git 源码快照。dated forensic repair 使用
`${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260821_schemeA/scripts/11_clusters.py`
（SHA-256 `d1d907f34a5c1fbe17aee538a3ca087de8e809f5c03286b6384a22347d994d2c`）；服务器根目录
`scripts/11_clusters.py` 是另一版本。后续运行应从 dated `deploy/<run_id>/scripts/` 或明确绑定的
Git commit 启动，并将源码包、完整环境、GTDB metadata/taxonomy/tree、全部 HMM 与完整命令写入
同一 manifest。旧文档中的 manifest SHA `c4e8c73b...` 已过时。

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

## 7. 运行清单（run_manifest.json）与完整命令清单

本地新版 `run_pipeline.sh` 默认将每次执行写入独立的 `runs/<run_id>/`，并在
启动阶段生成同目录的 `input_contract.json`。最终 manifest 额外绑定 `run_id`、
`run_root` 和该 contract；GTDB taxonomy/metadata/tree 未提供时保持 `pending`，
不得用其他文件冒充。复用仓库根目录必须显式传入 `--legacy-root-results`。

`run_pipeline.sh` 每步退出时向 `results/run_manifest.jsonl` 追加一行
`{step, exit_code, started, ended, note, command}`，最后调用 `run_manifest.py finalize`
固化为 `results/run_manifest.json`。本地新版在 `--strict-provenance` 下还要求并写入：

- 源码文件哈希与确定性的 `source_bundle_sha256`；
- Python/Conda/平台/cwd 环境信息；
- GTDB 输入、实际 HMM 输入（包括 06 的 v2 与 08c 的 legacy core 模型）哈希；
- 每个成功步骤的完整 argv 命令。

清单同时含：

- 元数据：生成时间、主机、git commit（当前 T141 主运行实际为 `null`，因此不能视为源码绑定）；
- 每步 exit code（fail-closed：任一步非零即中止）；
- 输入/输出文件 SHA-256（`--inputs` / `--outputs` 列表）。

若任何声明的输入或输出缺失或为空，`finalize` 返回非零且不会写出成功的运行清单。

**完整命令清单**（服务器按序执行，见 `pipeline/README_HPC.md`）：

```bash
conda activate phb_gtdb
cd ~/PHB_gtdb-ds
bash pipeline/scripts/05_predict_proteins.sh --threads 70        # 蛋白预测
bash pipeline/scripts/06a_filter_shards.sh                        # 过滤超长序列
bash pipeline/scripts/06_screen.sh --threads 70 --eval 1e-5      # 全库筛选(+06b 覆盖度)
conda run -n phb_gtdb python pipeline/scripts/07_process_hits.py \
    --hits data/screen/hits_all.tsv --shards data/proteins/shards_filt
conda run -n phb_gtdb python pipeline/scripts/08_validate.py
bash pipeline/scripts/08c_tier_rescore.sh
bash pipeline/scripts/09b_tier1_phylogeny.sh --threads 40        # 或 09_phylogeny.sh
conda run -n phb_gtdb python pipeline/scripts/10_distribution.py \
    --hits data/screen/genome_hits.tsv
conda run -n phb_gtdb python pipeline/scripts/11_clusters.py \
    --hits data/screen/hits_filtered.tsv                          # 必须含 locus
python pipeline/scripts/09i_tree_manifest.py                      # 树登记（仅在明确执行建树后）
```

> 命令块是本地 git 布局的示例。T141 为扁平 `scripts/`/`data/` 布局，不能直接把本地
> `pipeline/scripts/` 路径套到服务器；服务器运行只允许使用 dated `deploy/<run_id>/scripts/` 或
> 明确绑定的 Git 源码快照。当前治理状态下建树默认暂停；如已获得明确授权，使用
> `--run-phylogeny` 显式开启建树和树登记。

## 8. 小样本端到端测试与 CI

- 冒烟测试：`pipeline/tests/test_smoke.sh` —— 用 3–5 条合成蛋白 + 2 个迷你 HMM
  跑通 06b→07→08 的输入/输出契约（不依赖 GTDB 全库），验证 fail-closed 与列 schema。
  服务器执行：`bash pipeline/tests/test_smoke.sh`
- CI（GitHub Actions）：`.github/workflows/ci.yml` 静态检查所有 Python、bash `-n`、
  Python 冒烟测试与审计门禁回归测试；不依赖 HMMER。
- 发布门槛（Zenodo DOI 前）：同一次运行的 `run_manifest.json` 无 failed step、Git/源码包和 GTDB
  校验和已固化、环境与完整命令已登记、冒烟测试和审计测试通过，且 `tree_manifest.tsv` 无
  `EMPTY_OR_PARSE_FAIL`、`stale_input` 或 `input_not_registered` 状态。当前清单仍有 OH
  `stale_input`（2 条）和 ePhaZ CD-HIT `input_not_registered`（1 条），不满足发布门槛。
