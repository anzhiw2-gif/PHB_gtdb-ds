# STATUS — PHB_gtdb-ds 项目单一状态页

> 本页是项目**唯一的状态事实源**（单一状态页）。README、结果报告、审核报告等
> 只引用本页，不各自重复状态判断；口径冲突以本页为准。
> 更新日期：2026-08-24 ｜ 仓库：github.com/anzhiw2-gif/PHB_gtdb-ds

---

## 2026-08-24 本地运行隔离改造

- `run_pipeline.sh` 默认使用 `runs/<run_id>/`，拒绝覆盖已有运行目录；支持 `--run-id`、
  `--run-root`/`--run-dir` 和显式 `--legacy-root-results`。
- 每次启动生成 `run_context.env` 与 `input_contract.json`，区分 `verified`、`pending`、
  `missing` 输入状态。
- 最终 manifest 绑定 `run_id`、`run_root`、input contract、源码 bundle、命令和哈希。
- 本轮仅完成本地代码、测试和文档；未连接或重算 T141，未安装软件、删除结果、提交或推送。

## 0.1 2026-08-24 audit correction and authority map

This section supersedes conflicting numbers or provenance claims in older report
paragraphs. The canonical Figure 3 source table records `Pseudomonadota=26,850`
genomes (`results/figures/scheme_a/source_data/figure3_phylum_totals.tsv`). Any
older `26,855` value is stale and must not be reused. Figure 5 uses an
all-hit-neighborhood denominator (`candidate_loci`); it is not a tier1 count and
not the archaeal patatin subset. The patatin result remains 1,372 loci / 620
genomes, while the 112,926 tier1 records are a broad patatin-fold candidate set.

Authority is currently split and must be named explicitly: the local working
tree is an uncommitted audit-repair workspace; GitHub `main` remains at
`9a7d02d` and is not this state; the T141 dated deploy is the only server-side
source eligible for a future rerun. The audited T141 main manifest currently
hashes to `e991c3bd10a48f8faf9c450f0c17a5a3fb1f0315c256018f0772c5f64f71b2a3`,
which differs from the older `c4e8...` value printed elsewhere in this file.
The existing server manifest also lacks the new strict source-bundle and full
command provenance; a new dated run is required before claiming that contract.

The tree gate is unchanged: OH remains `stale_input` (recorded 1,465 versus
current 1,429), ePhaZ CD-HIT remains `input_not_registered`, and full ePhaZ/iPhaZ
trees plus HGT remain incomplete. These historical files are retained for
forensics and are not current scheme-A tree evidence.

## 1. 结论边界（措辞规范，全局适用）

本项目产出的是**基于 HMM 同源性 + 基序 + 邻域证据的候选同源蛋白**。SignalP 细分已完成，
但尚未经酶活、遗传或表型实验验证。因此：

- ✅ 允许表述：**“候选同源蛋白 / 功能潜力 / 推定（putative）”**。
- ❌ 禁止表述：**“PHB 降解基因确实存在”“功能基因”“首次证实/实证发现”**。

具体到古菌：只能说“古菌谱系检出 PHB 降解相关**候选同源蛋白**”，并附
“patatin 为广谱磷脂酶折叠、PhaJ 才是动员主通路”的 caveat（见 §5）。

---

## 2. 项目状态

- **阶段**：方案 A 候选集、生态分布与位点级基因邻域已在 T141 的隔离运行目录
  `/home/data/haoyu/PHB_gtdb-ds/runs/20260821_schemeA_03` 完成并固化；建树仍按用户决定暂停。
- **可复算的核心统计**：四家族去重值 **44,814 基因组 / 22.416%**（方案 A 的 OH `min-cov 0.6` 口径；已从本地最终 tier 表复算）。
- **运行验收**：T141 隔离运行的 `run_manifest.json` 含 11 个成功步骤；18 个输入和 18 个输出均非空且
  SHA-256 校验一致。当前文件实际 SHA-256 为
  `e991c3bd10a48f8faf9c450f0c17a5a3fb1f0315c256018f0772c5f64f71b2a3`；此前记录的
  `c4e8c73bc81e3085cfbd67a0f2ef5f153fbfc6f35a8508239b582cb06c833a7f` 已标记为过时，不能继续作为当前 manifest 哈希。
- **服务器**：T141（10.16.1.141），80 核 / 1TB 内存；当前无运行中任务。
- **线程上限**：计算任务合计 **≤70**（留 10 核余量）。已同步到
  `pipeline/config/params.yaml`、`05/06/run_pipeline` 脚本默认值。

**版本与运行 provenance（2026-08-24 核对）**：本地 `D:\PHB_gtdb-ds` 的 `HEAD` 与 GitHub
`main` 均为 `9a7d02d`，但本地工作树含未提交修改和未跟踪审计/结果文件；因此本地工作树不等同于
GitHub 已发布快照。T141 项目根目录没有 `.git`，主 manifest 的 `git_commit` 为 `null`，不能把该
运行宣称为某个 Git commit 的产物。dated forensic repair 使用
`/home/data/haoyu/PHB_gtdb-ds/deploy/20260821_schemeA/scripts/11_clusters.py`
（SHA-256 `d1d907f34a5c1fbe17aee538a3ca087de8e809f5c03286b6384a22347d994d2c`）；服务器项目根
`scripts/11_clusters.py` 是另一版本。以后只允许从 dated `deploy/<run_id>/scripts/` 或明确绑定的
Git commit 执行，不能直接运行服务器根目录旧脚本。

---

## 3. 数据流契约（路径 + 顺序 + fail-closed）

**路径契约（两套结构，务必区分）**：

| 环境 | 脚本目录 | 工作区/仓库根推导 | 种子/配置位置 |
|------|---------|------------------|--------------|
| 本地 git | `pipeline/scripts/` | `ROOT = SCRIPT_DIR/../..` | `pipeline/seeds/`、`pipeline/config/` |
| 服务器 T141 | `scripts/`（扁平） | `ROOT = SCRIPT_DIR/..`（或硬编码 `/home/data/haoyu/PHB_gtdb-ds`） | `data/seeds/`、`config/` |

- **本地 git 脚本**（`pipeline/scripts/`）用 `SCRIPT_DIR/../..` = 仓库根；`data/`、`results/` 在仓库根。
- **服务器执行脚本**（scp 到 `scripts/` 跑）用 `SCRIPT_DIR/..` = 工作区根，或硬编码
  `/home/data/haoyu/PHB_gtdb-ds`（如 `signalp_ePhaZ.sh`、`cdhit_tree.sh`）。**禁止**把
  `pipeline/` 当项目根，也**禁止**把本地 `../..` 约定套到服务器扁平结构上。
- 服务器无 `pipeline/` 目录；`data/seeds/` 对应本地 `pipeline/seeds/`。

**脚本顺序（主编排 run_pipeline.sh 强制）**：

```
05_predict_proteins.sh   → data/proteins/shards/         （Pyrodigal meta）
06a_filter_shards.sh     → data/proteins/shards_filt/    （剔除 >100K aa）
06_screen.sh (+06b)      → data/screen/hits_all.tsv      （hmmsearch + cov 列）
07_process_hits.py       → hits_filtered.tsv / genome_hits.tsv / family_seqs
08_validate.py           → validation.tsv + *_validated.faa
08c_tier_rescore.sh      → tiers/{fam}_tier{1,2}.faa
09/09b/09g/09i           → 系统发育 + 树登记
10_distribution.py       → results/tables/* 生态/分类分布
11_clusters.py           → cluster_context.tsv / cluster_summary.tsv
```

**fail-closed 约定**：任何一步非零退出立即中止（`run_pipeline.sh` 不再“记警告后
继续”）；`06_screen.sh` 的 hmmsearch/parallel、`08c` 的 hmmsearch 均拒绝吞失败。

`06a_filter_shards.sh`、`06_screen.sh` 与 `08c_tier_rescore.sh` 均先在带 PID 的临时构建目录完成
全量任务和 provenance/计数校验，再归档旧结果并原子发布；缺失输入、HMM、任务输出或计数不一致时不会
覆盖已有结果。

**关键输入契约（11 步）**：`11_clusters.py` 需要 locus 级输入，主编排传
`data/screen/hits_filtered.tsv`（含 locus），**不得传 `genome_hits.tsv`**（仅
基因组×家族矩阵，无 locus）。脚本内已 fail-closed 校验必需列。

**覆盖度（min-cov）**：`06b_aggregate_hits.py` 从 domtblout 计算每命中的 HMM
覆盖度并写入 `hits_all.tsv` 的 `cov` 列；`07_process_hits.py` 的 `--min-cov`
现在真正生效。**默认 0.0（不按覆盖度过滤，与已提交结果的旧口径一致）**；
非零阈值需用固定正负对照集校准后再启用（见 §6 待办）。

---

## 4. 核心结果数字（最终修正值）

| 家族 | tier1 序列（方案 A 重建后：OH 加 min-cov 0.6） |
|------|:---:|
| ePhaZ | 38,692 |
| iPhaZ | 32,926 |
| OH | 1,429（原 1,465，−36） |
| ArchPhaZ_hydrolase | 1,292 |

- **去重后核心解聚酶：44,814 基因组 / 22.416%**（不含广谱 patatin；原 44,821/22.419%，
  OH 加 min-cov 0.6 排除尼龙水解酶后净 −7 基因组）。
- ✅ **数据同步已修复**（2026-08-19）：本地 `data/screen/tiers/*.faa` 已从服务器
  同步为 Glu 修正后版本（ePhaZ 38,692 / iPhaZ 32,926 / OH 1,465 / hydrolase 1,292 /
  patatin 112,926 条广谱候选蛋白记录；古菌 patatin 讨论子集另为 1,372 个 locus / 620 个 genome）。同步期间曾误写 tier1.faa，已用 `rebuild_tier1_faa.py` 从
  `*_tier1.ids` + `*_validated.faa` 幂等重建（见 CHANGELOG §五）。
- ✅ **方案 A 已执行**（2026-08-19，`rerun_candidates.sh`）：06b→07(家族 min-cov)→07b→08→
  08c→09a 全链重跑，OH 加 `--min-cov 0.6`。tier1 变化仅 OH：1,465→1,429 序列 /
  1,444→1,410 基因组；四家族去重 44,821→44,814。

---

## 5. 古菌 patatin 结论（审慎口径）

- patatin 折叠是广谱磷脂酶/酯酶结构域，**不是 PHB 解聚酶专属**。
- 古菌 patatin 候选集：**1,372 位点 / 620 基因组** → 表述为“**古菌 patatin 候选集**”。
- 方案 A 的位点级邻域已以最终 `hits_filtered.tsv` 重跑并合并 80 个批次，输出
  `cluster_context.tsv`（212,240 行）、`cluster_locus_audit.tsv`（1,036,243 个审计位点）和
  新 schema 的 `cluster_summary.tsv`（30 个家族×标记组合）。这些文件目前保存在上述 T141 隔离运行目录。
- 新 schema 明确区分 `marker_hits`、`supporting_loci`、`supporting_genomes`。本次实际可用 marker
  为 `PhaC,PhaE,PhaJ,BdhA,phasin,PHA_gran_rgn`；服务器没有 `PhaA/PhaB/PhaP/PhaR` HMM，不能将其表述为完整 marker panel。
- Figure 5 的 `results/figures/scheme_a/source_data/figure5_neighborhood_rate.tsv` 使用**全命中邻域层**
  的 `candidate_loci` 作为 denominator（例如 patatin 行为 315,065 个 candidate loci），不是
  tier1 核心四家族，也不是古菌 patatin 的 1,372 loci / 620 genomes 子集。Figure 5 的 support rate
  只能描述该 all-hit 邻域层，不能回写成 patatin 子集的比例或已证明的代谢通路。
- 生物学 caveat：PhaZh1 型**体内角色有限**（敲除不影响动员），**PhaJ（烯酰-CoA
  水合酶）才是古菌动员主通路**。
- 历史 `cluster_summary.tsv` 保留以供追溯，但不再作为最终结论来源；新 schema 的输出明确区分
  `marker_hits`（出现次数）、`supporting_loci`（唯一命中位点）和 `supporting_genomes`（唯一基因组）。

---

## 6. 系统发育树状态（以 tree_manifest.tsv 为准）

审计清单：`results/trees_tier1/tree_manifest.tsv`（由 `09i_tree_manifest.py`
生成，含 kind/工具/叶数/树 SHA-256、记录时输入与当前输入 SHA-256、抽样名单和状态）。
该清单会保留历史输入哈希，防止输入变更后把旧树误标为当前树。

| 家族 | kind | 工具 | 叶数 | 说明 |
|------|------|------|:---:|------|
| OH | iqtree | IQ-TREE2 | 1,465 | **旧输入全量树，`stale_input`**；当前 OH tier1 为 1,429，必须重建后才能用于方案 A |
| ArchPhaZ_hydrolase | iqtree | IQ-TREE2 | 1,292 | **全量树**（1292 ≤ 2000，未抽样） |
| ePhaZ | cdhit | IQ-TREE2 | 10,648 | CD-HIT 去冗余后建树（非全量） |
| ePhaZ | iqtree | IQ-TREE2 | 1,942 | 抽样树（seed=42, N=2000） |
| iPhaZ | iqtree | IQ-TREE2 | 2,000 | 抽样树（seed=42, N=2000） |
| 4 家族 | fasttree | FastTree | 1,000 | 抽样树（seed=42, N=1000） |

- ePhaZ/iPhaZ 的 IQ-TREE 树仍为**抽样树**（非全量）；ePhaZ 另有 CD-HIT 去冗余
  （10,648 叶）近似树。**全量树（ePhaZ 38,692 / iPhaZ 32,926）规模过大，未建**。
- OH 的两棵树（FastTree/IQ-TREE2）记录输入为 1,465 条，而当前 tier1 为 1,429 条，均为 `stale_input`；不能用于方案 A 的推断。ePhaZ 的 CD-HIT 树缺少其服务器端去冗余 FASTA 的登记，状态为 `input_not_registered`。
- **建树已暂停（用户决定，2026-08-19）**：完整树（ePhaZ/iPhaZ 大族）与 HGT 检测暂不进行，现有树仅作记录，不继续扩展。
- 当前 `tree_manifest.tsv` 汇总为 `ok=6`、`stale_input=2`、`input_not_registered=1`；其中 OH 两棵树
  属于 `stale_input`，ePhaZ CD-HIT 树属于 `input_not_registered`。这三类状态都不能作为当前
  方案 A 的正式树证据。

---

## 7. 待办（按优先级）

**高（投稿前）**
- [x] 从服务器同步 Glu 修正后的 tier1.faa（已完成，2026-08-19；数字与报告一致）。
- [x] 用固定正负对照集校准 + 方案 A 重建候选集（OH 加 min-cov 0.6；见 §9；新四家族去重 44,814/22.416%）。
- [x] ePhaZ SignalP 细分（已完成：**有信号肽 56.5%**（Sec/SPI 42.8% + Lipo/SPII 13.4% + Tat 0.3%），无信号肽 OTHER 43.5%；N 端截断仅 ~2.5%，非主因，见 §10）。
- [x] 在 T141 用最终方案 A 输入重跑生态分布和位点级基因邻域，生成新 schema 的 `cluster_context.tsv`/`cluster_summary.tsv`（80 批全部成功，2026-08-23）。
- [ ] 重建 OH 树并登记实际输入；登记 ePhaZ CD-HIT 去冗余 FASTA 的来源与 SHA-256。
- [x] 生成 T141 隔离运行的 `run_manifest.json`，将输入、步骤 exit code 和输出哈希绑定到同一次运行（2026-08-23；见 §10）。
- [ ] 将源码包/明确 Git commit、完整环境锁定、GTDB taxonomy/metadata/tree 与 HMM 哈希补入同一 provenance 合同；当前 manifest 的 `git_commit=null`。

**低 / 发布前**
- [x] 种子 manifest 重建（reviewed 55/23 已修正 + evidence + split，见 §8）。
- [ ] 整理成论文；发布 HMM profiles + 命中表（GitHub Release + Zenodo DOI）。

> **建树已暂停（用户决定）**：iPhaZ CD-HIT 建树已停止，后续建树任务暂缓。
> 历史树保留于 `results/trees_tier1/` 以供追溯；OH 树为 `stale_input`，ePhaZ CD-HIT 树为
> `input_not_registered`，不能作为最终方案 A 的树证据。

### 已执行的服务器脚本

**方案 A 受审计重算记录（2026-08-21 至 2026-08-23）**：隔离运行目录完成 06b → 07
（OH `min-cov 0.6`）→ 07b → 08 → 08c → 09a → 10 → 11。cluster 首次并行运行中
`batch_001` 的空 FAA 已作为 `invalid_fasta_record` 审计并恢复；原失败 JSONL 保留为
`results/run_manifest.before_cluster_recovery.jsonl`，最终 manifest 只接受 11 个成功步骤。

**论文 + Zenodo（待办⑦，人工写作/发布）**：交付物清单见 `docs/reproducibility.md`
§6-8；HMM profiles + tier1 命中表走 GitHub Release + Zenodo DOI。

---

## 8. 种子 manifest 与训练/验证划分

- `02_collect_seeds.py` 已修复 `reviewed` 判定（不再用 `"reviewed" in entryType`，
  改为 `reviewed` 字段（返回 entryType）/ 前缀匹配），并新增：`--min-reviewed-only`、
  `retrieval_date`、`evidence`（PMID/DOI）、`split`（train/validation，80/20 按
  accession md5 取模）。
- ✅ **manifest 已重建**（2026-08-19，`rebuild_seeds_manifest.py` 用 UniProt 实时查询）：
  **reviewed=55 / unreviewed=23**（旧版误标为全 true），train/val=64/14，每条含
  evidence(PMID/DOI) 与 retrieval_date。78 条中的 18 条同源物标注见
  `pipeline/seeds/seeds_annotation.md`。

---

## 9. HMM 阈值校准（初步，固定正负对照集）

对照集（`pipeline/seeds/controls/`，`build_control_sets.py` 生成）：正对照 55 条
（真 PHB/PHA 解聚酶）、负对照 18 条（13 真核 BDH + 5 尼龙水解酶）。
方法：`calibrate_hmms.py` 用 curated HMM 对正负对照 hmmsearch，在阈值网格算
TP/FP/F1/MCC。结果：`results/tables/calibration_summary.tsv` + `calibration_hits.tsv`。

| 家族 | E-value 分离 | 关键发现 |
|------|:---:|------|
| ePhaZ | ✅ | TP=10/FP=0/FN=0，E=1e-2 即完美分离（默认 1e-5 更保守，无需改） |
| iPhaZ | ✅ | TP=17/FP=0/FN=0，同上 |
| OH | ⚠️ | TP=28/FP=5/FN=0；**5 条尼龙水解酶（nyl）被 OH HMM 强命中（E=1e-49~74），E-value 无法区分**；但 **cov 可区分**：nyl cov=0.31~0.53，真 OH cov≈0.99（唯一例外 Q4W8C9=phaZc, cov=0.203） |

**结论**：OH 家族建议启用 `--min-cov ≈ 0.6` 以排除尼龙水解酶假阳性（代价：误丢
1 条 Q4W8C9，其为 C. necator phaZc，与主流 3HB 寡聚体水解酶结构差异大）。
ePhaZ/iPhaZ 维持 E=1e-5、min-cov=0 不变。**此校准基于 78 条种子（小样本），
仅作阈值参考，重建候选集需在服务器用最终阈值全库重筛。**
（注：校准结果同步于服务器 results/tables/，本地副本见同名文件。）

### 9.1 SignalP 分析（ePhaZ N 端截断核查，`analyze_signalp.py`）

- SignalP6 类型分布（38,692 条）：**有信号肽 56.5%** = Sec/SPI 16,545(42.8%) +
  脂蛋白 Lipo/SPII 5,170(13.4%) + Tat 122(0.3%) + TatLipo 19；**无信号肽 OTHER
  16,836(43.5%)**。（注：早期只数 Sec/SPI 得 42.8%，漏了 Lipo/SPII 13.4%，已修正。）
- N 端截断核查（无信号肽 OTHER 组）：不以 Met 开头仅 2.8%（vs 有信号肽组 0.3%）、
  无内部终止符、长度中位数 364 vs 344 aa 相近 → **截断不是"无信号肽"的主因**。
- 结论：约 43.5% 无信号肽主要是"催化域与 ePhaZ 同源但非分泌"的变体（胞内型/
  非 Sec 分泌/革兰氏阳性胞外酶等），"胞外解聚酶"标签需审慎。
- **高置信胞外子集已导出**：`data/screen/tiers/ePhaZ_tier1_signalpeptide.faa`
  （21,856 条，header 带 SignalP 类型）+ `results/tables/ePhaZ_signalp_subset.tsv`，
  并按类型拆为 `ePhaZ_tier1_signalpeptide_{SP,LIPO,TAT,TATLIPO}.faa`。
- **门分布**（`phylum_dist.py`）：有信号肽 21,856 条 → Actinomycetota 30.3% /
  Pseudomonadota 25.5% / Bacteroidota 15.3%；SP-only 16,545 条、LIPO-only 5,170 条
  分别见 `results/tables/ePhaZ_SP_phylum.tsv` / `ePhaZ_LIPO_phylum.tsv`。

## 10. 复现性与审计

- 运行清单：T141 隔离运行目录的 `results/run_manifest.json`（由 `run_manifest.py`
  生成，含每步 exit code、起止时间、输入/输出 SHA-256）；其 SHA-256 见 §2。运行目录不覆盖仓库
  现有 `results/`，也不应把服务器机器本地产物误标为已提交的 Git 文件。
- 命令清单与数据溯源：`docs/reproducibility.md`、`pipeline/README_HPC.md`。
- 小样本端到端测试：`pipeline/tests/test_smoke.sh`（见 `docs/reproducibility.md` §7）。
