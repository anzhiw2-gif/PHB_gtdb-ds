# PHB_gtdb-ds 项目综合审计报告

> **历史快照（2026-09-01）**：本文记录当日 run-12 尚在运行时的状态，不能代表当前运行结论。当前权威状态见 [CURRENT_STATUS_20260902.md](CURRENT_STATUS_20260902.md)；原始证据保留用于追溯。

**审计截止时间：** 2026-09-01（Asia/Shanghai）
**项目：** `PHB_gtdb-ds`  ︱ **远端主机：** `<SERVER_USER>@<SERVER_HOST>`  ︱ **GitHub：** `anzhiw2-gif/PHB_gtdb-ds`

> 本报告是基于当前本地工作树、Git 引用、T141 dated deploy/run、输入契约、模型注册表和已有科研状态文档形成的审计快照。大规模 GTDB 原始蛋白数据保留在 T141，不复制进 Git；其完整性以分片 SHA-256 清单为证据。报告中的 HMM、domain、SignalP、邻域和树结果均表示候选同源或功能潜力，不等同于已验证 PHB 降解表型。

## 1. 审计结论摘要

| 审计面 | 当前状态 | 结论 |
|---|---|---|
| 代码与测试 | `changed-and-verified` | 本地测试 `118 passed, 1 skipped`；`compileall` 和 `git diff --check` 通过。 |
| 模型注册 | `verified-current` | 正式注册表有 10 个冻结模型，阈值为 `E=1e-5`；`OH` 使用 `min_cov=0.6`。 |
| ePhaZ 模型 | `verified-current` | 已拆分为 `ePhaZ_curated_core` 与 `ePhaZ_broad_discovery`；broad 仅用于探索，不作 PHB 特异性判定器。 |
| T141 预检 | `verified-current` | 100 个蛋白分片、GTDB taxonomy/metadata/tree 和 10 个 HMM 均已哈希锁定。 |
| T141 正式扫描 | `running` | 新 run `20260831_formal_frozen_scan_12` 正在运行，尚未生成最终命中表或完成 manifest。 |
| 旧正式扫描 | `failed-retained` | run 11 因 HMMER 超过 `100000 aa` 的工具限制失败，已终止但完整保留失败证据。 |
| 生物学结论 | `candidate-only` | 当前结果只能支持候选同源/功能潜力，不能宣称已发现或验证 PHB 降解基因。 |
| GitHub 同步 | `out-of-sync` | 本地 `HEAD=e0f073f`，GitHub `origin/main=aa47d50`；本地领先 4 个提交且存在未提交/未跟踪科研改动。 |
| 发布状态 | `pending` | 正式扫描未完成，README/STATUS 仍含旧状态描述；不应发布最终统计或 DOI 包。 |

**总判断：** 项目治理和可审计性已经明显改善，冻结模型和正式扫描入口具备可追溯基础；但当前仍处于“正式计算运行中、结果未验收、文档未完全同步”的阶段，不具备最终生物学结论或正式发布条件。

## 2. 权威关系与复现边界

### 2.1 当前权威来源

1. **本地工作树：** `<LOCAL_WORKSPACE>`。当前包含尚未提交的科研改动和运行归档，不能视为 GitHub 已发布版本。
2. **GitHub：** `origin/main` 指向 `aa47d50d045af3249f38d758bb9cd1fd68d4a384`，只代表已推送治理基线，不包含当前本地全部科研改动。
3. **服务器执行源：** 只能使用 dated `deploy/<run_id>/`；当前正式扫描使用 `deploy/20260831_formal_frozen_scan_12/`。
4. **服务器运行证据：** `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260831_formal_frozen_scan_12/`。原始 GTDB 蛋白分片不进入 Git。

### 2.2 运行目录约束

每次运行使用唯一 `runs/<run_id>/`，并保存 `logs/`、`inputs/`、`results/` 和 `input_contract.json`。不得覆盖历史 `results/` 或服务器历史 run。正式扫描不复用预检目录：

- 预检：`20260831_formal_frozen_preflight_10`，状态 `planned_not_run`；
- 正式扫描：`20260831_formal_frozen_scan_12`，状态 `running`；
- 失败证据：`20260831_formal_frozen_scan_11`，状态 `failed-retained`。

## 3. 代码、配置与测试

### 3.1 代码结构

- 主流程脚本：`pipeline/scripts/01-11` 及 `run_pipeline.sh`；
- 模型与 ePhaZ 治理脚本：`split_ephaz_seeds.py`、`clean_seed_library.py`、`calibrate_ephaz_layers.py`、`calibrate_ephaz_external_panel.py`、`audit_ephaz_iphaZ_competition.py`；
- 正式扫描：`formal_frozen_screen.sh`；
- 超长序列处理：`filter_hmmsearch_shard.py`；
- provenance：`run_context.py`、`run_manifest.py`；
- 测试：`pipeline/tests/`，当前 25 个测试文件。

### 3.2 环境与主要参数

- Conda 环境：`phb_gtdb`；Python 3.12；HMMER 3.4；DIAMOND 2.2.1；Pyrodigal 3.7.1；Biopython/Pandas/NumPy 等；
- GTDB：Release `R232`，GTDB taxonomy、metadata、tree 均在输入契约中记录；
- 预测：Pyrodigal meta 模式，配置线程上限 `70`，每个分片约 2,000 个基因组；
- 正式筛选：HMMER `E=1e-5`；OH 额外要求 HMM coverage `>=0.6`；
- HMMER 超长目标限制：超过 `100000 aa` 的记录不送入 HMMER，写入 `overlength_exclusions.tsv`，作为工具限制证据而不是阴性。

### 3.3 测试证据

```text
python -m pytest pipeline/tests -q
118 passed, 1 skipped, 6 subtests passed

python -m compileall -q pipeline/scripts pipeline/tests
PASS

git diff --check
PASS
```

测试涵盖运行隔离、输入契约、模型注册、ePhaZ 分层、外部校准、SignalP/iPhaZ 竞争审计、结构复核和聚合命中处理。正式扫描脚本的 Bash 语法已在 T141 dated deploy 上通过 `bash -n`。

## 4. 冻结模型注册表

权威注册表：`pipeline/config/formal_scan_models.tsv`；服务器运行副本：`runs/20260831_formal_frozen_scan_12/inputs/formal_scan_models.tsv`。

| 模型 | 来源 | E-value | min coverage | 报告组 | HMM SHA-256 |
|---|---|---:|---:|---|---|
| `ePhaZ_curated_core` | `split_run` | `1e-5` | 0.0 | ePhaZ | `96668851b42fba67d5fa8903e4ac25527769abad8b8e8d7366b3a62db1a9043f` |
| `ePhaZ_broad_discovery` | `split_run` | `1e-5` | 0.0 | ePhaZ | `3f981c5a12b7a7ec9389eab8e329f03e220875bdbce94f95849821e32dc127bc` |
| `iPhaZ` | `frozen_data_root` | `1e-5` | 0.0 | iPhaZ | `4a1d15d37e8f55bea89fbc6bc33b1e1c45b166c80aad3976c73c087937f7d311` |
| `OH` | `seedclean_run` | `1e-5` | 0.6 | OH | `037fc34b135b720e483710a397bbd120075e541d69b0e3352986e07a23a6875a` |
| `BdhA` | `seedclean_run` | `1e-5` | 0.0 | BdhA | `cb59917b36e0febf976e853071a85d5fe7e0229da81718d00b4f214fd7504aa7` |
| `ArchPhaZ_patatin` | `seedclean_run` | `1e-5` | 0.0 | ArchPhaZ_patatin | `c25085668d22ca39e654179a6f9d59d5e69dc8a6fd3083bdfcd5a63e7c274483` |
| `ArchPhaZ_hydrolase` | `seedclean_run` | `1e-5` | 0.0 | ArchPhaZ_hydrolase | `d667072ec2549f405c62930c57e844dafe418ec96a5edf23be28123366756ce8` |
| `PhaJ` | `seedclean_run` | `1e-5` | 0.0 | PhaJ | `69e0e9e604df2af1785dceba14130c5d7dda538c397c8dd07dc3ef22747bbc50` |
| `phasin` | `seedclean_run` | `1e-5` | 0.0 | phasin | `34ee7541a263de5c4508009efa1b0838ed937e6a6ae66c9b05bc413a6cfef826` |
| `PhaC` | `seedclean_run` | `1e-5` | 0.0 | PhaC | `7ce60b555be67dffbacbd01679e662e9b822d9df28c0807e8a12a5d9f773c973` |

### 4.1 ePhaZ 分层

- `ePhaZ_curated_core`：实验支持、完整长度、典型结构域架构；当前 4 条核心实验阳性全部命中，18 条负对照未命中。
- `ePhaZ_broad_discovery`：远缘、注释型及全部 61 条 `<200 aa` 序列；不因长度直接删除，短序列逐条记录在审查表中。
- broad 历史恢复文件：
  - HMM SHA-256：`3f981c5a12b7a7ec9389eab8e329f03e220875bdbce94f95849821e32dc127bc`；
  - alignment SHA-256：`13788bf99b3a07a8045718ebe884dd1f1e5e6c7e78196556559ec225e72951bc`；
  - seed FASTA SHA-256：`8acee9dec14a15dfc5e8d5f485c501e7b41eb461afa54af38f921151b88f8d08`。

## 5. 校准与独立面板

### 5.1 外部敏感性分析

证据：`runs/20260830_ephaz_external_panel_sensitivity_03/results/calibration/panel_sensitivity.tsv`。

- ePhaZ 模型在 `E=1e-5`、coverage 0.0/0.4 时 PHB 实验阳性召回为 `5/5`；coverage 0.6 降至 `2/5`，coverage 0.8 降至 `0/5`；
- `Q84C08` 是实验支持的 MCL-PHA depolymerase，明确不水解 PHB，但被两个 ePhaZ 模型命中；
- 因此 ePhaZ 仍冻结为 `E=1e-5`、coverage 0.0，coverage 只作下游分层证据，不作为通用硬过滤器；
- PHB、MCL-PHA、胞内 PHB、MCL-PHA 非 PHB、注释型和 challenge 分母必须分开报告，不允许合并成单一“阳性率”。

### 5.2 ePhaZ/iPhaZ 无信号肽竞争

证据：`runs/20260830_ephaz_competition_audit_05/results/competition/no_signal_competition.tsv`。

- SignalP `OTHER`：16,836 条；
- 当前竞争审计分层：ePhaZ-like 16,809、iPhaZ-like 9、ambiguous 18；
- 这些数字描述模型竞争关系，不是 PHB 表型分类；无信号肽不能直接解释为胞内或胞外功能。

### 5.3 570 条人工结构复核

证据：`runs/20260829_ephaz_ambiguous_structure_review_02/results/structural_review/ambiguous_structural_review.tsv`。

| 类别 | 数量 |
|---|---:|
| 样本总数 | 570 |
| 以 M 开始且无内部终止 | 536 |
| 可能 N 端截断 | 34 |
| `iPhaZ_consistent` | 511 |
| `partial_ePhaZ_signal` | 47 |
| `mixed_cross_family` | 1 |
| `insufficient_structural_support` | 11 |
| `provisional_iPhaZ_challenge` | 489 |
| `pending_manual` | 81 |

`pending_manual`、`partial_ePhaZ_signal` 和 `mixed_cross_family` 不能进入 ePhaZ 正对照分母，也不能直接支持 PHB 降解结论。

### 5.4 独立实验与负对照面板

证据：`runs/20260830_ephaz_independent_positive_negative_01/`。

- PHB 实验阳性：5 条；
- MCL-PHA 实验阳性：1 条，独立统计；
- 胞内 PHB 非 ePhaZ 负对照：3 条；
- MCL-PHA 非 PHB 负对照：1 条；
- 注释型近邻负对照：5 条；
- fragment/incomplete challenge：2 条，不进入正式阴性分母；
- 6 条 `iPhaZ-like` challenge 已从 ePhaZ broad 正对照分母移除，保留为跨家族挑战集；不升级为 core，不加入 bridge。

## 6. T141 正式扫描状态

### 6.1 预检 run

证据：`runs/20260831_formal_frozen_preflight_10/`。

- manifest：`results/preflight_manifest.json`；
- 状态：`planned_not_run`；
- 分片：100；
- 模型：10；
- `hmmsearch_outputs_created=false`；
- GTDB taxonomy SHA-256：`12de91fa3b1267cc2aeec843d6b9382aa98175b1e1d64804095eb09db80b3c4f`；
- GTDB metadata SHA-256：`1650d3164666e5839c20ee15d82511909a5bbd5269035a7b564b9048dd777893`；
- GTDB tree SHA-256：`9034e52f25ed0caead4e2153def93ee80ba0f513a0a2f65030127839a26fe02d`；
- 注册表 SHA-256：`8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`。

### 6.2 当前正式 run

**Run：** `20260831_formal_frozen_scan_12`
**Deploy：** `deploy/20260831_formal_frozen_scan_12/`
**服务器目录：** `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260831_formal_frozen_scan_12/`
**查询时间：** 2026-09-01 00:11（T141 本地时间）

当前只读快照：

| 指标 | 当前值 |
|---|---:|
| 模型×分片任务总数 | 1,000 |
| 已生成 `.tbl` | 251 |
| 已生成 `.dom` | 251 |
| 已过滤分片 | 25 / 100 |
| 失败任务 | 0 |
| 超长序列排除 | 3 条（ledger 含表头共 4 行） |
| `hits_all.tsv` | 尚未生成 |
| `scan_manifest.json` | 尚未生成 |

主进程 PID 为 `250321`，当前子进程为 `ePhaZ_curated_core × shard_0026`，CPU 约 153%，说明正在实际执行 HMMER，而非停死。单线程模型×分片串行方式导致每个约 2 GB 分片耗时较长。

### 6.3 失败 run 11

证据：`docs/T141_20260831_formal_frozen_scan_11_failed.md`。

HMMER 3.4 对超过 `100000 aa` 的目标序列直接中止，run 11 记录 5 个失败任务后终止。该失败是工具输入限制，不是生物学阴性；相关 `stderr`、`failed_tasks.tsv` 和部分 `hmmsearch.build` 均应保留用于取证。

## 7. 历史核心结果与解释边界

以下数字来自此前已完成的方案 A / P06-P07 结果，不应与当前尚未完成的 run 12 混用：

- 四个核心家族去重后候选基因组数：44,814，约 22.416%；
- tier1 代表性统计：ePhaZ 38,692，iPhaZ 32,926，OH 1,429，ArchPhaZ_hydrolase 1,292；
- SignalP ePhaZ tier1：有信号肽 21,856（56.5%），OTHER 16,836（43.5%）；
- SignalP 有信号肽组成：Sec/SPI 16,545、Lipo/SPII 5,170、Tat 122、TatLipo 19；
- 古菌核心解聚酶候选分布：Halobacteriota 112、Thermoproteota 32、Thermoplasmatota 154；
- 古菌 patatin 广谱折叠候选集：1,372 loci / 620 genomes；该集合不能等同于 PHB 专一 depolymerase。

关键解释边界：

1. `PhaJ` 是古菌 PHB 动员的主要已知通路线索；patatin 是广谱脂质酶/磷脂酶折叠，不能直接称为 PHB 特异性解聚酶。
2. Figure 5 使用 all-hit `candidate_loci` 分母，不是 tier1，也不是 1,372 loci 的古菌 patatin 分母。
3. 旧文档中的 `Pseudomonadota=26,855` 已被 `26,850` 取代，当前 Figure 3 源表为准。
4. 所有 HMM 命中、SignalP 分类、结构域覆盖、邻域共现和系统树均需经过分层审计后才能支持候选排序，不能直接转写成表型。

## 8. 系统发育与邻域分析状态

依据 `results/trees_tier1/tree_manifest.tsv` 和 `docs/STATUS.md` 的最后核对：

- `ArchPhaZ_hydrolase` IQ-TREE2 全量树：1,292 条，可作为当前输入记录；
- OH 树使用旧 1,465 条输入，而当前 OH tier1 为 1,429，状态 `stale_input`；
- ePhaZ CD-HIT 树状态 `input_not_registered`；
- ePhaZ/iPhaZ 全量树和 HGT 分析仍未完成；
- 已有抽样树和历史树只能作为追溯材料，不能宣称为当前全量系统发育证据；
- 基因邻域分析的正式输入必须使用带 locus 的 `hits_filtered.tsv`，不能误用只含 genome/family 的 `genome_hits.tsv`。

## 9. Git、GitHub 与工作树审计

### 9.1 版本关系

| 位置 | 当前引用 | 状态 |
|---|---|---|
| 本地分支 `main` | `e0f073f7d81212af34928eac85317addd8d68634` | 本地最新，含科研改动归档提交 |
| GitHub `origin/main` | `aa47d50d045af3249f38d758bb9cd1fd68d4a384` | 落后本地 4 个提交 |
| T141 正式 deploy | `deploy/20260831_formal_frozen_scan_12/` | 服务器执行权威 |
| T141 项目根 | 无 `.git` | 不能声明为某个 Git commit 的直接产物 |

### 9.2 工作树残留

当前工作树仍有未提交或未跟踪内容，包括科研脚本、测试、运行归档、图像、SignalP 结果和 deploy 目录。另有未跟踪文件 `--amino`，应在发布前确认来源和用途；不得在未确认前删除。当前状态不能被描述为“工作树干净”或“GitHub 已同步全部最新成果”。

### 9.3 README/STATUS 文档漂移

README 和 `docs/STATUS.md` 中仍存在历史语句，例如：

- “ePhaZ broad discovery 尚未找到实体文件”已被 T141 历史 run 恢复和当前 registry 取代；
- “尚未启动 268G 全量扫描”已被 run 12 的正式运行状态取代；
- README 仍将部分科研治理事项列为未完成，需在 run 12 完成后按证据更新；
- `docs/STATUS.md` 的更新时间为 2026-08-27，不能覆盖 2026-09-01 的运行快照。

建议将本报告作为审计快照，待正式扫描完成并审核后，再统一更新 `docs/STATUS.md`、README 和最终结果报告，避免多个文档各自维护不同数字。

## 10. 当前阻断项与风险分级

### 高风险 / 发布阻断

1. run 12 未完成，无法验收 1,000 个模型×分片任务，也没有最终 `hits_all.tsv` 和 scan manifest。
2. 本地、GitHub 和 T141 deploy 不是同一 Git 快照；当前科研改动尚未全部推送。
3. OH 树和 ePhaZ CD-HIT 树仍有 stale/unregistered 输入，不得作为当前方案 A 的正式树证据。
4. ePhaZ broad 命中包含明确的 MCL-PHA 非 PHB 例子，不能用于 PHB 特异性结论。

### 中风险 / 需在发布前处理

1. README、`docs/STATUS.md`、历史审计报告之间存在状态漂移。
2. 正式扫描脚本按模型×分片串行执行，运行时间长；当前有失败记录但缺少自动恢复/断点级任务清单。
3. 超长序列虽已从 HMMER 输入中隔离，但最终报告必须明确其数量、accession 和长度，不得放入阴性分母。
4. `pending_manual=81` 的结构复核尚未完全关闭；不能将其重新纳入 ePhaZ 阳性面板。
5. 未跟踪 `--amino` 和大量运行/结果文件需在提交前逐项归类，不能直接批量删除。

### 低风险 / 已控制

1. 运行目录覆盖保护、dated deploy、输入契约、HMM SHA-256 和测试门禁已建立。
2. 旧失败 run 和历史输出未被覆盖，取证链仍完整。

## 11. 建议的后续顺序

1. 继续监控 run 12，直到 1,000 个任务全部完成；确认 `failed_tasks.tsv` 为空、超长 ledger 完整、`hits_all.tsv` 非空、最终 scan manifest 生成。
2. 对最终输出执行模型数量、分片数量、文件非空、SHA-256、聚合行数和重复 accession 检查。
3. 将正式命中按 PHB、MCL-PHA、胞内、注释型和 challenge 分层统计；不合并阳性率。
4. 以 run 12 的 manifest 为唯一当前扫描证据，更新 `docs/STATUS.md` 和 README，明确旧数字、旧 pending 语句和新状态的替换关系。
5. 先整理本地工作树并拆分“脚本、测试、结果、文档”提交，再决定是否推送 GitHub；不得在扫描未验收前发布最终 HMM 命中表或 DOI 数据包。
6. 重新注册 OH 树和 ePhaZ CD-HIT 树输入后，再考虑继续系统发育/HGT；在此之前保留 `stale_input`/`input_not_registered` 标记。

## 12. 主要证据索引

- 项目规则：[AGENTS.md](../AGENTS.md)
- 旧单一状态页：[STATUS.md](STATUS.md)（更新时间早于本报告，需更新）
- 冻结模型注册表：[formal_scan_models.tsv](../pipeline/config/formal_scan_models.tsv)
- 正式扫描脚本：[formal_frozen_screen.sh](../pipeline/scripts/formal_frozen_screen.sh)
- 超长序列处理：[filter_hmmsearch_shard.py](../pipeline/scripts/filter_hmmsearch_shard.py)
- 预检 manifest：[preflight_manifest.json](../runs/20260831_formal_frozen_preflight_10/results/preflight_manifest.json)
- 正式扫描状态：[T141_20260831_formal_frozen_scan_12_status.md](T141_20260831_formal_frozen_scan_12_status.md)
- 失败证据：[T141_20260831_formal_frozen_scan_11_failed.md](T141_20260831_formal_frozen_scan_11_failed.md)
- ePhaZ 分层：[T141_20260828_ephaz_split_01_status.md](T141_20260828_ephaz_split_01_status.md)
- 外部敏感性：[T141_20260830_ephaz_external_panel_sensitivity_03_status.md](T141_20260830_ephaz_external_panel_sensitivity_03_status.md)
- 竞争审计：[T141_20260830_ephaz_competition_audit_05_status.md](T141_20260830_ephaz_competition_audit_05_status.md)
- 570 条结构复核：[T141_20260829_ephaz_ambiguous_structure_review_02_status.md](T141_20260829_ephaz_ambiguous_structure_review_02_status.md)
- 独立面板：[T141_20260830_ephaz_independent_panels_01_status.md](T141_20260830_ephaz_independent_panels_01_status.md)
- 完成测试：[pipeline/tests](../pipeline/tests)

**审计结论状态：** `running / release-blocked / candidate-only`。
**本报告不构成 PHB 降解表型验证，也不替代正式实验验证。**
