# T141 — Cys 发现层 HMM 有界定向召回：状态文档

> **⚠️ 数字更正声明（2026-09-28 追加；正文一字未改）**
>
> 本文正文把全库输入空间写作 **约 2.92 亿蛋白**，并据此计算命中率。该值是「单分片 2.92M × 100」的**外推**，已被实测取代：冻结的 100 个分片逐片计数之和为 **615,969,589**（`runs/20260928_phaded_scan13_z_scale_reconciliation_01/results/shard_record_counts.tsv`，100 片无 pending），**旧的 ~2.92×10⁸ 偏低约 2.1 倍**。
>
> **受影响的具体数字**：本文 L77 的「全蛋白空间命中率 0.012%（35,504/2.92 亿）」按真实 Z 应为 **≈0.0058%**（35,504/615,969,589）。该段的**结论不变**（发现层高召回、低特异性，只能用于召回），但引用其比率前必须先按实测 Z 重算。
>
> 现行权威：`docs/T141_20260928_phaded_evidence_model_redesign_status.md`（含实测 Z、逐片表与尺度修正的完整后果）。

**Run：** `runs/20260917_phaded_cys_targeted_recall_01`（本地与服务器同名）+ `deploy/20260917_phaded_cys_targeted_recall_01/`
**日期：** 2026-09-17
**状态：** `completed_candidate_only`
**证据层：** `discovery_hmm_uncalibrated`（发现层，可召回不可筛选）
**服务器：** 已使用；每任务 `--cpu 10`、并发 4、**总线程 40**（上限 40，启动时 load1=8.19，余量约 32 核）；未使用 GPU

---

## 0. 一句话结论

**不需要为召回而重扫全库。** 旧 10 模型 registry 对 Cys 型家族的召回率实测为 **99.94%**（35,484 / 35,504）；对完整蛋白空间用新模型重扫的净增益只有 **20 条**候选。
**真正被丢掉的是分层漏斗：4,301 条 Cys 型候选已被 scan 13 捞到，却在 Scheme-A 构造 109,087 候选池时被排除**——这批可以从既有 `hits_all.tsv` 直接恢复，**不需要任何重扫**。

---

## 1. 设计（预注册于打分之前）

| 项 | 值 |
|---|---|
| 模型 | `runs/20260917_phaded_cys_discovery_hmm_01/results/cys_discovery.hmm`（SHA-256 前缀 `70e1b347d8f8d12d`） |
| 输入空间 | `runs/20260901_formal_frozen_scan_13/inputs/scan_shards/shard_0001..0100.faa`，**267 GB**，约 **2.92 亿**条蛋白（单分片实测 292 万条） |
| 阈值 | `E < 1e-5`（预注册，**未调整**） |
| 线程政策 | 每任务 `--cpu 10` × 并发 4 = 总 40（上限 40，保留 ≥10 核余量；启动 load1=8.19） |
| 运行时长 | 全量 100 分片约 **3 分钟**（分片命中页缓存；单分片 8–9 秒） |
| 预注册文件 | `results/preregistration.json`（打分前写入） |

**为什么这是有界操作而非全库重扫**：分片是 scan 13 已存在的中间产物，本轮只对其**打分**，没有重跑 Pyrodigal/QC，没有重建 registry，没有修改 `formal_scan_models.tsv`。

---

## 2. 结果（全部实测，集合运算已在服务器用 `comm` 独立复核）

| 指标 | 数值 |
|---|---:|
| 全空间 Cys 命中（唯一蛋白，E<1e-5） | **35,504** |
| 其中已被旧 10 模型命中 | **35,484** |
| **旧 registry 对 Cys 型的召回率** | **99.9437%** |
| 真新召回（旧模型完全未命中） | **20** |
| 扫描命中中位于 109,087 候选池内 | 31,183 |
| 扫描命中中不在候选池内 | **4,321** |
| └ 其中**已在旧原始命中表**（=分层漏斗丢弃） | **4,301** |
| └ 其中不在旧原始命中表（真新） | **20** |
| 候选池是否 ⊆ 旧原始命中表 | **109,034 / 109,087（99.95%）** |

**4,301 条被丢弃者是由哪些旧模型捞到的**：iPhaZ **4,244**、PhaC **1,970**（部分同时命中）。
**全部 35,504 条 Cys 命中的旧模型来源**：iPhaZ **35,427（99.8%）**、PhaC 22,995、ePhaZ_broad_discovery 7 —— 印证 Cys 型家族与旧 `iPhaZ` 模型高度重叠。

**分片分布**：100/100 分片均有命中；每片 30–778 条（相差约 26 倍），说明分片按基因组 accession 切分、命中在分类学上不完全均匀。

---

## 3. 可执行结论（按价值排序）

1. **恢复 4,301 条被分层漏斗丢弃的 Cys 型候选** —— 从既有 `hits_all.tsv`（6,595,162 唯一蛋白）取出，无需重扫、无需新授权。这是本轮最有价值的发现。
2. **20 条真新召回**（`results/hits_NOT_in_old_registry.txt`）单列为优先复查队列。
3. **不需要 GTDB 全库重扫**；若将来仍要做，须走单独的 registry decision + 授权。
4. **修正上游文档错误**：`docs/T141_20260909_project_handoff.md` 的 Figure 1 漏斗把流程写成 `Pyrodigal+FAA/GFF QC → DIAMOND blastp (6,769,772) → HMMER 3.4`。实测 scan 13 全 run 内 `grep -ril diamond` **零命中**，`scan_manifest.json` 显示为纯 hmmsearch（100 分片 × 10 模型 = 1,000 任务）；`6,769,772` 实为 `results/logs/rerun_candidates.log` 中早先一次 900-tbl 聚合的命中行数。本 run 不改历史文档，仅记录。

---

## 4. 本轮自查发现的过程缺陷（如实留证）

| ID | 问题 | 证据 / 修正 |
|---|---|---|
| `recall_launch_bug_01` | 首次启动用 `xargs -n 1` 且脚本内含字面 `{}` 占位符 → 参数字位错乱（脚本收到文件名 `{}`），未产出任何 tblout | `logs/failed_attempt_01.txt`；改为把 shard 作为 xargs 追加的末位参数（签名 `OUTDIR HMM CPU SHARD`），单分片自检通过后才跑全量 |
| `pilot_extrapolation_error` | 试点选了**最小**分片（30 条命中，恰为全库最小值），外推全库约 3,000 条，实测 35,504，**低估 11.8 倍** | 教训：分片命中数相差 26 倍，试点须用中位分片，或直接跑全量（仅需数分钟） |
| `comm_crlf_error` | 池 ID 经 PowerShell `Set-Content` 上传为 CRLF，`comm` 首次比较给出**错误的 0 交集**和"池内 109,087 条全不在旧命中表"的假结论 | 加 `tr -d '\r'` 后重做，得 109,034/109,087 与 31,183 的正确值；已在 `recall_summary.json.validation.method_note` 记录 |
| `pilot_search_order` | （承接操作者指出的问题）在服务器上以 1 TB `runs/` 开头做单线程 `grep -r`，为找 0.6 GB 小目录里的匹配先扫 1 TB | 改用**元数据优先**：先 `find -name/-size`，内容 grep 只在小目录、并行 + `head` 短路、大目录加 `-maxdepth` |

---

## 5. 边界声明

- 本轮为**发现层候选召回**：不产生 family 判定、不写 registry、不改 `pipeline/config/formal_scan_models.tsv`、对现有 `subtype_call` 影响 **0 行**。
- 该发现层 HMM **高召回、低特异性**：全蛋白空间命中率 0.012%（35,504/2.92 亿），但候选池内命中 29.25%、563 条 x₁ 脂酶混淆集命中 **78.69%** —— 因此**只能用于召回，不得用于筛选、删除或降级任何候选**（与 AGENTS.md 一致）。
- 所有命中均为候选层同源/功能潜力证据，**不等同于已验证 PHB/PHA 降解表型**。

---

## 6. 产物

| 产物 | 位置 |
|---|---|
| 预注册 | 服务器与本地 `results/preregistration.json` |
| 100 个 tblout | 服务器 `results/shard_0001..0100.tbl` |
| 全空间命中清单（35,504） | `results/recall_unique_proteins.txt` |
| 真新召回（20） | `results/hits_NOT_in_old_registry.txt` |
| 已被旧模型命中（35,484） | `results/hits_already_in_old_registry.txt` |
| 分片命中分布 | `results/per_shard_hits.txt` |
| 汇总与验证 | `results/recall_summary.json` |
| 池 ID 清单（109,087） | `inputs/pool_accessions_109087.txt` |
| 失败留证 | `logs/failed_attempt_01.txt`、`logs/scan.log`、`logs/hmmsearch.err` |

---

## 7. 四路独立核查后的更正（2026-09-17 晚，`runs/20260917_phaded_task_completion_check_01/`）

核心数字**全部独立复算通过**（35,504 / 35,484 / 20 / 31,183 / 4,321 / 4,301 / 109,034 / 阈值 / 预注册时序 / -Z 机理）。以下为本状态文档与 run 产物中需要更正的表述（按核查报告 P1–P8）：

| # | 更正 | 原文（错误） | 实测（正确） |
|---|---|---|---|
| P1 | 运行时长 | 「约 3 分钟」「单片 8–9 秒」 | 墙钟 **499 s ≈ 8 分 19 秒**；单片 11–37 s，均值 20.1 s |
| P2 | 预注册时序 | （暗示试点也在预注册后） | 试点 tbl 早于预注册 **78 s**；预注册仍先于全部 100 个全量 tblout，无阈值篡改 |
| P3 | 最小分片 | 「shard_0051 恰好是最小值」 | **shard_0050 同为 30 条**，两者并列最小 |
| P4 | 输入字节数 | `287,146,000,000` | 实测 `286,194,724,100`（差 0.33%，267 GB 量级结论不变） |
| P5 | 证据保全 | （无） | xargs 驱动命令未留档；`--cpu 10` 只能由并发重建佐证（unverifiable） |
| P6 | /tmp 卫生 | （未声明） | /tmp 遗留 19 个本轮文件（含 CRLF 误配失败产物），已只报告未删除 |
| P7 | 可复算性 | （未声明） | 池打分 tblout 未留档 → **31,906 与 723 差异无法独立重算**（机理已双实验证实） |
| P8 | 产物表 | 未区分本地/服务器 | `recall_summary.json` 仅本地；服务器 `inputs/` 仅 `shard_list.txt` |

**已修复项：**

- **孤儿测试**（曾被取消的工作流留下 `test_parse_phaded_cys_targeted_recall.py` 而无实现 → 套件曾红 30 errors）：已按其 spec 补写 `pipeline/scripts/parse_phaded_cys_targeted_recall.py`，**全套件 600 tests OK (skipped=1)**。
- **run 目录缺 `input_contract.json`**：已按 run_context 结构补写（4,022 B；gtdb 三槽 pending；5 个输入实测哈希）。
- **线程上限 40 未写入 AGENTS.md**：已补一行（动态取 `min(40, 核数−负载−10)`、保留约 10 核余量）。
- **根目录 0 字节残留 `s01_recon.sh.b64`**：已移入本 run 的 `logs/server_scripts/`（未删除）。
- 文档路径精确化：`rerun_candidates.log` 位于项目根 `results/logs/`（不在 run13/run12 内）。
