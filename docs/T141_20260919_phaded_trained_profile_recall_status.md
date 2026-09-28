# T141 — trained profile 全库分片扫描（nPHAMCL/dPHAMCL 池外召回）：状态文档

**Run：** `runs/20260919_phaded_trained_profile_recall_01`（服务器 + deploy 绑定）
**日期：** 2026-09-19
**状态：** `completed_candidate_only`
**模型：** 9 个 trained profile（5 家族 + 4 超家族，来自 `profile_build_01/profiles/`，E<1e-5）
**执行：** 分片外层循环（每分片读一次），9 HMM 并行、每任务 `--cpu 4`（总 36 线程，load≈33/80，未超载）；约 **40 分钟** 跑完 900 个 (HMM×分片) 任务，0 错误。

---

## 1. 总体结果

| 指标 | 数值 |
|---|---:|
| 唯一命中蛋白 | 501,008 |
| 池内 | 70,861 |
| **池外候选** | **430,147** |

## 2. 池外候选按超家族分布（trained profile 覆盖的 4 个超家族）

| 超家族 | 池外命中（E<1e-5） | 强命中（<1e-30） | 解读 |
|---|---:|---:|---|
| intracellular **nPHAMCL** | **285,423** | **5,938** | ⚠️ 3 条训练序列的 HMM 过召：73% 是临界 E 值噪声（1e-10~1e-5），强命中才是可靠信号 |
| extracellular dPHASCL **type 2** | 85,671 | 11,677 | 偏宽 |
| extracellular dPHASCL **type 1** | 58,721 | 23,848 | 偏宽 |
| extracellular **dPHAMCL** | **332** | **460** | ✅ 特异（714 总命中中 460 强命中） |

## 3. 关键发现

1. **dPHAMCL 是唯一"池外召回干净"的 trained 超家族**：全库只有 714 条命中（池内 292 + 池外 332 + 重合），其中 460 条 E<1e-30 强命中。这个家族在全蛋白空间确实很小、很特异。
2. **nPHAMCL 的 285k 池外命中是过召**：hfam_4 的 HMM 只有 3 条训练序列，在 E<1e-5 阈值下匹配了大量通用 α/β-水解酶背景（73% 为临界值）。**这再次印证"任何 profile HMM 扫全蛋白空间都会过召；必须用强阈值 + 证据层过滤"**。
3. **与发现层扫描的结论完全一致**：无论 trained 还是 discovery HMM，"E<1e-5 全库命中数"都不是候选数；只有强命中（<1e-30）+ 后续证据层（motif/SignalP/InterPro）才有意义。

## 4. 八个超家族池外召回全景（合并发现层 + trained 两次扫描）

| 超家族 | 模型 | 池外命中 | 可信度 |
|---|---|---|---:|
| type 1 | 发现层 14 家族 + trained | ~22.9 万 | 宽（强命中 ~2.4 万） |
| type 2 | 发现层 7 家族 + trained | ~8.9 万 | 宽（强命中 ~1.2 万） |
| Cys（无 lipase box） | 发现层 9 家族 | 5,114 | 相对特异 |
| with lipase box | 发现层 hfam_2 | 1,238,812 | 噪声（通用 lipase-box） |
| PhaZ7-like | 发现层 hfam_7 | 2,089 | 特异 |
| periplasmic | 发现层 hfam_3 | 5 | 极特异 |
| **nPHAMCL** | trained | 285,423 | 过召（强命中 5,938） |
| **dPHAMCL** | trained | **332** | **特异** |

## 5. 产物

- 服务器 `results/shard_scores/`：900 个 tblout
- 服务器 `results/trained_pool_external.tsv`：池外候选清单（43 万行）
- `logs/scan.log`、`logs/hmmsearch.err`（0 错误）

## 6. 边界

所有命中为 candidate-only 同源/功能潜力证据；池外"命中数"≠候选数；不等同于已验证 PHB/PHA 降解表型。
