# T141 — 33 个 reference-only 家族发现层 HMM 全库分类（快路径）：状态文档

> **⚠️ 数字更正声明（2026-09-28 追加；正文一字未改）**
>
> 本文正文把全库输入空间写作 **约 2.92 亿蛋白**，并据此计算占比。该值是「单分片 2.92M × 100」的**外推**，已被实测取代：冻结的 100 个分片逐片计数之和为 **615,969,589**（`runs/20260928_phaded_scan13_z_scale_reconciliation_01/results/shard_record_counts.tsv`，100 片无 pending），**旧的 ~2.92×10⁸ 偏低约 2.1 倍**。
>
> **受影响的具体数字**：本文 L107 的「唯一命中蛋白 1,525,893（占 2.92 亿蛋白的 0.52%）」按真实 Z 应为 **≈0.248%**。唯一命中蛋白**绝对数不变**，只有占比需要按实测 Z 重算。
>
> 现行权威：`docs/T141_20260928_phaded_evidence_model_redesign_status.md`（含实测 Z、逐片表与尺度修正的完整后果）。

**Run：** `runs/20260918_phaded_family_discovery_01`（本地 + 服务器同名）+ `deploy/20260918_phaded_family_discovery_01/`
**日期：** 2026-09-18
**状态：** `completed_candidate_only`（快路径：候选池分类；全库 267 GB 召回尚未做）
**证据层：** `discovery_hmm_uncalibrated`（发现层，只召回不筛选）
**服务器：** 已使用；线程策略 `T = min(40, nproc − load1 − 10)`，实测 **T=12**（启动时 load≈58）；未用 GPU

---

## 0. 一句话结论

**33 个发现层 HMM 把 109,087 条候选中的 107,924 条（98.9%）认领进家族结构；其中 31,160 条是此前 `unassigned_PhaDED_like`（占 31,632 的 98.5%）。** 家族级"最佳家族"因 HMM 重叠而**不可靠**；**超家族级分类才有意义**（78.6% 候选只命中 1 个超家族）。

---

## 1. 做了什么

1. **训练集**：为 33 个 reference-only 家族各提取 DED 成员序列（合计 723 条来源里对应部分），每个家族 1–83 条（均值 ~22）。
2. **provenance 标记**：`inputs/family_training_provenance.tsv` 逐家族记录 `n_experimental_positive`、阳性 accession、PMID、DOI。
3. **构建**：33 个 MAFFT 比对 + hmmbuild，NAME = `{family}_discovery_uncalibrated`，层名 `discovery_hmm_uncalibrated`。
4. **打分**：对 109,087 候选池（44 MB）+ 723 参考面板（自检）各打一遍。
5. **合并**：`merge_phaded_family_discovery_scores.py` 产出三表（长表 / 最佳家族 / 家族汇总）。

## 2. 核心结果

| 指标 | 数值 |
|---|---:|
| 候选宇宙 | 109,087 |
| 被家族 HMM 认领 | **107,924（98.9%）** |
| 未认领 | 1,163（1.1%） |
| 原始命中（多对多，去重前） | 1,002,293 |
| 其中之前 `unassigned_PhaDED_like` 被认领 | **31,160 / 31,632（98.5%）** |
| 其中之前 `ambiguous_family` 被认领 | 9,411 / 9,664 |
| 其中之前 `assigned` | 67,315 / 67,749 |

### 2.1 超家族级分布（有意义的一级分类，含歧义档与实验阳性标记）

分类按**三档置信度**输出（`superfamily_classification.tsv`）：

| 置信度 | 候选数 | 占比 |
|---|---:|---:|
| `unique`（唯一超家族） | 84,822 | 78.6% |
| `ambiguous_2`（2 个超家族） | 20,383 | 18.9% |
| `ambiguous_3`（3 个超家族） | 2,719 | 2.5% |

按"最佳超家族"聚合（`superfamily_classification_summary.tsv`，含实验阳性标记）：

| 超家族 | 认领候选 | 实验阳性（仅 reference-only 家族） |
|---|---:|---:|
| extracellular dPHASCL type 1 | 67,456 | 3 |
| intracellular nPHASCL without lipase box | 31,971 | 1 |
| intracellular nPHASCL with lipase box | 5,866 | 1 |
| extracellular dPHASCL type 2 | 2,625 | 0 |
| periplasmic PHA depolymerases | 6 | 1 |

> **分类口径（本次修改后）**：`superfamily_claim` 为主分类（unique 时为超家族名，歧义时显式写 `ambiguous:sf1|sf2|sf3`）；家族命中作为**透明清单**（`families_hit` / `families_in_best_superfamily`），**不再硬写"最佳家族"**。示例：`GCA_000017645.1|CP000781.1_3394` = `ambiguous:type 1|type 2`（17 个家族横跨 2 个超家族，是真实的边界案例）。


### 2.3 实验阳性来源标记（你要求的关键）

33 个家族中**只有 6 个**有实验阳性来源，其余 27 个全为 annotation-only：

| 家族 | 实验阳性数 | 阳性 accession | 被认领候选数 |
|---|---:|---|---:|
| `DED_hfam_53` | 2 | ACG63775.1; O82950 | 331 |
| `DED_hfam_65` | 1 | CAJ92291.1（= PhaZ1） | 7,114 |
| `DED_hfam_2` | 1 | EAO52570.1 | 5,866 |
| `DED_hfam_58` | 1 | BAF35850.1 | 2,210 |
| `DED_hfam_3` | 1 | AAL30107.1（periplasmic） | 6 |
| `DED_hfam_7` | 1 | AAK07742.1（PhaZ7 型） | **0** |

> `hfam_7`（PhaZ7 型）候选池 **0 命中**、参考自检 4 命中——与之前 AHSMG 结论一致：候选池不含 PhaZ7 型。

## 3. 必须诚实声明的重要边界

1. **发现层 HMM 是"召回"工具，不是"分类器"**：它们由 annotation 序列训练、无判别力（参照 Cys 发现层对 563 条 x1 混淆集命中 78.69%）。
2. **家族级"最佳家族"不可靠**：33 个家族 HMM 高度重叠——只有 8,144 条候选恰命中 1 个家族，其余命中 2–22 个家族；`ambiguous_within_1e3` 占 11.8% 只是下限。**同一超家族内的家族（如 type 1 的 14 个家族）序列高度相似，"最佳家族"是最小 E 值的任意裁决，不应作为家族级结论引用。**
3. **超家族级才是有意义的分类**：78.6% 候选只命中 1 个超家族；跨超家族歧义（2–3 个超家族）占 21.4%。
4. **本轮只做了快路径**（候选池），**未做全库 267 GB 召回**——即"池外漏掉的候选"仍未知（参照 Cys 那次池外占 12%）。
5. 所有结果仍是 candidate-only，不等同于已验证 PHB/PHA 降解表型；不产生 family 判定、不进 registry、不改 `formal_scan_models.tsv`。

## 4. 产物

| 产物 | 位置 |
|---|---|
| 训练集（33 个） | `inputs/training/DED_hfam_*.faa` |
| provenance 清单 | `inputs/family_training_provenance.tsv` |
| HMM（33 个） | `results/hmms/DED_hfam_*.hmm` |
| 候选池 tblout（33 个） | `results/pool_scores/*.tbl` |
| 参考面板 tblout（33 个） | `results/ref_scores/*.tbl` |
| 长表（多对多认领） | `results/family_pool_claims_long.tsv` |
| **超家族分类表（主分类）** | `results/superfamily_classification.tsv` |
| **超家族汇总（含实验阳性）** | `results/superfamily_classification_summary.tsv` |
| 家族汇总（含实验阳性） | `results/family_discovery_summary.tsv` |
| 构建/打分日志 | `logs/build_score.log`、`logs/{mafft,hmmbuild,hmmsearch}.err` |
| 脚本 | `pipeline/scripts/build_phaded_family_discovery_sets.py`、`pipeline/scripts/merge_phaded_family_discovery_scores.py` |

## 5. 下一步

全库 267 GB（100 分片）召回**已完成**（2026-09-19 04:31 跑完）。结果与诚实解读见 §6。

---

## 6. 全库召回结果（2026-09-19，`results/shard_pool_external.tsv`，98 MB）

**规模**：3300 个 (HMM×分片) tblout，原始命中 2,577,233 行，去重后**唯一命中蛋白 1,525,893**（占 2.92 亿蛋白的 0.52%），其中**池外候选 1,419,842**。

**但必须诚实解读——这 1.42M 大部分不是 PHB 降解候选**。池外候选按"最佳家族所属超家族"分布极不均匀：

| 超家族 | 池外候选 | 占比 | 解读 |
|---|---:|---:|---|
| intracellular nPHASCL **with lipase box**（=仅 hfam_2） | **1,238,812** | 87.2% | ⚠️ hfam_2 的 HMM（29 条序列）极宽，实质是"通用 lipase-box α/β-水解酶检测器"，1.24M 命中是**背景噪声**，不是 PHB 候选 |
| extracellular dPHASCL type 1 | 170,397 | 12.0% | 仍偏宽（14 个家族 HMM 重叠） |
| intracellular nPHASCL without lipase box（Cys） | 5,114 | 0.36% | 相对特异 |
| extracellular dPHASCL type 2 | 3,425 | 0.24% | 相对特异 |
| extracellular native-SCL/PhaZ7-like | **2,089** | 0.15% | ⭐ **候选池完全漏掉了这一族**（池内 0 命中、池外 2,089） |
| periplasmic | 5 | — | 极特异 |

跨超家族歧义率仅 **0.4%**（5,378 / 1,419,842）——超家族级信号是干净的。

**三个结论**：

1. **发现层家族 HMM 的特异性差异极大**：hfam_2（有 lipase box）极宽（1.24M = 通用 lipase-box 背景），而 Cys 型（Cys-His-Asp 特异）、type 2、PhaZ7、periplasmic 相对特异（合计约 1.1 万条池外候选）。这再次印证"发现层只能召回、不能筛选"——**家族级 HMM 的召回数不能当候选数用**。
2. **候选池确实漏掉了 PhaZ7 型**：池内 0 命中、全库池外 2,089 命中。这是"旧 registry + 分层漏斗"漏掉一个整族的实证。
3. **真正可操作的池外候选**是特异家族：Cys 5,114 + type2 3,425 + PhaZ7 2,089 + periplasmic 5 ≈ **1.1 万条**；type 1 的 17 万条需进一步证据过滤；hfam_2 的 1.24M 是噪声，不应作为候选。
