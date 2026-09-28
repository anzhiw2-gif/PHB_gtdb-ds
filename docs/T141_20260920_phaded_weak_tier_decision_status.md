# T141 — 池外 weak/mid 命中层决策（审查报告第 5 项）：状态文档

**日期：** 2026-09-20
**性质：** 口径决策（**无新计算、不改任何冻结产物、不删除文件**）
**触发：** 项目全面审查报告「接下来该做什么」第 5 项
**边界：** candidate-only。

---

## 0. 决策

**池外 `weak`（E≥1e-10，3,423 条）与 `mid`（1e-30≤E<1e-10，1,503 条）**：

> **不升级**（不进入证据管线）、**不进任何候选/高可信度计数**、**保留在 score-only 层并显式标注 `score_tier=weak|mid`**。**不删除、不改写。**

与操作者推荐口径一致；下 §2 为与 AGENTS.md 的一致性核验。

## 1. 实测分布（本 run 独立逐行重算）

来源：`runs/20260919_phaded_high_confidence_relabel_01/results/pool_external_profile_hits.tsv`（35,558 行）

| score_tier | 条数 | 所属超家族分布 |
|---|---:|---|
| `strong`（E<1e-30） | **30,632** | type1 19,765 / type2 7,644 / Cys 2,049 / nPHAMCL 617 / PhaZ7-like 301 / dPHAMCL 254 / periplasmic 2 |
| `mid`（1e-30≤E<1e-10） | **1,503** | PhaZ7-like 859 / Cys 611 / dPHAMCL 30 / periplasmic 3 |
| `weak`（E≥1e-10） | **3,423** | Cys 2,454 / PhaZ7-like 929 / dPHAMCL 40 |
| **合计** | **35,558** | 全部 `evidence_layer = score_only`（35,558/35,558） |

**对账**：30,632（strong，已进证据管线）+ 4,926（mid+weak，未进）= 35,558 ✓。
**口径来源**：`score_tier` 由 `pipeline/scripts/build_phaded_pool_external_profile_layer.py` 按有效 E 值分层（strong <1e-30 / mid / weak ≥1e-10）；该 run 的 `results/pool_external_summary.json` 记 `per_score_tier = {strong: 30632, mid: 1503, weak: 3423}`。

**注意**：weak/mid 层**未曾进入**证据管线——池外证据筛选只跑 strong 层 30,632 条（见 `runs/20260919_phaded_pool_external_evidence_01/input_contract.json` 的 `result.strong_tier_total = 30632`），故 weak/mid **不在** `results/filter/high_confidence_candidates.tsv` 中。

## 2. 与 AGENTS.md 的一致性核验

**规则 1（发现层只召回不筛选）——原文逐字引用：**
> 「允许用 DED 家族序列（含 `annotation_only`）训练**发现层 HMM**，层名固定为 `discovery_hmm_uncalibrated` 且每份输出必须带该标注；发现层只用于**召回**，**不得**用于筛选、删除或降级任何候选（实测其对 563 条 x₁ 混淆嫌疑命中 443 条 = 78.69%，无判别力）。」

- ✅ **不删除**：weak/mid 全部原地保留在 score-only 层，未删任何一行。
- ✅ **不构成「降级」**：它们本就处于 weak/mid 层，**不升级 ≠ 降级**；本决策只是维持其现有层级，未把它们从任何更高层级打下来。
- ✅ **不用于「筛选」排除它们**：本决策不把它们当作被筛除的阴性，而是作为「未评估的证据层」保留。

**规则 2（池外 deferred 层不计入计数）——原文逐字引用（同构适用）：**
> 「with-lipase（`DED_hfam_2`）池外命中是已知噪声（2026-09-19 实测 1,206,655 条），必须保留为 **deferred 结构验证层**…而**禁止删除**；**该层不进入任何候选/高可信度计数**，后续只能经结构证据（PDB 8YNV Foldseek 或结构预测比对）召回，通过者再走证据层与校准 gate。」

- ✅ **不进任何计数**：weak/mid 不计入候选宇宙（109,087）、不计入池内高可信度（36,611）、不计入池外高可信度（4,131）。任何汇总里必须显式分开。
- ✅ **显式标注**：每行带 `score_tier=weak|mid` 与 `evidence_layer=score_only`，可随时区分。

**规则 3（发现层不得注册/不得产生 family 判定）：** weak/mid **不产生**任何 family 判定，**不进** `pipeline/config/formal_scan_models.tsv`，**不触发**任何 GTDB 重扫。✓

## 3. 允许 / 禁止的引用方式

**允许：**
- 作为「score-only 弱命中层」单列，写明 `score_tier` 与所属超家族；
- 用于**后续可选**的弱证据分析（例如补算证据以判断是否值得升级），但须**单独授权**且另行开 run。

**禁止：**
- 把 weak/mid 与 strong 层或任何高可信度集合**合并计数**；
- 在图/表中与高可信度候选**并列呈现而不加层级标注**（此前 §11.9.3 缺口 2 已就此提出）；
- 把 weak/mid 说成「阴性」「非候选」或删除它们。

## 4. 边界

- 本次为口径决策，**无文件改动**（分层标注在 `pool_external_profile_hits.tsv` 中已存在）。
- **未删除任何行**；历史产物原样保留。
- 若未来要升级 weak/mid：需**单独授权**的补算（SignalP + InterPro + motif 等），且升级后仍须走同一 fail-closed 过滤与校准 gate。
