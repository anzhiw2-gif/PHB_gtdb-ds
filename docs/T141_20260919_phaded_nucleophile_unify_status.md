# T141 — 池内亲核体类型统一 + 矛盾行标记 hold：状态文档

**Run：** `runs/20260919_phaded_nucleophile_unify_01`
**日期：** 2026-09-19
**状态：** `completed_candidate_only`
**性质：** **加性丰富**（additive enrichment）。过滤器的通过/失败逻辑**逐字不变**，只新增三样东西：`nucleophile_type` / `nucleophile_family` 两列 + `hold_candidates.tsv` 显式 hold 表。

---

## 0. 一句话结论

池内不同超家族的序列证据**确实不一样**；其中架构特异差异是**对的**（by design），而亲核体证据不对称（Cys 无直接标记、His/Asp 普遍不可判定）是真实缺口。本 run 把亲核体类型显式化（Ser 用 GxSxG、Cys 用家族归属），并把所有矛盾行显式标 hold，**高可信度 36,611 条不变**，但 71,860 条被排除候选现在有明确的 `hold_reason`。

## 1. 逐超家族证据覆盖审计（实测，n=109,087）

| 超家族 | n | lipase box | PF06850 | lid | 氧阴离子孔/type | His 可判定 | Asp 可判定 | SignalP 分泌 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| type 1 | 67,706 | 99.96% | 0.04% | 0.02% | 24.7% / type1 16,529 | 10.4% | 1.9% | 28.6% |
| 无 lipase box（Cys） | 30,977 | 2.2%(684) | 96.5% | 0.02% | 无 | 0.8% | 0% | 0% |
| nPHAMCL | 6,959 | 84.8%(1,059 缺) | 9.9% | 41.5% | 无 | 0.16% | 0.01% | 0% |
| type 2 | 2,658 | 100% | 0 | 0 | 82.6% / type2 2,362 | 0.45% | 0.68% | 89.6% |
| dPHAMCL | 292 | 100% | 0 | 0 | 无 | 1.4% | 0.3% | 32.9% |
| 有 lipase box | 22 | 3/19 | 3 | 0 | 无 | 4.5% | 0% | 0% |

**两类"不一样"**：① 架构特异的（PF06850= Cys 标记、lid = nPHAMCL 标记、氧阴离子孔/type = type1/2 标记）——正确；② 普遍结构性缺口——His/Asp 在几乎所有家族 `not_assignable` 达 89–99.9%（根因 `candidate_ded_alignment_members=0`，不可靠比对补）。

## 2. 本次落地产物

- `nucleophile_type`：`ser`（GxSxG）/ `cys`（Cys 家族 + 无 GxSxG）/ `ahsmg`（PhaZ7）/ `conflict`（矛盾）/ `undetermined`。
- `nucleophile_family`：直接亲核体标记——Cys 行为其家族 HMM（`families_in_best_superfamily`），Ser 行为 `GxSxG`，PhaZ7 行为 `AHSMG`。
- `hold_candidates.tsv`：71,860 条被排除候选，每条带 `hold_reason`。

## 3. 结果（实测，非转录）

- 高可信度 **36,611 条不变**；其中 `nucleophile_type` = `cys` 28,534 / `ser` 8,077（与超家族完全一致，无 pass 行带 conflict）。
- hold 分解：`localization_conflict` **48,488**（type1 等胞外标签无分泌信号）/ `architecture_criteria` 19,562 / `common_criteria` 2,538 / `nucleophile_conflict` **1,272**（其中 cys 带 GxSxG 233 + Ser 家族缺 GxSxG 1,039）。
- 对账：36,611 + 71,860 + 616（未分配）= 109,087 ✓。

## 4. 治理

- test-first：新增 4 个测试（`test_filter_phaded_nucleophile.py`），红灯→绿灯；全套 **624 tests OK (skipped=1)**；`compileall` exit 0；`git diff --check` exit 0。
- 新 run 结构完整：`input_contract.json`（gtdb 三槽 pending）、`logs/`（run + verify + acceptance + 契约生成）、`results/`（high_confidence + hold + summary）。
- 未 commit/push；未覆盖/删除任何历史产物；判据逻辑未改，仅加列加表。

## 5. 边界（不变）

- 全部候选 candidate-only；亲核体类型是**家族级标记证据**（Ser 经 GxSxG、Cys 经家族归属），不是逐条催化残基生化验证，也不是已验证 PHB/PHA 降解表型。
- His/Asp 的"逐条坐标验证"因结构断点**不可补**；本 run 未伪造该证据。
