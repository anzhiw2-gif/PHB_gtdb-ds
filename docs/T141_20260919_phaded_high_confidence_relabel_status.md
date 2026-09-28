# T141 — 高可信度过滤器可靠化（§11.9.4 三条）：状态文档

**Run：** `runs/20260919_phaded_high_confidence_relabel_01`
**日期：** 2026-09-19
**状态：** `completed_candidate_only`
**性质：** **加性诚实标注**（additive labelling）。池内过滤器的通过/失败判据**逐字不变**；本 run 只新增三处诚实标注，对应交接文档 §11.9.4 的三条可靠化建议。

---

## 0. 一句话结论

对 `20260919_phaded_high_confidence_filter_01` 的高可信度产物做**只加标注、不改判据**的可靠化：池内 Cys 候选被显式标注为"推定 Cys 型（催化残基未逐条验证）"，池外候选按 E 值分层并**改名去"高可信度"**。池内 36,611 条、池外 35,558 条与冻结产物**逐条一致**，仅新增标注列。

## 1. 三条可靠化建议的落地产物

| §11.9.4 | 落地 | 证据 |
|---|---|---|
| ① Cys 催化残基未逐条验证 | 池内新增列 `catalytic_residue_verification`（Cys=`inferred_cys_not_verified`，其余=`motif_level`）+ `filter_summary.json` 的 `honesty_notes` | `results/relabel_verification.json` |
| ② 池外 E 值分层 | 池外新增列 `score_tier`（strong/mid/weak） | 同上 |
| ③ 池外整层分隔、去"高可信度" | 原 `pool_external_high_confidence.tsv` → `pool_external_profile_hits.tsv`；每行 `evidence_layer=score_only`、架构/定位/InterPro 三列 `pending` | 同上 |

## 2. 关键验证（本 run 独立复算，非转录）

**池内：**
- 重跑 `filter_phaded_high_confidence.py` → **36,611 条**，与冻结产物一致；14 个既有列与冻结 `high_confidence_candidates.tsv` **逐字节一致**（`shared_columns_byte_identical=true`），仅新增 1 列。
- Cys 高可信度 **28,534** 条全部 `inferred_cys_not_verified`；其余 **8,077** 条 `motif_level`。

**池外（首次代码化，provenance 缺口已闭合）：**
- 新增 `build_phaded_pool_external_profile_layer.py`，复用冻结的 `classify_pool_external` 归属规则，复现已发布 **35,558** 条锚点（逐超家族一致：Cys 5,114 / type1 19,765 / type2 7,644 / PhaZ7 2,089 / nPHAMCL 617 / dPHAMCL 324 / periplasmic 5 / with-lipase 0）。
- 分层：strong **30,632** / mid **1,503** / weak **3,423**。其中 **Cys 弱命中 2,454（48%）、strong 2,049（40%）、mid 611** —— 与交接文档 §11.9.3 缺口 2 的实测数字一致。

## 3. 治理

- 新 run 结构完整：`logs/`（含 `run_filter.log`、`run_pool_external.log`、`verify_relabel.log`、`acceptance_*.log` + `acceptance_evidence_note.json`）、`inputs/`（`relabel_supplement.md`）、`results/`、`input_contract.json`（gtdb 三槽 `pending`）。
- 验收：`python -m unittest discover -s pipeline/tests` → **Ran 620 tests … OK (skipped=1)**（基线 600 + 本 run 新增 11 测试 + 并发任务新增 9，全部绿）；`compileall` exit 0；`git diff --check` exit 0。
- 未 commit / push；未覆盖、未删除任何历史产物；冻结的 `runs/20260919_phaded_high_confidence_filter_01`（含原 `pool_external_high_confidence.tsv`）原样保留。
- 代码改动 test-first：先写 11 个失败测试（红灯：1 failure + 10 errors），再实现至绿灯。

## 4. 边界（不变）

- 全部产物仍为 candidate-only；"高可信度"= 候选证据最强筛选层，池外层 = 仅 profile 得分层（score_only），均**不等同于已验证 PHB/PHA 降解表型**。
- 判据缺失一律不通过（fail-closed）；不把 pending 写成阴性。
- 池外 weak（E≥1e-10）命中**未删除、未降级为"非候选"**，只是显式标注为"score-only 弱"，供下游自行决定是否用于弱证据分析。
