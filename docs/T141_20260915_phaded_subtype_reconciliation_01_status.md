# T141 PhaDED 亚型证据 reconciliation amendment

日期：2026-09-15  
Run：`20260915_phaded_subtype_reconciliation_01`  
边界：本轮只做本地、只读、candidate-only 合并；没有启动服务器 HMM/结构任务，没有修改 registry、正式扫描结果或历史 run。

## 本轮完成

- 将全库 InterPro 按 accession 合并到最终台账。每条记录的 `interpro_status`（同时保留 `interpro_evidence_status`）仅使用：
  - `interpro_supported`：108,735
  - `interpro_not_reported`：1
  - `interpro_tool_input_excluded`：351
- 将 profile、domain/Pfam、motif、localization、structure、phylogeny 六层证据写入同一 accession 行，并生成 `subtype_call`、`subtype_confidence`。
- motif 字段保持独立：`motif_lipase_box`、`motif_catalytic_ser_cys`、`motif_his`、`motif_asp`、`motif_oxyanion_hole`、`motif_sbd`、`motif_linker`、`motif_lid`。当前全库只有 `partial_lipase_box_only=77,302` 和 `not_tested_full_library=31,785`；没有把 `GxSxG` 当作完整催化证据。
- 对 19 条人工重点候选完成 accession-bound 汇总：14 条 `retain_candidate_only`，3 条 `hold_architecture_conflict`，2 条 `hold_gene_model_or_structure`。
- 对冲突/模糊记录保留重叠关系：83 条 architecture conflict；42 条 superfamily ambiguity；9,664 条 family ambiguity。冲突或 ambiguity 不被自动升级为阳性或生物学阴性。
- 未训练 profile 缺口保持 `planned_not_run`：46 个 profile 中 9 个 trained、37 个 reference-only；reference-only 没有直接 calibrated score，不能从 parent-superfamily pool 推导 family call。

## 主要产物

- `runs/20260915_phaded_subtype_reconciliation_01/results/phaded_final_subtype_evidence_ledger.tsv`
- `runs/20260915_phaded_subtype_reconciliation_01/results/interpro_accession_summary.tsv`
- `runs/20260915_phaded_subtype_reconciliation_01/results/subtype_evidence_layer_counts.tsv`
- `runs/20260915_phaded_subtype_reconciliation_01/results/conflict_ambiguity_summary.tsv`
- `runs/20260915_phaded_subtype_reconciliation_01/results/priority_19_candidate_summary.tsv`
- `runs/20260915_phaded_subtype_reconciliation_01/results/profile_gap_plan.tsv`
- `runs/20260915_phaded_subtype_reconciliation_01/results/reconciliation_report.md`
- `runs/20260915_phaded_subtype_reconciliation_01/input_contract.json`
- `runs/20260915_phaded_subtype_reconciliation_01/run_manifest.json`

## 亚型主结果

| `subtype_call` | 数量 | 解释 |
|---|---:|---|
| `dPHASCL1_like_candidate` | 67,342 | 计算候选，不是表型确认 |
| `dPHASCL2_like_candidate` | 356 | 计算候选；仍需完整架构/motif复核 |
| `ambiguous_family_within_superfamily` | 9,631 | family 未解析，不能强行归类 |
| `ambiguous_superfamily` | 41 | superfamily 未解析 |
| `hold_architecture_conflict` | 83 | Pfam/架构冲突，保持 hold |
| `hold_gene_model_or_structure` | 2 | 来源完整性/结构 hold |
| `unassigned_PhaDED_like` | 31,632 | profile 覆盖缺口，保持 unresolved |

## 旧台账核对

旧 `phaded_full_evidence_ledger.tsv` 的逐行 `evidence_grade` 计数与旧 `evidence_tier_manifest.json` 一致：`HOLD_architecture_conflict=83`、`HOLD_sequence_integrity=1,951`、`L1_profile_plus_partial_architecture=65,316`、`L0_profile_only=10,467`、`PENDING_profile_unassigned=31,270`。此前中文状态文档中的另一组计数不再作为权威来源；新 amendment 绑定 TSV 和 manifest 的 SHA-256。

## 证据边界

HMM/profile、Pfam/InterPro、motif、定位、结构和系统发育只表示序列同源性或功能潜力。`not_tested`、`pending`、`ambiguous` 和 `hold` 都不是生物学阴性；任何等级均不等同于已验证 PHB/PHA 降解表型。

