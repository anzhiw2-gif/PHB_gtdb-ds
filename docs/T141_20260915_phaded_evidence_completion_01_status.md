# T141 PhaDED 证据完成 amendment

日期：2026-09-15  
主运行：`20260915_phaded_evidence_completion_01`  
范围：序列、Pfam/InterPro、motif、定位、结构和系统发育的 candidate-only 汇总；不构成 PHB/PHA 降解表型确认。

## 本轮完成

1. 建立 `20260915_phaded_motif_reconciliation_01`，对 `109,087` 条 accession 逐条检查 lipase box、catalytic Ser/Cys、His、Asp、oxyanion-hole，以及 SBD/linker/lid 状态，并核对 723 条参考序列和 38 个 family alignment 文件的哈希；其中 29 个含实际序列行，9 个为 header-only，均保留为审计状态。
2. 建立 `20260915_phaded_conflict_motif_review_01`，把 83 条 architecture conflict 和 42 条 superfamily ambiguity 与 motif 状态绑定；83 条全部保持 `hold_architecture_conflict`，其余 41 条保持 `ambiguous_superfamily_unresolved`。
3. 建立 `20260915_phaded_ambiguity_stratified_review_01`，把 9,664 条 family ambiguity 分成 13 个证据分层，按 accession SHA-256 确定性抽取 70 条复核样本；未降低阈值、未将 ambiguity 视为阴性。
4. 合并为最终候选等级和证据覆盖统计，保留 37 个 reference-only profile 的 `planned_not_run` 边界。

## 关键数量

| 层 | 状态 | 数量 |
|---|---|---:|
| motif | `partial_catalytic_pattern_panel` | 76,929 |
| motif | `conflict_relative_position` | 373 |
| motif | `not_tested_full_library` | 31,785 |
| architecture | `domain_conflict` | 83 |
| profile | `profile_ambiguous_superfamily` | 42 |
| profile | `profile_ambiguous_family` | 9,664 |
| profile | `profile_unassigned` | 31,632 |
| profile | `reference-only` | 37 |

## 主要输出

- `runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv`
- `runs/20260915_phaded_motif_reconciliation_01/results/reference_motif_panel.tsv`
- `runs/20260915_phaded_motif_reconciliation_01/results/reference_alignment_manifest.tsv`
- `runs/20260915_phaded_motif_reconciliation_01/results/phaded_motif_reconciled_ledger.tsv`
- `runs/20260915_phaded_conflict_motif_review_01/results/conflict_superfamily_focus.tsv`
- `runs/20260915_phaded_ambiguity_stratified_review_01/results/family_ambiguity_stratified_sample.tsv`
- `runs/20260915_phaded_evidence_completion_01/results/phaded_evidence_completion_report.md`
- `runs/20260915_phaded_evidence_completion_01/results/phaded_evidence_coverage.tsv`

## 证据边界

- `GxSxG`、His/Asp/oxyanion-hole 模式和 Pfam 域状态是序列/结构域模式证据，不是完整催化机制证明。
- architecture conflict 的优先级高于 profile、InterPro、motif、结构和系统发育支持，不能自动解除 hold。
- 37 个 reference-only profile 当前缺少足够独立阳性、family-resolved 阴性、challenge set 和留出校准，因此保持 `planned_not_run`，不得降低阈值强行归类。
- 所有等级和 subtype call 均为 candidate-only，不等同于已验证 PHB/PHA 降解菌。

## 运行边界

- 服务器预检已记录：`<SERVER_HOST>` 为 80 logical CPU、约 967 GiB 可用内存、2 × RTX 4090（0% 使用率）。
- 本轮 motif、冲突和分层统计均为确定性本地只读计算，未启动服务器 HMM/结构/GPU 任务；服务器单任务线程上限仍按 `AGENTS.md` 的 40 执行。
- 新源码均已放入对应 dated `deploy/<run_id>/inputs/`，历史 `runs/`、registry 和正式扫描结果未覆盖。
