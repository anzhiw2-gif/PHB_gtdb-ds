# T141 PhaDED 全库证据 amendment

日期：2026-09-16  
Run：`20260916_phaded_full_library_evidence_amendment_01`

## 本轮目的

将已有的 accession 级 InterPro、Pfam architecture、全库 motif panel、SignalP 状态、profile assignment、结构/系统发育代表性复核和 19 条人工复核记录合并为一份新的全库 candidate-only 矩阵。本轮不修改 formal scan registry，不新增 family 判定，也不把任何计算证据解释为 PHB/PHA 降解表型。

## 本轮完成

- 覆盖 `109,087` 条候选 accession。
- 使用 `20260915_phaded_motif_reconciliation_01` 的全库 motif panel 重新构建 subtype matrix。
- 保留 InterPro 三态：
  - `interpro_supported`: `108,735`
  - `interpro_tool_input_excluded`: `351`
  - `interpro_not_reported`: `1`
- 保留 profile gap：
  - `profile_trained_hit`: `67,749`
  - `profile_ambiguous_family`: `9,664`
  - `profile_ambiguous_superfamily`: `42`
  - `profile_unassigned`: `31,632`
- 保留 motif 证据边界：
  - `partial_catalytic_pattern_panel`: `76,929`
  - `conflict_relative_position`: `373`
  - `not_tested_full_library`: `31,785`
- 保留结构和系统发育的覆盖边界：
  - 结构 `not_tested_full_library`: `109,068`
  - 系统发育 `not_tested_full_library`: `109,045`
- 保留人工复核边界：19 条重点候选，其中 14 条 `retain_candidate_only`，5 条 hold。

## 解释

本轮完成的是“全库逐条纳入并显式标注证据状态”，不是把 109,087 条都升级为充分证据。完整的独立 profile、family-resolved negative、challenge set 和 held-out calibration 仍是 37 个 reference-only profile 的前置条件；在这些条件满足前，不产生新的 family 判定。

同样，`partial_catalytic_pattern_panel` 只表示 lipase-box/相关 pattern panel 的序列级支持，不等同于完整催化三联体、oxyanion-hole、SBD/linker/lid 或实验活性。`not_tested_full_library` 表示尚未运行或未绑定到该层，不表示阴性。

## 主要输出

- `runs/20260916_phaded_full_library_evidence_amendment_01/results/phaded_full_library_subtype_matrix.tsv`
- `runs/20260916_phaded_full_library_evidence_amendment_01/results/subtype_evidence_layer_counts.tsv`
- `runs/20260916_phaded_full_library_evidence_amendment_01/results/priority_19_candidate_summary.tsv`
- `runs/20260916_phaded_full_library_evidence_amendment_01/results/profile_gap_plan.tsv`
- `runs/20260916_phaded_full_library_evidence_amendment_01/input_contract.json`
- `runs/20260916_phaded_full_library_evidence_amendment_01/run_manifest.json`

## 服务器预检

2026-09-16 12:43（Asia/Shanghai）检查 `<SERVER_HOST>`：

- 80 logical CPU；
- 约 955 GiB 可用内存；
- GPU0：RTX 4090，利用率 100%，由其他结构任务占用；
- GPU1：RTX 4090，当前空闲；
- 本轮没有抢占服务器 GPU，也没有启动新的服务器计算任务。

## 下一阶段

1. 在服务器空闲资源允许且不超过单任务 40 线程的条件下，补 SignalP/定位全库运行。
2. 为 37 个 reference-only profile 获取独立阳性、family-resolved 阴性和 challenge set；完成前保持 `planned_not_run`。
3. 仅在校准后再生成新的 family-level call。
4. 结构和系统发育继续按冲突记录、模糊记录和分层代表序列扩展，不将未测试状态解释为阴性。
