# T141 PhaDED accession、结构域冲突与分层统计

日期：2026-09-13  
运行：`20260913_phaded_accession_structure_stratification_01`  
状态：已完成（只读重分析；未执行湿实验、GTDB 全库重扫或正式 registry 修改）

## 结论

1. 主分类输入为 `109,087` 条候选蛋白、`62,666` 个去重 genome。分层统计见 `runs/20260913_phaded_accession_structure_stratification_01/results/phaded_stratified_statistics.tsv`。
2. PhaDED 最佳 superfamily：extracellular dPHASCL type 1 为 `67,544`（61.92%），intracellular nPHAMCL 为 `6,959`（6.38%），extracellular dPHASCL type 2 为 `2,658`（2.44%），extracellular dPHAMCL 为 `292`（0.27%），无最佳 superfamily 为 `31,634`（29.00%）。
3. assignment status：`assigned` `67,749`，`ambiguous_family` `9,664`，`ambiguous_superfamily` `42`，`unassigned_PhaDED_like` `31,632`。`unassigned_PhaDED_like` 是未分配状态，不是生物学阴性。
4. 19 条人工重点候选中，3 条为 `hold_architecture_conflict`，2 条为 `hold_gene_model_or_structure`，共 5 条排除出保留候选；14 条保留为 `candidate_only`。逐条理由见 `architecture_conflict_exclusion.tsv`。
5. accession 文献桥接中，唯一 exact Tier A 记录为 `NC_008313.1` ↔ `RS_GCF_000009285.1`（*Cupriavidus necator* H16，PMID `12775684`）。该 genome 在本次候选集无 exact accession 重叠，因此候选集的 exact 文献支持数为 `0`。`O05527` 仅产生 taxonomy-only 上下文，不可升级为 Tier A。

## 证据边界

- HMM/PhaDED、Pfam/InterPro、SignalP、结构和系统发育结果表示候选同源或功能潜力，不等同于已验证 PHB 降解表型。
- 本轮结构域冲突排除只对已有的 19 条人工复核记录执行。主 assignment evidence 中全库 architecture 字段仍为 `pending`，所以不能据此宣称全库已完成结构域冲突清除。
- 9 个 profile 为正式训练 HMM，其余 profile 为 `reference_only`；family 层级结果须结合 `profile_manifest.tsv` 的 model status 解读。

## 产物

- `runs/20260913_phaded_accession_structure_stratification_01/input_contract.json`
- `runs/20260913_phaded_accession_structure_stratification_01/results/run_manifest.json`
- `runs/20260913_phaded_accession_structure_stratification_01/results/phaded_stratified_statistics.tsv`
- `runs/20260913_phaded_accession_structure_stratification_01/results/architecture_conflict_exclusion.tsv`
- `runs/20260913_phaded_accession_structure_stratification_01/results/accession_literature_bridge.tsv`

所有输入输出路径、文件大小和 SHA-256 已写入 `run_manifest.json`。输入来自主分类 run `20260909_phaded_candidate_mapping_01`、19 条人工复核 run `20260913_phaded_19_candidate_manual_review_01` 和完整 GTDB R232 文献桥接 run `20260913_gtdb_phb_degrader_literature_reconciliation_02`。
