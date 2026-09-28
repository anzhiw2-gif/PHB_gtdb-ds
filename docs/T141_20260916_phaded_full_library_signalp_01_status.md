# T141 PhaDED 全库 SignalP 定位证据

日期：2026-09-16  
Run：`20260916_phaded_full_library_signalp_01`

## 当前状态

SignalP 6.0 fast 模式已在 `<SERVER_HOST>` 完成，使用 dated deploy：

- 输入：`108,736` 条可运行序列；
- 工具输入排除：`351` 条，沿用 InterPro exclusion manifest；
- 参数：`torch_num_threads=32`、`write_procs=8`，合计不超过项目规定的 40；
- GPU：不使用；
- 预测结果：`108,736` 条，返回码 `0`；
- 预测类别：`SP=33,920`、`LIPO=11,434`、`TAT=329`、`TATLIPO=68`、`OTHER=62,985`；
- 加上 `351` 条工具输入排除后，完整覆盖 `109,087` 条候选。
- 本次实际墙钟时间约为 `2 小时 41 分钟`（约 13:05 至 15:46，2026-09-16，Asia/Shanghai）。

## 目的

补全全库定位层，逐条生成：

- `signalp_supported`：SignalP 已产生预测；
- `signalp_not_reported`：在输入中但结果未报告；
- `signalp_tool_input_excluded`：因工具输入限制未运行。

这些状态只用于定位证据，不代表分泌、底物特异性或 PHB/PHA 降解表型。

## 后处理已完成

1. accession 集合和重复检查通过；
2. 与全库 motif panel 合并；
3. 重新生成 subtype evidence matrix；
4. 更新 `localization_evidence_status`：
   - `export_signal_supported=45,751`
   - `localization_unassigned=62,985`
   - `localization_tool_input_excluded=351`
5. 保留 37 个 reference-only profile 的 `planned_not_run`，没有新增 family 判定。

## 当前亚型结论

SignalP 补全了定位层，但没有自动改变 profile/domain/motif 冲突记录的 subtype call。当前最终矩阵仍为 candidate-only：

- `dPHASCL1_like_candidate=67,342`
- `dPHASCL2_like_candidate=356`
- `ambiguous_family_within_superfamily=9,631`
- `ambiguous_superfamily=41`
- `unassigned_PhaDED_like=31,632`
- `hold_architecture_conflict=83`
- `hold_gene_model_or_structure=2`

完整催化 motif、结构和系统发育仍未覆盖全库，因此不能把这些亚型标签升级为实验确认或完整机制判定。
