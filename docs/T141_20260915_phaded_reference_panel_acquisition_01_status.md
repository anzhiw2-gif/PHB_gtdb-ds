# T141 PhaDED reference-only panel acquisition and leave-out calibration

日期：2026-09-15  
Run：`20260915_phaded_reference_panel_acquisition_01`

## 结果

本轮针对 profile manifest 中的 37 个 `reference_only` profile 建立了 111 行 acquisition 计划，每个 profile 各有 `independent_positive`、`family_resolved_negative` 和 `challenge_control` 三个角色。计划清单中的 accession 留空，表示尚未取得可绑定该 profile 的记录；没有用 parent-superfamily 结果推断 family。

复用并绑定了既有外部 panel 的 17 条 accession 记录（5 条独立实验 PHB 阳性、1 条 MCL-PHA 阳性/底物对照、4 条定位或底物阴性候选、8 条 challenge/annotation/fragment 控制），并复制了 33 个 DED family 的 132 个 source snapshot 文件。外部 panel 的 response 文件、FASTA、序列哈希、PMID/DOI 和旧 run provenance 均保存在本 run 的 `inputs/` 下。

严格 accession 绑定后：

- 绑定到 37 个 `reference_only` profile：`0`
- 未解析或仅属于既有 trained family：`17`
- `reference_only_insufficient_panel`：`33` 个 family profile
- `planned_not_run_superfamily_requires_family_resolution`：`4` 个 superfamily profile
- `calibrated_candidate_model`：`0`
- 新增 family 判定：`false`

## 校准门槛

只有同时满足以下条件才允许进入 family-specific 留出校准：至少 3 个独立阳性 accession、至少 1 个 held-out 阳性、至少 1 个 family-resolved negative、至少 1 个 challenge control、无未解释 negative/challenge 命中，并通过 accession/序列哈希和同源泄漏检查。本轮没有任何 profile 满足全部条件，因此没有执行新 HMM 训练、全库扫描或 family registry 更新。

## 服务器检查

执行前检查 `<SERVER_HOST>`：80 logical CPU；2 张 NVIDIA GeForce RTX 4090，GPU 利用率 0%。本轮不需要服务器计算，未启动服务器任务；项目硬规则中的单任务线程上限 40 仍适用。

## 主要产物

- [reference_panel_acquisition.tsv](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_01/results/reference_panel_acquisition.tsv)
- [panel_binding.tsv](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_01/results/panel_binding.tsv)
- [profile_calibration_readiness.tsv](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_01/results/profile_calibration_readiness.tsv)
- [leaveout_calibration_decisions.tsv](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_01/results/calibration/leaveout_calibration_decisions.tsv)
- [leaveout_calibration_report.json](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_01/results/calibration/leaveout_calibration_report.json)
- [input_contract.json](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_01/input_contract.json)

## 结论边界

当前 37 个 profile 仍是 `reference_only`，不能用于新增 family 判定。DED/InterPro/domain/motif/structure/phylogeny 和文献桥接只表示候选同源或功能潜力；本 run 不声明任何 GTDB 蛋白已验证 PHB/PHA 降解表型。
