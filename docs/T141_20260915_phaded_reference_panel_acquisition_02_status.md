# T141 PhaDED reference-only panel amendment 02

日期：2026-09-15  
Run：`20260915_phaded_reference_panel_acquisition_02`

## 本轮完成

本轮把已有 DED ledger 中的 586 条 accession 与 NCBI Protein GenPept 快照进行 accession 精确/无版本绑定，并将结果写入新的 dated run。全部 586 条记录均绑定到既有 ledger family；没有从 parent-superfamily 推断 family，也没有创建新 family 判定。

GenPept 输入为 4 个批次、共 522 个 accession 快照，包含 accession/version、definition 和可解析的 PMID 元数据。annotation-only 记录保留为挑战/背景证据，明确标记为 `not_eligible_annotation_only`，不转写为 formal negative。已有实验阳性如来自训练 accession 或同一 publication lineage，也不计作独立 held-out 阳性。

## 校准结果

- reference-only profile：37（33 family、4 superfamily）
- 面板记录：586；已绑定既有 family：586
- formal family-resolved negative：0
- 独立留出阳性：0
- `calibrated_candidate_model`：0
- `reference_only_insufficient_panel`：33
- `planned_not_run_superfamily_requires_family_resolution`：4
- `new_family_call_created`：`false`

校准门槛仍要求每个 family 至少 3 个独立阳性、1 个 held-out 阳性、1 个 family-resolved negative、1 个 challenge control，并通过泄漏/未解释命中检查。本轮仅有部分 challenge/background 记录，缺少 formal negative 和独立留出阳性，因此所有 profile 继续保持阻断状态。

## 产物

- [panel_evidence_amendment.tsv](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/panel_evidence_amendment.tsv)
- [profile_calibration_readiness.tsv](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/profile_calibration_readiness.tsv)
- [panel_acquisition_amendment_report.json](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/panel_acquisition_amendment_report.json)
- [leaveout_calibration_decisions.tsv](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_02/results/calibration/leaveout_calibration_decisions.tsv)
- [leaveout_calibration_report.json](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_02/results/calibration/leaveout_calibration_report.json)
- [input_contract.json](/<REPO_ROOT>/runs/20260915_phaded_reference_panel_acquisition_02/input_contract.json)

## 服务器与边界

SSH 预检确认 `<SERVER_HOST>` 有 80 logical CPU 和 2 张 NVIDIA GeForce RTX 4090，预检时 GPU 利用率为 0%。本轮是 accession/文献 provenance 整理，不需要重型计算，因此没有启动服务器任务，使用线程数为 0；项目单任务上限 40 仍然有效。

这些结果只支持序列层面的候选证据和面板充分性判断，不等于 PHB/PHA 降解实验阳性。37 个 reference-only profile 在独立面板和留出校准完成前不得用于新增 family 判定、训练新 HMM 或更新 family registry。
