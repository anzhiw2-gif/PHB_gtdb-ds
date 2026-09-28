# PhaDED 亚型证据矩阵与 InterPro 合并实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 在不覆盖历史运行和不宣称实验表型的前提下，生成 accession 级 InterPro 状态、完整 motif 检查状态和 PhaDED 亚型证据矩阵。

**Architecture:** 新建 dated run `20260915_phaded_subtype_reconciliation_01`。以现有 assignment、校正 motif/SignalP、Pfam architecture、全库 InterPro、19 条结构审查和 42 条系统发育审查为只读输入，输出一条 accession 一行的证据矩阵、冲突摘要、profile gap 报告和 manifest。缺失或未校准证据保留为 `pending/not_tested`，不升级为阴性。

**Tech Stack:** Python 标准库、TSV/JSON、HMMER 3.4 既有产物；服务器任务仅使用 dated deploy，单个任务不超过项目规则上限。

## Global Constraints

- 不修改历史 `runs/`、registry、HMM 或正式 GTDB 扫描结果。
- 所有输出保持 `candidate-only`，不把域、motif、定位、结构或树证据写成实验阳性。
- 每个新 run 保存 `logs/`、`inputs/`、`results/` 和 `input_contract.json`。
- InterPro 状态只能是 `interpro_supported`、`interpro_not_reported` 或 `interpro_tool_input_excluded`。
- 未训练 profile 只能生成 gap/`planned_not_run` 记录，不能降低阈值强行归类。

### Task 1: InterPro accession 汇总

**Files:**
- Create: `pipeline/scripts/merge_phaded_interpro_evidence.py`
- Test: `pipeline/tests/test_merge_phaded_interpro_evidence.py`

输出每个 accession 一行，聚合分析库、signature、InterPro ID、命中行数，并显式区分工具排除和未报告。

### Task 2: Motif 与亚型证据矩阵

**Files:**
- Create: `pipeline/scripts/build_phaded_subtype_matrix.py`
- Test: `pipeline/tests/test_build_phaded_subtype_matrix.py`

保留 lipase box、catalytic Ser/Cys、His/Asp、oxyanion hole、SBD/linker/lid 的独立状态；当前未执行的字段写入 `not_tested_full_library` 或 `pending_not_implemented`。合并 profile、domain、motif、localization、structure、phylogeny 状态，生成 `subtype_call` 与 `subtype_confidence`。

### Task 3: 新 dated run 与 profile gap

创建 `runs/20260915_phaded_subtype_reconciliation_01/`，生成输入契约、结果哈希、冲突/重点候选摘要和未训练 profile gap。profile 补建本轮只生成 `planned_not_run` 计划，直到有独立阳性、近邻阴性和 challenge set 后再执行校准。

### Task 4: 验证与报告

运行相关测试、全套测试、`compileall` 和 `git diff --check`，更新中文状态报告并明确候选证据边界。
