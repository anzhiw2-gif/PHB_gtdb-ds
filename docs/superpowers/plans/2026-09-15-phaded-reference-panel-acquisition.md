# PhaDED Reference Panel Acquisition and Leave-Out Calibration Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 为 37 个 `reference_only` PhaDED profile 建立可审计的独立阳性、family-resolved 阴性和 challenge panel，并仅对满足独立性、泄漏控制和留出条件的 profile 运行校准。

**Architecture:** 新建 dated run `20260915_phaded_reference_panel_acquisition_01`，先由 profile manifest 生成三角色 acquisition manifest，再通过 UniProt/NCBI/Europe PMC/PubMed/PhaDED 页面补充有明确 accession 和文献证据的记录。所有无法绑定到现有 family 定义的记录保持 `unresolved`；校准脚本只消费已绑定、去重、无同源泄漏的 panel。

**Tech Stack:** Python 标准库、TSV/JSON、NCBI E-utilities、Europe PMC REST、UniProt REST、PhaDED/DED source pages、HMMER/MAFFT（仅在 panel gate 通过后）。

## Global Constraints

- 新输出只能写入 `runs/20260915_phaded_reference_panel_acquisition_01/`，不得覆盖历史运行。
- 服务器仅执行 dated `deploy/<run_id>/` 中绑定源码；单任务线程不超过 `AGENTS.md` 的 40。
- 37 个 reference-only profile 在 panel gate 前保持 `reference_only`/`planned_not_run`，不得新增 family 判定。
- Homology、domain、motif、structure 和 phylogeny 只能作为 candidate evidence，不等同于 PHB/PHA 表型证明。
- 每条外部记录保存 accession、序列哈希、来源、版本、检索日期和 DOI/PMID/PMCID；缺失项写 `pending`。

### Task 1: Acquisition manifest

Create `pipeline/scripts/prepare_phaded_reference_panel.py` and its test. Read the authoritative profile manifest, select exactly 37 reference-only profiles, and emit one planned row for each role (`independent_positive`, `family_resolved_negative`, `challenge_control`) with no invented accession.

### Task 2: External evidence acquisition

Create a dated run with `logs/`, `inputs/`, `results/`, and `input_contract.json`. Query authoritative protein and literature APIs using bounded requests. Import existing independently sourced PHA/PhaZ controls only when the source record has explicit assay/result and preserve `family_binding_status=unresolved` unless the accession is already bound in the reference ledger.

### Task 3: Panel integrity and leakage audit

Implement validation for required roles, accession and sequence-hash deduplication, source independence, family binding, training/held-out separation, and explicit negative/challenge distinction. Report unresolved and insufficient rows rather than treating missing evidence as negative.

### Task 4: Leave-out calibration

Run family-specific calibration only for profiles with at least three independent positive accessions, at least one held-out positive, family-resolved negatives and challenge controls with zero unexplained negative hits. Emit `calibrated_candidate_model` only when all gates pass; otherwise emit `reference_only_insufficient_panel` or `planned_not_run` and leave the authoritative family registry unchanged.

### Task 5: Verification and handoff

Run focused tests, full test discovery, `compileall`, and `git diff --check`. Write a Chinese status report with counts, provenance hashes, unresolved gaps, server resource checks, and the explicit boundary that no new family calls were made.
