# ePhaZ Ambiguous Stratified Review Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 对 6,991 条 ePhaZ-vs-iPhaZ `ambiguous` 候选补算精确 domain 覆盖度，并按 bitscore 差值、长度、覆盖度和邻域证据生成可复核的固定种子分层样本。

**Architecture:** 使用独立 dated run 保存候选输入、两模型 domtblout、邻域快照、分层全表和样本表。脚本只读取竞争结果，解析精确 HMM 坐标覆盖度，按固定生物学边界形成互斥 strata，并在每个 strata 内使用确定性 SHA-256 排序抽样，避免随机状态漂移。

**Tech Stack:** Python 标准库、HMMER `hmmsearch`、TSV/FASTA、SHA-256 provenance。

## Global Constraints

- 运行必须使用新的 `runs/<run_id>/`，不得覆盖历史 `results/` 或服务器历史目录。
- 服务器只执行 dated `deploy/<run_id>/` 中绑定源码。
- 缺失 domain 或邻域证据必须写为 `pending`/`no_record`，不得推断为阴性。
- 该审计只表示候选同源、模型覆盖和邻域证据，不等同于 PHB 降解表型。
- 完成后运行相关测试、`compileall` 和 `git diff --check`。

### Task 1: Stratification script and tests

**Files:**
- Create: `pipeline/scripts/stratify_ephaz_ambiguous.py`
- Test: `pipeline/tests/test_stratify_ephaz_ambiguous.py`

- [x] Write failing tests for delta/length/coverage/neighborhood bins, deterministic per-stratum sampling, domain coverage parsing, and fail-closed missing columns.
- [x] Run the focused tests and verify the expected implementation failures.
- [x] Implement TSV readers, domtblout coverage parser, accession-to-locus neighborhood join, stratum labels, stable SHA-256 sampling, and provenance outputs.
- [x] Run focused tests until green.

### Task 2: Create and bind the T141 audit run

**Files:**
- Create: `runs/20260829_ephaz_ambiguous_sampling_01/`
- Create: `deploy/20260829_ephaz_ambiguous_sampling_01/`

- [x] Copy the exact competition TSV, ambiguous FASTA, two HMMs, SignalP reference, and neighborhood snapshots into the new run.
- [x] Copy the implementation and `run_context.py` into dated deploy.
- [x] Write `input_contract.json` with SHA-256 and leave unrelated GTDB inputs `pending`.

### Task 3: Small exact-model coverage run and stratified sampling

- [x] Run `hmmsearch` only on the 6,991 ambiguous FASTA records with the exact ePhaZ core and iPhaZ HMMs, writing domtblout under the new run.
- [x] Run the stratifier with fixed seed `20260829` and per-stratum cap `10`.
- [x] Validate counts, hashes, and no full-screening process/marker.

### Task 4: Documentation and verification

- [x] Create `docs/T141_20260829_ephaz_ambiguous_sampling_01_status.md` with strata definitions, sample counts, missing-evidence status, and interpretation boundary.
- [x] Run full tests, `compileall`, and `git diff --check`.
- [x] Do not commit or push unless separately authorized.
