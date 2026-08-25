# Run Isolation and Provenance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为本地 pipeline 增加 dated run context、GTDB 输入契约和项目规则入口，避免新的 manifest/log 覆盖历史记录。

**Architecture:** 新增纯 Python `run_context.py`，负责生成/校验 `runs/<run_id>/` 下的 `manifest.jsonl`、`run_manifest.json`、日志目录和输入契约；`run_pipeline.sh` 接受 `--run-id/--run-dir`，默认创建 UTC dated 目录，并把现有根目录 data 作为显式只读输入。完整数据目录迁移暂不在本轮强行重构，后续通过部署适配器逐步完成。

**Tech Stack:** Bash、Python 3.12、unittest、JSON/TSV SHA-256 manifests。

## Global Constraints

- 不启动 T141 重算，不安装软件，不删除历史结果，不提交或推送。
- 服务器旧脚本与本地脚本继续分开；不得把本地路径契约直接套到 T141 扁平目录。
- GTDB taxonomy/metadata/tree 缺失时必须标记 `pending`，不能伪造哈希或写成已绑定。
- 生产代码先有失败测试，再做最小实现。

### Task 1: Project Rules and Input Contract

**Files:** `AGENTS.md`, `pipeline/scripts/run_context.py`, `pipeline/tests/test_run_context.py`

- [ ] 创建 60 行以内项目规则入口，固定运行、科学措辞、服务器、清理和提交边界。
- [ ] 测试 run id 校验、目录创建、必需输入哈希和缺失 GTDB 输入的 `pending` 状态。
- [ ] 实现 JSON input contract，区分 `verified`、`pending`、`missing`。

### Task 2: Dated Run Manifest and Logs

**Files:** `pipeline/scripts/run_pipeline.sh`, `pipeline/scripts/run_manifest.py`, `pipeline/tests/test_audit_guards.py`

- [ ] 增加 `--run-id`、`--run-dir`，默认 `runs/<UTC timestamp>_<short git>/`。
- [ ] 将 JSONL、最终 manifest、主日志和输入契约写入 run directory，旧 `results/run_manifest*` 不再被覆盖。
- [ ] 严格 provenance 记录 run id、run directory、GTDB input contract 和 source bundle。
- [ ] 保留兼容的显式 `--legacy-root-results`，默认关闭且仅供审计旧流程。

### Task 3: Documentation and Handoff

**Files:** `docs/STATUS.md`, `docs/reproducibility.md`, `README.md`

- [ ] 记录新的 run directory 命令、目录结构和当前“根 data 只读、派生输出渐进迁移”的边界。
- [ ] 将 GTDB taxonomy/metadata/tree 未绑定标为 `pending`。
- [ ] 记录本轮未执行服务器部署和重算。

### Task 4: Verification

- [ ] 运行新增测试和全套 `python -m unittest discover -s pipeline/tests -v`。
- [ ] 运行 `compileall`、`git diff --check`；Linux CI 负责 `bash -n`。
