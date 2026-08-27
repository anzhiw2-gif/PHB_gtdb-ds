# 审计建议落地 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将 2026-08-24 审计建议落地为可验证的本地项目基线：运行清单绑定源码/环境/输入，执行链保持 fail-closed，结果文档统一统计口径并明确树与功能证据边界。

**Architecture:** 保留现有 `pipeline/scripts` Bash/Python 结构，在 `run_manifest.py` 增加可选但默认可生成的 provenance 区块；脚本域只修复失败传播和入口约束；文档域只引用已复算的本地结果并标记服务器部署/旧版本差异。服务器历史运行目录不写入、不删除。

**Tech Stack:** Bash、Python 3.12、`unittest`、现有审计 guard、PowerShell 只读核验。

## Global Constraints

- 不执行安装、删除、配置修改、提交或推送；保留工作树中已有修改和服务器历史目录。
- HMM、domain、SignalP、邻域和树只能表述为候选同源/功能潜力，不表述为已验证 PHB 降解表型。
- OH 1,429、patatin 1,372/620、patatin 112,926 记录与四家族去重 44,814/199,923 必须分层表述。
- 新增生产行为先写回归测试并观察测试失败，再实现最小修复。

### Task 1: Manifest Provenance

**Files:** `pipeline/scripts/run_manifest.py`, `pipeline/scripts/run_pipeline.sh`, `pipeline/tests/test_audit_guards.py`

- [x] 为 manifest 增加源码包/环境/GTDB 输入/完整命令 provenance 字段及哈希 helper。
- [x] 为缺少 provenance 的严格模式增加失败测试；保持旧调用可通过显式兼容模式。
- [x] 将主编排传入实际命令、环境入口和声明的 HMM/GTDB 文件（含 06 的 v2 与 08c 的 legacy core 两套实际模型）。
- [x] 运行 manifest 相关 unittest。

### Task 2: Fail-Closed Execution Chain

**Files:** `pipeline/scripts/05_predict_proteins.sh`, `pipeline/scripts/06_screen.sh`, `pipeline/scripts/08c_tier_rescore.sh`, `pipeline/scripts/run_pipeline.sh`, focused tests

- [x] 删除会吞掉并行/工具失败的 `|| true` 和 warning-then-continue 分支。
- [x] 缺少声明输入、空输出、任务返回非零时立即失败且保留可审计日志。
- [x] 增加脚本静态回归测试并运行现有全套测试。

### Task 3: Scientific Status and Documentation

**Files:** `docs/STATUS.md`, `docs/项目最新审计报告_20260824.md`, `docs/项目详细审查报告.md`, `README.md`, relevant figure/table source notes

- [x] 统一 26,850/26,855 差异，明确 Figure 5 denominator 与 tier1/patatin 层级。
- [x] 将 OH stale input、ePhaZ CD-HIT input_not_registered、全量树/HGT 未完成写入唯一状态页。
- [x] 明确 GitHub `main`、本地工作树、T141 deploy/run 的权威关系和当前 server manifest SHA 漂移。
- [x] 检查中文科学措辞不越过 phenotype evidence boundary。

### Task 4: Integration Verification

- [x] 运行 `python -m unittest discover -s pipeline/tests -v`。
- [x] 运行 Python compile/AST 检查和 `git diff --check`。
- [x] 检查 agent 修改是否冲突，输出本地部署清单；不写入 T141。
