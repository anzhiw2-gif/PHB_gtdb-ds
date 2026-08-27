# 审计一致性修复 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复 PHB_gtdb-ds 的失败传播、运行/树审计、结果 schema、文档状态和公开仓库卫生，使本地版本成为可验证的发布候选版本。

**Architecture:** 保留现有 Bash 编排和 Python 结果脚本，只增加集中式门禁与验证函数。运行清单记录每步与输入输出哈希；树清单根据当前 tier FASTA 判定历史输入；文档只引用经过验证的结果表。服务器部署不在本计划执行。

**Tech Stack:** Bash、Python 3.12、`unittest`/现有 smoke tests、GitHub Actions、PowerShell 只读核验。

## Global Constraints

- 不使用 `git reset`、`git checkout` 或清理命令覆盖已有未提交修改。
- 不启动 T141 长任务，不修改服务器原始数据。
- 不提交、不推送，除非用户另行明确授权。
- HMM、SignalP、基序、邻域和树只能表述为候选/功能潜力，不能表述为已证实表型。
- 所有新增行为先写失败测试，再写最小实现。

---

### Task 1: 运行门禁测试

**Files:**
- Create: `pipeline/tests/test_audit_guards.py`
- Modify: none initially

**Interfaces:**
- Tests will import small pure helpers from `run_manifest.py` and `09i_tree_manifest.py` after those helpers are introduced.

- [ ] 写测试：缺失输入/输出返回清晰错误；树输入哈希变化返回 `stale_input`；cluster summary schema 必须包含五列。
- [ ] 运行 `python -m unittest pipeline.tests.test_audit_guards -v`，确认因 helper 尚不存在而失败。

### Task 2: 修正 run manifest

**Files:**
- Modify: `pipeline/scripts/run_manifest.py`
- Modify: `pipeline/scripts/run_pipeline.sh`
- Test: `pipeline/tests/test_audit_guards.py`

**Interfaces:**
- `run_manifest.validate_paths(paths, label)` raises `ManifestError` on missing/empty paths.
- `run_manifest.build_manifest(...)` returns a JSON-serializable dictionary.

- [ ] 实现最小路径校验和 manifest builder。
- [ ] 让 finalize 对必需输入/输出缺失返回非零。
- [ ] 将主管线改为 `set -Eeuo pipefail`，步骤失败立即退出，并在 finalize 成功后才输出完成日志。
- [ ] 运行对应 unittest 和现有 smoke test。

### Task 3: 修正 HMM 输入门禁

**Files:**
- Modify: `pipeline/scripts/06_screen.sh`
- Modify: `pipeline/scripts/06b_aggregate_hits.py`
- Test: `pipeline/tests/test_audit_guards.py`

**Interfaces:**
- `06b_aggregate_hits.validate_pairs(tbl_files)` rejects missing `.dom` companions.

- [ ] 写 fixture 测试缺 HMM/缺 dom 时失败。
- [ ] 在 `06_screen.sh` 中对声明家族的 HMM 缺失直接退出。
- [ ] 在聚合器中校验配对文件和非空输入，保留明确的无命中状态。
- [ ] 运行测试。

### Task 4: 修正树清单与 cluster schema

**Files:**
- Modify: `pipeline/scripts/09i_tree_manifest.py`
- Modify: `pipeline/scripts/11_clusters.py`
- Modify: `results/trees_tier1/tree_manifest.tsv`
- Modify: `results/tables/cluster_summary.tsv`
- Test: `pipeline/tests/test_audit_guards.py`

**Interfaces:**
- `09i_tree_manifest.classify_status(n_leaves, tier_count, input_sha, current_sha)` returns `ok`, `stale_input`, or `EMPTY_OR_PARSE_FAIL`.
- `11_clusters.write_summary(rows, path)` writes `marker_hits`, `supporting_loci`, `supporting_genomes` columns.

- [ ] 让旧 OH 树被标记 `stale_input`，不删除历史树。
- [ ] 统一 cluster summary writer 和表头。
- [ ] 测试最终 tier 计数和树清单状态。

### Task 5: 文档和仓库卫生

**Files:**
- Modify: `.gitignore`
- Modify: `README.md`, `docs/STATUS.md`, `docs/项目审核报告.md`, `docs/项目详细审查报告.md`, `docs/reproducibility.md`
- Modify: `.github/workflows/ci.yml`

- [ ] 把 `.playwright-cli/` 加入忽略规则，并保留本地文件不删除。
- [ ] 明确报告是 local working-tree 状态，未发布内容不得写成 GitHub 已发布。
- [ ] 修正 SignalP checkbox、树旧输入状态、cluster 计数定义和核心共存数字。
- [ ] CI 增加审计测试。

### Task 6: 本地验证与交付边界

- [ ] 运行 Python AST/compile、unittest、smoke、`git diff --check`。
- [ ] 检查 GitHub HEAD 与本地差异，生成服务器部署清单，不执行服务器写入。
- [ ] 请求独立代码审查，记录 Critical/Important findings。
