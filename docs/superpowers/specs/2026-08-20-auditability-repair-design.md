# PHB_gtdb-ds 审计一致性修复设计

## 目标

把本地 `PHB_gtdb-ds` 工作树整理成一个可验证的发布候选版本：主管线失败必须失败退出，运行清单和树清单必须真实存在且绑定输入，报告数字必须与结果表一致，公开仓库不再携带浏览器/会话残留。

本阶段不启动 T141 长任务，不覆盖服务器原始结果，不提交或推送。

## 范围

### 运行门禁

- `pipeline/scripts/run_pipeline.sh` 使用 fail-closed 错误传播。
- 08、09、10、11、manifest 固化的非零退出、空结果和缺少关键产物都必须中止。
- `06_screen.sh` 缺少声明的 HMM 时失败；`06b_aggregate_hits.py` 校验每个 tbl/dom 配对及命中文件完整性。

### 审计一致性

- `run_manifest.py` 记录命令、git commit、输入/输出 SHA-256，并对缺失输入/输出返回非零。
- `09i_tree_manifest.py` 记录实际叶数、tier FASTA 计数和哈希；当树输入不匹配当前 tier1 时标记 `stale_input`。
- 树清单保留历史树，但不得把历史输入描述为当前最终树。

### 结果与文档

- `cluster_summary.tsv` 统一为 `marker_hits/supporting_loci/supporting_genomes` schema。
- 文档区分 `local-uncommitted`、`server-runtime` 和 `github-published` 状态。
- 明确 `44,814` 为四个核心家族去重基因组数；含 patatin 的 `30,058` 多家族共存不得冒充核心家族数字。

### 仓库卫生

- `.playwright-cli/` 进入忽略规则；会话文件从版本控制中移除但保留本地文件。
- CI 覆盖 Python 编译、Bash 语法和 smoke tests。
- 不在本阶段生成或伪造 `environment.lock.yml`、GTDB 校验和或服务器运行清单。

## 数据流

```text
inputs + scripts + parameters
        |
        v
run_pipeline.sh -- each step checked
        |
        +--> run_manifest.jsonl -- finalize --> run_manifest.json
        |
        +--> final tier FASTA --> tree tools --> tree_manifest.tsv
        |
        +--> cluster_context.tsv --> cluster_summary.tsv
```

任何箭头上的必需输入缺失、输出为空、schema 不符合或命令失败，都停止后续流程。

## 错误处理

- 管线使用 `set -Eeuo pipefail`，每个步骤通过统一包装器记录退出码。
- 允许的“无命中”只能是显式、已记录的分析结果；不能把缺文件或解析失败当作无命中。
- manifest 固化失败时不得打印完成日志。
- 树清单对输入 FASTA 变化使用 `stale_input` 状态，不自动重写树文件。

## 测试策略

- 先为每个修复写最小失败测试：命令失败传播、缺 HMM、缺 manifest 输入、stale tree、cluster schema。
- 使用现有无 HMMER 的 smoke test，补充纯 Python fixture 测试。
- 本地执行 Python AST/compile、smoke；Bash 语法在有 Bash 的环境或 CI 执行。

## 后续边界

本阶段完成后，另行准备 T141 部署包和重算命令；服务器执行需再次核对脚本 SHA-256、输入清单和输出状态，不能直接把本地文档声明视为服务器证据。
