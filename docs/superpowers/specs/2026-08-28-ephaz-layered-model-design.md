# ePhaZ 分层模型设计

## 目标

将现有宽 ePhaZ 种子库拆成两个互斥、可追溯的层：

- `ePhaZ_curated_core`：有明确实验支持、长度/完整性通过且结构域架构典型的序列。
- `ePhaZ_broad_discovery`：远缘、注释型、架构待确认或短序列，用于探索性搜索。

所有当前 `<200 aa` 的 61 条序列必须逐条出现在审查表中；长度本身不能触发删除。

## 数据流

`inputs/seeds_clean/ePhaZ.faa` 与 v2 manifest 进入 `split_ephaz_seeds.py`。脚本按 accession 唯一化并执行 fail-closed 校验，输出两份 FASTA、分层 manifest、61 条短序列审查 TSV、输入/输出 SHA-256。

只有 manifest 中有显式实验支持字段、长度不少于 200 aa、完整性与架构状态均为 `pass` 的 accession 才能进入 core；其他合法序列进入 broad，理由分别记录为 `short_sequence`、`evidence_missing`、`architecture_pending` 或 `remote_annotation`。

HMM builder 接受显式 family 列表，保留原九家族兼容性；新 run 显式传入两个 ePhaZ 层并分别输出 alignment 与 HMM。旧 `ePhaZ.hmm` 不参与新层校准。

校准器按层读取各自 HMM。core 使用 core 正对照，broad 使用全部 ePhaZ 正对照并保留每条对照的 evidence level；两个结果写入同一张表但 family 名称严格区分。

## 证据边界与错误处理

- `reviewed=true` 不是实验表型证据；实验支持必须来自 manifest 的 evidence/文献字段。
- HMM、domain coverage 和 SignalP 只代表候选同源或功能潜力。
- 重复 accession、非法氨基酸、空序列、manifest 缺项、缺失 HMM 或缺失对照均 fail-closed。
- 不删除历史输入或结果；新输出只能位于 dated `runs/<run_id>/`。

## 验收标准

1. 两层 accession 互不相交，且并集覆盖原始 ePhaZ accession。
2. 61 条 `<200 aa` 全部出现在 broad 或 `review_required`，无一自动进入 core。
3. 两个 HMM 可由 `hmmstat` 读取，alignment/HMM 均有 SHA-256。
4. 校准输出含 `ePhaZ_curated_core` 和 `ePhaZ_broad_discovery`，记录 TP/FP/TN/FN、覆盖度、阈值、HMM/输入哈希和工具版本。
5. 相关单元测试、`compileall`、`git diff --check` 通过；不启动全量 268G 筛选。
