# PhaDED 全库证据补齐实施计划

> **For agentic workers:** Use the existing candidate-only workflow and create a new dated run for every server execution. Do not modify historical runs or registry files.

**Goal:** 在不进行湿实验和正式 GTDB 重扫的前提下，补齐 109,087 条 PhaDED 候选的 InterPro/序列证据、分层统计和可追溯证据等级。

**Architecture:** 当前 InterProScan 全库任务独立完成域注释；完成后以 accession 为唯一连接键，将 InterPro、Pfam、序列完整性、SignalP、结构/系统发育覆盖状态和 exact accession 文献状态合并到新的候选台账。结构和系统发育只对冻结后的代表序列/重点候选执行，不对全库盲目预测。

**Tech Stack:** InterProScan 5.78-109.0、Python 标准库、TSV/JSON manifest、T141 dated `deploy/<run_id>/`。

## Global Constraints

- 所有结果保持 `candidate-only`，不把域、结构、树或同源性当作 PHB 降解表型。
- 不修改 `pipeline/config/`、历史 HMM、历史 `runs/` 或正式 registry。
- 服务器执行只能来自匹配的 dated `deploy/<run_id>/`，单任务 CPU 按当前项目约束不超过 40。
- 末端 `*` 只作为 `terminal_stop_only` 记录；内部异常字符才进入序列完整性 hold。
- 缺少 InterPro、结构、系统发育或文献证据必须标记 `pending/not_tested/unresolved`，不能解释为阴性。
- 每个输出保存路径、版本、大小和 SHA-256；不清理历史失败或工具排除记录。

## Task 1: Monitor the active InterProScan run

**Inputs:** `runs/20260914_phaded_interpro_full_01/input_contract.json`, server PID `2617348`.

**Outputs:** `results/interpro_full.tsv`, `results/interpro_full_run_manifest.json`, stdout/stderr logs.

- [ ] 每次检查记录 PID、Java 子进程、临时目录大小、最终 TSV 是否出现、stderr 是否非空、服务器 CPU/内存/GPU。
- [ ] 只要主进程和 Java 子进程存活，不重启、不重复提交。
- [ ] 完成后校验输出存在、输入 SHA-256 匹配、序列访问号覆盖数和 manifest SHA-256。

## Task 2: Prepare the InterPro parser locally while Task 1 runs

**Files:**
- Create: `pipeline/scripts/merge_phaded_interpro_evidence.py`
- Test: `pipeline/tests/test_merge_phaded_interpro_evidence.py`

**Consumes:** InterPro TSV schema、PhaDED assignment、Pfam feature table。

**Produces:** 每个 accession 一行的 InterPro 状态表；不命中的候选明确为 `not_reported`，工具排除明确为 `tool_input_excluded`。

- [ ] 先写失败测试：重复 accession、缺失 query、同一候选多个 InterPro 行聚合、空输出和工具排除状态。
- [ ] 实现按 accession 聚合 InterPro member/database、InterPro entry、GO/IPR lookup 字段；不从自由文本推断 PHB 活性。
- [ ] 运行该测试、`compileall` 和 `git diff --check`。

## Task 3: Validate and freeze the completed InterPro output

**Inputs:** Task 1 output、`candidate_accessions.txt`、`interpro_excluded.tsv`。

**Outputs:** `runs/20260914_phaded_interpro_full_01/results/interpro_validation.json`。

- [ ] 核对规范化 FASTA 的 `108,736` 条输入、`351` 条工具排除与原始 `109,087` 条候选总数完全对账。
- [ ] 检查 TSV 中 accession 是否重复、是否出现输入集合外 query、是否有截断行或异常列数。
- [ ] 将缺失结果标为 `not_reported`，不能标为 `negative`。
- [ ] 只有验证通过后才允许进入 Task 4。

## Task 4: Merge the complete candidate evidence ledger

**Inputs:** Task 3 validated InterPro table、`phaded_assignment.tsv`、`phaded_feature_consistency.tsv`、exact accession bridge、规范化序列完整性结果。

**Outputs:** 新 dated run 中的 `phaded_full_evidence_ledger.tsv` 和 manifest。

- [ ] accession 集合必须完全一致；重复或缺失直接失败。
- [ ] 每条记录保留 PhaDED family/superfamily、Pfam 架构、InterPro 状态、SignalP、motif、序列完整性、文献状态。
- [ ] 增加 `structure_evidence_status`、`phylogeny_evidence_status`、`interpro_evidence_status`。
- [ ] 证据等级建议保持分层：`HOLD_architecture_conflict`、`HOLD_sequence_integrity`、`L2_profile_plus_architecture_consistent`、`L1_profile_plus_partial_architecture`、`L0_profile_only`、`PENDING_profile_unassigned`。
- [ ] 不删除任何候选；hold 是审查状态，不是生物学阴性。

## Task 5: Generate stratified statistics and candidate exports

**Consumes:** Task 4 frozen ledger、GTDB R232 taxonomy/metadata。

**Produces:** protein-level、genome-level、family/superfamily-level统计，及按证据等级的候选清单。

- [ ] 分开统计 protein、genome、genome × family、genome × superfamily 分母。
- [ ] 分别报告 assigned、ambiguous_family、ambiguous_superfamily、unassigned_PhaDED_like。
- [ ] 按 evidence grade、Pfam architecture、InterPro state、SignalP、序列完整性和 GTDB 门/纲分层。
- [ ] 所有表格写明“candidate-only”和分母定义。

## Task 6: Run independent evidence branches after Task 4

这些分支可并行，但每个分支必须使用自己的 dated run/输出目录，不能同时写 Task 4 的权威台账。

- [ ] **文献 accession bridge：** 使用完整 GTDB R232 metadata/protein accession 表复核 exact genome/protein/gene accession；taxonomy-only 保持非 Tier A。
- [ ] **结构分支：** 仅对 19 条人工重点候选和冻结的代表序列做结构/Foldseek；3 条架构冲突和 2 条 gene-model hold 不能升级为阳性。
- [ ] **系统发育分支：** 按 PhaDED family/superfamily 建代表序列树，包含参考、challenge、ambiguous 和跨分类代表；不对 109,087 条全量建树。
- [ ] **生态/分类分支：** 依据 Task 5 统计表生成 GTDB 分类分布和候选优先级。

服务器资源策略：结构和系统发育若同时运行，每个任务不超过 20 CPU；默认优先串行运行两个重任务。GPU0 当前有其他用户作业，不能假定可用。

## Task 7: Final reconciliation and reporting

- [ ] 合并 Task 6 的分支结果时只增加证据字段，不覆盖 PhaDED 原始 assignment。
- [ ] 更新 dated run manifest、哈希清单和中文状态报告。
- [ ] 运行 `python -m unittest discover -s pipeline/tests -v`、`python -m compileall -q pipeline` 和 `git diff --check`。
- [ ] 报告中明确：序列/域/结构/树/文献桥接是证据等级，不是实验阳性；缺失证据不等于阴性。

## Parallelism Summary

| 阶段 | 可并行任务 | 必须串行的边界 |
|---|---|---|
| InterProScan 运行中 | 监控、解析器测试、文献表准备、文档准备 | 不启动第二个全库 InterPro/Pfam/HMM 任务 |
| InterPro 完成后 | 输出校验可与文献 bridge 准备并行 | 校验通过前不能合并总台账 |
| 总台账冻结后 | 分类统计、文献 bridge、代表结构、代表树可分支 | 每个分支独立输出；最终合并串行 |
| 最终收口 | 仅报告、manifest、哈希和测试 | 不再修改历史 run 或 registry |
