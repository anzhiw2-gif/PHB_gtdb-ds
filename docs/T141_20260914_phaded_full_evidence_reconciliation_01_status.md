# T141 PhaDED 全库证据分级对账

日期：2026-09-14  
运行：`20260914_phaded_full_evidence_reconciliation_01`  
状态：已完成（本地只读合并；未启动服务器计算、湿实验或 registry 修改）

## 本轮完成

- 将 `109,087` 条 PhaDED assignment 与全库 Pfam/InterPro 架构、序列完整性、SignalP 字段合并。
- 将 exact accession 文献桥接状态并入；本候选集 exact Tier A 文献支持为 `0`。
- 为每条候选生成 `evidence_grade`、`evidence_basis`、`structure_evidence_status`、`phylogeny_evidence_status` 和 `literature_evidence_tier`。
- 结构、系统发育未覆盖全库的记录显式标记 `not_tested_full_library`，不解释为阴性。

## 分级统计

| 证据等级 | 数量 | 含义 |
|---|---:|---|
| `HOLD_sequence_integrity` | 107,617 | 序列完整性警告；质量 hold，不是生物学阴性 |
| `L1_profile_plus_partial_architecture` | 779 | PhaDED profile + 部分架构证据 |
| `L0_profile_only` | 137 | 仅 profile 证据，独立架构仍 pending |
| `HOLD_architecture_conflict` | 83 | 明确 Pfam 架构冲突，保留但暂不升级 |
| `PENDING_profile_unassigned` | 471 | 当前阈值下未分配，仍保留为 unresolved candidate |

总数：`109,087`。所有记录仍为 candidate-only；没有任何等级等同于已验证 PHB 降解表型。

## 产物

- `runs/20260914_phaded_full_evidence_reconciliation_01/results/phaded_full_evidence_ledger.tsv`
- `runs/20260914_phaded_full_evidence_reconciliation_01/results/evidence_tier_manifest.json`
- `runs/20260914_phaded_full_evidence_reconciliation_01/results/evidence_grade_counts.tsv`
- `runs/20260914_phaded_full_evidence_reconciliation_01/logs/server_resource_check_20260914.txt`

## 服务器前置检查

`<SERVER_USER>@<SERVER_HOST>`（`grius`）：80 logical CPU、约 982 GiB 可用内存、2× RTX 4090（0% GPU，约 1 MiB 显存/卡）、负载约 1.26、磁盘可用约 100 TiB。由于本轮为本地合并，未在服务器启动新任务。

## 解释边界

Pfam/InterPro、HMM、SignalP、结构和系统发育证据表示同源性或功能潜力。`HOLD_sequence_integrity` 不代表生物学淘汰；`PENDING_profile_unassigned` 不代表阴性。下一步可在资源检查后，对需要的全库 InterPro/结构/系统发育覆盖任务创建 dated deploy，但不得直接运行服务器旧脚本。
