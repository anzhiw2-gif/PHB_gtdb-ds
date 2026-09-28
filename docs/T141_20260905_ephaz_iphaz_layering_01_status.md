# T141 20260905 ePhaZ/iPhaZ candidate-only 分层状态

## 状态

- run_id: `20260905_ephaz_iphaz_layering_01`
- 状态: `completed_candidate_only`
- 服务器: `<SERVER_USER>@<SERVER_HOST>`
- 执行解释器: `${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb/bin/python`
- 正式扫描: 未启动
- `pipeline/config/formal_scan_models.tsv`: 未修改

## 输入与源码

本轮固定使用 run 输入中的 ePhaZ/iPhaZ tier1、tier2 FASTA、ePhaZ SignalP 子集和竞争审计表。服务器执行的修复脚本 SHA-256 为
`7410183b08893909fb866f239e8d3f98b73c6c6b846dcfab16fafb4df7fdade8`。初始脚本快照 SHA-256 为
`79e06afdfbc58bf6f1f3c5b7258a960ed4a9ecceba5986bd980123423bb3f63d`，两者均保留。

## 输出

`results/layers_retry_01/` 保留为失败证据：原实现对 tier1/tier2 重叠 accession 重复输出，180,124 行中仅 109,087 个唯一 accession。

`results/layers_retry_02/` 为接受的 candidate-only 输出。修复逻辑以 `(family, accession)` 去重，优先级为 `competition_review > curated/tier1 > tier2 > broad/pending`。

| 层 | 蛋白数 | 基因组数 |
|---|---:|---:|
| `ePhaZ_competition_review` | 18 | 18 |
| `ePhaZ_curated_core_nonsecreted` | 16,818 | 13,213 |
| `ePhaZ_curated_core_secreted` | 21,856 | 17,223 |
| `ePhaZ_tier2_review` | 31,990 | 25,763 |
| `iPhaZ_tier1` | 32,926 | 25,920 |
| `iPhaZ_tier2_review` | 5,479 | 5,153 |
| 合计 | 109,087 | 76,141 `(genome,family)` |

去重前 180,124 行，移除 71,037 条重复来源；去重后 `(family, accession)` 重复数为 0。蛋白表与基因组表的键集合、层集合和 `protein_count` 已独立复算一致，格式错误数为 0。

## 解释边界

这些层只表示 HMM tier、SignalP 和竞争审计形成的候选同源/功能潜力。它们不等同于实验确认的 PHB 降解表型，也不应直接作为物种生态或底物特异性结论。ePhaZ curated secreted 与 nonsecreted 仍需分开报告；iPhaZ tier1 与 tier2 review 不能合并为同等证据层。

## 结果哈希

- `protein_layers.tsv`: `f4a9d99039ef5dc46fbea5b90dd98fd61adc2c65ae1da700c53976cb5b7df08c`
- `genome_layers.tsv`: `6eeee505d76c47127ca64cae9ce5d5b96028c068ddb60d705063a6974fa46fe5`
- `stratification_metadata.json`: `c108d480096e1cae8e3ce7fe2c58e9af1dabbd89375334137bbe5fb1b355961e`

## 后续边界

当前结果可用于候选分层、抽查和下游统计；在没有额外授权前，不注册新模型、不覆盖正式 registry、不启动新的 GTDB 扫描，也不删除 retry-01 失败证据。
