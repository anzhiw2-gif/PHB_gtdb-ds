# T141 PhaDED 竞争集架构与系统发育复核

日期：2026-09-12  
运行：`runs/20260911_phaded_competition_phylo_review_01/`  
状态：`completed_candidate_only`

## 输入与执行

已从权威 `priority_review_03/interpro_priority_v4` 和 `conflicting_review_dossier` 冻结 42 条去重候选：

| 面板 | 候选数 | 参考数 | 树叶数 |
|---|---:|---:|---:|
| `extracellular_explicit` | 4 | 71 | 75 |
| `phaC_vs_phaZ` | 12 | 52 | 64 |
| `structure_anomaly` | 8 | 301 | 309 |
| `ephaz_competition` | 18 | 284 | 302 |

每个面板均使用候选序列、对应 PhaDED 参考序列和独立的 challenge/参考角色标记。PhaC-vs-PhaZ 序列未被当作 PhaZ 阳性参考。T141 使用 dated deploy：`deploy/20260911_phaded_competition_phylo_review_01/`；服务器线程为 40；MAFFT `v7.525`，IQ-TREE `3.1.2`，蛋白模型 `LG+G4`，UFBoot `1000`。trimAl 已执行，但其版本命令没有产生可审计输出，manifest 保留为 `pending`。

## 竞争结论

| 结论 | 数量 | 含义 |
|---|---:|---|
| `retain_extracellular_dPHASCL2_like_candidate` | 3 | 最近 dPHASCL2 参考的 MRCA 支持度 ≥90；保留为候选，不是活性结论 |
| `unresolved_extracellular_placement` | 1 | `GCA_025458915.1|JALORL010000164.1_1` 最近参考虽为 DED_hfam_70，但 MRCA 支持缺失，暂不定型 |
| `PhaC_vs_PhaZ_unresolved` | 12 | 12 条均保留为 PhaC-like challenge；树与 nPHAMCL 参考的接近不能证明 PhaZ 解聚酶功能 |
| `cross_family_architecture_unresolved` | 8 | 结构/assignment 异常候选与最近参考家族不一致，系统发育不能消除架构冲突 |
| `retain_dPHASCL1_like_candidate` | 16 | 最近 dPHASCL1 参考距离 ≤1.5 且 MRCA 支持度 ≥90；仍是候选同源证据 |
| `ePhaZ_competition_unresolved` | 2 | `GCA_003021705.1|PXRJ01000007.1_55` 和 `GCA_963829385.1|CAWHHV010000004.1_298` 支持度/距离不足，暂不定型 |

完整的逐候选最近参考、距离、MRCA 支持度和决定见：
`runs/20260911_phaded_competition_phylo_review_01/results/competition_phylo_review_decisions.tsv`。

## 解释边界与后续

- `structure_anomaly` 和 `ephaz_competition` 的 IQ-TREE 日志出现 UFBoot 收敛提示；树已完成，但相关支持值应作为谨慎的候选证据，不应当作确定分类。
- 18 条 `ePhaZ_competition_review` 的 InterPro 字段在当前权威来源中为 `pending`，本轮没有伪造二次注释。
- PhaDED、Pfam、InterPro、motif、SignalP、结构和系统发育均只能支持候选同源/架构潜力，不能证明 PHB/PHA 降解表型。
- 下一步优先对 3 条 dPHASCL2-like、16 条 dPHASCL1-like 以及 1 条 unresolved extracellular 候选做结构域边界复核和结构预测；12 条 PhaC-vs-PhaZ 必须保持 challenge 集，不能进入 ePhaZ curated positive 训练集。
