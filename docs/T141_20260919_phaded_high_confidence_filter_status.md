# T141 — PhaDED 高可信度候选过滤（多证据筛选）：状态文档

**Run：** `runs/20260919_phaded_high_confidence_filter_01`
**日期：** 2026-09-19
**状态：** `completed_candidate_only`
**判据：** `inputs/filter_criteria.md`（预注册决策表）
**实现：** `pipeline/scripts/filter_phaded_high_confidence.py`

---

## 1. 结果：池内（109,087）高可信度候选 **36,611 条**

| 超家族 | 候选总数 | **高可信度** | 主要淘汰原因 |
|---|---:|---:|---|
| intracellular nPHASCL without lipase box（Cys） | 30,855 | **28,534** | 少部分无 PF06850 或强得分 |
| extracellular dPHASCL type 1 | 67,720 | **4,884** | **48,038 条无分泌信号**（胞外期望不满足）+ 14,433 条无 type1 几何 |
| extracellular dPHASCL type 2 | 2,653 | **2,110** | 少量无分泌信号/几何 |
| intracellular nPHAMCL | 6,950 | **987** | 5,807 条无 lid（胞内 mcl 必需判据） |
| extracellular dPHAMCL | 280 | **96** | 184 条无分泌信号 |
| intracellular nPHASCL with lipase box | 13 | **0** | 仅 13 条归属（trained 路径不覆盖该超家族），全部证据弱 |
| periplasmic | 0 | 0 | 池内无候选 |
| PhaZ7-like | 0 | 0 | 池内 AHSMG 零命中（已知） |
| **合计** | **108,471** | **36,611** | — |

## 2. 筛选后的核心发现

1. **type 1 是最大"水货"来源**：67,720 条 `dPHASCL1_like` 中，**48,038 条（71%）没有分泌信号**——它们顶着"extracellular"的标签却没有胞外定位证据。这是与 Knoll 2009 记录的"注释胞内实为胞外/反之亦然"同一类问题的全库量化。**只有 4,884 条同时满足分泌信号 + type1 几何 + 疏水 x₁。**
2. **Cys 型是最大高可信度池**：28,534 条满足"无 GxSxG + PF06850 命中 + 无分泌信号 + 强 profile"——这批是当前最有把握的胞内 Cys 型候选。
3. **nPHAMCL 的 lid 判据砍掉了 84%**：6,950 条 nPHAMCL 候选中只有 989 条全库有 lid 支持 → 987 条通过。lid 是胞内 mcl 与胞外 mcl 的判别标志（PhaZKT 判据），严格是对的。
4. **with-lipase 高可信度为零**：该超家族（hfam_2）没有 trained profile，发现层 HMM 又极宽（1.24M 池外噪声），因此池内无法构成高可信度集合——这如实反映了该家族的证据现状。

## 3. 池外候选（只有 profile 得分，架构/定位层 pending）

按预注册：池外候选单独产出"profile 强得分层"（`best_evalue < 1e-30`），**不与池内高可信度合并计数**。数量见 `runs/20260919_phaded_trained_profile_recall_status.md` 与 `runs/20260918_phaded_family_discovery_status.md`（强命中合计约 4.6 万条：trained 强命中 ~4.2 万 + 发现层特异家族 ~4 千）。

### 3.1 with-lipase 排除集的 deferred 归档（2026-09-19 决策）

池外 with-lipase 命中 **1,206,655 条**（trained 改判 32,157 条后）被排除在高可信度之外，
但**不删除**：已保留为 deferred 结构验证层，见
`runs/20260919_phaded_with_lipase_deferred_archive_01/`（含池内 13 行 filter 口径 +
5,711 行宽口径；构建脚本 `pipeline/scripts/build_phaded_with_lipase_deferred_tier.py`）。
该层不进入任何候选计数，后续仅经结构证据（PDB 8YNV Foldseek / 结构预测比对）召回。

## 4. 产物

- `results/high_confidence_candidates.tsv`（36,611 行，含全部判据列）
- `results/filter_summary.json`（各超家族通过/失败分解）
- `inputs/filter_criteria.md`（预注册判据表）
- `pipeline/scripts/filter_phaded_high_confidence.py`

## 5. 边界

- 高可信度 ≠ 已验证表型：只是"候选证据最强"的筛选层，仍是 candidate-only。
- 判据缺失一律不通过（fail-closed），未把缺失当阴性。
- 该集合用于后续分析（如实验验证优先级排序、结构预测候选池）。
