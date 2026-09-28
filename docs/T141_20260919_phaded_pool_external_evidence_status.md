# T141 — 池外 strong-tier 全证据筛选：状态文档

**Run：** `runs/20260919_phaded_pool_external_evidence_01`
**日期：** 2026-09-20
**状态：** `completed_candidate_only`
**目标：** 对池外 strong-tier（E<1e-30，30,632 条）施加与池内相同的多证据筛选。

---

## 0. 一句话结论

池外 strong-tier **30,632** 条经与池内相同的多证据筛选（SignalP 定位 + InterPro 域 + motif/催化域架构 + 唯一超家族归属）后，**高可信度 4,131 条（13.5%）**，其余 26,501 条被显式 hold（带 `hold_reason`）。

## 1. 结果

| 超家族 | 高可信度 | 占该族 % |
|---|---:|---:|
| extracellular dPHASCL type 2 | **2,238** | 29.3% |
| intracellular nPHASCL without lipase box（Cys） | **1,492** | 72.8% |
| extracellular dPHAMCL | 184 | 72.4% |
| extracellular dPHASCL type 1 | 134 | 0.68% |
| extracellular native-SCL/PhaZ7-like | 82 | 27.2% |
| periplasmic | 1 | 50% |
| intracellular nPHAMCL | 0 | 0% |

亲核体类型：ser 2,557 / cys 1,492 / ahsmg 82。hold 分解：architecture_criteria 15,509 / localization_conflict 7,650 / nucleophile_conflict 3,120 / common_criteria 222。

## 2. 关键发现

- **池外 type1 通过率极低（0.68%）**：19,765 条 type1 候选里仅 134 条通过——因为池外是「仅 profile 得分」层，绝大多数缺 `type1_verified` 几何（氧阴离子孔 + lipase box）。这与池内 type1（7.2% 通过）形成对比，说明发现层 HMM 的 type1 归属远比 trained 层宽。
- **Cys 型最稳**：2,049 条 Cys 候选中 1,492 条通过（72.8%），因 PF06850（检出 1,518 条）+ 无 GxSxG + 无分泌信号的判据相对特异。
- 结构断点仍成立：池外 `candidate_ded_alignment_members = 0`（与池内相同），His/Asp 等不可靠比对传递。

## 3. 证据层来源

- SignalP 6（`--torch 30 --write 6`）→ 30,632 条 `signalp_class`（SP 13,936 / OTHER 10,284 / LIPO 6,281 / TAT 75 / TATLIPO 21 / PILIN 35）。
- InterProScan 5.78（`--cpu 60`，重启）→ 427,479 行命中 / 30,410 序列；PF06850 检出 **1,518** 条、PF10503 7,299 条。
- motif 补全（Task 7 脚本，参考层输入冻结复用）→ lipase_box supported 23,583 / oxyanion supported 4,041 / PF06850 detected 1,518。
- 催化域 type（几何判定 helper）→ type1 580 / type2 3,001 / undetermined 27,051。

## 4. 治理与事故

- test-first 未新增测试（复用已验证脚本 + 3 个自检过的 helper）；过滤器本身含 nucleophile/hold 逻辑（此前已测）。
- **事故 A**：服务器 `strong_pool_external.noasterisk.faa` 曾被内联 ssh 引号 bug 截断为 2 字节 → 用完好原始 FASTA（`fc7a0c3d…`）本地重建（`dff1f83b…`）并重传。
- **事故 B**：首次 InterPro `--cpu 24` 跑 ~2h 后按操作者要求停掉、改 `--cpu 60` 重启（损失 ~2h，最终 06:13 完成）。
- `input_contract.json` 已写（gtdb 三槽 pending）；未覆盖/删除任何历史产物。

## 5. 边界

全部 candidate-only；「高可信度」= 候选证据最强筛选层，不等于已验证 PHB/PHA 降解表型。
