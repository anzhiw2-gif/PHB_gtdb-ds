# 家族特异性 HMM 筛选标准与下一步

日期：2026-09-03
适用范围：`PHB_gtdb-ds`、run 13 及其下游候选报告

## 标准

不同家族不强制使用相同的 coverage 阈值。统一要求是：阈值、校准依据、用途和证据边界必须逐家族记录。HMM 命中表示同源/功能潜力候选，不等于实验确认的 PHB 或 MCL-PHA 降解表型。

| 家族 | 正式筛选/计数 | 候选分层或用途 | 解释边界 |
|---|---|---|---|
| `OH` | `E<=1e-5` 且 `min_cov=0.6` | tier1 中 `coverage>=0.8` 可作为高可信优先复核层；`0.6-0.8` 保留为 review | `0.6` 有尼龙水解酶假阳性对照依据；不能据此宣称实验活性 |
| `ArchPhaZ_hydrolase` | 维持当前 registry 规则；不把全部 tier1 作为高可信古菌解聚酶 | `coverage>=0.8`：high-confidence；`0.6-0.8`：review；`<0.6`：exploratory | 低 coverage 可能是片段化/远缘同源；必须避免“古菌扩张”表述 |
| `ePhaZ_curated_core` / `ePhaZ_broad_discovery` | `E<=1e-5`，暂不设通用 coverage 硬门槛 | curated 用于严格层；broad 保持 discovery 层 | 远缘 MCL-PHA 候选仍需独立验证，不直接注册生产模型 |
| `PhaJ` | 不计入直接解聚酶阳性 | PHB 动员通路背景 | 不是直接解聚酶 |
| `BdhA` | 不计入直接解聚酶阳性 | 背景代谢 | 普遍存在，特异性不足 |
| `PhaC` | 不计入直接解聚酶阳性 | PHA 合成/颗粒背景或共现特征 | 不能单独证明降解 |
| `phasin` | 不设阳性/阴性硬门槛 | 可选辅助特征和模型灵敏度警示 | 当前命中数过低，不能解释为生物学缺失 |

## 当前推荐的下一步

1. 冻结上述解释口径，不再为统一阈值而重建所有 HMM。
2. 输出 OH 与 `ArchPhaZ_hydrolase` 的分层候选表，保留 accession、coverage、E-value、序列完整性和分类信息。
3. 对 OH 新增门类和 `0.6-0.8` 边界命中做小规模跨门抽查。
4. 对 `ArchPhaZ_hydrolase` 的 497 个 `coverage>=0.8` 基因组优先核查，尤其区分古菌与细菌来源。
5. 将 `PhaJ`、`BdhA`、`PhaC`、`phasin` 仅作为背景/辅助字段，不并入直接解聚酶阳性计数。
6. 只有出现新的、可复现的失败模式并获得独立实验/外部验证集时，才启动新的模型校准或正式扫描。

## 结果状态

- 正式 registry：保持不变。
- run 13 raw HMMER 输出：保持不变。
- `runs/20260902_oh_arch_aux_audit_01/`：作为 candidate-only 审计记录。
- 当前结果：可进入分层汇总、有限抽查和论文/报告准备；不能直接作为实验表型结论。
