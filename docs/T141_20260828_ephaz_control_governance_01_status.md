# T141 ePhaZ 控制治理与竞争审计状态

日期：2026-08-28
run：`20260828_ephaz_control_governance_01`
服务器：`<SERVER_USER>@<SERVER_HOST>`

## 本轮决策

- 六条复核序列统一标记为 `iPhaZ_like_challenge`，其注释理由保留
  `cross_family_confounded` 和 `annotation-only`。
- 六条序列从 ePhaZ broad 正对照分母移除，仅保留在
  `inputs/control_panels/iPhaZ_like_challenge.faa` 挑战集。
- 没有序列被升级为 `ePhaZ_curated_core` 或 ePhaZ bridge；当前
  `ePhaZ_architecture_remote` 数量为 0。
- 未启动 268G 全量筛选。

## 小规模分组校准

输入为 4 条实验核心、6 条 iPhaZ-like challenge 和 18 条负对照。ePhaZ
模型使用 `ePhaZ_curated_core.hmm`（SHA-256
`96668851b42fba67d5fa8903e4ac25527769abad8b8e8d7366b3a62db1a9043f`），
iPhaZ 模型使用 `iPhaZ.v2.hmm`（SHA-256
`4a1d15d37e8f55bea89fbc6bc33b1e1c45b166c80aad3976c73c087937f7d311`）。

| 模型 | 对照类别 | 测试数 | 命中数 | 解释 |
|---|---:|---:|---:|---|
| ePhaZ | `ePhaZ_curated_core` | 4 | 4 | 核心正对照全部命中 |
| ePhaZ | `iPhaZ_like_challenge` | 6 | 0 | 未被 ePhaZ 核心模型混入 |
| ePhaZ | `negative` | 18 | 0 | 本轮阈值下无假阳性 |
| iPhaZ | `ePhaZ_curated_core` | 4 | 0 | 无跨家族命中 |
| iPhaZ | `iPhaZ_like_challenge` | 6 | 6 | challenge 的家族归属信号一致 |
| iPhaZ | `negative` | 18 | 0 | 本轮阈值下无假阳性 |

完整表格和命中明细位于：
`runs/20260828_ephaz_control_governance_01/results/calibration_control_classes/`。

## 无信号肽竞争审计

对现有 ePhaZ tier1 FASTA 中 SignalP=OTHER 的 16,836 条蛋白，以 20 bits
差值作为竞争判定边界：

| 分配 | 数量 | 比例 |
|---|---:|---:|
| `ePhaZ_like` | 3,744 | 22.24% |
| `iPhaZ_like` | 4,865 | 28.90% |
| `ambiguous` | 6,991 | 41.52% |
| `no_reportable_hit` | 1,236 | 7.34% |

结果位于：
`runs/20260828_ephaz_control_governance_01/results/no_signal_competition/`。

`ambiguous` 占 41.52%，说明 ePhaZ 与 iPhaZ 模型在无信号肽子集中存在明显
竞争重叠。该结果只能说明模型同源/功能潜力分类，不能直接证明 PHB 降解
表型，也不足以授权全量重跑。

## 溯源与环境边界

- 本地合同：`runs/20260828_ephaz_control_governance_01/input_contract.json`。
- T141 原始合同副本：`runs/20260828_ephaz_control_governance_01/input_contract.t141.json`。
- 服务器只执行了 dated deploy：
  `deploy/20260828_ephaz_control_governance_01/scripts/`。
- 服务器合同中的 GTDB taxonomy、metadata、tree 保持 `pending`，没有伪造哈希。
- tier1 FASTA SHA-256：
  `4b09c57c42b2f59106c77cc214adbdcde2e08a442f82ab36237fdb8c4f8f3574`。
- SignalP 预测 SHA-256：
  `0c46b86ed942f92c23314d809cfde68d490b1aad16cbe1d3818a747541ca8b22`。

## 下一步建议

1. 对 6,991 条 `ambiguous` 做结构域覆盖、长度和邻域分层抽样复核。
2. 在确认可解释的 ePhaZ 架构证据前，不调整 broad 正对照分母，不执行
   268G 全量筛选。
3. 将本轮结果作为治理证据单独提交；不要覆盖历史 `results/` 或旧 run。
