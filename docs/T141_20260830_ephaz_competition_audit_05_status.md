# T141 ePhaZ/iPhaZ 无信号肽竞争审计状态

运行目录：`runs/20260830_ephaz_competition_audit_05/`

对象为 ePhaZ tier1 中 SignalP `OTHER` 的 16,836 条序列。以 bitscore 差值 20 为分类边界：

| 分类 | 数量 |
|---|---:|
| `ePhaZ_like` | 16,809 |
| `iPhaZ_like` | 9 |
| `ambiguous` | 18 |

逐 accession 表 `results/competition/no_signal_competition.tsv` 同时记录 ePhaZ/iPhaZ 的 E-value、bitscore、HMM 起止坐标、目标序列起止坐标、HMM coverage 和 target coverage。

这些分类只反映模型竞争关系，不是功能或表型注释。`iPhaZ_like` 与 `ambiguous` 序列在后续候选解释中必须单列，不能仅凭 ePhaZ 命中归入 PHB 胞外解聚酶。

本轮仅审计 ePhaZ tier1 子集，没有启动 268G 全量扫描。
