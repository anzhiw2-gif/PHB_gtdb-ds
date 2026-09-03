# T141 ePhaZ ambiguous 570 条结构复核

日期：2026-08-29
run：`20260829_ephaz_ambiguous_structure_review_02`

## 范围

本轮对分层抽样得到的 570 条 `ambiguous` 候选逐条复核。复核输入为固定的样本 FASTA/TSV、ePhaZ curated-core 与 iPhaZ HMMER `domtblout`、SignalP 结果和邻域上下文。HMM 覆盖度使用 domtblout 的 HMM 坐标（0-based 字段 15-16），同一蛋白的非重叠区段先求并集再除以模型长度；旧抽样脚本的 envelope 坐标和“只取最长区段”问题已修正，历史 run 保留不覆盖。

## 结果

| 指标 | 数量 |
|---|---:|
| 样本总数 | 570 |
| 序列完整（以 M 开始、无内部终止） | 536 |
| 可能 N 端截断（不以 M 开始） | 34 |
| `iPhaZ_consistent` | 511 |
| `partial_ePhaZ_signal` | 47 |
| `mixed_cross_family` | 1 |
| `insufficient_structural_support` | 11 |
| `provisional_iPhaZ_challenge` | 489 |
| `pending_manual` | 81 |

全部样本 SignalP 预测为 `OTHER`。邻域上下文精确匹配 226 条并含 marker，344 条无记录；`no_record` 不解释为“无邻域”。

## 解释边界

`provisional_iPhaZ_challenge` 仅表示序列完整性通过且 iPhaZ HMM 中央区段证据强于 ePhaZ 双区段架构证据。`partial_ePhaZ_signal`、`mixed_cross_family` 和 `pending_manual` 不得进入 ePhaZ 正对照分母，也不能直接用于 PHB 降解表型结论。任何 HMM、domain、SignalP 或邻域证据均只代表候选同源/功能潜力。

## 文件与 SHA-256

- 逐条复核表：`runs/20260829_ephaz_ambiguous_structure_review_02/results/structural_review/ambiguous_structural_review.tsv`
  - `b34695c72528cc9faf73ab17f04071da242964342aaaf0335ce9c18f3321b5b6`
- 复核元数据：`runs/20260829_ephaz_ambiguous_structure_review_02/results/structural_review/structural_review_metadata.json`
  - `3f6f6fc594de94e9188caa9d4e483ffb6512a422b19ac44ac86cd0125d7aff1e`
- 样本 FASTA 输入 SHA-256：`569728ce38929f696e0852ba5849f6e6554507187b79904baf0c49afadd34d3d`
- 样本 TSV 输入 SHA-256：`b610c63022c78026b43a4b3c01b32225f4fc50814499f4b9e8eadebf7b95bc34`

## 后续

优先人工查看 81 条 `pending_manual`，然后查看 47 条 `partial_ePhaZ_signal` 与唯一的 `mixed_cross_family`。建议补充 InterPro/Pfam、跨膜/信号肽和基因组邻域的人工证据字段；在这些证据完成前，不建立 `ePhaZ_architecture_remote`，不调整 broad 正对照分母，不启动 268G 全量重跑。
