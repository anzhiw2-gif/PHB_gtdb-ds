# T141 ePhaZ 双层 run 状态（2026-08-28）

运行目录：`runs/20260828_ephaz_split_01/`

本 run 从 `20260828_seedclean_07` 的 ePhaZ 输入建立两个互斥层：

- `ePhaZ_curated_core`：4 条；均有明确实验支持、典型架构、完整长度且不少于 200 aa。
- `ePhaZ_broad_discovery`：4,454 条；包含远缘/注释型序列及全部 61 条 `<200 aa` 序列。

61 条短序列逐条写入 `inputs/ephaz_layers/ePhaZ_short_sequence_review.tsv`，未因长度删除；短序列均保留在 broad 层。

## HMM 与校准

| family | HMM SHA-256 | alignment SHA-256 |
|---|---|---|
| `ePhaZ_curated_core` | `96668851b42fba67d5fa8903e4ac25527769abad8b8e8d7366b3a62db1a9043f` | `467966fc5a98b84f5e8ad37ee7bb836ab5eff0d980e149b88122b6afebb77f98` |
| `ePhaZ_broad_discovery` | `3f981c5a12b7a7ec9389eab8e329f03e220875bdbce94f95849821e32dc127bc` | `13788bf99b3a07a8045718ebe884dd1f1e5e6c7e78196556559ec225e72951bc` |

校准结果位于 `runs/20260828_ephaz_split_01/results/calibration/`：

- core：4/4 正对照命中，18 个负对照均未命中（所有测试阈值至 `1e-20` 均保持）。
- broad：10 条正对照中 4 条命中，6 条远缘/注释型对照未命中；18 个负对照均未命中。

broad 层结果是模型覆盖度证据，不是对未命中序列的生物学否定，也不等同于已验证 PHB 降解表型。当前不启动 268G 全量筛选；下一步应先对 6 条 broad 正对照进行结构域/架构复核，并决定是否扩充 broad seed 后再校准。
