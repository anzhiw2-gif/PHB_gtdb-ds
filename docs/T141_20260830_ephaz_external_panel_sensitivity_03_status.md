# T141 ePhaZ 外部面板敏感性校准状态

运行目录：`runs/20260830_ephaz_external_panel_sensitivity_03/`

本轮使用独立外部面板，对 `ePhaZ.v2`、`ePhaZ_curated_core` 和 `iPhaZ.v2` 做 E-value `1e-5/1e-10/1e-20` × HMM coverage `0.0/0.4/0.6/0.8` 网格校准。PHB、MCL-PHA、胞内、MCL-PHA 非 PHB、annotation-only 和 challenge 始终分开统计。

## 主要结果

- `ePhaZ.v2`：在 `1e-5/0.4` 下 PHB 实验阳性为 `5/5`；coverage `0.6` 时为 `2/5`，coverage `0.8` 时为 `0/5`。
- `ePhaZ_curated_core`：在较高 coverage 下同样损失部分 PHB 实验阳性召回，不适合直接采用 `0.6` 或 `0.8` 作为通用门槛。
- `iPhaZ.v2`：胞内 PHB 对照在宽松 coverage 下为 `3/3`，支持将其作为竞争审计模型，而非 ePhaZ 的替代判定器。
- `Q84C08`（明确不水解 PHB 的 MCL-PHA 解聚酶）被两个 ePhaZ 模型命中，说明 ePhaZ HMM 覆盖相近胞外 PHA 解聚结构空间，不具备单独的 PHB 底物特异性证明力。

完整网格位于 `runs/20260830_ephaz_external_panel_sensitivity_03/results/calibration/panel_sensitivity.tsv`。

## 冻结结论

现有 registry 的 E-value `1e-5`、coverage `0.0` 暂时保持不变；coverage 仅作为后处理分层证据，不升级为 ePhaZ 通用硬过滤。此次校准不授权 GTDB 全量重跑，也不改变任何 HMM 文件。

`ePhaZ_broad_discovery.hmm` 仍为 `pending`：本地、Git 历史、deploy 和 T141 目录均未找到实际文件，历史 SHA-256 不能替代实体文件。
