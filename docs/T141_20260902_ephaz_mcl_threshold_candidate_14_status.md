# ePhaZ MCL threshold sensitivity candidate status

日期：2026-09-02
运行：`20260902_ephaz_mcl_threshold_candidate_14`

本轮固定 `20260902_ephaz_mcl_subfamily_candidate_12` 的三个 profile 和 HMMER 输出，测试 4 个 E-value 阈值与 5 个 HMM coverage 阈值，共 20 组组合；未重建 HMM、未连接服务器、未启动正式扫描。

在保留的 `E<=1e-5` 下：

- coverage `0.0` 或 `0.6`：10/10 核心阳性，但 `J7K890` 产生 1 个 challenge 命中。
- coverage `0.8`、`0.9`、`0.95`：均为 10/10 核心阳性、0 challenge、0 negative。

因此建议下一轮候选参数采用 `E<=1e-5`、`min_domain_coverage=0.8`。这是当前小面板上的最小有效完整性门槛，不代表独立验证阈值。两个 divergent profile 仍为 singleton exploratory models，正式 GTDB 扫描和注册表更新继续 blocked。
