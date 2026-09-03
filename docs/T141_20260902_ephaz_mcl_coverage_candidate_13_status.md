# T141 ePhaZ MCL coverage-aware candidate status

日期：2026-09-02
运行：`20260902_ephaz_mcl_coverage_candidate_13`

本轮未重新构建 HMM，也未启动服务器或正式 GTDB 扫描；复用了
`20260902_ephaz_mcl_subfamily_candidate_12` 已核对的 `domtblout`，测试分类器新增的显式
`min_domain_coverage=0.9` 参数。

结果：

- `WHU94860.1`：`MCL_lipase_associated`，coverage `1.0`。
- `AZSS01000334.1:12616-13485(-)`：`MCL_streptomyces`，coverage `1.0`。
- `J7K890`：`no_hit`；E-value `1.3e-115` 但 coverage `0.637`。
- `O87189`、`Q7WT48`、`Q7WT49`：均 `no_hit`。

代码变化集中在 `pipeline/scripts/classify_ephaz_dual_profiles.py`：增加 domtblout 区间合并覆盖度解析，coverage 过滤保持默认关闭，只有显式传参才启用；同时修复 `mcl_` profile 名称的重复前缀。现有单 profile 接口行为保持兼容。

本轮仍是 candidate-only 软件校准，不提供 PHB/MCL-PHA 表型证明。两个 divergent profile 仍为 singleton，正式扫描与模型注册仍 blocked。
