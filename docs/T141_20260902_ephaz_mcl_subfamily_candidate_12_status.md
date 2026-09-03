# T141 ePhaZ MCL 子家族候选优化状态

日期：2026-09-02
运行：`20260902_ephaz_mcl_subfamily_candidate_12`
服务器：`<SERVER_USER>@<SERVER_HOST>`
软件：`MAFFT 7.525`、`HMMER 3.4`，`--thread/--cpu 60`

## 结果

- `mcl_classical`：8/8 个经典核心阳性恢复。
- `mcl_lipase_associated`：恢复 `WHU94860.1`，但该 profile 只有 1 条实验阳性序列，属于探索性 singleton。
- `mcl_streptomyces`：恢复 `AZSS01000334.1:12616-13485(-)`，但该 profile 只有 1 条实验阳性序列，属于探索性 singleton。
- 3 条 negative：0 个 accepted hit。
- challenge：0 个 accepted hit。
- `J7K890` 对 Streptomyces profile 的 E-value 为 `1.2e-115`，但 HMM 覆盖度仅 `0.637`，按 `E<=1e-5` 且 coverage `>=0.90` 被拒绝。

## 边界

这是 candidate-only 校准，不是独立验证，也不是 PHB/MCL-PHA 表型证明。两个 divergent profile 仍缺少跨属、独立、实验确认的额外阳性，不能直接用于正式 GTDB 扫描。`pipeline/config/formal_scan_models.tsv` 未修改，正式扫描未启动，release 仍 blocked。生产分类器仍保持原行为；本 run 的 coverage 规则只用于候选评估，后续需先完成代码级 coverage 接口和独立验证，再申请集成。

证据文件：

- `runs/20260902_ephaz_mcl_subfamily_candidate_12/results/SUBFAMILY_CANDIDATE_DECISION.md`
- `runs/20260902_ephaz_mcl_subfamily_candidate_12/results/subfamily_calibration.tsv`
- `runs/20260902_ephaz_mcl_subfamily_candidate_12/results/sha256sums.txt`
