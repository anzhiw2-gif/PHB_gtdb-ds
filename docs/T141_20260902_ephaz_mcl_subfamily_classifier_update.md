# ePhaZ MCL-PHA 子家族分类器更新

## 范围

本次只修改 candidate-only 软件逻辑，不启动 run 13，不修改
`pipeline/config/formal_scan_models.tsv`，也不改写既有运行目录或服务器历史
deploy。HMM 命中仍只表示同源/功能潜力，不等同于实验确认的 PHB 降解表型。

## 变更

- `pipeline/scripts/classify_ephaz_dual_profiles.py` 保留旧的 `--mcl-hmm` 单 profile
  入口，并新增可重复的 `--mcl-profile NAME=PATH` 入口。
- 多 profile 输出增加 `best_mcl_subfamily`，以及每个 profile 的独立 E-value 和
  bit score 列；显式名称 `classical`、`lipase_associated`、`Streptomyces` 分别可
  输出为 `MCL_classical`、`MCL_lipase_associated`、`MCL_Streptomyces`。
- `sp|...|`、`tr|...|`、`ref|...|`、`gi|...|` FASTA/HMMER target 统一解析为真实
  accession；带坐标的完整候选 ID 原样保留。
- `pipeline/scripts/ephaz_bridge_loo_calibration.py` 使用同一 accession 规则，避免
  LOO 命中表与 FASTA 控制集因 `sp`/`tr` 前缀而错配。

## 证据边界

当前 6 条 MCL 候选中，`Q51718`、`Q6UFW4`、`Q6MH49` 可组成 classical 分支的
candidate seed；`WHU94860.1`（lipase-associated）和
`AZSS01000334.1:12616-13485(-)`（Streptomyces）各只有 1 条实验确认完整序列，
暂不构建这两个正式 HMM。需继续补充跨属、实验确认且序列完整的阳性后，才可进行
子家族 HMM 的 candidate-only 校准，并重新评估是否申请正式扫描授权。

## 验证

已通过 `python -m unittest discover -s pipeline/tests -v`（123 tests，1 个平台权限
相关 skip）、`python -m compileall -q pipeline/scripts pipeline/tests` 和
`git diff --check`。
