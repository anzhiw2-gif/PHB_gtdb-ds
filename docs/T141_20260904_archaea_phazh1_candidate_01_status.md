# T141 20260904 古菌 PhaZh1 candidate-only 重建状态

## 结论

方案 B 已完成候选阶段：固定 8 条古菌 patatin-like 参考序列，构建一次 `archaea_PhaZh1_like_candidate.hmm`，并在固定支持面板/冲突负面板上校准。模型只能作为召回探针，不能独立证明 PhaZh1 或 PHB 解聚活性。

## 输入与来源

- 核心实验背景：`I3RBH0`，*Haloferax mediterranei* PhaZh1，Liu 2015，DOI `10.1128/AEM.04269-14`。
- 关键古菌条目：`M1XPT2`；其余为跨古菌类群的 UniProt patatin-like 支持同源，详见 `runs/20260904_archaea_phazh1_candidate_01/inputs/reference_manifest.tsv`。
- 负面板：4 条 UniProt patatin/phospholipase 注释序列，明确标记为广谱/冲突对照，不代表实验阴性。
- 原始 UniProt JSON 保存在 `runs/20260904_archaea_phazh1_candidate_01/inputs/uniprot_raw/`，并已写入 `input_contract.json`。

## 构建结果

- alignment：`runs/20260904_archaea_phazh1_candidate_01/results/archaea_PhaZh1_like_candidate.sto`
- HMM：`runs/20260904_archaea_phazh1_candidate_01/results/archaea_PhaZh1_like_candidate.hmm`
- HMM SHA-256：`f1794d6dc4a57cde82e0ac793e93bf36457389d8c9e39ed091900a6b39581c1a`
- alignment SHA-256：`a363c98cb5f3c9040f8c77eae48c36a1e64915dc6903b1c0fe956e786c22c98f`
- 远端环境：MAFFT/HMMER 位于 `${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb/bin/`；`hmmbuild` 报告版本 `7.525 (2024/Mar/13)`。

## 校准观察

| 面板 | 条目数 | 命中数 | 解释 |
|---|---:|---:|---|
| 支持面板 | 8 | 8 | 召回良好；仅部分有直接实验背景 |
| 广谱/冲突对照 | 4 | 4 | HMM 独立特异性不足 |
| 负面板中 `E<=1e-5` | 4 | 4 | 不能建立零假阳性阈值 |

因此不把该 HMM 注册到 `formal_scan_models.tsv`，也不启动 GTDB 正式扫描。下一步只对已有古菌候选执行联合审计：`coverage + 完整性 + Archaea taxonomy + BdhA/PhaJ 邻域 + PhaC/phasin 冲突`，输出 `PhaZh1_like_high`、`PhaZh1_like_review`、`archaeal_patatin_exploratory` 和 `non_PhaZh1_patatin` 四层。

联合审计实际结果：`PhaZh1_like_high=197`、`PhaZh1_like_review=7`、`archaeal_patatin_exploratory=2027`、`non_PhaZh1_patatin=76`。该分层基于同一 run 的 probe `domtblout` 和邻域上下文；当前尚未将独立 GTDB taxonomy TSV 的细粒度字段逐条回写，故不能用于正式生物学扩张统计。

## 运行边界

- run：`runs/20260904_archaea_phazh1_candidate_01/`
- deploy：`deploy/20260904_archaea_phazh1_candidate_01/`
- `formal_scan_models.tsv` 未修改。
- 未覆盖历史 `results/`、历史 run 或服务器旧脚本目录。
- 未启动 GTDB 正式扫描；所有结论保持 candidate-only。
