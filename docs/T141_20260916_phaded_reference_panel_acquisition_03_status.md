# T141 PhaDED reference panel acquisition 03

日期：2026-09-16  
Run：`20260916_phaded_reference_panel_acquisition_03`

## 获取结果

本轮从 PubMed/文献页面定位 PHB/PHA depolymerase 的实验研究，并从 NCBI Protein 获取 GenPept 快照。共形成 18 条外部面板候选，18 条均有完整序列 SHA-256；其中 11 条标记为文献实验阳性候选，1 条 family-resolved-negative 候选，6 条 challenge control 候选。

与既有 ledger 做 accession 约束绑定后，4 条记录精确绑定到已有 family，14 条保持 `unresolved_external_candidate`。4 条中 3 条（`AAA87070.1`、`BAF35850.1`、`AAL30107.1`）是训练 accession，不能作为 held-out；`AAF86381.1` 属于已训练 family，也不贡献 37 个 reference-only profile 的校准。因而本轮对 reference-only profile 的新独立留出阳性仍为 0。

## 当前门槛状态

- 新独立留出阳性：0
- 新 family-resolved formal negative：0
- challenge 候选：6（均未绑定到 reference-only family）
- reference-only profile 可校准：0
- 新 family 判定：`false`

未绑定的实验蛋白并未丢弃，而是保留为后续 profile/domain/结构证据的外部候选；未绑定的 negative/challenge 不能用于 family-specific 校准。annotation-only 记录不转写为 formal negative。

## 数据与来源

PubMed 检索/摘要涉及 PMID：`9406404`、`9872779`、`15489436`、`16199568`、`21948827`、`23951224`。原始检索 JSON、NCBI GenPept 快照、候选 FASTA、绑定表和哈希均保存在本 run 的 `inputs/` 与 `results/` 下。

服务器 `<SERVER_HOST>` 预检为 80 logical CPU、2× RTX 4090、GPU 利用率 0%；本轮主要是联网检索和本地 provenance 整理，未启动服务器计算，使用线程数 0，单任务上限 40 仍有效。

所有结论仅用于 sequence-level candidate evidence，不等同于 GTDB 菌株或蛋白已经获得 PHB/PHA 降解实验验证。
