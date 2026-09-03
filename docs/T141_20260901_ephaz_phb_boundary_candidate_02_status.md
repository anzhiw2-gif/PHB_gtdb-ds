# T141 ePhaZ PHB boundary candidate 02

本次仅执行 run-scoped 小规模校准，不修改正式模型注册表、run 13 或全库扫描。

- Run：`runs/20260901_ephaz_phb_boundary_candidate_02`
- Deploy：`deploy/20260901_ephaz_phb_boundary_candidate_02`
- 新增边界：`Q4W8C9`（经审阅的胞内 3HB-oligomer hydrolase，仅人工无定形 PHB 活性，天然/半结晶 PHB 阴性，PMID `16030206`）；`Q5Y152` 仅为未审阅的注释型 MCL-PHA-like 探索对照，不宣称有 accession 级实验验证。
- 两条序列均为完整 UniProt 记录，仅加入 `intracellular_non_ephaz_negative`，不加入 seed。Q5Y152 与 *P. putida* KT2442 论文 accession 不同，不能直接绑定 DOI `10.1074/jbc.M608119200`。

在 `E<=1e-5`、coverage=0 时，baseline、`PHB-focused_no_Q51718` 和
`PHB-focused_plus_independent` 的胞内/状态边界均为 `0/5` 命中。PHB
实验阳性仍为 `5/5`；去除 Q51718 的候选仍将 MCL-PHA 非 PHB 对照降为
`0/1`，但 MCL-PHA 阳性为 `0/1`，因此不能称为 PHB 特异性分类器。

服务器使用 HMMER 3.4、MAFFT 7.525、8 CPU；完整 SHA-256 见
`runs/20260901_ephaz_phb_boundary_candidate_02/results/output_sha256.tsv`
和 `input_contract.json`。所有命中只表示序列同源/功能潜力。
