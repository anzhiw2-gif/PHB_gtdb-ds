# T141 ePhaZ 独立实验阳性与结构近邻负对照面板

运行目录：`runs/20260830_ephaz_independent_positive_negative_01/`

本轮目标是建立一个不含 bridge 训练蛋白（`Q51871`、`Q5SLU4`）及当前 ePhaZ 核心种子（`B2NHN2`、`O05527`、`P12625`、`Q51718`）的外部验证面板。冻结 HMM、阈值和 GTDB 全量扫描均未改变或启动。

## 纳入标准

- 独立 PHB 阳性：原始论文明确报告具体蛋白/基因的克隆或测序，并有纯化、重组表达、清除圈或 PHB 降解活性证据；序列完整且 accession 稳定。
- MCL-PHA 阳性：有具体蛋白的实验支持，但底物主要为中链 PHA，单独统计，不能充当 PHB 阳性。
- 胞内负对照：实验或明确注释为胞内 PHB 解聚酶，缺少胞外分泌架构。
- MCL-PHA 非 PHB 负对照：有胞外结构和实验酶活，但论文明确不水解 PHB。
- 注释近邻负对照：具有 AB hydrolase/PhaZ/PHB-depolymerase 注释或相近结构，但缺少直接蛋白实验支持。
- 片段负对照：明确为 fragment 或不完整序列；保留用于挑战，不进入阳性分母。

所有原始 UniProt/NCBI/PubMed 响应均保存在 `inputs/`，输出记录及响应哈希见 `results/raw_response_index.tsv` 和 `results/ephaz_panel_evidence.tsv`。

## 当前面板

| 面板 | 条数 | 代表 accession | 解释 |
|---|---:|---|---|
| `independent_experimental_positive` | 5 | `AAB40611.1`, `O24719`, `A0A8W8`, `Q9LBN6`, `Q5YEW3` | 直接实验支持的 PHB 胞外解聚酶 |
| `mcl_pha_experimental_positive` | 1 | `Q6UFW4` | MCL-PHA 胞外解聚酶；底物范围控制 |
| `intracellular_non_ephaz_negative` | 3 | `O87189`, `Q7WT48`, `Q7WT49` | 胞内 PHB 解聚酶/非分泌架构 |
| `mcl_pha_non_phb_negative` | 1 | `Q84C08` | 实验明确不水解 PHB 的结构近邻 |
| `annotation_only_near_neighbor_negative` | 5 | `P26495`, `Q5Q138`, `Q7X5S3`, `A0A375HYL0`, `A0A1C9W3H4` | 仅注释或序列线索，不当作实验阴性 |
| `fragment_or_incomplete_negative` | 2 | `Q9AGB6`, `J7K890` | 片段/不完整，保留作挑战 |
| 近邻负对照合计 | 11 |  | 仅用于 ePhaZ 结构近邻压力测试 |

## 排除但保留挑战

- `Q71KW6`：与现有 `iPhaZ` challenge 序列完全相同，不能重复进入新负对照分母。
- `Q939Q9`：Paucimonas `PhaZ7` 有直接 PHB 实验证据，但 UniProt 标记为 fragment，且原论文描述为独立 family 9；不纳入完整阳性分母。

两条序列及排除理由保存在 `inputs/excluded_candidates.tsv`，序列保存在 `results/ephaz_excluded_challenge.faa`，不参与本轮校准分母。

## 产物

- `results/independent_experimental_positive.faa`
- `results/mcl_pha_experimental_positive.faa`
- `results/intracellular_non_ephaz_negative.faa`
- `results/mcl_pha_non_phb_negative.faa`
- `results/annotation_only_near_neighbor_negative.faa`
- `results/fragment_or_incomplete_negative.faa`
- `results/ephaz_excluded_challenge.faa`
- `results/ephaz_panel_evidence.tsv`
- `results/raw_response_index.tsv`
- `input_contract.json`

## 证据边界

这些面板用于独立验证、底物范围和结构近邻压力测试。HMM 命中仍只能解释为候选同源或功能潜力；不能据此宣称 GTDB 候选已验证 PHB 降解表型。下一步应在不改冻结模型的前提下做小规模外部校准，并分别报告 PHB、MCL-PHA、胞内和注释型分母。
