# T141 ePhaZ ambiguous 分层抽样复核

日期：2026-08-29
run：`20260829_ephaz_ambiguous_sampling_01`
服务器：`<SERVER_USER>@<SERVER_HOST>`

## 范围与方法

本轮只处理前一轮竞争审计中判定为 `ambiguous` 的 6,991 条蛋白，不改变
竞争阈值、不扩展候选集合，也未启动 268G 全量筛选。对这 6,991 条候选使用
本轮治理所绑定的两个模型重新输出 domtblout：

- ePhaZ：`ePhaZ_curated_core.hmm`，SHA-256
  `96668851b42fba67d5fa8903e4ac25527769abad8b8e8d7366b3a62db1a9043f`；
- iPhaZ：`iPhaZ.v2.hmm`，SHA-256
  `4a1d15d37e8f55bea89fbc6bc33b1e1c45b166c80aad3976c73c087937f7d311`。

domain 覆盖度按 domtblout 的 HMM 坐标 `(hmm_to-hmm_from+1)/qlen` 计算，
取每个模型的最大 domain 覆盖度。邻域只接受候选完整 ID
`genome|locus` 与上下文表 `genome + hit_locus` 的精确匹配；没有记录写为
`no_record`，不解释为“没有邻域”。

## 分层边界与总体分布

四个维度组成互斥 strata：

- bitscore 差值：`0-<2`、`2-<5`、`5-<10`、`10-<15`、`15-<20`；
- 长度：`<250`、`250-399`、`400-599`、`>=600 aa`；
- domain 覆盖：`<0.5`、`0.5-<0.8`、`>=0.8`；
- 邻域：`marker_present`、`no_record`。

总体边际分布：

| 维度 | 分布 |
|---|---|
| bitscore 差值 | 0-<2: 834；2-<5: 1,270；5-<10: 1,848；10-<15: 1,629；15-<20: 1,410 |
| 长度 | <250: 99；250-399: 5,595；400-599: 1,250；>=600: 47 |
| domain 覆盖 | <0.5: 810；0.5-<0.8: 6,089；>=0.8: 92 |
| 邻域 | marker_present: 1,750；no_record: 5,241 |

邻域与覆盖度交叉分布为：`marker_present` 中 `<0.5` 198、`0.5-<0.8`
1,528、`>=0.8` 24；`no_record` 中 `<0.5` 612、`0.5-<0.8` 4,561、
`>=0.8` 68。邻域快照共有 212,240 行，其中 ePhaZ context 82,380 行。

共形成 86 个非空 strata。每个 strata 使用固定种子 `20260829` 按
`SHA256(seed|stratum|accession)` 排序，最多抽取 10 条，得到 570 条样本。
样本 FASTA 与 TSV accession 一一相等，均为 570 条。

## 复核材料

- 全量分层表：`runs/20260829_ephaz_ambiguous_sampling_01/results/stratified_review_v2/ambiguous_stratified.tsv`
- 570 条样本表：`runs/20260829_ephaz_ambiguous_sampling_01/results/stratified_review_v2/ambiguous_sample.tsv`
- 570 条样本 FASTA：`runs/20260829_ephaz_ambiguous_sampling_01/results/stratified_review_v2/ambiguous_sample.faa`
- strata 汇总：`runs/20260829_ephaz_ambiguous_sampling_01/results/stratified_review_v2/strata_summary.tsv`
- 采样 provenance：`runs/20260829_ephaz_ambiguous_sampling_01/results/stratified_review_v2/sampling_metadata.json`
- ePhaZ domtblout：`runs/20260829_ephaz_ambiguous_sampling_01/results/domain_search/ePhaZ_curated_core.dom`
- iPhaZ domtblout：`runs/20260829_ephaz_ambiguous_sampling_01/results/domain_search/iPhaZ.dom`

关键输出 SHA-256：

| 文件 | SHA-256 |
|---|---|
| `ambiguous_sample.faa` | `569728ce38929f696e0852ba5849f6e6554507187b79904baf0c49afadd34d3d` |
| `ambiguous_sample.tsv` | `b610c63022c78026b43a4b3c01b32225f4fc50814499f4b9e8eadebf7b95bc34` |
| `ambiguous_stratified.tsv` | `216b09db1d227dee9499faa12d5436652f14e047979f74e909b713c929e81f42` |
| `sampling_metadata.json` | `87adb8c821990c3f9f0ed37e4f0fc3e23632d231c616fc829bbdfde6401e1775` |

## 解释边界与建议

`ambiguous` 中 6,089 条的最大模型 domain 覆盖落在 `0.5-<0.8`，而 41.52%
的候选此前已显示模型竞争接近。邻域表中 5,241 条为 `no_record`，这是当前
上下文输入未覆盖或未匹配的证据状态，不是邻域阴性。

因此，本轮抽样适合人工复核模型覆盖、蛋白长度、结构域架构和邻域记录是否
支持继续保留；不应把任何 strata 或抽样比例直接转成 PHB 降解表型结论，
也不足以授权 268G 全量重跑。优先人工查看：

1. `0-<5` bitscore 差值且 domain 覆盖 `<0.5` 的短序列；
2. `15-<20` bitscore 差值且 `>=0.8` 覆盖的长序列；
3. `marker_present` 与 `no_record` 成对 strata 的样本，确认邻域表是否存在系统性缺口。

固定样本中，低差值（`0-<2` 或 `2-<5`）且低覆盖（`<0.5`）有 82 条，
高差值（`15-<20`）且高覆盖（`>=0.8`）有 12 条；这两组应优先人工查看，
分别用于识别模型重叠造成的弱域命中，以及接近竞争边界但结构域完整的候选。
