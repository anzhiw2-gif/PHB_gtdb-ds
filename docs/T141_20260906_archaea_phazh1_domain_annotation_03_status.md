# T141 古菌 PhaZh1-like domain annotation 03

## 状态

- Run：`20260906_archaea_phazh1_domain_annotation_03`
- 状态：`completed_candidate_only`
- T141 run：`${PHB_REMOTE_ROOT}/PHB_gtdb-ds/runs/20260906_archaea_phazh1_domain_annotation_03`
- 绑定 deploy：`${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260906_archaea_phazh1_domain_annotation_03`
- 执行环境：`${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb`，`HMMER 3.4 (Aug 2023)`
- 不修改 `pipeline/config/formal_scan_models.tsv`，未启动正式 GTDB scan。

## 输入和交叉注释

输入为 audit_04 中的 204 条古菌记录：`PhaZh1_like_high=197`、`PhaZh1_like_review=7`。候选 FASTA、正控和负控均有路径、大小与 SHA-256 绑定在 `input_contract.json` 和 deploy manifest 中。

本次执行 18 个 HMMER 搜索（6 个模型乘 candidate、positive_control、negative_control），使用 `i-Evalue <= 1e-5`。逻辑模型为 `ArchPhaZ_patatin`、`ArchPhaZ_hydrolase`、`PhaJ`、`BdhA`、`PhaC` 与 `phasin`。

| 项目 | 结果 |
|---|---:|
| candidate domain annotation | 204/204 records |
| `domain_supported` | 204 |
| `domain_conflict` | 0 |
| `no_registered_domain_hit` | 0 |
| positive-control hits | 8 |
| negative-control hits | 5 |

所有 candidate 仅命中 `ArchPhaZ_patatin`；没有命中本轮登记的 `ArchPhaZ_hydrolase`、`PhaJ`、`BdhA`、`PhaC` 或 `phasin` 模型。正控 8/8 命中，但负控也有 5 条 `ArchPhaZ_patatin` 命中，因此该面板不是零交叉反应的阴性证明，不能用于提高候选为 PHB 表型。

独立 Pfam/InterPro 数据库或工具在 T141 上未登记，字段保持 `pending_no_registered_pfam_or_interpro_database`。

## 跨门人工复核

`cross_phylum_manual_review.tsv` 包含 9 条小规模人工复核子集：Halobacteriota review 2 条、Thermoproteota high 2 条和 review 5 条。每条保留 genome/locus、门/纲、coverage、E-value、N 端 `GxSxG`、模型命中、contig/坐标/链方向、两来源邻域标记和 `candidate_only` 处置。

所有 9 条的已记录邻域标记均为 `BdhA`，无已登记 `PhaC` 或 phasin 冲突命中。这是邻域/同源背景信息，不是降解表型证据。

## 结果与保留的失败证据

- 完成 manifest：`runs/20260906_archaea_phazh1_domain_annotation_03/results/completion_manifest.json`
- 主结果 SHA-256：`866cea0bc5c9b98828ee9e9c5538d181af6dbfffe5412babcc922fa453979a58`
- calibration SHA-256：`af5070da22793dda1474bfcdde4c60508b2cd402ae63f28733aa7f6228fe9537`
- 人工复核表 SHA-256：`869f3544d51dd668af16b9168ef664b8d94982ccf005db12a0870be3d8cc4328`

`20260906_archaea_phazh1_domain_annotation_01` 在 BOM header 检查阶段失败，`_02` 在 `phasin` SHA-256 检查阶段失败；两者均在 HMMER 前 fail-closed，服务器目录和日志保持原样。

## 解释边界

本结果只说明 PhaZh1-like candidate 与本项目已登记的 patatin 模型存在同源/结构域层面的一致性。HMM、domain、motif、邻域和分类信息都不等同于已验证的 PHB 降解表型。
