# T141 最终方案 A 重算交接清单

## 目的与边界

本清单只定义下一次 T141 运行需要验证和执行的工作，不包含远程写入操作。它用于使生态分布、
位点级基因邻域和 OH 系统树与同一最终方案 A 输入绑定。所有产物仅支持候选同源蛋白或功能潜力，
不构成 PHB 降解表型证据。

## 执行前门禁

1. 在 `/home/data/haoyu/PHB_gtdb-ds` 确认实际脚本为服务器的扁平 `scripts/` 结构；不要把本地
   `pipeline/scripts/` 的根路径规则直接复制过去。
2. 记录服务器 `scripts/`、`data/hmms/v2/`、最终 `data/screen/hits_filtered.tsv`、四个 tier1 FASTA
   以及 GTDB taxonomy/metadata 文件的 SHA-256、路径、大小和修改时间。
3. 确认 `hits_filtered.tsv` 含 `family`、`genome`、`locus` 列；确认 OH tier1 为 1,429 条。
4. 将本地修复以显式版本包部署到服务器后，重新核验脚本 SHA-256；不要假定同名脚本内容一致。
5. 以新的空运行目录执行，保留旧 `results/` 只读，不覆盖历史表或树。

## 必须重跑

1. 以最终方案 A 的位点表执行 `11_clusters.py`，输出新 `cluster_context.tsv` 和五列
   `cluster_summary.tsv`：`hit_family`、`marker_family`、`marker_hits`、`supporting_loci`、
   `supporting_genomes`。
2. 以同一最终 `genome_hits.tsv` 执行 `10_distribution.py`，并把输入 SHA-256 写入运行清单。
3. 从当前 OH tier1 FASTA（1,429 条）重建 OH 树；记录完整命令、抽样规则、输入 SHA-256、
   tree SHA-256 和叶数。旧 1,465 叶树不得覆盖或重命名为新树。
4. 为 ePhaZ CD-HIT 树登记去冗余输入 FASTA 的生成命令、路径、SHA-256、序列数与树 SHA-256；
   若无法获得这些证据，该树继续保持 `input_not_registered`。
5. 用失败即停止的编排生成新的 `run_manifest.jsonl` 和 `run_manifest.json`。清单应绑定步骤退出码、
   输入/输出 SHA-256；任一声明文件缺失或为空必须失败。

## 验收条件

- `run_manifest.json` 中所有步骤成功，且没有缺失输入或输出。
- 新 `cluster_summary.tsv` 为五列新 schema；任何统计同时报告出现次数、唯一位点和唯一基因组，
  不把三者混用。
- `tree_manifest.tsv` 中 OH 不再是 `stale_input`；ePhaZ CD-HIT 的状态可由登记证据决定。
- 生态、cluster 和树的输入哈希均可追溯至同一方案 A 候选集。
- 完成后才更新 `docs/STATUS.md`、两份审核报告和最终结果报告中的历史标记。
