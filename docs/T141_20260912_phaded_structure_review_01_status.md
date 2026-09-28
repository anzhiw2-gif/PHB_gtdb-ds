# T141 PhaDED 结构审查输入冻结

日期：2026-09-12  
运行：`runs/20260912_phaded_structure_review_01/`  
状态：`planned_not_run`

## 已完成

根据四面板系统发育结果，冻结 19 条高信息候选：

- 16 条 `retain_dPHASCL1_like_candidate`
- 3 条 `retain_extracellular_dPHASCL2_like_candidate`

每条候选均绑定了最近 PhaDED 参考序列，共冻结 19 条 reference controls。候选 FASTA、选择表和参考 FASTA 均已记录 SHA-256，合同见 `runs/20260912_phaded_structure_review_01/input_contract.json`。

## 工具检查

绑定的 T141 `phb_gtdb` 环境中未发现以下可执行文件；Python 包可用性尚未绑定或验证：

- ColabFold
- AlphaFold
- ESMFold
- Boltz
- Foldseek

因此本 run 没有启动结构预测、结构检索或结构比较，`structure_predictor.status=pending_tool`。没有安装软件、修改环境或连接外部结构服务。

## 解释边界

后续获得的预测结构、AlphaFold 模型或同源结构只能用于检查折叠、催化位点几何关系、信号肽后成熟区和结构域边界；它们不能单独证明 PHB/PHA 降解表型，也不能直接修改 `ePhaZ_curated_core` 或正式扫描模型。

## 下一步授权点

需要单独决定结构来源：

1. 在 T141 安装并绑定可复现的结构预测环境；或
2. 使用外部结构服务/数据库，并为每条候选记录查询日期、模型版本、结构文件和 SHA-256。

获得结构后，先审查 19 条高信息候选；12 条 PhaC-vs-PhaZ 和 8 条结构异常仍保持 challenge/hold，不得并入 ePhaZ 阳性训练集。
