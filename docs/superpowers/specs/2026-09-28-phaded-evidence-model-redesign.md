# PhaDED 证据模型与候选目录重构设计

## 目标

在不覆盖 2026-09-21 冻结结果、不直接启动 GTDB 全库重扫的前提下，把项目从“实验阳性决定能否建模、多个证据压成高可信二元标签”调整为“序列发现、序列家族归类、功能校准彼此独立”的证据体系。最终产品是一个尽可能完整、可追溯、按可靠度分层的 PhaDED 同源序列目录，而不是未经实验验证的降解表型目录。

## 已批准的设计假设

用户要求以序列发现的完整性和可靠性为主要目标，实验阳性的数量要求可以降低，但实验记录本身必须有更完整的证据。由此采用以下原则：

1. 实验阳性不是构建发现 HMM 或序列家族 HMM 的必要条件。
2. annotation-only 序列可进入发现或序列家族模型，但必须通过序列完整性、去重、结构域一致性和训练集合溯源检查。
3. 只有功能校准与功能特异性声明需要独立实验阳性、held-out、近邻阴性和 challenge 闭环。
4. 旧的 40,690、36,559、4,131 和相关 run 保持冻结；新规则产生 v2 输出和新 run，不回写历史结果。
5. 任何 GTDB 全库重扫、SignalP/结构/gRodon 重算、安装、服务器配置、删除和 push 均保留为显式授权动作。

## 四层模型体系

| 层级 | 训练来源 | 允许的结论 | 禁止的结论 |
|---|---|---|---|
| `reference_query_only` | 1–2 条参考或无法形成稳定比对的家族 | 参考相似性、人工复核入口 | 稳定家族 HMM、功能判定 |
| `discovery_hmm_uncalibrated` | 通过质量控制的实验与 annotation-only 参考 | 高召回同源候选、原始得分 | 删除/降级候选、family 判定、表型 |
| `sequence_family_hmm_validated` | 通过独立 held-out 和近邻混淆测试的序列家族模型 | sequence-defined family 归属 | PHB/PHA 降解功能确认 |
| `calibrated_candidate_model` | 高质量实验阳性、held-out、功能阴性和 challenge | 功能候选排序和有边界的模型验收 | 将每条命中写成实验表型 |

`discovery_hmm_uncalibrated` 名称和只召回约束保持不变。新增的 `sequence_family_hmm_validated` 是独立层，不进入 `pipeline/config/formal_scan_models.tsv`，除非以后另行通过发布决策。

## 证据字段拆分

候选记录不再把历史分类、转运预测、定位、基序和催化残基压成单一标签。v2 至少包含：

- `historical_superfamily`：Knoll/PhaDED 历史分类；
- `sequence_family_call` 与 `sequence_family_confidence`：序列家族归类；
- `catalytic_domain_type`：type 1/type 2/其他/未定，仅依据催化域几何和可靠参考比对；
- `accessory_domain_architecture`：SBD、linker、lid 等独立记录；
- `transport_signal_prediction`：SP/LIPO/TAT/TATLIPO/OTHER/未测试/输入不合格；
- `localization_evidence`：实验确认、其他证据、仅历史标签、未知；
- `nucleophile_identity`：Ser/Cys/unresolved；
- `motif_class`：GxSxG/AHSMG/Cys-associated/other/unresolved；
- `experimental_evidence_grade`：E3/E2/E1/A；
- `model_layer` 与 `functional_calibration_status`：序列模型状态和功能校准状态分开；
- `primary_disposition`：每条蛋白唯一最终处置；
- `evidence_flags`：允许多个非互斥支持或失败原因。

AHSMG 是含 Ser 的 motif class，不再作为第三种亲核体。SignalP 结果只描述预测转运状态，不自动等于实验定位。SBD 不参与 type 1/type 2 的定义。

## 候选目录

v2 主目录保留全部已进入分析范围的蛋白，每条蛋白只有一个 `primary_disposition`：

1. `core_sequence_homolog`：序列家族模型独立验证通过，覆盖和归属差值满足预注册标准；
2. `probable_sequence_homolog`：多项序列/结构域证据支持，但定位或唯一归属未完全闭合；
3. `remote_homolog_candidate`：仅发现层支持；
4. `function_unresolved`：同源关系可能成立，但无法排除竞争功能，例如 nPHAMCL-like；
5. `deferred_structure_review`：with-lipase 等需结构竞争面板；
6. `excluded_input_quality`：仅限序列缺失、明显截断或输入错误，不能把科学不确定性写成输入错误。

旧的 high-confidence 表仅作为 v1 比较输入。v2 不用一个二元“高可信”列混合不同家族的证据标准。

## 实验证据分级

- E3：纯化蛋白直接底物转化、产物或动力学证据，最好包含催化突变；
- E2：敲除/回补、细胞聚合物变化或强体内因果证据；
- E1：表达、相关表型或其他间接实验关联；
- A：数据库注释、序列推断或历史家族继承。

一条记录必须绑定 accession/version、序列 SHA-256、菌株/基因、论文、实验单位、底物、实验体系、结果、对照和独立性分组。重复 accession、同一基因的多个数据库记录以及同一实验的重复报道不能重复计为独立阳性。

一个 E2/E3 锚点足以让家族获得 `experimentally_anchored` 说明，但不足以自动获得 `calibrated_candidate_model`。现有 ≥3 独立阳性规则保留为功能校准治理门槛，并补充跨属/独立性字段的机器检查；它不再阻止发现型模型建设。

## 数据流

```text
冻结 reference ledger + FASTA
  -> 证据分级与训练资格审计
  -> reference_query_only / discovery_hmm_uncalibrated
  -> 在参考、held-out、近邻混淆集和既有候选上评估
  -> sequence_family_hmm_validated（仅通过者）
  -> v2 候选证据表
  -> 唯一 primary_disposition + 多 evidence_flags
  -> v1/v2 影响分析
  -> 经授权后才执行 SignalP、结构、gRodon 或全库重算
```

## 特殊家族

- `hfam_70`：状态改写为 `candidate_gate_passed_not_promoted`，保留未 finalize、未进 registry、未重扫边界。
- `hfam_52`：可作为序列模型候选，弱实验记录单独分级。
- `hfam_4` / nPHAMCL：降为 `function_unresolved`；现有 987 条不得与功能判别较强的集合混合展示。
- `hfam_2` / with-lipase：1,206,655 条改称“百万级待判别命中”，继续 `deferred_structure_review`，不得进入主候选计数，也不得删除。
- 只有 1–2 条参考的家族：`reference_query_only`，不把先验主导的单/双序列 profile 表述为稳定家族 HMM。

## HMMER 分片尺度

正式入口必须记录完整目标序列数 `Z`。如果使用分片搜索，则每个 shard 使用相同的全库 `-Z`，并在 manifest 中保存 `database_size_Z`、shard 序列数和命令。对 scan 13 先做只读命令和日志审计：若缺少 `-Z` 且原始 bitscore/shard 计数足够，则创建新 reconciliation run 做可证明的尺度统一；若不足，则停止并请求全库重扫授权，不能把未经校正的 shard E-value 当全库统一尺度。

## 下游统计

gRodon 的组名改为“候选基因携带”与“在规定搜索/质量条件下未检出候选”。必须断言最终阳性基因组与对照集合零重叠，验证 `growth_rate_per_h` 的公式和单位，并把属水平平均差设为主估计对象。旧的 713 条复用预测只有在工具、参数、输入和公式可比时进入主分析，否则只进入敏感性分析。

## 验收边界

设计完成的最低标准是：

1. 旧 run 和旧结果没有被覆盖；
2. 每条 v2 候选有一个且只有一个最终处置；
3. AHSMG/Ser、SignalP/定位、SBD/type、序列模型/功能模型已经拆分；
4. 52 条 demotion 在规则层可重复产生；
5. 616 条未解释流转、4,301/6,214 集合关系和 gRodon 66 条差额有机器可读解释；
6. HMM/profile 绑定训练集合、alignment、HMM、工具版本和 SHA-256；
7. 发现层不筛除候选、不产生 family call；
8. 所有新运行使用新 `runs/<run_id>/` 和 dated `deploy/<run_id>/`，线程遵守动态上限；
9. 测试、`compileall` 和 `git diff --check` 通过；
10. 正式发布前完成本地、GitHub、deploy 和 run manifest 权威核对。
