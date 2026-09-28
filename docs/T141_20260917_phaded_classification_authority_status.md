# T141 PhaDED 分类权威口径（层级错配修复）

日期：2026-09-17
Run：`runs/20260917_phaded_classification_authority_01`
计划：`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` §3 Task 1
边界：本轮只做本地、只读、candidate-only 口径固化；没有启动服务器任务，没有修改 `pipeline/config/formal_scan_models.tsv`、历史 HMM、registry、正式扫描结果或任何历史 run，没有删除任何文件。

> **candidate-only 声明：** 本项目的 HMM、profile、domain、motif、SignalP、结构与系统发育结果只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。本 Task 建立的分类权威表同样只是层级口径，不构成任何表型结论。

---

## 1. 问题：层级错配（layer mismatch）

项目此前把 PhaDED 的两层实体当作同一类东西处理：

| 层级 | 实体 | 来源性质 |
|---|---|---|
| 超家族 | 8 个 | **功能/定位先验**（硬框架） |
| 同源家族 | 38 个 | **2009 年序列相似性 + 系统发育聚类**（次级细化层） |

两者混用会导致一个具体风险：把 2009 年的相似性聚类当作与功能先验同等的证据，从而让一个纯聚类层获得 registry 资格。本 Task 用一张显式权威表把两层拆开，并把"聚类层不得进入 registry"写成机器可检验的不变量。

---

## 2. 文献依据（逐条可复核）

来源：Knoll M, Hamm TM, Wagner F, Martinez V, Pleiss J. *The PHA Depolymerase Engineering Database.* BMC Bioinformatics 2009;10:89. PMCID PMC2666664；原文全文已在库内：`research/europepmc/comp_ft_knoll_phaded.txt`。

1. **8 个超家族来自 28 条实验验证种子序列。** 原文正文（“Database construction”）："As a first step, **28 seed sequences of proteins with experimentally validated depolymerase activity** (Table 1[15-33]) were stored in the database and annotated. These seed sequences were assigned to 6 previously described superfamilies **based on their function** [34]. Additionally the families of intracellular nPHASCL depolymerases (lipase box) … and the family of periplasmatic PHA depolymerases … were introduced. **Thus, a total of 8 superfamilies were introduced**."
   → 超家族层是**功能先验**构成的硬框架，不是聚类产物。

2. **38 个同源家族来自序列相似性加系统发育聚类。** 原文："**Superfamilies were subdivided into homologous families, which were introduced based on sequence similarity and phylogenetic analysis (Fig. 1).**" 数据库总体量为 "735 sequence entries which code for 587 different proteins … assigned to **8 superfamilies and 38 homologous families**"。
   → 家族层是**次级细化层**，其证据类型是聚类，与超家族层不同类。

3. **DED 自 2009 年起未再更新。** Database Commons 记录（ID 3446，<https://ngdc.cncb.ac.cn/databasecommons/database/id/3446>）：`Last update: 2009`，版本 v1.1；该记录经库内既有调研转录：`research/web/web_survey_report.md` L31、`docs/literature_survey_report.md` L333。本轮**未重新联网核对**该页面，故此处按"库内既有记录 + 计划 §1.2 基线"引用，不声称本轮实测。
   → 家族层绑定的是**冻结快照** `Knoll2009_v1.1_frozen`，不是持续更新的数据库。

### 2.1 三处注释与序列归属矛盾 —— "不得以注释定 family"的原始依据

Knoll 原文明确记录了三处 **GenBank 注释与 DED 序列归属相反**的情形，这是本项目"注释不能作为分类依据"的直接文献依据：

| # | 序列 | GenBank 注释 | DED 归属（按序列相似性） | 原文要点 |
|---|---|---|---|---|
| 1 | `gi:194292521`（*Cupriavidus taiwanensis*） | "intracellular PHA depolymerase" | **extracellular dPHASCL type 1** | "…were assigned to the family of extracellular dPHASCL depolymerases (catalytic domain type 1) **due to their sequence similarity**" |
| 2 | `gi:74267419`（*Ralstonia eutropha* H16） | "intracellular PHA depolymerase" | **extracellular dPHASCL type 1** | 同上；且该酶"highly active against artificial amorphous PHB granules, and is lacking a signal peptide, a linker domain, and a substrate binding domain" |
| 3 | `gi:34452171`（*Pseudomonas* sp.） | "**extracellular** PHA depolymerase" | **intracellular nPHAMCL** | "Another exception is the PHA depolymerase from Pseudomonas sp. which is annotated as 'extracellular PHA depolymerase' in the GenBank but was assigned to the family of intracellular nPHAMCL depolymerases in the DED" |

**定位例外（第 4 处，非注释矛盾但同类性质）：** PhaZ7（*Paucimonas lemoignei*）是分类上的定位例外——胞外酶但**只对天然颗粒有活性**（原文："One exception of this classification is an extracellular nPHASCL depolymerase from Paucimonas lemoignei which is active only against native PHA granules [6]"）。
→ 其所属超家族 `extracellular native-SCL/PhaZ7-like` 在权威表中的 family 行已显式带上该例外说明。

**结论：** 注释（含 GenBank 与自动注释）既可能把胞外酶标成胞内，也可能把胞内酶标成胞外。因此 `sequence_clustering_2009` 层的归属只能由序列相似性 + 系统发育支撑，**注释不得上升为 family 判定依据**，更不得赋予 registry 资格。

### 2.2 诚实记录的文献内不一致（本轮新发现）

Knoll 2009 **摘要**写的是 "587 PHA depolymerases, which were assigned to 8 superfamilies and 38 homologous families **based on their sequence similarity**"，与正文的分层叙述（超家族=功能先验；家族=聚类细分）措辞不同。

- 本 Task 的权威表**采用正文口径**，因为：正文给出了可核对的构建步骤（28 条种子 → 6 个按功能划分的既有超家族 + 新增 2 个 = 8 个超家族；随后"superfamilies were subdivided into homologous families"）。
- 摘要的宽松措辞**不得**被用来论证"家族层也是功能先验"。若要推翻本 Task 的分层，需要新的原始构建记录，而不是摘要措辞。

---

## 3. 对 109,087 条候选的数值影响：**0 行**

> **本 Task 对 109,087 条候选为 0 行数值影响。**

不是估计，而是有三条可复核的理由：

1. 8 个超家族标签的 `registry_eligible` **全部保持 `true`**，因此没有任何超家族范围的判定因本 Task 降级。表中 `registry_eligible=true` 计数 = 8。
2. `subtype_call` **只由 `phaded_superfamily_best` 决定**（`pipeline/scripts/build_phaded_subtype_matrix.py:_base_subtype`，L153–165 只读 `phaded_superfamily_best`）。`phaded_family_best` 不参与任何候选调用；把 family 层设为 `registry_eligible=false` **在计算路径上无法改变任何一行**。
3. 在本 Task 之前，`registry_eligible` 在整个仓库（排除本 Task 新增文件）**只有 1 处出现，且是计划文档本身**（`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` L105 的 Step 1.1 定义）；`pipeline/` 下出现次数为 **0**。即该字段此前**没有任何运行时消费者**，本表目前是治理产物，无计算副作用。

候选宇宙实测（本轮实测，非引用计划）：

| 项 | 值 |
|---|---|
| `runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv` 数据行数 | **109,087**（含表头 109,088 行） |
| 该文件 size | 64,498,737 |
| 该文件 SHA-256 | `796cb0a32cc7fde1c768b2c0f2f8e30542c7aa37e773f2b8004a482ac878569a` |
| 本 Task 改变的行数 | **0** |

---

## 4. 权威表

路径：`runs/20260917_phaded_classification_authority_01/inputs/phaded_classification_authority.tsv`
字段（顺序即声明顺序）：`layer_id, layer_kind, source_evidence_type, source_version, registry_eligible, gate_profile_id, notes`

- `layer_id`：superfamily 行用 DED 超家族标签原文；family 行用 `phaded_family_id`（`DED_hfam_<整数>`）。
- `gate_profile_id`：绑定 `profile_manifest.tsv` 中的真实 profile id（family/superfamily profile），可反查训练与校准状态。

### 4.1 实测计数

| 口径 | 实测值 |
|---|---:|
| superfamily 层行数 | **8** |
| family 层行数 | **38** |
| 总行数 | **46** |
| `source_evidence_type=functional_prior` | 8 |
| `source_evidence_type=sequence_clustering_2009` | 38 |
| `registry_eligible=true` | **8** |
| `registry_eligible=false` | **38** |
| `source_version=Knoll2009_v1.1_frozen` | **46（全部）** |
| `gate_profile_id` 解析成功 / 总数 | **46 / 46（pending=0）** |

`gate_profile_id` 覆盖：46 个 profile 中 **9 个 trained + 37 个 reference-only**（33 family + 4 superfamily），本轮从 `profile_manifest.tsv` 实测：family trained=5 / family reference_only=33 / superfamily trained=4 / superfamily reference_only=4，与计划 §1.1 基线一致。

### 4.2 核心不变量（19 项，全部通过）

| 不变量 | 内容 |
|---|---|
| `superfamily_count_is_8` | 超家族层恰好 8 行 |
| `family_count_is_38` | 家族层恰好 38 行 |
| `layer_ids_are_unique_and_non_empty` | `layer_id` 是唯一非空主键 |
| `superfamilies_are_functional_prior` | 8 行全部 `functional_prior` |
| `superfamilies_are_registry_eligible` | 8 行全部 `registry_eligible=true` |
| `families_are_sequence_clustering_2009` | 38 行全部 `sequence_clustering_2009` |
| `families_are_bound_to_the_frozen_knoll2009_v1_1_snapshot` | 38 行 `source_version` 全部为 `Knoll2009_v1.1_frozen` |
| **`families_are_not_registry_eligible`** | **38 行全部 `registry_eligible=false`（核心不变量）** |
| `evidence_type_and_registry_eligibility_are_never_crossed` | 只允许 `(functional_prior,true)` 与 `(sequence_clustering_2009,false)`；交叉即报错 |
| `every_layer_is_bound_to_the_frozen_ded_snapshot` | 46 行同一冻结快照 |
| `every_layer_carries_its_provenance_rationale` | 每行 `notes` 非空 |
| `superfamily_layer_ids_match_ded_import` | 与 `phaded_family_definitions.tsv` 无新增无缺失 |
| `family_layer_ids_match_ded_import` | 同上（38 个 family id 一一对应） |
| `superfamily_layer_ids_match_hard_framework` | 等于 8 个 DED 超家族标签原文 |
| `family_layer_ids_use_ded_hfam_format` | 全部符合 `DED_hfam_<整数>` |
| `gate_profiles_resolve_in_profile_manifest` | 46/46 解析成功，pending 即 fail |
| `gate_profiles_are_unambiguous` | 无 profile 归属歧义 |
| `family_gate_profiles_confirm_superfamily_parentage` | 每个 family 的 gate profile 指向与 DED 导入相同的父超家族 |
| `candidate_universe_row_count_matches_plan_fact_baseline` | 候选宇宙实测 = 109,087 |

### 4.3 源文件哈希（全部实测）

| 角色 | size | SHA-256 | 状态 |
|---|---:|---|---|
| `phaded_family_definitions` | 16,879 | `9c9a21a778390873b6377903f90bb29a400747351081aac570f60b4fcbeea62a` | verified |
| `profile_manifest` | 13,721 | `50f584eab55516c7874dbedaed459b366f6fbed9b8be63a1844aabbc1a4b42c6` | verified |
| `phaded_classification_authority`（本 Task 产出） | 37,670 | `a978369172a070e73b7b469766710ad503e4a7f60bef62e1e213e11e0843ee31` | verified |
| `candidate_evidence` | 64,498,737 | `796cb0a32cc7fde1c768b2c0f2f8e30542c7aa37e773f2b8004a482ac878569a` | verified |

**`phaded_family_definitions.tsv` 的 SHA-256 与计划中给定的期望值 `9c9a21a778390873b6377903f90bb29a400747351081aac570f60b4fcbeea62a` 逐字符一致（实测，match=True）。** 该文件的两个上游副本（`runs/20260909_phaded_reference_evidence_01/inputs/`、`runs/20260915_phaded_profile_gap_audit_01/inputs/`）以及本 run 内的副本哈希全部等于该值。

`profile_manifest.tsv` 在仓库中共有 6 份副本，**6 份 SHA-256 全部等于 `50f584eab55516c7874dbedaed459b366f6fbed9b8be63a1844aabbc1a4b42c6`**（测试逐份核对），因此不存在"用错那份 manifest"的风险。

---

## 5. 报告口径统一（Step 1.6）

以下三条自本 Task 起对**所有新输出**生效；本 Task 自身产出已遵守。

### 5.1 版本标注

所有引用 PhaDED 分类的新输出必须标注：

```
PhaDED (Knoll 2009, v1.1, frozen snapshot)
```

并附**导入哈希**（本次为 `phaded_family_definitions.tsv` = `9c9a21a7…ea62a`、`profile_manifest.tsv` = `50f584ea…b42c6`）。禁止使用不带快照哈希的裸 "PhaDED" 或 "DED" 表述，避免与任何"已更新版本"混淆。

### 5.2 无 SBD 家族：Pfam PF06850 缺失是结构性预期，不是阴性

涉及层级（权威表中已写入 `notes`）：

- `extracellular dPHAMCL`（family `DED_hfam_8`）
- `extracellular native-SCL/PhaZ7-like`（family `DED_hfam_7`）

声明：这两层**本来就没有底物结合域（SBD）**（Knoll 2009：*dPHAMCL depolymerases possess no substrate binding domain*；PhaZ7 型同样无 SBD），因此 **Pfam PF06850（PHB 解聚酶 C 端 SBD）缺失是结构性预期，不得计为阴性结果**，也不得解释为"证据不足"。凡是把 PF06850 缺失当作阴性证据的表述都必须改写。

### 5.3 `dPHASCL1_like_candidate` 与 `type1_verified` 必须分开表述

- `dPHASCL1_like_candidate`：**超家族候选标签**（由 `phaded_superfamily_best` 派生），只表示"与 extracellular dPHASCL type 1 超家族相似"。本轮在 `runs/20260915_phaded_subtype_reconciliation_01/results/phaded_final_subtype_evidence_ledger.tsv` 上**实测**为 67,342 / 109,087 行（该台账总行数 109,087，与候选宇宙一致）。
- `type1_verified`：**几何验证过**的 type 1 归属（需同时具备 lipase box 与 oxyanion hole 坐标，且 oxyanion 位置相对 lipase box 满足 type 1 关系）。

两者**不得互相替换、不得合并计数**。特别是：不得把 67,342 写成"已验证 type 1"，也不得在未做几何验证时使用 `type1_verified` 这一措辞。该量化由计划 Task 4 负责，本 Task 只固化口径。

---

## 6. 本轮完成清单

- ✅ Step 1.1 失败测试先行：`pipeline/tests/test_validate_phaded_classification_authority.py`（28 tests）
- ✅ Step 1.2 确认失败：`FileNotFoundError: pipeline/scripts/validate_phaded_classification_authority.py`；`Ran 28 tests … FAILED (failures=11, errors=14)`
- ✅ Step 1.3 实现 fail-closed 校验器：`pipeline/scripts/validate_phaded_classification_authority.py`（输出计数 + 全部源文件哈希；任一不变量失败即 `status=blocked` 且退出码 1；`pending` gate profile 视为失败）
- ✅ Step 1.4 生成权威表并让测试通过：`Ran 28 tests … OK`（exit 0）
- ✅ Step 1.5 本状态文档（文献依据 + 三处注释矛盾 + 0 行影响声明）
- ✅ Step 1.6 报告口径统一（§5）

**测试同时覆盖 fail-closed 行为**（不是只测"文件存在"）：family 被提权为 registry eligible、superfamily 被降级为聚类层、`source_version` 被改成非冻结值、多出/缺失 layer id、缺列、空 `registry_eligible` —— 7 种篡改输入全部必须被拒绝。

## 7. 主要产物

- `runs/20260917_phaded_classification_authority_01/inputs/phaded_classification_authority.tsv`
- `runs/20260917_phaded_classification_authority_01/inputs/phaded_family_definitions.tsv`（上游副本，哈希一致）
- `runs/20260917_phaded_classification_authority_01/inputs/profile_manifest.tsv`（上游副本，哈希一致）
- `runs/20260917_phaded_classification_authority_01/results/authority_validation.json`
- `runs/20260917_phaded_classification_authority_01/input_contract.json`
- `runs/20260917_phaded_classification_authority_01/logs/commands.txt`
- `pipeline/scripts/validate_phaded_classification_authority.py`
- `pipeline/tests/test_validate_phaded_classification_authority.py`
- 本文件

### 7.1 input_contract.json 的 `status` 说明（避免误读）

契约顶层 `status` 为 `pending`，**唯一原因**是三个 GTDB 槽位（taxonomy / metadata / tree）未解析——本 Task 不消费任何 GTDB 文件。契约内 `run_inputs_status=verified`，且已声明的 3 个 run 输入（family definitions、profile manifest、authority table）全部 `status=verified` 并记录了 size 与 SHA-256。此处**没有**把缺失项伪造成已完成。

## 8. 限制与未做事项（诚实记录）

1. **无运行时消费者。** 权威表目前不被任何 pipeline 脚本读取；它是治理与口径产物。接入下游属于后续任务，且接入时必须另开 run。
2. **文献内措辞不一致已如实记录**（§2.2），未隐去摘要与正文的差异。
3. **28 vs 30 不得混淆。** Knoll 2009 的 28 是**文献中 Table 1 的实验验证种子序列数**（DED 建库用）；本项目 `phaded_reference_ledger.tsv` 的 `experimental_positive` 记录为**30 条**（723 行 = 30 experimental_positive + 693 annotation_only），超家族 profile 绑定的阳性 accession 合计也是 30。两者是不同量，本 Task **不**声称项目超家族层由那 28 条文献种子构成。
4. **Database Commons 的 "Last update: 2009, v1.1" 未在本轮联网复核**，按库内既有调研记录与计划基线引用（§2 第 3 条）。
5. 未修改 `pipeline/config/formal_scan_models.tsv`、历史 HMM、历史 run；未执行任何服务器任务；未删除任何文件。

## 9. 验证命令与结果

```
python -m unittest pipeline.tests.test_validate_phaded_classification_authority -v
  -> Ran 28 tests in 0.150s ... OK   (exit 0)
python -m unittest discover -s pipeline/tests
  -> Ran 346 tests in 8.258s ... OK (skipped=1)   (exit 0,连续两次一致)
python -m compileall -q pipeline
  -> exit 0（无输出）
git diff --check
  -> exit 0（无空白错误）
```

完整套件实测 **346 tests OK (skipped=1)**（本 Task 后，连续两次运行结果一致）；基线（计划 §1.1）为 277，差值来自本 Task 新增 28 项与并行 agent 的新增测试。**本 Task 自身的 28 项测试全部通过。**

**本轮观察到的一次瞬时失败（如实记录，未修改）**：在一条较早的完整套件运行中，曾出现 `FAILED (failures=5, errors=1, skipped=1)`，6 条诊断全部落在 `test_assign_phaded_catalytic_domain_type.py`（计划 Task 4，**属于其他并行 agent 的测试文件**，落盘时间 2026-09-17 1:03:46，当时仍在写入中）。按分工要求**只报告、不修改**；随后两次完整运行该模块均已通过，因此判定为并行写入过程中的瞬时状态，与本 Task 无关。

**历史产物未改动核对**：`git status --porcelain` 中唯一的已跟踪改动是 ` M AGENTS.md`（1 行删除），该改动**在本会话开始前即已存在**，非本 Task 产生；`runs/`、`results/`、`pipeline/config/` 下无任何已跟踪文件被修改或删除。受保护文件哈希本轮末尾复测：`phaded_family_definitions.tsv`（上游）= `9c9a21a7…ea62a`（未变）、`pipeline/config/formal_scan_models.tsv` = `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`（未变）、`motif_candidate_evidence.tsv` = `796cb0a3…8569a`（未变）。

---

*本状态文档为 candidate-only 研究记录；所有层级口径表示候选同源或功能潜力，不等同于已验证 PHB/PHA 降解表型。*
