# PhaDED `intracellular nPHASCL (no lipase box)`（Cys 型）判据 decision record

**任务：** 计划 v2 Task 6（`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` §3 Task 6）
**Run：** `runs/20260917_phaded_inphascl_cys_decision_01`（新建，未触碰任何历史 run / `results/` / `pipeline/config/` / 历史 HMM）
**日期：** 2026-09-17
**状态：** decision record 完成；推荐方案为**方案 A + 候选层注释证据层**，该证据层已实现并落盘（`status=completed_candidate_only`）
**服务器：** **未使用**。本轮仅在 2026-09-17T01:46:02+08:00 做了一次只读负载核对（见 §5.3），未提交任何服务器任务，占用 0 线程

> **证据边界（贯穿全文）：** 本记录与本次 run 的全部产出只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。本任务**不产生任何新的 family 判定**，不修改 `pipeline/config/`，对 109,087 条候选的 `subtype_call` 影响为 **0 行**。

---

## 1. 问题与范围

`intracellular nPHASCL without lipase box` 是 PhaDED 中最大的胞内超家族（参考台账 276 行）。其判据在文献中被描述为"极硬"：

- 催化三联体为 **Cys-His-Asp**（催化亲核体是 Cys 而不是 Ser）；
- **没有** `Gx1Sx2G` lipase box；
- 催化 Cys 前一位（cysteine-1）几乎全为疏水残基，且几乎全为 **Val**。

本轮需要回答的是：**是否应当、以及是否可能**为它建立一个正式的 profile / registry 条目，或者是否存在另一条不违反项目硬规则的证据路径。要求比较三个方案（保持现 registry / 新增窄范围 Cys profile 并保留历史模型 / 版本化替换旧模型），预注册证伪条件，给出明确推荐，并在推荐为"不需要新 HMM"时实现候选层注释证据层。

---

## 2. 事实基线（全部本轮独立复核，可复核）

### 2.1 参考台账事实 —— 含对任务前提的一处更正

从冻结的 `runs/20260915_phaded_motif_reconciliation_01/results/reference_motif_panel.tsv`（723 行）与本 run 独立重算：

| 家族 | 参考数 | evidence_status |
|---|---:|---|
| `DED_hfam_61` | 42 | annotation_only 42 |
| `DED_hfam_62` | 5 | annotation_only 5 |
| `DED_hfam_63` | 64 | annotation_only 64 |
| `DED_hfam_64` | 11 | annotation_only 11 |
| **`DED_hfam_65`** | **73** | **annotation_only 72 + experimental_positive 1** |
| `DED_hfam_66` | 2 | annotation_only 2 |
| `DED_hfam_67` | 45 | annotation_only 45 |
| `DED_hfam_68` | 33 | annotation_only 33 |
| `DED_hfam_69` | 1 | annotation_only 1 |
| **Cys 超家族合计** | **276** | **annotation_only 275 + experimental_positive 1** |

> **前提更正（必须记录）：** 任务书写作"hfam_61 到 hfam_69 的参考全部是 annotation_only，实验阳性为 0"。**实测该表述对 9 个家族中的 8 个成立，对 `DED_hfam_65` 不成立。**
> `DED_hfam_65_0001` = **`CAJ92291.1`**，`evidence_status=experimental_positive`、`positive_negative_control=positive_seed`、`organism=Cupriavidus necator`、`reported_localization=intracellular`、`substrate_class=SCL`、`experimental_assay="Knoll 2009 Table 1 seed designation"`、`notes` 记 `Knoll 2009 Table 1 seed GI 3641686`（来源文件 `runs/20260911_phaded_pfam_architecture_01/inputs/phaded_reference_ledger.tsv` 第 124 行）。
> 这正是任务书中点名的模型酶 H16 PhaZ1 / GI 3641686 —— 也就是说任务书的两处陈述互相矛盾，冻结台账给出了裁决：**Cys 超家族恰有 1 条实验阳性**。
> 该更正**强化**而不是削弱下列结论：可训练的严格记录只有 1 条，仍远低于 gate 的 `positive_count >= 3`。同时须诚实记录：这条阳性的 `experimental_assay` 是 Knoll 2009 Table 1 的种子指定，即**文献策展结果**，不是本项目新做的实验。

### 2.2 硬规则与 gate（代码级依据）

- **禁训规则：** `pipeline/scripts/validate_phaded_reference_ledger.py` L172–173
  `if strict_training and status in {"annotation_only", "pending_review"}: raise ValueError(...)`；
  由 `pipeline/tests/test_validate_phaded_reference_ledger.py::test_rejects_annotation_only_record_from_strict_training` 强制。
  → 276 条 Cys 参考中 **275 条（99.6%）在严格训练中不可用**。
- **校准 gate：** `pipeline/scripts/finalize_phaded_reference_panel_calibration.py`（`gate = kind=="family" and positive_count >= 3 and heldout >= 1 and negative >= 1 and challenge >= 1 and family_resolved_negative_status == "sufficient" and challenge_status == "sufficient"`，fail-closed）。
  → 即使把 1 条实验阳性算作 `positive_count`，仍为 **1 < 3**；且 held-out = 0。
- **阴性面：** Phase 2 验证报告 §6 + Task 5 产物：960 条挑战记录中 `formal_negative_eligible` = **0**；冻结台账 723 行里没有任何 `experimental_negative`（693 annotation_only + 30 experimental_positive）。"台账绑定 + 显式实验阴性"这一组合在现有证据下**不可能存在**。
- **注册权威层级：** `runs/20260917_phaded_classification_authority_01/inputs/phaded_classification_authority.tsv`：8 个 superfamily = `functional_prior` + `registry_eligible=true`；38 个 family = `sequence_clustering_2009` + `registry_eligible=false`。Cys 超家族的 9 个 hfam 全部落在**不可注册**的那一层。

### 2.3 候选层事实（决定性）

- 109,087 条候选中 **0 条**出现在 723 行参考台账里，**0 条**属于任何冻结 DED 比对（Task 7 `candidate_ded_alignment_members = 0`）。
- 候选层从来没有任何一条带 Cys 先验：`motif_candidate_evidence.tsv` 的 `motif_expected_catalytic_residue` **109,087/109,087 全为 `Ser`**；`catalytic_ser_cys_state` 只有 `supported`(77,302) 与 `not_detected_pattern`(31,785) 两态。
  → **候选层无法用比对列传递任何催化 Cys 坐标**（Task 2 遗留 blocker：Cys 型亲核体坐标不可得）。
- Task 7 全库层：lipase box `supported 77,302 / not_detected_pattern 31,785`；PF06850 原始结合 `detected 30,621 / not_detected_in_tested_pfam 67,798 / not_tested_no_pfam_evidence 10,668`。

### 2.4 DEFECT_D 与它打开的那条路

DEFECT_D：在本参考集中 **PF06850 不是胞外 SBD**，而是胞内 Cys 型超家族的标记。本 run 独立重算 `reference_mapping_manifest.tsv` 的 `pfam_accessions`：

| 超家族 | 参考数 | 携带 PF06850 |
|---|---:|---:|
| intracellular nPHASCL without lipase box | 276 | **274（99.28%）** |
| extracellular dPHASCL type 1 | 284 | **0** |
| extracellular dPHASCL type 2 | 71 | **0** |
| 其余 5 个超家族合计 | 92 | **0** |

→ PF06850 在本参考集内是**该超家族的专属标记**（274/276，其余 447 条全为 0）。这提示存在一条**不需要训练新 HMM** 的候选层注释路径：`PF06850 命中 + 无 lipase box + Cys 上下文序列特征`。

---

## 3. 三个方案

### 方案 A：保持现 registry 不变，仅发布 candidate-only 说明（+ 本轮的候选层注释证据层）

| 维度 | 评估 |
|---|---|
| **治理可行性** | ✅ **完全合规**。不改 registry、不改 `pipeline/config/`、不训练模型、不启动服务器任务；新增产出只写新 dated run。候选层注释证据表**只新增证据列**，`subtype_call` 影响 0 行。 |
| **证据充分性** | 对"注释证据"这一用途：**充分**。三个判据都可在冻结证据上测量，且 PF06850 的参考层对比极强（274/276 vs 0/447）。对"family 判定"这一用途：**不充分**，且候选层没有任何独立标签可用来验证命中集。 |
| **预期收益** | 把 Cys 型候选从"完全没有证据层"变成"有一张可复核、可引用、可证伪的注释证据表"：**29,949 条**同时满足 PF06850 命中且无 lipase box，其中 **25,061 条**进一步携带 VCQ。为下一轮（若有人工实验资源）提供明确的候选优先级清单，且不动摇任何既有结论。 |
| **风险** | 主要风险是**误读**：把候选层注释证据当成 family 判定或表型证据。缓解措施：每一行都带 `evidence_layer / family_call_made=false / hmm_score_present=false / subtype_call_impact=none / new_family_call_made=false` 与 phenotype boundary 文本；manifest 显式声明无 HMM、无 gate 变动；status 文档单列一节说明二者区别。 |
| **可证伪条件** | ① PF06850 在本参考集外被证明不是 Cys 型标记；② 得到独立非 Cys 阴性集后，判别器在其中的命中率不可接受；③ lipase box 缺失项被证明不可判定（miss 不可归因）。任一成立即回退到"只有 PF06850 原始结合、不发布判别器"。 |

### 方案 B：新增窄范围 Cys 型 profile，同时完整保留历史模型

| 维度 | 评估 |
|---|---|
| **治理可行性** | ❌ **被项目硬规则禁止，即使获得人类授权也不应执行。** 训练一个 fitted profile 的严格输入只能是 `experimental_positive` / `experimental_negative` 记录：276 条 Cys 参考中 275 条是 `annotation_only`，被 `validate_phaded_reference_ledger.py` L172–173 直接拒绝；唯一可用的严格记录是 **1 条**（`CAJ92291.1`）。以 1 条序列拟合 profile 在方法与治理上都不成立。 |
| **证据充分性** | ❌ 不足以支撑。除禁训规则外，还有三重缺口：held-out = 0；family-resolved negative = 0；Cys 型亲核体坐标 `pending_reference_annotation`（无法做残基级校准）。gate 明确不满足。 |
| **预期收益** | 名义上可给出"窄范围、只覆盖 Cys 型"的新 profile；但由于最小输入是 1 条序列，任何打分结果都无法与本轮已发布的 PF06850 证据层在证据强度上区分。收益不能成立。 |
| **风险** | 极高：① 直接违反项目硬规则；② 单序列 profile 会给出近乎任意的高分，制造大量假阳性；③ 会污染 `formal_scan_models.tsv` 的版本权威；④ 一旦发布，无法与历史模型做口径分离。 |
| **可证伪条件** | 若将来台账中出现 ≥3 条 Cys 超家族的 `experimental_positive`（且各自带 assay/result、held-out 与 family-resolved negative 齐备），则本方案重新变成可评估对象。在此之前不应执行。 |

### 方案 C：用版本化新模型替换旧模型

| 维度 | 评估 |
|---|---|
| **治理可行性** | ❌ **双重禁止。** 继承方案 B 的全部禁训违规；此外"替换"意味着改写历史模型/registry 的权威关系，与 AGENTS.md"不得覆盖历史 `results/`"、"发布前核对本地/GitHub/服务器 deploy/run manifest 的权威关系"直接冲突。 |
| **证据充分性** | ❌ 最低。替换需要比"新增"更强的证据（需证明旧模型在新口径下失效），而现有证据连"新增"都支撑不了。 |
| **预期收益** | 无。Cys 型超家族目前没有任何已注册 profile 可被"改进"，替换的对象不存在；真实的替换只会发生在其他 4 个已注册超家族上，而那正是本轮明确不做的部分。 |
| **风险** | 最高：不可逆地破坏历史可比性，并使 109,087 条候选的既有 `subtype_call` 失去基线。 |
| **可证伪条件** | 只有在方案 B 的门槛全部满足、且新模型在**独立**挑战集上相对旧模型有一致优势时，才可讨论"版本化并存"（而不是替换）。 |

---

## 4. 预注册证伪条件（本轮已执行，输入见 run 的 `inputs/falsification_preregistration.tsv`）

下列条件在**运行之前**写入 run 的 `inputs/`，由实现脚本 fail-closed 评估（`evaluate_falsification_gates`）：任一条 `failed` → manifest `status=falsification_triggered`、全部 `cys_type_candidate_evidence` 记为 `withheld_falsification_triggered`、进程返回非 0。

| ID | 预注册条件 | 期望观测 | 违反时的回退 | **本轮实测** |
|---|---|---|---|---|
| F1 | PF06850 不得结合 Cys 超家族以外的任何参考 | 超家族外命中 = 0 | 判定失败，撤回全部候选判别，非 0 退出 | ✅ **passed**，超家族外命中 **0** |
| F2 | PF06850 在 Cys 超家族的参考覆盖率不得低于 0.90 | ≥ 0.90 | 判定失败；PF06850 不是可用标记 | ✅ **passed**，**274/276 = 0.992754** |
| F3 | 已知胞外参考（**284 条 type 1 + 71 条 type 2**）中不得出现未被解释的判别命中 | 两超家族命中 = 0 | 判定失败并回退到方案 A（registry 不变 + 仅 candidate-only 说明） | ✅ **passed**，两超家族命中 **0**（type 1: 0/284，type 2: 0/71） |
| F4 | lipase box 判据的 miss 必须可归因：`miss_is_assignable=true` 且全部带坐标的文献锚点被复现 | 4/4 锚点复现 | 判定失败；"无 lipase box"这一项不可解释，证据层不得发布 | ✅ **passed**，`anchors_with_coordinate=4`、`anchors_reproduced=4`、`miss_is_assignable=true`（引自冻结的 Task 7 summary `ee1ce845…`） |
| F5 | 本任务不得引入任何 fitted HMM/profile | 实现不拟合任何模型、不写模型文件 | 任何模型文件或训练后的 profile 出现即整体作废 | ✅ **passed**，`hmm_fitted=false`，run 内无模型文件 |
| F6 | 候选层命中没有独立验证标签，只能是注释证据、绝不可是 family 判定 | 每行带 boundary 列且 `family_call_made=false` | 任何把候选命中升格为 family 判定或表型的表述必须撤回 | ✅ **passed**（列级强制，见 §6） |

**未预注册、事后发现的描述符（必须如此标注，不得当作 gate）**
在 F1–F6 全部通过、预注册产物落盘之后，本轮又做了一次**额外的探索性 pass**：由 Knoll 2009 的 cysteine-1 规则出发观察候选层 cysteine+1，发现三肽 **`VCQ`**（Val-Cys-Gln）在冻结参考台账上的**样本内**对比远优于单看 cysteine-1：

| 描述符 | Cys 超家族（276） | 其他 7 超家族（447） |
|---|---|---|
| cysteine-1 为 Val（`V-C` 二肽） | 256（92.75%） | 93（**20.81%**） |
| `VCQ` 三肽 | **251（90.94%）** | **1（0.22%）** |

- 唯一命中的非 Cys 参考是 `EDL63834.1`（`intracellular nPHASCL with lipase box` / `DED_hfam_2` / annotation_only），属相邻超家族的预期交叉，可解释。
- 模型酶 `CAJ92291.1` 恰好携带 1 个 `VCQ`。
- **诚实限定：** `VCQ` 没有文献依据（Knoll 2009 只写了 cysteine-1），是**从冻结参考台账本身data-derived出来的**；因此它的敏感度/假阳性率是**样本内**数字，工作区内**不存在 held-out 集**（`in_sample=true`、`held_out_set=none`）。它被明确标为 `exploratory_in_sample_descriptor_discovery`、`load_bearing=false`，**不参与任何判别、不作为 gate、不进入 registry 论证**。任何把 `VCQ` 用作过滤器的做法必须先在一个独立的非 Cys 阴性集上测量。
- 由此新增的下一轮可证伪条件 **F7（本轮未预注册，登记为下一轮条件）：** 若在一个独立的非 Cys α/β-水解酶阴性集上 `VCQ` 命中率 > 5%，则该描述符必须降级为纯描述列，不得出现在任何候选优先级论证中。

---

## 5. 推荐

### 5.1 明确推荐

**推荐方案 A**：registry 与 `pipeline/config/` 保持完全不变，**不训练任何新 HMM、不新增任何 profile**，改为发布一层**候选层注释证据**（已实现，见 §6）。

补充理由（按重要性）：

1. **方案 B/C 违反硬规则，不是"暂缓"而是"不应执行"。** 275/276 的 Cys 参考是 `annotation_only`，被 `validate_phaded_reference_ledger.py` L172–173 拒绝进严格训练；唯一可训练的 1 条序列不足以拟合 profile，也不满足 gate 的 `positive_count >= 3`。因此即使人类操作者本轮已授权"完成所有任务"，正确的执行方式仍是**不执行训练**并如实说明原因。
2. **任务书的前提经核实需要更正**（§2.1）：Cys 超家族并非"实验阳性为 0"，而是"恰有 1 条"（`CAJ92291.1`）。这条更正值恰好就是任务书点名的模型酶，但它仍然不足以支撑任何 profile，所以结论方向不变、依据更准确。
3. **DEFECT_D 提供了一条真实、可证伪、零新增模型的替代路径**，且它的参考层对比强到足以支撑"注释证据"这一用途：PF06850 在 276 条 Cys 参考中命中 274 条，在其余 447 条参考中命中 0 条。
4. **失败模式是可控的。** 这条路径没有"训练→过拟合→污染 registry"的单向不可逆风险；它只是新增证据列，`subtype_call` 影响 0 行，且 **预注册条件 F1–F6 共 6 条**（其中 4 条数值化、全部 fail-closed 可回退；F5/F6 为结构性声明）**全部通过**。

> **更正痕迹（`20260917_phaded_final_consistency_fix_01`，2026-09-17）：** 本条原文写「4 条预注册 gate 全部 fail-closed 可回退」，缩略了条数。实测 `inputs/falsification_preregistration.tsv` 注册 **F1–F6 共 6 条**；§4 的表已列全 6 条，其中 4 条被数值化评估（`cys_decision_manifest.json` 的 `falsification_gates` 4 个键，`outcome` 全 `passed`），**F5**（不得引入任何 fitted HMM/profile）与 **F6**（候选命中只能是注释证据）为**结构性声明**，分别由 `hmm_fitted=false` 与列级 `family_call_made=false` 强制。另：§6 产物表原写测试文件「18 个测试，先写后实现」，已改为 **15 个先写后实现（有红灯日志）+ 3 个后补写的探索性 pass 测试（红灯未落盘）**，与同表上一行的红灯证据（`Ran 15 tests`）保持一致。
5. **另外三个缺失超家族仍然保持 `reference_only`**，本轮不作任何处理：periplasmic（全库仅 1 成员且 2025 综述承认生理角色未明）、i-nPHASCL（有 lipase box，参考台账仅 29 条）、native-SCL/PhaZ7-like（AHSMG 与 Bacillus 脂酶 abH18.01 共享，单独不充分）。

### 5.2 推荐方案是否需要服务器计算

**不需要。** 本轮实现只做"读冻结证据 + 读冻结 FASTA + 序列字符串扫描"，全部本地完成（两次 pass 合计约 2 分钟）。不需要 profile 比对、不需要 hmmscan、不需要任何服务器算力。

之所以**不做**"对 DED hfam_61–hfam_69 参考序列做比对打分"这类服务器任务，原因不是成本，而是治理与证据：① 用 1 条可训练记录拟合 profile 违反 annotation_only 禁训规则；② 冻结的 38 个 `.aln` 是原始 ClustalW 块，不含 Knoll 2009 声称的人工标注列（Task 2 已记录该前提本地不可验证），因此也无法用来读取权威的催化 Cys 列；③ 候选层与这些比对零重叠，比对列的结论传不到 109,087 条候选上。

本轮确实探查过"用 Cys 型比对的不变 Cys 列反推亲核体"这条纯本地思路，结论是**不可判定**并已如实放弃：9 个 Cys 家族比对中不变 Cys 列不唯一（如 `aln61` 有 237/413/425 三列、`aln65` 有 120/218/411 三列），仅靠保守性无法唯一指认催化 Cys。该失败探查记录在此，不计入任何产物。

### 5.3 服务器负载读数（只读核对，未提交任何任务）

`2026-09-17T01:46:02+08:00`，`<SERVER_USER>@<SERVER_HOST>`：
`load average: 4.50, 4.77, 4.81`（80 逻辑核）；`Mem: 1007 GB total, 956 GB available`；两张 RTX 4090 均 `0 %` / `1 MiB`。
**本次 run 使用的服务器线程数 = 0**，未创建 dated `deploy/` 目录，未执行任何远程源码。

---

## 6. 已实施的产物（方案 A 的候选层证据层）

Run：`runs/20260917_phaded_inphascl_cys_decision_01`（含 `logs/`、`inputs/`、`results/`、`input_contract.json`；GTDB taxonomy/metadata/tree 均写 `pending`，未伪造哈希）

| 产物 | 大小 | SHA-256 | 说明 |
|---|---:|---|---|
| `results/cys_candidate_evidence.tsv` | 42,397,093 | `3cf1edaeb73a3fd00c585105d66b70d237a1707c0e65944a39ad040c76e47c4d` | 109,087 行 × 20 列候选层证据 |
| `results/cys_reference_contrast.tsv` | 909 | `2a14e7b6d0d190feb4786b964730617ad021c4f5e2d89492f3cd14a72bc4121f` | 8 个超家族的参考层对比 |
| `results/cys_decision_manifest.json` | 6,945 | （见该文件 `outputs`） | 预注册 pass 的 manifest |
| `results/cys_typing_descriptor_evidence.tsv` | 34,288,167 | `ad0017c26a56250fce874fd577db90ec6576ca2b41e3b2d795d7d83235bf89ff` | 探索性 `V-C` / `VCQ` 描述符 pass（第二次 pass，只写新文件名） |
| `results/cys_typing_descriptor_discovery.json` | 16,371 | `e019f7e7487d85678d62c80d4e45640a675adc7e8f68575a9a6425981302e9d4` | 样本内对比报告，含 `in_sample=true` / `held_out_set=none` |
| `results/cys_typing_descriptor_manifest.json` | — | — | 记录**预注册产物未被改写**的哈希 |
| `inputs/cys_criterion_definitions.tsv` | 954 | — | 三个判据的可观察定义、文献依据、是否 load-bearing |
| `inputs/falsification_preregistration.tsv` | 1,827 | — | F1–F6 预注册条件 + 回退动作 |
| `inputs/input_manifest.tsv` | 1,068 | — | 每个输入的 path/size/SHA-256；GTDB 三项 `pending` |
| `logs/step1_failing_test_red.log` | 34,006 | — | **先写失败测试**的红灯证据（15 个测试全部因实现文件不存在而 error） |
| `logs/annotation_stdout.log`、`logs/descriptor_pass_stdout.log`、`logs/commands.txt`、`logs/descriptor_pass_commands.txt` | — | — | 两次 pass 的完整 stdout 与执行命令 |
| `pipeline/scripts/annotate_phaded_inphascl_cys_evidence.py` | — | — | 实现（不拟合任何模型） |
| `pipeline/tests/test_annotate_phaded_inphascl_cys_evidence.py` | — | — | 18 个测试：**15 个先写后实现**（红灯日志见上一行），另 **3 个为第二次探索性 descriptor pass 后补写**（红灯未落盘，证据缺口见状态文档 §11） |

**候选层结果（109,087 条全覆盖，无遗漏行）**

| 判据 / 组合 | 计数 |
|---|---:|
| PF06850 `detected` | 30,621 |
| PF06850 `not_detected_in_tested_pfam` | 67,798 |
| PF06850 `not_tested_no_pfam_evidence` | 10,668 |
| lipase box `supported` | 77,302 |
| lipase box `not_detected_pattern` | 31,785 |
| **`PF06850 detected` ∧ `无 lipase box`（声明的判别证据）** | **29,949** |
| └ 其中无超家族先验（`unassigned_PhaDED_like`） | 29,271 |
| └ 其中先验为 `intracellular nPHAMCL` | 654 |
| └ 其中先验为 `extracellular dPHASCL type 1` | 24 |
| 上述再叠加 `V-C` | 26,023 |
| 上述再叠加 `VCQ`（探索性） | 25,061 |

**强制边界列（109,087 行每行都带）**
`evidence_layer=candidate_layer_annotation_evidence`、`family_call_made=false`、`hmm_score_present=false`、`subtype_call_impact=none`、`new_family_call_made=false`、`evidence_boundary="candidate_only; candidate-layer annotation evidence of candidate homology or functional potential and not a validated PHB/PHA degradation phenotype"`。

**候选层与 fitted HMM 证据的区别（必须分清）**

| | 本轮的候选层注释证据 | 一个 fitted HMM/profile 会给出的东西 |
|---|---|---|
| 证据类型 | 域注释（Pfam 绑定）+ 成熟序列模式（GxSxG 及其缺失）+ 序列上下文（cysteine±1） | 概率模型的对数似然打分与 E-value |
| 是否有模型 | **无**（`hmm_fitted=false`，run 内无模型文件） | 有，且需版本化绑定 |
| 能否进 registry | **不能**，也不试图进 | 可被考虑，但必须先通过 gate |
| 参考层对比强度 | PF06850 274/276 vs 0/447（强） | 由 held-out 集决定（当前为 0） |
| 已校准的阴性/挑战面 | **无**（family-resolved negative = 0） | 同样为 0，故仍不可解锁 |
| 可否作为表型证据 | **不可以** | 同样不可以（candidate-only） |
| 可否删除候选 | **不可以** | 不可以 |

---

## 7. 未解决项与交接

1. **候选层判别集没有任何独立验证标签**（109,087 条中 0 条属 DED 比对、0 条属参考台账）。因此 29,949 / 25,061 这两个数字**只能**作为候选优先级，不能被写成敏感度/特异度。
2. **Cys 型催化亲核体坐标仍然 `pending_reference_annotation`。** 本轮的 cysteine±1 只是**上下文描述符**，不是残基级判定；候选层的 "expected catalytic residue" 依然全是 `Ser`。
3. **`VCQ` 是样本内发现**，需要独立非 Cys 阴性集才能评估（F7）。在 F7 通过之前不得用作过滤器。
4. **gate 未变、也未达成。** 本任务不降低校准门槛，37 个 profile 的阴性/held-out 仍为 0。
5. **上游 gate 计数器口径**（Phase 2 报告 Task 5 指出）只统计绑定到 profile 自身 family 的记录，跨家族轴对照贡献恒为 0；该缺陷不在本任务范围，仅报告不修。
6. **其他三个缺失超家族**（periplasmic、i-nPHASCL with lipase box、native-SCL/PhaZ7-like）继续保持 `reference_only`。

---

## 8. 引用

- Knoll M, Hamm TM, Wagner F, Martinez V, Pleiss J. **The PHA Depolymerase Engineering Database.** BMC Bioinformatics 2009;10:89. [PMC2666664](https://europepmc.org/article/PMC/PMC2666664)
  - 无 lipase box 家族使用催化 Cys："A few intracellular nPHASCL depolymerases have no lipase box, but have a catalytic triad consisting of cysteine, histidine, and aspartic acid."
  - cysteine-1 规则："all family members of the family of intracellular nPHASCL depolymerases (no lipase box) also have a hydrophobic residue (almost all valine) at position cysteine-1."
- 项目内部冻结证据：`runs/20260915_phaded_motif_reconciliation_01/results/reference_motif_panel.tsv`、`runs/20260917_phaded_motif_completion_full_01/results/motif_completion_summary.json`、`runs/20260917_phaded_reference_residue_mapping_01/results/reference_mapping_manifest.tsv`、`runs/20260911_phaded_pfam_architecture_01/inputs/phaded_reference_ledger.tsv`、`runs/20260917_phaded_phase2_verification_01/results/phase2_verification_report.md`
- 相关状态文档：`docs/T141_20260917_phaded_inphascl_cys_decision_status.md`

---

## 9. 2026-09-17 人类授权修订（`runs/20260917_phaded_cys_discovery_hmm_01`）

> **本节为纯追加内容。** 上面 §1–§8 的原文一字未改（包括 §5.1「推荐方案 A / 不训练任何新 HMM」与 §5.2「不需要服务器计算」的结论）；本节记录人类操作者随后给出的授权如何**在其自身范围内**改变执行范围。两处结论的适用口径不同，**不构成互相推翻**：§5 讨论的是**校准层 / registry 层**的拟合 profile，本节授权的是**发现层、未校准**的 HMM。

**修订时间：** `2026-09-17T14:04:52+08:00`（= `2026-09-17T06:04:52Z`），本地工作区 `<REPO_ROOT>`，分支 `main`，`HEAD = 2e44aafe2c91263447f37bac625f97c56ddafc05`（未提交、未推送）。
**授权 run：** `runs/20260917_phaded_cys_discovery_hmm_01`（本地与服务器同名；本轮新建）。
**服务器：** `<SERVER_USER>@<SERVER_HOST>`（grius），`${PHB_REMOTE_ROOT}/PHB_gtdb-ds`。本轮**使用**服务器 CPU（HMMER/MAFFT，单任务线程 ≤ 40、无 GPU 需求），这与 §5.2 的"不需要服务器计算"是**不同轮次的不同结论**。

### 9.1 人类授权的内容（逐条）

1. **授权训练 HMM。** 人类操作者于 2026-09-17 明确授权"可以训练 hmm"，并要求执行本轮 dated run。
2. **范围选择：** (a) 目标家族 = **Cys 型家族**（`intracellular nPHASCL without lipase box`，冻结台账 276 条参考）；(b) **对现有 109,087 条候选池打分**（`runs/20260909_phaded_candidate_mapping_01/inputs/candidate_union.faa`）。
3. **训练数据可以包含 DED 家族序列（含 `annotation_only`）。** 依据：**DED 自身就是这么做的** —— Knoll 2009 明说 PHA Depolymerase Engineering Database 为每个 family 提供 **HMMER profile** 用于 in silico 鉴定。因此把 DED 家族成员用作比对/建 profile 的输入，是**复现文献方法**，而不是把 `annotation_only` 记录升格为实验证据。
4. **新模型层命名：`discovery_hmm_uncalibrated`**（发现层 HMM，未校准）。该名字是本轮新引入的层级标签，必须在 HMM 文件、打分表、claims 表、闸门结果与状态文档中逐处标注。

### 9.2 授权**不**改变的东西（治理红线，逐条保留）

| 红线 | 本轮事实 |
|---|---|
| **训练 ≠ 校准** | HMM 是 `discovery_hmm_uncalibrated`；未做任何阈值校准、未做 held-out 校准、未声明敏感度/特异度。 |
| **训练 ≠ 进 registry** | `pipeline/config/formal_scan_models.tsv` **绝不修改**（前后 SHA-256 逐字节比对）；不新增 registry 条目；本轮模型文件只存在于本轮 dated run。 |
| **命中仍是 candidate-only** | 命中只表示候选同源或功能潜力，**不产生 family 判定**（`family_call_made=false`、`new_family_call_made=false`）；**对 `subtype_call` 影响 0 行**。 |
| **`annotation_only` 禁训规则仍适用于严格/校准层** | `validate_phaded_reference_ledger.py` L172–173 继续有效；本轮**不申请**严格训练资格，只是以 DED 家族成员复现文献的 profile 构建方式。 |
| **校准 gate 原样保留** | 每家族 ≥3 独立阳性 + ≥1 held-out + ≥1 family-resolved 阴性 + ≥1 challenge + 零未解释命中，全部**未触碰、仍未达成**；§5.1 关于"gate 不满足"的判断继续有效。 |
| **历史 run 只读** | `runs/20260901_*`–`runs/20260917_*` 既有 run 全部只读，不覆盖、不修改、不删除。 |

### 9.3 §5.1 / §5.2 与本节的适用口径对照（避免误读）

| §5 的结论 | 适用对象 | 本节是否推翻 |
|---|---|---|
| §5.1「不应执行训练」 | **校准层 / registry 层** fitted profile（需通过 gate 才能进 registry） | **不推翻** —— 该类 profile 本轮**仍然不建**，gate 仍然不满足 |
| §5.2「不需要服务器计算」 | **方案 A**（仅本地候选层注释证据） | **不推翻** —— 方案 A 自身确实不需要；本节是**新增的一条并行发现层路线**，不是对方案 A 的否定 |
| 「下一会话不应再尝试训练该 HMM」 | 同上（特指"用 1 条可训练记录拟合 registry 级 profile"） | **在该措辞的对象上仍然成立**；但**不适用于**本节授权的 `discovery_hmm_uncalibrated` 发现层（其设计输入为 270 条 DED 家族成员、共 1 张 MSA，且被明确标注为未校准、不进 registry） |

### 9.4 本轮同时登记的证据更正（与授权相关的关键事实）

见 `runs/20260917_phaded_cys_discovery_hmm_01/results/reference_evidence_audit.tsv`：任务书点名的 **`O87189` 与冻结台账 `DED_hfam_65_0001`（`CAJ92291.1`）序列逐字节相同**（同 419 aa），`Q7WT48`/`Q7WT49` 分别与 `AAP74581.1`/`AAP74580.1` 逐字节相同 —— 即 Task 5 面板中被记为"胞内定位阴性"的三条记录**其实都已经是冻结台账成员**，其中 `O87189` 就是台账**唯一**那条 `experimental_positive` 的别名。该事实**不改动冻结台账**（台账按规则只读），只在本轮 run 内记录，并直接影响本轮训练集与闸门的边界声明。

---

*本 decision record 为 candidate-only 记录；所有 profile、domain、motif 与序列上下文证据只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
