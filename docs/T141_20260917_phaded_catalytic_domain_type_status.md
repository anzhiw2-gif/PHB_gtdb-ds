# T141 PhaDED 催化域 type 1/type 2 归属量化 + AHSMG 事实固化（计划 v2 · Task 4）

日期：2026-09-17
Run：`runs/20260917_phaded_catalytic_domain_type_01`
分支/HEAD：`main` / `2e44aaf`（本 Task 未触碰 git 历史）
计划依据：`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` 第 3 节 Task 4
判据文献：Knoll M, et al. **The PHA Depolymerase Engineering Database.** BMC Bioinformatics 2009;10:89（PMC2666664）

边界（不可突破）：本轮为**本地、只读、candidate-only 纯分析**。没有 ssh、没有服务器任务、没有启动任何 HMM/结构/系统发育重算，没有修改 `pipeline/config/formal_scan_models.tsv`、历史 HMM、registry、历史 run 或任何既有 `results/`。催化域 type 归属是**候选同源/功能潜力的序列几何证据，不等同于已验证的 PHB/PHA 降解表型**。

---

## 1. 一句话结论

项目当前以超家族 profile 命中把 **67,342** 条候选标为 `dPHASCL1_like_candidate`；其中真正通过 **type 1 几何验证**（oxyanion hole 严格位于 lipase box N 端）的只有 **16,511 条 = 24.5181%**，另外 **50,458 条（74.93%）连 oxyanion hole 坐标证据都没有**，**373 条（0.55%）相对位置与标签冲突**。因此 `dPHASCL1_like_candidate` 是一个**标签**，`type1_verified` 是一次**几何验证**，二者绝不可互换表述。

---

## 2. 判据与比较口径

### 2.1 文献原文（本仓库内可复核）

`research/europepmc/comp_ft_knoll_phaded.txt` L17 逐字记录 Knoll 2009 原文：

> "Within the sequences of **type 1** catalytic domains, the oxyanion hole can be found **N-terminal to the lipase box**, similar to lipases. Within the sequences of **type 2** catalytic domains, the oxyanion hole is found **C-terminal to the catalytic triad**."

`research/europepmc/comp_classification.md` L60–61 同口径复述（type 1 = oxyanion hole 位于 lipase box N 端；type 2 = 位于催化三元组 C 端）。

### 2.2 实现口径（复用既有代码，不另立第二套口径）

`pipeline/scripts/audit_phaded_motifs.py` L177–180：

```python
if lipase and oxy and "extracellular dPHASCL type 1" in subtype and oxy.start() > lipase.start():
    oxy_state = "conflict_relative_position"
if lipase and oxy and "extracellular dPHASCL type 2" in subtype and oxy.start() < lipase.start():
    oxy_state = "conflict_relative_position"
```

即：**type 1 的 oxyanion hole 起点必须小于 lipase box 起点（N 端）**；起点大于则被判为相对位置冲突。本 run 完全复用该口径，并以测试 `test_convention_matches_audit_phaded_motifs_conflict_semantics` 把两个模块钉在同一几何上（同一序列在 audit 模块得到 `supported` / `conflict_relative_position`，在本模块必须得到 `type1_verified` / `undetermined_position_conflict`）。

### 2.3 与任务书 Step 4.1 一句表述的差异（必须记录，不得静默）

任务书 Step 4.1 写"type 1 要求 oxyanion 位置**严格大于** lipase 位置"。该表述与**三处独立证据**相反：

1. Knoll 2009 原文（§2.1）：type 1 = oxyanion hole 位于 lipase box **N 端**，即坐标更小；
2. `audit_phaded_motifs.py` L177–180（任务书同时要求"复用其口径"）：type 1 分支把 `oxy.start() > lipase.start()` 记为 **conflict**，说明 type 1 期望的是 `oxy.start() < lipase.start()`；
3. 候选实测：`motif_reference_subtype` 为 type 2 的 2,362 条记录**全部**是 `oxyanion 起点 > lipase 起点`（2,362/2,362 = 100%），若采用"type 1 = 大于"的读法，则这 2,362 条会被全部判成冲突、而 100% 的 type 2 证据被清零。

因此本 run 采用**文献与既有实现一致的口径**（type 1 = oxyanion 严格 N 端 / 坐标更小），并把被否证的读法量化为翻转灵敏度（见 §3.4），以公开方式证明该读法不自洽，而不是静默二选一。

### 2.4 fail-closed 规则（Step 4.1 的核心不变量，已全部测试）

| 条件 | 判定 | basis |
|---|---|---|
| 缺 `oxyanion_hole_coordinates` 或 `lipase_box_coordinates`（空/缺失） | `undetermined_no_oxyanion` | `oxyanion_coordinates_missing` |
| 坐标不可解析（非 `start-end`、含 0/负、end<start） | `undetermined_no_oxyanion` | `invalid_coordinate_format` |
| `oxyanion_hole_state = not_detected_pattern` | `undetermined_no_oxyanion`（任何情况下都不得获得 type） | `oxyanion_signature_not_supported` |
| `oxyanion_hole_state = conflict_relative_position` | `undetermined_position_conflict`（不得升级为已验证） | `oxyanion_relative_position_conflict` |
| 两者起点相等 | `undetermined_position_conflict` | `coordinate_start_positions_equal` |
| oxyanion 起点 < lipase 起点（state=supported） | `type1_verified` | `oxyanion_n_terminal_of_lipase_box` |
| oxyanion 起点 > lipase 起点（state=supported） | `type2_verified` | `oxyanion_c_terminal_of_lipase_box` |

`catalytic_domain_type` 取值被**封闭**在 4 个值内：`type1_verified` / `type2_verified` / `undetermined_no_oxyanion` / `undetermined_position_conflict`。

**标签不参与判定**：`subtype_call` 只在判定完成后用于交叉列表（测试 `test_label_is_never_used_to_assign_a_type` 证明同几何在任一标签下结果相同），从而避免 Task 4 自身引入循环论证。

---

## 3. Step 4.4 量化结果（全部实测，可复核）

### 3.1 全库 type 归属

输入 `runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv`（109,087 行，SHA-256 `796cb0a3…569a`）。

| `catalytic_domain_type` | 数量 | 占 109,087 | 含义 |
|---|---:|---:|---|
| `type1_verified` | **16,530** | 15.1530% | oxyanion hole 严格位于 lipase box N 端 |
| `type2_verified` | **2,362** | 2.1652% | oxyanion hole 严格位于 lipase box C 端 |
| `undetermined_no_oxyanion` | **89,822** | 82.3398% | 无可用 oxyanion 坐标证据（**不是阴性**） |
| `undetermined_position_conflict` | **373** | 0.3419% | 上游已标记相对位置冲突 |
| 合计 | 109,087 | 100% | 自洽（`totals_self_consistent = true`） |

对上游客观分布的核对（完全对齐，无一处改写）：

- `oxyanion_hole_state`：`not_detected_pattern` 89,812 + `supported` 18,902 + `conflict_relative_position` 373 = 109,087；
- 89,822 = 89,812（`not_detected_pattern`）+ **10**（state=`supported` 但 `lipase_box_coordinates` 为空 → fail-closed，10 条均为 `unassigned_PhaDED_like`）；
- 18,892 = 16,530 + 2,362 = 同时具备两个坐标字段且 state=`supported` 的行数 → **有能力做几何判定的仅占全库 17.3181%**。

### 3.2 关键量化：67,342 条 `dPHASCL1_like_candidate` 的 type 1 归属

| 分组 | 计数 | 占该标签 |
|---|---:|---:|
| `dPHASCL1_like_candidate` 标签总数 | **67,342** | 100% |
| 其中 `type1_verified` | **16,511** | **24.5181%** |
| 其中 `type2_verified` | 0 | 0% |
| 其中 `undetermined_no_oxyanion` | 50,458 | 74.9295% |
| 其中 `undetermined_position_conflict` | 373 | 0.5539% |
| `dPHASCL2_like_candidate` 标签总数 | **356** | 100% |
| 其中 `type2_verified` | **323** | **90.7303%** |
| 其中 `undetermined_no_oxyanion` | 33 | 9.2697% |
| 其中 `undetermined_position_conflict` | 0 | 0% |

即：**67,342 → 16,511（24.52%）**，与预期方向一致（远小于 67,342，因为全库只有 18,902 条有 oxyanion 模式命中）。**75.48% 的 `dPHASCL1_like_candidate` 目前没有任何几何支持**，其中绝大部分（74.93%）是"证据缺失"而非"几何不符"。

额外的诚实观察：**373 条相对位置冲突全部落在 `dPHASCL1_like_candidate` 标签内**（`label_contradicts_geometry = 373`），`dPHASCL2_like_candidate` 为 0 条。

### 3.3 标签 × 几何交叉表（逐行实测）

| `subtype_call` | 标签数 | `type1_verified` | `type2_verified` | `undetermined_no_oxyanion` | `undetermined_position_conflict` |
|---|---:|---:|---:|---:|---:|
| `dPHASCL1_like_candidate` | 67,342 | 16,511 | 0 | 50,458 | 373 |
| `dPHASCL2_like_candidate` | 356 | 0 | 323 | 33 | 0 |
| `ambiguous_family_within_superfamily` | 9,631 | 19 | **2,020** | 7,592 | 0 |
| `ambiguous_superfamily` | 41 | 0 | 0 | 41 | 0 |
| `hold_architecture_conflict` | 83 | 0 | 19 | 64 | 0 |
| `hold_gene_model_or_structure` | 2 | 0 | 0 | 2 | 0 |
| `unassigned_PhaDED_like` | 31,632 | 0 | 0 | 31,632 | 0 |

两个副产品发现（本轮新得到，均为候选层几何证据，非表型结论）：

1. 在 family 未解析的 `ambiguous_family_within_superfamily`（9,631 条）中，几何判定解出了 **2,039 条（21.17%）**：19 条 type 1 + **2,020 条 type 2**。几何证据在 family 层无法解析处提供了独立于 profile 的判别信息；
2. 在 `hold_architecture_conflict`（83 条）中，19 条具备 type 2 几何（保持 hold，不因本轮的 type 判定而解除 hold）。

跨标签来源的一致性核对：0916 三份全库 subtype matrix（signalp / evidence-amendment / foldseek-amendment）的 `subtype_call` 在 **109,087 行上逐行完全一致**。该核对在运行内执行并**落盘留证**（`run_manifest.json → label_source_consistency`：`checked_matrices` 3 份、`checked_accessions` 109,087、`accession_sets_identical=true`、`subtype_call_identical=true`、`disagreeing_accession_count=0`）；任一行分歧即 fail-closed 拒绝输出（测试 `test_disagreeing_matrices_are_rejected`、`test_accession_set_mismatch_is_rejected`）。因此 67,342 这一标签数不是单一文件的偶然产物。

### 3.4 被否证读法的翻转灵敏度（公开量化，仅用于证伪，不作结果）

若把 §2.3 那句表述当作口径（type 1 = oxyanion 坐标更大），同一份数据会变成：`type1_verified` 2,362 / `type2_verified` 16,530 / 冲突 373 —— 即 100% 的 type 2 证据反转、`dPHASCL2_like_candidate` 的 323 条 type 2 验证全部失效。该读法与 Knoll 2009 原文和 `audit_phaded_motifs.py` 实现同时冲突，故不采用。数值保留在 `results/type_summary.json → orientation_sensitivity_if_inverted`，供独立复核者自行判定。

---

## 4. Step 4.5 AHSMG 事实固化

### 4.1 参考面板阳性对照（已通过，可复核）

`runs/20260915_phaded_motif_reconciliation_01/results/reference_motif_panel.tsv`（723 条参考，SHA-256 `230ca43f…91a0`）：

| 维度 | 实测 |
|---|---|
| `ahsmg_state = supported` | **2** / 723 |
| `not_detected_pattern` | 721 / 723 |

| 参考 accession | AHSMG 坐标 | PhaDED family | 参考亚型 | 面板 `evidence_status` |
|---|---|---|---|---|
| `AAK07742.1` | 172-176 | `DED_hfam_7` | `extracellular native-SCL/PhaZ7-like` | `experimental_positive` |
| `2VTVA` | 134-138 | `DED_hfam_7` | `extracellular native-SCL/PhaZ7-like` | `annotation_only`（PDB 结构条目） |

`2VTVA` 是 PhaZ7 的 PDB 结构 **2VTV 的链 A**（计划 v2 §1.3 记录；本条为结构注释参考，其面板证据等级为 `annotation_only`，而 `AAK07742.1` 为 `experimental_positive`——两者共同构成阳性对照，其中 1 条带实验阳性证据等级）。

**本轮独立复现**：脚本用与 `audit_phaded_motifs.AHSMG_PATTERN` 完全相同的模式 `AHSMG`，在参考序列集 `runs/20260911_phaded_pfam_architecture_01/inputs/phaded_reference.faa`（723 条，SHA-256 `00403a57…ed5f`）上重跑，得到**恰好同 2 条命中、坐标完全相同**：`AAK07742.1` 172-176（序列长 380）、`2VTVA` 134-138（序列长 342），其余 721 条无命中。因此"阳性对照通过"不是转述结论，而是本轮可复现的观测。

### 4.2 候选层零检出

候选全库 109,087 条 `ahsmg_state` **全部为 `not_detected_pattern`（109,087 / 109,087，命中 0）**。

### 4.3 结论表述（唯一允许的写法）

> **以已通过阳性对照的 AHSMG 模式，在 109,087 条候选池中未检出 PhaZ7 型（native-SCL/PhaZ7-like）基序。**

**绝对不得**写成（违反即为口径错误）：

- "PhaZ7 型家族不存在" / "该族不存在"；
- "自然界不存在该酶" / "该酶不存在"；
- "候选全库阴性"（零检出是**未检出**，不是生物学阴性）。

**候选池来源可能不含该族的说明（必须与上句同时出现）**：本项目 109,087 条候选宇宙来自 GTDB 全库 profile HMM 命中集合（见 `runs/20260911_phaded_pfam_architecture_01/inputs/phaded_candidate_union.faa` 的构建链路）。这个集合**可能根本不含 PhaZ7 型家族成员**，因此零结果首先说明的是**候选池的覆盖范围**，而不是该基序或该酶在自然界的存在性。该口径已作为 `AHSMG_CONCLUSION_BOUNDARY` / `CANDIDATE_POOL_SOURCE_NOTE` 写入脚本常量，并随 `results/type_summary.json` 与 `run_manifest.json` 一并落盘；测试 `test_absence_statement_boundary_forbids_family_absence_claims` 断言该文本包含 scoping 且**不含**上述禁止表述。

---

## 5. Step 4.6 标签与验证严格分开（可发表口径模板）

| 术语 | 是什么 | 不是什么 |
|---|---|---|
| `dPHASCL1_like_candidate`（67,342） | superfamily **profile 命中标签**（`build_phaded_subtype_matrix.py` L153–156 由 `phaded_superfamily_best` 决定） | **不是** type 1 验证，**不是**表型确认 |
| `dPHASCL2_like_candidate`（356） | superfamily **profile 命中标签**（同一函数 L157–158） | 同上；其 type 2 归属同样需要几何验证 |
| `type1_verified`（16,530） | 由 oxyanion hole 与 lipase box 坐标**严格相对位置**得出的几何验证（本轮产出） | **不是**家族归属、**不是**完整催化三联体、**不是**降解活性 |
| `type2_verified`（2,362） | 同上，方向相反 | 同上 |

必须分开书写的理由（本轮实测支撑）：标签 67,342 条与几何验证 16,511 条相差 **50,831 条**（4.08 倍）；若把两者当同一件事，会把 74.93% 无 oxyanion 证据的候选误表述为"已确认 type 1"。同理，`dPHASCL2_like_candidate` 356 条与 `type2_verified` 323 条也不是同一个集合，且 2,362 条 `type2_verified` 中只有 323 条带该标签——**反向也不成立**。

---

## 6. 已知限制（不得回避）

1. **模式启发式**：oxyanion 判定来自 `H[A-Z]?G[A-Z]?C[A-Z]?Q`（HGCXQ-like 模式）、lipase box 来自 `GxSxG` 的**首个**匹配。坐标是模式命中位置，**不是**结构上的 oxyanion hole 残基；本轮的"几何"是**序列级相对位置**，与 Knoll 2009 的序列级定义一致，但不等于三维结构验证。
2. **覆盖率有限**：仅 17.3181% 的候选（18,892 / 109,087）可做几何判定。82.34% 的 `undetermined_no_oxyanion` **不构成 type 1/type 2 的否定证据**，任何"多数候选并非 type 1"的表述都不被本轮数据支持。
3. **面板不完整**：本轮只消费 oxyanion + lipase box 两个 motif 字段。His 模式 723/723 从未触发、Asp 仅对 type 1 生效、lid/linker/SBD 未绑定（计划 v2 §1.3）——**催化三联体完整性不在本 Task 的可辩护范围内**（属 Task 2）。因此 `type1_verified` 只表示"oxyanion 相对位置的序列几何与 type 1 定义一致"。
4. **输入未绑定 GTDB taxonomy/metadata/tree**：本分析不消费这三项，故 `input_contract.json` 的 `status` 如实为 `pending`，而 6 个实际输入全部为 `verified`（含 SHA-256）。
5. **不解除任何 hold**：`hold_architecture_conflict` 等既有 hold 与 ambiguity 分类一律保留，type 判定不得作为解除依据。

---

## 7. 未做 / 未修改（授权记录）

- `server_execution_started: false`、`ssh_started: false`、`formal_scan_started: false`；
- `formal_registry_modified: false`、`formal_scan_models_modified: false`、`historical_run_modified: false`；
- 未删除任何文件；未执行 git commit / push / clean / reset；未修改其他 agent 的产物；
- 本轮**只读**消费 `runs/20260915_*`、`runs/20260916_*`、`runs/20260911_*` 中的既有证据表，未回写任何历史 run。

## 8. 主要产物

- `pipeline/scripts/assign_phaded_catalytic_domain_type.py`（新；SHA-256 `154cf45e…dfa9`）
- `pipeline/tests/test_assign_phaded_catalytic_domain_type.py`（新；24 个测试；SHA-256 `631db33f…02db8`）
- `runs/20260917_phaded_catalytic_domain_type_01/results/catalytic_domain_type.tsv`（109,087 行；SHA-256 `cb97645b…ec77`）
- `runs/20260917_phaded_catalytic_domain_type_01/results/type_summary.json`（SHA-256 `09ed6a5d…b88c`）
- `runs/20260917_phaded_catalytic_domain_type_01/inputs/source_manifest.tsv`（6 个输入的 path/size/SHA-256；SHA-256 `b175ba41…030f`）
- `runs/20260917_phaded_catalytic_domain_type_01/logs/catalytic_domain_type.log`（SHA-256 `33954c37…c268e`）
- `runs/20260917_phaded_catalytic_domain_type_01/run_manifest.json`（SHA-256 `1d2ff9de…1e12`）
- `runs/20260917_phaded_catalytic_domain_type_01/input_contract.json`（SHA-256 `dc399dc9…9321`）

输入 SHA-256（与 `inputs/source_manifest.tsv` 一致）：

| 输入 | size | SHA-256 |
|---|---:|---|
| `motif_candidate_evidence.tsv` | 64,498,737 | `796cb0a32cc7fde1c768b2c0f2f8e30542c7aa37e773f2b8004a482ac878569a` |
| `phaded_full_library_localization_subtype_matrix.tsv`（signalp） | 156,562,006 | `6f7f74b12a71f2fc5d2109b638d474bec837110087189ec259bf4e1c411e6517` |
| `phaded_full_library_subtype_matrix.tsv`（evidence amendment） | 156,558,847 | `15d8eb37ce9fef81de724b063e90d952b5a16efd02c279bb1d7b858f2cc10318` |
| `phaded_full_library_foldseek_subtype_matrix.tsv` | 171,507,779 | `691805823029b03225858d1632ed0f11204b9350db326353efeb1ae57627868b` |
| `reference_motif_panel.tsv` | 394,009 | `230ca43fb2ad069b6a0a43acae23688116838c8dd5ba4f0951fd3d979f7a91a0` |
| `phaded_reference.faa` | 300,777 | `00403a576ae492e7f42113cd9fc895e0d8ae6cdac48a7cdbaaeca06ad75fed5f` |

## 9. 复现命令

```powershell
python pipeline/scripts/assign_phaded_catalytic_domain_type.py `
  --motif-evidence runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv `
  --subtype-matrix runs/20260916_phaded_full_library_signalp_01/results/phaded_full_library_localization_subtype_matrix.tsv `
  --subtype-matrix runs/20260916_phaded_full_library_evidence_amendment_01/results/phaded_full_library_subtype_matrix.tsv `
  --subtype-matrix runs/20260916_phaded_full_library_foldseek_amendment_01/results/phaded_full_library_foldseek_subtype_matrix.tsv `
  --reference-panel runs/20260915_phaded_motif_reconciliation_01/results/reference_motif_panel.tsv `
  --reference-fasta runs/20260911_phaded_pfam_architecture_01/inputs/phaded_reference.faa `
  --run-dir runs/20260917_phaded_catalytic_domain_type_01

python -m unittest discover -s pipeline/tests
python -m compileall -q pipeline
git diff --check
```

`--expected-total-candidates`（默认 109,087）与 `--expected-dphascl1-labeled`（默认 67,342）在运行中做 fail-closed 校验：候选宇宙规模或标签数与预期不符时**拒绝写出结果**，避免在未预期的候选宇宙上静默产出。

## 10. 验证记录（本次实测）

| 检查 | 命令 | 结果 |
|---|---|---|
| 本 Task 聚焦测试 | `python -m unittest discover -s pipeline/tests -p "test_assign_phaded_catalytic_domain_type.py"` | **24 tests OK** |
| 全量测试套件 | `python -m unittest discover -s pipeline/tests` | 354 tests，**3 failures + 3 errors（全部来自 `test_complete_phaded_reference_residue_mapping.py`，见下）**，skipped=1 |
| 全量测试套件（排除并行 agent 红阶段文件后） | 见下说明 | **349 tests，0 failures，0 errors，skipped=1 → OK** |
| 语法编译 | `python -m compileall -q pipeline` | exit 0 |
| 空白/补丁检查 | `git diff --check` | exit 0 |
| 工作区改动审计 | `git status --porcelain` | 仅新增未跟踪文件；唯一 `M` 项为 `AGENTS.md`（**进入本 Task 前即已存在**，本轮未触碰） |

关于上述 6 个失败/错误的归属说明（只报告、不修改）：它们全部来自 `pipeline/tests/test_complete_phaded_reference_residue_mapping.py`（mtime `2026-09-17 01:05`），报错内容统一为 `required file is missing: pipeline/scripts/complete_phaded_reference_residue_mapping.py`，即该文件是**并行执行的 Task 2 agent 的失败测试先行（红阶段）产物**，其实现脚本尚未落地。该文件不属于本 Task 产物，本轮一律未修改、未删除。

历史输入只读核对：`motif_candidate_evidence.tsv` 的 `LastWriteTime` 仍为 `2026-09-15 02:25:52`、size 仍为 64,498,737，未被本轮改写。

---

*本文件为 candidate-only 研究产物；所有 type 判定表示候选同源或功能潜力，不等同于已验证 PHB/PHA 降解表型。*
