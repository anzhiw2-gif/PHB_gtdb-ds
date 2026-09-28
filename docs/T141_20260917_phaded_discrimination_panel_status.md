# T141 按判别轴重设计 challenge panel（计划 v2 · Task 5）

日期：2026-09-17
Run：`runs/20260917_phaded_discrimination_panel_01`
分支/HEAD：`main` / `2e44aaf`（本 Task 未 commit、未 push、未触碰 git 历史）
计划依据：`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` 第 3 节 Task 5（Step 5.1–5.3）
判据文献：Knoll M, et al. **The PHA Depolymerase Engineering Database.** BMC Bioinformatics 2009;10:89（PMC2666664）；lid 依据 Appl Microbiol Biotechnol 2025（PMC12504323）

边界（不可突破）：本轮为**本地、只读消费既有证据**的设计阶段。没有 ssh、没有服务器任务、没有新增任何网络检索、没有训练 HMM、没有全库重扫、没有修改 `pipeline/config/formal_scan_models.tsv`、历史 HMM、registry、历史 run 或任何既有 `results/`。所有产出为 candidate-only：HMM/profile/domain/motif/结构证据只表示**候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型**。

---

## 1. 一句话结论

本轮为 4 条判别轴建立了**可审计、fail-closed、绑定唯一的轴锚点面板**，共 **960 条**记录（正端锚点 291、负端显式实验阴性 4、challenge 665）；但**不降低门槛**，也不产生任何新 family 判定：**4 条轴上的 family-resolved 阴性（可计入 gate 的阴性槽位）仍然是 0**，37 个 reference-only profile 仍然**全部阻断**（`profiles_meeting_full_gate = 0`）。本轮的实质进展是：**阴性对照从此有明确的锚定依据、证据来源与"为什么落在该轴该端"的可复核理由**，同时把"缺什么才能补齐"量化到单条记录级别。

---

## 2. 交付物与哈希

| 产物 | 路径 | size | SHA-256 |
|---|---|---:|---|
| 轴锚定 challenge panel | `runs/20260917_phaded_discrimination_panel_01/results/challenge_panel_by_axis.tsv` | 860,439 | `0431a9edf17685b353ac835e45469bed9f5f3bba995c47fd895332b7237d0755` |
| 逐轴可达性 + 逐 profile 缺口 | `runs/20260917_phaded_discrimination_panel_01/results/axis_reachability.json` | 176,725 | `555c68e5b1793bbbf83dadbfdc425751fa0b9a1974a42682fdff0110e05eb016` |
| 设计脚本 | `pipeline/scripts/design_phaded_discrimination_panel.py` | 61,558 | `e7a47de93d07c9ea407ac5761f8c8a0adc12a102845cd0598d332537ee9b47ed` |
| 失败测试先行 + 不变量测试 | `pipeline/tests/test_design_phaded_discrimination_panel.py` | 30,069 | `7b9edce932a2704b4ebdf80ce318b001353b4820852f75b0e8f6ede60e15a499` |
| 版本化声明输入（轴定义） | `runs/.../inputs/discrimination_axis_definitions.tsv` | 3,138 | `59ed44431956915a…`（全长见 `run_manifest.json`） |
| 版本化声明输入（绑定规则） | `runs/.../inputs/discrimination_axis_assignment_rules.tsv` | 2,545 | `3f16262bf1da1dc7…` |
| 版本化声明输入（外部锚点声明） | `runs/.../inputs/external_axis_anchor_declarations.tsv` | 10,741 | `b38ec65d7387d578…` |
| 输入契约 / run manifest | `runs/.../input_contract.json`、`run_manifest.json` | — | 15 条输入全部 `verified`；GTDB taxonomy/metadata/tree 如实写 `pending` |

两个结果文件与生成脚本的一致性用 `logs/check_determinism.py` 验证：把设计阶段重跑到临时目录，panel TSV **逐字节相同**（`0431a9ed…`），`axis_reachability.json` 除 `generated_at` 时间戳外内容完全相同。

新 run 目录含 `logs/`、`inputs/`、`results/`、`input_contract.json`（由 `pipeline/scripts/run_context.py` 生成，契约 `run_id` 与目录名一致）。

---

## 3. 四条判别轴与两端锚点状态

轴、端状态、允许的证据种类（可锚定 / 仅可作 challenge）与绑定优先级规则都写在 `inputs/discrimination_axis_definitions.tsv` 与 `inputs/discrimination_axis_assignment_rules.tsv` 中，脚本逐条校验；证据种类若未在该轴该端声明，脚本直接 fail-closed 报错，不允许"顺手放到某个端"。

| 轴 | 正端状态 | 负端状态 | 正端锚点状态 | 负端锚点状态 |
|---|---|---|---|---|
| `axis_1_substrate_chain_length` | SCL/PHB 底物（实验阳性） | MCL 非 PHB 底物 | 有实验阳性锚点（3 条台账 + 3 条外部），family 全部 unresolved | 有显式实验阴性锚点（Q84C08），family unresolved |
| `axis_2_particle_state` | denatured-only 胞外 | native-active PhaZ7 型 | 有实验阳性锚点（20 条 Knoll 2009 种子，family-resolved） | **无显式实验阴性记录**（4 条仅 annotation_only / 阳性） |
| `axis_3_localization` | 胞外 + 信号肽 + SBD | 胞内 + 无信号肽 + 检出 lid | 仅外部 acquisition 记录（5 条实验阳性，family unresolved） | 有显式实验阴性锚点（O87189、Q7WT48、Q7WT49），family unresolved |
| `axis_4_fold_neighborhood` | PF10503 域阳性 | 同 fold 脂酶/酯酶 + 极性 x1 | 域证据锚点（260 条 PF10503 命中，evidence_status 仍为 annotation_only） | **无显式实验阴性记录**（563 条极性 x1 + 34 条 α/β 水解酶域） |

### 3.1 轴定义的硬约束

- `negative_anchor` **只能**由 `experimental_negative_assay` 产生；任何标签、域状态、模式状态落在负端都只能是 `challenge`；
- 反之，`annotation_only` 记录被测试与实现双重禁止获得 `formal_negative_eligible=true` 或 `gate_credit_slot=family_resolved_negative`；
- lid 模式证据（`supported_pattern_state`）只在 `axis_3_localization:negative` 端声明，放到任何其他轴/端直接报错；
- **PF06850 绝不作为任何轴的锚点**（脚本对 `anchor_evidence_basis` 含 `pf06850` 直接报错），落实 Task 2 的 DEFECT_D：PF06850 只命中 intracellular nPHASCL without lipase box（274/276），在 284 条胞外 dPHASCL type 1 与 71 条 type 2 参考中 0 命中，因此它不是胞外 SBD 标记；
- 候选宇宙（Task 3 的 563 条）**只能**使用 `candidate_polar_x1_pattern_state` 一种证据：本 run 的输入契约（Task 2 的参考层 lid 结果 + Task 3 的 563 条）不含候选层全库 motif 补全层。该补全层在本 run 生成期间由**并行 agent 的另一个 run** 产出（`runs/20260917_phaded_motif_completion_full_01`，`results/motif_completion_full.tsv` 291,698,708 B；lid 独立检出 989 supported + 1,910 partial_loop1 + 11 partial_loop2 = 2,910 条任一环检出，SBD 判据仍为 `not_assignable`，与 DEFECT_D 一致）。本 run **未消费**该 run 的任何文件，也未修改它；据此把候选层接入 `axis_3` 负端应另开 run，且仍只能是 `challenge`（候选 accession 无 family 绑定）。

---

## 4. Panel 组成（独立重算，`logs/audit_panel_outputs.log`）

| 维度 | 计数 |
|---|---|
| panel 记录总数 / 唯一 accession | 960 / 960 |
| role：positive_anchor / negative_anchor / challenge | 291 / 4 / 665 |
| axis_end：positive / negative | 293 / 667 |
| evidence_status：experimental_positive / experimental_negative / annotation_only | 40 / 4 / 916 |
| family_binding_status：resolved_existing_family / unresolved | 379 / 581 |
| 逐轴记录数 | axis_1 = 18、axis_2 = 24、axis_3 = 60、axis_4 = 858 |
| **formal_negative_eligible = true 的记录** | **0** |
| 未锚定记录（如实列出，不硬塞进轴） | 344 条台账行 + 10 条外部记录 |

外部 acquisition 的 28 条 accession 逐条裁决（`axis_reachability.json → external_records_resolution`）：18 条作为轴锚点纳入；10 条不纳入，理由分别为
`fragment_below_completeness_threshold`（J7K890、Q9AGB6）、
`represented_by_ledger_row`（P26495、AAA87070.1、BAF35850.1、AAL30107.1、AAB40611.1）、
`evidence_conflict`（AAF86381.1：台账 `annotation_only` vs acquisition `experimental_positive`；Q939Q9：同一 run 内既有 fragment 排除又有阳性面板行）、
`ambiguous_axis_evidence`（PCH05226.1：定义行同时含 Esterase 与 PHB depolymerase，且无 PMID/DOI）。

344 条未锚定台账行的构成是：273 条 intracellular nPHASCL without lipase box（仅 PF06850 命中的 Cys 型）、59 条 extracellular dPHASCL type 2（无 PF10503）、12 条 extracellular dPHASCL type 1（无 PF10503）。它们**没有任何一条已声明轴端上的 supported 状态**，其轴归属若成立只能靠 `reported_localization`/`substrate_class` 标签——而 Knoll 2009 自己记录了三处注释与序列归属矛盾，计划 §1.2 明确"注释不可作分类依据"，因此本 run 把它们如实列为未锚定，而不是用标签硬填。

---

## 5. 逐轴可达性：当前能提供多少个锚定明确轴的 family-resolved 阴性

| 轴 | 负端锚定记录 | 其中显式实验阴性 | 其中 **family-resolved 且可计入门槛** | 距"每 profile ≥1 个 family-resolved 阴性"的缺口 |
|---|---:|---:|---:|---|
| `axis_1_substrate_chain_length` | 12 | 1（Q84C08） | **0** | 1 |
| `axis_2_particle_state` | 4 | 0 | **0** | 1 |
| `axis_3_localization` | 53 | 3（O87189、Q7WT48、Q7WT49） | **0** | 1 |
| `axis_4_fold_neighborhood` | 598 | 0 | **0** | 1 |
| **合计** | **667** | **4** | **0** | 37 个 profile 各缺 1（合计 37） |

为什么这 4 条显式实验阴性**都不能计入**：

1. 它们**没有任何一条绑定到既有 family**（`family_binding_status = unresolved`，`bound_family` 为空）——Q84C08、O87189、Q7WT48、Q7WT49 都不在 723 条冻结台账中，也没有唯一无版本号匹配；
2. 冻结台账本身**一个实验阴性记录都没有**（723 行的 `evidence_status` 只有 `annotation_only` 693 + `experimental_positive` 30，`positive_negative_control` 只有 `not_assessed` 693 + `positive_seed` 30），因此任何"台账绑定 + 显式实验阴性"的组合在现有证据里**不可能存在**；
3. 667 条负端记录中有 35 条确实绑定到 reference-only profile，但**全部是 annotation_only**（4 条 hfam_7 的 native-active 参考 + 31 条 fold 邻域记录），因此只能作 challenge；且它们只覆盖 **4 个** profile（hfam_2、hfam_7、hfam_67、hfam_71）；
4. MCL 端（hfam_4/hfam_8）与 lid 端（hfam_4）的锚点**全部属于已训练 family**，不是 33 个 reference-only 家族中的任何一个，所以对 37 个 profile 天然无效。

---

## 6. 距离校正门槛还差多少

门槛逐字取自 `pipeline/scripts/finalize_phaded_reference_panel_calibration.py` L75（fail-closed）：**每 family ≥3 独立阳性 + ≥1 held-out + ≥1 family-resolved 阴性 + ≥1 challenge + 零未解释命中**。本 run 的 `axis_reachability.json → gate` 原样复制该定义，并显式写入 `gate_relaxation_applied = false`、`gate_relaxation_supported = false`。

| 槽位 | 37 个 reference-only profile 的现状 | 缺口合计 |
|---|---|---:|
| 阳性（≥3） | 6 个 profile 有 1–2 条（最多 `DED_hfam_53` 的 2 条），31 个为 0；本 panel 新增可用阳性 **0** | 104 |
| held-out 阳性（≥1） | 全部 0 | 37 |
| **family-resolved 阴性（≥1）** | **全部 0**（本 panel 提供 0） | **37** |
| challenge（≥1） | 32 个 profile 已有（来自既有 acquisition amendment），5 个为 0（`DED_hfam_3` 与 4 个超家族 profile） | 5 |
| 零未解释命中 | **无法评估**：37 个 reference-only profile 连 HMM 都没有（`hmm_sha256` 为空、`hmm_path` 为空、`alignment_sha256` 为空），本 panel 无法对它们打分 | — |
| 通过全部槽位 | **0 / 37** | — |

因此：**本 panel 目前对 gate 的净贡献是 0 个槽位**，37 个 profile 全部仍然 `reference_only_insufficient_panel`。这不是本 run 的执行缺陷，而是证据本身的缺口，必须由新的 acquisition/建模步骤补齐（见 §7），**不得通过放宽任何阈值解决**。

补充的两个结构性事实（本 run 新识别，供后续决策）：

1. **跨家族轴对照当前不被任何计数器承认。** `reconcile_phaded_reference_panel.py` L135–147 与 `acquire_phaded_reference_panel_amendment.py` L184–191 只统计"绑定到该 profile 自身 family"的记录，因此"用 MCL 酶去挑战胞外 dPHASCL profile"这类**跨家族**轴对照，在现行 gate 计数下贡献恒为 0。轴化 panel 要真正生效，需要的是**一个新的、经评审的计数口径决策**，而不是本 run 自行改写 gate。
2. **37 个 reference-only profile 从未构建 HMM**，所以"零未解释命中"这一条在它们身上是**不可评估**而非"通过"；本 run 把它写成 `not_assessable_no_profile_hmm`。

---

## 7. 补齐缺口需要什么（明确列出，本 run 不做）

| 缺口 | 需要什么 | 归属 |
|---|---|---|
| 4 条轴上的 family-resolved 阴性 | 带**显式实验阴性**陈述的 accession，且能唯一绑定到 reference-only family；现有 4 条阴性都必须先解决 family 归属（它们不在 723 条台账内，需要新的 family 归属证据或新的 acquisition run） | 新 acquisition run（需记录检索日期与哈希） |
| axis_2 负端（native-active） | native-active 端的显式阴性实验记录；现有 PhaZ7 参考只有 annotation_only 或"能降解 native PHA"的阳性 | 文献 acquisition |
| axis_4 负端（脂酶/酯酶） | 脂酶/酯酶**功能性阴性**记录（而非极性 x1 模式）；563 条极性 x1 中仅 200 条（35.5%）另有 α/β 水解酶 fold 域证据，358 条只有 PF06850（Cys 型标记，不是 fold 证据），5 条无任何 tested Pfam | 文献 acquisition + Pfam 证据补充 |
| 344 条未锚定 reference-only 行 | 经独立验证的胞内/胞外判别证据（本 run 拒绝用 `reported_localization` 标签锚定） | Task 2 遗留 blocker（胞外 SBD 判据、linker、Cys 型亲核体坐标） |
| 候选宇宙的 lid/轴证据 | 候选层全库 motif 补全层已由并行 run `20260917_phaded_motif_completion_full_01` 产出（lid 989 supported + 1,921 partial），但**未被本 run 消费**；把它接入 `axis_3` 负端需要新 run，且候选无 family 绑定，只能作 challenge | 新 run（Task 7 产物已就绪） |
| 零未解释命中 | 37 个 profile 的 HMM 必须先构建（≥3 独立阳性是前置条件） | Task 6/后续授权 |

---

## 8. 本 run 顺带发现、只报告不修改的问题

1. **另一个 run 的字段错位（P 级，实测）**：`runs/20260916_phaded_reference_panel_acquisition_03/inputs/external_panel_candidates.tsv` 表头 9 列、数据行 10 列，导致 `bound_family` 之后**整体右移一列**：`pmid` 为空、`doi` 列存放的是 PMID、`source_database` 列存放的是 DOI。其结果文件 `results/panel_candidates.tsv` 同样受影响（`sequence_sha256` 列在该文件中为 `NCBI Protein`）；`results/panel_candidates_bound.tsv` 的 `sequence_sha256` 已被后续 enrichment 覆盖为正确哈希，但 pmid/doi/source_database 仍是错位值。**本 run 未修改该 run 的任何文件**；引用 PMID/DOI 时改用未错位的 `external_panel_evidence_prior.tsv` 与 acquisition 报告的 `literature_pmids`。
2. **证据冲突（未裁决）**：`AAF86381.1` 在台账为 `annotation_only`，在 acquisition 面板为 `experimental_positive`；`Q939Q9` 在同一 run 内既被 `ephaz_excluded_challenge.faa` 标为 fragment 排除、又出现在阳性面板行；`P26495`（无版本号）在 ephaz 上下文被当作 annotation_only near-neighbor negative，而台账 `P26495.1` 是 `DED_hfam_4` 的 `experimental_positive` 种子。本 run 一律**不做升级也不做降级**，保持不纳入并在 `axis_reachability.json` 记录理由。
3. **轴端共享是设计固有属性**：`axis_3` 的胞外正端与 `axis_4` 的 PF10503 正端描述的是同一批胞外 dPHASCL 记录。为满足"每条记录恰好绑定一个轴"，本 run 用**声明过的优先级规则**做唯一指派，并另设 `secondary_axis_capability` 列如实记录该记录同时满足的其他轴端（例如 `P12625.1` 记录为 `axis_4_fold_neighborhood:positive;axis_3_localization:positive_label_only`），不隐藏共享。
4. **正交性声明**：本 run 未使用 Task 2 已知自相矛盾的 His 共有模式（`GMXHXXPXXG` 10 残基 vs 文献例 `GMGHAWSGG` 9 残基），也未把 AHSMG 零检出写成"该族不存在"。

---

## 9. 复现与验收命令

```powershell
# 1) 生成声明输入（幂等重写本 run 自己的 inputs/）
python runs/20260917_phaded_discrimination_panel_01/logs/author_discrimination_panel_inputs.py
# 2) 设计 panel（结果文件已存在时拒绝覆盖）
python pipeline/scripts/design_phaded_discrimination_panel.py
# 3) 聚焦测试
python -m unittest pipeline.tests.test_design_phaded_discrimination_panel -v
# 4) 独立重算（不 import 设计脚本）
python runs/20260917_phaded_discrimination_panel_01/logs/audit_panel_outputs.py
# 5) 结果文件与脚本的可复现性（重跑到临时目录后逐字节比对）
python runs/20260917_phaded_discrimination_panel_01/logs/check_determinism.py
```

| 命令 | 结果 |
|---|---|
| 聚焦测试 `test_design_phaded_discrimination_panel` | 先红：`FileNotFoundError: pipeline/scripts/design_phaded_discrimination_panel.py`（`logs/step5.1_failing_test_red.log`）；后绿：**`Ran 38 tests … OK`（无 skip）**（`logs/step5.2_focused_tests_green.log`） |
| `python -m unittest discover -s pipeline/tests` | **`Ran 450 tests in 10.828s` / `OK (skipped=1)`**（`logs/full_test_suite.log`）。450 = 本阶段基线 375 + 本 run 新增 38 + 并行 agent 在同期新增的测试；无失败、无 error |
| `python -m compileall -q pipeline` | exit 0（`logs/compileall.log`） |
| `git diff --check` | exit 0（`logs/git_diff_check.log`） |
| 独立重算 | `logs/audit_panel_outputs.log`：960 条 / 唯一 960；负端 667；显式实验阴性 4；family-resolved 可计入 0；HMM 已构建的 reference-only profile 0/37 |
| 可复现性 | `logs/check_determinism.log`：`deterministic` |

设计阶段的 fail-closed 不变量（全部有测试）：每条记录**必须且只能**绑定一个已声明轴与轴端；`anchor_justification`、`anchor_evidence_basis`、`source_doi_or_pmid` 不得为空；`negative_anchor` 必须携带 `experimental_negative`；`annotation_only` 永不得为 formal negative 或 gate-eligible；未绑定既有 family 的记录必须保持 `unresolved` 且不得携带 `bound_family`；`ambiguous_existing_family` 不得被升级；证据种类必须在该轴该端被声明；候选宇宙不得使用 lid 证据；PF06850 不得作为锚点；输出存在时拒绝覆盖。门槛保护由三层测试共同断言：(a) `GATE` 常量必须与 `finalize_phaded_reference_panel_calibration.py` L75 的 `positive_count >= 3` / `heldout >= 1` / `negative >= 1` / `challenge >= 1` 逐条一致（并直接在该源文件中搜索这些片段）；(b) panel 字段名、模块属性名中不得出现任何放宽门槛的词汇；(c) `axis_reachability.json` 必须**显式**声明 `gate_relaxation_applied = false` 与 `gate_relaxation_supported = false`，且除这两个显式声明字段外，整份 JSON 中不得出现相关词汇。

---

## 10. 未 commit / 未越界核对

- `git rev-parse HEAD` = `2e44aafe2c91263447f37bac625f97c56ddafc05`（与计划基线一致），本 run **未 commit / push / clean / reset**。
- `pipeline/config/formal_scan_models.tsv` = `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`（与 Phase 1 验证报告记录值一致，未变）。
- `AGENTS.md` = `e42adbe3cb7d5b6e6cd713c5c3a797060ec638f838844a8811df48b8a5c71f2e`，mtime 2026-09-12 18:32:37（**早于本阶段**，非本 run 改动）。
- 扫描 `runs/20260901_*`–`runs/20260916_*` 全部历史 run 目录：**没有任何文件在 2026-09-17 当天或之后被修改**（0 违规）；本 run 只新建了 `runs/20260917_phaded_discrimination_panel_01/`，并新增 `pipeline/scripts/design_phaded_discrimination_panel.py`、`pipeline/tests/test_design_phaded_discrimination_panel.py`、本状态文档。
- 未修改其他 agent 的产物（含并行 run `20260917_phaded_motif_completion_full_01`、`20260917_phaded_p2_p4_remediation_01`）：对这些 run 只读引用并报告观察。

---

## 11. 边界复述

本 run 全部产出是 candidate-only 的序列/域/模式与文献锚点证据组织结果。每条 challenge panel 记录仍然只表示"该序列在某个判别轴的某一端具备某类证据"，**不表示** PHB/PHA 降解表型、不代表任何 family 判定，也没有改变任何校准阈值：`gate_relaxation_applied = false`、`new_family_call_created = false`、对 109,087 条候选的分类结果为 **0 行影响**。
