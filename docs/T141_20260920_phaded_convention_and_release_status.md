# T141 — 口径统一（C1–C10）与发布准备（审查报告第 6 项）：状态文档

**日期：** 2026-09-20
**性质：** 口径落版 + 发布**准备**（**不执行任何推送 / 不创建 DOI**）
**触发：** 项目全面审查报告「接下来该做什么」第 6 项
**边界：** candidate-only。**每条口径的数值均在本 run 于冻结产物上独立实测**，逐条给出产物路径；未凭空裁决任何数值，未编造任何哈希。

---

## A. C1–C10 最终引用口径

> 来源：交接文档 §5（`docs/T141_20260917_project_handoff.md`，C1–C10）。本节把每条冲突**落定为单一引用口径**，并在冻结产物上逐条复核。**未编辑交接文档**（由操作者统一更新）。

### A.1 汇总表

| 条 | 最终引用口径（一句话） | 关键数值 / 分母 | 复核产物 |
|---|---|---|---|
| C1 | 状态计数与真实传递数**分开写**；禁止单独说「48/48/1 来自锚定列传递」 | 状态计数 48(his)/48(asp)/1(oxyanion)；真实传递 45/45/0；分母 723 条参考 | `runs/20260917_phaded_reference_residue_mapping_01/results/reference_mapping_manifest.tsv`；`runs/20260917_phaded_p2_p4_remediation_01/results/p2_p4_remediation_report.md` |
| C2 | **65 = 比对行数；64 = 面板内真实带锚行数**，定义不同不可互换 | 65（=46+13+6）；64；分母 723 | `runs/20260917_phaded_motif_completion_full_01/results/motif_completion_summary.json` |
| C3 | **外部面板声明 6 → 锚入面板 4 → 绑定 family 0** | 6 / 5 / 4 | 见 A.3（含路径更正） |
| C4 | **必须写明是哪一层**：0915/计划/Task 4 层 = 373；Task 7 层 = 无此状态 | 373；分母 109,087 | `runs/20260915_phaded_motif_reconciliation_01/results/subtype_matrix_manifest.json`；`runs/20260917_phaded_catalytic_domain_type_01/results/type_summary.json` |
| C5 | 读冻结 `motif_completion_summary.json` 的 `criteria.*.before` 时**必须同时读旁证** field_notes | 旧字段名 + lid/linker/sbd 硬编码 0（应读 `not_available`） | 同名 run 的 `results/motif_completion_summary_field_notes.json` |
| C6 | **不写死数值**；计划 v2 是活文档，以磁盘实测为准，快照只在 manifest | 历史序列（均非当前值）：23,378 / 24,158 / 29,852 / 31,225 / 33,695 B | `runs/20260917_phaded_final_consistency_fix_01/results/final_consistency_fix_manifest.json`；`runs/20260917_phaded_closeout_01/results/closeout_manifest.json` |
| C7 | 引用时必须写明「**来自计划 v2 事实基线，未经本批独立重算**」 | 62,666 基因组（未重算；本地 gtdb = `pending`） | 计划 v2；交接文档 §6.4 |
| C9 | **必须写明是哪一层**：31,632 = subtype 判定层；31,634 = 先验层 | 31,632 / 31,634（差 2 行已定位） | `runs/20260915_phaded_motif_reconciliation_01/results/subtype_matrix_manifest.json`；`runs/20260917_phaded_inphascl_cys_decision_01/results/cys_decision_manifest.json` |
| C10 | **ahsmg 的 0 是实测值；lid/linker/sbd 的 0 是硬编码**，不可混谈 | ahsmg=0（measured）；lid/linker/sbd=0（`not_available`） | 同 C5 两文件 |

### A.2 C1 复核细节（本 run 实测）

`reference_mapping_manifest.tsv`（723 行）实测：

| 列 | `supported_anchor_column` | `supported_literature_pattern` | 其它 |
|---|---:|---:|---|
| `his_state` | **48** | 42 | pending 633 |
| `asp_state` | **48** | 23 | not_detected_at_mapped_column 2；pending 650 |
| `oxyanion_hole_state` | **1** | 255 | pending 467 |

→ 确认 **48/48/1 是 `supported_anchor_column` 的状态计数**；真实传递数 45/45/0 出自 `p2_p4_remediation_report.md`（Phase 2 独立重算）。**引用规范：写状态计数时标 `supported_anchor_column`；写传递数时标「去掉锚行自身文献坐标后的传递」。**

### A.3 C3 复核细节 + **一处路径更正（本 run 发现）**

| 数值 | 含义 | 实测产物 |
|---:|---|---|
| **6** | 外部面板**声明**的显式实验阴性 | `runs/20260917_phaded_discrimination_panel_01/inputs/external_axis_anchor_declarations.tsv`（28 行中 `evidence_status=experimental_negative` = **6**；其中 `include=true` = **4**：Q84C08、O87189、Q7WT48、Q7WT49） |
| **5** | 绑定进面板候选表的条数 | `runs/20260916_phaded_reference_panel_acquisition_03/results/panel_candidates_bound.tsv`（**5** = 上述 4 条 + `J7K890`） |
| **4** | 锚入判别面板的条数 | `runs/20260917_phaded_discrimination_panel_01/results/challenge_panel_by_axis.tsv`（**4**） |

> **⚠️ 路径更正（本 run 实测发现）：** 交接文档 §5 C3 把第三个来源写作 `runs/20260916_phaded_reference_panel_acquisition_03/inputs/external_axis_anchor_declarations.tsv`——**该路径不存在**。实际文件位于 `runs/20260917_phaded_discrimination_panel_01/inputs/external_axis_anchor_declarations.tsv`。**本 run 只报告该路径差异，不改写交接文档**（由操作者统一更新）。

**最终口径：** 「外部面板声明 **6** 条显式实验阴性 → 绑定进面板 **5** 条 → 锚入判别面板 **4** 条 → 其中 **0** 条绑定 family（均 `family_binding_status=unresolved`、不可计入 gate 的 family-resolved 阴性槽位）。」

### A.4 C4 复核细节（本 run 实测）

| 层 | 状态 | 实测产物 |
|---|---|---|
| 0915 层 | `motif_evidence_status.conflict_relative_position` = **373** | `runs/20260915_phaded_motif_reconciliation_01/results/subtype_matrix_manifest.json`（实测 `{'partial_catalytic_pattern_panel': 76929, 'conflict_relative_position': 373, 'not_tested_full_library': 31785}`） |
| Task 4 层 | `type_counts.undetermined_position_conflict` = **373** | `runs/20260917_phaded_catalytic_domain_type_01/results/type_summary.json`（实测 `{'type1_verified': 16530, 'type2_verified': 2362, 'undetermined_no_oxyanion': 89822, 'undetermined_position_conflict': 373}`） |
| Task 7 层 | **无该状态**（`oxyanion_hole_state` 仅两态） | `runs/20260917_phaded_motif_completion_full_01/results/motif_completion_full.tsv` |

**引用规范：** 说「373 条 conflict」时必须写明「**0915/Task 4 口径**」；讨论 Task 7 层时不得声称有 conflict 状态。**另须标注** Task 7 的 oxyanion supported **18,939** vs Task 4 的 **18,902**（+37 差异，两 run supported 口径不完全同一）。

### A.5 C5 / C10 引用规范

- 冻结的 `total_library/motif_completion_summary.json` 的 `criteria.*.before` **仍**使用旧字段名 `independent_detection_supported`，且 lid/linker/sbd 的 before **仍**为硬编码 `0`。**冻结产物按规则不可改写**，因此：
  - **必须**同时读 `runs/20260917_phaded_motif_completion_full_01/results/motif_completion_summary_field_notes.json`（本 run 实测存在）——其中把 ahsmg 标为 `measured`、lid/linker/sbd 标为 `not_available`。
  - **禁止**把 ahsmg 的 0 与 lid/linker/sbd 的 0 当作同一性质引用（C10）。

### A.6 C6 引用规范

计划 v2 是**活文档**：`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md`。**任何引用都不得写死大小/SHA-256**；某一时点的快照值只写在 `runs/20260917_phaded_final_consistency_fix_01/results/final_consistency_fix_manifest.json`（实测存在）与 `runs/20260917_phaded_closeout_01/results/closeout_manifest.json`（实测存在）中。历史序列 23,378 / 24,158 / 29,852 / 31,225 / 33,695 B **全部只是历史快照**。

### A.7 C7 引用规范

「62,666 基因组」**来自计划 v2 事实基线，本批未独立重算**（需 GTDB metadata，本地 `gtdb` 输入为 `pending`）。任何引用必须带此限定语，不得表述为实测值。

### A.8 C9 引用规范（本 run 实测）

| 层 | 指标 | 值 | 产物 |
|---|---|---:|---|
| subtype **判定**层 | `counts.subtype_call.unassigned_PhaDED_like` | **31,632** | `runs/20260915_phaded_motif_reconciliation_01/results/subtype_matrix_manifest.json` |
| **先验**层 | `candidate_prior_superfamily_counts["(no_superfamily_prior)"]` | **31,634** | `runs/20260917_phaded_inphascl_cys_decision_01/results/cys_decision_manifest.json`（实测完整分布：`(no_superfamily_prior) 31634 / extracellular dPHAMCL 292 / extracellular dPHASCL type 1 67544 / extracellular dPHASCL type 2 2658 / intracellular nPHAMCL 6959`） |

**差 2 行已定位**：`GCA_964247575.1|CAXTVD010000242.1_5`、`GCA_964589695.1|CAZPUI010000267.1_2`（`candidate_assignment_status=ambiguous_superfamily` 但 `candidate_prior_superfamily` 为空）。

---

## B. 发布准备（**只准备，不执行**）

> 本节仅为清单与步骤草案。**本 run 未 push GitHub、未创建 Zenodo DOI、未做任何网络写操作。**

### B.1 发布物清单（本 run 逐一实测存在性）

**(a) HMM / profile**

| 类别 | 内容 | 路径 | 实测 |
|---|---|---|---|
| 正式 registry | 10 行模型：`ePhaZ_curated_core` / `ePhaZ_broad_discovery` / `iPhaZ` / `OH` / `BdhA` / **`ArchPhaZ_patatin`** / **`ArchPhaZ_hydrolase`** / `PhaJ` / `phasin` / `PhaC`（阈值均 e-5） | `pipeline/config/formal_scan_models.tsv` | 10 行已读取确认；**禁改** |
| 已入库 HMM | 21 个 HMM 文件 | `data/hmms/`（含 `data/hmms/v2/`） | `git ls-files data/hmms` = 21 |
| PhaDED trained profile | **9** 个 profile HMM | `runs/20260909_phaded_reference_evidence_01/results/profile_build_01/profiles/` | 9 个 `.hmm` |
| 发现层 HMM | **33** 个家族 HMM，层名固定 `discovery_hmm_uncalibrated` | `runs/20260918_phaded_family_discovery_01/results/hmms/` | 33 个 `.hmm` |

> **发布时必须区分三层**：① 正式 registry（10 模型，含古菌**召回探针**）；② PhaDED trained profile（9，但**未通过校准 gate**）；③ 发现层 33（**只召回，禁筛选**）。**三者不得混为「本项目 HMM 集」一句话。**

**(b) 轻量命中表**

| 表 | 行数 | 路径 | 实测 |
|---|---:|---|---|
| 池内高可信度 | 36,611 | `runs/20260919_phaded_high_confidence_relabel_01/results/high_confidence_candidates.tsv` | 存在 |
| 池外高可信度 | 4,131 | `runs/20260919_phaded_pool_external_evidence_01/results/filter/high_confidence_candidates.tsv` | 存在 |
| 池外 hold（含 `hold_reason`） | 26,501 | `runs/20260919_phaded_pool_external_evidence_01/results/filter/hold_candidates.tsv` | 存在 |
| 池外 score-only 层（含 `score_tier`） | 35,558（strong 30,632 / mid 1,503 / weak 3,423） | `runs/20260919_phaded_high_confidence_relabel_01/results/pool_external_profile_hits.tsv` | 存在 |

### B.2 发布物**必须附带**的声明（逐条，缺一不可）

1. **candidate-only 声明**：所有 HMM/domain/SignalP/motif/结构/定位/系统发育证据只表示**候选同源或功能潜力**，**不等同于已验证的 PHB/PHA 降解表型**。
2. **结构断点声明**：`candidate_ded_alignment_members = 0`——109,087 条候选**没有任何一条**属于任何冻结 DED 比对，故候选层的催化残基（His/Asp/Ser/Cys 身份）**无法**靠比对列传递获得。
3. **校准 gate 未解锁声明**：**没有任何家族达到 ≥3 条独立实验阳性**（Cys 型最多 2 条：PhaZ1 + PhaZ2），因此**未发布任何校准 profile**；`finalize_phaded_reference_panel_calibration.py` 的 gate 为 fail-closed，37 个 reference-only profile 全部 `profiles_meeting_full_gate = 0`。
4. **发现层声明**：33 个家族 HMM 为 `discovery_hmm_uncalibrated`，**只用于召回**，无判别力（实测对 563 条 x₁ 混淆嫌疑命中 443 条 = 78.69%）；不得用于筛选/删除/降级任何候选。
5. **with-lipase deferred 声明**：`intracellular nPHASCL with lipase box` 池外命中 1,206,655 条为已知噪声，保留为 deferred 结构验证层（`runs/20260919_phaded_with_lipase_deferred_archive_01`），**不进任何候选/高可信度计数**。
6. **古菌声明**：`ArchPhaZ_patatin` / `ArchPhaZ_hydrolase` 为「**patatin 折叠召回探针**」，其命中**不得**被读作「古菌 PHB 解聚酶」（见 `docs/T141_20260920_phaded_archaeal_hmm_convergence_status.md`）。
7. **池外 weak/mid 声明**：池外 `weak`（3,423）/ `mid`（1,503）为 score-only 层，**不升级、不进任何计数**（见 `docs/T141_20260920_phaded_weak_tier_decision_status.md`）。
8. **Cys 推定声明**：Cys 型高可信度候选标注 `catalytic_residue_verification=inferred_cys_not_verified`（**推定** Cys 型，催化残基未逐条验证，未排除无 GxSxG 的 Ser 型混入）。
9. **结构验证层 + 序列层口径声明（2026-09-20 追加，源自 with-lipase 结构验证试点）**：deferred 层（1,206,655）的结构召回路径**已实测可行但尚未执行**——PDB **8YNV**（B. thuringiensis PhaZ，链 A 298 aa，含 `GWSMG`/`DRDYVV`）可下载、Foldseek 可用（自测 TM 1.00）、ColabFold 1.6.1 可在 GPU 运行，但**单序列 AF2 预测不可用于结构比对**（同序列 TM 仅 0.20、pLDDT 36），全量 1.2M 结构预测在算力上不可行（~324 GPU-日）。因此：**未产生任何「结构通过」条目，结构结果不入任何计数**；若日后执行，须先按预注册判据（Foldseek **qtmscore ≥0.5 且 E ≤1e-3**、预测 **pLDDT ≥70** 质量门、fail-closed：无结构/低置信/无命中一律记 `not_assessed`，**不得记为阴性**）并在**独立 `structural_validation` 层**发布。**同时必须声明序列层口径**：对**单条**参考（8YNV 链 A）的 blastp，其 **E 值按单序列搜索空间计算，比 GTDB 规模宽松约 8.4×10⁵ 倍**，**不可**与 HMM E 值或任何全库 E 值比较；该类搜索**不构成过滤器**，**无命中不得读作阴性**（依据 `docs/T141_20260920_phaded_structure_verification_status.md` §2.2 更正记录——其中记录了「E 值列被误当 bitscore 列」这一**在两个独立 agent 上各复发一次**的错误，及其机器校验脚本 `runs/20260919_phaded_pool_external_evidence_01/plan/verify_pilot_blastp_columns.py`）。

### B.3 GitHub Release + Zenodo DOI 操作步骤草案（**待操作者明确授权后执行**）

> **前置授权要求（AGENTS.md）：** 「push 到 GitHub 仍需操作者明确授权」「重算、安装、配置变更、删除与推送仍须明确授权」。以下步骤**本 run 一律未执行**。

**阶段 0 — 预检（只读，可先做）**
1. `git status --porcelain` 确认工作树；记录 `git rev-parse HEAD`。
2. 跑 `pipeline/tests/test_public_repo_safety.py`（禁服务器身份串入库）。
3. 跑 `python -m unittest discover -s pipeline/tests`、`python -m compileall -q pipeline`、`git diff --check`，三者 exit 0。

**阶段 1 — 打包发布物（本地，不推送）**
4. 建发布目录（如 `release/<version>/`），放入：HMM（三层**分目录**）、轻量命中表、`B.2` 的 **9** 条声明（写成 `README_DISCLOSURE.md`）、以及每个文件的 SHA-256 清单。
5. 为每份产物生成 `SHA256SUMS`；**哈希以实测为准**。

**阶段 2 — GitHub Release（需授权）**
6. 打 annotated tag（如 `v<version>-candidate-only`）并 push tag。
7. 创建 Release，上传阶段 1 的压缩包与 `SHA256SUMS`；Release 正文**必须内嵌 B.2 全部 9 条声明**。

**阶段 3 — Zenodo DOI（需授权）**
8. 在 Zenodo 关联该 GitHub 仓库并开启 release 归档（或手动上传同一压缩包）。
9. **DOI 仅用于引用**；Zenodo 描述里同样内嵌 B.2 的 9 条声明与 candidate-only 边界。
10. 把 DOI 写回项目文档（`README` / 交接文档）——**该写入由操作者执行**。

**阶段 4 — 发布后核对（AGENTS.md：发布前核对权威关系）**
11. 核对本地 / GitHub / 服务器 `deploy/` / run manifest 四处的版本权威关系，避免混用。

---

## C. 本次发现并报告的口径/路径问题（不擅自改写交接文档）

| # | 问题 | 处理 |
|---|---|---|
| 1 | 交接文档 §5 C3 的第三个来源路径 `20260916_phaded_reference_panel_acquisition_03/inputs/external_axis_anchor_declarations.tsv` **不存在**；实际在 `20260917_phaded_discrimination_panel_01/inputs/` | 本文档 A.3 如实记录；**未改写交接文档** |
| 2 | `ArchPhaZ_patatin`（103 种子）与 `ArchPhaZ_hydrolase`（12 种子）**确实已在** `pipeline/config/formal_scan_models.tsv` 注册（10 行中的第 6/7 行）——「未注册」只适用于 8 种子的 `archaea_PhaZh1_like_candidate.hmm` | 已在 `docs/T141_20260920_phaded_archaeal_hmm_convergence_status.md` 第 2/2b 条精确区分 |

## D. 边界与未做事项

- 本 run **未改**任何历史 run 的 `results/`、未改 `pipeline/config/`、未改冻结台账、**未删除任何文件**。
- 本 run **未编辑** `docs/T141_20260917_project_handoff.md`（由操作者统一更新）。
- 本 run **未** `git commit` / `git push` / 创建 Zenodo DOI / 做任何网络写操作。
- **待操作者授权**：B.3 的全部执行阶段，以及是否把 A 节口径写回交接文档。
