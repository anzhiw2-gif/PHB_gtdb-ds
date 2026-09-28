# PhaDED 分类权威口径与架构判据验证计划（v2）

> **For agentic workers:** 按 Task/Step 执行，代码改动一律**先写失败测试，再实现**。每个 Task 使用新的 dated run（如需服务器则用 dated deploy）。本文件**取代**同日 v1 草稿（v1 未被执行，差异见文末"修订记录"）。

**制定日期：** 2026-09-17（v2）
**工作区：** `<REPO_ROOT>`（分支 `main`，HEAD `2e44aaf`，自 2026-09-09 未变）
**服务器：** `<SERVER_USER>@<SERVER_HOST>:${PHB_REMOTE_ROOT}/PHB_gtdb-ds`（SSH 用 `C:\Windows\System32\OpenSSH\ssh.exe`）
**授权范围：** 本地仅改动 `<REPO_ROOT>`；服务器仅改动 `${PHB_REMOTE_ROOT}/PHB_gtdb-ds`；服务器其他路径只读，确需改动时先复制进 PHB_gtdb-ds 再改。

---

## 0. 计划目的：三条主线

| 主线 | 内容 | 对应任务 |
|---|---|---|
| **A. 口径修复** | 超家族（功能先验，可进 registry）与家族（2009 聚类细化层，禁止 registry 资格）错配问题 | Task 1 |
| **B. 判据补全与验证** | 把文献四判据（x₁ 疏水 / 催化域 type 1-2 / lid / AHSMG）在现有数据上落地、验证、去循环化 | Task 2–4、7 |
| **C. 解锁校准** | 按判别轴重设计 challenge panel；唯一可合理绕开门槛的 Cys 型 profile 论证 | Task 5–6 |

外加 **Task 0**（收尾两件已知小事）与贯穿性的口径统一。

**边界（不可突破）：** 所有产出仍为 candidate-only，表示候选同源或功能潜力，**不等同于已验证 PHB/PHA 降解表型**。

---

## 1. 事实基线（v2 重新核查，全部可复核）

### 1.1 项目状态（2026-09-17 00:52 复核，与 v1 时无变化）

- Git：分支 `main`，HEAD `2e44aaf`；自 09-17 00:41 以来工作区唯一新增文件是 v1 计划本身。
- 候选宇宙：109,087 蛋白 / 62,666 基因组；46 profiles = 9 trained + 37 reference-only（33 family + 4 superfamily）；27 个 family 零实验阳性。
- 门槛：≥3 独立阳性 + ≥1 held-out + ≥1 family-resolved 阴性 + ≥1 challenge + 零未解释命中（`finalize_phaded_reference_panel_calibration.py` L75，fail-closed）；37 个 profile 的阴性/held-out 全为 0 → **降门槛无法解锁任何 profile**。
- subtype 判定只由 `phaded_superfamily_best` 决定（`build_phaded_subtype_matrix.py` L153–165）；`phaded_family_best` 中从未出现 reference-only family → **合并/删除 27 个家族对候选为 0 行影响**。
- 测试基线：277 tests OK（skipped=1）（受限沙箱下 138 个 error 均为 tempfile 权限，非代码缺陷）。

### 1.2 文献依据（Knoll 2009 等）

- PhaDED = 8 超家族（28 条实验验证种子的**功能/定位先验**）+ 38 同源家族（**序列相似性 + 系统发育聚类**）[PMC2666664](https://europepmc.org/article/PMC/PMC2666664)。
- DED 自 2009 年未更新（Database Commons：`Last update: 2009`, v1.1）[链接](https://ngdc.cncb.ac.cn/databasecommons/database/id/3446)；但 2025 综述仍沿用该框架 [PMC11893044](https://europepmc.org/article/PMC/PMC11893044)。
- 四判据：**x₁ 疏水性**（脂酶/酯酶 x₁ 极性 vs PHA 解聚酶 x₁ 疏水）；**type 1/2** = oxyanion hole 相对 lipase box 的位置（N 端/C 端）；**lid**（i-nPHAMCL 有、e-dPHAMCL 无，[PMC12504323](https://europepmc.org/article/PMC/PMC12504323)）；**AHSMG**（PhaZ7 型替代基序，与 Bacillus 脂酶 abH18.01 共享，单独不充分）。
- Knoll 原文记录三处注释与序列归属矛盾（gi:194292521、gi:74267419 归入 e-dPHASCL type 1；gi:34452171 归入 i-nPHAMCL；PhaZ7 定位例外）→ **注释不可作分类依据**。

### 1.3 判据数据现状（v2 新核查：候选层 + 参考层对照）

**候选层**（`runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv`，109,087 行）：

| 列 | 分布 |
|---|---|
| `lipase_box_x1` | 未测 31,785；疏水（A,C,F,I,L,M,V,W,Y）76,739（99.27%）；**非疏水 563**（D 89/R 84/K 82/G 80/E 62/T 47/Q 36/S 31/P 30/H 12/N 10） |
| `oxyanion_hole_state` | not_detected 89,812 · supported 18,902 · conflict_relative_position 373 |
| `ahsmg_state` | not_detected_pattern **109,087（全部）** |
| `lid_state` | not_expected_for_subtype 101,836 · pending_reference_annotation 7,251（**由 subtype 先验赋值，循环论证**） |

**参考层**（`reference_motif_panel.tsv`，723 条参考；v2 新增核查）：

| 列 | 分布 | 结论 |
|---|---|---|
| `lipase_box_state` | supported 441 | ✅ 可用 |
| `oxyanion_hole_state` | supported 263 | ✅ 部分可用 |
| `asp_state` | supported 仅 25（全是 e-dPHASCL type 1，坐标约 238–248） | ⚠️ 只对 type 1 生效 |
| **`his_state`** | **not_detected_pattern 723（全部，含本应命中 His 的 P12625.1 等）** | ❌ **His 模式从未触发**，三联体检测实质缺失 |
| `ahsmg_state` | **supported 2**（`AAK07742.1`、`2VTVA`，均为 PhaZ7 家族参考） | ✅ **阳性对照已存在且通过** → 候选全库 0 命中是真实结果：候选池不含 PhaZ7 型 |
| `sbd_state` | pending_tool_input 717 | ⚠️ 从未运行（可由既有 Pfam PF06850 数据补齐） |
| `linker_state` | pending_reference_annotation 717 | ⚠️ 从未绑定 profile |
| `lid_state` | pending 58（i-nPHAMCL 52 + e-dPHAMCL 6） | ⚠️ 循环赋值，需独立检测 |

**v2 关键结论（改变计划形态的两点）：**

1. **AHSMG 问题已解决**：阳性对照（2VTVA = PhaZ7 的 PDB 结构 2VTV 链 A）已通过。v1 中"先做阳性对照再下结论"的 Task 降级为**事实固化与报告**（Task 4 的一个 step）。结论表述：候选池以已验证的模式**未检出** PhaZ7 型，**不得**表述为"该族不存在"。
2. **真正的技术缺口是 His/Asp/lid/linker/SBD 参考映射未完成**：`audit_phaded_motifs.py` 如实记录了 `cys_his_asp_reference_mapping_pending` / `lid_profile_not_bound` / `linker_profile_not_bound` / `sbd_reference_or_pfam_binding_pending`。当前 motif 层实质是"**lipase box + oxyanion hole**"面板，不是完整催化三联体面板。**这是比 x₁ 分层更优先的补全项**（Task 2），因为：type 1/2 几何验证、催化完整性声明、lid 判别、challenge panel 的可操作性全部依赖它。

### 1.4 已知待收尾（Task 0 素材）

- provenance：SignalP 实际在 T141 跑完（返回码 0），但三个 09-16 contract 仍写 `server_execution_started:false`、最新 reconciliation manifest 仍写 `execution.mode=local_read_only`；`amend_phaded_provenance.py` + 测试已写好且通过，未执行。
- 计数不一致：`subtype_call` 的 `ambiguous_family_within_superfamily` 9,631 / `ambiguous_superfamily` 41 与 profile 层 9,664 / 42 的差值（33、1）未解释。

---

## 2. 全局约束

- 新输出只写新的 `runs/<run_id>/`；**不得覆盖**历史 run、`results/`、registry、正式扫描结果。`run_id` 可审计、无路径遍历；目录含 `logs/`、`inputs/`、`results/`、`input_contract.json`。
- 输入记录路径、版本、大小、SHA-256；缺失写 `pending`，不得伪造哈希。
- 服务器只执行 dated `deploy/<run_id>/` 中绑定源码；单任务线程 ≤ 40。
- 代码改动先写失败测试；完成后跑相关测试、`python -m compileall -q pipeline`、`git diff --check`。
- 不得修改 `pipeline/config/formal_scan_models.tsv`、历史 HMM、`ePhaZ_curated_core`、历史 run。
- **不得启动 GTDB 全量重扫**，除非 registry decision record 评审通过并获得新授权。
- 保留所有失败重试、hold、challenge、negative、superseded 证据；不删除。

---

## 3. 任务清单

### Task 0（P0，微）：收尾两件已知小事

**Files:** 已存在 `pipeline/scripts/amend_phaded_provenance.py` + 测试；新建 run `runs/20260917_phaded_housekeeping_amendments_01/`

- [x] **Step 0.1** 执行 provenance amendment（沿用 v1 Task 0 的 6 步：基线测试 → dated run → 执行 → 校验 `server_dated_deploy`/`returncode=0`/`no_historical_run_overwritten=true` → 记录哈希 → 状态文档，声明原 contract 不可追溯修改）。
- [x] **Step 0.2** 计数不一致 amendment：只读核对 `reconciliation_manifest.json` 中 9,631/41 与 9,664/42 的来源行，写出差值解释（预期为 33 条与 1 条在 subtype 层与 profile 层的归类差异）；**不覆盖**历史 manifest。
- [x] **Step 0.3** 状态文档 `docs/T141_20260917_phaded_housekeeping_amendments_status.md`；跑完整测试 + `compileall` + `git diff --check`。

> **执行结果（Task 0）：** `runs/20260917_phaded_housekeeping_amendments_01` — **completed**。34 条 reclassified（33 `profile_ambiguous_family` + 1 `profile_ambiguous_superfamily`，全部 `subtype_call=hold_architecture_conflict`，`resolution=hold_branch_precedence`），`unexplained_reclassified_count=0`，差值 33/1 解释闭合；Phase 1 验证报告 P1（contract 不合规）已在同 run 修复（`status` 域、`gtdb` 三输入 `pending`），pre-fix 字节由 `results/input_contract_superseded_pre_p1_fix.json`（2141 B / `91b4b3ee…`）保留；根因未闭合见 Phase 2 验证报告 V2。

### Task 1（P0）：分类权威口径（层级错配修复，零数值影响）

**Files:** 新建 `pipeline/scripts/validate_phaded_classification_authority.py` + 测试；run `runs/20260917_phaded_classification_authority_01/`；`inputs/phaded_classification_authority.tsv`；`results/authority_validation.json`；状态文档

- [x] **Step 1.1（失败测试先行）** 断言：8 个 superfamily 全部 `source_evidence_type=functional_prior` 且 `registry_eligible=true`；38 个 family 全部 `sequence_clustering_2009` + `source_version="Knoll2009_v1.1_frozen"` 且 **`registry_eligible=false`**（核心不变量）；层 ID 与 `phaded_family_definitions.tsv` 一一对应。
- [x] **Step 1.2** 运行确认失败 → 实现 fail-closed 校验 → 生成权威表 → 测试通过。
- [x] **Step 1.3** 状态文档引用文献：8 超家族来自 28 条实验验证种子；38 家族来自聚类；DED 2009 后未更新；三处注释矛盾。并**显式写入**：本 Task 对 109,087 条候选为 0 行影响。
- [x] **Step 1.4** 报告口径统一并入本 Task 文档：所有新输出标注 `PhaDED (Knoll 2009, v1.1, frozen snapshot)` + 导入哈希；对无 SBD 家族（e-dPHAMCL、PhaZ7 型）显式声明 **PF06850 缺失是结构性预期而非阴性**；`dPHASCL1_like_candidate`（超家族候选标签）与 `type1_verified`（几何验证）必须分开表述。

> **执行结果（Task 1）：** `runs/20260917_phaded_classification_authority_01` — **completed**。权威表 46 行 = 8 superfamily（`functional_prior` / `registry_eligible=true`）+ 38 family（`sequence_clustering_2009` / `registry_eligible=false`），46/46 `source_version=Knoll2009_v1.1_frozen`，38 个 `layer_id` 与 `phaded_family_definitions.tsv` 双向零差集，`gate_profile_resolution` 46/46 resolved / pending 0；`candidate_universe_impact.rows_changed=0`。

### Task 2（P0，技术核心）：参考残基/判别域映射补全 + lid 去循环化

**目的：** 把 motif 层从"lipase box + oxyanion hole"面板升级为**文献可辩护的催化三联体 + 判别域面板**，并修正 lid 循环论证。

**Files:** 新建 `pipeline/scripts/complete_phaded_reference_residue_mapping.py` + 测试；run `runs/20260917_phaded_reference_residue_mapping_01/`；`inputs/literature_residue_reference.tsv`（版本化证据表）；`results/reference_mapping_manifest.tsv`

**Consumes:** DED 标注比对（38 个 alignment 文件，29 个含序列行、9 个 header-only，哈希已记录于 `reference_alignment_manifest.tsv`）；文献坐标：LtPHBase（PDB **8DAJ**：Ser121-His270-Asp197、oxyanion Cys40、四基序 [PMC9601781](https://europepmc.org/article/PMC/PMC9601781)）；P. funiculosum（PDB **2D80**，type 2 结构）；PhaZKT（Ser102-Asp221-His248 + lid 环 FNGIG aa34–38 / YYWQLF aa190–195）；PhaZGK13（无 lid 对照）[PMC12504323](https://europepmc.org/article/PMC/PMC12504323)

- [x] **Step 2.1（失败测试先行，核心不变量）**
  - His 模式必须在 **8DAJ 对齐的 LtPHBase** 参考上命中 His270、Asp 模式命中 Asp197（当前实现 His=0/723，**测试必失败**）；
  - Asp 模式必须在 type 2 结构（2D80 对应参考）上命中；
  - lid 检测器输入**不得含 `subtype_call` 或任何 superfamily/family 标签**（防循环论证的硬性测试）；且在 PhaZKT 对照命中 lid、在 PhaZGK13 对照不命中；
  - AHSMG 阳性对照（`2VTVA`/`AAK07742.1`）必须继续通过（回归保护）。
- [x] **Step 2.2** 运行确认失败 → 从 DED 标注比对 + 文献坐标构建逐超家族的残基映射表（缺失写 `pending`，不得推断）。
- [ ] **Step 2.3** 实现补全逻辑，使 723 条参考的 His/Asp/lid/linker 状态有据可查；SBD 用既有 Pfam **PF06850** 绑定补齐（`pending_tool_input` 转 `supported/not_detected`，不新写 motif）。—— **未完全达成**：His 90、Asp 71、lid 39、sbd 274（PF06850 绑定）已落地，但 **linker 仍 0 条可判定**（无 profile 绑定，41 条仅 `delimited_candidate_region_not_profiled`）；且 PF06850 的实际身份被后续证据推翻（见 DEFECT_D）。
- [x] **Step 2.4** 写 `literature_residue_reference.tsv`（accession、坐标来源、DOI/PMID、检索日期、哈希），并记录每列在参考面板上的命中率作为**实现校准报告**。
- [x] **Step 2.5** 状态文档：显式记录**缺陷 A 修正**（原 `lid_state` 由 subtype 先验派生，构成循环论证）与本任务后 motif 面板的可辩护范围。
- [x] **Step 2.6** 完整测试 + `compileall` + `git diff --check`。

> **执行结果（Task 2）：** `runs/20260917_phaded_reference_residue_mapping_01` — **partial**（自述 partial）。参考层（723 行 × 53 列，独立重算确认）：lipase_box 444、catalytic_ser_cys 444、his 90（48 `supported_anchor_column` + 42 `supported_literature_pattern`）、asp 71（48 + 23 + 2 `not_detected_at_mapped_column`）、oxyanion 256（255 模式 + 1 锚列）、sbd 274、lid 39、ahsmg 2、linker 0；其中真正由比对列**传递**得到的是 45/45/0（48/48/1 是含锚行自身文献坐标的状态计数）。遗留 blocker：Cys 型亲核体坐标不可得；intracellular nPHASCL（有/无 lipase box）、periplasmic、extracellular nPHASCL 四个超家族无带坐标的文献锚点；type 2 环状置换列保守度 0.0769；linker 与胞外 SBD 无检测器（**DEFECT_D**：PF06850 在本参考集中标记 intracellular nPHASCL without lipase box，274/276；在 284 条 extracellular dPHASCL type 1 与 71 条 type 2 中 0 命中）。

**边界：** 本 Task 不产生任何新的 family 判定；只完善参考映射与检测实现。

### Task 3（P1，纯分析，零新计算）：x₁ 疏水分层

**Files:** 新建 `pipeline/scripts/stratify_phaded_lipase_box_x1.py` + 测试；run `runs/20260917_phaded_x1_hydrophobicity_01/`

- [x] **Step 3.1（失败测试先行）** 断言 20 种氨基酸疏水判定与 Knoll 口径一致（疏水集 A,C,F,I,L,M,V,W,Y；G/P 单列非疏水）；非疏水 x₁ 的 `recommended_confidence_cap` 不得为 `high_candidate_only`；空值 → `not_tested`，不得计为阴性。
- [x] **Step 3.2** 实现 → 运行 → 核对总数自洽（疏水 + 非疏水 + 未测 = 109,087；非疏水预期 **563**）。
- [x] **Step 3.3** 产出 `results/confounder_candidates.tsv`（563 条）+ 状态文档（可发表口径：已测 99.27% 符合疏水 x₁ 特征；0.73% 为脂酶/酯酶混淆嫌疑，此前未识别）。**非疏水 x₁ 不是生物学阴性，不得删除候选。**

> **执行结果（Task 3）：** `runs/20260917_phaded_x1_hydrophobicity_01` — **completed**。独立重算 `x1_stratification.tsv` 109,087 行：hydrophobic 76,739 + non_hydrophobic 563（151 acidic + 178 basic + 80 glycine + 30 proline + 124 polar_uncharged）+ not_tested 31,785 = 109,087；`confounder_candidates.tsv` 恰 563 行；`negative_call_count=0`（非疏水只是混淆嫌疑，不是阴性）。

### Task 4（P1，纯分析/报告）：type 1/2 归属量化 + AHSMG 事实固化

**Files:** 新建 `pipeline/scripts/assign_phaded_catalytic_domain_type.py` + 测试；run `runs/20260917_phaded_catalytic_domain_type_01/`

- [x] **Step 4.1（失败测试先行）** 断言：仅当同时有 `oxyanion_hole_coordinates` 与 `lipase_box_coordinates` 才可判 type，否则 `undetermined_*`（fail-closed）；**type 1 = oxyanion hole 位于 lipase box 的 N 端（即坐标更小）**，type 2 反之；`not_detected_pattern` 不得获得 type 判定。
  > **勘误（2026-09-17，Phase 1 执行后更正）**：本 Step 初稿误写为「type 1 要求 oxyanion 位置 > lipase 位置」。该写法与三处独立证据矛盾：Knoll 2009 原文（`research/europepmc/comp_ft_knoll_phaded.txt` L17）明确 type 1 的 oxyanion hole 位于 lipase box 的 **N 端**；上游 `pipeline/scripts/audit_phaded_motifs.py` L177–180 把 `oxy.start() > lipase.start()` 对 type-1 标签序列标为 `conflict_relative_position`；且全部 2,362 条 `extracellular dPHASCL type 2` 标签行都是 oxyanion start > lipase start，按初稿字面规则会使 100% 的 type 2 证据失效。实现按文献口径执行（见 `docs/T141_20260917_phaded_catalytic_domain_type_status.md` §2.3）。
- [x] **Step 4.2** 实现 → 运行 → 量化：67,342 条 `dPHASCL1_like` 中真正通过 type 1 几何验证的比例（预期远小于 67,342，因全库仅 18,902 条有 oxyanion 证据）。
- [x] **Step 4.3** AHSMG 事实固化（v2 新增，替代 v1 Task 4）：在结果文档中写明——AHSMG 检测器**已在 2 条 PhaZ7 参考上通过阳性对照**（`AAK07742.1`、`2VTVA`），候选全库 109,087 条 0 命中 → 结论是"候选池未检出 PhaZ7 型"，**不得**写"该族不存在"。
- [x] **Step 4.4** 状态文档分开表述 `dPHASCL1_like_candidate`（标签）与 `type1_verified`（验证）。

> **执行结果（Task 4）：** `runs/20260917_phaded_catalytic_domain_type_01` — **completed**。独立重算 109,087 行：type1_verified 16,530 / type2_verified 2,362 / undetermined_no_oxyanion 89,822 / undetermined_position_conflict 373（合计 109,087）；`dPHASCL1_like` 67,342 中 type1_verified **16,511**（24.5181%），`dPHASCL2_like` 356 中 type2_verified 323（90.7303%）；AHSMG 参考层 2 supported / 721 not_detected，候选层 0/109,087（`conclusion_boundary` 明写"候选池未检出"而非"该族不存在"）；Step 4.1 计划文本与实现的 type 1 方向冲突已在计划内以勘误块更正（Phase 1 验证报告 P3 已标识，产物侧口径正确）。

### Task 5（P2，解锁校准）：按判别轴重设计 challenge panel

**Files:** 新建 `pipeline/scripts/design_phaded_discrimination_panel.py` + 测试；run `runs/20260917_phaded_discrimination_panel_01/`

四个判别轴（正端 / 负端）：

| 轴 | 正端 | 负端/challenge |
|---|---|---|
| 底物链长 | SCL 酶（PHB） | MCL 酶（已有 **Q84C08**，实验明确不水解 PHB） |
| 颗粒状态 | denatured-only（胞外） | native-active（PhaZ7 型，参考 `2VTVA`） |
| 定位 | 胞外（信号肽 + SBD） | 胞内（无信号肽 + 有 lid，依赖 Task 2 结果） |
| fold 邻域 | PF10503 阳性 | 脂酶/酯酶（同 fold、**极性 x₁**，对接 Task 3 的 563 条） |

- [x] **Step 5.1（失败测试先行）** 断言：每个 challenge 记录必须绑定一个明确判别轴且注明锚定依据；annotation-only 不得转写为 formal negative；未绑定既有 family 的保持 `unresolved`。
- [x] **Step 5.2** 实现（只读消费既有证据，如需新序列则另开 acquisition run 记录检索日期与哈希）→ 输出每轴可达性评估：当前能提供多少个**锚定明确轴**的 family-resolved negative。
- [x] **Step 5.3** 状态文档诚实报告距 gate 的差距。**不降低门槛。**

> **执行结果（Task 5）：** `runs/20260917_phaded_discrimination_panel_01` — **completed**。`challenge_panel_by_axis.tsv` 960 行 × 18 列（axis_1 18 / axis_2 24 / axis_3 60 / axis_4 858），960/960 `gate_credit_slot=none`、`formal_negative_eligible=false`；916 annotation_only + 4 experimental_negative + 40 experimental_positive；family 绑定 379 resolved / 581 unresolved。可达性：positive_anchor 291 + negative_anchor 4 + challenge 665 = 960；**37 个 reference-only profile 全部阻断**（`profiles_meeting_full_gate=0`），gap：positive 104 / held-out 37 / family-resolved negative 37 / challenge 5；`gate_relaxation_applied=false`。**证据缺口（V1）**：本 run 未留 `compileall` / `git diff --check` 日志。

### Task 6（P3，需授权 + 服务器）：`i-nPHASCL (no lipase box)` Cys 型判据方案

**为什么是它：** DED 最大家族（224 条，38%）；判据极硬（Cys-His-Asp + 无 `GxSxG` + Cys 前疏水 Val）；模型酶 *Cupriavidus necator*（原 *Ralstonia eutropha*）H16 PhaZ1（`CAJ92291.1`）功能已验证 → 是 4 个缺失 superfamily 中唯一能**合理绕开"≥3 条新实验阳性"门槛**的对象（判据来自功能已验证的模型酶，而非聚类）。其余三个保持 `reference_only`（periplasmic 全库仅 1 成员且生理角色未明；i-nPHASCL-lipase-box 仅 1 成员；native-SCL/PhaZ7-like 的 AHSMG 与 Bacillus 脂酶共享，单独不充分）。

- [x] **Step 6.1a（由并发 agent 完成，非本会话）** 建立 Cys 型判据定义表与**预注册证伪条件**，并在冻结证据上跑出候选层证据：`runs/20260917_phaded_inphascl_cys_decision_01/inputs/cys_criterion_definitions.tsv`、`inputs/falsification_preregistration.tsv`（**F1–F6 共 6 条**）、`results/cys_decision_manifest.json`（其中 **4 条数值化闸门全部 passed**；F5/F6 为结构性声明）。
- [x] **Step 6.1b** 写 decision record `docs/superpowers/plans/2026-09-17-phaded-inphascl-cys-decision-record.md`（三方案：保持 registry / 新增窄范围 Cys profile 并保留历史模型 / 版本化替换），预注册证伪条件（在脂酶/酯酶 challenge 集上出现未解释命中即回退）。—— **已完成：该文件存在**（落盘 mtime `2026-09-17T01:51:47` CST，早于 Phase 3 交接文档的落盘；**该文件是活文档：活文档的大小与 SHA-256 会因其自身更新而变化，请以磁盘实测为准；某一时点的快照值记录在 `runs/20260917_phaded_final_consistency_fix_01/results/final_consistency_fix_manifest.json` 与本轮新增的收口 run（`runs/20260917_phaded_closeout_01/results/closeout_manifest.json`）中**），§5.1 **推荐方案 A**，§5.2 **不需要服务器**。
  > **更正痕迹（`20260917_phaded_closeout_01`，2026-09-17）—— 残留自引用过期（P2）：** 本行原文把该文件的哈希与大小写死为「（**24,971 B / `5a5aa62611ad540f6cda82cfaf67aa35a0d7478a4705cc449aa1d7c6a17fca93`** / mtime `2026-09-17T01:51:47` CST，早于 Phase 3 交接文档的落盘）」。**该组数值是上一轮 fix run 动手之前的修复前快照，不是当前值**（上一轮为更正 P6/P8 已改写该文件，其改后快照见 `runs/20260917_phaded_final_consistency_fix_01/results/final_consistency_fix_manifest.json` 的 `changed_files[3].after`），故已按固定说明替换；**修复前快照事实予以保留：该文件在 01:51:47 已存在，且早于交接文档的落盘**。
- [ ] **Step 6.2** 提交评审，**等待明确授权**；此前不建 run、不改 config、不启动服务器任务。—— **未完成，且按治理规则不是"待办"**：方案 B/C 的 profile 训练路径被 `annotation_only` 禁训规则 + gate 的 `positive_count >= 3`（实测可训练记录 1 条）**判定不可行**，故"授权"不再是本 Task 的前置条件，方案 A 也不需要新 profile。
- [ ] **Step 6.3（获授权后）** 建 profile → dated deploy → T141 ≤40 线程打分 → candidate-only 输出 → 测试 + `compileall` + `git diff --check`。—— **未完成，已判定不可行**（同上）：**下一会话不应再尝试训练该 HMM**。

> **执行结果（Task 6）：** **partial — 进度由并发 agent 在 Phase 3 收尾期间推进，非本会话所做。**（**状态口径**：Step 6.1a/6.1b 已完成，Step 6.2/6.3 未完成但**已被治理规则判死**，故为 `partial` 而非 `completed`，也非"待授权"。）证据层已就绪：`runs/20260917_phaded_inphascl_cys_decision_01`（`status=completed_candidate_only`）产出 `results/cys_decision_manifest.json`（`cys_type_candidate_evidence_count = 29,949`）、`results/cys_candidate_evidence.tsv`（42,397,093 B）、`results/cys_reference_contrast.tsv`，并**预先注册证伪条件 F1–F6 共 6 条，其中 4 条数值化闸门全部 passed**（`pf06850_reference_specificity` observed_hits=0、`pf06850_reference_sensitivity` 274/276 = 0.992754、`extracellular_reference_no_false_positive` observed_hits=0、`lipase_box_miss_assignable` anchors_reproduced 4/4），**F5 / F6 为结构性声明**（不引入任何 fitted HMM/profile；候选命中只能是注释证据，由 `hmm_fitted=false` 与列级 `family_call_made=false` 强制）；另产出 `results/cys_typing_descriptor_discovery.json` 与 `results/cys_typing_descriptor_evidence.tsv`。**边界如实**：`hmm_fitted=false`、`new_family_call_made=false`、`subtype_call_rows_changed=0`、`gate_note` 明写校准 gate 未变且仍不满足。Step 6.1 的 decision record **已存在**（见 Step 6.1b 行），无任何授权证据、未建 profile、未 dated deploy、未启动服务器任务（`hmm_fitted=false` 佐证）—— 原因**不是**授权缺失，而是**方案 B/C 被硬规则禁止**（276 条参考中 275 条 `annotation_only` 禁训，唯一可严格训练记录 1 条 < gate 的 3 条）。**故本 Task 判 partial。**

> **更正痕迹（`20260917_phaded_final_consistency_fix_01`，2026-09-17）：** 本 Task 原文本（原文逐一保留如下引号内）有三处需更正：① **Step 6.1b 原文**「- [ ] **Step 6.1b** …—— **该文件经 glob 检查仍不存在**。」→ 该文件**存在**（见上；**该组数值 24,971 B / `5a5aa626…` / mtime 01:51:47 是修复前快照，不是当前值** —— 见本行上方 `20260917_phaded_closeout_01` 的更正痕迹），计划最后一次写入为 01:51:05（比记录早 42 秒），故起草时为真但**未被回填**；复选框已改为 `[x]`。② **原执行结果行**「**仍缺 Step 6.1 的 decision record**（… 经 glob 检查**不存在**）… **故本 Task 不判 completed。**」→ decision record 已存在；`partial` 的理由改为 Step 6.2/6.3 被治理规则判死。③ **原文**「预先注册 4 个证伪闸门且全部 passed」→ 实测预注册 **F1–F6 共 6 条**，其中 4 条被数值化评估（全部 passed），F5/F6 为结构性声明（依据：`inputs/falsification_preregistration.tsv` 6 个数据行；`cys_decision_manifest.json` 的 `falsification_gates` 4 个键，`outcome` 全为 `passed`）。核对证据见 `runs/20260917_phaded_final_consistency_fix_01/logs/verify_p1_p8.py` 的 `P1/P2/P8` 条目。

### Task 7（P2，依赖 Task 2）：全库 motif 补全层重跑（本地）

**Files:** 新建 run `runs/20260917_phaded_motif_completion_full_01/`

- [x] **Step 7.1** 用 Task 2 完成后的映射/检测器，对 109,087 条重算 His/Asp/lid/linker/SBD 状态，写入**新** dated run（不覆盖 0915 run）。
- [x] **Step 7.2** 输出新的 motif 层分布并与 0915 对比（His/Asp 从 0/25 → 参考校准后的实际检出率；lid 从循环赋值 → 独立检出）。
- [x] **Step 7.3** 状态文档：更新候选侧"催化三联体完整性"可辩护口径，明确 `partial_catalytic_pattern_panel` 与完整三联体证据的区别。
- [x] **Step 7.4** 若任一 motif 全库运行需要服务器（预计不需要：均为序列模式/profile 比对，本地即可），必须先建 dated deploy 并获授权。

> **执行结果（Task 7）：** `runs/20260917_phaded_motif_completion_full_01` — **partial**（自述 partial）。候选层 after（独立重算，`assignable + not_assignable = 109,087` 8 条判据全部成立）：lipase_box 77,302 / 31,785；his 7,317 supported（低于 0915 的 265 → **+7,052**），not_assignable 101,770（**93.29%**；更正：原文写 93.27%，实测 `101,770 / 109,087 = 93.2925%`，见 `20260917_phaded_final_consistency_fix_01` 报告 N7）；asp 1,290 support（**+0**，模式与 0915 字节相同）；oxyanion 18,939（**+37**，0915 的 373 条 `conflict_relative_position` 在新层并入 not_detected）；lid 989 supported + 1,921 partial（双环/单环）= `any_loop_detected` 2,910，106,177 not_detected（0915 的 lid 列 0 条独立检出，全部来自 subtype 先验）；ahsmg 0/109,087；linker 与 sbd 各 0 条可判定。**决定性事实**：`candidate_ded_alignment_members = 0`（109,087 条中无一条属于任何冻结 DED 比对），故 Task 2 的锚定列传递对候选层覆盖为 0，候选层只能是序列模式层。**口径更正（premise_corrections）**：Phase 1 简报引用的 His "0 命中" / Asp "约 25" 是参考层数字，候选层 0915 原值为 his 265 / asp 1,290。

---

## 4. 执行顺序与依赖

| 批次 | 任务 | 依赖 | 服务器 | 产出价值 |
|---|---|---|---|---|
| **P0** | Task 0（收尾） | 无 | 否 | 修正执行事实与计数口径 |
| **P0** | Task 1（权威口径） | 无 | 否 | 解决层级错配；0 行影响 |
| **P0** | Task 2（参考映射补全 + lid 去循环） | 无 | 否 | **motif 层从"部分"变"可辩护"；是后续一切判别的基础** |
| **P1** | Task 3（x₁ 分层） | 无 | 否 | 563 条混淆嫌疑（全新发现） |
| **P1** | Task 4（type 1/2 + AHSMG 固化） | 无（可并行） | 否 | 量化 67,342 标签证据强度 |
| **P2** | Task 5（challenge panel） | Task 2、3 | 否 | 解锁 family-resolved negative |
| **P2** | Task 7（全库补全重跑） | Task 2 | 否（预计） | 升级全库 motif 层 |
| **P3** | Task 6（Cys 型 record） | Task 1、2、5 + **授权** | **是** | 突破 4 个缺失 superfamily 的僵局 |

**建议起步（第一批，全本地、零风险）：Task 0 + Task 1 + Task 2。** 其中 Task 2 是技术核心：它不产出新分类，但决定 x₁ 之外所有判据是否可信。

---

## 5. 不做什么（prohibited）

- 不降低校准 gate（已证明无效：37 个 profile 的 negative/held-out 全为 0）。
- 不合并/删除 27 个零证据家族（已证明 0 行影响）。
- 不修改 `pipeline/config/formal_scan_models.tsv`、历史 HMM、`ePhaZ_curated_core`、历史 run。
- 不启动 GTDB 全量重扫；未获授权不启动任何服务器计算。
- 不把 AHSMG 零检出写成"该族不存在"（只写"以已通过阳性对照的模式未检出"）。
- 不用 subtype 派生的 `lid_state` 检验 subtype。
- 不覆盖 `runs/20260907_progress_report_01/`、`runs/20260915_*`、`runs/20260916_*`。
- 不删除失败重试、hold、challenge、negative 证据。

---

## 6. 验证命令（每个 Task 完成后）

```powershell
python -m unittest discover -s pipeline/tests -v
python -m compileall -q pipeline
git diff --check
```

服务器资源核对（执行前）：`& C:\Windows\System32\OpenSSH\ssh.exe -o BatchMode=yes <SERVER_USER>@<SERVER_HOST> 'uptime; nproc; free -g | head -2'`

---

## 7. 交付物清单（预期）

| 产物 | 路径 |
|---|---|
| 分类权威表 | `runs/20260917_phaded_classification_authority_01/inputs/phaded_classification_authority.tsv` |
| 文献残基参考表 | `runs/20260917_phaded_reference_residue_mapping_01/inputs/literature_residue_reference.tsv` |
| 参考映射校准报告 | `runs/20260917_phaded_reference_residue_mapping_01/results/reference_mapping_manifest.tsv` |
| x₁ 混淆清单 | `runs/20260917_phaded_x1_hydrophobicity_01/results/confounder_candidates.tsv` |
| type 归属 | `runs/20260917_phaded_catalytic_domain_type_01/results/catalytic_domain_type.tsv` |
| challenge panel | `runs/20260917_phaded_discrimination_panel_01/results/challenge_panel_by_axis.tsv` |
| Cys 型 decision record | `docs/superpowers/plans/2026-09-17-phaded-inphascl-cys-decision-record.md` |
| 收尾 amendment | `runs/20260917_phaded_housekeeping_amendments_01/` |
| 全库补全层 | `runs/20260917_phaded_motif_completion_full_01/` |
| 各 Task 状态文档 | `docs/T141_20260917_*_status.md` |

---

## 8. 文献依据

- Knoll M, et al. **The PHA Depolymerase Engineering Database.** BMC Bioinformatics 2009;10:89. [PMC2666664](https://europepmc.org/article/PMC/PMC2666664)
- **Biodegradation of PHA: current state and future prospects.** Front Microbiol 2025. [PMC11893044](https://europepmc.org/article/PMC/PMC11893044)
- **DED — Database Commons**（Last update 2009, v1.1）[链接](https://ngdc.cncb.ac.cn/databasecommons/database/id/3446)
- **α/β Hydrolases: Toward Unraveling Entangled Classification.** Proteins 2025. [PMC11878206](https://europepmc.org/article/PMC/PMC11878206)
- **Lid in mclPHA intracellular depolymerase.** Appl Microbiol Biotechnol 2025. [PMC12504323](https://europepmc.org/article/PMC/PMC12504323)
- **LtPHBase.** Protein Sci 2022. [PMC9601781](https://europepmc.org/article/PMC/PMC9601781)（PDB 8DAJ 坐标）
- **Brigham et al. R. eutropha PhaZ.** PLoS ONE 2012. [PMC3430594](https://europepmc.org/article/PMC/PMC3430594)
- 项目内部：`research/europepmc/comp_classification.md`、`comp_ft_knoll_phaded.txt`、`comp_ft_thermophile_phb.txt`

---

## 9. 修订记录（v1 → v2）

1. **Task 结构重组**：8 → 7 个任务；"报告口径统一"并入 Task 1 与 Task 4；"计数不一致"并入 Task 0。
2. **AHSMG 任务降级**：v2 核查发现阳性对照（`2VTVA`/`AAK07742.1`）已在参考面板通过 → v1 的"实现审计"改为 Task 4 的"事实固化"步骤，不再需要联网取序列。
3. **新增 Task 2（技术核心）**：v2 核查发现 His 模式 723/723 从未触发、Asp 仅对 type 1 生效、SBD/linker 从未运行、lid 循环赋值——这些 `pending` 是真实缺口，且是所有判据可信度的前置条件。
4. **优先级调整**：Task 2 与 Task 1 并列 P0；x₁ 分层（v1 Task 2）降为 P1 纯分析。
5. **新增 Task 7**：Task 2 完成后对全库 109,087 条重跑补全 motif 层（新 dated run，不覆盖 0915）。
6. 事实基线补入 1.2（文献）与 1.3（参考层/候选层对照表），并记录本次状态复核结论（09-17 00:52 无新增变更）。

---

*本计划为 candidate-only 研究计划；所有产出表示候选同源或功能潜力，不等同于已验证 PHB/PHA 降解表型。*
