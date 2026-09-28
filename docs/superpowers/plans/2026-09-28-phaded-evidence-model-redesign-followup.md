# PhaDED 证据模型重构 —— 未完成工作后续计划（2026-09-28）

> **主计划（已完成的 Task 1–13）**：`docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign.md`
> **设计（不变）**：`docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md`
> **现行状态/权威表**：`docs/T141_20260928_phaded_evidence_model_redesign_status.md`
> **B 类运行的授权申请包**：`docs/T141_20260928_phaded_evidence_run_authorization_packets.md`
> **本计划起点**：HEAD `0fa4cd1`（3 个本地 commit，**未 push**）；全库 **1,243 项测试 OK**、`compileall` 0、`git diff --check` 0；v1 冻结产物 SHA-256 与开工时逐位相同。
> **本文档性质**：candidate-only 计划件。所有 HMM/profile/domain/motif/SignalP/结构/系统发育/分类与生长速率证据只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。

---

## 0. 未完成工作的三类划分（先划边界，再排任务）

主计划 Task 1–13 已交付**代码、测试、配置、文档与 CI 门**，但刻意停在授权 checkpoint 之前。剩余工作按**是否需要新的计算/运行授权**分三类：

| 类 | 含义 | 可否立即开工 |
|---|---|---|
| **A 类：本地代码与只读分析** | 只读冻结产物、或只改 `pipeline/` 下源码与测试；**不新建** `runs/<run_id>/` 或 `deploy/<run_id>/`，不调用任何计算工具 | ✅ 只需你确认开工（无需新的计算授权） |
| **B 类：证据采集与重算运行** | 需新建 dated `runs/<run_id>/` + `deploy/<run_id>/`，或调用 SignalP / HMMER / Foldseek / 结构预测 / gRodon / GTDB 全库 | ❌ **每个计算类别单独授权**，一类批准不等于另一类 |
| **C 类：发布链** | 研究快照入版本控制、push / tag / GitHub Release / Zenodo DOI | ❌ 需明确授权 |

**总览表**（任务明细见 §2–§4；状态列为 2026-09-28 执行后的实测状态）：

> **执行已完成并推送（2026-09-28）**：F1–F9、F13–F17 已执行；8 个 commit 已 push 到 `origin/main`（`3117287..f8b0cf8`），tag `phaded-evidence-model-v2-20260928` 已推送。最终门禁：**1759 项测试 OK**（skipped 1，环境性）、`compileall` 0、`git diff --check` 0、`runs/results/deploy` 零改动、v1 三件权威件 SHA-256 逐位不变。
> **本表的权威记录是 `docs/T141_20260917_project_handoff.md` §12.11**（含四个科学结果、两处被推翻的计划前提、以及仍未完成项）。
> **F10（P3）、F11（P4）、F12（P5）仍未执行**——它们各自卡在 §5 的待决策项上（#7 是否新建 MAFFT 比对、#8 竞争面板 PDB 闭集、#9 抽样规模与种子）。这些是**科学设计决策**，编造它们会违反项目规则，故等待操作者裁决。
> **GitHub Release 与 Zenodo DOI 未创建**——`gh` 未登录、无 Zenodo API token；push 与 tag 已用 Git Credential Manager 的既有凭据完成。

| # | 类 | 任务 | 状态 | 授权 |
|---|---|---|---|---|
| F1 | A | 616 条候选流转残留：归因 + 处置规则提案 | ✅ **完成**：472 `profile_unassigned_no_discovery_claim` + 140 `discovery_claim_not_unique` + 4 `score_tier_not_upgraded`，`Σ=616`、`unexplained=0`（双向完整性等式证明）；全部映射 `function_unresolved` | 否 |
| F2 | A | gRodon 66 条差额分桶 | ✅ **完成（前提被推翻）**：66 = 23 实测 + 43 不可解；原「本地可闭合」前提错误，见 §F2 留痕 | 否 |
| F3 | A | 三个下游脚本改读 `model_layer` | ✅ **完成**：3 处判据点（计划只列了 2 处）；缺列 manifest 回退 `reference_query_only` 并声明 | 否 |
| F4 | A | `formal_scan13_tier_processing.sh` 绑定全库 `-Z` | ✅ **完成**：fail-closed（缺 Z 即 exit 1）；manifest 记尺度字段；两处缺口测试已翻正 | 否 |
| F5 | A | `prepare_phaded_structure_review.py` 解硬编码与 `pending_tool` 误述 | ✅ **完成**：19 字面量消失；工具可用性改为**实测**（区分 `unavailable` 与 `pending_tool`），并区分「可用」与「已执行预测」 | 否 |
| F6 | A | 发现层/deferred 层下游误用审计 | ⏳ 进行中 | 否 |
| F7 | B | scan-13 `-Z` 尺度对账 run（**无 HMMER**） | ⏳ **Step A 运行中**（服务器 40 jobs，100 分片 / 267 GB）；deploy + run + input_contract 已绑定；`hits_all.tsv` SHA-256 与冻结 manifest 逐位一致 | ✅ 已授权 |
| F8 | B | P1：723 参考层 SignalP 补测 | ⬜ 未开始 | ✅ 已授权 |
| F9 | B | P2：11,610 条 SP-less type-1 配件域核查 | ⬜ 未开始 | ✅ 已授权 |
| F10 | B | P3：29,974 条 Cys 参考锚定残基映射 | ⬜ 未开始 | ✅ 已授权 |
| F11 | B | P4：987 条 nPHAMCL 结构竞争面板 | ⬜ 未开始（F5 已完成，前置解除） | ✅ 已授权 |
| F12 | B | P5：with-lipase 分层抽样结构验证 | ⬜ 未开始 | ✅ 已授权 |
| F13 | B | 723 行台账证据分级 curation | ✅ **完成**：E3=1 / E2=26 / E1=0 / A=3（30 条阳性）；唯一独立性分组计数 **24**；`families_meeting_minimum = [hfam_52(6), hfam_70(7)]`；strict 校验通过；冻结台账 SHA 未变 | ✅ 已授权 |
| F14 | B | v2 候选目录真实重算 + 对账 + v1/v2 影响报告 | ⏳ 进行中 | ✅ 已授权 |
| F15 | B | `hfam_70` promotion 与 finalize | ⏳ 进行中（操作者已授权提升；hfam_52 按设计刻意不提升） | ✅ 已授权 |
| F16 | C | 研究快照入版本控制 | ⬜ 未开始 | ✅ 已授权 |
| F17 | C | push / tag / GitHub Release / Zenodo DOI | ⬜ 未开始（远端 `github.com/anzhiw2-gif/PHB_gtdb-ds`，`gh` 2.95.0 可用） | ✅ 已授权 |

---

## 1. B 类共享执行契约（每个 B 类任务都必须遵守，不再逐条重复）

1. **新 run、新 deploy**：`runs/<run_id>/{logs,inputs,results}` + `input_contract.json`；服务器只执行 dated `deploy/<run_id>/` 中已绑定源码。`run_id` 可审计、无路径遍历（建议 `^[0-9]{8}_[a-z0-9_]+$`）。
2. **动态线程**：启动前紧邻实测 `min(40, nproc − 1 分钟 loadavg − 10)`，写进资源记录；**40 是上限不是启动值**；结果低于 1 则中止。
3. **溯源**：GTDB taxonomy/metadata/tree、HMM、比对、源码、环境的路径 + 版本 + 大小 + SHA-256 全部落 manifest；缺失写 `pending`，**禁止伪造哈希**。
4. **candidate-only**：任何输出不得删除候选；结果只能改 `evidence_flags` 或 `primary_disposition`。
5. **E 值尺度**：任何 HMMER 分片搜索必须传**同一个全库 `-Z`**，并在 manifest 记 `database_size_Z`、逐片序列数、命令模板；`06_screen.sh` 缺 `--database-size-z` 时 fail-closed。
6. **完成判据**：相关测试 + `python -m compileall -q pipeline/scripts pipeline/tests` + `git diff --check` 全过；交接文档留痕。

---

## 2. A 类：可立即执行的本地任务

### F1：616 条候选流转残留 —— 归因与处置规则提案

**背景（实测）**：109,087 − 36,611 − 71,860 = **616**。这 616 条既不在池内高可信层、也不在 hold 层；与池外 strong 层交集为 **0**；换 demote 后口径（36,559）残差为 **668 = 616 + 52**，故它们**不是**那 52 条 SignalP-export demote。`reconcile_phaded_candidate_flow.check_hold_residual` 目前持续报 `mismatch`，这是**有意的响亮失败**。

**文件**：
- 只读：`runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv`、`runs/20260919_phaded_nucleophile_unify_01/results/{high_confidence_candidates.tsv,hold_candidates.tsv,filter_summary.json}`、`runs/20260919_phaded_pool_external_evidence_01/`
- 产出（新 run）：`runs/<run_id>/results/residual_616.tsv`、`residual_616_attribution.tsv`、`residual_616_disposition_proposal.md`、manifest

**步骤**：
- [ ] Step 1 用 `check_hold_residual(...)` 导出 616 条 accession 名单（唯一来源，不做集合相减的近似）。
- [ ] Step 2 与 `protein_layers.tsv` 对照：逐条取 `layer`、`major_superfamily`、`score_tier`、E 值档、是否 registry 模型命中。
- [ ] Step 3 归因：判定 616 条是在**哪一阶段**落到两表之外（tier 阈值 / registry 家族范围裁剪 / 池外分桶边界 / 表构建 bug），每条给一个机器可读 `attribution` 值（如 `tier_below_threshold`、`outside_registry_family_range`、`unexplained`）。
- [ ] Step 4 提出**处置规则**（不是逐条手改）：为每个归因给出 `primary_disposition` 映射，写进 `disposition_proposal.md`；无法归因的必须列 `unexplained` 并计数，**不得默认塞进某个桶**。
- [ ] Step 5 把 proposal 作为 F14 的输入之一；`check_hold_residual` 的期望值保持「未解释只数」直到 F14 落地。

**验收**：`residual_616.tsv` 恰 616 行且 accession 唯一；`Σ attribution 计数 = 616`；proposal 中每个归因都有对应处置与理由；输出只新增 run，**不改** v1 冻结表。

**授权**：不需要新计算授权（只读分析 + 新 run 目录）——但**新建 run 目录**本身按 AGENTS.md 属「新运行」，需你在开工时确认；若你希望连 run 目录都不建，可退化为「输出到临时目录 + 文档记录」。

### F2：gRodon 66 条差额 —— 按 accession 集合分桶闭合

> **⚠️ 本节前提已被 F2 执行结果推翻（2026-09-28，留痕不掩盖）**
>
> 本节原写「已实测本地可闭合」，并假定 `4,507 = 30,623 − 26,116` 是可 accession 的集合运算。**实测证明该假定错误**：
>
> 1. `26,116`（冻结 deploy 的 `excluded_pos`）**不指向任何 accession**——它是 `build_grodon_manifest.py` 在 6 处累加的**标量**，同时混合了两种不同原因：整个属无对照臂（L177/185/187/192）**与**阳性被 1:1 截断挤出（L195）。
> 2. 本节原建议的运算方向也错：`degrader_genomes.tsv (30,623) − pool_genomes_exclusion.txt (62,666) = 905`，**不是 4,507**；且 30,623 条中有 **29,718 条本身就在该排除文件里**（该文件是对照臂的候选池账本，不是"无对照"账本）。
> 3. **正确闭合式**（全部实测）：
>    ```text
>    4,507（文档口径）= 4,441（manifest 阳性，实测）
>                     +    23（预测失败的阳性，实测，≥ accession 级集合）
>                     +    43（残差：冻结证据中无任何输入能为这些记录命名 → pending）
>    ```
>    接手文档原记的「`4,441 + 23 failed ≠ 4,507` 是悖论」**不是悖论**：那 54 条失败记录是**先被选进 manifest、随后预测失败**，故 `failed_positive ⊂ manifest_positive`；朴素相加把 23 条重复计入。
> 4. 真实集合差是 `len(degrader_input − manifest_positive) = 26,182`（**不是 66**）；`identity_check.passes=false` 是**故意的**——没有任何桶可以认领这 26,182 条，执行方拒绝伪造 eligible 名单去凑一个"通过"。
>
> **已产出的 F2 结果**：新 run `runs/20260928_phaded_grodon_66_reconciliation_01/`（`results/grodon_66_buckets.tsv`、`grodon_66_reconciliation.json`），新脚本 `pipeline/scripts/reconcile_phaded_grodon_manifest_delta.py` + 44 项测试。
>
> **永久闭合的最省路径（尚未执行，需单独授权）**：给 `build_grodon_manifest.py` 加 `--emit-excluded-ledger`，把被丢弃的阳性 accession 连同原因逐条写出。这是**纯选择回放**，不需要 HMMER 或 gRodon；该单一产物即可让 43 条残差变成可 accession，从而把文档口径的 66 彻底闭合为 `23 实测 + 43 有据`。

**背景（本会话实测，全部来自冻结产物）**：

| 量 | 值 | 来源 |
|---|---:|---|
| `degrader_input_genomes` | 30,623 | `results/manifest_stats.json` |
| `degrader_genomes_excluded_no_control` | 26,116（**标量，非 accession 集合**） | 同上 |
| 文档口径 eligible（30,623 − 26,116） | 4,507 | 与接手文档 §8.1 一致 ✅ |
| `manifest_positive` | 4,441 | 同上 |
| **差额** | **66 = 23 实测 + 43 不可解** | F2 实测 |
| `manifest_rows` | 8,882（4,441 + 4,441） | 同上 |
| 预测 ok → 平衡 | 8,828 → 8,776（4,388 对） | 接手文档 §8.1 |
| failed | 54 = **23 阳性 + 31 对照**（表内有 `phaZ_status`） | `results/grodon_failed_genomes_newdeg40k.tsv` |
| `missing_fasta` | **0** | `manifest_stats.json` |

**关键**：`missing_fasta = 0` 意味着 Task 12 合成测试里「20 条缺 FASTA」那一桶在**真实数据中为 0**；真实分桶必须重新推导，不得照抄合成比例。

**文件**：
- 只读：`runs/20260920_phaded_grodon_growth_01/inputs/{degrader_genomes.tsv,pool_genomes_exclusion.txt}`、`results/{grodon_growth_manifest_newdeg40k.tsv,grodon_failed_genomes_newdeg40k.tsv,grodon_growth_balance_audit_newdeg40k.tsv,grodon_growth_balanced_genus_selected_counts_newdeg40k.tsv,manifest_stats.json}`
- 已产出（新 run）：`runs/20260928_phaded_grodon_66_reconciliation_01/{results/grodon_66_buckets.tsv,results/grodon_66_reconciliation.json,inputs/,logs/,input_contract.json}`

**步骤（实际执行）**：
- [x] Step 1 以 accession 集合推导可闭合部分；**拒绝**了本节原建议的错误运算（已在上方留痕）。
- [x] Step 2 用 `grodon_reanalysis_v2.explain_manifest_difference` 分桶；实测桶：`prediction_failed_positive_arm = 23`、`prediction_failed_control_arm = 31`、`manifest_positive_arm = 4,441`、`absent_from_manifest = 26,182`、`residue_unresolvable_from_the_frozen_evidence = 43`（0 accession，pending）。
- [x] Step 3 拒绝以占位 accession 把 `identity_check` 凑成 true（该"能通过"的情形已被实测并**标注为陷阱**）。
- [x] Step 4 结论：`check_grodon_manifest_delta` 在文档口径上持续 `mismatch` 是**真阳性**，不是 bug；要变成可复现 `ok`，唯一缺的输入是**per-accession eligible 名单**（等价于被命名的 26,116 条）。

**验收**：实测桶计数之和可解释 `4,507 = 4,441 + 23 + 43`；`grodon_66_buckets.tsv` 恰好覆盖全部 30,623 条输入且 accession 唯一；43 条残差以 `pending` 明确列出并写明缺失输入（GTDB `bac120_taxonomy_r232.tsv` + 冻结 deploy 的 seed-42 shuffle 状态，二者本地均不存在）。

**授权**：不需要计算授权（纯只读 + 本地集合运算）。

### F3：三个下游脚本由 `model_status` 改读 `model_layer`

**背景**：Task 5 把分层语义改为 `model_layer`，但 `model_status` 词表为向后兼容未改，于是**发现层 HMM 在旧判据下也显示为 `trained`**，下游可能把「发现层」当「判别层」使用——这正是 AGENTS.md 明令禁止的。

**实测命中点**：
- `pipeline/scripts/run_phaded_candidate_hmmsearch.py:35`
- `pipeline/scripts/build_phaded_with_lipase_deferred_tier.py:18`
- `pipeline/scripts/build_phaded_family_discovery_sets.py`（需在实现时逐处确认判据）

**步骤**：
- [ ] Step 1 先写失败测试：构造 `model_layer="discovery_hmm_uncalibrated"` 但 `model_status="trained"` 的 profile manifest，断言这些脚本**不得**把它当判别层（不得用于筛选/降级/家族判定）。
- [ ] Step 2 改为读 `model_layer`；对 `reference_query_only` 必须拒绝评分（无 HMM）。
- [ ] Step 3 保持 `model_status` 列本身不变（历史 manifest 兼容），只在**判据**处切换，并在 docstring 说明。
- [ ] Step 4 验证：相关测试 + `compileall` + `git diff --check`；确认冻结产物未被读写。

**授权**：不需要计算授权（本地代码 + 合成 fixture）。

### F4：`formal_scan13_tier_processing.sh` 绑定同一个全库 `-Z`

**背景**：Task 9 已让 `06_screen.sh` fail-closed 要求 `--database-size-z`，但 **tier 处理入口自身仍未绑定尺度**——它的 `08c_tier_rescore.py` 会重跑 HMMER。当前只由测试钉住缺口（`test_tier_entrypoint_does_not_yet_pin_a_database_size`、`test_tier_manifest_schema_still_lacks_scale_fields`）。

**文件**：`pipeline/scripts/formal_scan13_tier_processing.sh`、`pipeline/scripts/08c_tier_rescore.py`、`pipeline/tests/test_hmmer_database_size.py`（把两个「缺口断言」改为正向断言）

**步骤**：
- [ ] Step 1 入口接受 `--database-size-z` + `--database-size-basis`，正整数校验，缺失即 fail-closed（与 `06_screen.sh` 同规）。
- [ ] Step 2 `08c_tier_rescore.py` 的所有 `hmmsearch` 调用经 `pipeline/scripts/hmmer_command.py::build_hmmsearch_command` 构造，传入同一 Z；禁止 `--domZ`。
- [ ] Step 3 `tier_processing_manifest.json` 增加 `database_size_Z`、`database_size_basis`、逐片序列数、命令模板。
- [ ] Step 4 把两个缺口测试改成「必须存在 `--database-size-z`」的正向断言。
- [ ] Step 5 验证：`test_hmmer_database_size`、`test_formal_scan13_tier_run` + `compileall` + `git diff --check`。

**授权**：不需要计算授权。**注意**：本任务不重跑 tier 处理，只改代码与测试。

### F5：`prepare_phaded_structure_review.py` 解除硬编码（P4 前置）

**背景（实测）**：`pipeline/scripts/prepare_phaded_structure_review.py:54-55` 硬编码 `if len(rows) != 19: raise ...`；L91 的 manifest 声明 `structure_predictor.status = "pending_tool"`（理由「无可用预测/比对可执行文件」），与实测「Foldseek 可用（自测 TM 1.00）」**直接冲突**。后果：P4 的 987 条输入会 **fail-closed 中止**。

**步骤**：
- [ ] Step 1 失败测试：目标集为 987 条（或任意 N）时必须正常构建；`len(rows)` 与声明期望不符时报错，但期望值由 `--expected-candidates` / 输入契约给出，**不再写死 19**。
- [ ] Step 2 工具可用性改为**实测驱动**：运行时探测 Foldseek/结构预测可执行文件，可用则记版本 + 路径 + SHA-256 并置 `status=available`，不可用才 `pending_tool`；不得预写结论。
- [ ] Step 3 manifest 的 `candidates_checked` 与 `reason` 由实测填充；`pending` 只用于真正未探测项。
- [ ] Step 4 验证 `pipeline/tests/test_prepare_phaded_structure_review.py` + `compileall` + `git diff --check`。

**授权**：不需要计算授权（只改代码与测试；探测逻辑在**运行时**才执行，届时属 P4 授权范围）。

### F6：发现层/deferred 层下游误用审计

**背景**：AGENTS.md 规定发现层只召回、不进 `formal_scan_models.tsv`、不产生 family 判定；with-lipase deferred 层禁删、禁入计数。除 F3 的两个命中点外，可能还有调用点把两者当判别层/当主计数。

**步骤**：
- [ ] Step 1 只读扫描 `pipeline/scripts/*.py`：所有读取 `profile_manifest`、`formal_scan_models.tsv`、deferred TSV 的脚本，逐个判定其用途（召回 vs 判别 vs 计数）。
- [ ] Step 2 输出审计表 `discovery_layer_usage_audit.tsv`：脚本 / 行号 / 用途 / 是否违规 / 依据。
- [ ] Step 3 违规项写入后续任务清单或就地修复（若属 A 类）；**不得**删除任何 deferred 层文件。
- [ ] Step 4 把审计表纳入 `pipeline/tests/` 的一条静态断言测试（防止未来回归）。

**授权**：不需要计算授权。

---

## 3. B 类：需要授权的证据采集与重算

> 每个任务的**完整口径**（决策问题 / 目标群体 / 抽样依据 / 工具版本 / 假阳性竞争者 / 输出 schema / 结果如何只改证据标记 / 线程资源 / 预估规模 / 不属于本次授权的动作）已在 `docs/T141_20260928_phaded_evidence_run_authorization_packets.md` 中成文；本节只列**开工前必须补齐的前置**与**验收判据**。

### F7：scan-13 `-Z` 尺度对账 run（**无 HMMER**，优先级最高的 B 类任务）

**依据（Task 9 只读审计，实测）**：scan-13 执行的 hmmsearch **没有 `-Z`、没有 `--domZ`**（证据：冻结 run 自记的 driver SHA-256 与本地 `formal_frozen_screen_parallel.sh` 逐位相同，其 L115 即命令模板；独立复核见 `runs/20260917_phaded_task_completion_check_01/01_recall/README.md` §8）。因此**每个分片的 E 值用的是该分片自己的序列数**（shard_0001 = 4,694,121 vs shard_0051 = 2,923,820，差 1.61 倍）。E 随 Z 严格线性已被两次独立实验证实，故**不需要重跑任何 HMMER**。

**步骤**：
- [ ] Step A 对冻结的 100 个分片跑 `grep -c '^>'`，得逐片 `Z_shard` 与**实测全库 Z**（闭合长期 `pending` 的 ~2.92×10⁸ 近似值）。
- [ ] Step B 用既有 `pipeline/scripts/parse_phaded_cys_targeted_recall.py::rescale_evalue`（`E_full = E_shard × Z_total / Z_shard`）重标 6,743,198 行 → 新表 + 新 manifest（含 `database_size_Z`、逐片计数、命令模板）。
- [ ] Step C 在新列上重新施加 `-E 1e-5`，报告与冻结 6,743,198 行的**命中数差**；**新旧 E 值永不混用**，冻结表一字不改。

**验收**：逐片计数之和 = 实测全库 Z；重标后的表与其 manifest 自洽；差异报告给出方向与量级；冻结 `hits_all.tsv` 未改。

**授权**：需 ✅（新建 run + dated deploy + 服务器只读访问冻结 scan-13 目录）。**不需要** GTDB 重扫。

### F8–F12：P1–P5 五个定向证据运行

| 任务 | 目标 | 开工前必须补齐的前置 |
|---|---|---|
| **F8 P1** | 723 参考层 SignalP 补测，按实验定位/功能证据/继承历史标签分层，标定「已知胞外参考的预测信号肽检出率」 | 无（目标集可精确枚举） |
| **F9 P2** | 11,610 条 SP-less type-1 的配件域核查（FN3/Ig-like/TSP3/CHB，独立于信号肽） | taxonomy join 源与版本（§5 #6）；已确认目标集可由**单键**复现 |
| **F10 P3** | 29,974 条 Cys 的参考锚定残基映射，显式区分 substitution / truncation / uncertain | 是否新建 MAFFT 比对（若需 → 触发 Checkpoint 3，**另一授权类别**） |
| **F11 P4** | 987 条 nPHAMCL 结构**竞争面板**（非单一 8YNV） | **F5 完成** + 竞争面板成员 PDB 闭集（§5 #8） |
| **F12 P5** | with-lipase 1,206,655 条的**分层代表抽样**结构验证 | 抽样规模、种子、六个分层键中四个的 join 源（§5 #9） |

**共同验收判据**（每包单独）：决策问题有明确答案或明确「不可判定」；目标集行数与来源可复现；工具版本与路径落 manifest；**输出只改 `evidence_flags` / `primary_disposition`，不删除任何候选**；deferred 层计数一律不进主候选。

**授权**：SignalP、HMM 构建、Foldseek/结构预测、GTDB 搜索是**四类分别的批准**；批准其中一类不自动授权其他类。

### F13：723 行真实台账证据分级 curation

**背景**：Task 4 已交付校验器与 15 列新字段，但**真实台账仍为未 curation 状态**（compat 模式报告 15 列 `pending`、`uncurated_experimental_positive_rows = 30`、`functional_positive_count = 0`）。这是文献/证据采集工作，不是纯工程。

**步骤**：
- [ ] Step 1 逐条为 30 条 `experimental_positive` 与 693 条 `annotation_only` 填 `experimental_evidence_grade`（E3/E2/E1/A），**每条必须绑定** `study_id` / `experiment_unit_id` / 底物 / 体系 / 对照 / `independence_group`。
- [ ] Step 2 用 `summarize_independence` 报告「按唯一 `independence_group` 计数」的功能阳性数（预期显著低于 30，因为存在同基因重复报道）。
- [ ] Step 3 走 `strict` 模式校验；未取得的值写 `pending`，**不得推断**。
- [ ] Step 4 产出新的台账输入（新 run），**不覆盖**冻结台账。

**授权**：✅ 需（新建 run + 文献获取）。**产出直接决定 F14 与 F15 的门槛是否可满足**。

### F14：v2 候选目录真实重算 + 对账 + v1/v2 影响报告

**步骤**：
- [ ] Step 1 用 `build_phaded_candidate_catalog_v2.py` 生成 v2 目录（每条蛋白唯一 `primary_disposition` + 多 `evidence_flags`）。
- [ ] Step 2 用 `reconcile_phaded_candidate_flow.py` 闭合流转（含 F1 的 616 处置），要求 `unaccounted = 0`、`multiply_disposed = 0`。
- [ ] Step 3 用 `compare_phaded_catalog_versions.py` 出三场景影响报告，**76:24 只能作为「具名规则下的候选构成」表述**（护栏已内置）。
- [ ] Step 4 冻结 v1 表保持只读；v2 全部写新 run。

**授权**：✅ 需（计划 Authorization Checkpoint 5）。

### F15：`hfam_70` promotion 决策与 finalize

**背景**：六腿 gate 已通过，但按新规范状态为 **`candidate_gate_passed_not_promoted`**；`calibrated_candidate_model` **只能**由显式授权动作产生。

**步骤**：
- [ ] Step 1 你裁决是否提升（§5 #1）；若提升，运行 `finalize_phaded_reference_panel_calibration.py --promote-calibrated-model`。
- [ ] Step 2 记录：仍未进 `formal_scan_models.tsv`、仍未重扫（全库召回已在 `20260919_phaded_trained_profile_recall_01` 完成，99.94% 召回）。
- [ ] Step 3 明确边界：即使提升，也不代表任何单条命中是实验表型。

**授权**：✅ 需（finalize 是显式授权动作）。

---

## 4. C 类：发布链

### F16：研究快照入版本控制（Checkpoint 1）

范围：`pipeline/scripts` 下约 95 个历史脚本、`pipeline/tests` 下约 100+ 历史测试、`docs` 下约 85 份历史文档。**禁止 `git add -A`**，须逐目录审阅后分批提交，并在交接文档记录每个 commit。**建议先在临时克隆上试提交**，确认无身份信息泄露（`test_public_repo_safety` 已是门禁）。

### F17：push / tag / GitHub Release / Zenodo DOI

**前置**：F16 完成 + 发布前四路核对（本地 / GitHub / 服务器 deploy / run manifest 权威一致）+ 交接文档齐备。**push 必须操作者明确授权**，且授权是针对**具体 commit 范围**的，不是泛化的「以后都可以推」。

---

## 5. 开工前待决策清单（9 项，按阻塞面排序）

| # | 决策项 | 阻塞 | 建议默认 |
|---|---|---|---|
| 1 | `hfam_70` 是否提升为 `calibrated_candidate_model` | F15 | 只做序列层验证、暂不提升（gate 通过 ≠ 功能已校准） |
| 2 | 616 条的处置规则（F1 会给提案） | F14 | 按归因批量映射，`unexplained` 保持显式计数 |
| 3 | 4,301 条 Cys 漏斗丢弃是否恢复（接手文档 §9 #3） | F10 | 先出归因再决定，不默认恢复 |
| 4 | 池外 mid/weak 4,926 条是否升级补算 | 无（当前决策：不升级） | 维持不升级 |
| 5 | compat 模式下 `E1+direct` / `A+非annotation` 是否升级为致命 | F13 | 维持「报告不致命」，strict 模式致命 |
| 6 | P2 分层所需的 taxonomy join 源与版本 | F9 | 用 `input_contract.json` 已绑定的 GTDB 版本 |
| 7 | P3 是否新建 MAFFT 比对 | F10 | 优先复用既有比对；新建则另走 Checkpoint 3 |
| 8 | P4 竞争面板的 PDB 成员闭集 | F11 | 由你定；不得只留 8YNV 一条 |
| 9 | P5 抽样规模 / 种子 / 分层键 join 源 | F12 | 规模与种子须在申请包中预注册后再抽样 |

另有两项低风险可选项：把 `test_public_repo_safety` 与 `test_finalize_phaded_subtype_reconciliation` 纳入 CI 门（在真实仓库中均通过，仅在无数据沙箱失败）；跨属最小数是否由「仅报告」升级为「硬性强制」。

---

## 6. 依赖与建议顺序

```text
立即可做（A 类，互不阻塞，可并行）
  F3 ── F6        （下游分层语义修正 + 误用审计）
  F4 ── F7        （tier 绑定 Z → 才能开 Z 对账 run）
  F5 ── F11       （解除硬编码 → 才能开 P4）
  F1 ── F14       （616 处置规则 → v2 重算）
  F2              （gRodon 66 闭环，独立）

授权后（B 类）
  F13 ── F14 ── F15      （台账 curation → v2 目录 → 提升决策）
  F7（独立，优先级最高：不重跑 HMMER 即可统一 6.74M 行的 E 值尺度）
  F8 / F9 / F10 / F12     （四类计算分别授权，可并行但各需前置决策 §5 #6/#7/#9）

最后（C 类）
  F16 ── F17
```

**建议的第一批**：A 类全部（F1–F6），因为它们**不需要任何计算授权**，且 F4/F5 是 B 类两个高价值运行的直接前置。

---

## 7. 授权检查点（逐项，保持与主计划一致）

代码、测试、文档、本地合成 fixture 可继续推进。以下动作**每一类**都需要你的明确授权，且**一类批准不等于另一类**：

1. 把 2026-09-21 未跟踪研究快照纳入版本控制 → **F16**
2. 创建任何真实分析 run 或 dated deploy → **F1–F15**
3. 从 723 参考集重建 MAFFT 比对或 HMM → **F10（若需）**
4. 运行 SignalP、Foldseek、结构预测、HMMER against GTDB shards、gRodon → **F7–F12**（四类分别批准）
5. 用真实项目数据重算候选目录 → **F14**
6. 安装/变更环境、改服务器配置、删除或归档数据、push 到 GitHub → **F16–F17**

---

## 8. 验收边界

1. 旧 run 与旧结果**未被覆盖**（v1 三件权威件 SHA-256 逐位不变）。
2. 每条 v2 候选有且仅有一个 `primary_disposition`（F14 验收）。
3. AHSMG/Ser、SignalP/定位、SBD/type、序列模型/功能模型已拆分且**由测试钉死**。
4. 52 条 demotion 在规则层可重复产生（已完成）。
5. 616 条、4,301/6,214、gRodon 66 条差额**有机器可读解释**（F1、F2 收口）。
6. HMM/profile 绑定训练集合、alignment、HMM、工具版本与 SHA-256（已完成）。
7. 发现层不筛除候选、不产生 family call（F3、F6 收口）。
8. 所有新运行使用新 `runs/<run_id>/` 与 dated `deploy/<run_id>/`，线程遵守动态上限。
9. 测试、`compileall`、`git diff --check` 通过。
10. 发布前完成本地、GitHub、deploy 与 run manifest 的权威核对。

---

*本文档为 candidate-only 计划件。所有 profile/HMM/domain/motif/SignalP/结构/系统发育/分类与生长速率证据只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
