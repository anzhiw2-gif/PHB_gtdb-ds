# T141 2026-09-17 Task 0（收尾两件已知小事）状态文档

- **run_id：** `20260917_phaded_housekeeping_amendments_01`
- **run 目录：** `runs/20260917_phaded_housekeeping_amendments_01/`（用 `pipeline/scripts/run_context.py::create_run_layout` 新建，含 `logs/`、`inputs/`、`results/`、`input_contract.json`）
- **依据：** `docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` v2 第 3 节 Task 0（Step 0.1 / 0.2 / 0.3）
- **工作区：** `<REPO_ROOT>`，分支 `main`，HEAD `2e44aaf`
- **阶段性质：** 全本地；未使用 SSH，未连服务器，未提交/启动任何服务器任务
- **边界（candidate-only）：** 本 Task 只做事实修正与计数口径说明。所有 profile、domain、motif、SignalP、结构、系统发育与 Foldseek 证据仍只表示候选同源或功能潜力，**不等同于已验证 PHB/PHA 降解表型**。本 Task 不新增任何 family 判定、不解除任何 architecture conflict、不动任何候选分类结果。

---

## 1. 核心声明：原始 contract 不可追溯修改，amendment 才是当前解释

> **原始 `input_contract.json`（09-16 三份）与 09-16 reconciliation manifest 一律不可追溯修改（not retroactively editable）。**
> 它们在当时如实记录了写盘时刻的已知状态，是历史证据，必须原样保留。
> **自本 amendment 起，后续工作必须以下列 amendment 产物为当前解释（current interpretation）：**
> - `runs/20260917_phaded_housekeeping_amendments_01/results/provenance_amendment.json`（执行 provenance 的当前解释）
> - `runs/20260917_phaded_housekeeping_amendments_01/results/historical_contract_authorization_snapshot.json`（三份历史 contract 的只读快照与 `server_execution_started: false` 原文留证）
> - `runs/20260917_phaded_housekeeping_amendments_01/results/count_discrepancy_amendment.json`（subtype 层与 profile 层计数差值的当前解释）
>
> 引用 09-16 contract 时，**必须**同时引用本 amendment；单独引用历史 contract 会得到已被本 amendment 修正的过时结论。

三份历史 contract 的 `authorization` 字段实测值（本次只读读取，未做任何写入）：

| run | 文件 SHA-256 | size | `server_execution_started` | 读时校验 |
|---|---|---|---|---|
| `20260916_phaded_full_library_signalp_01` | `41915306795ef2e0218576b0a7ccfc3574405e01c458683badbaed603e9594e2` | 4766 | `false` | `read_only_verified: true` |
| `20260916_phaded_full_library_evidence_amendment_01` | `9ffbcbf7c737ca78af583cd3088d7959386a153c3dcaaad04427c88930a7d5c9` | 3816 | `false` | `read_only_verified: true` |
| `20260916_phaded_full_library_reconciliation_amendment_02` | `2a2497bb210eb595f759b5e55d4d9b90f4b60eba49c94ed406d496cd13d6b5ac` | 3332 | `false` | `read_only_verified: true` |

四字段原文（三者一致）：`candidate_only_execution: true` / `formal_scan_authorized: false` / `formal_registry_modified: false` / `server_execution_started: false`。

**全量不改动校验（本 Task 开始前记录、结束后复算）：** 7 个被读取的历史文件（上述 3 份 contract + SignalP 摘要 + Foldseek manifest + 最终矩阵 + reconciliation manifest）SHA-256 与 size 全部 `UNCHANGED`，结果 `ALL_HISTORICAL_UNCHANGED = True`。记录见 `inputs/pre_run_historical_hashes.json` 与 `logs/artifact_sha256_manifest.txt`。

---

## 2. Step 0.1 — provenance amendment（已执行并校验）

### 2.1 SignalP 输入来源：既有摘要，**未重建、未伪造**

计划允许在摘要缺失时以最小合法 JSON 重建。实测**摘要已存在**，因此**未新建任何替代文件**：

- `runs/20260916_phaded_full_library_signalp_01/inputs/signalp_run_summary.json`
  size 641，SHA-256 `960aa328f27ef9a47704446d88532102b1ad5f3a3fdd43d3406d2b2fd607b079`
- 该哈希与 `runs/20260916_phaded_full_library_signalp_01/input_contract.json` 中 `inputs.signalp_run_summary.sha256` 的记录**完全一致** → 摘要与历史 contract 自洽，可直接作为 amendment 输入。
- 摘要内容：`status: completed`、`returncode: 0`、`candidate_only: true`，并含完整 signalp6 命令（含 `--mode fast`、`--torch_num_threads 32`、`--write_procs 8`），命令中的 `--fastafile` 路径为
  `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260916_phaded_full_library_signalp_01/inputs/candidate_interpro_normalized.faa`。

### 2.2 实际执行事实（amendment 记录的当前解释）

| 字段 | 值 | 依据 |
|---|---|---|
| `execution.mode` | `server_dated_deploy` | SignalP 在 T141 dated deploy 上完成 |
| `execution.server` / `account` | `<SERVER_HOST>` / `<SERVER_USER>` | 计划 v2 第 1 行服务器定义 |
| `execution.signalp_returncode` | `0` | `signalp_run_summary.json` |
| `execution.signalp_threads` | `32` | 摘要命令 `--torch_num_threads 32`（同命令另含 `--write_procs 8`，两者语义不同，均原样留在 `signalp_command` 中） |
| `execution.deploy_path` | `${PHB_REMOTE_ROOT}/PHB_gtdb-ds/deploy/20260916_phaded_full_library_signalp_01` | **由摘要命令的 `--fastafile` 路径推导**，非独立观测 |
| `execution.logical_cpu` / `available_memory` | `pending` / `pending` | 本 Task 未接触服务器，**未测即写 pending，未填 0** |
| `execution.gpu_execution` / `gpu_observation` / `server_observed_at` | `not_recorded` | 摘要无 GPU 记录；未执行服务器核对命令 |
| `correction.no_historical_run_overwritten` | `true` | 历史 run 与 contract 均未改写 |
| `profile_boundary` | 本 amendment 不授权任何新 family 判定或 profile 提升 | — |

### 2.3 校验结果（脚本 fail-loud 断言，不满足即抛错）

`run_phaded_housekeeping_amendments.verify_provenance_assertions()` 实测通过：

- `execution.mode == "server_dated_deploy"` ✅
- `execution.signalp_returncode == 0` ✅
- `correction.no_historical_run_overwritten is True` ✅

### 2.4 输入绑定（path / size / SHA-256）

| 输入 | path（简写） | size | SHA-256 |
|---|---|---|---|
| SignalP 运行摘要 | `runs/20260916_phaded_full_library_signalp_01/inputs/signalp_run_summary.json` | 641 | `960aa328f27ef9a47704446d88532102b1ad5f3a3fdd43d3406d2b2fd607b079` |
| Foldseek 合并 manifest | `runs/20260916_phaded_full_library_foldseek_amendment_01/results/foldseek_merge_manifest.json` | 1404 | `cd741ce47a6a51f8ddba45e837ca2f4aef875969b34bedec0bafc58effaebdfa` |
| 最终 Foldseek 矩阵 | `runs/20260916_phaded_full_library_foldseek_amendment_01/results/phaded_full_library_foldseek_subtype_matrix.tsv` | 171507779 | `691805823029b03225858d1632ed0f11204b9350db326353efeb1ae57627868b` |

Foldseek manifest 状态 `completed_candidate_only`（脚本校验通过），其自身记录的 output 哈希 `69180582…` 与实际矩阵文件哈希**一致**。

### 2.5 已知遗留缺陷（只报告，未修）

`pipeline/scripts/amend_phaded_provenance.py` 先对 `provenance_amendment.json` 计算 SHA-256 并把该值写进文件，**随后又重写该文件**（L101→L126），因此文件内嵌的自身哈希对应的是**重写前版本**，与最终文件不符：

- 内嵌值：`4312c1df1bc5c31ababd9a4a0ccb2a479fc94b7fce6ee914f795721cf1c4d900`（size 2867，重写前版本）
- 最终文件：`8b87ed4bfa56d4f30822df3f212537a5473b122686368907631791a2463a021a`（size 3078）

处理方式：**不改动该脚本**（属其他 agent 的既有产物，其测试已通过），改由 `results/provenance_outputs.json` 记录权威的最终 SHA-256，并显式标注 `matches_final_file: false`，避免下游误用陈旧哈希。`provenance_readme.md` 只写一次，其内嵌哈希与最终文件一致（`matches_final_file: true`）。

---

## 3. Step 0.2 — 计数不一致 amendment（已解释，差值全部归因）

### 3.1 结论：`resolved_layer_precedence`，差值 33 + 1 全部找到来源，未解释 = 0

实测从 171,507,779 字节矩阵流式重算 109,087 行（`row_count = 109087`），与 09-16 manifest 记录值**逐键完全一致**（`manifest_matrix_count_mismatches = {}`，不一致条目数 0）。

| profile 层状态 | profile 计数 | subtype call | subtype 计数 | 差值 | 被改判行数 | 是否完全解释 |
|---|---|---|---|---|---|---|
| `profile_ambiguous_family` | 9664 | `ambiguous_family_within_superfamily` | 9631 | **33** | 33 | ✅ |
| `profile_ambiguous_superfamily` | 42 | `ambiguous_superfamily` | 41 | **1** | 1 | ✅ |

### 3.2 来源：判定分支优先级（不是数据错误，也不是少算）

`pipeline/scripts/build_phaded_subtype_matrix.py` L210–L219 的分支顺序为：

1. `manual_review_decision == "hold_architecture_conflict"` **或** `domain_evidence == "domain_conflict"` → `hold_architecture_conflict`
2. `manual_review_decision == "hold_gene_model_or_structure"` → `hold_gene_model_or_structure`
3. `assignment_status == "unassigned_PhaDED_like"` → `unassigned_PhaDED_like`
4. `assignment_status == "ambiguous_superfamily"` → `ambiguous_superfamily`
5. `assignment_status == "ambiguous_family"` → `ambiguous_family_within_superfamily`

hold 分支（1、2）在 ambiguous 分支（4、5）**之前**求值，因此一个在 profile 层被判为 ambiguous、同时在 domain 层为 `domain_conflict` 的候选，其 `subtype_call` 记为 `hold_architecture_conflict`，不再计入 subtype 层 ambiguous 计数。**两层回答的是不同问题（profile 层问"profile 归属是否唯一"，subtype 层问"最终候选归类"），都没有算错**；subtype 层的当前顺序是更保守的一种。

### 3.3 逐条证据（34 行，全部为 `hold_architecture_conflict` + `domain_conflict`）

| 触发分支组合 | 行数 |
|---|---|
| `domain_conflict`（19 条人工复核之外的 33 行） | 33（含 32 行 ambiguous_family + 1 行 ambiguous_superfamily） |
| `manual_hold_architecture_conflict + domain_conflict` | 2（均为 ambiguous_family；这 2 行同时命中分支 1 的两个条件） |
| **合计** | **34** |

- 33 行 `profile_ambiguous_family` 全部改判为 `hold_architecture_conflict`（31 行仅 domain 层 conflict，2 行同时为人工 hold）。
- 1 行 `profile_ambiguous_superfamily`（`GCA_005787565.1|SYAP01000036.1_14`）同样改判为 `hold_architecture_conflict`。
- 逐行证据：`results/reclassified_ambiguous_hold_rows.tsv`（34 行 + 判定字段），摘要见 `results/count_discrepancy_report.md` 第 3 节。

**自洽核对：** `hold_architecture_conflict` 总数 83 = 49（profile_trained_hit）+ 33 + 1；`hold_gene_model_or_structure` 2（均 profile_trained_hit）；subtype 层合计 9631+41+67342+356+83+2+31632 = 109087，profile 层合计 9664+42+67749+31632 = 109087，两层总数与矩阵行数三者一致。hold 行按触发条件的完整分布：`domain_conflict` 80、`manual_hold_architecture_conflict + domain_conflict` 3、`manual_hold_gene_model_or_structure` 2（合计 85）。

### 3.4 反事实与处置

若把 ambiguous 分支提到 hold 分支之前，`subtype_call_counts` 会读作 9664 / 42——但那会把"因 domain 层结构冲突而被 hold"的行重新计入 ambiguous，属重复计数。**故不修改历史 manifest 的数值**：09-16 `reconciliation_manifest.json` 保持原样（SHA-256 `1f55f4e2f1d414f8648ab01f986d92054836033c8a6803d51a19c017fa1c8189`，实测 `historical_manifest_unmodified = true`，读取前后哈希与 mtime 均不变），由本 amendment 提供当前解释。

---

## 4. 产物清单与 SHA-256（实测）

`logs/artifact_sha256_manifest.txt`（其自身 SHA-256 `c446552814c08846ade884f9c8f9da47d46cbd39f5440d4b108979a7adb28975`，size 1423）列出本 run 目录下除自身以外的全部文件：

| 相对路径 | size | SHA-256 |
|---|---|---|
| `input_contract.json` | 2141 | `91b4b3eefc73d0de265c076dd9097b399e00d747fa7032a0cab6490f29f44d55` |
| `inputs/input_manifest.json` | 2561 | `94b4a9d7291e9b31c827cbbc2b768363e3a298be2032da10cc2d2658007c7ebd` |
| `inputs/pre_run_historical_hashes.json` | 1418 | `a8d78f2a3bcd4f8dc37c6037d722d94a97b1196abf2397064f32ecd7164395dc` |
| `logs/artifact_sha256_manifest.txt` | 1423 | `c446552814c08846ade884f9c8f9da47d46cbd39f5440d4b108979a7adb28975` |
| `results/count_discrepancy_amendment.json` | 25988 | `a03e929e6dba08adefdfcbb9d390a1b5a9df6de98baa0224bf14542ef374c272` |
| `results/count_discrepancy_outputs.json` | 1123 | `36a2372e903d07d85957c8e90f591e540c222deaa840186e1dd59f4f6638bc51` |
| `results/count_discrepancy_report.md` | 9655 | `013493baa4dca5a5f68ee79707699a371a9bde3cfb0beb740dccfeb94d4a9f74` |
| `results/historical_contract_authorization_snapshot.json` | 2572 | `4542d6df904893f544d55afff6679be8319c33c34c61ad78aa58b7f8ab566912` |
| `results/provenance_amendment.json` | 3078 | `8b87ed4bfa56d4f30822df3f212537a5473b122686368907631791a2463a021a` |
| `results/provenance_outputs.json` | 1191 | `a41097afff751d84c61bb09dabf8e5f3bc177097362fb9d4d9ae8e750e94f532` |
| `results/provenance_readme.md` | 657 | `c00f3ca063170393133efd5e97ef50b26208e41cdcd0e1993a265a62c7dc230e` |
| `results/reclassified_ambiguous_hold_rows.tsv` | 8033 | `92ad5c73a5397c500494785cda2a4d04611abcd55f676841cff4d135b0cbf992` |

代码（本 Task 新建，先写失败测试后实现）：

| 文件 | size | SHA-256 |
|---|---|---|
| `pipeline/scripts/run_phaded_housekeeping_amendments.py` | 30540 | `19c7531415be67d9b2486931e98dae61c8cd970ed9f74fde7144d8d4a8eef13f` |
| `pipeline/tests/test_phaded_housekeeping_amendments.py` | 15389 | `619d28319e304a5f62001a9987cbfc66d50704549b4930f680425afbd369d77c` |

说明：`input_contract.json` 由既有脚本 `amend_phaded_provenance.py` 写出（其设计行为是写到 `output_dir.parent`），内容为该脚本的 provenance 契约；完整输入绑定另见 `inputs/input_manifest.json`（由 `run_context.build_input_contract` 生成，其中 GTDB taxonomy/metadata/tree 因本 Task 不消费而如实为 `pending`）。

---

## 5. 验证（实测结果与时间点）

| 命令 | 结果 |
|---|---|
| `python -m unittest pipeline.tests.test_amend_phaded_provenance -v` | **OK**（Ran 1 test）— 基线确认，先于任何改动 |
| `python -m unittest pipeline.tests.test_phaded_housekeeping_amendments -v` | 实现前：`ModuleNotFoundError`（失败先行）；实现后：**OK**（Ran 10 tests） |
| 聚焦合计（两模块） | **OK**（Ran 11 tests，exit 0） |
| `python -m unittest discover -s pipeline/tests` | 01:01 快照：**Ran 325 tests, OK (skipped=1)**；01:02:2x 快照：`Ran 346 tests, FAILED (errors=21, skipped=1)`；01:03 快照：`Ran 346 tests, FAILED (failures=3, errors=1, skipped=1)`（4 项**全部**位于 `test_assign_phaded_catalytic_domain_type.RunIntegrationTests`） |
| 排除他人模块后运行全部 94 个测试模块 | **Ran 325 tests, 0 errors, 0 failures, skipped=1 → OK** |
| `python -m compileall -q pipeline` | 01:02:2x：**exit 0** |
| `git diff --check` | **exit 0**（无输出；未执行任何 git 写操作） |

**关于 `discover` 出现非零失败数的诚实说明（非本 Task 引入，只报告不修）：** 全部失败项**只**位于 `pipeline/tests/test_assign_phaded_catalytic_domain_type.py`（Task 4 模块，`git status` 显示为 `??` 未跟踪，属并行 agent），且**从未**出现在 `test_phaded_housekeeping_amendments` 或 `test_amend_phaded_provenance`。该模块在我执行期间被反复写入（实测 mtime 由 01:02:13 → 01:02:18），其中一次 `compileall` 抓到中间态语法错误（该文件 L349 `SyntaxError: expected ':'`），下一次即恢复；其失败项数由 21 errors 变为 3 failures + 1 error，呈持续变动，属其"失败测试先行 → 实现"的正常中间态。**我未修改该文件，也未修改其任何代码**。排除该模块后套件全绿（94 modules，325 tests，0 error，0 failure，skipped=1 → OK）。

---

## 6. 本 Task 未做的事（合规声明）

- 未 SSH、未连服务器、未提交或启动任何服务器任务（`logical_cpu` / `available_memory` 因此如实为 `pending`）。
- 未修改、未覆盖、未删除任何历史 `runs/`、`results/`、三份 09-16 contract、09-16 reconciliation manifest、`AGENTS.md`、历史 HMM。
- 未修改 `pipeline/config/formal_scan_models.tsv`。
- 未执行 `git commit` / `push` / `clean` / `reset`；未删除任何文件。
- 未修改其他 agent 的产物（含 `amend_phaded_provenance.py` 的内嵌自哈希缺陷，仅在第 2.5 节报告）。
- 未降低校准 gate、未合并/删除任何家族、未新增 family 判定、未改动任何候选的 `subtype_call`。

---

## 7. 后续引用口径（必须遵守）

1. 引用 SignalP 执行事实时，写 **`server_dated_deploy`（T141 dated deploy，returncode 0）**，并附 amendment 路径；**不要**据 09-16 contract 的 `server_execution_started: false` 推断"SignalP 未在服务器执行"。
2. 引用 ambiguous 相关计数时，**必须点明所在的层**：profile 层 `profile_ambiguous_family = 9664` / `profile_ambiguous_superfamily = 42`；subtype 层 `ambiguous_family_within_superfamily = 9631` / `ambiguous_superfamily = 41`；差值 33 / 1 由 hold 分支优先级解释，见第 3 节。
3. 引用历史 contract 原文时，**同时**引用 `historical_contract_authorization_snapshot.json` 作为"原文 + 只读校验"留证。
4. 所有结论维持 candidate-only：计算证据不构成 PHB/PHA 降解表型的验证。

---

## 8. P1 修复：`input_contract.json` 契约一致化（2026-09-17，Phase 1 验证报告 P1）

> **本节为 P1 修复任务追加，第 1–7 节原文一字未改。** 本修复只改本 run 的契约登记，不改任何计算结果、不新增任何 family 判定。

### 8.1 问题（引用验证报告 P1）

来源：`runs/20260917_phaded_phase1_verification_01/results/phase1_verification_report.md`，size 14526，SHA-256 `a91737a16263af02f83a0187c9c729e33eb71b9e964d511b47291ca738f305f9`，问题清单 **P1（合规，中等）**。原文逐字引用见 `results/p1_contract_conformance_fix.json` 的 `problem.verbatim_finding`。

| 缺陷 | 修复前实测 | run_context 1.0 / AGENTS.md 要求 | 修复后实测 |
|---|---|---|---|
| `run_dir` | 缺失 | 契约须含 run 目录 | `<REPO_ROOT>\runs\20260917_phaded_housekeeping_amendments_01` |
| `generated_at` | 缺失 | 契约须含生成时间戳 | `2026-09-16T17:23:35+00:00`（UTC，= 本地 2026-09-17 01:23:35） |
| `gtdb` 块 | **缺失**（taxonomy/metadata/tree 连 `pending` 都没有） | AGENTS.md：GTDB taxonomy/metadata/tree 必须记录，缺失写 `pending` | 三项齐备，`path: null` + `status: pending` + `sha256: null` + `size: null` |
| `inputs[*].status` | 3 条全部无该字段 | run_context 1.0：每条输入须有 status | 3 条全部 `verified`（size 与 SHA-256 实测重算一致） |
| 顶层 `status` | `completed_candidate_only`（不在取值域内） | run_context 取值域 `{verified, pending}` | `pending`（唯一原因是三条 GTDB 槽为 pending），原值迁至 `run_completion_status` |

**根因（只报告，未修）**：`pipeline/scripts/amend_phaded_provenance.py` L127–L142 用自带字典写本 run 的 `input_contract.json`，未走 `run_context.write_input_contract`，故缺 `run_dir`/`generated_at`/`gtdb`，`inputs` 无 `status`，顶层 `status` 直接取自 Foldseek manifest 的 `completed_candidate_only`。该脚本属其他任务的既有产物（其自身测试通过），本修复不动它。

### 8.2 修复动作（信息零丢失）

1. **先留证**：原契约按字节复制为 `results/input_contract_superseded_pre_p1_fix.json`（2141 B，SHA-256 `91b4b3ee…`），与改写前实测哈希逐位一致（脚本 fail-closed 校验：不符即拒改）。
2. **取骨架**：用 `pipeline/scripts/run_context.py::build_input_contract` 生成契约骨架，GTDB taxonomy/metadata/tree 显式声明为 `None` → 如实落为 `pending`。
3. **重算绑定**：3 条输入逐条重算 size 与 SHA-256 并与原声明比对（**全部一致**）后，才补 `status: verified`；未新增、未改动任何路径与哈希。
4. **保留原文**：`execution` 与 `phenotype_boundary` 逐字保留；原 `status` 值迁到 `run_completion_status: "completed_candidate_only"`。
5. **顶层 status**：`pending`，含义与其余 4 个 run 完全相同（仅因 GTDB 三槽 pending），并写入 `status_note` 说明。
6. **未直接用 `write_input_contract`**：该函数只写自身字典，会丢掉本 run 已有的 `execution` 与 `phenotype_boundary`；故用 `build_input_contract` 取骨架后回填保留字段（记录于 `remediation.why_not_write_input_contract`）。

### 8.3 逐项差异（`diff` 字段机器可读版本见修复记录）

| 项 | 变化 |
|---|---|
| 顶层键新增 | `gtdb`、`run_dir`、`generated_at`、`run_completion_status`、`run_inputs_status`、`superseded_contract`、`status_note` |
| 顶层键删除 | 无 |
| 顶层值变更 | `status`：`completed_candidate_only` → `pending`（原值保留在 `run_completion_status`） |
| `inputs` 条目数 | 3 → 3（`signalp_run_summary`、`foldseek_merge_manifest`、`final_foldseek_matrix`），顺序与路径、size、SHA-256 **全部不变** |
| 每条 `inputs` 新增键 | `status`（值 `verified`），`values_changed = {}` |
| 字段顺序 | 条目内由 `path,size,sha256` 变为 `path,status,sha256,size`（run_context 口径），内容不变 |

### 8.4 产物与 SHA-256（实测）

| 相对路径 | size | SHA-256 |
|---|---|---|
| `input_contract.json`（修复后） | 3496 | `6164a46830c82662df0da9eed1e0a55b680fa726b3787a4912c9140132167317` |
| `results/input_contract_superseded_pre_p1_fix.json`（修复前字节） | 2141 | `91b4b3eefc73d0de265c076dd9097b399e00d747fa7032a0cab6490f29f44d55` |
| `results/p1_contract_conformance_fix.json` | 15552 | `2e9aad6373182e92b8952d48c63e60a1da3c988a5eb1f7771a66c9a58fc3b6ed` |
| `logs/p1_contract_conformance_check.py` | 14301 | `2d0042fdb515e00e6582ab938874aef0fe1fc6fee459522a6e6f4cf9c87fde75` |
| `logs/repair_p1_input_contract.py` | 13080 | `19136b936e5ced8af4a9be17e99bf847015439c5351b0f3128448a5371391ba0` |
| `logs/write_p1_fix_record.py` | 12476 | `891036c1a6c4bc421df2a0b4a3daf4a306469b76f019ae18563a1560225f5e66` |
| `logs/p1_contract_conformance_check_pre_fix.log` / `.json` | 5456 / 5439 | `8a9151df5ca0b5cf72efedda7841401b67e8e2017a423ce9e54969e87b676a37` / `73d79a8f6fff6a8ac6d50a99e2a8d0b19a2bb8e035a66eaa1adf87280342716a` |
| `logs/p1_contract_conformance_check_post_fix.log` / `.json` | 7966 / 8302 | `d6c7b85d4df44923439fa2f2bf0f26ba60d464a9408efa25b5013a1de4bd994f` / `ee5e4938437a4b2d983aa27c17e6d466656a9429bb83f9d392c751843b357791` |
| `logs/p1_contract_repair.log` | 1432 | `4336e0259b88e6d2c79d8d2efd05735fb2e967d45a1b246861bc2c6304a95599` |
| `logs/p1_contract_conformance_repair_summary.json` | 7867 | `98b78957d2fb4258da8fec598cd257589e8f18f5ac4bbf654341802d23d61a0e` |
| `logs/p1_acceptance_commands.log` | 634 | `5acf5108d147da3c7ea651d5245ab906e64952362071389cc5e86195e538fdb0` |
| `logs/p1_full_suite_snapshot.log` | 192 | `ce99e6c919ad8a90aae99bf2fdc972b786952ba01af49abf5fb8bb18745a24c7` |
| `logs/p1_contract_conformance_fix_record_write.log` | 476 | `29b1b35d522b609a189158f4f915c5c9035387e205ad08d03d27b7057ac59261` |
| `logs/p1_integrity_crosscheck.py` | 7984 | `ab49826f6520ffa87e2c9b1ee1f7fb1d9cbfb68a6224a3394c6011435ef5017f` |
| `logs/p1_integrity_crosscheck.log` | 2432 | `6132e74cf53a502f84c1c0dc8a1726c8d22a0b088b1b1e1620d64e9368bffa0d` |
| `logs/write_p1_post_fix_manifest.py` | 3479 | `ddd7c39202c478e0b030f5629dc72aef2c564b4fb39e64f74af27bc4987e345c` |
| `logs/p1_post_fix_artifact_sha256_manifest.txt` | 4162 | `358857900845ed4e66a7e87f8398e2515df58cf6266bc92a1cf96fb2a0ac94ca` |

该 manifest 列出本 run 目录下除自身以外的全部 29 个文件（run 目录共 30 个文件），已用 `logs/p1_integrity_crosscheck.py` 逐条复算：**28/28 条目与磁盘一致**（自身日志因在本脚本读取 manifest 后重写而显式排除），修复前 11 条历史产物中 10 条哈希不变、唯一变化即 by-design 的 `input_contract.json`。

**过程披露（诚实记录）：** 该 manifest 的第一版由 PowerShell `Set-Content -Encoding utf8` 生成，实测**带 UTF-8 BOM**（首个条目的解析因此偏移），且当时文件清单尚未包含后续新增的交叉校验脚本。已改为由 `logs/write_p1_post_fix_manifest.py` 以 Python 生成（UTF-8 无 BOM、LF），并迭代到不动点：连续两次重新生成的 size/sha256 **完全相同**（4162 / `35885790…`），证明清单生成幂等、内容不再自指振荡（清单内不内嵌自身哈希，也不回显交叉校验日志自身的哈希）。

`logs/artifact_sha256_manifest.txt`（01:01:58 写就）**未改动**；其中 `input_contract.json` 条目（`91b4b3ee…` / 2141）现已过时，权威清单改为 `logs/p1_post_fix_artifact_sha256_manifest.txt`（逐文件列出当前哈希，并在头部显式标注该条目的新旧哈希）。修复记录本身**不内嵌自哈希**（自哈希内嵌正是本 run 第 2.5 节记录的缺陷模式）。

### 8.5 校验输出（本 run `logs/`）

- 修复前：`logs/p1_contract_conformance_check_pre_fix.log` → `SUMMARY passed=5 failed=8`、`RESULT: FAIL`、exit 1（复现 P1：缺 `run_dir`/`generated_at`/`gtdb`、`inputs` 无 `status`、顶层 status 越域）。
- 修复后：`logs/p1_contract_conformance_check_post_fix.log` → `SUMMARY passed=16 failed=0`、`RESULT: PASS`、exit 0；含与其余 4 个 run 契约的关键字段集合比对、GTDB 三槽齐备且为 `pending`、3 条输入 size+SHA-256 重算一致、`run_id`/`run_dir` 与目录一致。
- 参照契约（只读）：`20260917_phaded_classification_authority_01`、`20260917_phaded_reference_residue_mapping_01`、`20260917_phaded_x1_hydrophobicity_01`、`20260917_phaded_catalytic_domain_type_01`。
- 完整性交叉校验：`logs/p1_integrity_crosscheck.log` → `SUMMARY failures=0`、`RESULT: PASS`（exit 0）。校验项：`20260901_*`–`20260916_*` 中修改时间 ≥ 2026-09-17 的文件数 **0**；4 个参照契约 mtime 均停留在 Phase 1 时段（01:01–01:13，未因本修复变动）；Phase 1 验证报告 SHA-256 仍为 `a91737a1…`，与修复记录所载一致；6 个 JSON 产物全部可解析；post-fix manifest 28/28 条目与磁盘一致。

### 8.6 验收命令（01:24 快照，`logs/p1_acceptance_commands.log`）

| 命令 | 结果 |
|---|---|
| `python -m unittest pipeline.tests.test_phaded_housekeeping_amendments pipeline.tests.test_amend_phaded_provenance` | **OK**（Ran 11 tests），exit 0 |
| `python -m unittest discover -s pipeline/tests` | **Ran 383 tests … OK (skipped=1)**，exit 0 |
| `python -m compileall -q pipeline` | exit 0 |
| `git diff --check` | exit 0（已跟踪改动仍仅 ` M AGENTS.md`，与本修复无关） |

**诚实说明：** Phase 1 基线为 375 tests；当前 383 tests 是并行 agent 新增（+8），本修复未改任何 `pipeline/` 源码。本修复过程中有一次中间态全套运行报 `failures=1`/`errors=4`，全部位于 `pipeline/tests/test_phaded_calibration_field_semantics.py`（另一 agent 的 P4 模块，字段语义改名，属"失败测试先行"中间态），随后快照即为全绿；**我未修改该模块，也未据其失败做任何改动**。

**01:26 复核（同一套命令，第二次独立快照）：** `Ran 400 tests … FAILED (failures=3, errors=1, skipped=1)`，失败项**全部**位于 `pipeline/tests/test_complete_phaded_motif_completion_full.py`（另一 agent 的 Task 7 模块）。证据表明该模块正处于"失败测试先行 → 实现"的写入中间态：测试文件 mtime `01:25:35`，其实现 `pipeline/scripts/complete_phaded_motif_completion_full.py` mtime `01:26:45`（即在我复核的同一分钟内被写入），且相邻两次复核的失败数由 `failures=4, errors=1` 变为 `failures=3, errors=1`。**该失败与本修复无关：本修复未修改任何 `pipeline/` 源码（`git status --porcelain` 的已跟踪改动始终仅 ` M AGENTS.md`，且该改动早于本阶段），只改了本 run 的 `input_contract.json`、在本 run 内新增 logs/results 产物、并追加本文档第 8 节。** 因此本修复的权威验收证据是 01:24 快照；01:26 的全套失败应记在并行 agent 的 Task 7 模块名下（只报告，不修）。

### 8.7 未做的事与边界

- 未修改 `pipeline/scripts/amend_phaded_provenance.py`（根因只报告）。
- 未改动其余 4 个 0917 run、Phase 1 验证 run、`20260901_*`–`20260916_*` 历史 run、`results/`、`AGENTS.md`、`pipeline/config/formal_scan_models.tsv`、历史 HMM。
- 未执行 `git commit` / `push` / `clean` / `reset`；未删除任何文件（原契约以字节快照完整保留）。
- 未新增任何 family 判定、未改动任何候选 `subtype_call`、未降低任何 gate。所有证据仍为 candidate-only，**不等同于已验证 PHB/PHA 降解表型**。

### 8.8 引用口径更新

引用本 Task 的输入契约时，以修复后的 `input_contract.json`（3496 B，`6164a468…`）为准；需要 2026-09-17 01:01 时刻的原文时，引用 `results/input_contract_superseded_pre_p1_fix.json`（2141 B，`91b4b3ee…`，与修复前字节一致）。P1 已修复；验证报告的 P2–P5 不在本次修复范围内。
