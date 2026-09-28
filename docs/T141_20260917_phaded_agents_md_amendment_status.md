# T141 — AGENTS.md 修正（发现层与模型可复现性条款 + commit/push 授权明确化）：状态文档

**任务：** 事项 3 —— 把发现层模型规则与模型可复现性规则写进 `AGENTS.md`，并明确 commit/push 授权条款
**Run：** `runs/20260917_phaded_agents_md_amendment_01`（本地；**未使用服务器计算**，仅只读查询资源快照）
**修改的项目文件：** `AGENTS.md`（已跟踪；本轮明确授权范围）；**新增** `pipeline/tests/test_agents_md_rules.py`
**日期：** 2026-09-17
**状态：** `completed_with_pending_operator_confirmation` —— 修正已实施并通过全部验收；**commit/push 条款为提案，待操作者一行确认**

> **边界声明：** 本 run 只改**规则文件**、**新增**一个规则回归测试与本文档，**未**改任何 pipeline 逻辑、`pipeline/config/`、历史 HMM、历史 run 或 `results/`。本文中所有哈希、计数与判据均由本 run 的脚本从磁盘重算；凡引用而非重算的数字均已逐处标注。修正后的 `AGENTS.md` 仍然声明：本项目 HMM/profile/domain/motif/定位/系统发育证据只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。

---

## 0. 一句话结论

工作区里 `AGENTS.md` 相对已跟踪版本**少了那行授权条款**（990 B vs 1,116 B，逐字差异见 §1），而该行是全项目**唯一**约束 commit/push 的常驻规则 —— 全库扫描证实它**零机器强制**（无测试、无脚本读取 `AGENTS.md`），删除后自动加载的规则里**没有任何一句**限制提交或推送。本轮按任务要求把它替换为**更明确、可执行的版本**（本地 commit 允许但每次须记录；push 仍需操作者授权），并基于本会话实测证据新增**发现层模型**、**模型可复现性**、**失败证据归档**三组规则，共新增 6 条规则；**既有 8 条规则一条未删**（逐条字符串比对，缺失列表为空）。修正后 2,190 B / 16 行 / 14 条规则（增长 2.15×，规则上限 30 行）。

---

## 1. 修改前后（Step 1 + Step 2）

| 对象 | 字节 | 行 | 规则条数 | SHA-256 |
|---|---:|---:|---:|---|
| `git show HEAD:AGENTS.md`（已跟踪，blob `cade387f…`） | 1,116 | 11 | 9 | `5b9f11ad5f0decbd4a2a4b4d3d798e99137e4b07b5bb956515e918267dc225b5` |
| 工作区**修正前**（快照 `inputs/agents_md_before.md`） | 990 | 10 | 8 | `e42adbe3cb7d5b6e6cd713c5c3a797060ec638f838844a8811df48b8a5c71f2e` |
| 工作区**修正后**（快照 `inputs/agents_md_after.md`） | 2,190 | 16 | 14 | `82230bee098d36ecde8876eec8b2964b021fabb33b55482058e013ba4df9ff30` |

- 行尾全为 **LF**；`inputs/agents_md_after.md` 与工作区 `AGENTS.md` **逐字节一致**（`acceptance_summary.json` → `agents_md.after_matches_worktree = true`）。
- **diff 语义须注意：** 操作者的删除**未暂存**，索引仍是 HEAD，因此 `git diff -- AGENTS.md` 是 **HEAD → 修正后**（`1 file changed, 6 insertions(+), 1 deletion(-)`），**不是**"修正前 → 修正后"。三种视图分别存档：
  - `inputs/agents_md_worktree_vs_head.diff`（HEAD → 修正前；唯一差异 = 删除那行）
  - `inputs/agents_md_index_to_worktree.diff`（HEAD → 修正后）
  - `inputs/agents_md_before_vs_after.diff`（快照 before → 快照 after）
- **回退方式：** `inputs/agents_md_before.md` 即修正前字节，回退不需要任何 git 操作。

### 1.1 删除那行的实际影响（实测，非推断）

| 项 | 实测 |
|---|---|
| 全库扫描命中（排除 `.git`/`.uv-cache`/`research` 第三方树） | **77 行 / 32 个文件**（`results/agents_md_clause_dependency_scan.{tsv,json}`） |
| 逐字复述该条的**其它**文件 | **0 处** —— 该句现在**只存在于 git 历史**（HEAD blob） |
| 机器强制 | **零**：`grep -r AGENTS pipeline/` 仅 1 处命中（`finalize_phaded_subtype_reconciliation.py:189` 把它记为**溯源输入**，不读规则）；`pipeline/tests/` 在本轮之前**无任何**测试引用它（本轮的 `test_agents_md_rules.py` 是第一个） |
| 残留约束位置 | 交接文档 §9/§10 第 5 条、`docs/T141_20260906_conversation_process_and_decision_log.md:27,287`、`docs/T141_20260912_project_handoff.md:36` —— 均为**时点性交接记录**，非常驻规则 |
| 为何后果真实 | `AGENTS.md` 是 agent 的**自动加载指令面**（本会话实测：修改后系统立即推送新内容并标注 "This file changed after it was loaded"）。删除后该面上**再无任何一句**约束提交/推送 |
| 当下的风险暴露度 | `git status --porcelain` = **1 条已跟踪改动（即 `AGENTS.md`） + 288 个未跟踪文件**；最后一次非 `AGENTS.md` 提交 `2e44aaf`（**2026-09-03**），即 **14 天**跟踪层成果只有本地一份；`.gitignore` L67 排除 **`runs/`**（本地 15.61 GB、服务器 1,016 GB 产物按设计不在版本控制内） |

---

## 2. 新增/修订的 6 条规则（逐条证据见 `results/amendment_rationale.md`）

1. **发现层 HMM 允许训练，但只能召回** —— 允许用 DED 家族序列（含 `annotation_only`）训练发现层 HMM，层名固定 `discovery_hmm_uncalibrated` 且每份输出带标注；只用于**召回**，**不得**用于筛选、删除或降级任何候选。证据：HMM `70e1b347…`（210,746 B）、`LENG 451`/`NSEQ 270`（本 run 读 HMM 文件头实测）；**对 563 条 x₁ 混淆嫌疑命中 443 条 = 78.69%，无判别力**（引用 `discovery_manifest.json`，本 run 复核其为 32,644 B / `982a70a4…`）。
2. **发现层不得进 registry、不得产生 family 判定；只能对既有中间产物打分** —— 证据：G3 记 `formal_scan_models.tsv` 前后同为 `8a1c0508…b7a3`（本 run 实测当前仍 450 B / 同哈希）；G4 记 `subtype_call_rows_changed = 0`、`call_columns = []`；manifest 记 `registry_modified = false`、`new_family_call_made = false`。
3. **模型可复现性：绑定比对与输入哈希，禁止声称逐位可复现** —— 证据：MAFFT **7.525** 同命令两次得到不同比对；本 HMM 绑定 `2084c4ce…`（本 run 实测该比对 271,045 B，且该 run 目录下 `*.aln` **仅此一个**）。第一次比对的哈希**转引自** Cys 状态文档（其字节未落盘，本 run 无法独立复核 —— 已在理由书中标注）。
4. **失败证据归档：可归档、禁止删除** —— 证据：服务器 `du -sh runs/` = **1016G**（两个时点 `19:38:46` / `19:41:16` 一致，`logs/server_resources_snapshot.txt`），53 个 run 目录，项目总量 1.7T，`/home/data` **仍余 100T**（故属卫生规则而非救火；如实标注）。规则同时**逐字保留**原句"清理前先确认精确路径与可恢复性"。
5. **本地 `git commit` 允许但每次须记录；push 到 GitHub 仍需操作者明确授权（待一行确认）** —— 理由（两条，均有实测支撑）：14 天未跟踪成果只有本地一份（§1.1 末行）；push 涉及公开面脱敏判断，不可机械判定。**这是唯一一处有意不同于被删原句的措辞**，已在规则文本内标注"待操作者一行确认"。
6. **通用授权条不再含"提交"二字** —— 改为"重算、安装、配置变更、删除与推送仍须明确授权；只读审计应保持只读"，以免与第 5 条自相矛盾；测试断言钉住该分工（原句不得逐字回归，通用条必须仍在）。

> **既有规则零删除的实测保证：** `results/amendment_manifest.json` → `rule_preservation.pre_existing_bullets_absent_from_after = []`（字段 `all_pre_existing_rules_preserved = true`）。该检查在本轮**真实抓到过一次回归**：首版把「保留历史运行残留和失败证据；清理前先确认精确路径与可恢复性」整句改写为归档条，检查立即报缺失；改为**原句保留 + 归档条另起一条**后才通过。

---

## 3. 交叉一致性（Step 4）

| 关系 | 判定 |
|---|---|
| 发现层规则 vs「新 run / 不覆盖历史 results」 | 无冲突（发现层产物只存在于 dated run；对既有候选池只读打分） |
| 发现层规则 vs「服务器只执行 dated deploy 中绑定源码」 | 无冲突（该 run 的脚本已绑定并逐文件哈希一致） |
| 发现层规则 vs「不得启动 GTDB 全量重扫」（交接文档 §10 第 6 条） | **已显式划界**：规则正文写"**只能对既有中间产物打分**"+"**全库重扫仍需单独授权**"，消除"未禁即允许"的读法 |
| 可复现性规则 vs「必须记录路径/版本/大小/SHA-256」 | **细化**：规定"工具不可位复现时记**哪一个**哈希并显式声明"，不推翻原条 |
| 归档规则 vs「保留历史运行残留和失败证据；清理前先确认…」 | 无冲突（原句**逐字保留**；归档条另起且含"**禁止删除**"） |
| commit 条 vs 通用授权条 | 已消除冲突（通用条删去"提交"，由专项条治理；测试钉住） |
| commit 条 vs 交接文档 §10 第 5 条（"不得执行 git commit/push/clean/reset"） | **存在需操作者裁决的差异，本 run 不擅自裁决** —— 已列为待确认第 1 项 |

---

## 4. 验收证据（全部实测）

| 检查 | 结果 | 证据 |
|---|---|---|
| 新测试**先红**（test-first） | `Ran 10 tests … FAILED (failures=9)` | `logs/02_test_first_red.log`（59 行 `AssertionError` 明细） |
| 新测试**后绿** | `Ran 10 tests … OK` | `logs/03_test_first_green.log` |
| **剔除并发 agent 在写模块**后的完整套件 | **`Ran 534 tests … OK (skipped=1)`，exit 0** | `logs/full_test_suite_excluding_inflight.log` |
| 套件计数归因 | **534 = 522 文档基线 + 本 run 10 + 并发 agent 2**（该 2 个落在其 19:45:13 改写的 `test_acquire_phaded_reference_panel_amendment.py` 内；被剔除的两个在写模块另含 36 个测试） | `results/acceptance_summary.json` → `full_suite_excluding_inflight.attribution` |
| 完整套件（**不过滤**，如实记录） | 首跑 `Ran 568 tests … FAILED (failures=2, errors=6, skipped=1)`；复跑 `Ran 570 tests … FAILED (errors=1)` | `logs/full_test_suite.log`、`logs/full_test_suite_rerun.log` |
| **失败归因** | 全部落在**并发 agent 19:42–19:45 正在写的模块**：`test_finalize_phaded_superfamily_gate.py`（其实现脚本 `finalize_phaded_superfamily_gate.py` **不存在** → test-first 红灯）、`test_credit_phaded_axis_anchored_controls.py`（其脚本 19:43:20 刚落盘，`reconcile_phaded_reference_panel.py` 19:43:36 被改写）。**无一涉及 `AGENTS.md` 或本 run 的测试** | 上述两份日志 + 模块 mtime 实测 |
| `python -m compileall -q pipeline` | exit 0 | `logs/compileall.log` |
| `git diff --check` | exit 0（无输出） | `logs/git_diff_check.log` |
| **HEAD 未变、无 commit/push** | `2e44aafe2c91263447f37bac625f97c56ddafc05`（= 交接文档基线口径），`head_matches_documented_baseline = true`；`commit_or_push_performed = false` | `results/acceptance_summary.json` → `git`；`logs/git_head_check.txt` |
| run 内产物清单复验 | **47 条** path/bytes/sha256，`missing 0 / size 0 / sha256 0 / volatile 0 / 未收录 0` → **`MATCH`**；外部产物（`AGENTS.md`、本状态文档、新测试）另列并复验 → **`MATCH`** | `logs/artifact_sha256_manifest.txt`、`logs/artifact_manifest_verification.json` |
| 保护文件未变 | `pipeline/config/formal_scan_models.tsv` = **450 B / `8a1c0508…b7a3`**；`runs/20260917_phaded_motif_completion_full_01/results/motif_completion_full.tsv` = **291,698,708 B / `4f182243…`**（与 Cys manifest 记录逐字一致） | `results/acceptance_summary.json` → `protected_files_unchanged` |
| 输入契约 | `input_contract.json` **15 项输入**全部 `status=verified`（含 path/size/SHA-256）；GTDB taxonomy/metadata/tree **三项均 `pending`**（未伪造哈希）；总体 `status = pending`（因 GTDB 项缺失，符合 `run_context` 语义） | `input_contract.json` |
| 服务器 | **只读查询**，未写入任何文件、未建 deploy、未跑计算；快照 `load 5.22/4.80/4.89`、80 核、可用内存 969 GB、`runs/` 1016G、`/home/data` 余 100T | `logs/server_resources_snapshot.txt` |

> **日志编码披露：** 本机 `Tee-Object` 写出的日志为 **UTF-16LE + BOM**，直接 grep 会失败。按项目既有惯例（保留原始日志不覆盖），本 run 另存 **UTF-8 镜像** `logs/utf8_mirrors/`（12 个文件，编码转换、文本等价），两侧 SHA-256 均记录在 `results/acceptance_summary.json` → `utf16_log_mirrors`，便于后续会话直接检索。
> **脚本自纠记录（如实保留）：** 本轮修了 3 处**本 run 自身**的缺陷 —— (a) `00_init_run.py` 因编辑器已建 `logs/` 而触发 fail-closed 复用告警（已改为"仅当 `logs/` 外无产物时才容忍"）；(b) `04_build_amendment_manifest.py` 把 `bytes` 当 `Path` 用；(c) `07_write_acceptance_summary.py` 的行尾/引号解析（CRLF 未归一、Python repr 非 JSON）。三处均为**证据读取/生成**缺陷，不涉及任何科学判据或数据。

---

## 5. 待操作者确认（3 项，详见 `results/pending_operator_confirmation.md`）

1. **commit/push 条款（最重要，须一行回复）** —— 选定 A（本 run 推荐：本地 commit 允许+记录、push 需授权）/ B（完全恢复被删原句：提交与推送都需授权）/ C（本地与远端均允许）。
2. **发现层边界表述** —— 是否同意把该能力从一次性授权升格为常设能力；是否要求 HMM 文件 `NAME` 行也逐字为 `discovery_hmm_uncalibrated`（**当前实测为 `cys_discovery_uncalibrated`，若要求则未满足，需另行授权重建**）。
3. **归档策略** —— 是否允许"移动而非删除"；归档目录方案未定（本 run **未创建**任何归档目录），且移动会令以原路径为输入的 manifest 失效，是否要求同时产出"原路径 → 归档路径"映射清单（本 run 建议需要）。

---

## 6. 明确不做 / 未解决项

- **未**执行 `git commit` / `push` / `clean` / `reset`；**未**触碰任何历史 run、`results/`、`pipeline/config/`、历史 HMM；**未**在服务器写入任何文件。
- **未**重算 78.69% 混淆命中率（**引用 + 文件哈希复核**，非重跑 `hmmsearch`）；如需重算须单独授权服务器任务。
- **未**测定两张 MAFFT 比对的生物学等价性（故规则只要求"显式声明 + 绑定具体哈希"，不设容差阈值）。
- **未**评估工作区的公开脱敏边界（这正是 push 需操作者确认的理由之一）。
- **未**拆分统计服务器 `runs/` 1,016 GB 的逐 run 明细。
- **并发写入如实披露：** 本 run 执行期间另有一 agent 在本工作区写 `pipeline/scripts/` 与 `pipeline/tests/`（19:42:45–19:45:13 实测 mtime），导致**不过滤的完整套件当前为红**。该失败**不属本 run 的回归**（归因见 §4）。本 run 未修改也未回滚其任何文件。

---

## 7. 复现命令

```powershell
# 1) 修正前状态与条款依赖扫描
python runs\20260917_phaded_agents_md_amendment_01\logs\01_scan_clause_dependencies.py
# 2) 规则回归测试（对当前 AGENTS.md 应为 Ran 10 tests OK）
python -m unittest pipeline.tests.test_agents_md_rules -v
# 3) 剔除并发在写模块后的完整套件（应为 Ran 534 tests OK, skipped=1）
python runs\20260917_phaded_agents_md_amendment_01\logs\06_run_suite_excluding_inflight.py
# 4) 验收汇总 + 产物清单复验
python runs\20260917_phaded_agents_md_amendment_01\logs\07_write_acceptance_summary.py
python runs\20260917_phaded_agents_md_amendment_01\logs\05_verify_artifacts.py
# 5) 既有验收命令
python -m unittest discover -s pipeline/tests
python -m compileall -q pipeline
git diff --check
git rev-parse HEAD
```

> 第 2 步的计数会随**工作区被再次修改**而变化；`04_build_amendment_manifest.py` 与 `05_verify_artifacts.py` 会重写本 run 自身产物（不改历史 run）。按项目惯例，**本文档不写自身 SHA-256**（其当前值记录在 `logs/artifact_manifest_verification.json` 的 `external_artifacts` 块中）。

---

*本文档为 candidate-only 记录。本 run 只改规则与测试，不产生任何生物学结论；`AGENTS.md` 修正后的全部条款仍要求：HMM/profile/domain/motif/定位/系统发育证据只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
