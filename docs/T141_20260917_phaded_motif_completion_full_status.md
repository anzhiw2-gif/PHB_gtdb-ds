# T141 20260917 全库 motif 补全层重跑（Task 7）状态文档

**Run：** `runs/20260917_phaded_motif_completion_full_01/`
**日期：** 2026-09-17
**任务：** 计划 v2 Task 7（Step 7.1–7.4，本地，无服务器）
**依赖：** Task 2 `runs/20260917_phaded_reference_residue_mapping_01/`（已完成的参考映射与检测器）
**输入（只读，未修改）：** `runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv`（109,087 行，64,498,737 B，SHA-256 `796cb0a3…8569a`，实测一致）、`runs/20260911_phaded_pfam_architecture_01/inputs/candidate_union.faa`（109,087 条，46,119,635 B，SHA-256 `64895bc7…b4cc`，与 `runs/20260909_phaded_candidate_mapping_01/inputs/candidate_union.faa` 实测同文件）
**边界：** 本轮仍为 candidate-only。所有序列模式、结构域结合与参考映射只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。本轮未 ssh、未连接服务器、未提交服务器任务、未 commit/push、未修改 `pipeline/config/formal_scan_models.tsv`（实测 SHA-256 仍为 `8a1c0508…c5b7a3`）、未修改任何历史 run 或历史 HMM。

**结论一句话：** 重算本身完成且可审计；但补全后**只有 `lipase_box`、`oxyanion_hole`、`ahsmg`、`lid` 四条判据在 109,087 条全库上完全可判定**，`his`/`asp` 只对模式命中的少数行可判定，`linker`/`sbd` **在候选层一条都不可判定**。这不是执行失败，而是 Task 2 参考映射边界在全库上的真实投影。

---

## 1. 为什么候选层只能是"序列模式层"

Task 2 用 **DED 比对列传递**（alignment column transfer）把文献坐标传给参考面板，这是它把 His 从 0/723 提到 90/723 的唯一机制。该机制的有效范围是**比对成员资格**。实测（本轮脚本 `candidate_alignment_members()`，只读解析 `runs/20260909_phaded_reference_evidence_01/inputs/raw/ded/` 下全部 `*.aln`，同目录另有 39 个 `.html` 与 46 个 `.txt` 抓取残页，**不当作比对解析**）：

| 实测项 | 数值 |
|---|---:|
| 候选序列（`candidate_union.faa`） | 109,087 |
| 与 723 条参考 ledger 的 accession 交集 | **0** |
| 解析的冻结 DED 比对文件（`*.aln`） | 38（其中 29 个含序列行、9 个 header-only，与 Task 2 报告一致） |
| 与全部比对行名的交集 | **0** |
| Task 2 报告中取得锚定列传递的参考数 | **65**（参考面板内部；任务书写的 64 为旧值） |

**因此参考锚定坐标一条也落不到候选行。** 本轮的全部候选层证据只有三类：①文献校准的序列模式；②已冻结的 Pfam 结构域结合（PF06850/PF10503）；③记录在 0915 面板里的 superfamily 先验（**只用于判定"该判据在此行是否可判定"，绝不作为检测器输入**）。

## 2. 方法：可判定性（assignability）是唯一新增的语义层

每条候选 × 每条判据输出 `state` + `assignability` + `assignability_reason` + `assignability_basis`。

**规则（全部由实测数据推导，脚本内无硬编码阈值）：**

1. **检测器阳性 ⇒ 一律 assignable。**
2. **检测器阴性/miss ⇒ 只有在该模式对"所有带坐标的文献锚点"都可复现时才 assignable。** 实测（`criterion_miss_assignability_basis`，脚本重扫文献锚点序列得出）：

| 判据 | 该角色带坐标的文献锚点数 | 模式可复现数 | 参考面板阳性对照数 | miss 可判定？ |
|---|---:|---:|---:|---|
| `lipase_box`（Gx₁Sx₂G） | 4 | **4** | 444 | 是 |
| `oxyanion_hole`（HGCXQ） | 2 | **2** | 256 | 是 |
| `his`（GMxH） | 4 | **1** | 90 | **否** |
| `asp`（GxxDYTV） | 4 | **1** | 71 | **否** |
| `ahsmg`（AHSMG 字面五肽） | 0（无坐标） | – | 2 | 是（有已验证阳性对照） |
| `lid`（FNGIG+YYWQLF 双环 + 间距 156±8） | 0 | – | 39 | 是（有阳性与阴性对照） |

   His 只复现 1/4、Asp 只复现 1/4 —— 因为 PhaZKT His248（DDGHL）、PhaZGK13 His260（GGHEW）、2D80 His135（AVHTF）与两者的 Asp 上下文都不含该模式。**模式不命中不代表该残基不存在，只代表该代理模式覆盖不了这个上下文**，所以这些行的 His/Asp 必须写 `not_assignable_*`，绝不能写 `not_detected_pattern`。
3. **`lid` 对每一行都 assignable**：检测器只读序列（`detect_lid(sequence, signature=None)`，源码中不含 subtype/superfamily/family 字样，已由测试钉住），先验单独放在 `lid_expectation_from_prior` 列，可以拿检测结果去对先验，而不是复述先验。
4. **`linker`、`sbd` 对每一行都 not_assignable**：linker 无任何已绑定 profile（只能"划定候选区间"，那不是 linker 判定）；PF06850 在本面板不是胞外 SBD（Task 2 DEFECT_D：274/276 命中 intracellular nPHASCL without lipase box，在 284 条胞外 type 1 与 71 条 type 2 参考中 0 命中）。PF06850 的原始结合结果仍按行发布。

## 3. 结果：每条判据的可判定 / 不可判定分解（109,087 条全库）

| 判据 | 0915 候选层 before | 本轮 after | assignable | not_assignable | 判定 |
|---|---|---:|---:|---:|---|
| `lipase_box` | supported 77,302 / not_detected_pattern 31,785 | supported 77,302 / not_detected_pattern 31,785 | **109,087** | 0 | 完全可判定（模式未变） |
| `oxyanion_hole` | supported 18,902 / not_detected 89,812 / conflict 373 | supported 18,939 / not_detected_pattern 90,148 | **109,087** | 0 | 完全可判定 |
| `ahsmg` | not_detected_pattern 109,087 | not_detected_pattern 109,087 | **109,087** | 0 | 完全可判定（0 阳性） |
| `lid` | not_expected_for_subtype 101,836 / pending 7,251；**独立检出 0** | supported 989 / loop1-only 1,910 / loop2-only 11 / not_detected_pattern 106,177 | **109,087** | 0 | 完全可判定（2,910 条至少检出 1 个环） |
| `his` | supported **265** / not_detected_pattern 108,822 | supported **7,317** / not_assignable_* 101,770 | **7,317** | 101,770 | 部分可判定 |
| `asp` | supported **1,290** / not_detected_pattern 107,797 | supported **1,290** / not_assignable_* 107,797 | **1,290** | 107,797 | 部分可判定 |
| `linker` | pending 108,795 / not_expected 292 | not_assignable_no_linker_profile_bound 109,087 | **0** | 109,087 | **全库不可判定** |
| `sbd` | supported 30,621 / not_detected 67,652 / pending_tool_input 10,522 / not_expected 292 | not_assignable_pf06850_is_not_the_extracellular_sbd_in_this_panel 109,087 | **0** | 109,087 | **全库不可判定** |

每条判据的 assignable + not_assignable 均严格等于 **109,087**（测试逐条核对，并从 TSV 独立重算一遍，不信任 summary）。

### 不可判定原因分解（`results/coverage_decomposition.tsv`）

| 判据 | 原因代码 | 行数 | share |
|---|---|---:|---:|
| his | `not_assignable_candidate_not_in_any_ded_alignment` | 7,236 | 6.63% |
| his | `not_assignable_circularly_permuted_fold_column_not_transferable` | 2,646 | 2.43% |
| his | `not_assignable_no_literature_anchor_for_prior_superfamily` | 60,495 | 55.46% |
| his | `not_assignable_no_superfamily_prior_and_no_alignment_membership` | 31,393 | 28.78% |
| asp | 同上四类 | 7,249 / 2,640 / 66,274 / 31,634 | – |
| linker | `not_assignable_no_linker_profile_bound` | 109,087 | 100% |
| sbd | `not_assignable_pf06850_is_not_the_extracellular_sbd_in_this_panel` | 109,087 | 100% |

四类原因的含义：①先验有文献锚点（extracellular dPHAMCL、intracellular nPHAMCL），但候选不是任何 DED 比对成员，列传递到不了；②先验有锚点但锚定列该折叠不可传递（type 2 环状置换，Task 2 实测保守度 0.0769）；③先验根本没有带坐标的文献锚点（extracellular dPHASCL type 1 = 16 个家族，以及 intracellular nPHASCL 有/无 lipase box、periplasmic、native-SCL/PhaZ7-like）；④记录的先验为空（0915 的 `unassigned_PhaDED_like`）。

## 4. 关键量化：补全 His 与 Asp 模式后，全库到底多了多少可判定条目？

- **His：265 → 7,317 assignable（+7,052，×27.6）**。0915 的 265 条阳性**全部保留**（GM[A-Z]H[A-Z]{2}P[A-Z]{2}G 是 GMxH 的超集，实测 265/265 仍是 supported）。但 101,770 条（93.3%）仍不可判定。
- **Asp：1,290 → 1,290 assignable（+0）**。Task 2 的 Asp 正则 `G[A-Z]{2}DYTV` **与 0915 逐字节相同**，所以候选层一条也没多。Task 2 参考层的 25 → 71 全部来自锚定列传递，而该机制到不了任何候选（第 1 节）。**这是本轮最重要的负面结论：参考层的补全不会自动转化为候选层的覆盖率。**
- **His + Asp 合计：1,555 → 8,607 assignable（+7,052，占全库 7.89%）**；至少一条不可判定的行 100,480 条（92.11%）。两条都可判定者仅 1,066 条。
- **lid：独立检出 0 → 989 条双环 supported（+1,910 loop1-only、+11 loop2-only，共 2,910 条至少一个环）**，且 109,087 条全部获得可判定状态。与先验对照：detected_and_expected 989 / detected_but_not_expected 0 / not_detected_but_expected 5,970 / not_detected_and_not_expected 102,128。
- **linker / sbd：可判定条目 0 → 0**，本轮没有、也不可能产生任何一条可判定的 linker 或胞外 SBD 结果。
- **`oxyanion_hole`：18,902 → 18,939（+37）**；0915 的 373 条 `conflict_relative_position` 在本轮**不再产生**，因为那个几何检查用 subtype 先验当输入，几何归属属 Task 4 run（`type1_verified` 16,530 / `type2_verified` 2,362 / conflict 373）。

### 需要显式纠正的三处前提（任务书/计划文本 vs 实测）

1. **"his 全库 0 命中"不成立于候选层。** 0/723 是**参考层** `reference_motif_panel.tsv` 的数字；0915 **候选层** recorded `his_state=supported` 共 **265** 条（被 0915 那个自相矛盾的模式 GM[A-Z]H[A-Z]{2}P[A-Z]{2}G 命中）。
2. **"asp 全库约 25 命中"同样是参考层数字。** 0915 **候选层** `asp_state=supported` 为 **1,290** 条。
3. **"仅 64/723 条参考可传递坐标"**：Task 2 报告实测为 **65**；且它是参考层数字，候选层的对应值是 **0**。

## 5. "部分催化模式面板"与"完整催化三联体证据"的区别（必须分开表述）

- 本轮整库的 `motif_panel_status` 分布：`partial_catalytic_pattern_panel` 21,512 / `partial_lipase_box_only` 56,025 / `no_catalytic_residue_pattern_detected` 31,550。**没有任何一行达到 `complete_*`。**
- `partial_catalytic_pattern_panel` 只表示"**在序列上同时命中了亲核体模式与至少一条其他催化模式**"，即**部分催化模式面板**。它**不**等于"催化三联体已确认"，因为：(a) Ser/His/Asp 三者从未在同一条候选上以**带坐标的参考锚定**被确认（候选层没有锚定传递）；(b) His/Asp 的阳性仅是代理模式的命中，且这两个模式各自只覆盖 4 个文献锚点中的 1 个；(c) AHSMG 全库 0 命中，Cys 型亲核体坐标不可得。
- 因此正确表述是"**候选层序列模式面板（partial/partial-catalytic）**"，而不是"催化三联体完整"。任何下游若需要"完整催化三联体证据"，必须另开任务（例如把候选并入 DED 家族比对再传递锚定列，或引入结构证据），**本轮不提供**。
- 与 0915 词汇表的差异（已刻意标注，不是静默改名）：① 0915 的 `not_tested_full_library` 桶被替换为 `no_catalytic_residue_pattern_detected` —— 0915 说"全库未测"，本轮确实测了全库，继续用那个词会失真；② 0915 的 `conflict_relative_position` 移交 Task 4（几何轴）；③ 新增的 `not_assignable_*` 家族是本轮唯一新增状态词，**只用于不可判定的行**，与 `supported` / `not_detected_pattern` / `pending_reference_annotation` 互斥且不重叠。

## 6. 不可判定不等于阴性（本轮硬不变量）

- 全部 109,087 × 8 条判据逐格检查：`assignability=not_assignable` 的行，其 `state` 必须以 `not_assignable` 开头，且**不得**含 `not_detected_pattern` / `not_detected_in_tested_pfam` / `not_detected_at_mapped_column` / `not_expected_for_subtype` / `not_expected_for_superfamily`（测试逐行强制，并从产物 TSV 独立复扫一遍）。
- 换言之：**His/Asp 有 101,770 / 107,797 行、linker 与 sbd 各有 109,087 行是"我们不知道"，不是"没有"。** 这些行**不得**被任何下游当成阴性、不得用于计数分母、不得据此删候选。
- `sbd` 的 30,621 → 0 是**降级**（把 PF06850 原始结合撤回为候选级证据），**不是证据丢失**：原始结合在本轮 `sbd_pf06850_binding_state` 中仍是 detected 30,621，与 0915 完全一致。同理，0915 `pending_tool_input` 10,522 与本轮 `not_tested_no_pfam_evidence` 10,668 的差 146，来自 0915 把 292 条 extracellular dPHAMCL 强行写成 `not_expected_for_subtype`（其中 146 条其实没有任何 Pfam 结果），本轮按数据如实发布。

## 7. 本任务不产生任何新的 family 判定

- 每行写入 `new_family_call_made=false`，摘要写入 `new_family_or_superfamily_call_made=false`，测试强制。所有 superfamily 标签均来自 0915 已冻结的 `motif_reference_subtype`，只用于第 2 节的可判定性分支，**不构成新判定、不改动任何 family 归属、不改动任何 profile 或 registry 资格**。
- 一个**观察（不是判定）**：0915 的 `unassigned_PhaDED_like` 共 31,634 条，其中 29,906 条带 PF06850、30,700 条无 lipase box、仅 1 条带催化域 PF10503；而 67,544 条 `extracellular dPHASCL type 1` 中只有 27 条带 PF06850。这与 Task 2 的 DEFECT_D（PF06850 标记的是胞内 Cys 型而非胞外 SBD）方向一致，提示该群体值得一个**单独授权的**任务去核对（例如 Task 6 的 Cys 型 decision record）。**本轮不把它写成任何 family 或 superfamily 判定。**

## 8. 交付物与校验

| 产物 | 说明 |
|---|---|
| `results/motif_completion_full.tsv` | 109,087 行；每行 8 条判据的 state / assignability / reason / basis，加 prior 与 PF06850 原始证据列 |
| `results/motif_completion_summary.json` | 每条判据 before/after 分布、assignable/not_assignable、原因分解、前提纠正、prior 与锚点索引、PF06850 普查、输入哈希 |
| `results/coverage_decomposition.tsv` | 每判据 × 桶（assignable / 每个 not_assignable 原因）的行数与占比 |
| `logs/step7.1_failing_test_red.log` | 失败测试先行证据（模块缺失即红） |
| `logs/step7.3_focused_tests.log` | 聚焦测试通过日志 |
| `logs/full_test_suite.log` | 完整测试套件日志 |
| `logs/artifact_sha256_manifest.txt` | 10 个产物（3 个 results、input_manifest、runner、commands、input_contract、脚本、测试、本文档）的 size/SHA-256 |
| `inputs/input_manifest.tsv` | 全部外部输入的 path/size/SHA-256；GTDB taxonomy/metadata/tree 三项如实写 `pending` |
| `logs/commands.txt` | 完整命令行（本地执行，无服务器） |
| `input_contract.json` | run_context schema；GTDB taxonomy/metadata/tree 三项如实写 `pending` |

**验收命令（本地，最终一次合并运行）：**

| 命令 | 结果 |
|---|---|
| `python -m unittest discover -s pipeline/tests` | **Ran 450 tests … OK (skipped=1)**，exit 0 |
| 本任务聚焦测试 `pipeline.tests.test_complete_phaded_motif_completion_full` | **Ran 29 tests … OK**（failures 0 / errors 0） |
| `python -m compileall -q pipeline` | exit 0 |
| `git diff --check` | exit 0 |

450 = Phase 1 基线 375 + 本任务 29 + 其他并行任务的测试；本任务相关的 29 条全绿。

## 9. 并行工作与遗留（诚实披露）

- 本轮第一次跑完整套件时（约 01:2x），`pipeline/tests/test_design_phaded_discrimination_panel.py`（**属于并行工作的 Task 5 agent**，其目标脚本 `pipeline/scripts/design_phaded_discrimination_panel.py` 当时尚未创建）报出 31 个 `FileNotFoundError`；随后该 agent 补齐脚本，短暂出现 5 error + 1 failure。**这些都是他人的进行时产物，我只报告、未修改。** 最终一次合并运行（第 8 节）已全绿，上述中间态未对本任务任何产物造成影响。
- **仍然不可判定的判据：** `his`（101,770 行）、`asp`（107,797 行）、`linker`（109,087 行）、`sbd`（109,087 行）。阻塞条件与 Task 2 相同：候选不是 DED 比对成员、四个超家族无带坐标的文献锚点、type 2 折叠环状置换、linker profile 未绑定、胞外 SBD 判据未解决、Cys 型亲核体坐标不可得。
- **本轮状态判定：** 重算本身**完成**（109,087 条逐行、可复算、可审计），但五条目标判据中有两条（linker、sbd）在候选层**根本无法判定**、两条（his、asp）仅部分可判定。按计划 Step 7.1–7.4 的字面要求本轮为 **partial**：不是执行失败，而是"整库重算完成 + 部分判据结构性无法执行"的如实结果。
- **未做（越界）：** 未执行 GTDB 全量重扫、未改 profile/registry、未做任何 family 合并或删除、未把任何 `not_assignable` 行当阴性使用、未 ssh、未触碰服务器。

---

*本文件为 candidate-only 记录：所有 motif、结构域、参考映射与几何结果只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
