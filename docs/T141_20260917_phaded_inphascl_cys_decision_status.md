# T141 — `i-nPHASCL (no lipase box)` Cys 型判据 decision record 与候选层证据：状态文档

**任务：** 计划 v2 Task 6（`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` §3 Task 6）
**Run：** `runs/20260917_phaded_inphascl_cys_decision_01`（新建；未触碰任何历史 run、`results/`、`pipeline/config/`、AGENTS.md、历史 HMM）
**Decision record：** `docs/superpowers/plans/2026-09-17-phaded-inphascl-cys-decision-record.md`
**日期：** 2026-09-17
**状态：** `completed_candidate_only`
**服务器：** 未使用（0 线程、0 任务、未建 dated `deploy/`）

---

## 0. 一句话结论

推荐并已执行**方案 A**：registry 与 `pipeline/config/` **完全不变**、**不训练任何新 HMM/profile**，改为发布一层可复核、可证伪的**候选层注释证据**。为 Cys 型超家族建立 fitted profile（方案 B）**被项目硬规则禁止**（276 条参考中 275 条为 `annotation_only`，唯一可严格训练的记录 1 条 < gate 要求 3），因此**即使已获授权也不应执行**。

---

## 1. 本任务不产生任何新的 family 判定（显式声明）

| 声明 | 事实依据 |
|---|---|
| **不产生新的 family 判定** | 实现脚本不写任何 family/superfamily 列；109,087 行每行 `family_call_made=false`、`new_family_call_made=false`；manifest `new_family_call_made=false`、`family_call_made=false` |
| **不训练任何模型** | `hmm_fitted=false`、`hmm_score_present=false`；run 内无任何模型文件（`results/` 只有 2 张 TSV + 3 个 JSON + 1 张探索性 TSV） |
| **不改校准 gate** | gate 常量与判据未触碰；`positive_count>=3`、`held-out>=1`、`family-resolved negative>=1`、`challenge>=1` 全部保持原样且**仍未达成** |
| **不改 registry 权威层级** | 8 个 superfamily 仍 `functional_prior / registry_eligible=true`；38 个 family 仍 `sequence_clustering_2009 / registry_eligible=false` |
| **不启动 GTDB 全量重扫** | 未运行任何 hmmsearch/pfam_scan/profile 扫描 |

---

## 2. 对 109,087 条候选的 `subtype_call` 影响为 **0 行**

- 本任务**只读**消费冻结证据（Task 7 全库 motif 完成表、Task 2 参考映射表、候选 FASTA、参考 FASTA），**未写入**任何既有文件（实现脚本对输入只读；测试 `test_end_to_end_outputs_are_self_consistent_and_inputs_untouched` 断言四个输入的字节在运行前后完全相同）。
- 产出表 `results/cys_candidate_evidence.tsv` 与 `results/cys_typing_descriptor_evidence.tsv` 都是**只新增证据列**的表，不含 `subtype_call` 列，且每行显式写 `subtype_call_impact=none`。
- `results/cys_decision_manifest.json` 显式记录 `subtype_call_rows_changed: 0`；探索性 pass 的 manifest 同样记录 `subtype_call_rows_changed: 0`。
- 依据（计划 §1.1）：subtype 判定只由 `phaded_superfamily_best` 决定（`build_phaded_subtype_matrix.py` L153–165），本任务既不产出也不消费任何 `*_best` 字段作为输入。

---

## 3. 候选层注释证据 vs fitted HMM 证据：必须分清的区别

| 维度 | 本轮发布的东西 | 一个 fitted HMM/profile 会给出的东西 |
|---|---|---|
| 证据类型 | Pfam 域绑定（PF06850 原始结合）+ 成熟序列模式（`Gx1Sx2G` 及其缺失）+ 序列上下文（cysteine−1/+1） | 概率模型的 log-odds 打分与 E-value |
| 是否有模型参数 | **无** | 有，且必须版本化绑定 |
| 参考层对比强度 | PF06850：276 条 Cys 参考命中 274，其余 447 条命中 0（比值 0.9928 vs 0.0000） | 由 held-out 集决定 → 当前 0 |
| 能否进 registry | **不能**，本任务也不试图进 | 可被考虑，但必须先通过 gate |
| 能否支撑"催化完整性"表述 | **不能**。他的/Asp 在候选层仍是模式层（His 可判定 7,317 / 不可判定 101,770；Asp 可判定 1,290 / 不可判定 107,797），不是三联体证据 | 同样不能单独支撑 |
| 能否作为表型证据 | **不能** | 不能 |
| 能否据此删除候选 | **不能** | 不能 |
| 失败可逆性 | 高：删掉新 run 即可完全回退，无任何下游依赖 | 低：一旦进 registry 会污染版本权威 |

**边界表述（可直接引用）：** 本轮产物是"候选层注释证据"，表示候选同源或功能潜力，**不是 family 判定、不是 HMM 打分、不是已验证的 PHB/PHA 降解表型**。

---

## 4. 交付物

| 产物 | 大小 (B) | SHA-256 |
|---|---:|---|
| `results/cys_candidate_evidence.tsv` | 42,397,093 | `3cf1edaeb73a3fd00c585105d66b70d237a1707c0e65944a39ad040c76e47c4d` |
| `results/cys_reference_contrast.tsv` | 909 | `2a14e7b6d0d190feb4786b964730617ad021c4f5e2d89492f3cd14a72bc4121f` |
| `results/cys_decision_manifest.json` | 6,945 | `11f6607d7564df85f5bb58cabf173c03c1ec7b3a299b5634cc946ad33c4edfcc` |
| `results/cys_typing_descriptor_evidence.tsv` | 34,288,167 | `ad0017c26a56250fce874fd577db90ec6576ca2b41e3b2d795d7d83235bf89ff` |
| `results/cys_typing_descriptor_discovery.json` | 16,371 | `e019f7e7487d85678d62c80d4e45640a675adc7e8f68575a9a6425981302e9d4` |
| `results/cys_typing_descriptor_manifest.json` | 8,435 | `2427da482444447a695c7bc882697710416bc99c33412b30b94b4c547282a6f0` |
| `input_contract.json` | 3,421 | `a9943ce389a14f97965a0455d263db3279ce8cd5dcf80b459831054df6694fa3` |
| `inputs/cys_criterion_definitions.tsv` | 954 | `0c7ae80b291aa4f23680924e625da503fe3e26c083d573dbf7648645eff8f147` |
| `inputs/falsification_preregistration.tsv` | 1,827 | `6ce352be129cc8574d9cbb0d433288bf4a07eddb269991b4dccc14b3760feb2b` |
| `inputs/input_manifest.tsv` | 1,068 | `89e62a7eb173d2ef22c1a8c4d0688b40feef269a0f3b8ab80e1fddbd4e7a94b3` |
| `logs/artifact_sha256_manifest.txt` | — | 覆盖 run 内全部 22 个文件（本文件自身除外） |
| `pipeline/scripts/annotate_phaded_inphascl_cys_evidence.py` | 新文件 | 实现（不拟合任何模型） |
| `pipeline/tests/test_annotate_phaded_inphascl_cys_evidence.py` | 新文件 | **18 个测试**：其中 **15 个属"先写后实现"的红灯测试**（有落盘日志）；**另 3 个属第二次探索性 descriptor pass 后补写的测试**（红态在会话内被观察到但**未落盘日志**，证据缺口见 §11） |

run 目录含 `logs/`、`inputs/`、`results/` 与 `input_contract.json`；`run_id` 经 `pipeline/scripts/run_context.py::validate_run_id` 校验，契约由 `run_context.write_input_contract` 写出；GTDB `taxonomy` / `metadata` / `tree` 三项均写 `pending`（**未伪造哈希**）。`inputs/input_manifest.tsv` 记录全部 5 个输入与 3 项 GTDB 的 path/size/SHA-256 或 `pending`。

---

## 5. 证据与计数

### 5.1 参考层（723 条冻结参考；独立重算，不调用上游脚本）

| 超家族 | 参考数 | PF06850 命中 | 无 lipase box | `V-C` | `VCQ` |
|---|---:|---:|---:|---:|---:|
| **intracellular nPHASCL without lipase box** | 276 | **274 (0.9928)** | 273 | 256 (0.9275) | **251 (0.9094)** |
| extracellular dPHASCL type 1 | 284 | **0** | 0 | 65 | **0** |
| extracellular dPHASCL type 2 | 71 | **0** | 2 | 18 | **0** |
| 其余 5 个超家族（dPHAMCL / PhaZ7-like / nPHAMCL / with lipase box / periplasmic） | 92 | **0** | 4 | 10 | **1** |

- 唯一 1 条携带 `VCQ` 的非 Cys 参考是 `EDL63834.1`（`intracellular nPHASCL with lipase box` / `DED_hfam_2` / `annotation_only`），属相邻超家族预期交叉。因此全部 7 个非 Cys 超家族合计 `VCQ` = **1/447（0.0022）**。

### 5.2 候选层（109,087 条，100% 覆盖，0 行遗漏）

| 项 | 计数 |
|---|---:|
| PF06850 `detected` / `not_detected_in_tested_pfam` / `not_tested_no_pfam_evidence` | 30,621 / 67,798 / 10,668 |
| lipase box `supported` / `not_detected_pattern` | 77,302 / 31,785 |
| **`PF06850 detected` ∧ 无 lipase box（声明的判别证据）** | **29,949** |
| └ 叠加 `V-C` | 26,023 |
| └ 叠加 `VCQ`（探索性） | 25,061 |
| 无超家族先验中满足判别证据 | 29,271 |
| `cys_evidence_class` 12 类计数之和 | 109,087 ✅ |

### 5.3 关键治理事实（本轮经核实）

- **Cys 超家族恰有 1 条 `experimental_positive`**：`CAJ92291.1`（`DED_hfam_65`，`Cupriavidus necator`，`positive_seed`，Knoll 2009 Table 1 seed GI 3641686）—— 即任务书点名的模型酶。任务书"hfam_61 到 hfam_69 实验阳性为 0"的说法因此**需要更正**（详见 decision record §2.1）。该记录是文献策展结果，不是本项目新做的实验。
- 其余 275 条参考均为 `annotation_only` → 被 `validate_phaded_reference_ledger.py` L172–173 拒绝进入严格训练。
- 候选层与参考层零重叠：0/109,087 条候选属任何冻结 DED 比对，0 条出现在 723 行台账中（Task 7 已固化，本任务沿用）。

---

## 6. 三个方案与推荐（摘要；完整论证见 decision record §3–§5）

| 方案 | 治理可行性 | 证据充分性 | 结论 |
|---|---|---|---|
| **A. registry 不变 + candidate-only 说明（+ 本轮候选层证据层）** | ✅ 完全合规 | 对"注释证据"充分；对"family 判定"不充分 | **推荐并已执行** |
| B. 新增窄范围 Cys profile 且保留历史模型 | ❌ 违反 annotation_only 禁训规则；gate 不满足 | ❌ 仅 1 条可训练序列 | **不执行，即使获授权** |
| C. 版本化替换旧模型 | ❌ 双重禁止（继承 B + 破坏历史权威） | ❌ 最低，且无替换对象 | **不执行** |

**不需要服务器计算**（§9）。

---

## 7. 预注册证伪条件与实测结果

条件在运行前写入 `inputs/falsification_preregistration.tsv`，由 `evaluate_falsification_gates` 做 fail-closed 评估：任一条 `failed` → `status=falsification_triggered`、全部判别置为 `withheld_falsification_triggered`、进程返回 3。

| ID | 条件 | 实测 | 结果 |
|---|---|---|---|
| F1 | PF06850 不得结合 Cys 超家族以外的参考 | 超家族外命中 **0** | ✅ passed |
| F2 | PF06850 在 Cys 超家族覆盖率 ≥ 0.90 | **274/276 = 0.992754** | ✅ passed |
| F3 | **284 条 type 1 + 71 条 type 2 胞外参考中不得出现未被解释的判别命中** | **0**（type 1: 0/284；type 2: 0/71） | ✅ passed |
| F4 | lipase box miss 必须可归因（4/4 锚点复现、`miss_is_assignable=true`） | `anchors_with_coordinate=4`、`anchors_reproduced=4` | ✅ passed |
| F5 | 不得引入任何 fitted HMM/profile | `hmm_fitted=false`，无模型文件 | ✅ passed |
| F6 | 候选层命中只能是注释证据、不可升格为 family 判定 | 列级强制（见 §1） | ✅ passed |

**F7（本轮事后登记，未预注册，属下一轮条件）：** 若在独立非 Cys α/β-水解酶阴性集上 `VCQ` 命中率 > 5%，则该描述符必须降级为纯描述列。`VCQ` **不得**在 F7 通过前被用作过滤器。

---

## 8. 探索性描述符与诚实限定

在预注册 pass 全部通过、产物落盘之后，本轮又做了一次**明确标注为探索性**的第二次 pass（只写新文件名，不覆盖任何既有产物）：

- 发现：由 Knoll 2009 的 cysteine−1 规则出发，观察 cysteine+1，得到三肽 **`VCQ`**。样本内对比 251/276（Cys 超家族）vs 1/447（其余超家族）。
- **诚实限定：** `VCQ` **没有文献依据**（Knoll 2009 只写了 cysteine−1），是**从冻结参考台账本身 data-derived 出来的**，其敏感度/假阳性率是**样本内**数字；工作区内**不存在 held-out 集**（`in_sample=true`、`held_out_set=none`）。因此它被标为 `exploratory_in_sample_descriptor_discovery`、`load_bearing=false`，**不参与判别、不作为 gate、不进入 registry 论证**。
- 第二个 pass 的 manifest 显式记录**预注册产物 `cys_candidate_evidence.tsv` 的 SHA-256**（`3cf1eda…`），使审阅者可验证"预注册结果未在事后被调整"。`input_contract.json` 由第一次 pass 写出后**刻意未被第二次 pass 改写**（第二次 pass 的 provenance 记在自己的 manifest 与 `logs/descriptor_pass_commands.txt` 中）。

---

## 9. 服务器

**本任务未使用服务器。** 推荐方案不需要比对打分、不需要 hmmscan、不需要任何服务器算力（两次 pass 全本地，合计约 2 分钟）。

只读负载核对（未提交任何远程任务，占用 **0 线程**），`2026-09-17T01:46:02+08:00`，`<SERVER_USER>@<SERVER_HOST>`：

```
2026-09-17T01:46:02+08:00
 01:46:02 up 6 days, 10:26,  6 users,  load average: 4.50, 4.77, 4.81
80
               total        used        free      shared  buff/cache   available
Mem:            1007          50          77           0         885         956
0, 0 %, 1 MiB
1, 0 %, 1 MiB
```

未创建 `deploy/20260917_phaded_inphascl_cys_decision_01/`，因为**没有需要执行的服务器源码**。

---

## 10. 验收命令

| 命令 | 结果 | 日志 |
|---|---|---|
| `python -m unittest discover -s pipeline/tests` | **Ran 484 tests — OK (skipped=1)**，exit 0（两次独立运行结果一致） | `logs/full_test_suite.log`、`logs/full_test_suite_final.log` |
| `python -m compileall -q pipeline` | exit 0 | `logs/compileall.log` |
| `git diff --check` | exit 0（无空白错误） | `logs/git_diff_check.log` |
| `python -m unittest pipeline.tests.test_annotate_phaded_inphascl_cys_evidence` | **Ran 18 tests — OK** | `logs/step2_focused_tests_green.log` |
| 未执行 `git commit` / `push` / `clean` / `reset` | 是 | — |

**收尾核验（本任务自查）：**

| 核验项 | 结果 |
|---|---|
| `pipeline/config/formal_scan_models.tsv` | 450 B / `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3` → **未变** |
| 历史 run（`20260901_*`–`20260916_*`） | 扫描 **3,177** 个文件，mtime ≥ 2026-09-17 的 **0** 个 |
| 本次消费的 5 个冻结输入 | 运行后逐条重算 SHA-256：**5/5 MATCH**（证明只读） |
| `cys_candidate_evidence.tsv` 行数 / 唯一 accession | 109,087 / 109,087（0 重复，0 遗漏；与 Task 7 源表 accession 集合完全相等） |
| 原始证据列与 Task 7 源表比对 | `sbd_pf06850_binding_state` 与 `lipase_box_state` 逐行 **0 处不一致** |
| boundary 列全行检查 | 109,087/109,087 行满足 `evidence_layer` / `family_call_made=false` / `hmm_score_present=false` / `subtype_call_impact=none` / `new_family_call_made=false` |

**测试计数口径说明：** Phase 2 结束时的基线是 450 tests。当前 484 中有 **18 个属于本任务**（`pipeline/tests/test_annotate_phaded_inphascl_cys_evidence.py`）。其余增量来自**并行工作的其他 agent**，本任务未修改、也未评估它们的产物。`git status --porcelain` 中唯一的已跟踪改动仍是 ` M AGENTS.md`（Phase 2 已记录的既有状态，本任务未触碰）。

---

## 11. Test-first 证据与缺口披露

- **已保留的红灯证据：** `logs/step1_failing_test_red.log`（34,006 B，UTF-16LE）—— 实现文件尚不存在时运行测试模块的原始输出：`Ran 15 tests`、`FAILED (errors=15)`，全部 error 为 `FileNotFoundError: ... pipeline\scripts\annotate_phaded_inphascl_cys_evidence.py`。即 **test-first 成立**（这 15 个测试先存在且先红），15 个测试名与当前 18 个方法中的 15 个逐一相同。
- **缺口披露（如实记录）：** 后续新增的探索性描述符测试（3 个：`test_descriptor_pass_fails_closed_on_missing_sequence`、`test_descriptor_pass_requires_the_preregistered_artifact_and_leaves_it_byte_identical`、`test_vcq_positions_are_one_based_and_glutamine_specific`）的**红状态在会话内被观察到但原始日志未落盘**。观察到的报错为：
  `AttributeError: module 'annotate_phaded_inphascl_cys_evidence' has no attribute 'annotate_typing_descriptor'`（2 次）与 `... has no attribute 'val_cys_gln_cysteine_positions'`（1 次），
  随后实现补齐、`Ran 18 tests ... OK`（`logs/step2_focused_tests_green.log`）。本任务**不伪造**该红灯日志文件，只在此如实登记。
- **口径更正（`20260917_phaded_final_consistency_fix_01`，2026-09-17）：** 本文档 §4 交付物表原先把 18 个测试整体描述为「18 个测试，**先写后实现**」，其中"先写后实现"仅对**首轮 15 个**有落盘红灯证据成立（`Ran 15 tests … FAILED (errors=15)`），**后加 3 个无 test-first 红灯日志**；该表已改为按 15 + 3 分述。本 Task 的 decision record（§6 产物表）同处自述也已同步更正：它写的是「18 个测试，先写后实现」，而同一张表的上一行本就记「红灯证据（15 个测试全部因实现文件不存在而 error）」——两行原本互相矛盾。**该缺口不影响任何产物**：3 个后加测试属**未被预注册**的探索性 descriptor pass（`exploratory_in_sample_descriptor_discovery` / `load_bearing=false`），不参与 F1–F6 任何闸门。

---

## 12. 明确不做的事

- ❌ 不训练任何新 HMM/profile；不为 Cys 超家族新增 registry 条目。
- ❌ 不修改 `pipeline/config/`（尤其 `formal_scan_models.tsv`）、历史 HMM、历史 run、`results/`、AGENTS.md。
- ❌ 不降低校准 gate；不把零检出写成"不存在"。
- ❌ 不把候选层注释证据表述为 family 判定或表型证据；不据此删除任何候选。
- ❌ 不启动 GTDB 全量重扫；不提交任何服务器计算任务。
- ❌ 不删除任何文件（探索性 pass 因此只写新文件名，而非改写预注册产物）。

---

## 13. 未解决项

1. 候选层判别集**没有任何独立验证标签**（0/109,087 属 DED 比对、0 条属参考台账）→ 29,949 / 25,061 只能作候选优先级，**不得**表述为敏感度/特异度。
2. Cys 型催化亲核体坐标仍为 `pending_reference_annotation`；候选层"expected catalytic residue"依然全是 `Ser`。cysteine±1 只是**上下文描述符**，不是残基级判定。
3. `VCQ` 需独立非 Cys 阴性集评估（F7）。
4. 另外三个缺失超家族（periplasmic、i-nPHASCL with lipase box、native-SCL/PhaZ7-like）继续保持 `reference_only`。
5. 上游 gate 计数器只统计绑定到 profile 自身 family 的记录（Phase 2 报告 Task 5 指出，跨家族轴对照恒为 0）—— 不在本任务范围，**仅报告、不修改**（属其他 agent 的产物）。
6. 本轮尝试过"用 Cys 型 DED 比对的不变 Cys 列反推亲核体"这一纯本地思路，结论为**不可判定**（`aln61` 有 237/413/425 三列、`aln65` 有 120/218/411 三列，保守性无法唯一指认），已放弃且未进入任何产物；失败探查在此登记。

---

*本状态文档为 candidate-only 记录；所有 profile、domain、motif、定位与序列上下文证据只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
