# T141 — Cys 催化残基验证 / 显式降级（方案 a）：状态文档

**Run：** `runs/20260920_phaded_cys_verification_01`
**日期：** 2026-09-20
**状态：** `completed_candidate_only`
**任务来源：** 审查报告「接下来该做什么」第 2 项（Cys 催化残基验证或显式降级）

---

## 0. 一句话结论

**判定为方案 (a)：可校准的 Cys 序列级检测器存在，已实现为新的证据列。**

关键依据：Cys 型催化亲核体**确实有文献锚点**，且项目此前的冻结台账把它记为 `pending` 而错过了它——**PhaZ1 Cys183**（Ralstonia eutropha H16；定点突变显示 C183A / D355A / H388Q 失活，而 **C183S 不失活**；[PMID 16233560](https://pubmed.ncbi.nlm.nih.gov/16233560/)）。以该锚点校准的 `[疏水]-C-Q` 模式在 PhaZ1 上**只命中位置 183**，并在参考面板与候选池上都产生了清晰分离。

> **边界（必须随本列一起引用）**：这是**序列上下文匹配**，不是催化活性验证、不是 family 判定、**不得**读作"已验证 Cys"。

---

## 1. 可行性判定（为什么是 (a) 而不是 (b)/(c)）

| 判据（任务书对 (a) 的要求） | 实测 | 结论 |
|---|---|---|
| 存在文献锚点 | PhaZ1 **C183**（三联体 C183-D355-H388），定点突变证据 | ✅ |
| 模式在 723 条参考上**复现锚点** | PhaZ1 上模式命中位置 = **[183]**，三肽 = `VCQ`（唯一命中；该蛋白另 3 个 Cys 87/370/386 均不命中） | ✅ |
| miss 可判定 | 该角色**唯一**带坐标的锚点被完整复现（1/1）；按项目既有规则（"复现每一个带坐标的锚点"）miss 可判定 | ✅ |
| 模式有文献依据（非纯 data-derived） | **-1 位疏水**来自 Knoll 2009 的 `cysteine-1` 规则；**+1 位 Gln** 来自 PhaZ1 锚点上下文 | ✅ |
| 特异性 | 参考面板：Cys 271/276 = **0.9819** vs 非 Cys 22/447 = **0.0492**（精确 `VCQ` 层：0.9094 vs **0.0022**） | ✅ 可用 |

**为什么不是 (b)（家族 HMM 强绑定）**：Cys 候选的超家族归属本身就来自发现层家族 HMM，用同一来源再标一次"家族级 Cys 证据"是循环论证，不增加独立信息。

**为什么不是 (c)（显式降级）**：文献锚点真实存在且可复现，因此不必降到纯文档声明；(c) 的核心保护（不得写成"已验证"）已**同时**实现（见 §4）。

> **与 09-17 旧结论的关系**：09-17 把 `VCQ` 标为 `exploratory_in_sample_descriptor_discovery` / `load_bearing=false` 并禁止其作为过滤器，理由正是"**无文献依据**"。本 run 修正了该前提——`VCQ` 是 PhaZ1 C183 锚点的上下文，**有文献依据**。但 09-17 的另一条限定仍然成立并保留：**特异性只在样本内测过、无 held-out 集**，因此本列**仍未接入过滤器**。

---

## 2. 校准（模式 + 面板分离）

- 模式：`[疏水]-C-Q`（`疏水 = ACFILMVWY`）。
  - `-1` 疏水：Knoll 2009 `cysteine-1` 规则；
  - `+1` Gln：PhaZ1 锚点 `V182-C183-Q184`。
- 锚点复现（fail-closed）：模式必须在锚点序列的**报告位置**上命中，否则进程以 `SystemExit` 退出、**不写任何候选行**。

| 面板 | n | 模式命中 | 精确 `VCQ` |
|---|---:|---:|---:|
| Cys 超家族参考 | 276 | **271 (0.9819)** | 251 (0.9094) |
| 其余 7 超家族参考 | 447 | 22 (0.0492) | **1 (0.0022)** |

> `[疏水]-C`（只有 Knoll 2009 规则、不含 +1）在非 Cys 参考上是 **330/447 = 0.7383**——几乎无特异性。**特异性主要由 +1 的 Gln 提供**，这也是把锚点上下文（而非仅规则）写进模式的原因。

---

## 3. 候选层分层结果（139,719 条 = 池内 109,087 + 池外 30,632）

| 分组 | n | `anchor_exact_vcq` | `anchor_context_hydrophobic_c_q` | 模式命中率 |
|---|---:|---:|---:|---:|
| **Cys 超家族（全部）** | 32,904 | 26,506 | 3,298 | **0.9058** |
| ├ Cys **高可信度** | 28,534 | 24,885 | 2,832 | **0.9714** |
| └ Cys 被 hold | 4,370 | 1,621 | 466 | 0.4776 |
| extracellular dPHASCL type 1 | 87,485 | 543 | 4,625 | 0.0591 |
| intracellular nPHAMCL | 7,567 | 80 | 675 | 0.0998 |
| extracellular dPHASCL type 2 | 10,297 | 84 | 511 | 0.0578 |
| extracellular dPHAMCL | 534 | 1 | 14 | 0.0281 |
| PhaZ7-like / with lipase box / periplasmic | 316 | 11 | 22 | 0.104 |

**两点值得注意**：
1. **Cys 高可信度组的模式命中率（97.1%）显著高于同超家族被 hold 组（47.8%）**——这是与高可信度过滤相互独立的内部一致性证据。
2. 非 Cys 超家族的模式命中率仅 **2.8%–10%**，与 Cys 超家族的 90.6% 分离清晰。

---

## 4. 诚实边界（写进产物、测试钉死）

| 声明 | 实现方式 |
|---|---|
| **不得读作"已验证 Cys"** | 每行 `catalytic_activity_verified = false`；校准 JSON 里 `phenotype_boundary = "candidate-only; this column must not be read as a verified Cys"` |
| 不是 family 判定 | `family_call_made = false`；不产出任何 family/superfamily 列 |
| 不训练模型 | `hmm_fitted = false` |
| 未接入过滤器 | 本列**未**被 `filter_phaded_high_confidence.py` 消费；`filter_wired=false` |
| 独立验证未完成 | 每行 `independent_validation = "pending"`；校准 JSON `scope = in_sample_reference_panel`、`held_out_set = none` |
| 单一锚点 | `evidence_basis = sequence_context_match_to_single_literature_anchor_phaZ1_C183_PMID16233560` |

> **9-17 的 F7 条件仍然待办**：在独立非 Cys α/β-水解酶阴性集上评估；命中率 > 5% 则该描述符必须降级为纯描述列。本 run 只把 `VCQ` 从"无文献依据"升级为"有文献锚点"，**没有**解除 F7。

---

## 5. Test-first 证据

| 步骤 | 证据 |
|---|---|
| 先写 15 个测试（实现文件不存在） | `logs/step1_failing_test_red.log`：`Ran 15 tests … FAILED (errors=15)`，全部 `FileNotFoundError: … annotate_phaded_cys_nucleophile_motif.py` |
| 实现后转绿 | `logs/step3_focused_tests_green.log`：`Ran 18 tests … OK` |
| 实测暴露 2 个 bug → 先写测试钉住 → 再修 | 见 §6 |

**实测发现并修复的 2 个 bug（均由"真实数据跑一遍"暴露，而非单测）：**

1. **CLI 解析 bug**：`--candidate-fasta` 误加 `type=Path`，导致 `POOL=PATH` 解析抛 `TypeError: argument of type 'WindowsPath' is not iterable`。→ 先加 `CliTests.test_cli_accepts_pool_equals_path`（红灯：同 `TypeError`），再去掉 `type=Path`（绿灯）。
2. **accession 截断 bug（更严重）**：`accession_of()` 对候选 header `GCA_1|CTX_1_5` 取 `|` 之后的部分，产出 `CTX_1_5`，**丢失基因组前缀** → 下游任何 join 全部落空（首次分析跑出全部 `(unassigned)`）。→ 先加 `test_candidate_accession_keeps_full_header`（红灯：`'CTX_1_5' != 'GCA_1|CTX_1_5'`）与 `test_reference_lookup_still_resolves_piped_accession`，再改成"候选 header 原样输出、参考表用 `|` 之后部分解析"（绿灯）。

**验收命令（本 run 实测）**：

| 命令 | 结果 | 日志 |
|---|---|---|
| `python -m unittest pipeline.tests.test_annotate_phaded_cys_nucleophile_motif` | **Ran 18 tests — OK**，exit 0 | `logs/step3_focused_tests_green.log` |
| `python -m unittest discover -s pipeline/tests` | **Ran 642 tests — OK (skipped=1)**，exit 0 | `logs/full_test_suite.log` |
| `python -m compileall -q pipeline` | exit 0 | `logs/compileall.log` |
| `git diff --check` | exit 0 | `logs/git_diff_check.log` |

（642 = 本任务前的 624 + 本任务新增 18。）

---

## 6. 产物

| 产物 | 大小 (B) |
|---|---:|
| `results/cys_nucleophile_motif_evidence.tsv` | 24,160,768（139,719 行） |
| `results/cys_nucleophile_motif_calibration.json` | 2,311 |
| `results/cys_motif_tier_by_superfamily.json` | 1,990 |
| `input_contract.json` | 6,224 |
| `pipeline/scripts/annotate_phaded_cys_nucleophile_motif.py` | 新文件（实现） |
| `pipeline/tests/test_annotate_phaded_cys_nucleophile_motif.py` | 新文件（18 个测试） |

`input_contract.json` 记录全部输入/输出/源码的 path + size + SHA-256，GTDB 三槽 `pending`（未伪造哈希）。

---

## 7. 未做 / 未解决（如实登记）

1. **未接入过滤器**：本列是新增证据，不参与高可信度判定；是否可作为判据需先过 F7（独立阴性集）。
2. **无 held-out 集**：特异性数字全部是样本内。
3. **单一文献锚点**：Cys 型只有 PhaZ1 一条带坐标的实验证据（PhaZ2 有实验阳性但无催化残基坐标；PhaZd 是 Ser 型）。
4. **`C183S` 不失活这一事实本身**提示 Ser 可替代 Cys——因此"未命中该模式"**不得**读作"非 Cys 型"的强证据（参考面板上仍有 5/276 漏检）。
5. **冻结台账未改**：`runs/20260909_phaded_reference_evidence_01` 里 276 条 Cys 参考的 `catalytic_residues` 仍为 `pending`；本 run 以**独立 run** 提供锚点，未回写台账。
6. 未修改共享交接文档 `docs/T141_20260917_project_handoff.md`（按任务约束）。

---

*本状态文档为 candidate-only 记录；序列上下文匹配只表示候选同源或功能潜力，不等同于已验证的催化活性或 PHB/PHA 降解表型。*
