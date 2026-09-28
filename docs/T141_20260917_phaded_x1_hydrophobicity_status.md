# T141 x₁（lipase box）疏水分层状态

日期：2026-09-17
Run：`20260917_phaded_x1_hydrophobicity_01`
任务：计划 v2 Task 3（P1，纯分析，零新计算）
状态：`completed_candidate_only`

## 结论摘要

对候选宇宙 `109,087` 条蛋白的 `lipase_box_x1`（lipase box `GxSxG` 的 x₁ 位残基）做疏水性分层：

| 分层 | 条数 | 占全体 | 占已测（`77,302`） |
|---|---|---|---|
| `hydrophobic`（A,C,F,I,L,M,V,W,Y） | `76,739` | `70.35%` | **`99.27%`** |
| `non_hydrophobic`（D,E,H,K,N,Q,R,S,T,G,P） | **`563`** | `0.52%` | **`0.73%`** |
| `not_tested`（`lipase_box_x1` 为空） | `31,785` | `29.14%` | — |
| 合计 | **`109,087`** | `100%` | — |

- 自洽性：`76,739 + 563 + 31,785 = 109,087` ✅（manifest `totals_self_consistent=true`）；
- 与计划第 1.3 节记录的既有实测分布**逐项一致**（含各残基明细），本 Task 未改动任何数字；
- 非疏水 x₁ 的 563 条构成**脂酶/酯酶混淆嫌疑集**，此前项目从未识别过这批候选；
- x₁ 轴发出的**生物学阴性判定为 `0` 条**（`negative_call_count=0`）。

## 判据与口径

依据 Knoll M, et al. *The PHA Depolymerase Engineering Database*. BMC Bioinformatics 2009;10:89
（PhaDED；PMC2666664）：脂酶与酯酶的 `Gx₁Sx₂G` 基序中 x₁ 通常是**极性**残基，而 PHA 解聚酶几乎全部为**疏水**残基。
因此本 Task 把 x₁ 落成三分层：

- `hydrophobic` = `A,C,F,I,L,M,V,W,Y`；
- `non_hydrophobic` = `D,E,H,K,N,Q,R,S,T`（极性/带电）+ 单列标注的 `G`（glycine，无侧链）与 `P`（proline，环状仲胺）。
  G 与 P 计入非疏水（因此不能满足 PHA 解聚酶的疏水 x₁ 特征），但不与极性残基合并，逐个残基单独标注；
- `not_tested` = x₁ 为空值 = **未测**，不是阴性。

`C`（cysteine，8 条）按项目既定口径计入疏水集。判据表已随 run 版本化落盘：
`inputs/x1_residue_policy.tsv`（20 种标准氨基酸 + `not_tested` 行，含引用）。

实现为 **fail-closed**：出现非标准残基（如 `X`/`B`/`Z`/`J`/`U`/`O`/`*`/小写字母/多字符）时直接报错终止，
不把未知值悄悄归入任何一类。本次运行输入仅含 20 种标准大写残基与空值，未触发该分支。

## 可发表口径

在 `109,087` 条候选中，有 `77,302` 条（`70.86%`）测得了 lipase box 的 x₁ 残基；其中
**`76,739` 条（`99.27%`）的 x₁ 为疏水残基，符合 PHA 解聚酶的 x₁ 特征**，而
**`563` 条（`0.73%`）为疏水 x₁ 之外残基，构成脂酶/酯酶混淆嫌疑集**——该嫌疑集在本次分析之前从未被识别。

`99.27%` 的分母是**已测**候选（`77,302`），不是全体候选；`31,785` 条未测候选既不支持也不削弱该结论。

## 混淆嫌疑集（563 条）的构成

`results/confounder_candidates.tsv`（563 条 + 表头）：

| 维度 | 分布 |
|---|---|
| x₁ 分组 | `basic` 178（R 84、K 82、H 12）· `acidic` 151（D 89、E 62）· `polar_uncharged` 124（T 47、Q 36、S 31、N 10）· `glycine` 80 · `proline` 30 |
| 残基 | D 89 · R 84 · K 82 · G 80 · E 62 · T 47 · Q 36 · S 31 · P 30 · H 12 · N 10 |
| `motif_panel_status` | `partial_catalytic_pattern_panel` 563（全部） |
| 基因组 | 分布于 `560` 个基因组；单基因组最多 2 条 |
| `recommended_confidence_cap` | `moderate_candidate_only` 563（全部） |

（`basic` 分组按惯例含 His；分组列仅为描述性标注，实际分层只由疏水集是否命中决定。）

## 置信度上限语义

`results/x1_stratification.tsv` 的 `recommended_confidence_cap` 是**本单轴**给出的上限：

| x1_class | `recommended_confidence_cap` | 含义 |
|---|---|---|
| `hydrophobic` | `high_candidate_only` | x₁ 轴不额外下调（是否达到 high 仍取决于其他判据） |
| `non_hydrophobic` | `moderate_candidate_only` | x₁ 轴把上限压到 high 以下 |
| `not_tested` | `not_assessed_by_x1_axis` | x₁ 轴不提供信息，既不通过也不否决 |

候选的**实际**上限是各判别轴上限的最小值；本 Task 不发布跨轴合成结论（合成属 Task 5 及以后）。
硬性不变量：**非疏水 x₁ 行的 `recommended_confidence_cap` 永不为 `high_candidate_only`**，已由测试与运行后核验确认（563/563）。

## 边界声明（必须随结论引用）

- **非疏水 x₁ 不是生物学阴性。** 它只是提示该候选的归属需要更强的判别证据（催化域 type 1/2 几何、
  lid、AHSMG、SBD/linker、结构、系统发育等），既不证明它是脂酶/酯酶，也不证明它不降解 PHA。
- **未测（`not_tested`）也绝不是阴性**：`31,785` 条只是没有测到 lipase box 的 x₁，属于证据缺失。
- x₁ 轴在本 Task 中**删除 0 条候选、重判 0 条候选**：563 条全部保留在候选宇宙内，只被标注为混淆嫌疑。
- 本项目全部 HMM、profile、domain、motif、SignalP、结构与系统发育证据都只表示**候选同源或功能潜力**，
  **不等同于已验证的 PHB/PHA 降解表型**。

## 产出与可复现性

Run 目录：`runs/20260917_phaded_x1_hydrophobicity_01/`

| 文件 | 内容 |
|---|---|
| `results/x1_stratification.tsv` | 109,087 行，按 `accession` 升序；列：`accession, lipase_box_x1, x1_class, confounder_flag, recommended_confidence_cap, x1_group` |
| `results/confounder_candidates.tsv` | 563 行嫌疑集；列：`accession, genome, lipase_box_x1, x1_group, confounder_flag, recommended_confidence_cap, motif_panel_status` |
| `results/x1_summary.tsv` | 9 个指标行（含 `negative_call_count 0`） |
| `results/x1_residue_breakdown.tsv` | 每残基计数与占比（20 残基 + `not_tested`） |
| `results/x1_manifest.json` | 判据、计数、输入哈希、输出哈希、边界声明 |
| `inputs/x1_residue_policy.tsv` | 版本化判据表（含 PMC2666664 引用） |
| `logs/x1_stratification.log` | 执行日志（含 `server_execution_started: false`） |
| `input_contract.json` | 由 `pipeline/scripts/run_context.py` 生成 |
| `run_manifest.json` | 计数、授权标志、源码快照哈希、非改动声明 |

输入（只读消费，未修改）：

- `runs/20260915_phaded_motif_reconciliation_01/results/motif_candidate_evidence.tsv`
  size `64,498,737`，SHA-256 `796cb0a32cc7fde1c768b2c0f2f8e30542c7aa37e773f2b8004a482ac878569a`
  —— 与既有历史 contract 中记录的同文件哈希一致，可确认输入未被改写。

确定性核验：在临时目录重跑同一输入，四个结果 TSV 的 SHA-256 与 run 内产物**逐字节一致**；
`x1_stratification.tsv` 亦已核验按 `accession` 严格升序且无重复 accession。

测试：`pipeline/tests/test_stratify_phaded_lipase_box_x1.py`，10 个测试，覆盖 20 种氨基酸口径、
G/P 单列标注、空值 → `not_tested` 非阴性、非疏水上限不为 `high_candidate_only`、
端到端输出的排序/自洽/嫌疑集过滤、未知残基与重复 accession 的 fail-closed、以及输入文件字节不变。

## 已知口径说明（不是缺陷）

1. `input_contract.json` 的 `status` 为 `pending`：合同模板按项目规则固定声明 GTDB
   `taxonomy`/`metadata`/`tree` 三项，而本 Task 为候选 motif 表的纯分析，**未消费任何 GTDB 参考文件**，
   故这三项 `path=null`。本 Task 实际消费的两个输入均为 `verified`（path/size/SHA-256 均已记录）。
2. `input_contract.json` 的 `generated_at` 为 UTC（`2026-09-16T17:01:22+00:00`），与本机本地日期
   `2026-09-17` 相差时区，属 `run_context.py` 的既定行为，非日期不一致。
3. 相对计划要求的最小列集，`x1_stratification.tsv` 增加了末列 `x1_group`，用于把 G/P 与极性残基
   区分开来（计划 Step 3.1 要求 G 与 P「单列标注」）；另额外产出 `x1_residue_breakdown.tsv`。

## 非改动与非授权声明

- 只新建 `runs/20260917_phaded_x1_hydrophobicity_01/` 与 `pipeline/scripts/stratify_phaded_lipase_box_x1.py`
  及其测试、本状态文档；**未覆盖、未修改、未删除**任何历史 `runs/`、`results/`、`pipeline/config/`、
  `AGENTS.md` 或历史 HMM；未修改 `pipeline/config/formal_scan_models.tsv`。
- 全程本地执行：**未启用 SSH、未连接服务器、未提交任何服务器任务**；
  `server_execution_started=false`、`formal_scan_started=false`、`formal_registry_modified=false`。
- 未执行任何 `git commit` / `push` / `clean` / `reset`，未删除任何文件。
- 本 Task 不降低校准 gate，不合并/删除任何家族，不改变任何 profile registry 状态。

*本文件为 candidate-only 研究产出；所有结论表示候选同源或功能潜力，不等同于已验证 PHB/PHA 降解表型。*
