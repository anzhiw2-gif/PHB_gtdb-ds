# T141 PhaDED 定向证据运行授权包（Task 11 准备件）

> **计划**：`docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign.md`（Task 11，Step 1–5）
> **设计**：`docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md`（四层模型体系、验收边界）
> **状态文档**：`docs/T141_20260928_phaded_evidence_model_redesign_status.md`（v1 权威表、授权边界表）
> **接手文档（v1 冻结口径）**：`docs/T141_20260921_phaded_comprehensive_review.md`
> **本文档性质**：**准备件（preparation artifact）**。所有计数、路径、工具与先例均逐条标注来源；来源未给出的数字一律写 `pending`，不推断、不编造哈希或时长。

---

## 0. 性质与授权边界声明

本文档**只是一份准备件，用于支持操作者的授权决策**：它列出五个定向证据运行包（按科学优先级排序）的决策问题、目标群体与精确 accession 列表推导、抽样或全集依据、工具与版本、假阳性竞争者、输出 schema、证据标记变更范围、线程与资源规则、预估规模与时长，以及**明确不属于本次授权的动作**。**截至本文档编写时，这五个运行一个都没有执行**：没有创建任何 `runs/<run_id>/` 或 `deploy/<run_id>/`，没有调用 SignalP、HMMER、MAFFT、Foldseek、ESMFold、AlphaFold 或任何计算工具，没有安装或变更任何环境，没有触碰 `runs/`、`results/`、`deploy/` 下的任何已有文件（只读访问除外），历史结果与冻结权威表保持原样。**每个计算类别需要分别、显式的授权**：按计划 Task 11 Step 3，**SignalP、HMM building、Foldseek/结构预测、GTDB searching 是四项彼此独立的批准，批准其中任何一项都不授权其余任何一项**，也不授权创建运行目录、安装软件、修改服务器配置、删除或归档数据、push 到 GitHub。本文档本身也不构成对任何运行的批准。

---

## 1. 总览：五个授权包（按科学优先级排序）

| 优先级 | 运行名（计划给定 / 建议） | 计算类别 | 目标条数 | 需要授权 | 状态 |
|---:|---|---|---:|---|---|
| P1 | `20260928_phaded_reference_signalp_v2_01`（计划给定） | SignalP 6.0h 序列层预测 | 723 | SignalP 执行 + 创建 dated run/deploy | `pending-authorization` |
| P2 | `20260928_phaded_type1_accessory_domain_v2_01`（建议名，计划未给定，见 §3.2.4） | HMMER 3.4 / Pfam-A 域扫描（`hmmscan`） | 11,610 | HMMER 执行 + 创建 dated run/deploy | `pending-authorization` |
| P3 | `20260928_phaded_cys_mapping_v2_01`（计划给定） | 参考锚定残基映射（项目脚本；比对重建属另一类） | 29,974 | 创建 dated run/deploy；如需重建 MAFFT 比对则另需「建比对/HMM」授权 | `pending-authorization` |
| P4 | `20260928_phaded_structure_review_v2_01`（计划给定） | Foldseek 结构比对 +（可选）结构预测，**竞争面板** | 987 | Foldseek/结构预测执行 + 创建 dated run/deploy | `pending-authorization` |
| P5 | `20260928_phaded_with_lipase_deferred_sample_v2_01`（建议名，计划未给定，见 §3.5.4） | 分层抽样（本地只读）+ 抽样子集的结构召回（Foldseek/结构预测） | 抽样规模 `pending`（池 1,206,655，**禁删、禁入任何主候选计数**） | 创建 dated run/deploy；结构召回另需授权 | `pending-authorization` |

**优先级依据**：顺序即计划 Task 11 Step 1 的 1–5 项顺序，也即接手文档 §5.3 与 §9 给出的「建议的裁决动作」顺序——先标定信号肽门严格度（P1，最大证据标准缺口），再据此裁决被信号肽门单独卡住的最大候选集（P2），再处理亲核体层最大的推断风险（P3），再处理序列层已定论不可判别的家族（P4），最后处理体量最大但已被隔离的噪声层（P5）。

---

## 2. 五个包共用的通用条款

### 2.1 命名与 run 契约（授权后才创建）

每个真实运行必须同时具备：

```text
deploy/<run_id>/                                  # 服务器只允许执行 dated deploy 中已绑定源码
runs/<run_id>/logs/
runs/<run_id>/inputs/
runs/<run_id>/results/
runs/<run_id>/input_contract.json
```

- 依据：`AGENTS.md` L2–L6；计划 Global Constraints L13–L16。
- **这些目录一律由协调者/操作者在授权之后创建**，本文件不创建、也不预置任何内容。
- `run_id` 必须可审计且**不含路径遍历**：本文件建议 `^[0-9]{8}_[a-z0-9_]+$`，含 `..`、`/`、`\`、盘符、绝对路径或空白字符者一律拒绝；创建前必须打印并人工核对最终绝对路径，确认落在 `D:\PHB_gtdb-ds\runs\`（本地）或服务器既有 `runs/` 根之下。
- 计划给定的三个 run 名（P1/P3/P4）沿用计划 Task 11 的文件列表原文；**日期前缀必须替换为实际启动日**（`AGENTS.md` 要求 dated deploy），因此上面的 `20260928` 是计划提议值，不是固定 ID。
- `input_contract.json` 必须绑定：源码快照、参考 FASTA/台账、HMM/比对、环境与工具版本、命令、线程数、文件 size 与 SHA-256；**缺失项写 `pending`，不得伪造哈希**（`AGENTS.md` L5、L12）。
- 五个包都**不得覆盖**历史 `results/`、历史 `runs/` 或服务器历史目录（`AGENTS.md` L2）。

### 2.2 线程与资源（动态，禁止固定值）

- 服务器计算单任务线程上限 **40**；启动前必须按实际空闲动态取 `min(40, 总核数 − 负载 − 10)`，保留约 10 核余量；总核数以 `nproc`、负载以 `/proc/loadavg` 1 分钟值**实测**为准。
- 依据：`AGENTS.md` L7；计划 Global Constraints L16（结果 <1 时必须中止）。
- **该值在启动时刻实测，永不预先固定**。本文档不给出任何具体线程数作为「本次将使用」的值。
- 历史先例中的超限值不可沿用：`docs/T141_20260920_phaded_three_gaps_status.md` L5、L17 记录 2026-09-20 的 SignalP 补测使用 60 线程（`--torch_num_threads 60`），那是操作者当次临时授权；`docs/T141_20260921_phaded_comprehensive_review.md` §8.1 记录 gRodon 侧 70 线程同样为临时授权。**P1 及所有包都不得复制这些数值**。
- 代码侧已固定：`pipeline/scripts/run_phaded_full_signalp.py` L12 `MAX_THREADS = 40`，L24–L25 在 `torch_threads + write_processes > 40` 时直接抛错。

### 2.3 数字溯源与 `pending` 规则

- 本文档引用的每个计数都标注了文件路径与章节；**凡计划/接手文档未给出数字处一律写 `pending`**。
- 标注为「本会话只读复核」的计数，是在**没有创建任何 run、没有调用任何计算工具**的前提下，用标准库 `csv` 对既有冻结 TSV 做的只读重算；运行授权后仍须在 `input_contract.json` 中重新计算并绑定 SHA-256，不得直接抄用本文档数字。

### 2.4 结果一律只改证据标记，不删候选

- 五个包的输出**都只能新增证据列/证据标记**，**任何候选在任何情况下都不被删除**：不删行、不删文件、不改写历史结果。
- 允许变化的字段仅限：`evidence_flags`（多值、非互斥）、以及由证据驱动的 `primary_disposition` 迁移；受控词表来自设计文档 §候选目录 与计划 Task 3 Step 2：
  `core_sequence_homolog` / `probable_sequence_homolog` / `remote_homolog_candidate` / `function_unresolved` / `deferred_structure_review` / `excluded_input_quality`。
- **不可变化**：任何包的阴性/未命中结果都**不得**产生 `excluded_input_quality`（设计文档 §候选目录 第 6 条：不能把科学不确定性写成输入错误）；`excluded_input_quality` 仅限序列缺失、明显截断或输入错误。
- **`deferred_structure_review` 的 with-lipase 层（P5 对象）在任何包中都不得进入主候选/高可信度计数**（`AGENTS.md` L11）。

### 2.5 通用的未授权动作（适用于全部五个包）

批准任何一个包，都**不**授权：创建任何未在批准文本中点名的 run/deploy；安装或升级任何软件/环境；修改服务器配置；删除、移动或归档任何数据；重建 MAFFT 比对或 HMM（计划 Authorization Checkpoint 3）；对 GTDB 分片做 HMMER 重扫（Checkpoint 4）；用真实项目数据重算 v2 候选目录（Checkpoint 5）；push 到 GitHub（Checkpoint 6）。

---

## 3. 五个授权包

## 3.1 P1 —— 723 参考层 SignalP 补测（标定信号肽门严格度）

### 3.1.1 决策问题

「当前候选层 type-1 门（必须预测到分泌信号，否则 hold）所依据的判据，在**已知胞外**的参考层上检出率是多少？」——即用一个参考层实测的检出率，为信号肽门设定可辩护的严格度，而不是继续让参考层不受该标准检验。

### 3.1.2 目标群体

- **规模：723 条**，即冻结参考台账全量，不做抽样（全集）。
- 权威来源表：`runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv`
  - 该表为 v1 冻结权威件（`docs/T141_20260928_phaded_evidence_model_redesign_status.md` §1.1 行「校准台账（723 参考）」）。
  - 实测口径（本会话只读复核）：**723 数据行 × 26 列**；表头列为 `reference_id / accession / sequence_sha256 / sequence_source / source_database / source_database_version / retrieval_date / organism / taxonomy_id / phaded_superfamily / phaded_family_id / reported_localization / substrate_class / substrate_detail / experimental_assay / experimental_result / evidence_status / positive_negative_control / catalytic_residues / lipase_box_or_ahsmg / oxyanion_hole_evidence / domain_architecture / primary_doi / pmid / pmcid / notes`。
- 组成：**30 `experimental_positive` + 693 `annotation_only`**（`docs/T141_20260921_phaded_comprehensive_review.md` §2.1）。
- 序列输入：`runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference.faa`（size 300,777 B，本会话只读实测）。
- **accession 列表推导（无需筛选，全集直取）**：`accession` 列 723 个值，去重后应仍为 723；重复或缺失即 fail-closed 中止，不补、不猜。
- 存在的证据缺口（本包要解决的正是它）：整个 723 台账**0 行提及 signal**，参考层**从未跑过 SignalP**；284 条 type-1 参考的 `domain_architecture` 全为 `pending`（接手文档 §5.2 缺口 #4）。
- 分层键（计划 Task 11 Step 1 第 1 项要求的三层）：
  1. **实验定位**：`reported_localization`（实验/来源记录的定位）；
  2. **功能证据**：`evidence_status` + `experimental_assay` + `experimental_result` + `primary_doi/pmid`（区分有功能证据的 30 条与纯注释的 693 条）；
  3. **继承历史标签**：`phaded_superfamily` + `phaded_family_id`（继承自 DED 超家族名，属**定义性**标签而非实测，接手文档 §5.2 缺口 #4 原文：「"extracellular" 标签是**定义性**的」）。
  - 注：计划 Task 4 新增的 15 个证据质量与独立性列（含 `localization_evidence`、`experimental_evidence_grade`、`independence_group`）**尚未落到真实台账**（状态文档 §1.3 末行：「真实台账 curation 未授权」）。因此本包的分层键**只能使用现有 26 列**；若操作者希望按 E3/E2/E1/A 分级分层，必须先把 Task 4 Step 6 的台账 curation 一并授权，否则该维度写 `pending`。

### 3.1.3 抽样或全集依据

- **全集**（723/723），不抽样：目标是**参考层检出率**这一比例量，723 条本身已是全量且规模小；抽样会引入与候选层不可比的抽样误差。
- 输出必须**按 §3.1.2 的三个分层键分组报告检出率**（每组给出分子/分母），并单独给出 284 条 type-1 参考、276 条 Cys 型参考、52 条 nPHAMCL 参考的分组结果——这些是后续门严格度讨论的分母。
- 若操作者要求分层抽样的替代方案，抽样设计、分层键权重与随机种子必须**预先书面固定**，种子值写 `pending`（本文档不预设）。

### 3.1.4 工具与版本

- 工具：**SignalP 6.0h**，`--organism other --mode fast`；执行入口为项目脚本 `pipeline/scripts/run_phaded_full_signalp.py`（L28–L44 构造命令），授权后按需修改（计划 Task 11 文件列表已列该脚本为「授权后修改」）。
- 版本记录位置（可核查的既有出处）：
  - `docs/T141_20260920_phaded_three_gaps_status.md` L17：「SignalP 6.0h（`--organism other --mode fast --torch_num_threads 60 --write_procs 8`，~31 min，returncode 0）」；
  - `docs/T141_20260916_phaded_full_library_signalp_01_status.md` L8：「SignalP 6.0 fast 模式已在 `<SERVER_HOST>` 完成，使用 dated deploy」；
  - 环境归属：conda env `signalp6`（`docs/项目详细审查报告.md` L60：「`signalp6`（SignalP 6.0）」）。
- 本次运行的**精确二进制路径与 `signalp6 --version` 输出**：`pending`（必须在授权后的 `input_contract.json` 中记录 path + size + SHA-256；不得复用历史 run 的哈希）。
- GPU：不使用（`docs/T141_20260916_phaded_full_library_signalp_01_status.md` L13）。

### 3.1.5 假阳性竞争者

本包要对抗的**不是序列假阳性，而是判据不对称**（接手文档 §5.3，审查重点）：

- **候选层 type-1 门**：必须预测到分泌信号，否则 hold；实测 **48,038 / 67,720（71%）** 因此被 hold（接手文档 §5.3 第 1 条；权威机器可读来源 `runs/20260919_phaded_nucleophile_unify_01/results/filter_summary.json` → `per_superfamily["extracellular dPHASCL type 1"].fail_localization = 48038`、`total = 67720`）。
- **参考层 type-1**：**从未被这条标准检验**（接手文档 §5.3 第 2 条）。
- **文献反例（真实假阴性）**：Knoll 2009 自录 `gi:74267419`（*R. eutropha* H16）被归入胞外 type-1，且原文写明其 "lacking a signal peptide, a linker domain, and a substrate binding domain"；`gi:194292521`（*C. taiwanensis*）GenBank 注胞内、按序列相似性判胞外（接手文档 §5.3 第 3 条）。
- 因此本包的判别对象是：**「已知胞外参考中预测不到信号肽的比例」**。若该比例不可忽略，则「无分泌信号 ⇒ hold」这一门本身需要放宽或改写；若该比例接近 0，则 11,610 条（P2）更可能属折叠误判。**两种结局都只改证据口径，不改候选存在性。**

### 3.1.6 输出 schema

- 原始输出：SignalP 6.0h `--format txt` 的 `prediction_results.txt`（沿用既有布局，见 `runs/20260920_phaded_three_gaps_01/results/signalp_29521/prediction_results.txt`）。
- 逐参考注释表（建议文件名 `reference_signalp_annotation.tsv`），字段与受控词表：

| 字段 | 类型 | 受控词表 / 来源 |
|---|---|---|
| `accession` | string | 台账 `accession` 原样 |
| `reference_id` | string | 台账 |
| `phaded_superfamily` | enum | 台账既有 8 值（不加新值） |
| `phaded_family_id` | string | 台账 |
| `evidence_status` | enum | `experimental_positive` / `annotation_only` |
| `reported_localization` | string | 台账原样（**实验/来源定位，不是预测**） |
| `signalp_class` | enum | `SP` / `LIPO` / `TAT` / `TATLIPO` / `OTHER` / `PILIN` / `tool_input_excluded` / `pending` |
| `transport_signal_prediction` | enum | 同 `signalp_class`（设计文档要求与定位拆列） |
| `localization_evidence` | enum | `experimental` / `other_evidence` / `historical_family_label_only` / `unknown` |
| `permanence_note` | string | 固定串：`prediction_only_not_experimental_localization` |
| `tool_name` / `tool_version` | string | `signalp6` / `6.0h`（精确版本以运行内 `--version` 输出为准） |
| `run_id` / `source_sha256` | string | 运行内绑定 |

- **禁止**：把 `signalp_class` 直接写入候选层的 `localization_evidence=experimental`；设计文档与计划 Task 3 Step 2 的测试 `test_signalp_does_not_become_localization_truth` 已把这条钉死。

### 3.1.7 结果如何只改变证据标记/处置

- **不删除任何参考、任何候选**；723 参考台账**只读**，输出写在新 run 的 `results/`，台账本体不动（与 `docs/T141_20260920_phaded_three_gaps_status.md` L59「未改冻结台账（amendment 追加制）」同一原则）。
- 允许变化：参考层新增 `transport_signal_prediction` 证据列；候选层据参考层检出率**重新标定**信号肽门的严格度标记（写入 `evidence_flags`，例如 `localization_gate_calibrated_against_reference_detection_rate`）。
- 不允许变化：`primary_disposition` **不因本包单独迁移**——本包产出的是**门的标定值**，不是对任何具体候选的判定；把结果直接用于迁移某条候选的处置，必须走后续独立授权。
- 特别地：`excluded_input_quality` 在 P1 中**不可能被触发**。

### 3.1.8 线程与资源

- 动态规则同 §2.2：`min(40, nproc − load_1m − 10)`，实测于启动时刻；由 `torch_num_threads + write_processes ≤ 40` 约束（`pipeline/scripts/run_phaded_full_signalp.py` L24–L25）。
- 具体数值写 `pending`（启动前实测记录，随 `input_contract.json` 落盘）。
- 本包为 723 条小规模输入，**不得**复制历史 60 线程先例。

### 3.1.9 预估规模与运行时长

| 项 | 值 | 来源 |
|---|---|---|
| 输入规模 | 723 条 | 台账（本会话只读复核） |
| 本项目实测先例 A | **108,736 条 → 约 2 小时 41 分钟**（2026-09-16，13:05–15:46 Asia/Shanghai） | `docs/T141_20260916_phaded_full_library_signalp_01_status.md` L17 |
| 本项目实测先例 B | **29,521 条 → 约 31 分钟**（SignalP 6.0h fast，60 线程临时授权） | `docs/T141_20260920_phaded_three_gaps_status.md` L17 |
| 本次预估 | **`pending`** | 723 条规模本项目**未实测**；两个先例的输入规模差异近 3 个数量级，线性外推属编造，本文档不做 |

### 3.1.10 明确不属于本次授权的动作

批准 P1 **不**授权：重建 MAFFT 比对或 HMM（Checkpoint 3）；Foldseek 或任何结构预测（Checkpoint 4）；对 GTDB 分片做 HMMER 重扫（Checkpoint 4）；gRodon（Checkpoint 4）；用真实数据重算 v2 候选目录（Checkpoint 5）；安装/升级 SignalP 或任何环境、修改服务器配置、删除/归档数据、push（Checkpoint 6）；创建除 P1 之外的任何 run/deploy；把 SignalP 结果写成实验定位或表型证据。

---

## 3.2 P2 —— 11,610 条 SP-less type-1 候选的配件域核查

### 3.2.1 决策问题

「在 11,610 条因**未预测到分泌信号**而被 hold、但催化域几何判定为 `type1_verified` 的候选里，有多少条具有**独立于信号肽的**胞外配件域证据（FN3 / Ig-like / TSP3 / CHB），从而支持它们是真 type-1（Knoll 模式：信号肽未检出/序列缺失）而非折叠误判？」

### 3.2.2 目标群体

- **规模：11,610 条**（计划 Task 11 Step 1 第 2 项；接手文档 §5.3 第 4 条与 §9 第 2 项）。
- 权威来源表（v1 冻结，只读）：`runs/20260919_phaded_nucleophile_unify_01/results/hold_candidates.tsv`（`docs/T141_20260928_phaded_evidence_model_redesign_status.md` §1.1 冻结表 + §1.2；该表 11 列：`accession / genome / superfamily / nucleophile_type / hold_reason / lipase_box_state / sbd_pf06850_binding_state / signalp_class / catalytic_domain_type / interpro_status / profile_evidence_status`）。
- **精确推导（本会话只读复核，未创建 run）**——单键过滤即精确复现 11,610：

```text
hold_candidates.tsv
  WHERE superfamily = 'extracellular dPHASCL type 1'
    AND hold_reason = 'localization_conflict'
    AND catalytic_domain_type = 'type1_verified'
  => 11,610 行 / 11,610 个唯一 accession / 8,942 个基因组
```

  - 逐级计数（同一只读复核）：`hold_candidates.tsv` 全表 **71,860** 行；type-1 子集 **62,836** 行（`localization_conflict` 48,038 / `architecture_criteria` 14,433 / `common_criteria` 365）；type-1 ∩ `localization_conflict` = **48,038** 行；再加 `catalytic_domain_type = 'type1_verified'` = **11,610** 行。
  - 该 48,038 与 `runs/20260919_phaded_nucleophile_unify_01/results/filter_summary.json` 的 `fail_localization = 48038` 一致；71,860 与同文件 `total_held = 71860` 一致。
- **对交接文档措辞的一处诚实更正（重要，避免不可复现的过滤条件）**：接手文档 §5.3 把该集合描述为「同时具备 `type1_verified` 几何 + trained profile 命中 + InterPro 支持」。只读复核显示：
  - 这 11,610 条中 `profile_evidence_status` = `profile_trained_hit` **11,602** 条、`profile_ambiguous_family` **8** 条 → 若把「trained profile 命中」写成强制性过滤键，会**漏掉 8 条**；
  - 48,038 条的全部 `interpro_status` 均为 `interpro_supported` → 「InterPro 支持」在该文件里**不是区分性条件**。
  - 因此本包采用的**唯一可复现过滤器**是 §3.2.2 的三键式（其中 `superfamily` 与 `hold_reason` 界定母集，`catalytic_domain_type` 界定目标集）；「trained profile 命中」与「InterPro 支持」只作为**输出列的记录项**，不作为过滤键。
- 这 11,610 条的实测画像（同一次只读复核）：`signalp_class` 全为 `OTHER`（11,610/11,610，即全部**未预测到分泌信号**）；`lipase_box_state` 全为 `supported`；`nucleophile_type` 全为 `ser`；`sbd_pf06850_binding_state` 全为 `not_detected_in_tested_pfam`（该状态本身是**用错标记**的历史产物，见下）。
- 已有部分覆盖的证据层（v1 冻结，可直接复用、不需重算）：`runs/20260920_phaded_three_gaps_01/results/type1_sbd/type1_accessory_domain_layer.tsv` —— 覆盖 **5,018** 条 type-1 高可信度候选（池内 **4,884** 条已标注 `in_matrix=yes`；池外 **134** 条 `pending`，fail-closed）。**注意：这 5,018 条是「已通过」的 type-1，与本包 11,610 条「被信号肽门 hold」的集合不重叠**（来源：`type1_sbd_annotation_summary.md` §4–§6；接手文档 §5.3 亦明示该层「已有 5,018 条」是既有资产）。本包的增量是**把同一套配件域标注施加到 11,610 条 hold 集合上**。

### 3.2.3 抽样或全集依据

- **全集**（11,610/11,610），不抽样：这是**裁决性**集合——每一条的归属都会影响门严格度的解释；且规模对 `hmmscan` 属小批量。
- 分层键（用于分组报告与后续敏感性分析，**不用于抽样**）：`catalytic_domain_type`（本集合内恒为 `type1_verified`，仅作核对）、`profile_evidence_status`（`profile_trained_hit` 11,602 / `profile_ambiguous_family` 8）、按 `genome` 聚合的**属/科级分类**（8,942 个基因组 → 需 join GTDB taxonomy，join 源与版本 `pending`）、以及 `lipase_box_x1` 疏水性（该列在 hold 表中不存在，需从 v1 矩阵或 motif 层 join，源 `pending`）。
- 若操作者选择「先跑分层代表子集」，则抽样规模、分层键权重与随机种子必须预先书面固定；**种子与规模写 `pending`**。

### 3.2.4 工具与版本

- 工具：**HMMER 3.4（Aug 2023）的 `hmmscan`**，对 **Pfam-A** 库做域扫描；标记集沿用项目已固定的 12 个信号（`runs/20260920_phaded_three_gaps_01/results/type1_sbd/type1_sbd_annotation_summary.md` §6）：
  - FN3：`PF00041` / `SM00060` / `cd00063` / `SSF49265` / `PS50853` / `G3DSA:2.60.40.10`
  - Ig-like：`PF16403` / `PF17957`
  - TSP3：`PF02412`
  - CHB：`PF13290`
  - 催化域（核对用）：`PF10503`（`IPR010126`，GO extracellular）
  - **明确不用 `PF06850` 作为胞外 SBD 标记**（`type1_sbd_annotation_summary.md` §2、§5；接手文档 §6 错误 #1）
- 版本记录位置：
  - HMMER 3.4 (Aug 2023)：`docs/T141_20260911_phaded_pfam_architecture_01_status.md` L11（`${PHB_REMOTE_ROOT}/miniconda3/envs/phb_gtdb/bin/hmmscan`，`--cut_ga --cpu 40`）；同版本另见 `docs/T141_20260906_archaea_phazh1_domain_annotation_03_status.md` L9；
  - Pfam-A 库：`${PHB_REMOTE_ROOT}/GTDB/pfam/Pfam-A.hmm`，size **2,246,909,846 B**，SHA-256 `4b0da6399b97d2b23329de82190b2e5705bde1dea4cb2d028a52712d5966739e`（`docs/T141_20260911_phaded_pfam_architecture_01_status.md` L10）；
  - **Pfam-A 的发布版本字符串（如 Pfam 3x.x）：`pending`** —— 项目文档只记录了路径、size 与 SHA-256，未记录 release 号，本包不推测。
- 可选替代：**InterProScan 5.78-109.0**（`docs/T141_20260914_phaded_interpro_full_01_status.md` L29；`docs/superpowers/plans/2026-09-14-phaded-evidence-completion.md` L9 记为 `InterProScan 5.78-109.0`）。若操作者选择 InterProScan 路径，那是**同一计算类别内的另一种工具**，仍需在本包批准文本中显式点名，不得默认替代。
- 运行脚本：计划 Task 11 未点名本包的脚本；v1 中同类产物由 `runs/20260920_phaded_three_gaps_01/plan/build_type1_accessory_layer.py` 生成（`type1_sbd_annotation_summary.md` §6 原文）。授权后应把逻辑移植进**新的 dated deploy**（`deploy/<run_id>/`），**不得直接运行服务器根目录或历史 run 下的旧脚本**（`AGENTS.md` L6）。
- 建议 run 名：`20260928_phaded_type1_accessory_domain_v2_01`（**建议值，计划未给定**；日期前缀须改为实际启动日）。

### 3.2.5 假阳性竞争者

- **信号肽门单独卡住的 71%**：候选层 type-1 门要求预测到分泌信号，实测 **48,038 / 67,720（71%）** 因此被 hold（接手文档 §5.3 第 1 条；`filter_summary.json` `fail_localization = 48038`、`total = 67720`）。本包正是在这 71% 内部做裁决，因此**信号肽阳性不能作为本包的判别依据**——判别必须来自**独立于信号肽**的配件域证据。
- **折叠误判（generic α/β-hydrolase / 酯酶型背景）**：`PF10503` 是 `Esterase, PHB depolymerase` 催化域（`type1_sbd_annotation_summary.md` §3），单独出现即为通用酯酶样折叠信号；因此「只有催化域、无任何配件域」是**默认结局而非阳性证据**。既有实测：5,018 条已通过 type-1 中 **≥1 配件域仅 1,430 条（29.3%）**，其余约 **71% 仅催化域 `PF10503`**（`type1_sbd_annotation_summary.md` §6）→ 本包必须显式与之对照。
- **短型/截断与基因模型伪影**：合作方既有开放问题原文——「"配件域趋异未入 Pfam" vs "短型无配件域" 需 Foldseek/FN3·Ig HMMER 才能彻底区分」（`type1_sbd_annotation_summary.md` §5 第 3 条）。因此**配件域缺失不得写成阴性**，只能写 `pending`/`not_detected_in_tested_models`。
- **同源冗余**：8,942 个基因组承载 11,610 条候选，比例接近 1.3:1，存在同一基因簇多拷贝；计数必须同时给出「条数」与「去冗余后条数」两个口径，去冗余键 = 序列 SHA-256（v1 冻结台账已有 `sequence_sha256` 字段范式）。

### 3.2.6 输出 schema

沿用 v1 既有证据层列（可直接对照 `runs/20260920_phaded_three_gaps_01/results/type1_sbd/type1_accessory_domain_layer.tsv` 的 11 列：`accession / pool / in_matrix / has_fn3 / has_ig_like / has_tsp3 / has_chb / n_accessory_domains / accessory_domains / pfam_accessions / interpro_signatures`），并按下表扩展：

| 字段 | 类型 | 受控词表 / 说明 |
|---|---|---|
| `accession` | string | 候选 header **原样**（含 `GCA_...\|CTX_..._pos` 前缀；v1 曾出现 accession 截断 bug，见 `docs/T141_20260920_phaded_cys_verification_status.md` §6 第 2 条，本包测试须钉死） |
| `genome` | string | hold 表原样 |
| `source_layer` | enum | `hold_type1_sp_less_v2` |
| `catalytic_domain_state` | enum | `PF10503_supported` / `PF10503_not_detected` / `tool_input_excluded` / `pending` |
| `has_fn3` `has_ig_like` `has_tsp3` `has_chb` | enum | `yes` / `no` / `pending` |
| `n_accessory_domains` | int | ≥0；`pending` 仅当 `has_*` 含 `pending` |
| `accessory_domains` | string | `;` 连接的具体 signature 列表 |
| `pfam_accessions` | string | `hmmscan` 命中 accession 列表 |
| `interpro_signatures` | string | 如启用 InterProScan 路径 |
| `accessory_evidence_independence` | enum | `independent_of_signal_peptide`（固定串，声明本证据不依赖 SignalP） |
| `evidence_flags` | 多值 | 追加式，如 `accessory_domain_supported_independent_of_signalp` / `catalytic_domain_only` / `accessory_domain_not_assessed` |
| `tool_name` `tool_version` `database_sha256` | string | `hmmscan` / `HMMER 3.4 (Aug 2023)` / Pfam-A SHA-256（运行内实测） |
| `run_id` `source_sha256` | string | 运行内绑定 |

- **禁止**：出现任何形式的 `negative`/`absent` 字面量来表达「未检出配件域」；只能用 `not_detected_in_tested_models` + `pending` 语义（与本项目 fail-closed 口径一致）。

### 3.2.7 结果如何只改变证据标记/处置

- **不删除任何候选**：11,610 条 hold 记录在 v1 冻结表中**原样保留**；本包只在新 run 的 `results/` 产出注释层。
- 允许变化：`evidence_flags` 新增配件域支持/未评估标记；命中 ≥1 个**独立于信号肽**的配件域者，可支持把该候选的 `primary_disposition` 由 `remote_homolog_candidate` 上移为 `probable_sequence_homolog`（**仅在 v2 目录构建脚本被单独授权后生效**）。
- 不允许变化：**任何候选不得因配件域未检出而被降级或删除**；不得因本包结果解除 `architecture_criteria`/`common_criteria` 的既有 hold（那是另外的判据）；不得写入 `excluded_input_quality`（除非命中设计文档 §候选目录 第 6 条的三种输入质量情形，且该判定必须有独立输入证据）。

### 3.2.8 线程与资源

- 动态规则同 §2.2；具体线程数 `pending`（启动前实测）。
- 服务器单任务上限 40；`hmmscan --cpu N` 的 N 即为动态值，**不得**沿用历史 `--cpu 40` 作为固定值，也不得沿用 `--cpu 60` 的历史临时授权值（`docs/T141_20260919_phaded_pool_external_evidence_status.md` L37 记录过 60 线程的临时授权）。

### 3.2.9 预估规模与运行时长

| 项 | 值 | 来源 |
|---|---|---|
| 输入规模 | 11,610 条 | 本会话只读复核（§3.2.2） |
| 本项目实测先例 A | Cys 发现层 HMM 打分：候选集 **7.14 s** / 参考集 1.47 s / 混淆集 1.85 s（`hmmsearch`，`--cpu 40`） | `docs/T141_20260917_phaded_cys_discovery_hmm_status.md` L98 |
| 本项目实测先例 B | 全库 100 分片发现层召回墙钟 **499 s ≈ 8 分 19 秒**（单分片 11–37 s，均值 20.1 s） | `docs/T141_20260917_phaded_cys_targeted_recall_status.md` §1、§P1 更正行 |
| 本项目实测先例 C | Pfam-A 全库参考扫描与定向扫描的执行参数已记录（`--cut_ga --cpu 40`），**未记录墙钟时间** | `docs/T141_20260911_phaded_pfam_architecture_01_status.md` L11–L12 |
| 本次预估 | **`pending`** | 11,610 条对 Pfam-A（2.25 GB 库）的 `hmmscan` 本项目**无同规模实测**；先例 A/B 是不同库、不同输入规模，不可直接换算 |

### 3.2.9.1 与 P1 的依赖关系（供排程）

本包的**解释**依赖 P1：只有在参考层信号肽检出率已知的前提下，「配件域支持 vs 折叠误判」的比例才有可比的基准。因此建议排程为 **P1 先授权、P2 后授权**；但两个包的**执行**彼此独立（P2 不需要 P1 的输出作为输入）。

### 3.2.10 明确不属于本次授权的动作

批准 P2 **不**授权：SignalP 执行（P1 是独立批准）；结构预测或 Foldseek（Checkpoint 4）；对 GTDB 分片做 HMMER 重扫（Checkpoint 4，本包只对**已冻结的 11,610 条序列**做域扫描，**不接触 GTDB 分片**）；重建 MAFFT 比对或 HMM（Checkpoint 3）；安装或升级 HMMER/Pfam/InterProScan 环境（Checkpoint 6）；创建 P2 之外的 run/deploy；把配件域证据写成胞外定位或表型结论。

---

## 3.3 P3 —— 29,974 条 Cys 候选的参考锚定残基映射

### 3.3.1 决策问题

「29,974 条 Cys 型候选在文献锚点（PhaZ1 Cys183）所对应的比对列上，究竟是**命中、替换、截断，还是无法判定**？」——即把当前 `inferred_cys_not_verified` 的**整体推断**，升级为逐条的**锚定残基映射状态**，并且**不得把不确定性写成阴性**。

### 3.3.2 目标群体

- **规模：29,974 条**（计划 Task 11 Step 1 第 3 项；接手文档 §4.5 亲核体层）。
- 权威来源表（v1 冻结，只读）：`runs/20260920_phaded_high_confidence_pool_merge_01/results/merged_high_confidence_candidates.tsv`（18 列：`accession / pool_origin / genome / superfamily / signalp_class / profile_best_evalue / profile_evidence_status / lipase_box_state / lipase_box_x1 / catalytic_domain_type / catalytic_residue_verification / nucleophile_type / nucleophile_family / lid_state / sbd_pf06850_binding_state / ahsmg_state / interpro_status / high_confidence`）。
- **精确推导（本会话只读复核，未创建 run）**：

```text
merged_high_confidence_candidates.tsv        (40,690 行; SHA-256 已记录于状态文档 §1.2)
  WHERE nucleophile_type = 'cys'
  => 29,974 行 / 29,974 个唯一 accession / 23,995 个基因组
  （等价条件：superfamily = 'intracellular nPHASCL without lipase box'，两条件在本表中逐行同解）
  catalytic_residue_verification 全部 = 'inferred_cys_not_verified'（29,974/29,974）
  pool_origin 拆分：pool_in 28,482 + pool_external 1,492
```

  - 该 29,974 与接手文档 §4.5 `cys 29,974` 完全一致；`ser 10,634` / `ahsmg 82` 同表同源。
  - **口径交叉核对**：接手文档 §5.2 缺口 #5 与 `docs/T141_20260919_phaded_high_confidence_relabel_status.md` L26 记「Cys 高可信度 **28,534** 条全部 `inferred_cys_not_verified`」。28,534 是**池内 demote 前**口径；28,534 − 52（SignalP-export demote）= **28,482** = 本表 `pool_in` 数，与 `runs/20260920_phaded_three_gaps_01/results/cys_demotion/cys_demotion_summary.json` 的 `moved_to_hold = 52` 一致。三处数字互洽，**引用时必须写明层级**（状态文档 §1.2 口径提醒）。
- 更宽的分母（供分层参考，不属目标集）：Cys 超家族候选全层 **32,904** 条 = 高可信 28,534 + 被 hold 4,370（`docs/T141_20260920_phaded_cys_verification_status.md` §3）。
- 文献锚点（本包唯一锚点）：**PhaZ1 Cys183**（*Ralstonia eutropha* H16；定点突变 C183A / D355A / H388Q 失活，而 **C183S 不失活**；PMID 16233560）→ 模式 `[疏水]-C-Q`（`-1` 疏水来自 Knoll 2009 `cysteine-1` 规则，`+1` Gln 来自锚点上下文 `V182-C183-Q184`）（`docs/T141_20260920_phaded_cys_verification_status.md` §0–§2；接手文档 §2.5）。
- 既有证据层（v1 冻结，可复用）：`runs/20260920_phaded_cys_verification_01/results/cys_nucleophile_motif_evidence.tsv`（139,719 行 = 池内 109,087 + 池外 30,632，size 24,160,768 B），列含 `anchor_exact_vcq` / `anchor_context_hydrophobic_c_q`；该层**每行 `catalytic_activity_verified=false`、`independent_validation=pending`、`hmm_fitted=false`、`family_call_made=false`、`filter_wired=false`**（同文档 §4）。

### 3.3.3 抽样或全集依据

- **全集**（29,974/29,974），不抽样：这是**逐条状态赋值**任务，抽样无法给出每条的处置；且规模与 v1 既有同类运行（139,719 行）同量级，可行性已被先例覆盖。
- 分层键（用于分组报告）：`pool_origin`（`pool_in` 28,482 / `pool_external` 1,492）、`profile_evidence_status`、`superfamily`（本集合内单一值）、`genome`（23,995 个）、以及既有 Cys 模式的三种状态（`anchor_exact_vcq` / `anchor_context_hydrophobic_c_q` / 未命中）。
- 已知分层实测（可作为本包的对照基准，来源 `docs/T141_20260920_phaded_cys_verification_status.md` §3）：Cys 高可信度 28,534 条中 `anchor_exact_vcq` **24,885** + `anchor_context_hydrophobic_c_q` **2,832** → 模式命中率 **0.9714**；Cys 被 hold 4,370 条命中率 0.4776；非 Cys 超家族 2.8%–10%。
- 特异性边界（必须随结果一起引用）：参考面板 Cys 271/276 = **0.9819** vs 非 Cys 22/447 = **0.0492**；精确 `VCQ` 层 0.9094 vs **0.0022** —— **特异性只在样本内测过、无 held-out 集**；F7 条件（在独立非 Cys α/β-水解酶阴性集上评估，命中率 >5% 则降级为纯描述列）**未解除**（同文档 §1、§4 注）。

### 3.3.4 工具与版本

- 工具：项目脚本 `pipeline/scripts/complete_phaded_reference_residue_mapping.py`（计划 Task 11 列为「授权后修改」）+ 既有 Cys 锚定检测实现 `pipeline/scripts/annotate_phaded_cys_nucleophile_motif.py`。
  - 前者的映射机制已在参考层落地并被测试钉死：文献坐标表（`inputs/literature_residue_reference.tsv`，20 行 × 26 列）、**锚定列传递**（只在既有 DED 比对内部进行）、**保守度闸门 ≥0.5**（不达标一律 fail-closed 为 `pending_anchor_column_not_conserved`，**绝不**记为 `not_detected_*`），见 `docs/T141_20260917_phaded_reference_residue_mapping_status.md` §3–§4。
  - 后者提供 `[疏水]-C-Q` 模式与锚点复现的 fail-closed 断言（模式必须在锚点报告位置命中，否则 `SystemExit` 且**不写任何候选行**），见 `docs/T141_20260920_phaded_cys_verification_status.md` §2。
- 版本/出处记录位置：
  - 文献锚点坐标表：`runs/20260917_phaded_reference_residue_mapping_01/inputs/literature_residue_reference.tsv`（size 15,306；SHA-256 前缀 `387c2641b3d1…`，完整值见该 run 的 `input_contract.json`）；
  - 若需要**新建比对**来做候选 ↔ 参考的列传递，则涉及 **MAFFT 7.525**（`docs/T141_20260901_ephaz_phb_boundary_candidate_02_status.md` L15；`docs/T141_20260902_ephaz_mcl_subfamily_candidate_12_status.md` L6），而 MAFFT 输出**不可比特复现**，必须把**具体比对哈希**写进 manifest 并声明 `bit_reproducible=false`（`AGENTS.md` L12；计划 Global Constraints L21）。
  - **本包是否必须新建比对：`pending`** —— 取决于候选序列是否已存在于既有 DED 比对行中；该项必须在授权前的输入盘点中判定，不得默认。
- 既有 DED 比对资产：`runs/20260909_phaded_reference_evidence_01/inputs/raw/ded/aln*.aln`（38 个 family 比对；本会话只读实测该目录存在同名文件）。参考层实测「可用比对 29 个、其中含文献锚行者 3 个（`DED_hfam_4` / `DED_hfam_8` / `DED_hfam_70`）、覆盖 64 条参考」（`docs/T141_20260917_phaded_reference_residue_mapping_status.md` §4 第 2 点）。

### 3.3.5 假阳性竞争者

- **「未命中 = 非 Cys 型」这一错误推论本身**是本包最大的假阳性来源。直接依据：**C183S 不失活**，说明 Ser 可替代 Cys；因此"未命中该模式"**不得**读作"非 Cys 型"的强证据（`docs/T141_20260920_phaded_cys_verification_status.md` §7 第 4 条；接手文档 §10 第 4 问：**Cys 型 29,974 条的催化残基未逐条验证，且 C183S 不失活使"未命中"不能当强阴性**）。
- **无 GxSxG 的 Ser 型混入**：Cys 判型是推断（无 `GxSxG` + `PF06850` + 强 Cys 家族 profile），未排除该混入（接手文档 §5.2 缺口 #5；`filter_summary.json` `honesty_notes.cys_typing`）。
- **单一锚点的过拟合风险**：Cys 型只有 PhaZ1 一条带坐标的实验证据（PhaZ2 有实验阳性但无催化残基坐标；PhaZd 是 Ser 型）；参考面板仍有 **5/276 漏检**（同文档 §7 第 3、4 条）。
- **类泛素/其他 Cys 水解酶背景**：`[疏水]-C`（只有 Knoll 规则、不含 `+1` Gln）在非 Cys 参考上是 **330/447 = 0.7383**，几乎无特异性——特异性主要由 `+1` 的 Gln 提供（同文档 §2 注）。因此任何"含 Cys"的宽模式匹配**不得**作为判别证据。

### 3.3.6 输出 schema

输出必须**显式区分 substitution / truncation / uncertain 三类状态**，并给每行坐标与映射依据：

| 字段 | 类型 | 受控词表 / 说明 |
|---|---|---|
| `accession` | string | 候选 header 原样（含 `GCA_...\|CTX_..._pos`；v1 accession 截断 bug 见 §6 第 2 条） |
| `genome` | string | 合并表原样 |
| `nucleophile_identity` | enum | `cys` / `ser` / `unresolved`（**不得**出现 `ahsmg` 作为亲核体值；AHSMG 是 Ser 的 motif class，计划 Task 3 Step 2） |
| `motif_class` | enum | `GxSxG` / `AHSMG` / `Cys-associated` / `other` / `unresolved` |
| `anchor_reference` | string | `PhaZ1_C183_PMID16233560` |
| `anchor_column` | int | 比对列号；无则为空 |
| `anchor_coordinates` | string | 形如 `183-183`；缺失为空 |
| `mapping_status` | enum | **`match_at_anchor_column`** / **`substitution_at_anchor_column`** / **`truncation_before_anchor_column`** / **`uncertain_no_alignment_row`** / **`uncertain_no_anchor_for_this_family`** / **`uncertain_anchor_column_not_conserved`** / **`uncertain_gap_at_mapped_column`** / `tool_input_excluded` |
| `observed_residue` | char | 映射列上实测残基（`substitution` 时必填，如 `S`） |
| `expected_residue` | char | `C` |
| `mapping_basis` | string | `literature_anchor:…` / `anchor_column_transfer:…` / `literature_pattern:anchor_context_hydrophobic_c_q`；前缀词表沿用参考层实现（`pipeline/scripts/complete_phaded_reference_residue_mapping.py` L67–L68） |
| `conserved_fraction` | float | 该列在比对中的保守度（`>= 0.5` 才允许传递） |
| `catalytic_activity_verified` | bool | 恒 `false`（禁止读作已验证 Cys） |
| `independent_validation` | enum | `pending`（F7 未解除） |
| `evidence_flags` | 多值 | 如 `cys_anchor_substitution_ser_not_negative` / `cys_mapping_uncertain` |
| `run_id` `source_sha256` | string | 运行内绑定 |

- **强制规则**：`substitution_at_anchor_column` 与 `truncation_before_anchor_column` **都不是阴性**；前者必须同时写 `evidence_flags += cys_anchor_substitution_ser_not_negative`，且**不得**触发任何降级。所有 `uncertain_*` 一律为 `pending` 语义，**不得**折叠成 `not_detected`。

### 3.3.7 结果如何只改变证据标记/处置

- **不删除任何候选**：29,974 条在 v1 合并表中原样保留；本包只新增映射证据列。
- 允许变化：`evidence_flags` 新增 `cys_anchor_match` / `cys_anchor_substitution_ser` / `cys_mapping_uncertain` 等标记；`nucleophile_identity` 可由 `unresolved` 收敛为 `cys` 或 `ser`（**当且仅当**锚定列给出明确残基）；据此可支持 `primary_disposition` 由 `remote_homolog_candidate` 上移为 `probable_sequence_homolog`。
- 不允许变化：**不得**因未命中/替换/截断而把任何候选降级或置为 `excluded_input_quality`；**不得**接入高可信度过滤器（既有 `filter_wired=false` 状态不变，`docs/T141_20260920_phaded_cys_verification_status.md` §4）；**不得**据此产生 family 判定（`family_call_made=false` 保持不变）。

### 3.3.8 线程与资源

- 动态规则同 §2.2；具体线程数 `pending`。
- 本包以本地序列/比对处理为主；若含 `hmmsearch`/`hmmscan` 打分，线程按动态值；若需新建比对（MAFFT），则**该动作属另一授权类别**（Checkpoint 3），必须单独批准后再排程。

### 3.3.9 预估规模与运行时长

| 项 | 值 | 来源 |
|---|---|---|
| 输入规模 | 29,974 条 | 本会话只读复核（§3.3.2） |
| 本项目实测规模先例 | v1 同类运行一次处理 **139,719 行**，产出 `cys_nucleophile_motif_evidence.tsv` 24,160,768 B；**未记录墙钟时间** | `docs/T141_20260920_phaded_cys_verification_status.md` §3、§6 |
| 本项目实测参考层先例 | 参考层映射 run 覆盖 723 行 × 53 列输出；**未记录墙钟时间** | `docs/T141_20260917_phaded_reference_residue_mapping_status.md` §7 |
| 本次预估 | **`pending`** | 项目证据中没有本规模的耗时记录；不推算 |

### 3.3.10 明确不属于本次授权的动作

批准 P3 **不**授权：重建 MAFFT 比对或 HMM（Checkpoint 3）；SignalP（P1）；结构预测/Foldseek（Checkpoint 4）；GTDB 分片 HMMER 重扫（Checkpoint 4）；安装/升级环境（Checkpoint 6）；创建 P3 之外的 run/deploy；把映射结果接入过滤器或写成催化活性验证；把 `pending`/`uncertain` 写成阴性。

---

## 3.4 P4 —— 987 条 nPHAMCL-like 候选的结构竞争面板

### 3.4.1 决策问题

「987 条 nPHAMCL-like 候选，在**竞争面板**（多种 α/β-水解酶折叠参照物同时在场）中更接近解聚酶参照物、还是更接近通用脂肪酶/内酯酶参照物？」——即用结构证据回答序列层已定论**不可判别**的家族归属问题。

### 3.4.2 目标群体

- **规模：987 条**（计划 Task 11 Step 1 第 4 项；接手文档 §4.2 `intracellular nPHAMCL = 987（2.4%）`、§9 第 13 项）。
- 权威来源表（v1 冻结，只读）：`runs/20260920_phaded_high_confidence_pool_merge_01/results/merged_high_confidence_candidates.tsv`
  - **精确推导（本会话只读复核）**：`WHERE superfamily = 'intracellular nPHAMCL'` → **987 行 / 987 个唯一 accession / 985 个基因组**；`catalytic_domain_type` 全部为 `undetermined_no_oxyanion`（987/987）。
  - 与 `runs/20260919_phaded_nucleophile_unify_01/results/filter_summary.json` 的 `per_superfamily["intracellular nPHAMCL"].pass = 987`（`total = 6950`、`fail_architecture = 5807`、`fail_common = 156`）一致。
- 为什么必须用**竞争面板**而不是单一 8YNV：`hfam_4`（nPHAMCL）在序列层已被实测**定论不可判别**（接手文档 §5.4 与 §2.6）：
  - challenge 5 家族 × 10 挑战序列中，hfam_4 **3/10 假阳性**：`P24640`（*Moraxella* 脂肪酶，5.4e-12）、`Q02104`（*Psychrobacter* 脂肪酶，4.1e-11）、`Q88N36`（PcaD γ-内酯酶，6.6e-11）（接手文档 §2.6；`docs/T141_20260920_phaded_three_gaps_status.md` §④）；
  - 催化三联体（亲核 Ser 与肘部两侧 G）在 10 条挑战中 **10/10 保守** = 共享的、非判别；唯一判别位点是肘部 `G-V-S-W-G` 的 V/W（挑战 0/10 共享）；根因 `EFFN=0.378`（Dirichlet 先验把 profile 拉向背景）；
  - `--pnone` 重建后 3 假阳性清零但 `eff_nseq=0.00`（阳性 E=0、全挑战 E=999）= **精确匹配过拟合、零召回，假胜利**；
  - **定论**：hfam_4 无法做成判别性家族 HMM（序列信号 = 通用 `PF00561` 折叠），nPHAMCL 改走结构证据；操作者已决定不再修改 HMM。
  - 与 8YNV 的关系：`PDB 8YNV`（*B. thuringiensis* 胞内 PhaZ 首个晶体结构，1.42 Å，Wang 2024，PMID 39592048）是本家族**建族谱系**的参照物（`docs/T141_20260917_project_handoff.md` §578 预注册；`runs/20260919_phaded_with_lipase_deferred_archive_01/results/deferred_tier_notes.md` §3）。单一参照物只能回答"是否像 8YNV"，无法回答"是否**更**像脂肪酶"——这正是 3/10 假阳性的教训。

### 3.4.3 抽样或全集依据

- **全集**（987/987），不抽样：规模小、且每条都需要一个可归属的结构处置。
- 面板设计（**成员清单必须以运行内绑定为准，本文档不虚构 PDB ID**）：
  - 必含：解聚酶参照物 **8YNV**（链 A，298 aa；序列含规范 `GxSxG` `GWSMG` 与 Asp motif `DRDYVV`）与 **8YNW**（`docs/T141_20260917_project_handoff.md` §578 记录的 8YNV/8YNW 结构验证路径）；
  - 必含既有挑战集中的竞争折叠：`P24640`（脂肪酶）、`Q02104`（脂肪酶）、`Q88N36`（PcaD γ-内酯酶）（接手文档 §2.6）——其结构来源与 PDB 映射 **`pending`**；
  - 通用 α/β-水解酶 / 脂肪酶 / 内酯酶折叠的**额外**面板成员的 PDB ID 清单：**`pending`**（必须由操作者在授权前固定为闭集，且与 §3.4.5 的竞争者定义对应）。
  - **预注册判据（提议，需操作者确认后才执行）**：结构通过 = `Foldseek qtmscore ≥ 0.5` 且 `E ≤ 1e-3`；预测质量门 = `pLDDT ≥ 70` 才纳入比对；fail-closed = 无结构 / `pLDDT < 70` / 无 Foldseek 命中 → `not_assessed`，**不得记为阴性**（`docs/T141_20260920_phaded_structure_verification_status.md` §4）。
- 抽样（若操作者要求先跑子集）：规模、分层键（如按 `profile_best_evalue` 分位、按基因组分类分位）与随机种子必须预先固定；**种子与规模写 `pending`**。
- 需显式记录的既有资产与冲突：
  - 既有结构审查输入构造脚本 `pipeline/scripts/prepare_phaded_structure_review.py` **硬编码期望 19 条**（L54–L55：`if len(rows) != 19: raise ValueError`），且其 manifest 声明 `"structure_predictor": {"status": "pending_tool", … "No predictor or structure-comparison executable is available in the bound T141 phb_gtdb environment"}`（L91）；
  - 而项目其他文档实测 **Foldseek 可用**：`${PHB_REMOTE_ROOT}/tools/foldseek_20260912/foldseek/bin/foldseek`，自测 8YNV vs 8YNV → TM 1.00 / E 9.8e-70（`docs/T141_20260920_phaded_structure_verification_status.md` §1；接手文档 §6 错误 #5 明确"Foldseek 未安装"是漏查 `tools/` 的错误结论）。
  - **该冲突必须在授权前解决**：`prepare_phaded_structure_review.py` 的 19 条硬编码必须参数化，其 `pending_tool` 声明必须按实测更正；否则脚本会在 987 条输入上直接 fail-closed 中止。此项**不属本包自动授权范围**，需在批准文本中单独列出。

### 3.4.4 工具与版本

- 结构比对：**Foldseek**
  - 路径：`${PHB_REMOTE_ROOT}/tools/foldseek_20260912/foldseek/bin/foldseek`
  - 版本标识：`463739e0014a1549a527de589102cde98f802f37`；可执行文件与分发包哈希记录在 `runs/20260912_phaded_reference_structure_comparison_01/results/tool_manifest.json`（`docs/T141_20260912_phaded_reference_structure_comparison_01_status.md` L9：安装为 dated user-level binary，4 线程执行）。
- 结构预测（**仅当**需要把无实验结构的候选纳入比对时）：
  - **ColabFold 1.6.1（JAX 0.5.3，GPU0）**；实测 AF2 单序列 296 aa → **84 s**，产出 PDB（`docs/T141_20260920_phaded_structure_verification_status.md` §1）；
  - **硬约束**：AF2 单序列预测**不可用于结构比对** —— 用 8YNV 自身序列做 AF2 单序列预测再与晶体结构比对，**TM 仅 0.20（pLDDT 36.1）**（同文档 §1、§2 第 6 点；`docs/T141_20260920_phaded_convention_and_release_status.md` §9 同）；
  - ESMFold 环境**损坏**（torch `2.12.1+cpu`、缺 `openfold`），修复属"安装"操作，**需单独授权**（同文档 §1、§6 第 1 项）；
  - ColabFold 无 ESMFold 选项；可用模型：alphafold2 / alphafold2_ptm / multimer v1–v3 / deepfold_v1（同文档 §1）。
- 面板参照结构：PDB 8YNV 可下载（实测 `curl https://files.rcsb.org/download/8YNV.pdb` → HTTP 200，1,750,491 B，9,538 ATOM）；8YNV 为 4 链 A/B/C/D（298/296/297/296 aa）（同文档 §1）。
- **本文档不给出任何未在项目文档中出现过的 PDB ID、TM 值或版本号。**

### 3.4.5 假阳性竞争者

必须显式对抗的通用折叠（直接来自 hfam_4 challenge 实测与家族归因）：

| 竞争者类别 | 具名实例 | 实测证据 |
|---|---|---|
| 细菌**脂肪酶**（lipase） | `P24640`（*Moraxella* sp.，E=5.4e-12）、`Q02104`（*Psychrobacter* sp.，E=4.1e-11） | 接手文档 §2.6；`three_gaps_status.md` §④ |
| **γ-内酯酶**（lactonase，PcaD 型） | `Q88N36`（PcaD，E=6.6e-11） | 同上 |
| 通用 **α/β-水解酶折叠**（`PF00561`） | 无单条具名实例；为 hfam_4 的序列信号本体 | 接手文档 §5.4「真正的序列层不可分目前只有 hfam_4 一个实测案例」、§5.4 hfam_4 补充段「催化三联体（Ser + 两侧 G）10/10 挑战共享 = 通用 α/β-水解酶机器、非判别」 |
| **with-lipase** 侧同类问题（仅作机制参照，不属本包对象） | `DED_hfam_2` 池外 1,206,655 条 | 接手文档 §2.7「α/β-水解酶折叠序列不可分」；`AGENTS.md` L11 |

- 竞争面板的判别逻辑必须是**相对**的：候选须**更接近解聚酶参照物而非上述任一竞争折叠**；仅"与 8YNV 相似"不构成通过（3/10 假阳性正是单参照物的产物）。

### 3.4.6 输出 schema

| 字段 | 类型 | 受控词表 / 说明 |
|---|---|---|
| `accession` | string | 候选 header 原样 |
| `genome` | string | 合并表原样 |
| `query_type` | enum | `experimental_structure` / `predicted_structure` / `no_structure` |
| `query_source_id` | string | PDB ID 或预测作业 ID；无则空 |
| `query_plddt` | float | 预测结构的平均 pLDDT；实验结构为空 |
| `panel_member_id` | string | 面板成员 ID（PDB/UniProt，运行内绑定） |
| `panel_member_role` | enum | `depolymerase_reference` / `competitor_lipase` / `competitor_lactonase` / `competitor_generic_ab_hydrolase` |
| `qtmscore` `ttmscore` `alntmscore` | float | Foldseek 输出字段 |
| `evalue` | float | Foldseek E 值（**必须记录搜索空间尺度**） |
| `alignment_coverage` | float | 覆盖度 |
| `structure_call` | enum | `pass_depolymerase_closer` / `fail_competitor_closer` / `not_assessed` |
| `call_basis` | string | 预注册判据版本 + 阈值（如 `panel-v1:qtmscore>=0.5;evalue<=1e-3;plddt>=70`） |
| `evidence_flags` | 多值 | 如 `structure_competition_panel_support` / `structure_not_assessed` |
| `tool_name` `tool_version` | string | `foldseek` / `463739e0014a1549a527de589102cde98f802f37` |
| `run_id` `source_sha256` | string | 运行内绑定 |

- **强制规则**：`not_assessed` 是 fail-closed 的唯一默认值；**不得**出现 `negative`/`absent`。

### 3.4.7 结果如何只改变证据标记/处置

- **不删除任何候选**：987 条在 v1 合并表中原样保留。
- 允许变化：`evidence_flags` 新增结构竞争面板标记；`pass_depolymerase_closer` 者可支持把 `primary_disposition` 由 `function_unresolved` 上移为 `probable_sequence_homolog`（**987 条当前受 `function_unresolved` 保护：设计文档 §特殊家族「hfam_4 / nPHAMCL：降为 `function_unresolved`；现有 987 条不得与功能判别较强的集合混合展示」**——因此即使结构通过，也不得直接把 987 条并入功能判别较强的集合展示）。
- 不允许变化：**不得**因 `fail_competitor_closer` 或 `not_assessed` 而删除或降级；**不得**写入 `excluded_input_quality`；**不得**据此声称 PHB/PHA 降解表型或产生 family 判定（`docs/T141_20260920_phaded_convention_and_release_status.md` §9 要求结构结果只在独立 `structural_validation` 层发布，**不入任何计数**）。

### 3.4.8 线程与资源

- GPU 任务：建议**只用 GPU0**（GPU1 已被他人占用 23.8 GB / 62–100%），CPU 侧 ≤8 线程用于 IO/后处理；线程/GPU 规则同 §2.2（`min(40, nproc − load_1m − 10)`，保留 ~10 核余量）（`docs/T141_20260920_phaded_structure_verification_status.md` §5 注）。
- Foldseek 历史执行使用 **4 线程**（`docs/T141_20260912_phaded_reference_structure_comparison_01_status.md` L9，为历史绑定值，**不作为本次固定值**）。
- 具体线程数/GPU 占用核对结论：`pending`（启动前实测并记录）。

### 3.4.9 预估规模与运行时长

| 项 | 值 | 来源 |
|---|---|---|
| 输入规模 | 987 条 | 本会话只读复核（§3.4.2） |
| Foldseek 比对 | **分钟级**（项目表格原文） | `docs/T141_20260920_phaded_structure_verification_status.md` §5 |
| 结构预测单条实测 | ColabFold AF2 单序列 296 aa → **84 s** | 同文档 §1 |
| 文档内既有换算（作为方法先例，非本包数字） | 1,550 条 → ~36 GPU-h；4,013 条 → ~94 GPU-h（AF2 单序列）；ESMFold 若修好（~5 s/条）→ 1,550 条 ~2.2 GPU-h；AF2+MSA（2–5 min/条）→ 1,550 条 ~1–2 GPU-日；**全量 1.2M → ~324 GPU-日（不可行）** | 同文档 §5 |
| 本次预估（987 条） | **`pending`** | 项目**未实测** 987 条规模；本文档不据 84 s/条 自行换算成工时 |

### 3.4.10 明确不属于本次授权的动作

批准 P4 **不**授权：修复 ESMFold 环境（安装 `openfold` + GPU 版 torch，属 Checkpoint 6）；批准 AF2+MSA 联网模式；SignalP（P1）；HMMER/Pfam 域扫描（P2）；HMM 重建（Checkpoint 3）；GTDB 分片 HMMER 重扫（Checkpoint 4）；with-lipase 1.2M 层（P5 是独立批准）；创建 P4 之外的 run/deploy；把结构结果计入任何候选/高可信度计数；产生 family 判定或表型结论；修改 `pipeline/scripts/prepare_phaded_structure_review.py` 的 19 条硬编码（该项需在批准文本中单独点名）。

---

## 3.5 P5 —— with-lipase deferred 池的分层代表抽样

### 3.5.1 决策问题

「在 1,206,655 条 with-lipase deferred 命中中，能否用**分层代表抽样 + 结构召回**，以可外推的方式估计其中『真带 lipase box 的胞内解聚酶』的比例，并召回经结构证据支持的代表性成员？」——即在不删除、不计数的前提下，为这个百万级噪声层建立唯一既定的召回路径。

### 3.5.2 目标群体

- **池规模：1,206,655 条**；**禁删、禁入任何主候选计数**（`AGENTS.md` L11；接手文档 §2.7、§4.1；`docs/T141_20260921_phaded_comprehensive_review.md` §7.2 陷阱 #5：「用 1,206,655（trained 覆盖 discovery 规则后），不用 1,238,812（规则前桶）」）。
- 权威来源表（v1 冻结 deferred，只读）：`runs/20260919_phaded_with_lipase_deferred_archive_01/results/pool_external_with_lipase_deferred.tsv`
  - **本会话只读复核**：**1,206,655 数据行**（含表头 1,206,656 行；文件 size 183,156,801 B）；列为 `protein_id / genome / superfamily / discovery_best_family / discovery_families_hit / discovery_best_evalue / trained_best_model / trained_best_evalue / override_by_trained / model_layer`（10 列）。
  - 同 run 另两件：`pool_internal_with_lipase_discovery_unique.tsv`（5,711 条，池内发现层 unique 归属）、`pool_internal_with_lipase_filter_consistent.tsv`（13 条，与高可信过滤器同优先级口径，`fail_common = 13`、`pass = 0`）。
  - 计数器一致性：`results/deferred_summary.json` `with_lipase_pool_external_deferred = 1206655`；池外 8 超家族合计 1,512,873（`deferred_tier_notes.md` §1）。
  - 分配规则（复现锚点）：trained 改判**当且仅当** `trained_best_evalue < discovery_best_evalue`（严格小于）；hfam_2-best 1,238,812 行 − trained 更强 32,157 行 = **1,206,655**（`deferred_tier_notes.md` §2）。
- **序列取回配方**（deferred TSV 不含序列）：用 `protein_id`（`GCA_...|CTX_..._pos`）从 GTDB 分片或全库 FASTA 按头名提取；分片体系与 SHA-256 以其 `input_contract` 为准（`deferred_tier_notes.md` §3）。
- **分层键**（计划 Task 11 Step 1 第 5 项要求的六项）与其**当前可得性**（必须如实区分，不得假装六项都现成）：

| 分层键 | 在该 TSV 中的可得性 | 说明 |
|---|---|---|
| 序列簇 | **需新建**（`pending`） | TSV 无簇列；需按序列聚类（阈值与工具 `pending`） |
| 分类（taxonomy） | **需 join**（`pending`） | TSV 只有 `genome`；GTDB taxonomy 源路径/版本 `pending` |
| 长度 | **需 join**（`pending`） | 需按取回配方提取序列后计算 |
| 完整度 | **需 join**（`pending`） | 需基因模型/完整度元数据，源 `pending` |
| 得分 | **现成** | `discovery_best_evalue`、`trained_best_evalue` |
| 架构 | **部分现成** | `superfamily`、`discovery_best_family`、`discovery_families_hit`；配件域架构需另 join（`pending`） |

  → 因此本包的抽样设计必须**先固定**上述 join 的权威源与版本，否则分层不成立；未固定的部分写 `pending`。

### 3.5.3 抽样或全集依据

- **分层抽样**（非全集）：全集 1.2M 的结构预测已被实测判定**不可行**（~324 GPU-日，`docs/T141_20260920_phaded_structure_verification_status.md` §5）。
- 抽样设计必须预先书面固定并在 `input_contract.json` 中绑定：分层键与层权重、每层配额、随机种子、以及**抽样框 SHA-256**。**本文档不预设抽样规模与种子——两者均写 `pending`。**
- **外推性设计参照（已有实测先例，必须复用其验证方法）**：2026-09-20 的 1,001 条抽样采用 deferred TSV（按 `genome accession` 排序）上的**等步长系统抽样**——首行 index 0、末行 1,206,000、步长恒为 **1,206**（1,001 × 1,206 ≈ 全表，无跳区）；代表性经机器判定：样本按 `discovery_best_evalue` 的排名中位 **605,238（≈第 50.2 百分位）**，进入 top-1000 者 **1 条**（均匀期望 0.83 条），**Spearman(行序, E 值秩) = 0.0074**（行序与 E 值轴独立），`GCA_000` 前缀占比 样本 0.0030 / 群体 0.0026（135/194 个前缀覆盖）（`docs/T141_20260920_phaded_structure_verification_status.md` §2.2；校验脚本 `runs/20260919_phaded_pool_external_evidence_01/plan/verify_pilot_blastp_columns.py`，只读，全 PASS 才退出 0）。
- **预筛建议（文档既有方案，需操作者确认）**：三段漏斗 = Stage 1 严格 E 值预筛（**E<1e-50 → 1,550 条**，推荐保守口径；**E<1e-30 → 4,013 条**，宽口径备选）→ Stage 2 结构预测（仅 Stage 1 子集，GPU0）→ Stage 3 Foldseek（8YNV 链 A 为 query vs 预测结构库）→ Stage 4 通过者走证据层/校准 gate（仍 candidate-only，不入现有计数）（同文档 §3）。
- **序列层口径警告（必须随结果引用）**：`hfam_2` 对 8YNV 的得分是 **E 5.8e-165、bitscore 533.6**，而 1.2M 中 **99.6% 是 E>1e-30 的弱命中**；E 值分布桶：`<1e-100` 1,031 / `1e-100–1e-50` 519 / `1e-50–1e-30` 2,463 / `1e-30–1e-10` 676,148 / `1e-10–1e-5` 526,223（同文档 §2 第 4、5 点）。另：对**单条参考**（8YNV 链 A，298 aa）做 blastp 的 E 值按单序列搜索空间计算，**比 GTDB 规模宽松约 8.4×10⁵ 倍**（log10 = 5.92），**不可**与任何全库 E 值比较，**不构成过滤器，无命中不得读作阴性**（同文档 §2.2；换算 GTDB 规模：E<1e-5 ↔ 578/1,001，E<1e-10 ↔ 174/1,001）。

### 3.5.4 工具与版本

- 抽样与分层：**项目本地只读脚本**（无外部计算工具）。建议 run 名 `20260928_phaded_with_lipase_deferred_sample_v2_01`（**建议值，计划未给定**；日期前缀改为实际启动日）。
- 序列取回：按 `deferred_tier_notes.md` §3 配方，用 `protein_id` 从 GTDB 分片/全库 FASTA 按头名提取；分片路径与 SHA-256 以 `runs/20260901_formal_frozen_scan_13/inputs/...` 的 `input_contract` 为准（**服务器权威件，本地无副本**；状态文档 §1.1 记「服务器全库扫描（scan 13）… 本地无副本」）。
- 结构召回（Stage 2/3，**需单独授权**）：
  - Foldseek：路径与版本同 §3.4.4（`${PHB_REMOTE_ROOT}/tools/foldseek_20260912/foldseek/bin/foldseek`，`463739e0…`）；
  - 结构预测：ColabFold 1.6.1（GPU0）；ESMFold 环境损坏（修复需授权）；**AF2 单序列不可用于结构比对（同序列 TM 0.20 / pLDDT 36）**（`docs/T141_20260920_phaded_structure_verification_status.md` §1–§2）。
- 既有试点脚本（服务器 `runs/20260919_phaded_pool_external_evidence_01/plan/`）：`pilot_extract_deferred_sample.sh`、`pilot_analyze_deferred_sample.sh`、`pilot_foldseek_validate.sh`、`verify_pilot_blastp_columns.py` 等（同文档 §3）。**授权后必须移植进新的 dated deploy**，不得直接运行历史 run 下脚本（`AGENTS.md` L6）。
- **列映射不变式（防复发）**：凡 `blastp` 自定格式产物，必须先跑列不变式检查再解读——同一列错位在本项目中**已复发两次**（同文档 §2.2 更正记录与 §2.2「复发提示」）。

### 3.5.5 假阳性竞争者

- **层内主体噪声**：1.2M 中 **99.6% 为 E>1e-30 的弱命中**；`deferred_tier_notes.md` §4 原文：本层全部为 `discovery_hmm_uncalibrated` / recall-only 证据，实测该发现层对 563 条 x₁ 混淆嫌疑命中 **443 条（78.69%）**，**无判别力**（`AGENTS.md` L9 同）。
- **生物学竞争者**：细菌脂肪酶/酯酶——hfam_2 的判别特征「带 lipase box 的胞内 SCL-PHB 解聚酶」与细菌脂肪酶/酯酶的通用信号同源（同为 α/β-水解酶折叠 + 同样肘部基序），**一维序列层面没有可判别的家族特异指纹**；任何序列 HMM 只有两个结局：先验主导成"近背景"（海量假阳性）或收紧成精确匹配（零召回，参见 hfam_4 `--pnone` 教训）（接手文档 §2.7 第 1 点）。
- **实验锚点不足**：本家族实验阳性只有 **1 条**（`EAO52570.1`，出自 Knoll 2009 单一出处 PMID 19296857）；台账 29 条参考中 28 条 `annotation_only`；远低于 gate 的「≥3 独立阳性」，也无跨属、跨论文的活性表征（接手文档 §2.7 第 3 点）。
- **高可信度 0 通过**：池内 13 条全部 `common_criteria` hold、池外 0 通过 → **0 条进入 40,690**（唯一一个 0 通过的胞内超家族）（接手文档 §2.7 第 4 点；`deferred_tier_notes.md` §1）。
- **单参照物 blastp 的伪证据**：抽样实测 987/1,001 在 E≤1e-3 有命中（单序列空间），其中 E<1e-5 = **954**、最小 E = **1.16e-122**、pident 中位 24.254%、最大 54.138% → 结论是「对单条参考做 blastp **不具判别力、不能用作过滤器**」，而**不是**"命中即阳性"（同文档 §2.2 更正后的科学解读）。

### 3.5.6 输出 schema

- 抽样清单（`deferred_stratified_sample.tsv`；**只列抽样成员，绝不导出或复制整层**）：

| 字段 | 类型 | 受控词表 / 说明 |
|---|---|---|
| `protein_id` | string | deferred TSV 原样（`GCA_...\|CTX_..._pos`） |
| `genome` | string | 原样 |
| `superfamily` | string | 原样（本层主体为 `intracellular nPHASCL with lipase box`） |
| `discovery_best_family` `discovery_families_hit` | string | 原样 |
| `discovery_best_evalue` `trained_best_evalue` | float | 原样 |
| `override_by_trained` `model_layer` | string | 原样 |
| `stratum_id` | string | 层标识（由预注册分层键生成） |
| `stratum_keys` | string | 六键取值：`cluster=…;taxonomy=…;length=…;completeness=…;score_bin=…;architecture=…`（未可得键写 `pending`） |
| `selection_method` | enum | `systematic_step` / `stratified_random` / `census` |
| `selection_seed` | int | 抽样种子（`pending` 至操作者固定） |
| `selection_index` | int | 在抽样框中的序号（系统抽样的可审计锚点） |
| `sequence_sha256` | string | 取回序列的 SHA-256（未取回写 `pending`） |
| `deferred_layer` | string | 固定串，声明 `not_in_main_counts` |
| `run_id` `frame_sha256` | string | 运行内绑定（抽样框 SHA-256） |

- 结构召回输出：沿用 §3.4.6 的字段（`qtmscore` / `ttmscore` / `alntmscore` / `evalue` / `plddt` / `structure_call` ∈ {`pass` / `fail` / `not_assessed`}），并额外记录 `search_space_scale`（声明单序列 vs 全库尺度）与 `criteria_version`。
- **语义强制**：抽样是为**估计比例**与**召回**，不是为筛选；`not_assessed` 为 fail-closed 默认值；**层内任何计数都不得进入主候选计数**。

### 3.5.7 结果如何只改变证据标记/处置

- **不删除任何一条**：1,206,655 条在 deferred 归档中原样保留（`AGENTS.md` L11「禁止删除」；接手文档 §9 第 6 项）。
- 允许变化：被结构召回并预注册判据通过者，可在**独立的 `structural_validation` 层**新增条目与 `evidence_flags`（`docs/T141_20260920_phaded_convention_and_release_status.md` §9：「若日后执行…在**独立 `structural_validation` 层**发布」）；通过者再走证据层与校准 gate（`AGENTS.md` L11）。
- 不允许变化：**本层不进入任何候选/高可信度计数**；`primary_disposition` 的 `deferred_structure_review` **不因本包自动撤销**——通过者必须先经证据层与校准 gate，且该迁移需另行授权；不得把 `not_assessed` 或未命中写成阴性；不得据此产生 family 判定。

### 3.5.8 线程与资源

- 动态规则同 §2.2；抽样本身为本地 IO 任务（线程数 `pending`，启动前实测）。
- 结构阶段：**只用 GPU0**（GPU1 被他人占用）；CPU 侧 ≤8 线程；服务器 GPU/CPU 占用与线程值 `pending`（启动前记录 `nproc`、`/proc/loadavg`、计算出的线程数与 dated deploy 身份，计划 Task 11 Step 5）。

### 3.5.9 预估规模与运行时长

| 项 | 值 | 来源 |
|---|---|---|
| 池规模 | 1,206,655 条（183,156,801 B） | 本会话只读复核 + `deferred_summary.json` |
| 序列提取实测先例 | 抽样 1,001 条、扫 100 分片（268 GB）→ **293 s**；文档注明"成本在扫分片，故全量 1.2M 提取耗时相同" | `docs/T141_20260920_phaded_structure_verification_status.md` §1–§2 |
| 结构预测先例 | AF2 单序列 **84 s/条**；1,550 条 ~36 GPU-h；4,013 条 ~94 GPU-h；全量 1.2M ~324 GPU-日（不可行）；ESMFold 若修好 ~5 s/条 → 1,550 条 ~2.2 GPU-h | 同文档 §5 |
| Foldseek 比对先例 | **分钟级** | 同文档 §5 |
| 抽样与分层本身 | **`pending`** | 项目无该步耗时记录 |
| 本次抽样规模与结构阶段时长 | **`pending`** | 抽样规模未定（§3.5.3）；结构阶段时长取决于抽样规模与预筛口径，本文档不自行换算 |

### 3.5.10 明确不属于本次授权的动作

批准 P5 **不**授权：删除、移动、归档或改写 deferred 层任何一行（**永久禁止**）；把本层任何计数写入候选/高可信度计数或图表候选数；SignalP（P1）；HMMER/Pfam 域扫描（P2）；HMM 重建（Checkpoint 3）；GTDB 分片 HMMER 重扫（Checkpoint 4，含"序列层数据库规模搜索"——文档已明确该形式若日后要加，必须**先预注册阈值并单独授权**）；修复 ESMFold 环境或批准 AF2+MSA 联网（Checkpoint 6）；创建 P5 之外的 run/deploy；用单条参考 blastp 作为过滤器或把无命中读作阴性。

---

## 4. 排程与相互依赖（供操作者决定批准顺序）

| 依赖 | 说明 |
|---|---|
| P1 → P2（解释层） | P2 的「配件域支持 vs 折叠误判」比例需要 P1 给出的参考层信号肽检出率作为基准；两者**执行**互不依赖，**解释**有先后 |
| P1/P2 → 候选层 type-1 门标定 | 任一结果都只改证据口径；门本身的改写需另一次显式决策 |
| P3 独立 | 与 P1/P2 无输入依赖；但若需新建 MAFFT 比对，触发 Checkpoint 3（另需授权） |
| P4 独立 | 依赖结构工具可用性；`prepare_phaded_structure_review.py` 的 19 条硬编码与 `pending_tool` 声明冲突**必须在授权前解决** |
| P5 独立且隔离 | 与 P4 共享 Foldseek/预测工具但**层与计数完全隔离**；P5 结构通过的条目**不得**回流 P4 的 987 条计数体系 |
| 全部 | 任何包都**不**授权 GTDB 全库重扫；任何包都**不**产生 family 判定；任何包都**不**删除候选 |

---

## 5. 逐数字溯源表（本文档引用的全部计数）

| 数字 | 来源（路径 / 章节） |
|---|---|
| 723 参考（30 `experimental_positive` + 693 `annotation_only`） | `docs/T141_20260921_phaded_comprehensive_review.md` §2.1；表 `runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv`（本会话只读复核：723 数据行 × 26 列） |
| 46 profile = 9 trained + 37 reference-only | 接手文档 §2.1；`docs/T141_20260928_phaded_evidence_model_redesign_status.md` §1.3 |
| 参考层 0 行提及 signal；284 条 type-1 参考 `domain_architecture` 全 `pending` | 接手文档 §5.2 缺口 #4 |
| 候选宇宙 109,087 = ePhaZ 70,682 + iPhaZ 38,405 | 接手文档 §3、§4.1；状态文档 §1.2 |
| 池内高可信 36,611（demote 前）/ 36,559（demote 后） | 状态文档 §1.2；接手文档 §2.4；`runs/20260920_phaded_three_gaps_01/results/cys_demotion/cys_demotion_summary.json` |
| 合并高可信 40,690（0 重叠；池内 36,559 + 池外 4,131） | 状态文档 §1.2（含 SHA-256）；接手文档 §2.4 |
| 亲核体 cys 29,974 / ser 10,634 / ahsmg 82 | 接手文档 §4.5；本会话只读复核 `merged_high_confidence_candidates.tsv` |
| type-1 hold 48,038 / 67,720（71%） | 接手文档 §5.3 第 1 条；`runs/20260919_phaded_nucleophile_unify_01/results/filter_summary.json` → `fail_localization = 48038`、`total = 67720` |
| 11,610（SP-less type-1，本包目标） | 接手文档 §5.3 第 4 条、§9 第 2 项；计划 Task 11 Step 1 第 2 项；**本会话只读复核单键复现**：`hold_candidates.tsv` 三键过滤 = 11,610 |
| hold 全表 71,860；type-1 子集 62,836；架构 hold 14,433；common 365 | `runs/20260919_phaded_nucleophile_unify_01/results/filter_summary.json`（`total_held`、`per_superfamily`）；本会话只读复核 |
| type-1 高可信 5,018（池内 4,884 + 池外 134） | `runs/20260920_phaded_three_gaps_01/results/type1_sbd/type1_sbd_annotation_summary.md` §4–§6；接手文档 §4.2 |
| 配件域：FN3 1,278 / Ig-like 344 / TSP3 74 / CHB 119 / ≥1 配件域 1,430（29.3%） | `type1_sbd_annotation_summary.md` §6 |
| Cys 高可信 28,534（池内 demote 前）；28,482（池内 demote 后）+ 1,492（池外）= 29,974 | 接手文档 §5.2 缺口 #5；`docs/T141_20260919_phaded_high_confidence_relabel_status.md` L26；`filter_summary.json` `pass = 28534`；本会话只读复核合并表 `pool_origin` 拆分 |
| Cys 超家族候选 32,904 = 28,534 + hold 4,370 | `docs/T141_20260920_phaded_cys_verification_status.md` §3 |
| Cys 模式命中：24,885 + 2,832（命中率 0.9714）；hold 组 0.4776 | 同文档 §3 |
| Cys 面板特异性 271/276 = 0.9819 vs 22/447 = 0.0492；精确 `VCQ` 0.9094 vs 0.0022；`[疏水]-C` 330/447 = 0.7383 | 同文档 §1–§2 |
| Cys 层规模先例：139,719 行；产物 24,160,768 B | 同文档 §3、§6 |
| `inferred_cys_not_verified`（Cys 高可信全量） | 同文档 §4；`filter_summary.json` `honesty_notes.cys_typing` |
| nPHAMCL 987（985 基因组；`catalytic_domain_type` 全 `undetermined_no_oxyanion`） | 接手文档 §4.2；`filter_summary.json` `pass = 987`、`total = 6950`；本会话只读复核 |
| hfam_4 challenge 3/10 假阳性：P24640 5.4e-12 / Q02104 4.1e-11 / Q88N36 6.6e-11 | 接手文档 §2.6；`docs/T141_20260920_phaded_three_gaps_status.md` §④ |
| hfam_4 EFFN 0.378；`--pnone` 后 `eff_nseq=0.00`（阳性 E=0、挑战全 E=999） | 接手文档 §2.6、§5.4 |
| 结构层覆盖 109,068/109,087 = 99.98% `not_tested` | 接手文档 §4.4、§5.2 缺口 #2 |
| 结构证据 19 条 review/support/conflict；系统发育审查 42 条 | 接手文档 §2.3、§4.4 |
| deferred 池 1,206,655（不含序列，10 列） | `runs/20260919_phaded_with_lipase_deferred_archive_01/results/deferred_summary.json`；`deferred_tier_notes.md` §1；**本会话只读复核 1,206,655 数据行** |
| 池内 with-lipase 5,711（discovery unique）/ 13（filter consistent）；池外合计 1,512,873 | `deferred_tier_notes.md` §1 |
| 规则前桶 1,238,812 − trained 更强 32,157 = 1,206,655 | `deferred_tier_notes.md` §2 |
| deferred E 值桶：1,031 / 519 / 2,463 / 676,148 / 526,223；E<1e-50 = 1,550；E<1e-30 = 4,013；99.6% E>1e-30 | `docs/T141_20260920_phaded_structure_verification_status.md` §2 第 4–5 点 |
| hfam_2 对 8YNV：E 5.8e-165、bitscore 533.6 | 同文档 §2 第 4 点 |
| 8YNV：链 A 298 aa，含 `GWSMG` / `DRDYVV`，1.42 Å，PMID 39592048；可下载 1,750,491 B / 9,538 ATOM | 同文档 §1；`docs/T141_20260917_project_handoff.md` §578 |
| Foldseek 版本 `463739e0014a1549a527de589102cde98f802f37`；自测 8YNV vs 8YNV TM 1.00 / E 9.8e-70 | `docs/T141_20260920_phaded_structure_verification_status.md` §1；`docs/T141_20260912_phaded_reference_structure_comparison_01_status.md` L9 |
| AF2 单序列 296 aa → 84 s；同序列 TM 0.20 / pLDDT 36.1；ESMFold 破损（torch 2.12.1+cpu、缺 openfold） | `docs/T141_20260920_phaded_structure_verification_status.md` §1–§2 |
| 3 段漏斗算力：1,550 → ~36 GPU-h（AF2）/ ~2.2 GPU-h（ESMFold）；4,013 → ~94 GPU-h；全量 1.2M → ~324 GPU-日；Foldseek 分钟级 | 同文档 §5 |
| 预注册判据：qtmscore ≥0.5 且 E ≤1e-3；pLDDT ≥70；fail-closed → `not_assessed` | 同文档 §4 |
| 抽样 1,001 条提取 293 s（扫 100 分片 / 268 GB）；系统抽样步长 1,206；Spearman 0.0074；top-1000 命中 1 条（期望 0.83） | 同文档 §1–§2.2 |
| blastp vs 8YNV 单序列：987/1,001 E≤1e-3；E<1e-5 = 954；最小 1.16e-122；pident 中位 24.254%、最大 54.138%；单序列空间宽松 ~8.4×10⁵ 倍 | 同文档 §2.2 |
| SignalP 6.0h：29,521 条 ~31 min；108,736 条 ~2 h 41 min | `docs/T141_20260920_phaded_three_gaps_status.md` L17；`docs/T141_20260916_phaded_full_library_signalp_01_status.md` L17 |
| SignalP 29,521 分类分布：OTHER 29,469 / SP 49 / LIPO 1 / TAT 2 | `docs/T141_20260920_phaded_three_gaps_status.md` §① |
| HMMER 3.4 (Aug 2023)；Pfam-A size 2,246,909,846 B + SHA-256 | `docs/T141_20260911_phaded_pfam_architecture_01_status.md` L10–L11 |
| MAFFT 7.525（不可位复现） | `docs/T141_20260901_ephaz_phb_boundary_candidate_02_status.md` L15；`docs/T141_20260902_ephaz_mcl_subfamily_candidate_12_status.md` L6；`AGENTS.md` L12 |
| InterProScan 5.78-109.0 | `docs/T141_20260914_phaded_interpro_full_01_status.md` L29；`docs/superpowers/plans/2026-09-14-phaded-evidence-completion.md` L9 |
| HMMER 全库先例：9 profile × 100 shard = 41 min；单 profile ≈5–10 min | 接手文档 §2.6 |
| Cys 发现层打分：7.14 s / 1.47 s / 1.85 s | `docs/T141_20260917_phaded_cys_discovery_hmm_status.md` L98 |
| 全库 100 分片召回墙钟 499 s（单分片 11–37 s，均值 20.1 s） | `docs/T141_20260917_phaded_cys_targeted_recall_status.md` §1、P1 行 |
| 发现层 x₁ 混淆命中 443/563 = 78.69% | `AGENTS.md` L9；`deferred_tier_notes.md` §4 |
| v2 目录 `phaded_candidate_catalog_v2.tsv` 未创建 | 状态文档 §1.3（`pending-authorization`） |
| 52 条 demote：36,611→36,559；40,742→40,690 | 接手文档 §2.4；`cys_demotion_summary.json`；`docs/T141_20260920_phaded_three_gaps_status.md` §① |
| 参考层映射先例：723 行 × 53 列；锚定列可用比对 29 个 / 含锚 3 个 / 覆盖 64 条参考；保守度闸门 0.5 | `docs/T141_20260917_phaded_reference_residue_mapping_status.md` §4、§7 |

---

## 6. 已登记的开放问题（必须在授权前解决或明确接受为 `pending`）

| # | 问题 | 影响包 | 处置 |
|---:|---|---|---|
| 1 | 参考层 15 个证据质量/独立性列（E3/E2/E1/A、`independence_group` 等）**尚未落到真实台账** | P1 | P1 只能用台账现有 26 列分层；若要求按证据分级分层，须同时授权 Task 4 Step 6 的台账 curation |
| 2 | Pfam-A **发布版本字符串**未记录（只有路径/size/SHA-256） | P2 | `pending`；授权后在 `input_contract.json` 记录实际库版本 |
| 3 | P3 是否需要**新建 MAFFT 比对**（触发 Checkpoint 3） | P3 | 授权前完成输入盘点后判定，不得默认 |
| 4 | `pipeline/scripts/prepare_phaded_structure_review.py` 硬编码 19 条（L54–L55）且 manifest 声明 `pending_tool`（L91），与"Foldseek 可用"的实测结论冲突 | P4 | 必须在批准文本中单独点名修改；否则 987 条输入直接 fail-closed 中止 |
| 5 | `blastp` 自定格式的列映射错位**已复发两次** | P5 | 任何序列层产物必须先跑列不变式检查（`verify_pilot_blastp_columns.py` 范式）再解读 |
| 6 | with-lipase 竞争面板的 PDB 成员闭集、抽样规模、抽样种子、六项分层键的 join 源与版本 | P4 / P5 | 全部 `pending`；须由操作者在授权前固定并写入批准文本与 `input_contract.json` |
| 7 | 计划仅给定 3 个 run 名（P1/P3/P4）；P2/P5 的 run 名与日期前缀 | P2 / P5 | 本文档给建议名；实际 ID 由协调者在授权后按启动日确定并核对路径 |
| 8 | 历史 60 线程（SignalP）与 70 线程（gRodon）为临时授权 | 全部 | 不得沿用；一律按 §2.2 动态实测 |

---

## 7. 未执行声明（诚实边界）

1. 本文档编写期间**未执行**任何计算工具：无 SignalP、无 HMMER/`hmmscan`/`hmmsearch`、无 MAFFT、无 Foldseek、无 ESMFold/AlphaFold/ColabFold、无 gRodon、无 GTDB 扫描。
2. **未创建**任何 `runs/<run_id>/` 或 `deploy/<run_id>/`；未写入、未移动、未删除 `runs/`、`results/`、`deploy/` 下任何文件。
3. 本文档对 `runs/` 下冻结 TSV 的访问**全部为只读**：读取表头、统计行数、按既有列做只读分组计数（用于验证 11,610 / 29,974 / 987 / 1,206,655 等目标群体推导）。这些只读复核结果**不构成运行产物**，授权后必须在新的 `runs/<run_id>/` 中重新计算并绑定 SHA-256。
4. 未安装或变更任何环境；未修改服务器配置；未 push；未 `git add`/`commit`/`reset`/`checkout`。
5. 本文档是 **candidate-only** 记录：所有 HMM/profile/domain/motif/SignalP/结构/系统发育/分类与生长速率证据只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**；结构通过者也只在独立 `structural_validation` 层发布，不入任何候选/高可信度计数。

---

*本文档为授权准备件（preparation artifact），不含任何已执行运行的结论。五个包的状态一律为 `pending-authorization`；每个计算类别（SignalP、HMM building、Foldseek/结构预测、GTDB searching）需要分别的显式批准。*
