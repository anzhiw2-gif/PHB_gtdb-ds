# T141 PhaDED 证据模型与候选目录重构 —— 实施状态文档

> **计划**：`docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign.md`
> **设计**：`docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md`
> **接手文档**：`docs/T141_20260921_phaded_comprehensive_review.md`（v1 冻结口径）
> **本文档性质**：candidate-only 实施记录。所有 HMM/profile/domain/motif/SignalP/结构/系统发育/分类与生长速率证据只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。

---

## 0. 一句话状态

2026-09-28 重构按计划实施中：v1 冻结产物 **未被读写改动**（哈希实测复核一致），v2 组件以「新文件 + 新测试」方式落地，任何真实重算（SignalP / 结构 / HMMER 全库 / gRodon / 建 HMM）**仍为未授权状态**，全部停在授权 checkpoint 之前。

---

## 1. Task 1：权威图谱冻结（本文件即产出）

### 1.1 不可变 v1 权威表

| Topic | Frozen authority | Status |
|---|---|---|
| candidate universe | `runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv` | frozen-v1 |
| pool-in high confidence | `runs/20260919_phaded_nucleophile_unify_01/results/high_confidence_candidates.tsv` | frozen-v1 |
| final v1 merge | `runs/20260920_phaded_high_confidence_pool_merge_01/results/merged_high_confidence_candidates.tsv` | frozen-v1 |
| with-lipase deferred | `runs/20260919_phaded_with_lipase_deferred_archive_01` | frozen-deferred |
| v2 candidate catalog | not created | pending-authorization |

补充冻结件（v1 口径引用时需要，同样只读）：

| Topic | Frozen authority | Status |
|---|---|---|
| 52 条 SignalP-export Cys demote | `runs/20260920_phaded_three_gaps_01/results/cys_demotion/` | frozen-v1 |
| 候选映射 + 亚型矩阵（109,087） | `runs/20260916_phaded_full_library_signalp_01/results/phaded_full_library_localization_subtype_matrix.tsv` | frozen-v1 |
| 池外 strong 层证据补齐 | `runs/20260919_phaded_pool_external_evidence_01/` | frozen-v1 |
| 校准台账（723 参考） | `runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv` | frozen-v1 |
| 分类权威层（46 profile） | `runs/20260917_phaded_classification_authority_01/` | frozen-v1 |
| 服务器全库扫描（scan 13） | 服务器 `runs/20260901_formal_frozen_scan_13/`（本地无副本） | frozen-v1（remote） |
| 生长速率分析 | `runs/20260920_phaded_grodon_growth_01/` | frozen-v1 |

### 1.2 v1 权威件实测哈希与行数（2026-09-28 本会话实测，非转录）

| 文件 | bytes | 数据行数 | SHA-256 |
|---|---:|---:|---|
| `runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv` | 10,957,989 | 109,087 | `F4A9D99039EF5DC46FBEA5B90DD98FD61ADC2C65AE1DA700C53976CB5B7DF08C` |
| `runs/20260919_phaded_nucleophile_unify_01/results/high_confidence_candidates.tsv` | 12,879,302 | 36,611 | `7772B34C6A9FFA79B5BEDFB2446E2301DCD46DB936D9267E5396F1E81E81101F` |
| `runs/20260920_phaded_high_confidence_pool_merge_01/results/merged_high_confidence_candidates.tsv` | 14,286,866 | 40,690 | `7B6F9208CFCC175932AA21DA5114DBA10281D6311966DADFAAA3EB2C99715BDB` |

**行数交叉核对（与接手文档 §2/§3 逐项一致）**：

- 候选宇宙 **109,087** = ePhaZ 70,682 + iPhaZ 38,405 ✅
- 池内高可信 **36,611**（demote 前口径；demote 后 36,559 见 `cys_demotion/`）✅
- 合并高可信 **40,690**（池内 36,559 + 池外 strong 4,131，0 重叠）✅

> **口径提醒**：`high_confidence_candidates.tsv` 是 **36,611** 的 demote 前快照；36,559 是 demote 修订版口径。引用时必须写明层级，避免 §7.2 第 3 条陷阱（30,026 vs 29,974 同类问题）。

### 1.3 v1 → v2 权威来源表（source of truth）

| 主题 | v1 权威（冻结，只读） | v2 权威（本次重构） | v2 状态 |
|---|---|---|---|
| 候选宇宙 | `protein_layers.tsv`（109,087） | `phaded_candidate_catalog_v2.tsv` | pending-authorization（构建脚本已落地） |
| 高可信二元标签 | `merged_high_confidence_candidates.tsv`（40,690） | 唯一 `primary_disposition` + 多 `evidence_flags` | 代码已落地，真实重算未授权 |
| 模型分层 | 9 trained + 37 reference-only | `reference_query_only` / `discovery_hmm_uncalibrated` / `sequence_family_hmm_validated` / `calibrated_candidate_model` | 代码已落地 |
| 功能校准 | 六腿合成 gate（仅 hfam_70 文档级通过） | `functional_calibration_status` 独立列；hfam_70 改写为 `candidate_gate_passed_not_promoted` | 代码已落地，未 finalize |
| 亲核体 | `nucleophile_type = cys/ser/ahsmg`（三值） | `nucleophile_identity ∈ {ser,cys,unresolved}` × `motif_class`（AHSMG 归 Ser motif） | 代码已落地 |
| 定位 vs 转运 | SignalP 结果与「胞内/胞外」混用 | `transport_signal_prediction` 与 `localization_evidence` 拆列 | 代码已落地 |
| with-lipase 1,206,655 | deferred 归档层 | 「百万级待判别命中」+ `deferred_structure_review` | 层级不变 |
| 流转对账 | 散落在多文档 | `candidate_flow.tsv` + 六项具名 reconciliation check | 代码已落地 |
| 生长速率组名 | degrader / non-degrader | 候选基因携带 / 在规定搜索与质量条件下未检出候选 | 代码已落地，重算未授权 |
| 参考台账 | 723 行 × 53 列（无实验分级列） | +15 列证据质量与独立性字段 | 校验器已落地；**真实台账 curation 未授权** |

---

## 2. 实施边界与授权状态（2026-09-28）

### 2.1 本次允许且已完成的动作类型

- 新增源码、配置、测试、文档；
- 本地只读审计（`git status`、`Get-FileHash`、读取冻结 TSV）；
- 本地 `git commit`（每次在 `docs/T141_20260917_project_handoff.md` 留痕）。

### 2.2 仍然禁止、需逐项显式授权的动作（计划 Authorization Checkpoints）

| # | 动作 | 状态 |
|---|---|---|
| 1 | 把 2026-09-21 未跟踪研究快照纳入版本控制 | **未授权**（仅盘点，未 `git add`） |
| 2 | 创建任何真实分析 run 或 dated deploy | **未授权** |
| 3 | 从 723 参考集重建 MAFFT 比对或 HMM | **未授权** |
| 4 | 运行 SignalP / Foldseek / 结构预测 / HMMER against GTDB shards / gRodon | **未授权** |
| 5 | 用真实项目数据重算候选目录 | **未授权** |
| 6 | 安装或变更环境、改服务器配置、删除/归档数据、push 到 GitHub | **未授权** |

> 上述状态在本轮结束时保持 **未授权**；任何 step 若需要上表动作，一律停止并请求授权，不以「代码已就绪」替代授权。

---

## 3. Task 2–13 实施记录

> 节内逐任务记录：产出文件、验证命令与结果、偏离计划之处。**未执行的授权动作在此明确标注为 pending，不写成已完成。**

### Task 2：动态服务器资源与 run 契约

- 状态：**完成**（代码、测试、配置、驱动脚本全部落地；未在服务器上执行任何任务）。
- 产出：
  - 新增 `pipeline/scripts/server_resources.py`：`compute_thread_limit(total_cores, load_1m, reserve=10, hard_cap=40)`、`resource_record(...)`、`parse_loadavg`、`count_cores`、`read_measurements`、`measure`、CLI（读 `nproc` / `/proc/loadavg`，输出 JSON 资源记录，`--out` 写文件，容量不足时非零退出）。
  - 新增 `pipeline/tests/test_server_resources.py`（22 项）。
  - 修改 `pipeline/scripts/run_pipeline.sh`：删除 `THREADS_PREDICT=70` / `THREADS_SCREEN=70`，新增 `measure_threads()`，**在每个并行 stage 启动前即时实测**并取返回值作为启动值。
  - 修改 `pipeline/scripts/05_predict_proteins.sh`、`pipeline/scripts/06_screen.sh`：默认 `THREADS=40`（上限，非启动值），新增 `--server`；`--server` 激活时实测 `min(40, nproc − loadavg − 10)`，请求值超过实测上限即**中止**（fail-closed），不超载主机。
  - 修改 `pipeline/config/params.yaml`：`prediction.threads` 与 `screening.threads` 70 → 40，并注明「上限（ceiling），非启动值」与 `min(40, nproc - loadavg - 10)`。
- 关键语义：**40 是上限而不是启动值**；启动值必须逐次实测，任何硬编码 70 都视为治理缺陷（旧值等于假设「80 核且完全空闲」）。
- 验证：`python -m unittest pipeline.tests.test_server_resources pipeline.tests.test_run_isolation` → **22 项 OK**；`bash -n` 三个 shell 脚本全部 exit 0；`python -m compileall -q pipeline/scripts pipeline/tests` exit 0；`git diff --check` exit 0。
- `pending`：服务器上的实测线程数与资源记录需在**获授权的真实运行**时产生，本轮不执行。

### Task 3：共享证据本体（evidence ontology）

- 状态：**完成**。
- 产出：新增 `pipeline/scripts/phaded_evidence_schema.py`（受控词表 + 校验器 + CLI）、`pipeline/tests/test_phaded_evidence_schema.py`、`pipeline/config/phaded_sequence_model_gates.tsv`（`sequence-gate-v1` 预注册行）。
- 关键语义（与设计一致）：
  - `MODEL_LAYERS` 四层与 `FUNCTIONAL_CALIBRATION_STATUSES` **互相独立**，`ALLOWED_MODEL_STATES` 禁止二者矛盾（`calibrated_candidate_model` 只出现在校准层）。
  - **AHSMG 是 Ser 的 motif class，不是第三种亲核体**：`normalize_nucleophile("extracellular native-SCL/PhaZ7-like", "", "supported")` → `{"nucleophile_identity": "ser", "motif_class": "AHSMG"}`；`nucleophile_identity="ahsmg"` 被拒绝。
  - **SignalP 不是定位真相**：`normalize_localization` 只回显预测与证据两列；行级校验要求 `experimental_confirmed` 必须给出非 SignalP-only 的 `localization_evidence_source`。
  - **发现层不得有 family 判定、不得排除候选**；`excluded_input_quality` 必须有受控的输入质量理由，科学不确定性（`function_unresolved` 等）一律拒绝。
  - 缺列一律按「最弱默认值 + 不产生任何声明」处理（present → strict，absent → not validated）。
- 验证：54 项测试 OK；`compileall` exit 0；`git diff --check` exit 0；只读集成检查对冻结台账 723 行判为 `valid`（不臆造新增列的取值）。

### Task 4：723 参考台账的证据分级与独立性

- 状态：**完成（校验器与配置）**；**真实台账 curation 未执行**（计划 Step 6，需单独授权）。
- 产出：扩展 `pipeline/scripts/validate_phaded_reference_ledger.py`（238 → 878 行，保留全部原有检查与报告键）、`pipeline/tests/test_validate_phaded_reference_ledger.py`（7 → 28 项，原 7 项未改且仍通过）、新增 `pipeline/config/phaded_reference_evidence_columns.tsv`（15 列）。
- 关键语义：
  - E3/E2 必须同时给出 `study_id`、`experiment_unit_id`、`experimental_substrate_class`、`experimental_system`、直接结果与来源（DOI/PMID/PMCID）；**E1 永不计为直接功能阳性**。
  - **功能阳性按唯一 `independence_group` 计数，绝不按行计数**；同一基因/同一实验的重复报道不重复计数。
  - `discovery_training_eligible=true` 仅当序列完整性 `complete` 且家族一致；`functional_calibration_eligible` 对 grade A 必须为 false。
  - 矛盾的直接实验阴性**阻断功能校准但不从审计中移除**（该行仍被校验、仍列入审计）。
  - 默认 `compat` 模式：15 个新列缺失只报告为 `pending`，**冻结台账的校验结果逐字不变**（723 行 / 30 positive / 693 annotation_only / 38 families，exit 0）；`strict` 模式（opt-in）要求全部 15 列齐备且已 curation。
- 验证：28 项测试 OK；冻结台账 compat 模式 exit 0；strict 模式按设计 exit 1 并列出缺失列；`compileall` exit 0；`git diff --check` exit 0。
- `pending`：真实台账的 15 列仍为 `pending`（报告为 `uncurated_experimental_positive_rows = 30`）；填入属独立证据采集 run，**需授权**。

### Task 5：分离 HMM 构建与功能校准

- 状态：**完成**（未构建任何真实 HMM）。
- 产出：修改 `pipeline/scripts/build_phaded_profiles.py`（341 → 558 行）、`pipeline/tests/test_build_phaded_profiles.py`（25 项）、`pipeline/scripts/finalize_phaded_reference_panel_calibration.py`（377 → 500 行）、`pipeline/tests/test_finalize_phaded_reference_panel_calibration.py`（13 项）。
- 关键语义：
  - **实验阳性不再是建模的必要条件**：≥3 条**唯一**合格序列（`sequence_sha256` 去重，两条同哈希只算一条）即建 `discovery_hmm_uncalibrated`；1–2 条 → `reference_query_only`。
  - 构建器**只**产出 `reference_query_only` / `discovery_hmm_uncalibrated`，**永不**产出 `sequence_family_hmm_validated` / `calibrated_candidate_model`。
  - manifest 新增 `model_layer`、`functional_calibration_status`、`training_sequence_count`、`experimental_anchor_count`、`annotation_only_count`、`training_set_sha256`、`bit_reproducible`；`training_set_sha256` 是训练 FASTA 真实 SHA-256，无训练集时为 `pending`；**`bit_reproducible` 恒为 `false`**（MAFFT 不可位复现，具体比对哈希才是权威）。
  - 新增 `profile_training_audit.tsv`：每行（含被排除行）都写明 `model_layer`、`discovery_training_eligible`、`training_selected`、`exclusion_reason`、`blocks_functional_calibration`、`block_reason`。
  - 功能校准 gate 改为 `qualified_e2_e3_positive_count ≥ 3` 且 `independent_genus_count ≥ 3` 且 heldout/negative/challenge ≥ 1 且两个 status 均 `sufficient`；**通过但未提升 → `candidate_gate_passed_not_promoted`**，`calibrated_candidate_model` **只能**由显式授权动作（`--promote-calibrated-model`，默认关闭）产生。
  - 矛盾实验阴性：`blocks_functional_calibration=true`、`functional_calibration_status=blocked_contradictory_experimental_negative`，序列仍在审计中可见。
- 本轮由协调者补齐的跨任务接口（重要）：把 `blocked_contradictory_experimental_negative` 加入 `phaded_evidence_schema.FUNCTIONAL_CALIBRATION_STATUSES` 与其 `ALLOWED_MODEL_STATES`（原先只允许三个状态，会让带矛盾阴性的家族被 Task 6 判为「未识别状态」而永久无法进入序列家族层）。
- 验证：`test_build_phaded_profiles` + `test_finalize_phaded_reference_panel_calibration` = 38 项 OK；与 `test_phaded_evidence_schema`、`test_finalize_phaded_superfamily_gate` 合并 109 项 OK；`compileall` / `git diff --check` exit 0。
- **同步更新的既有测试（语义已被本任务取代，留痕说明）**：`pipeline/tests/test_finalize_phaded_superfamily_gate.py` 的 2 处合成断言——`test_family_gate_still_passes_on_the_legacy_slots` 与 `test_superfamily_gate_is_reachable_when_every_slot_is_met`——已改为断言新语义（缺 qualified 列时 fail-closed 为 `reference_only_insufficient_panel`；gate 通过但未提升为 `candidate_gate_passed_not_promoted`；授权后才为 `calibrated_candidate_model`）。**冻结文件的字节级回归测试 `test_family_gate_is_byte_identical_on_frozen_readiness` 未改且仍通过**，历史决策未被改写。
- `pending`：任何真实 HMM 重建（需授权）；下游按 `model_status == "trained"` 判断的脚本（`run_phaded_candidate_hmmsearch.py`、`build_phaded_with_lipase_deferred_tier.py`、`build_phaded_family_discovery_sets.py`）应改为读 `model_layer`——**本轮未改**，列为待办。

### Task 6：独立校验序列家族模型

- 状态：**完成**（未构建任何真实模型）。
- 产出：新增 `pipeline/scripts/evaluate_phaded_sequence_models.py`（1002 行，纯标准库）、`pipeline/tests/test_evaluate_phaded_sequence_models.py`（61 项）。
- 关键语义：
  - **序列家族层由 held-out / 混淆集表现赢得，绝不由实验阳性条数赢得**：`experimental_anchor_count ≥ 3` 但 held-out 或混淆集不合格，仍然不能升到 `sequence_family_hmm_validated`。
  - **`EFFN < 1` 是诊断字段，不是自动失败**（`eff_nseq_below_one` 不在 `blocking_reasons` 中）——本项目 EFFN 体检显示 8/9 的 profile 先验主导，把 EFFN 当闸门会误杀。
  - **held-out 少于 5 条时必须全部召回**：`heldout_rule = all_recovered_below_five`（要求 `recovered == count`）vs `min_recall_at_or_above_five`（要求 `recall ≥ 0.8`）；规则显式记录，**不靠 0.8 的取整推断**。实测 4/4 过、4/3 失败、5/4 与 10/8 过、10/7 失败。
  - **本模块永不改动 `functional_calibration_status`**：三重证明——① 归一化只接受 `ALLOWED_MODEL_STATES["sequence_family_hmm_validated"]` 允许的状态（构造上排除已校准值）；② AST 静态测试证明没有以该字段为赋值目标的语句，且 `calibrated_candidate_model` 不出现在任何可执行字符串字面量中；③ 运行时 9 行矩阵断言输出与输入逐字相同，manifest 记 `functional_calibration_status_changed=false` / `upgraded=false`。
  - 每个输出带 `sequence_claim = "sequence-defined family homology"`、`phenotype_boundary = "does not establish PHB/PHA degradation"`、以及**原样复制**的 `functional_calibration_status`。
  - 20 个稳定 snake_case 阻塞原因码（如 `discovery_family_call_present`、`heldout_recall_below_minimum`、`alignment_sha256_invalid`）；缺失值一律 `pending` 且**阻止晋升**（`bit_reproducible=false` 本身不阻止——权威是具体的 `alignment_sha256`）。
  - `reference_query_only` 输入被判为不可晋升，且**不会被评估成已验证层**。
- 验证：61 项 OK；`compileall` exit 0；`git diff --check` exit 0。
- `pending`：真实 profile 的 held-out/混淆集评分需在获授权的运行中产生。

### Task 7：v2 候选目录（字段语义修正）

- 状态：**完成（代码）**；**用真实数据重算未授权**（计划 Checkpoint 5）。
- 产出：新增 `pipeline/scripts/build_phaded_candidate_catalog_v2.py`、`pipeline/tests/test_build_phaded_candidate_catalog_v2.py`（66 项）；`pipeline/scripts/filter_phaded_high_confidence.py` **仅加弃用/边界注释**（可执行 token 流逐字不变，1941 tokens 前后相同；原有 3 项冻结 v1 过滤器测试仍通过）。
- 关键语义：
  - **AHSMG = Ser 的 motif class**（`nucleophile_identity=ser` + `motif_class=AHSMG`）；**SBD 不定义 type 1/type 2**；**SignalP 不等于实验定位**。
  - 20 列目录，**每条蛋白唯一 `primary_disposition`**；16 个受控 `evidence_flags`（`|` 连接、固定顺序、非互斥、未知值报错）。
  - `excluded_input_quality` **只能**经有据可查的输入质量理由到达；把科学不确定性措辞写进输入质量字段会**报错**（`scientific_uncertainty_is_not_an_input_error`）而不是静默降级；两者并存的对抗性行：**有据可查的输入质量理由优先**（坏输入无法分类），但不确定性仍保留在 `evidence_flags` + `scientific_uncertainty_present`，可审计、可回退。
  - **功能视图不能改动主目录处置**：源码级测试证明 `classify_row` / `primary_disposition` 从不引用 `functional_view_eligible`；该列在全部处置定稿后才计算，且返回副本；行为上前后快照相同。
  - **52 条 demotion 精确计入一次**：声明的可覆盖常量 `EXPECTED_CYS_EXPORT_DEMOTIONS = 52`（含来源依据与被覆盖标记）；重复 accession、空 accession、条数不等于声明值、以及 demotion ID 不在候选宇宙内，**全部报错**（不静默补行）。
- 验证：123 项（含 schema 与冻结 v1 过滤器）OK；`compileall` / `git diff --check` exit 0；冻结合并表实测未变（40,690 行、mtime 2026-09-20、SHA-256 `7B6F9208…`）。
- `pending`：真实候选目录的生成（需授权）；`run_pipeline.sh` 未接入该脚本（有意）。

### Task 8：候选流转对账（关闭计数差）

- 状态：**完成（对账工具）**；真实数据仅**只读**驱动验证。
- 产出：新增 `pipeline/scripts/reconcile_phaded_candidate_flow.py`、`pipeline/tests/test_reconcile_phaded_candidate_flow.py`（99 项）。
- 输出：`candidate_flow.tsv`（每条 accession 唯一处置）、`candidate_flow_summary.json`、`overlap_matrix.tsv`、逐 check 的 `check_<name>.json`、`checks_summary.json`。退出码：**2 = 流转不闭合（`unaccounted` 或 `multiply_disposed` 非零）**，1 = 某 check 不匹配，0 = 干净。`check_*` **绝不因科学不一致而抛异常**，只返回 `status=ok|mismatch` + 解释 + 计数 + 违规 accession；结构性错误（缺列、坏表、路径不存在）立即抛 `ValueError`。
- 六项具名 check 的**真实只读结论**（实测，非转录）：

| check | 真实结论 |
|---|---|
| `check_hold_residual` | **616 条确实尚未解释**：池内高可信 36,611 + hold 71,860 之外仍有 616 条不在任何桶中；换 demote 后口径 36,559 时残差为 **668 = 616 + 52**，故正确表述是「616 条既不在池内高可信也不在 hold，且这 616 不是那 52 条 demote」；池外 strong 层与这 616 的交集为 **0** |
| `check_cys_demotion` | **52 条不并入 hold 集合**：它们带 `localization_conflict_signalp_predicted_export`，但**没有一条**出现在 71,860 行 hold 表里 → `join_hold=false`、`separate_disposition=true`，52 条各列出一次 |
| `check_cys_funnel_union` | **PhaC 1,970 + iPhaZ 4,244 = 6,214 行，并集 4,301 → 重叠 1,913**（1,970+4,244−4,301）；宇宙内部分 3,913 |
| `check_grodon_manifest_delta` | **66 条差额目前仍无法分桶**：磁盘上没有 eligible-positives 列表，且 4,441 + 23 failed ≠ 4,507 → check 持续报 `mismatch` 直到 Task 12 的处置表就位（**这是有意的响亮失败**） |
| `check_motif_coverage_vs_resolvability` | 覆盖率与可判定率只有在「至少一个判据可判定」的口径下才满值；逐判据真实可判定数是 his 7,317、asp 1,290、**linker / sbd = 0** → 「覆盖但不可判定」默认视为预期，「可判定但未覆盖」永远硬失败 |
| `check_profile_inventory` | **8 超家族 + 38 家族 = 46 定义；manifest 46 = 9 trained + 37 reference-only；并强制 `profile_id == gate_profile_id`** |

- 验证：99 项 OK；`compileall` / `git diff --check` exit 0；对冻结层只读驱动端到端（输出写入仓库外临时目录）：109,087 行、`unaccounted=0`、`multiply_disposed=0`、`demoted_from_v1=52`、hold 原因非互斥计数 `{localization_conflict 48,488; architecture_criteria 19,562; common_criteria 2,538; nucleophile_conflict 1,272}`、`hold_reason_union_count=71,860`，计数锚点 `cys_demoted_count=52` / `hold_residual_unexplained_count=616` / profile 全部匹配。
- `pending`：**616 条的裁定**（需人给它们一个 `primary_disposition`）；gRodon 66 条的分桶（依赖 Task 12 的处置表）；v2 目录与 profile manifest 齐备后的端到端真实运行。

### Task 10：v1/v2 影响与规则敏感性报告

- 状态：**完成**（未读取或重算任何真实项目数据）。
- 产出：新增 `pipeline/scripts/compare_phaded_catalog_versions.py`、`pipeline/tests/test_compare_phaded_catalog_versions.py`（66 项）。
- 三个**预注册**场景（定义为数据、集中一处、不可漂移；未知名字报错并列出合法名）：
  - `current_v1_rules`：8 超家族范围 → 共同判据（interpro、架构不冲突、强 profile、非混淆子）→ **定位 hold 规则**（外泌型需 SP/LIPO/TAT/TATLIPO，非外泌型需 OTHER/pending；正是扣下 48,038/67,720 = 71% type-1 的那条）→ 逐超家族架构判据（含 x1 疏水、pf06850）。
  - `localization_uncertain_not_excluded`：仅移除 `v1_localization_export_requirement`（非外泌侧保留）。
  - `sequence_validated_models_only`：纯按 `model_layer == "sequence_family_hmm_validated"` 限制。
- **76:24 护栏**：`COMPOSITION_CAVEAT` 是固定字符串，写明「仅指在该具名规则下的候选构成、不是生物学比例」，给出判据不对称（type-1 71% 卡在分泌信号）与研究偏向，并引用接手文档 §7.3；每个构成块带 `basis="candidate_composition_under_named_rule"`、`scenario`、`rule`、`denominator`、`caveat`；写盘前 `validate_report_guard` 拒绝缺少/错配情境、caveat 缺失、块被移到情境之外，以及任何键或标签含 `biological_ratio` / `true_ratio` / `in_vivo_ratio` 等被禁片段。测试断言整份报告 JSON 文本不含被禁片段。
- **漂移护栏**：测试在只读前提下导入冻结的 `filter_phaded_high_confidence.py`，断言本模块的 `V1_SUPERFAMILY_CRITERIA` / `EXPORT` / `NO_EXPORT` / `HYDROPHOBIC` / `v1_nucleophile_type` 与其相等。
- 验证：66 项 OK；`compileall` / `git diff --check` exit 0。
- `pending`：混淆子成员资格是唯一「行外」的 v1 规则——需 `--confounders` 或行内 `is_confounder`，否则判为 `pending`（不猜）；真实场景需要候选级证据列（含 profile/interpro/架构），只喂冻结的 40,690 输出表会让多数行合法地成为 `pending`，此时 `total_proteins_status` 记为 `lower_bound_<n>_rows_pending_not_decidable`。

### Task 11：按科学优先级准备定向证据运行

- 状态：**完成（准备文档）**；五个运行**全部未执行**。
- 产出：新增 `docs/T141_20260928_phaded_evidence_run_authorization_packets.md`（701 行）：§0 授权边界声明 → §1 五包总览表（状态全部 `pending-authorization`）→ §2 通用条款（run 契约、动态线程、溯源与 `pending`、只改证据标记不删候选）→ §3 五个包 P1–P5（每个含决策问题/目标群体/抽样或全集依据/工具版本/假阳性竞争者/输出 schema/结果如何只改证据标记/线程与资源/预估规模与运行时长/明确不属于本次授权的动作）→ §4 排程与依赖 → §5 逐数字溯源表（40+ 行）→ §6 开放问题 8 项 → §7 未执行声明。
- 五个包：P1 723 参考 SignalP 补测；P2 11,610 SP-less type-1 配件域核查；P3 29,974 Cys 参考锚定残基映射（显式区分 substitution / truncation / uncertain，且**「未命中」不能当强阴性**，因 C183S 不失活）；P4 987 nPHAMCL 结构**竞争面板**（非单一 8YNV）；P5 with-lipase 1,206,655 deferred 分层代表抽样（禁删、禁入主计数）。
- **两项可复现性更正（写入文档，可追溯）**：
  1. 11,610 目标集可由**单键**精确复现：`hold_candidates.tsv` 中 `superfamily='extracellular dPHASCL type 1'` ∧ `hold_reason='localization_conflict'`（48,038）∧ `catalytic_domain_type='type1_verified'` → **11,610 行 / 11,610 唯一 accession / 8,942 基因组**（11,602 `profile_trained_hit` + 8 `profile_ambiguous_family`）。接手文档 §5.3 的「trained profile 命中 + InterPro 支持」措辞会漏掉那 8 条，且 48,038 条的 `interpro_status` 全为 `interpro_supported`（非区分条件）。
  2. `pipeline/scripts/prepare_phaded_structure_review.py` **硬编码期望 19 条**（L54–L55）且 manifest 声明 `structure_predictor` 状态为 `pending_tool`（L91，理由「无可用预测/比对可执行文件」），与实测「Foldseek 可用（自测 TM 1.00）」**直接冲突** → 987 条输入会 fail-closed 中止。已在 §6 列为授权前必须解决项，并明确「修改该脚本」不属任何包的自动授权。
- 验证：`git diff --check` exit 0；0 个 U+FFFD；16 张表格式完好。
- `pending`（不得臆造，文档内逐条列出）：P1 723 条运行时长；P2 hmmscan 时长；Pfam-A **发布版本字符串**（项目只记路径/size/SHA-256）；P2 分层所需的 GTDB taxonomy join 源与版本；P3 是否需新建 MAFFT 比对（若需则触发 Checkpoint 3，属另一授权类别）与运行时长；P4 987 条结构阶段时长与竞争面板 PDB ID 闭集；P5 抽样规模、种子与六个分层键中四个的 join 源；P2/P5 的实际 run_id 与日期前缀；各包启动时实测的线程数。

### Task 12：gRodon 标签、单位、重叠与依赖性修正

- 状态：**完成（代码与测试）**；**未创建 dated deploy、未运行 gRodon / R**。
- 产出：新增 `pipeline/scripts/grodon_reanalysis_v2.py`、`pipeline/tests/test_grodon_reanalysis_v2.py`（77 项）；修改 `pipeline/tests/test_build_grodon_manifest.py`（6 旧 + 4 新）、`pipeline/tests/test_grodon_group_stats.py`（6 旧 + 5 新）、`pipeline/scripts/plot_phaded_grodon_growth_nature.py`、`pipeline/tests/test_plot_phaded_grodon_growth_nature.py`（5 旧 + 15 新）。**旧测试无一被削弱或删除**（仍对着冻结 deploy 脚本运行并通过）。
- 关键语义：
  - **组名**：`candidate_gene_carrier`（候选基因携带）vs `candidate_not_detected_under_defined_search`（在规定搜索/质量条件下未检出候选）；`FORBIDDEN_LEGACY_LABELS` 含 `degrader` / `non-degrader`，图件模块被测试断言不再渲染它们。
  - **阳性/对照零重叠**：`build_manifest` 在三处抛含字面量 `positive/control overlap` 的 `ValueError`——(G1) 携带集 ∩ 排除集 ≠ ∅；(G2) 携带集 ∩ 声明的对照池 ≠ ∅；以及**每个返回 manifest 的后置条件**：`carrier_genomes ∩ control_genomes == ∅`。为让 (G2) 可真正触发（原调用形态在构造上不可能重叠，无法证明守卫有效），新增可选第 7 个关键字参数 `control_pool`；强制测试即原调用形态 + `control_pool={"G1"}`，并有对照测试证明同样输入在干净池下产出零重叠。
  - **单位与单次换算**：冻结 deploy 的 `08_run_grodon_one.R` 输出 `res$d` = 最小倍增时间（小时），驱动写入 `0.6931471805599453 / parts[11]`（即 ln2/d，单位 1/h）。`growth_rate_per_h(d) = ln2/d`，`d ≤ 0` 或非有限值抛 `ValueError("doubling time must be positive")`。`resolve_growth_rate`：速率字段存在则**原样返回（绝不二次换算）**并校验单位拼写；否则对倍增时间字段**恰好换算一次**；两者并存时以记录的速率为准、用 ln2/d 交叉校验，**相对偏差 > 5% 抛错并指名 `output_formula`**（可同时抓出二次换算与错单位）。
  - **属水平平均差为主估计对象**：bootstrap 与置换**按属整簇重采样**（`RESAMPLING_UNIT`），绝不以基因为单位；中位数、符号检验、逐基因组分析、胞内 vs 胞外一律标记为次要（`SECONDARY_ANALYSES`）。
  - **不使用「无差异」**：除非预注册等价边界，一律表述为「在当前设计下未检出差异」（`EQUIVALENCE_MARGIN_ABSENT` 提供有边界的等价陈述路径）。
  - **713 条复用预测审计**：`audit_reused_predictions` 要求工具版本、模型/参数、温度处理、CDS 准备、核糖体标记流程与输出公式一致；不兼容记录**排除出主分析并保留在敏感性表**（含指名不一致字段的理由）。
  - **66 条差额**：`explain_manifest_difference` 全部由 accession 级集合计算——`difference = len(eligible − manifest)`，每个桶计数是集合大小，**桶计数之和必须精确等于 `difference`**，桶之间不重叠、两个方向都不越界，未解释者按 accession 列出。锚定测试 30 + 20 + 14 + 2 = **66**；负对照证明「条数相同但成员不同」的 manifest 会报 `difference = 0` 同时暴露 28 条不合格与 28 条未解释（仅做聚合相减会把两者都掩盖）。
- 验证：118 项 OK（57.6 s）；`compileall` / `git diff --check` exit 0。
- `pending`：dated deploy 必须**同时**携带两个脚本（绘图模块经 `sys.path` 导入 `grodon_reanalysis_v2`）；合成渲染的碰撞审计报 `FIX BEFORE DELIVERY`（panel a 的 x 轴刻度标签 24 pt 间距、panel d 说明文字压到 panel e；panel 字母干净、对齐门 PASS），**需在真实渲染上复核**——标签比被替换的表型词更长属固有效应，测试不断言该结论；创建 deploy 与运行 gRodon 需授权。

### Task 9：审计并统一 HMMER 分片数据库尺度（`-Z`）

- 状态：Step 1–5 完成；Step 6 停在执行决策，**未执行任何 run / deploy / HMMER**。

**Step 1 审计结论（只读，仅依据冻结的已执行命令与日志，不依据当前源码）：**

```text
normalizable_from_frozen_scores
```

依据（逐条为冻结证据，非推断）：

| 项 | 实测 | 证据路径与命中原文 |
|---|---|---|
| 执行的 hmmsearch 命令模板 | `hmmsearch --tblout "$BUILD/$base.tbl" --domtblout "$BUILD/$base.dom" -E 1e-5 --cpu 1 "$hmm" "$filtered"` | `runs/20260917_phaded_cys_targeted_recall_01/logs/s02_provenance.raw.txt` L144 记 `deploy/20260901_formal_frozen_scan_13/scripts/formal_frozen_screen_parallel.sh` sha256 `b5e47723a783de278f97b6b49f4bdd30eae5658634d713408a21b4710bcfee53`；本地 `pipeline/scripts/formal_frozen_screen_parallel.sh` 实测 sha256 与该值**逐位相同**，其 L115 即上列模板 |
| `-Z` 是否存在 | **不存在** | 同上 L115；独立复核见 `runs/20260917_phaded_task_completion_check_01/01_recall/README.md` §8 L56：「全库扫描：score_one.sh 无 -Z（= HMMER 默认 Z=每分片序列数），-E 1e-5；100/100 片一致。」 |
| `--domZ` 是否存在 | **不存在** | 同上 L115；`runs/20260917_phaded_cys_targeted_recall_01/results/recall_summary.json` `results.pool_vs_fullscan_difference_explanation`：「全库扫描用 HMMER 默认 Z（≈每分片序列数）」 |
| 因此实际 Z | 每分片各自的序列数（**片间不同**） | 同上两处；`docs/T141_20260921_phaded_comprehensive_review.md` L29：「100 个蛋白分片（约 267 GB / 2.92 亿蛋白），10 模型 × 100 分片 = 1,000 个 hmmsearch 任务」 |
| 全库目标序列总数 | 约 `292000000`（2.92×10⁸，GTDB R232 199,923 基因组）；**近似值，精确实测总数为 pending** | `runs/20260917_phaded_cys_targeted_recall_01/results/recall_summary.json` `design.input_space_note`「267 GB 实测（du -sh）；单分片实测 292 万条蛋白」；`runs/20260921_phaded_stage_report_figures_02/results/source_data/fig1_funnel.tsv` 第 1 行 `is_approx=True`、note「single-shard measured 2.92 M × 100」；`runs/20260917_phaded_cys_targeted_recall_01/logs/s01_recon.raw.txt` L73–74 之 `du -sb` 实测 `286194724100` 字节（267 GB） |
| 分片序列数 | shard_0001 = **4,694,121** 条；shard_0051 = **2,923,820** 条（两者差 ≈1.61 倍，证明 Z 片间不同）。**完整 100 片逐片实测表 pending**（冻结 manifest 未记录任何 per-shard 计数） | 4,694,121：`runs/20260917_phaded_cys_targeted_recall_01/logs/s02_provenance.raw.txt` L200–201「=== shard_0001.faa record count ===」/「4694121」；2,923,820：`runs/20260917_phaded_task_completion_check_01/01_recall/evidence/raw_s09_s10_s11_output.txt` L18「shard_0051 records=2923820  -> expected ratio 2923820/109087 = 26.803」，其计数命令见同 run `tmp/s09_fixes.sh` L44 `grep -c '^>' "$R13/inputs/scan_shards/shard_0051.faa"` |
| 冻结 manifest 是否记录尺度 | **没有**，既无 `database_size_Z` 也无分片序列数 | `runs/20260917_phaded_cys_targeted_recall_01/logs/s02_provenance.raw.txt` L1–107 逐字复制了 run13 `scan_manifest.json`：键为 `schema_version/status/run_id/created_utc/parent_run/reused_tasks/threads/task_total/task_completed/registry_sha256/models/overlength_exclusions_sha256/hits_all_sha256`，**无任何尺度字段** |
| 原始 bitscore 是否保留 | 是（`score` 列） | `…/s02_provenance.raw.txt` L203–208 记录 `hits_all.tsv` 表头 `family shard protein tacc E-value score bias domE qname cov` 与 6,743,198 行命中；`pipeline/scripts/06b_aggregate_hits.py` L128 把 tblout 的 `p[5]`（bit score）写入 `score` 列 |
| 尺度可统一性 | **E 随 Z 严格线性**，两独立实验证实：同一目标 `-Z 2920000 → E=1.4e-08`、`-Z 109087 → E=5.1e-10`，比值 26.80（预测 26.803）／实测均值 26.9 | `runs/20260917_phaded_task_completion_check_01/01_recall/evidence/raw_s09_s10_s11_output.txt` L18–19、L28–30；`problems_and_unverifiable.md` U3「E 值随 Z 线性缩放（两独立实验证实，比值 26.8–27.5）」；`README.md` §8 L56–L61 |

判定为 `normalizable_from_frozen_scores` 的理由：原始 bitscore 与分片归属齐全（`hits_all.tsv` 的 `score`/`shard` 列），换算式 `E_full = E_shard × (Z_total / Z_shard)` 已在冻结证据上被两次独立实验证实，**不需要重跑任何 HMMER**。据此在 `pipeline/scripts/06_screen.sh`、`pipeline/scripts/hmmer_command.py`、`pipeline/scripts/06_validate_screen_manifest.py` 中把「同一全库 `-Z`」固化为入口契约，并在 `pipeline/tests/test_hmmer_database_size.py`、`pipeline/tests/test_formal_scan13_tier_run.py` 中加入判定性测试（49 项全部通过）。

**pending（不得臆造）：**① 全库 Z 的**实测总数**（现仅有 ~2.92×10⁸ 近似值与 `fig1_funnel.tsv` 的 `is_approx=True` 标注）；② **完整 100 片逐片序列数**（需在服务器对冻结 `runs/20260901_formal_frozen_scan_13/inputs/scan_shards/*.faa` 计数；该目录本地不存在）；③ `hits_all.tsv` 的 6,743,198 行未本地留存（`runs/` 内无任何 `scan_manifest.json` / `hits_all.tsv`）。

### Task 13：扩展 CI 并发布带版本的科学发布候选（release candidate）

- 状态：**完成（CI 门、当前口径语言、AGENTS.md 规则、本文档记录）**；**未 push、未打 tag、未创建 release、未创建任何 run/deploy**（本任务不执行 Step 6 的 commit，git 交由协调者）。
- 产出（严格限于本任务拥有的文件）：
  - `.github/workflows/ci.yml`：新增具名步骤 **`Run PhaDED governance and evidence tests`**，门禁 **12 个模块 / 479 项测试**（`test_run_context`、`test_run_isolation`、`test_server_resources`、`test_phaded_evidence_schema`、`test_validate_phaded_reference_ledger`、`test_build_phaded_profiles`、`test_evaluate_phaded_sequence_models`、`test_build_phaded_candidate_catalog_v2`、`test_reconcile_phaded_candidate_flow`、`test_hmmer_database_size`、`test_compare_phaded_catalog_versions`、`test_finalize_phaded_reference_panel_calibration`）。既有 4 个步骤（checkout/setup/`pip install -r requirements-ci.txt`/`compileall`/Bash 语法/空白检查）全部保留未改。
  - `README.md`：顶部加入三句绑定措辞；新增 §4「v2 证据模型（四层 `model_layer`）」；修正「高可信度」的读法（具名规则下的候选构成，非表型判定）。
  - `CHANGELOG.md`：**新增**「〇、当前口径（2026-09-28 追加）」条目（三句定稿 + 术语更正表 + 四层模型 + CI 门）；**历史章节一字未改**，只在 §六 处加一条留痕说明其「高可信度」按新表读法解释。
  - `AGENTS.md`：新增 6 条规则（见下），**全部既有规则逐条保留未改**（`test_agents_md_rules.py` 10 项断言全部仍通过）。
- **离线验证方法（可复现）**：把整个 `pipeline/` 树复制到**仓库外**临时目录（该目录**没有** `runs/`、`deploy/`、`results/`），在其中逐个模块运行；判据 = 数据无关才算通过。另在该树中跑全量 discovery 作对照。
  - 12 个门禁模块在无数据树中 **exit 0 / 479 项 OK**（`Ran 479 tests ... OK`，与在真实仓库中的 479 项结果一致）。
  - **全量 discovery 在无数据树中的实测对照**（说明「为什么 discovery 不能当发布门」）：
    `Ran 1168 tests in 61.220s` → **`FAILED (failures=14, errors=54, skipped=16)`**，exit 1；68 条 FAIL/ERROR 行全部落在 **13 个模块**上，其中 **11 个**确实依赖被忽略的冻结产物（如 `test_phaded_calibration_field_semantics`、`test_complete_phaded_motif_completion_full` 读 `runs/`；`test_build_grodon_manifest`、`test_dedupe_predictions`、`test_grodon_group_stats` 读 `deploy/`）。
  - **显式排除并有实测理由（排除 ≠ 削弱；持有冻结产物的机器上照常运行并通过）**：
    | 被排除模块 | 具体理由（实测错误，非推断） |
    |---|---|
    | `test_finalize_phaded_superfamily_gate` | 读冻结 readiness TSV `runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/profile_calibration_readiness.tsv`。`runs/` 被 `.gitignore` 忽略，`git ls-files runs/` 仅 3 个文件且都不是它 → 干净检出无此文件，4 项测试报 `ValueError: readiness input is not a regular file`。**该文件未被 git 跟踪，故按规则排除** |
    | `test_grodon_reanalysis_v2` | 读被忽略的 `deploy/` 树（`deploy/` 在 `.gitignore` 内、`git ls-files deploy/` 为 0）：`test_legacy_deploy_scripts_still_exist_and_stay_untouched` 在无 deploy 树时报 `AssertionError: False is not true`（77 项中 1 项失败）。其余 76 项为纯 fixture，但模块只有整体通过才可门禁 |
    | 其余 **116** 个被发现模块（128 = 12 门禁 + 116 未门禁） | 其中 **11 个**实测因冻结 `runs/`/`deploy/` 产物缺失而失败（见上条实测清单与「13 个模块」计数）；另外 **105 个**未门禁模块多为直接读冻结 run 输出的 audit/figure/run13 套件，本就不属本门禁范围。**没有任何测试被修改或删除**（新增/修改测试文件不在本任务拥有清单内） |
    | **2 个「非 fixture 问题」模块**：`test_public_repo_safety`、`test_finalize_phaded_subtype_reconciliation` | 二者在无数据树中失败**只因为该树不是 git 仓库 / 缺根目录 `AGENTS.md`**：`test_public_repo_safety` 的 `git ls-files -z` 报 `CalledProcessError ... exit status 128`；`test_finalize_phaded_subtype_reconciliation` 报 `FileNotFoundError: AGENTS.md`（`AGENTS.md` **是**被跟踪文件）。补齐根文件 + 仅补 `deploy/20260920_phaded_grodon_growth_01/` 后，这两个模块与 `test_grodon_reanalysis_v2` 合计 **89 项 OK**（在真实仓库中单跑这两者亦为 `Ran 3 tests ... OK`）。故它们**能**在 CI 通过，但本任务不擅自扩大门禁范围，列为**后续单独审阅的候选**。（另有第 13 个沙箱失败模块 `test_agents_md_rules`，其失败纯粹来自该沙箱的陈旧根文件副本，不属此类，在真实仓库中 10 项断言全部通过。） |
  - `test_finalize_phaded_reference_panel_calibration` **纳入门禁**：13 项全部为合成 fixture，在无数据树中通过（与其读 readiness 的同族模块相反）。
  - **discovery 步骤的决定**：`python -m unittest discover -s pipeline/tests -t . -p 'test_*.py'` 现在可用（`pipeline/__init__.py`、`pipeline/tests/__init__.py` 包标记 + `-t .` 顶层参数），本机实测 **`Ran 1243 tests ... OK (skipped=1)`**；但**故意不作为发布门**——无数据检出实测为 `Ran 1168 tests ... FAILED (failures=14, errors=54, skipped=16)`，其中上述 **11 个数据依赖模块纯因被忽略的冻结 fixture 缺失**而失败，把整个 discovery 当门会让构建为**环境**而不是代码回归失败。该决定与实测数字已写入 CI 注释。
- **三句定稿措辞（英文原文，已逐字写入 README 顶部与 CHANGELOG 新条目）**：
  ```text
  PhaDED outputs are a versioned sequence-homology candidate resource.
  Discovery and sequence-family HMMs do not prove PHB/PHA degradation.
  hfam_70 passed the project candidate gate but is not finalized, registered, or released as a function-calibrated production model.
  ```
- **旧措辞逐处更正（README/CHANGELOG 全域搜索 `degrader` / `non-degrader` / `false positive` / 作为表型断言的「高可信度」）**：实测搜索结果——
  1. **`"million false positives"` / `degrader` / `non-degrader` 在 README 与 CHANGELOG 的正文中原先 0 命中**（实测 `Contains` 检查）；`degrader` / `non-degrader` 在仓库内只出现在本文档 §1.3 与 §Task 12 记录，属 v1→v2 对照的**历史留痕**，本任务不改写本文档已有部分。
  2. **实际需更正的旧措辞 2 处**：README 原「"高可信度" = 候选证据最强的筛选层」（`高可信度` 被当作筛选层名，易读成表型断言）与 CHANGELOG §六「多证据高可信度筛选后池内 36,611 条」。
     - README：改写为「**在该具名筛选规则下**候选证据最强的层（`candidate composition under a named rule`，**不是**表型判定、也不是生物学比例）」。
     - CHANGELOG：§六 历史数字叙述**一字未改**（保留 36,611 与当时口径），改为在文首**新增**「〇、当前口径」条目，用更正映射表把旧措辞标注为「被取代的旧措辞」，并说明 §六 的「高可信度」按新表读法解释；同时给出全域搜索结论（除映射表右列外，这些旧措辞作为**主张**的出现次数为 0）。
  3. **新措辞落地**：「million-scale unresolved hits excluded from main counts pending discrimination」落在 CHANGELOG 更正映射表，并在 README §4 以「待判别的百万级命中（million-scale unresolved hits）」表述（1,206,655 条 = `deferred_structure_review`，不进主计数、不得删除）；「candidate carrier / not detected under defined search」落在 CHANGELOG 更正映射表（对应 Task 12 已落地的 `candidate_gene_carrier` / `candidate_not_detected_under_defined_search` 组名）。
- **`AGENTS.md` 新增规则（6 条，均在原有 bullet 风格内，既有规则一条未删）**：
  1. **40 是上限（ceiling），不是启动值**；启动值须紧邻启动前实测，任何硬编码线程数（含历史 `THREADS_*=70`）一律视为治理缺陷。
  2. 四个 `model_layer` 名称 + `functional_calibration_status` **独立**、不得由 `model_layer` 推导或反推；`calibrated_candidate_model` **只能**由显式授权提升（`--promote-calibrated-model`，默认关闭）产生。
  3. AHSMG 是 Ser 的 **motif class**、**不是**第三种亲核体身份；SignalP 只预测转运、**不是**定位真相；SBD **不定义** catalytic type 1/2。
  4. manifest 必须同时绑定 `training_set_sha256` / `alignment_sha256` / `hmm_sha256`；MAFFT 比对一律 `bit_reproducible=false`，**具体比对哈希**才是权威。
  5. HMMER 分片搜索必须传**同一个全库 `-Z`**，记 `database_size_Z` + 逐片序列数 + 命令；`06_screen.sh` 缺 `--database-size-z` 时 **fail-closed**。
  6. 候选流转：每条 accession 唯一 `primary_disposition`；`excluded_input_quality` **仅**用于有据可查的输入质量缺陷，**绝不**用于科学不确定性。
  - 既有「允许用 DED 家族序列（含 `annotation_only`）训练发现层」规则**保留**，措辞对齐为 `model_layer=discovery_hmm_uncalibrated`。
- **本任务未新增/未修改任何测试文件**（不在拥有清单内），因此不存在「削弱或删除测试」的情形；被排除的模块只是不进 CI 门，测试本体保持原样可运行。
- `pending`：① 协调者执行 Step 6 的 `git add` / `git commit`（本任务按治理规则不做）；② **push / tag / GitHub Release / Zenodo DOI 未授权**，发布候选的最终对外发布待操作者授权；③ Task 12 的 dated deploy 与真实 gRodon 运行仍未授权；④ 本文档 §Task 8 的 616 条处置、gRodon 66 条分桶仍 pending。
- 最终验收核对（本任务范围内）：CI 门绿 → 见 §4；三句措辞在本文件、README、CHANGELOG 一致；`runs/`、`results/`、`deploy/` 三个目录**零改动**（见 §4 的 `git status` 结果）。

---

## 4. 累计验证基线（本地可离线执行）

```powershell
python -m unittest discover -s pipeline/tests -t . -p 'test_*.py' -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
```

> `-t .`（top-level directory）是 `discover -s pipeline/tests` 的必需参数：包标记
> `pipeline/__init__.py` 与 `pipeline/tests/__init__.py` 落地后，模块名以 `pipeline.tests.*`
> 导入，缺少 `-t .` 会因顶层目录推断而导入失败。

**2026-09-28 Task 13 实测结果（真实输出，非转录）：**

| 检查 | 命令 | 实测结果 |
|---|---|---|
| 全量 discovery（本机，数据齐备） | `python -m unittest discover -s pipeline/tests -t . -p 'test_*.py' -v` | **`Ran 1243 tests in 73.498s` → `OK (skipped=1)`，exit 0** |
| 编译 | `python -m compileall -q pipeline/scripts pipeline/tests` | exit 0 |
| 空白/EOF | `git diff --check` | exit 0 |
| CI 门（数据无关子集） | 12 个模块（见 Task 13 记录） | **`Ran 479 tests ... OK`**（在无 `runs/`/`deploy/` 的临时树中实测同样 479 项 OK） |
| 冻结目录 | `git status --short` / `git diff --stat` | `runs/`、`results/`、`deploy/` **零条目**（只读边界保持） |

其中唯一 skip 为**环境相关**：`test_run_isolation` 的 `test_path_symlink_ancestor_fails_closed`
在 Windows 无符号链接特权时报
`skipped "symlink creation unavailable: [WinError 1314] …"`（模块自身 fail-closed 判定，非削弱）；
在允许创建符号链接的 Linux CI 上该项会实际执行。

**数据/工具依赖测试**：显式跳过并写明理由，**不得删除或削弱**。
（CI 的具体排除清单与逐条实测理由见 §3 Task 13；被排除者≠被削弱，它们在持有冻结
`runs/` / `deploy/` / `results/` 产物的机器上照常运行并通过。）

---

*本文档为 candidate-only 实施记录。所有 profile/HMM/domain/motif/SignalP/结构/系统发育/分类与生长速率证据只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
