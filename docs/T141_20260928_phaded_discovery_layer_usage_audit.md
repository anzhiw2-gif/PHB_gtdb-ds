# PhaDED 发现层 / with-lipase deferred 层下游误用审计（F6）

> **任务**：`docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign-followup.md` §2 **F6**（A 类，无需计算授权）
> **性质**：只读静态审计。本文件只报告，不修复任何脚本，不删除任何 deferred 层文件，不写 `runs/`、`results/`、`deploy/`。
> **审计器**：`pipeline/scripts/audit_discovery_layer_usage.py`（标准库 + `ast`）
> **回归测试**：`pipeline/tests/test_audit_discovery_layer_usage.py`（34 项，全过）
> **快照**：工作树 HEAD `5bb5943bd67cdd83196180727dfb0f775a8b6e69`，审计时刻 **2026-09-28 22:05（本地）**。
> **candidate-only 声明**：本审计只判定「模型层被如何使用」，不产生任何候选、家族或表型结论。

---

## 0. 结论速览

| 量 | 值 |
|---|---:|
| 扫描脚本数 | **192**（`*.py` 172 + `*.sh` 20；`audit_discovery_layer_usage.py` 自身按声明排除） |
| 审计表行数（findings） | **184** |
| `violation=true` | **2** |
| `violation=unreviewed` | **105** |
| `violation=false`（合规） | **77** |
| 有发现（任一类别）的脚本 | 39 |
| 完全未提及发现层/deferred 层的脚本 | 153 |
| 退出码 | **1**（有违规项，设计如此） |
| `pipeline/config/formal_scan_models.tsv` 是否列出发现层 | **否（absent）** |
| `pipeline/config/formal_scan_models.tsv` SHA-256 | `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3` |

> **重跑漂移（可复核）**：22:12 重跑时并发的其他任务新建了 `pipeline/scripts/build_phaded_v2_catalog_input.py`，于是 192→**193** 脚本、findings 184→**185**、unreviewed 105→**106**（多出的唯一一行是 `build_phaded_v2_catalog_input.py:805 discovery_decision/unreviewed`）。用集合比对确认：除该新增文件外，两次快照**逐行完全相同**，违规项仍为同样 2 处。行号会随并发落地漂移，修复时请以「路径＋问题描述」重新定位。

`use_class` 分布：`recall_only` 51、`unknown` 91、`thresholding` 8、`family_call` 8、`registry_write` 16、`deferred_layer_touch` 9、`counting` 1。

**一句话结论**：发现的 2 处违规都是 **R2（发现层不得产生 family 判定）**，位于 `pipeline/scripts/attribute_phaded_candidate_flow_residual.py:1139` 与 `:1141`（发现层的 `superfamily_claim` 被写进 `major_superfamily` 字段）。其余 105 项是**审计无法在静态源码上定论**的调用点，已逐条列出「需要什么才能判定」。R3（注册表）、R4（deferred 禁删）**未发现违规**：注册表 16 处引用全为只读，deferred 层 9 处引用中 1 处为新建写入（属"保留"要求）、0 处删除。

---

## 1. 审计口径（可复核）

### 1.1 被审计的治理规则

| 规则 | 内容（AGENTS.md） |
|---|---|
| `R1_recall_only_only` | 发现层（`model_layer="discovery_hmm_uncalibrated"`）只用于**召回**；不得用于筛选、删除或降级任何候选 |
| `R2_no_family_call` | 发现层不得产生 family 判定；每份输出必须带其层标注 |
| `R3_registry_exclusion` | 发现层不得写入 `pipeline/config/formal_scan_models.tsv` |
| `R4_deferred_preserved` | with-lipase（`DED_hfam_2`，1,206,655 条池外命中）必须保留为 deferred 结构验证层，**禁删** |
| `R5_deferred_out_of_counts` | deferred 层不进入任何候选 / 高可信度计数 |
| `R6_reference_query_only_unscoreable` | `reference_query_only` 无 HMM，不得被当作有 HMM 打分 |
| `R7_layer_status_independent` | `model_layer` 与 `functional_calibration_status` 独立，不得互相推导 |
| `NONE` | 该行不引用任何规则（仅记录，用于"无法判定是否触犯规则"的构造） |

### 1.2 `use_class`（用途分类）与 `violation` 语义

`use_class ∈ {recall_only, thresholding, family_call, counting, registry_write, deferred_layer_touch, unknown}`：

* `recall_only`：只打分 / 排序 / 标注 / 路由，不排除任何候选；
* `thresholding`：用分数决定保留或剔除（无验证闸门）；
* `family_call`：把层/分数写成 family / superfamily / subtype / confidence 判定；
* `counting`：治理层文本进入计数；
* `registry_write`：引用 `formal_scan_models.tsv`（再区分只读与写入）；
* `deferred_layer_touch`：读/写 with-lipase deferred 归档（`note` 里必写 `read_only=true|false`）；
* `unknown`：静态源码无法定论。

`violation ∈ {true, false, unreviewed}`，**没有第四种取值**：无法判定的构造一律 `unreviewed`，绝不静默写 `false`。规则 7 的误用形状与规则 2 同类（用层**制造判定**），因此按任务给定的封闭词表归入 `family_call`，精确构造写在 `signal` 与 `note` 中。

### 1.3 方法（为什么不是 grep）

* 每个 `*.py` 用 `ast` 解析；**注释与 docstring 对 AST 不可见**，所以"注释里写着 `if discovery_score > 0.5: drop(...)`"不会被判成违规。文本扫描（`scan_text_signals`）会把这类行记录为 `*_in_prose` / `prose_only_governance_claim` 的 `unreviewed` 行——**刻意记录、绝不判违规**，测试里用同一片段同时证明「文本扫描会命中」与「AST 不会误判」。
* 对 `formal_scan_models.tsv`、with-lipase deferred 归档这类**只以字面量出现的配置/路径引用**，使用文档化的文本扫描，并按「路径绑定名是否被写/删」判定 `read_only`。
* `*.sh` 无 AST，只用文档化文本启发式；除字面量路径引用以外一律 `unreviewed`。
* 污点分析有四种来源标记：`layer`（标识发现层）、`score`（发现层分数/E 值/命中/家族字段）、`family`、`status`（功能校准状态）。**污点按函数作用域传播**，不做文件级名字传播。
* 判定按被守卫块的**效果**给出：`drop`/`select`（保留剔除）、`assign`（写判定）、`count`（计数）、`rank`/`label`（仅打分标注）、`collect`（收集）、`emit_then_skip`（先写行再跳过＝标注而非过滤）、`suppress`/`weaken`（写空值/最弱值）、`raise`（fail-closed 执法）、以及无法判定。
* 审计器自身被排除（`--exclude` 默认值），因为它的 token 表与正则本身就包含这些字面量；该排除是**声明式**的，会写进 summary 的 `excluded` 字段。

### 1.4 已内置的假阳性抑制（本次实做中逐个验证过）

1. **注释/docstring 不算代码**（AST 不可见）；
2. `ePhaZ_broad_discovery` 是 ePhaZ 模型名，**不是** PhaDED 发现层——不匹配；
3. `cov` 曾误命中 `dis*covery*`、`core` 曾误命中 `select_s*core*able`：所有短提示词改为词界匹配；
4. **容器字面量不传播元素污点**（否则一张含发现层成员的表会污染该表的所有后续使用）；
5. 全字符串常量的 dict 视为**改名表/字段表**，不是数据写入；
6. 模块级 `SCREAMING_CASE` 常量只是**命名**一条规则（`REASON_DISCOVERY_FAMILY_CALL = "..."`），不是 family 判定；
7. 含空格的散文字符串不算发现层 provenance；
8. 路径的读/写判定只看**该行绑定的名字**（否则文件里任一写入名会把每一句涉及 deferred 的注释都判成写入）；
9. `model_status` 的关键词判定只在**包含式比较**（`==` / `in`）且产生打分/候选/家族后果时才算违规，纯报表计数不算。

---

## 2. 全脚本用途分类总表

> 39 个脚本有发现；其余 **153** 个被扫描脚本（受审计器自身排除规则影响的 1 个除外）在整个文件中**没有任何**发现层 / deferred 层引用，故不出现在表中。

| 脚本 | 行数 | V=true | U=unrev | C=false | use_class（行数） |
|---|---:|---:|---:|---:|---|
| `acquire_phaded_reference_panel_amendment.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `annotate_archaea_domains.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `annotate_phaded_inphascl_cys_evidence.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `assign_phaded_candidates.py` | 1 | 0 | 0 | 1 | recall_only(1) |
| `attribute_phaded_candidate_flow_residual.py` | 20 | **2** | 18 | 0 | family_call(4) unknown(16) |
| `audit_phaded_profile_gap.py` | 4 | 0 | 4 | 0 | thresholding(1) unknown(3) |
| `audit_phaded_profiles.py` | 2 | 0 | 2 | 0 | thresholding(1) unknown(1) |
| `build_archaea_manual_review.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `build_phaded_candidate_catalog_v2.py` | 20 | 0 | 10 | 10 | family_call(1) recall_only(10) unknown(9) |
| `build_phaded_cys_discovery_training_set.py` | 2 | 0 | 0 | 2 | recall_only(2) |
| `build_phaded_family_discovery_sets.py` | 5 | 0 | 3 | 2 | recall_only(2) thresholding(1) unknown(2) |
| `build_phaded_pool_external_profile_layer.py` | 4 | 0 | 1 | 3 | deferred_layer_touch(3) family_call(1) |
| `build_phaded_profiles.py` | 11 | 0 | 5 | 6 | recall_only(6) unknown(5) |
| `build_phaded_subtype_matrix.py` | 2 | 0 | 2 | 0 | unknown(2) |
| `build_phaded_with_lipase_deferred_tier.py` | 16 | 0 | 8 | 8 | deferred_layer_touch(3) family_call(2) recall_only(5) unknown(6) |
| `calibrate_ephaz_layers.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `compare_phaded_catalog_versions.py` | 4 | 0 | 4 | 0 | unknown(4) |
| `curate_phaded_reference_evidence.py` | 5 | 0 | 3 | 2 | deferred_layer_touch(1) recall_only(1) unknown(3) |
| `design_phaded_discrimination_panel.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `evaluate_phaded_sequence_models.py` | 5 | 0 | 2 | 3 | recall_only(3) unknown(2) |
| `finalize_ephaz_calibration.py` | 1 | 0 | 0 | 1 | registry_write(1) |
| `finalize_phaded_reference_panel_calibration.py` | 4 | 0 | 4 | 0 | unknown(4) |
| `formal_frozen_screen.sh` | 4 | 0 | 0 | 4 | registry_write(4) |
| `formal_frozen_screen_parallel.sh` | 4 | 0 | 0 | 4 | registry_write(4) |
| `formal_scan13_tier_processing.sh` | 3 | 0 | 0 | 3 | registry_write(3) |
| `merge_phaded_family_discovery_scores.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `monitor_formal_scan.sh` | 1 | 0 | 0 | 1 | registry_write(1) |
| `overall_pha_phb_archaea_audit.py` | 3 | 0 | 3 | 0 | unknown(3) |
| `parse_phaded_cys_discovery_hmm_scores.py` | 14 | 0 | 8 | 6 | recall_only(6) thresholding(3) unknown(5) |
| `parse_phaded_cys_targeted_recall.py` | 10 | 0 | 6 | 4 | recall_only(4) thresholding(1) unknown(5) |
| `phaded_evidence_schema.py` | 8 | 0 | 4 | 4 | recall_only(4) unknown(4) |
| `prepare_phaded_mapping_inputs.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `prepare_phaded_reference_panel.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `preregister_phaded_cys_discovery_gates.py` | 5 | 0 | 0 | 5 | recall_only(2) registry_write(3) |
| `reconcile_phaded_candidate_flow.py` | 3 | 0 | 0 | 3 | counting(1) deferred_layer_touch(2) |
| `reconcile_phaded_reference_panel.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `run_phaded_candidate_hmmsearch.py` | 8 | 0 | 4 | 4 | recall_only(4) thresholding(1) unknown(3) |
| `split_ephaz_seeds.py` | 1 | 0 | 1 | 0 | unknown(1) |
| `validate_phaded_reference_ledger.py` | 4 | 0 | 3 | 1 | recall_only(1) unknown(3) |

---

## 3. `violation=true`（2 项）

两条都属同一处构造（`discovery_family_assignment`＝判定点、`discovery_family_field_emitted`＝写出点），**规则 R2**。

| # | 位置 | signal | 证据（源码） | 为什么算违规 | 什么会推翻它 |
|---|---|---|---|---|---|
| 1 | `pipeline/scripts/attribute_phaded_candidate_flow_residual.py:1139` | `discovery_family_assignment` | `claim = _value(row, "discovery_superfamily_claim")` → `if claim:` 进入返回分支 | 发现层的 `superfamily_claim` 决定 `major_superfamily` 取值（发现层"产生"了家族/超家族判定） | 若调用方只把该值当"带来源标注的归因证据"而不作为家族判定使用（需人读调用方与下游目录） |
| 2 | `pipeline/scripts/attribute_phaded_candidate_flow_residual.py:1141` | `discovery_family_field_emitted` | `{"major_superfamily": claim, "major_superfamily_source": "superfamily_classification::superfamily_claim"}` | 该行把一个 **family/superfamily 判定字段**填成发现层值；行内只有"来源列"，**没有** `model_layer=discovery_hmm_uncalibrated` 标注 | 若补上显式 `model_layer` 层标注（如 `build_phaded_candidate_catalog_v2.resolve_family_call` 的做法：发现层行直接输出空 family call），则该构造会降为 `unreviewed`（与 with-lipase 归档同级处理） |

**对照（合规写法）**：`pipeline/scripts/build_phaded_candidate_catalog_v2.py:791-793` 对发现层行返回空 family call；`pipeline/scripts/phaded_evidence_schema.py:686-698` 对"发现层行携带 family call"直接 fail-closed 抛错。这两处被审计判为 `recall_only/false`（R2 执法）。

---

## 4. `violation=unreviewed`（105 项）——每项需要什么才能判定

### 4.1 发现层/deferred 值驱动的决策，效果不可判（62 项）

| 子类 | 条数 | 位置（文件:行） | 需要什么才能判定 |
|---|---:|---|---|
| `discovery_decision` + `NONE`：只比较 `model_layer`/校准字段（**没有点名发现层**），是否走到发现层分支取决于运行时 manifest | 35 | `annotate_archaea_domains.py:42`；`attribute_phaded_candidate_flow_residual.py:972,973,1069,1078`；`build_archaea_manual_review.py:22`；`build_phaded_candidate_catalog_v2.py:745,751,859,863,933,942,1018`；`build_phaded_family_discovery_sets.py:100`；`build_phaded_with_lipase_deferred_tier.py:115,188`；`calibrate_ephaz_layers.py:236`；`curate_phaded_reference_evidence.py:1418,1423`；`evaluate_phaded_sequence_models.py:604,606`；`finalize_phaded_reference_panel_calibration.py:258,264,268,270`；`overall_pha_phb_archaea_audit.py:72,85,87`；`phaded_evidence_schema.py:454,458`；`run_phaded_candidate_hmmsearch.py:126,195`；`split_ephaz_seeds.py:295`；`validate_phaded_reference_ledger.py:623,629` | 冻结 `profile_manifest.tsv` 的各分支取值直方图（哪些 profile 真为 `discovery_hmm_uncalibrated`）＋ 该分支输出表的消费方 |
| `discovery_decision` + `R1`：发现层 provenance 的值驱动了决策，但**没有与阈值比较**（真值判断/成员判断），既可能是标注也可能是筛选 | 27 | `annotate_phaded_inphascl_cys_evidence.py:895`；`attribute_phaded_candidate_flow_residual.py:629,760,789,791,802,805,828,834,862,878,1071,1675`；`build_phaded_profiles.py:307,348,406`；`build_phaded_with_lipase_deferred_tier.py:292,500`；`compare_phaded_catalog_versions.py:571,580`；`parse_phaded_cys_discovery_hmm_scores.py:159,195,201,537`；`parse_phaded_cys_targeted_recall.py:147`；`phaded_evidence_schema.py:676`；`validate_phaded_reference_ledger.py:632` | 该块写出/追加的目标集合**下游是否是候选集**（追到消费者）；只有追到"候选/高可信度"才升级为 R1 违规 |

### 4.2 用发现层分数做保留/剔除，但集合不是"候选集"（6 项，`thresholding`）

| 位置 | 形状 | 需要什么才能判定 |
|---|---|---|
| `build_phaded_family_discovery_sets.py:150` | `if row["model_layer"] not in DISCOVERY_SET_MODEL_LAYERS: continue` | 确认该表只用于"发现集/召回集"（层标签路由），不含候选流转语义 |
| `parse_phaded_cys_discovery_hmm_scores.py:257` | `_hits_in(...)`：`entry["best_e"] < DISCOVERY_EVALUE_THRESHOLD` → 收 hit | G1/G2 闸门的分母/命中集是否被当作候选集使用（当前看是召回/闸门统计） |
| `parse_phaded_cys_discovery_hmm_scores.py:311` | 汇总 `called = [...]`（发现层阈值内的训练 accession） | 同上（G2 灵敏度统计） |
| `parse_phaded_cys_discovery_hmm_scores.py:454` | 同一阈值比较的重复出现（另一处调用） | 同上 |
| `parse_phaded_cys_targeted_recall.py:201` | `if record["best_e"] >= DISCOVERY_EVALUE_THRESHOLD: continue` | 该脚本自述为 targeted **recall**；需确认输出表不充当候选集 |
| `run_phaded_candidate_hmmsearch.py:197` | `elif layer == DISCOVERY_LAYER:` 路由到独立 recall-only 表 | F3 已写明"发现层绝不进判别表"；需确认该表不回流候选 |

### 4.3 发现层 family/confidence 字段（4 项 + 1 项本地绑定）

| 位置 | 形状 | 需要什么才能判定 |
|---|---|---|
| `attribute_phaded_candidate_flow_residual.py:1056` | 字段 `"discovery_superfamily_confidence"` 明确指出是发现层置信度 | 判定该字段是否被下游当"家族判定"消费（字段名已自带 provenance） |
| `attribute_phaded_candidate_flow_residual.py:817` | `confidence = _value(row, "discovery_superfamily_confidence")` | 该局部是否流向 family call（见 §3 同一函数） |
| `build_phaded_pool_external_profile_layer.py:119` | 输出行含 `discovery_best_family` 字段 | 是否只作 provenance 记录（该脚本自述沿用 deferred-tier 口径） |
| `build_phaded_with_lipase_deferred_tier.py:311,351` | deferred 归档行内 `superfamily` 由发现层 `best_family` 或 `discovery_unique` 得出，行内带 `model_layer` 标注 | 由协调者裁决："带层标注的归档级超家族字段"是否仍算 R2 的 family 判定（本审计判 `unreviewed`，不判违规） |

### 4.4 `reference_query_only` 被包含式选中（8 项，R6）

`acquire_phaded_reference_panel_amendment.py:152`；`audit_phaded_profile_gap.py:160,205`；`audit_phaded_profiles.py:27`；`build_phaded_subtype_matrix.py:226`；`design_phaded_discrimination_panel.py:1128`；`prepare_phaded_reference_panel.py:172`；`reconcile_phaded_reference_panel.py:110`。

**需要什么才能判定**：这些选择点在**本函数内没有打分调用**（审计已核），但需跟随返回值确认调用方是否把 `reference_query_only` 送进 `hmmsearch`/tblout 流程（那将是 R6 违规：为没有 HMM 的 profile 编造分数）。按仓库自身文档，`reference_query_only` 只作"参考词汇表"保留，故预期为合规。

### 4.5 遗留 `model_status` 判据（5 项，F3 同类缺陷，但**不在** F3 的三个文件内）

`audit_phaded_profile_gap.py:101,159`；`audit_phaded_profiles.py:26`；`build_phaded_subtype_matrix.py:62`；`prepare_phaded_mapping_inputs.py:87`。

**为什么不是违规**：审计把这些位置判为"报表分区/标注"（集合不是候选/高可信度集，见 §5.2 的反例判据），因此停在 `unreviewed`。但判据仍是遗留列——`build_phaded_profiles.py:407-408` 会同时写 `model_status="trained"` 与 `model_layer="discovery_hmm_uncalibrated"`，所以任何以 `model_status == "trained"` 为判据的地方**在运行时无法区分发现层**。
**需要什么才能判定**：读冻结 `profile_manifest.tsv`，确认这些位置的分支是否真的会接收发现层 profile；若是，`build_phaded_subtype_matrix.py:62/222/226`（把发现层判成 `high_candidate_only`/`moderate_candidate_only` 置信度）风险最高。

### 4.6 `model_layer` 反向推导（1 项，R7）

`build_phaded_candidate_catalog_v2.py:555`：`status == "calibrated_candidate_model"` → `return "calibrated_candidate_model"`（**由状态推层**）。
**需要什么才能判定**：裁决"把行内已声明的校准状态解析为其层"是否属 R7 禁止的"反推"。注意同文件的 `resolve_functional_calibration_status`（`:569-581`）只读声明值，是合规写法；`build_phaded_profiles.py:402-405` 的状态由 `blocking_negatives` 独立得出，也合规。

### 4.7 仅出现在散文里的治理断言（13 项 prose claim + 5 项散文分数比较）

* `prose_only_governance_claim`（docstring/注释断言"不筛选/不删除/不产生 family call"）：`build_phaded_candidate_catalog_v2.py:152,780`；`build_phaded_family_discovery_sets.py:22`；`build_phaded_profiles.py:13,199`；`build_phaded_with_lipase_deferred_tier.py:16,17`；`compare_phaded_catalog_versions.py:78`；`curate_phaded_reference_evidence.py:150`；`merge_phaded_family_discovery_scores.py:2`；`parse_phaded_cys_targeted_recall.py:14`；`phaded_evidence_schema.py:619`；`run_phaded_candidate_hmmsearch.py:16`。
* `discovery_score_comparison_in_prose`（注释/散文字符串里出现发现层比较，**不是可执行代码**）：`compare_phaded_catalog_versions.py:44`；`parse_phaded_cys_discovery_hmm_scores.py:15`；`parse_phaded_cys_targeted_recall.py:11,149,197`。

**需要什么才能判定**：散文断言不可静态验证。审计**刻意记录**而不判违规（这与"注释里写了过滤"不得报违规的要求一致）。要闭环需人工逐条核对断言对应的代码路径，或在 CI 中改为可执行断言。

---

## 5. 合规证据（正面确认，77 项）

### 5.1 R3：注册表只被读，未被写（16 项）

`formal_frozen_screen.sh:11,30,103,170`；`formal_frozen_screen_parallel.sh:14,45,146,150`；`formal_scan13_tier_processing.sh:95,116,142`；`monitor_formal_scan.sh:23`；`finalize_ephaz_calibration.py:97`；`preregister_phaded_cys_discovery_gates.py:20,110,237`。全部 `read_only=true`（哈希校验、运行前复制进 `inputs/`、或声明"未修改"）。**0 处写入 config 本体**。
特别地，`preregister_phaded_cys_discovery_gates.py` 的 G3 闸门就是「运行前后 `formal_scan_models.tsv` SHA-256 不变」，与本次审计的注册表核对互相印证。

### 5.2 R1/R2：标注与边界执法（38 项）

* `discovery_layer_label_emitted` 21 项：输出行显式带 `model_layer=discovery_hmm_uncalibrated`（`build_phaded_profiles.py`、`parse_phaded_cys_*`、`build_phaded_with_lipase_deferred_tier.py` 等）——**强制 provenance 标注**，合规。
* `discovery_layer_gate_compliant` 17 项 + `layer_gate_on_model_layer` 9 项：发现层行被写空值/最弱值，或 fail-closed 抛错（`phaded_evidence_schema.py:686-698`；`build_phaded_candidate_catalog_v2.py:753,791`）。
* `reference_query_only` 标注/拒绝而非打分 4 项（含 `evaluate_phaded_sequence_models.py:608` 追加 `REASON_MODEL_LAYER_REFERENCE_ONLY`）。

### 5.3 R4/R5：deferred 层被保留，计数未进候选（10 项）

`build_phaded_pool_external_profile_layer.py:18,61,62`（只读引用）；`build_phaded_with_lipase_deferred_tier.py:42,57`（标签常量）、`:367`（**新建** deferred TSV，`read_only=false` —— 属 AGENTS.md 要求的"保留"动作）；`curate_phaded_reference_evidence.py:247`；`reconcile_phaded_candidate_flow.py:107,118`（只读）。**0 处删除**（`unlink`/`rmtree`/`rm` 均未出现在 deferred 归档上）。
`reconcile_phaded_candidate_flow.py:120` 的 `REVIEW_WITH_LIPASE_DEFERRED = 1206655` 被判为 `counting/false`：该常量记录的是 deferred 层**自身规模**，注释明确要求它"留在每个候选计数之外"；审计只能确认它本身不是候选计数（该常量当前未被任何计数消费——这点由静态分析确认）。

### 5.4 R6/R7：未发现违规

* 无任何脚本把 `reference_query_only` 送进打分调用（审计在 8 个包含式选择点所在函数内均未发现打分调用）；
* 无任何脚本把 `functional_calibration_status` **赋值**为层推导结果（R7 正向无违规）。

---

## 6. `pipeline/config/formal_scan_models.tsv` 缺失确认

审计只读取该文件（未修改；`git status` 中该文件无改动）：

* 路径：`pipeline/config/formal_scan_models.tsv`
* 大小：450 字节；列：`model, hmm_source, threshold, min_cov, report_group`；行数：10 个模型
* **SHA-256：`8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`**
* 模型清单：`ePhaZ_curated_core, ePhaZ_broad_discovery, iPhaZ, OH, BdhA, ArchPhaZ_patatin, ArchPhaZ_hydrolase, PhaJ, phasin, PhaC`
* **发现层是否存在：否（`discovery_layer_absent = true`）**。匹配的是层标识符（`discovery_hmm_uncalibrated` / `cys_discovery_uncalibrated` / `uncalibrated`），**不是**裸词 `discovery`——`ePhaZ_broad_discovery` 是 ePhaZ 家族模型，与 PhaDED 发现层无关，被显式排除（回归测试对此有专门断言）。

---

## 7. 本审计自身的限度（静态分析看不到什么）

1. **不执行**任何被扫描脚本；不 import、不 trace、不跑 HMMER。
2. **动态派发不可见**：`getattr`/`setattr`、`importlib`、拼字符串的模块名、字典驱动的可调用对象、插件注册表、运行时拼装的 `subprocess` argv——只经动态派发到达的误用会被报成 `recall_only` 或 `unknown`，而不是违规。
3. **数据相关分支不可见**：某张表里是否是发现层行、某个 manifest 行是否发现层、真实数据走哪个分支，都是运行时事实。provenance 只经数据到达（例如值来自路径不含层名的表）的构造一律 `unknown`/`unreviewed`。
4. **`*.sh` 只有文本启发式**：除字面量路径引用外，shell 行一律 `unreviewed`；shell 内嵌 heredoc 里用 Python 写注册表的情形可能被报成只读。
5. **中间层名字**：`disc_e`（`^disc_`）可识别，但 `best_evalue` 这类**从发现层表读出却没有层名**的字段不可识别——这是 §4.1 那 27 项存在的原因。
6. **无违规行 ≠ 合规证明**：只能说明审计在该构造上看不到误用。
7. **快照性**：本次审计期间有并发任务在改脚本（见 §8.5），行号与结论对应 2026-09-28 22:05 的工作树状态；重跑请以审计表 `script:line` 为准。
8. **散文断言不可验证**：docstring/注释里写的"不筛选/不删除"只是断言，审计只记录（13 项）。

---

## 8. 修复清单（Repair list）——本任务不修复，只点名

> **本审计不修改任何现有管线脚本**（治理要求：并发任务正在编辑相关文件，两个 agent 同时改同一文件会产生冲突）。以下按"是否已被并发任务覆盖"分组。

### 8.1 违规项（必须处理：R2）

| 位置 | 问题 | 建议动作 | 并发归属 |
|---|---|---|---|
| `pipeline/scripts/attribute_phaded_candidate_flow_residual.py:1139-1143` | 发现层 `superfamily_claim` 写入 `major_superfamily` | 二选一：① 发现层来源时输出空值（`build_phaded_candidate_catalog_v2.resolve_family_call` 的做法）；② 保留值但给该行补显式 `model_layer=discovery_hmm_uncalibrated` 标注，并声明它是**召回字段**而非 family call | **该文件为本次审计期间新建的未跟踪文件**（mtime 2026-09-28 21:51），**不在** F3 三个文件内。归属需协调者确认（疑为 F1 的归因脚本）；**不要与 F1 并行修改** |

### 8.2 F3 已覆盖（不要重复改）

| 位置 | 问题 | 状态 |
|---|---|---|
| `pipeline/scripts/run_phaded_candidate_hmmsearch.py` | 原 `model_status == "trained"` 判据（原 L35） | ✅ 已被 F3 改写：现读 `model_layer`，发现层只进 `phaded_profile_scores_recall_only.tsv`，`reference_query_only` 不可打分（审计判合规：`:197` recall-only 路由、`:250` 行内 `family_call=""`） |
| `pipeline/scripts/build_phaded_with_lipase_deferred_tier.py` | 原 `model_status` 判据（原 L18/L76） | 🔄 F3 编辑中（工作树已修改）：`read_trained_model_superfamily` 等判据在改；审计对其判 `unreviewed`（`:311,:351` 归档级超家族字段） |
| `pipeline/scripts/build_phaded_family_discovery_sets.py` | 原 `model_status == "reference_only"` 判据（原 L63） | 🔄 F3 编辑中；审计对其 `:150` 层标签路由判 `unreviewed` |

### 8.3 需要后续任务处理（F3 未覆盖）

| 优先级 | 位置 | 问题 | 需要的动作 |
|---|---|---|---|
| 高 | `pipeline/scripts/build_phaded_subtype_matrix.py:62,222,224,226` | 以 `model_status == "trained"` 决定 `subtype_confidence=high_candidate_only/moderate_candidate_only`——发现层在遗留列下也是 `trained` | 改用 `model_layer`；发现层/`reference_query_only` 最多 `low_candidate_only` |
| 中 | `pipeline/scripts/prepare_phaded_mapping_inputs.py:87` | 以 `model_status == "trained"` 选 profile | 改用 `model_layer` |
| 中 | `pipeline/scripts/audit_phaded_profiles.py:26`、`audit_phaded_profile_gap.py:101,159` | 报表脚本按遗留列分区 | 改读 `model_layer`（低风险，但同类判据） |
| 中 | `pipeline/scripts/build_phaded_candidate_catalog_v2.py:555` | 由 `functional_calibration_status` 反推 `model_layer`（R7） | 裁决是否合法；若否，改为要求显式 `model_layer` |
| 中 | `attribute_phaded_candidate_flow_residual.py:802,828,834,1056,1071` 等 | 发现层 provenance 值驱动的决策/字段（§4.1–4.3） | 逐点确认下游消费方 |
| 低 | `build_phaded_family_discovery_sets.py:150`、`parse_phaded_cys_*` 6 处、`run_phaded_candidate_hmmsearch.py:197` | 发现层阈值用于"召回/命中/闸门"集合 | 在 manifest 或表头显式声明"recall-only 表，不得作为候选输入"，把口径钉死 |
| 低 | 8 处 `reference_query_only` 包含式选择（§4.4） | R6 潜在面 | 确认调用方不打分；可在选择处直接加断言 |
| 低 | 13 项散文断言（§4.7） | 断言与代码可能不一致 | 人工核对；或把断言改成测试 |

### 8.4 无需处理（正面确认）

注册表只读（R3）、deferred 层无删除且计数未进候选（R4/R5）、发现层标注与 fail-closed 边界执法（R1/R2）均**未发现违规**，见 §5。

### 8.5 并发编辑快照（写作时）

本次审计期间工作树处于被并发任务修改的状态（`git status`）：`pipeline/scripts/08c_tier_rescore.py`、`pipeline/scripts/formal_scan13_tier_processing.sh`、`pipeline/scripts/build_phaded_with_lipase_deferred_tier.py` 为已修改；`run_phaded_candidate_hmmsearch.py`、`build_phaded_family_discovery_sets.py`、`attribute_phaded_candidate_flow_residual.py` 等为未跟踪新文件。因此**行号会随 F1/F3/F4 落地而漂移**，修复时请以本报告的路径＋问题描述为准重新定位。

---

## 9. 复现与验证

```powershell
# 审计表 + summary（写入新目录；目录必须为空或不存在；有违规时退出码 1）
python pipeline/scripts/audit_discovery_layer_usage.py `
  --scripts-dir pipeline/scripts `
  --registry pipeline/config/formal_scan_models.tsv `
  --out-dir <新的空目录>

# 回归测试 / 语法 / 空白
python -m unittest pipeline.tests.test_audit_discovery_layer_usage -v
python -m compileall -q pipeline/scripts pipeline/tests
git diff --check
python pipeline/scripts/audit_discovery_layer_usage.py --help
```

审计产物：`discovery_layer_usage_audit.tsv`（列：`script, line, evidence_kind, signal, use_class, violation, rule_cited, note`）＋ `discovery_layer_usage_audit_summary.json`。
**本次未把 TSV 落入仓库**：F6 的所有权限定为审计器、测试与本报告三个文件，且 `runs/`、`results/`、`deploy/` 为冻结证据；上表与 §3/§4 已完整转载其内容。若需要版本化该 TSV，请由协调者指定路径后再生成（本任务不擅自新增仓库文件）。

**实测数字（2026-09-28 22:05 快照）**：脚本 192（py 172 / sh 20）、findings 184、violation=true 2、unreviewed 105、compliant 77。
**22:12 重跑**：脚本 193（py 173 / sh 20）、findings 185、violation=true 2、unreviewed 106、compliant 77 —— 差异仅来自并发新增的 `build_phaded_v2_catalog_input.py`（见 §0 的重跑漂移说明）。
