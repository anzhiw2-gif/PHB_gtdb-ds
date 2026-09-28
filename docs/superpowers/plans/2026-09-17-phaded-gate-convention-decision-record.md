# Decision record — PhaDED 校准 gate 的两个口径缺陷（规则 #3 与 #4）

**Run：** `20260917_phaded_gate_convention_fix_01`
**日期：** 2026-09-17
**状态：** 已决策并已实现（决策先于实现落盘；实现与测试证据见本轮状态文档 `docs/T141_20260917_phaded_gate_convention_fix_status.md`）
**边界声明：** 本记录及本轮全部产物均为 **candidate-only**。所有 profile、domain、motif、锚点与判别轴证据只表示**候选同源或功能潜力**，**不等同于已验证的 PHB/PHA 降解表型**。本记录**不放宽任何数值门槛**（3/1/1/1 保持不变），只修正**计数口径**与**结构性排除**。

---

## 0. 一句话结论

本轮修正两处**口径缺陷**：(#3) 校准计数器只认"绑定到 profile 自身 family"的记录，使 `20260917_phaded_discrimination_panel_01` 的 960 条四轴判别面板在校准上被完全忽略；(#4) `finalize_phaded_reference_panel_calibration.py` L75 的 `kind == "family"` 使 `profile_kind == "superfamily"` **结构上永远无法通过 gate**，与 `20260917_phaded_classification_authority_01` 确立的"超家族为分类主框架、家族为细化层"相矛盾。修法为：**新增**轴锚定显式字段（不改旧字段语义）+ 为 superfamily 定义**独立且明确**的 gate（保留 family gate 逐字节不变）。**修法不改变任何证据，也不解锁任何 profile** —— 实测 37/37 仍全部不通过，瓶颈在实验阴性为 0 与 held-out 阳性为 0，而非口径。

---

## 1. 事实基线（逐项实测，来源可核）

| 事实 | 实测值 | 来源 |
|---|---|---|
| 参考层冻结台账 | 723 行；`evidence_status` = `annotation_only` 693 + `experimental_positive` 30 + **`experimental_negative` 0** | `runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv`（426,900 B / `0281911c…`，本轮独立重算） |
| 分层权威表 | 46 行 = 8 `superfamily`（`functional_prior`，`registry_eligible=true`）+ 38 `family`（`sequence_clustering_2009`，`registry_eligible=false`） | `runs/20260917_phaded_classification_authority_01/inputs/phaded_classification_authority.tsv`（本轮独立重算） |
| 37 个 reference-only profile 的校准 readiness | `external_bound_negative` 合计 **0**；`heldout_positive` 合计 **0**；`existing_experimental_positive`+`external_bound_positive` 合计 **7**（单 profile 最大 2）；`external_bound_challenge` 合计 579；`family_resolved_negative_status = missing` **37/37** | `runs/20260915_phaded_reference_panel_acquisition_02/results/amendment/profile_calibration_readiness.tsv`（6,112 B / `5a88a4c7…`，37 行，本轮独立重算） |
| 四轴判别面板 | 960 条 = 291 positive_anchor + 4 negative_anchor + 665 challenge；`evidence_status` = `annotation_only` 916 + `experimental_negative` 4 + `experimental_positive` 40；`family_binding_status` = `resolved_existing_family` 379 + `unresolved` 581 | `runs/20260917_phaded_discrimination_panel_01/results/challenge_panel_by_axis.tsv`（本轮独立重算，与同 run `axis_reachability.json` 的 `totals` 逐项一致） |
| 4 条显式实验阴性 | `Q84C08`(axis_1, 负端)、`O87189`/`Q7WT48`/`Q7WT49`(axis_3, 负端)；**4 条全部 `family_binding_status=unresolved`**、`formal_negative_eligible=false`、`gate_credit_slot=none` | 同上（本轮独立重算） |
| gate 现状 | `profiles_meeting_full_gate = 0`、`profiles_blocked = 37`、`gate_relaxation_applied = false`、`profiles_with_built_hmm = 0` | `runs/20260917_phaded_discrimination_panel_01/results/axis_reachability.json` |
| 旧口径输出的可复现性（回归锚点） | 用当前 `finalize` 从上述冻结 readiness 重算，输出与冻结产物**逐字节相同**：6,408 B / `4a04d5eabb9c1ed285349cfe669745f160102ba247c8ff97b620c5aabb228c93` | `runs/20260915_phaded_reference_panel_acquisition_02/results/calibration/leaveout_calibration_decisions.tsv` vs 本轮临时目录重算（未写回任何历史目录） |

---

## 2. 缺陷 #3：gate 计数器只认本家族绑定

### 2.1 缺陷的精确位置与机制（本轮逐行核对）

- `pipeline/scripts/reconcile_phaded_reference_panel.py`
  - L143：`existing_positive = sum(row.phaded_family_id == family and row.evidence_status == "experimental_positive" for row in ledger)` —— 只统计**绑定到该 profile 自身 family** 的台账行；
  - L135–138：`by_profile` 只收 `row["reference_only_profile_id"]` 非空（即 `profile_kind == family` 且 family 命中参考层）的记录；
  - L145–147：`positive` / `negative` / `challenge` 三个计数器**全部**只在这些 `bound_rows` 上求和。
- `pipeline/scripts/acquire_phaded_reference_panel_amendment.py`
  - L184–185：`by_profile[row["profile_id"]]`，而 records 在 L142–145 已被 `family not in target_by_family` 过滤；
  - L189–192：三个计数器同样只在本 profile 的行上求和。

**后果（可量化）：** 判别面板中**按判别轴而非按 family 锚定**的对照对本 profile 的 gate 贡献恒为 0。以本轮实测的 `challenge_panel_by_axis.tsv` 为例：960 条面板记录中 **581 条 `unresolved`**（含全部 4 条显式实验阴性）在旧口径下**不可能**进入任何 profile 的计数器；即使 Task 5 补齐了轴锚定阴性，旧计数器也读不到它们。`axis_reachability.json` 的 `blockers[3]` 已把这一条如实登记。

### 2.2 决策：新增显式轴锚定字段，**不改**旧字段语义

**决策 D3-1（新增字段）。** 在校准 readiness 层**新增**四个字段，与旧字段并列存在：

| 新字段 | 含义 |
|---|---|
| `axis_anchored_negative_count` | 该 profile 在**判别轴上可被归因**的轴锚定**实验阴性**记录数 |
| `axis_anchored_challenge_count` | 该 profile 在**判别轴上可被归因**的轴锚定 challenge 记录数 |
| `axis_anchor_basis_status` | `audited` / `not_audited`；只有 `audited` 的计数才可被任何 gate 读取 |
| `axis_anchored_*_credit_ids` | 可归因记录的 accession 列表（审计用，逐条可回查） |

**决策 D3-2（旧字段语义冻结）。** `external_bound_positive` / `external_bound_negative` / `external_bound_challenge` / `heldout_positive` / `family_resolved_negative_status` / `challenge_status` 的**定义、取值与计算路径一个字都不改**。新字段**只增不改**，因此历史 run 的 readiness 与 decisions 产物在语义上仍与当时一致，历史可比性得到保护（这一点由 §5 的逐字节回归测试强制）。

**决策 D3-3（fail-closed 的锚定依据）。** 新字段只有同时满足**全部**下列条件才计数，任一不满足即**拒绝计数**：

1. `axis_id` 非空**且**存在于冻结的轴定义表（`axis_reachability.json` 的 `axis_definitions`）。未知轴 → `ValueError`（fail-closed，不静默计 0）；
2. `anchor_justification` 非空（可审计的锚定依据文本）；
3. `source_doi_or_pmid` 非空（可回溯到文献标识）；
4. `anchor_evidence_kind` ∈ 该轴该端的**预注册合格证据种类**（取自轴定义表的 `*_qualifying_evidence_kinds`）。
   - 负端**实验阴性**通道只接受 `experimental_negative_assay`；
   - 特别地，`candidate_polar_x1_pattern_state`（563 条候选层极性 x₁ 记录）**结构上不可能**进入阴性通道 —— Task 3 已明确"非疏水 x₁ 不是生物学阴性"，因此它们只能作为 challenge，且必须由轴定义的 challenge 种类白名单放行；
5. **轴对端必须有 family-resolved 锚点**：记录 R 位于轴 a 的 e 端时，被归因的 profile P 必须在 a 的**对端**拥有 ≥1 条 `family_binding_status` 以 `resolved` 开头的面板记录。这一条是整个新口径的**科学内核**：轴锚定对照之所以能作为 P 的对照，是因为 P 与该对照处在**同一判别轴的两端**（即形成了一个真实的判别对）；同端记录、无对端锚点的记录一律不计。

**决策 D3-4（family 层与 superfamily 层分别归因）。**
- family profile P（`phaded_family_id = F`）：对端锚点必须是**绑定到 F 本身**的记录；
- superfamily profile S：对端锚点可以是**绑定到 S 的任一成员 family** 的记录（超家族 = 其成员家族的并集，见 §3.2）。
- 归因按**记录**去重计数，不按 family 重复计数。

### 2.3 为什么这是**口径修正**而不是"放宽标准"

1. **数值门槛一个没动。** gate 的 `3 / 1 / 1 / 1` 与"零未解释命中"四项阈值逐字未变（`finalize` 的 `report["gate"]` 与 `gate_convention_impact.json` 的 `thresholds` 两处同时记录）。
2. **旧字段没有被稀释。** 报告口径没有变成"把轴对照塞进 negative 槽位"，而是**新增一个独立字段**。family gate 读的仍然只有旧字段（见 §3.1），因此 family 层的判定**与修改前逐字节等价**（§5 回归测试）。
3. **新字段的可采信条件比旧字段更严，不是更松。** 旧字段只要求"绑定到自身 family"；新字段额外要求轴定义存在、锚定依据非空、文献标识非空、证据种类在白名单内、且存在对端 family-resolved 锚点 —— 五个条件**同时**成立才计 1。任何一条缺失都 fail-closed 归零。
4. **不产生新证据。** 新字段是对**既有冻结面板**的重新计数，不新增任何 accession、不改任何 `evidence_status`、不把 `annotation_only` 升级为实验证据（`693` 条 `annotation_only` 记录依旧不得进入任何严格证据槽位）。
5. **不改历史。** 见 §4。

### 2.4 证伪条件（可执行）

- **F3-a（污染检测）**：若修法后**任何既有 profile 在证据未变的情况下新通过 gate**，则该修法被污染，**必须回退**。本轮实测：新口径下 `calibrated_candidate_model = 0`，37/37 仍不通过。
- **F3-b（旧字段漂移）**：若 `external_bound_negative` / `external_bound_positive` / `external_bound_challenge` / `heldout_positive` 中任何一个的值与本轮之前对同一冻结输入的重算结果不一致，则该修法被污染，**必须回退**。由 `test_legacy_counters_are_not_changed_by_axis_credit` 强制。
- **F3-c（无锚定依据的计分）**：若任何缺少 `axis_id` / `anchor_justification` / `source_doi_or_pmid` 中任一字段的记录仍被计分（或未知轴被静默计 0 而非 raise），则 fail-closed 契约失效，**必须回退**。由 `test_missing_anchor_justification_is_rejected_fail_closed`、`test_missing_source_doi_or_pmid_is_rejected_fail_closed`、`test_unknown_axis_is_rejected_fail_closed` 强制。
- **F3-d（把非阴性证据当作阴性）**：若任何 `candidate_polar_x1_pattern_state` 或 `reported_label` 记录被计入 `axis_anchored_negative_count`，则该修法越界为放宽，**必须回退**。由 `test_polar_x1_pattern_never_counts_as_negative` 强制。

---

## 3. 缺陷 #4：superfamily 结构上永远无法校准

### 3.1 缺陷的精确位置

`pipeline/scripts/finalize_phaded_reference_panel_calibration.py` L75：

```python
gate = kind == "family" and positive_count >= 3 and heldout >= 1 and negative >= 1 and challenge >= 1 and row["family_resolved_negative_status"].strip() == "sufficient" and row["challenge_status"].strip() == "sufficient"
```

`kind == "family"` 使 4 个 reference-only 超家族（`extracellular native-SCL/PhaZ7-like`、`intracellular nPHASCL with lipase box`、`intracellular nPHASCL without lipase box`、`periplasmic PHA depolymerases`）**无论证据多充分都不可能通过**。这与分层权威表（§1：8 个超家族 `functional_prior` + `registry_eligible=true`，38 个家族 `sequence_clustering_2009` + `registry_eligible=false`）以及 Task 1 的结论"**超家族是分类主框架、家族是细化层**"**直接矛盾**：主框架层反而没有校准判据。

### 3.2 决策：为 superfamily 定义**独立且明确**的 gate，family gate 原样保留

**决策 D4-1（family gate 逐字节保留）。** L75 的合取式**一个字不改**，其输入仍是旧字段。任何对 family 行的输出（`decision` / `calibration_status`）必须与修改前**逐字节相同**。

**决策 D4-2（新的 superfamily gate，独立且明确）。** 对 `profile_kind == "superfamily"` 的行，判定为：

| 槽位 | 新 gate 的取值来源 | 阈值（**未变**） |
|---|---|---|
| positive | `superfamily_positive_count`（= Σ 成员 family 的 `existing_experimental_positive + external_bound_positive`）；列缺失时回退到 `positive_count`（对 superfamily 行为 0） | ≥ 3 |
| held-out positive | `superfamily_heldout_positive`（= Σ 成员 family 的 `heldout_positive`）；列缺失时回退到 0 | ≥ 1 |
| negative | `axis_anchored_negative_count`（§2） | ≥ 1 |
| challenge | `axis_anchored_challenge_count`（§2） | ≥ 1 |
| 零未解释命中 | `zero_unexplained_hits_status == "sufficient"`；列缺失或不等于 `sufficient` → 不通过 | 必须满足 |
| 锚定依据 | `axis_anchor_basis_status == "audited"` | 必须满足 |
| 新 family 判定 | `new_family_call_created` ∈ {false,0,no}，否则 **raise** | 必须为 false |

**为什么 positive / held-out 在超家族层用成员家族并集**：超家族在权威表中**就是**其成员 family 的并集，因此"成员家族里的独立实验阳性"在超家族分辨率下按定义就是超家族的阳性；这是**分辨率适配**，不是新增证据。为避免任何隐藏的放宽，本轮在输出里**同时**并排给出**原始 family-bound 计数**与**超家族口径计数**两列（`gate_detail`），使二者的差异**逐行可见**。
**为什么 negative / challenge 在超家族层读轴锚定字段**：超家族级不需要"绑定到某一个家族"这一更细的分辨率要求，只需要**轴对端**这一层级的判别对；这正是 §2.2 第 5 条的产物。

**决策 D4-3（分辨率必须在输出里显式标注）。** 超家族级校准的语义是**"超家族分辨率判定"**，比 family 级**粗**（同一超家族内的家族差异不被该判定区分）。因此：
- 新增输出 `leaveout_calibration_gate_detail.tsv` 携带 `gate_kind` 与 `gate_resolution`（`family` / `superfamily`）两列，**逐行显式标注分辨率**；
- 报告 JSON 的 `gate` 块拆为 `gate.family` 与 `gate.superfamily` 两个子块，各自声明 `resolution`；
- **禁止**把 `gate_resolution = superfamily` 的通过结果表述为 family 级判定，也**禁止**据此写任何 `subtype_call` 或表型结论。

**决策 D4-4（保留 fail-closed 断言）。** `new_family_call_created` 非 false 时 **raise**（对两种 kind 都生效，与修改前一致）。这一断言是整条评测链的兜底：即使某个 superfamily 新通过了新 gate，它也**不能**因此产生任何新的 family 判定。

### 3.3 为什么这是**结构性排除的移除**而不是"放宽标准"

1. 阈值 `3/1/1/1` 与"零未解释命中"**逐字未变**（两处记录，见 §2.3 第 1 条）。
2. **family 层完全未动**：family 行的 gate 表达式与输出逐字节等价（§5 用冻结输入与冻结产物的双回归强制）。
3. **superfamily 层的每个槽位都有明确的证据来源、明确的 fail-closed 回退值、且输出里逐行标注分辨率**；缺失任何输入列时**回退到阻塞值**（0 / `not_audited` / `not_assessable`），因此"缺少证据"只会导致**不通过**，不可能导致**通过**。
4. **本轮无任何 profile 新通过**（§5 实测）：4 个 reference-only 超家族的 positive 并集各为 1（< 3）、held-out 各为 0（< 1），即使 negative 槽位由 Q84C08 的轴锚定得到 1，gate 仍在 positive 与 held-out 两处失败。
5. 修法**只增加可评估性**：旧口径下超家族连"如何评估"都不存在（结构上不可能），新口径下 gate 可评估、可解释、可预期在证据到位时通过 —— 这正是移除"结构性排除"的定义。

### 3.4 证伪条件（可执行）

- **F4-a（污染检测，与 F3-a 同）**：若修法导致任何既有 profile 在证据未变的情况下新通过 gate，**必须回退**。本轮实测 `calibrated_candidate_model = 0`。
- **F4-b（family 漂移）**：若用冻结 readiness 重算得到的 `leaveout_calibration_decisions.tsv` 与冻结产物（6,408 B / `4a04d5ea…`）**不是逐字节相同**，则 family gate 被改动，**必须回退**。由 `test_family_gate_is_byte_identical_on_frozen_readiness` 强制。
- **F4-c（fail-closed 回退值）**：若任一缺失列导致 superfamily gate **通过**（而不是不通过），则 fail-closed 契约失效，**必须回退**。由 `test_superfamily_gate_fails_closed_when_columns_are_absent` 强制。
- **F4-d（分辨率标注缺失）**：若任何 gate detail 行缺少 `gate_resolution`，或其值不等于该行 `profile_kind` 对应的分辨率，则判定不可被正确解读，**必须回退**。由 `test_gate_detail_labels_resolution_for_every_row` 强制。
- **F4-e（断言被绕过）**：若 `new_family_call_created` 为 true 时不再 raise，则兜底断言失效，**必须回退**。

---

## 4. 追溯性声明（历史不可改写）

1. **历史 run 的既有判定不被追溯修改。** `runs/20260901_*`–`runs/20260917_*` 中任何既有产物（含 `leaveout_calibration_decisions.tsv`、`profile_calibration_readiness.tsv`、`axis_reachability.json`、`challenge_panel_by_axis.tsv`）**本轮一律只读**，其字节、判定与"当时口径"的含义**保持不变**。旧口径下"37/37 不通过"的历史结论**依然成立且依然可引用**。
2. **新口径只对新的 dated run 生效。** 本轮新建 `runs/20260917_phaded_gate_convention_fix_01/`（以及未来任何 dated run）在新口径下的重算结果，**不覆盖、不回填、不取代**任何历史 run 的产物；两套口径的差异在 `results/gate_convention_impact.json` 中**并列呈现**。
3. **阈值与判据的权威声明位置不变**：`pipeline/config/formal_scan_models.tsv`（本轮逐字节未变）、冻结台账、分层权威表均未改动。
4. **不得据本记录删除或降级任何候选**：本记录不产生任何 `subtype_call` 变更，也不构成候选删除依据（candidate-only）。

---

## 5. 实现与强制（Test-first）

新增/修改的代码与测试（先写失败测试，红灯后再实现）：

| 文件 | 作用 |
|---|---|
| `pipeline/scripts/credit_phaded_axis_anchored_controls.py`（新增） | §2.2 的轴锚定归因 + §3.2 的超家族口径并集；全部 fail-closed 分支在此 |
| `pipeline/scripts/finalize_phaded_reference_panel_calibration.py`（修改） | family gate 逐字保留；新增 superfamily gate；新增 `leaveout_calibration_gate_detail.tsv`；`axis_anchored_*` 无 `audited` 依据时 raise |
| `pipeline/scripts/reconcile_phaded_reference_panel.py`（修改） | **仅在显式提供轴面板输入时**追加轴锚定字段；未提供时输出列与数值逐字节不变 |
| `pipeline/scripts/acquire_phaded_reference_panel_amendment.py`（修改） | 同上 |
| `pipeline/tests/test_credit_phaded_axis_anchored_controls.py`（新增） | F3-a…F3-d |
| `pipeline/tests/test_finalize_phaded_superfamily_gate.py`（新增） | F4-a…F4-e + family 逐字节回归 |

**强制回归（本轮实测证据见状态文档）：** `python -m unittest discover -s pipeline/tests`、`python -m compileall -q pipeline`、`git diff --check`、以及用冻结 readiness 重算 decisions 的逐字节比对。

---

## 6. 明确的"不做"

❌ 不修改 `pipeline/config/`（尤其 `formal_scan_models.tsv`）❌ 不修改任何历史 run 的任何文件 ❌ 不放宽 `3/1/1/1` 或"零未解释命中"❌ 不把 `annotation_only` 升级为实验阴性（693 条依旧禁入严格槽位）❌ 不把非疏水 x₁ 或 `reported_label` 当作阴性 ❌ 不为任何 profile 生成 family 判定 ❌ 不产生 `subtype_call` 变更 ❌ 不据此删除或降级任何候选 ❌ 不执行 `git commit` / `push` / `clean` / `reset` ❌ 不删除任何文件 ❌ 不启动服务器计算。

---

*本决策记录为 candidate-only 文件。其全部结论只涉及"如何计数既有证据"，不涉及任何 PHB/PHA 降解表型的成立与否。*
