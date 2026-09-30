# PHB_gtdb-ds 项目审计快照（2026-09-29）

> **本文档身份**：面向审计的当前项目状态总览。所有数字均取自冻结产物或本次实测复核（复核命令见文末「如何核验」）；不来自记忆、不含推断、不含 DOI 占位符。
>
> **依据**：
> - 权威状态：`docs/T141_20260928_phaded_evidence_model_redesign_status.md`
> - 发布记录：`docs/T141_20260928_phaded_v2_release_notes_and_commands.md`
> - 计划：`docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign-followup.md`（F1–F17）
> - 逐轮交接：`docs/T141_20260917_project_handoff.md` §12（共 55 个小节，2,032 行）
> - 本快照生成时实测：`origin/main = 5e96c9b`、工作区 0 改动、frozen 目录（`runs/ results/ deploy/`）0 改动、HEAD 中服务器路径字面量 0 处。

---

## 0. 一句话总结

**F1–F17 已按计划依次执行完毕**：16 项任务完成且逐项经机械清册核验；F17 的 push / tag / GitHub Release 完成，**Zenodo DOI 由操作者于 2026-09-29 决定暂缓**（既非待办亦非阻塞）。**两条结构证据线（P4 普查、P5 尾部）均已跑完并得出互为补充的结论**：序列标签在**小边际（≤约 14 bits）无结构担保、大边际（≥约 157 bits）有很强担保**。尚待操作者决定的事项为 4 项（历史重写、472 条 disposition、`hfam_70` 传播、目录 schema），均已在文档中留痕并附默认状态。

---

## 1. 任务完成状态（F1–F17 逐项）

> 本表每一行都有机械核验支撑：`runs/20260928_phaded_f_task_inventory_01/verify_inventory.py` 检查 run 目录存在性与 commit 可解析性，判定 `INVENTORY COMPLETE`（17/17 任务可追溯，无一缺失）。

### A 类：基础修正

| 任务 | 结果 | 产物 / 证据 |
|---|---|---|
| **F1** 616 条候选流转残留归因 | **完成** | `runs/20260928_phaded_candidate_flow_616_attribution_01`、`runs/20260928_phaded_616_disposition_resolution_01`。**发现词表缺口**：613 条 `remote_homolog_candidate` 中 **472 条既无发现层行也无 profile 判定**，六项词表无适配项（见 §4 决策②） |
| **F2** gRodon 66 条差额分桶 | **完成** | `runs/20260928_phaded_grodon_66_reconciliation_01`。**闭合式 `4,507 = 4,441 + 23 + 43`**；此前「26,116 本地可闭合」的前提被实测推翻 |
| **F3** 三个下游脚本改读 `model_layer` | **完成** | commit `a5d2cc4`（序列发现与功能校准分离；`functional_calibration_status` 独立字段） |
| **F4** tier 处理绑定同一全库 `-Z` | **完成** | commit `417e838`（测试覆盖；入口无 `--database-size-z` 即 fail-closed） |
| **F5** 结构脚本解硬编码 | **完成** | commit `417e838`（`prepare_phaded_structure_review.py`） |
| **F6** 发现层 / deferred 层下游误用审计 | **完成** | `runs/20260928_phaded_catalog_schema_audit_01` |

### B 类：数据与证据运行

| 任务 | 结果 | 产物 / 证据 |
|---|---|---|
| **F7** scan-13 `-Z` 尺度对账 | **完成（本会话最重要的数字更正）** | `runs/20260928_phaded_scan13_z_scale_reconciliation_01`。**实测全库 Z = 615,969,589**（旧 `~2.92×10⁸` 系外推，偏低约 2.1 倍，已作废）。冻结 6,743,197 行按真实尺度重标后，通过 `1e-5` 的行由 **6,742,621 → 6,305,428**，即 **437,193 行（6.48%）新被拒**（逐家族最大落差 `ePhaZ_broad_discovery` −205,106；`phasin` 7 行全部出界）。**冻结表一字未改，新旧 E 值不混用** |
| **F8** P1 参考层 SignalP | **完成** | `runs/20260928_phaded_reference_signalp_v2_01`（732 个结果文件） |
| **F9** P2 type-1 SP-less 配件域 | **完成** | `runs/20260928_phaded_type1_sp_less_accessory_v2_01` |
| **F10** P3 Cys 四态 | **完成** | `runs/20260928_phaded_cys_anchor_states_v2_01`（substitution / truncation / uncertain 明确标注） |
| **F11** P4 nPHAMCL 竞争面板 | **完成（序列层 + 结构层均已收官）** | 序列层 `runs/20260928_phaded_nphamcl_competition_panel_v2_01`（987 条**全部**同时显著命中锚与竞争者，评分器已自带 `READ WITH CARE` 声明）；结构层 `runs/20260928_phaded_structural_survey_v2_01` —— **30/30 完成，结果见 §2** |
| **F12** P5 with-lipase 延迟池 | **完成（全链路 + 结构 tranche 已收官）** | 五个 run（分层 / 尾部全取 / 序列提取 / 序列层竞争 / 结构 tranche）。结构 tranche `runs/20260928_phaded_deferred_tail_structural_v2_01` —— **20/20 完成，结果见 §2** |
| **F13** 723 行台账证据分级 curation | **完成** | `runs/20260928_phaded_reference_evidence_curation_01` |
| **F14** v2 候选目录真实重算 | **完成** | `runs/20260928_phaded_v2_catalog_recompute_01` + `runs/20260928_phaded_v2_catalog_input_delta_verify_01`（**增量性证明**：adapter 修订为纯增量，0 值失配于约 530 万单元格比较）。交付目录为固定 **20 列** schema |
| **F15** `hfam_70` 提升 | **完成（提升 + 传播缺口均已量化）** | `runs/20260928_phaded_hfam70_promotion_01`（46 个 profile 中**恰好 1 个**被提升：`family_DED_hfam_70_9818be7f78e3`）+ `runs/20260928_phaded_promotion_propagation_audit_01`（见 §4 决策③） |

### C 类：发布

| 任务 | 状态 | 证据 |
|---|---|---|
| **F16** 研究快照入版本控制 | **完成** | commit `c1b4155`、`b2a01b4` |
| **F17** push | **完成** | `origin/main = 5e96c9b`，0 ahead / 0 behind |
| **F17** tag | **完成** | annotated tag `phaded-evidence-model-v2-20260928` |
| **F17** GitHub Release | **完成** | Release **id `398452852`**，公开、prerelease（与 tag 的 release-candidate 措辞一致），发布说明 = tag 注解逐字（1,802 字符），发布后经 API 取回复核 **0 身份违规**。URL：https://github.com/anzhiw2-gif/PHB_gtdb-ds/releases/tag/phaded-evidence-model-v2-20260928 |
| **F17** Zenodo DOI | **操作者决定暂缓（2026-09-29）** | **DOI 不存在，且任何地方未写占位符**（编造标识符比留空更糟）。此前实测：本机无 token；Zenodo 公开 API 该账号 0 条记录；发布未铸出 DOI ⇒ 集成未开启。**将来要铸的路径不变**（Zenodo 网页端开启该仓库的 GitHub 集成后重新发布）；本仓库内无代码路径可完成它 |

---

## 2. 两条结构证据线的最终结果

### 2.1 P4 结构普查（30/30，完整）

**产物**：`runs/20260928_phaded_structural_survey_v2_01/results/FINAL_survey_report.md`；评分摘要 `full_read_scores/survey_structural_summary.json`。Foldseek 210 行（30 × 7 面板成员）比对，exit 0；日志零错误；十分位 0–9 全覆盖。

| 量 | 值 |
|---|---|
| 判给竞争者 / 锚 | **30 / 0** |
| 竞争者构成 | `Q88N36` ×27、`Q02104` ×2、`P24640` ×1 |
| 锚 TM | **0.720 – 0.772** |
| 竞争者 TM | **0.773 – 0.804** |
| 结构边际（锚 − 竞争者） | −0.068 … −0.013（中位 −0.042；跨度 **0.0554 TM**） |
| 序列边际跨度 | **13.80 bits**（−14.2 … −0.4） |
| 越过预注册 `TM ≥ 0.5` | 30 / 30 |
| Spearman | +0.123；n=30 阈值 0.361 ⇒ **不显著** |

**要点**：两个 TM 区间**不重叠**（锚最好 0.772 < 竞争者最差 0.773）。13.80 bits 的序列偏好扫描只让结构裁决移动 0.042 TM，且一次未改变。⇒ **小边际下，`nphamcl_like_supported` 标签不携带结构担保。** 这是关于**标签**的陈述：987 条仍为 `function_unresolved`，无任何删除/降级/排除。

### 2.2 P5 尾部结构 tranche（20/20，完整）

**产物**：`runs/20260928_phaded_deferred_tail_structural_v2_01/results/FINAL_P5_report.md`；评分摘要 `tail_structural_scores/survey_structural_summary.json`。Foldseek 700 行（20 × 5 模型 × 7 靶标），sha256 `816e1eaa…`；完成标记 `P5_COMPLETE.json` `complete: true`；耗时 2h41m（17:05:32Z → 19:46:49Z）。

| 量 | 值 |
|---|---|
| 判给锚 / 竞争者 | **10 / 10** |
| 序列边际跨度 | **312.9 bits**（−233.0 … +79.9） |
| 结构边际跨度 | **0.3991 TM** |
| 越过预注册 `TM ≥ 0.5` | 19 / 20 |
| Spearman | **+0.872；显著且可解读** |

**要点**：关系单调——**锚偏好最强的 8 条（−233 至 −157 bits）全部 8 条判给锚**（结构边际 +0.198 … +0.232 TM，已独立验证）；中间带（−41 至 −6）混合、边际小；最偏竞争者的一条（+79.9）以 −0.167 TM 判给竞争者。**交叉验证**：普查的全模型文件（1050 行）与 rank_001 文件（210 行）评分**逐项相同** ⇒ P5 与普查的差异是**人群属性**而非文件形态。

### 2.3 两个结果合并后的正确结论

> **低于约 14 bits 边际时序列标签没有结构担保；高于约 157 bits 时它有很多。**

单看任一 tranche 都会得出错误的推广。两者**互斥抽样**（0 个共享 accession，`runs/20260928_phaded_tranche_disjointness_01` 已核验），这正是解读前先验互斥性的价值。

**边界（不变）**：折叠均为**预测**（ColabFold MSA 模式），非实验结构；P5 的 n=20 且 10/10 对分；P5 尾部属 **with-lipase 延迟池**（已知含噪层），**其成员不是候选，本结果不改变任何处置**；P5 有一条锚 TM 0.414 < 0.5，预注册阈值在此不是干净分隔符（分隔靠边际）。

---

## 3. 本会话对既有记录的关键修正（审计要点）

本会话除执行任务外，还**修正了若干此前已公开的错误记录**。审计应知悉这些修正的存在与理由：

| # | 修正 | 严重性 | 位置 |
|---|---|---|---|
| 1 | **公开 README 两处低报自身状态**：①「no family has reached ≥3 independent positives」——按唯一独立性分组重算 `DED_hfam_70` 达 7、`DED_hfam_52` 达 6；②「hfam_70 未被最终化」——F15 已在授权下提升为 `calibrated_candidate_model`。另把已完成的 Release 从待办中勾掉 | 高（公开门面） | README L5/L10/L185 |
| 2 | **`08c_tier_rescore.py` 用法示例含被推翻的旧 Z**（292,000,000）——照抄者会让 E 值过于宽松 2.1 倍，恰是 Z 契约要防止的失败模式 | **高（运维隐患）** | 已改为实测值并加警告；见 `4cd32a6` |
| 3 | 权威状态文档三项 `pending` 被本会话关闭（全库 Z、逐片计数、重标尺表本地镜像） | 中 | 状态文档 §L258–269 |
| 4 | 两份按旧 Z 计算比率的 dated 文档（0.012% → **0.0058%**；0.52% → **0.248%**）加更正横幅，正文未改 | 中 | `20260917_phaded_cys_targeted_recall_status.md`、`20260918_phaded_family_discovery_status.md` |
| 5 | **「结构层被硬约束阻塞」的说法被推翻**：该约束只对 `--msa-mode single_sequence` 成立；MSA 模式实测 pLDDT 92.6–93.4 | 高（影响资源决策） | 交接 §12.14；早先段落已加前向指针 |
| 6 | **两个评分器会输出易被误读的统计量**：结构评分器的 Spearman（小 n、窄带宽下无意义）与序列评分器的 `nphamcl_like_supported`（987/987 同时显著命中两侧时不构成判别）。已在工具内做守卫（`spearman_is_interpretable` / `spearman_is_significant` / `READ WITH CARE`） | 中 | `403c5ce`、`05157aa` |
| 7 | **身份门禁事故与修复**：commit `414b324` 曾把服务器路径字面量带入公开历史并被推送；已由 `538f331` 修复（HEAD 0 处），`653cfdd` 如实记录。经实测该字符串**自 `8162e38`（Initial commit）起即在公开历史、横跨 10 个 commit** ⇒ 不是新暴露类别 | 高（已修复） | 交接 §12.45/§12.46 |
| 8 | 「668」数字撞车消歧（F14 几何缺口所指的 accession 数 vs 流程残差 `668 = 616 + 52`） | 低 | 交接 §12.37(d) |
| 9 | 交接文档自身结构缺陷（两段孤儿文字）修复 + 结构/索引核验器落地 | 低 | §12.40 |

---

## 4. 待操作者决策的事项（4 项）

> 每一项都**已带证据记录在案，且保守默认（不做）已在生效**。未决 ≠ 阻塞；当前所有交付物在这些默认下自洽。

| # | 事项 | 证据要点 | 当前默认（生效中） |
|---|---|---|---|
| ① | **git 历史重写** | 公开历史中 **94 个**身份字符串违规 blob，最早来自 `8162e38`（Initial commit），横跨约 24 个 commit；**HEAD 与工作区 0 违规**（`runs/20260928_phaded_git_history_identity_audit_01`，只读实测）。重写需明确授权 | **不重写**。若决定做，应一次性处理全部 94 个而非个别 commit |
| ② | **472 条的 `primary_disposition`** | 613 条 `remote_homolog_candidate` 中 472 条既无发现层行也无 profile 判定，而该标签断言「仅发现层支持」——对它们不成立；六项词表**无适配项** | **保留原标签**（不把「词的缺口」用次优标签掩盖）；可选方向：为词表加第 7 项、或映射到 `function_unresolved`（F1 的提案） |
| ③ | **`hfam_70` 提升是否传播到 1,087 行** | 提升的是 1 个 profile；1,087 行目录行名它：仅 **132** 条为 `profile_trained_hit`/`assigned`，**955** 条为 `ambiguous`；`family_score_gap` 在 **1.0 bits** 处把它们完全分开（0.0–0.9 vs 1.0–3.3）。当前适配器**无提升传播输入**（结构性，非疏忽） | **不传播**（1,087 行保持 `discovery_hmm_uncalibrated`） |
| ④ | **目录 schema 是否承载命中强度** | 交付目录为固定 **20 列**；`evidence_flags` 词表纯分类 | **不改**；若改需另立版本并重做增量性验证 |
| ⑤ | **是否发 v2.1 Release**（§12.32 的 Option B） | Option A（发布既有 tag）已执行；v2.1 会在当前 main（含本会话全部修正）上新建版本标记。发布说明草案已备 | **不发**（无新版本决策需求） |

**已了结（不再需要决策）**：**Zenodo DOI —— 操作者 2026-09-29 决定暂缓**（§12.55；既非待办亦非阻塞；README 的项目级「论文 + profiles + DOI」规划保留不变，那是与本会话任务范围不同的另一件事）。

---

## 5. 治理与完整性状态（生成快照时实测）

| 项 | 实测值 |
|---|---|
| 发布门禁 | `python pipeline/scripts/run_release_gate.py` **5 项全过**（字面路径 / diff-check / compileall / 身份凭据门禁 / 全量单测）。全量单测 **1,912 项 OK（skipped 1）** + 门禁模块自身 11 项 |
| 交接核验器 | 索引核验：§12 共 **55** 小节，引用与实际**一一对应**；结构核验：小节连续升序、无重复、每节有正文、落款后无内容 —— 全绿 |
| frozen 目录 | `runs/`、`results/`、`deploy/` **0 改动**（本次实测） |
| 工作区 | **0 改动**；`origin/main` 0 ahead / 0 behind |
| 身份 | HEAD 中服务器路径字面量 **0 处**；公开门禁 7 项 OK |
| v1 权威件 | 会话中途（§12.33）实测三份 v1 权威件哈希与会话开始时**逐位一致**（`F4A9D99039EF…`、`7772B34C6A9F…`、`7B6F9208CFCC…`） |
| 冻结扫描注册表 | `pipeline/config/formal_scan_models.tsv` sha256 `8A1C0508…` 与第 13 轮审计记录一致；无 `DED_hfam_*` 条目 |

---

## 6. 本会话产出清单（供审计抽样）

- **新增脚本 13 个**（§12.36 索引，各配失败测试；测试模块 11 个，全部经存在性核验）
- **新增 run 33 个**（`runs/20260928_phaded_*`，含 F1–F17 的 25 个任务 run 与 8 个核验/审计 run）
- **关键冻结产物**：
  - `scale_delta_report.json`（z_total=615,969,589；delta=−437,193）
  - `FINAL_survey_report.md`（30/30 竞争者；TM 区间不重叠）
  - `FINAL_P5_report.md`（10/10；ρ=+0.872 显著；含全模型交叉验证）
  - `P5_COMPLETE.json`（20/20，complete: true）
- **门禁设施**：`pipeline/scripts/run_release_gate.py`（五项门禁、首个失败即停、退出码即决定、真实树上通过）

---

## 7. 如何核验本文件

1. **门禁**：`python pipeline/scripts/run_release_gate.py`（退出码为准）
2. **任务清册**：`python runs/20260928_phaded_f_task_inventory_01/verify_inventory.py` → `INVENTORY COMPLETE`
3. **交接一致性**：`python runs/20260928_phaded_handoff_index_verify_01/verify_index.py` 与 `python runs/20260928_phaded_handoff_structure_verify_01/verify_structure.py`
4. **两个结构结论**：分别读两个 `FINAL_*_report.md`，与同 run 的 `survey_structural_summary.json` 对照
5. **F7 数字**：`runs/20260928_phaded_scan13_z_scale_reconciliation_01/results/scale_delta_report.json`
6. **发布**：https://github.com/anzhiw2-gif/PHB_gtdb-ds/releases/tag/phaded-evidence-model-v2-20260928 （Release id 398452852）

*本文档为 candidate-only 审计记录。所有被引用的 profile、domain、motif、SignalP、结构、定位与系统发育证据仍只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。*
