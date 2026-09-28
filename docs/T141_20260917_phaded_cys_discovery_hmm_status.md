# T141 — Cys 型发现层 HMM（`discovery_hmm_uncalibrated`）训练与候选池打分：状态文档

**任务：** 计划 v2 Task 6 的后续授权轮 —— 人类操作者于 2026-09-17 授权"可以训练 hmm"，范围 = Cys 型家族（`intracellular nPHASCL without lipase box`）+ 对现有 109,087 条候选池打分
**Run：** `runs/20260917_phaded_cys_discovery_hmm_01`（本地与服务器同名；未触碰任何历史 run、`results/`、`pipeline/config/`、AGENTS.md、历史 HMM）
**Decision record：** `docs/superpowers/plans/2026-09-17-phaded-inphascl-cys-decision-record.md`（§9 为本轮追加的授权修订，§1–§8 原文未改）
**日期：** 2026-09-17
**服务器：** **已使用** —— 仅 CPU（MAFFT + HMMER 3.4），单任务线程 ≤ 40，未使用 GPU
**状态：** `completed_candidate_only` —— 预注册闸门 **G1–G4 全部 passed**（**但 G2 是样本内自洽度量，不是泛化估计**；见 §6）

> **边界声明（贯穿全文）：** 本轮产物是 **发现层、未校准** 的同源/功能潜力证据。命中**不是** family 判定，**不是**校准过的 profile，**不进 registry**，对 109,087 条候选的 `subtype_call` 影响 **0 行**，且**不等同于已验证的 PHB/PHA 降解表型**。

---

## 0. 一句话结论

在人类操作者授权下，用**冻结 DED 台账的 270 条 Cys 型成员**（Knoll 2009 式 profile 构建，允许含 `annotation_only`）做了 1 张 MAFFT 比对并建出 1 个**未校准发现层 HMM**（`NAME cys_discovery_uncalibrated`，451 个模型位点）。在**预注册且钉死**的 E < 1e-5（`-Z 109087`）阈值下：**109,087 条候选中 31,906 条命中**；723 条参考 + Q84C08 的 724 条打分集中**恰好只有 276 条 Cys 型参考命中、其余 448 条（含 365 条胞外参考与 Q84C08）一条也不命中**。**G1–G4 全部通过**，但必须与下列三条诚实限定一起读：**(a)** G2 的 270/270 = 1.00 是**样本内**召回（训练成员自身即分母），**不是**泛化能力；**(b)** x₁ 疏水混淆集 563 条中有 **443 条命中（78.69%）**，说明该发现层**很宽**，不具 x₁ 分辨力；**(c)** 校准 gate **未变、仍未达成**（本 run 的审计只把 Cys 家族实验阳性记录数从 1 提到 **2**，仍 < 3）。

---

## 1. 授权与治理设计

| 项 | 内容 |
|---|---|
| 授权 | 2026-09-17 人类操作者授权"可以训练 hmm"，范围 = Cys 型家族 + 对现有 109,087 条候选池打分 |
| 新模型层命名 | **`discovery_hmm_uncalibrated`**（发现层 HMM，未校准）—— 逐处标注于 HMM 文件（`NAME cys_discovery_uncalibrated`）、训练 FASTA 头、训练集 manifest、预注册文件、打分表、claims 表、`gate_evaluation.json`、`discovery_manifest.json` |
| 训练数据口径 | 允许使用 DED 家族序列（含 `annotation_only`）。依据：**DED 自身就是这么做的**（Knoll 2009 为每个 family 提供 HMMER profile 用于 in silico 鉴定），因此这是**复现文献方法**，不是把 `annotation_only` 升格为实验证据 |
| 训练 ≠ 校准 | 未做阈值校准、未做 held-out 校准，**不声明敏感度/特异度**；G2 明确标为 `in_sample_not_held_out` |
| 训练 ≠ 进 registry | `pipeline/config/formal_scan_models.tsv` **逐字节未变**（运行前/后均为 450 B / `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`，见 §6 G3）；模型文件只存在于本轮 dated run |
| 命中仍是 candidate-only | 打分表每行 `family_call_made=false`、`subtype_call_impact=none`；claims 表每行另带 `new_family_call_made=false`；**对 `subtype_call` 影响 0 行** |
| `annotation_only` 禁训规则 | `validate_phaded_reference_ledger.py` L172–173 **继续有效**；本轮**不申请**严格训练资格，只以 DED 家族成员复现文献的 profile 构建方式 |
| 校准 gate | 每家族 ≥3 独立阳性 + ≥1 held-out + ≥1 family-resolved 阴性 + ≥1 challenge + 零未解释命中 —— **未触碰、仍未达成** |

完整授权修订文本见 decision record **§9（纯追加，§1–§8 一字未改）**，其中 §9.3 逐条说明 §5.1/§5.2 与本节的适用口径（§5 讨论**校准层/registry 层**，本轮授权的是**发现层**，两者不互相推翻）。

---

## 2. 训练集构成与排除理由

来源：冻结台账 `runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv` 中 `phaded_superfamily == "intracellular nPHASCL without lipase box"` 的 **276 条**记录。**台账本身只读，未改动**（哈希见 §10）。

| 判据 | 规则 | 结果 |
|---|---|---|
| 定位 | `reported_localization` 必须为 `intracellular` | 276/276 通过，**排除 0** |
| SignalP | 冻结参考层**没有任何 SignalP 运行**（候选层 SignalP 表里参考 accession 命中数为 **0**）→ 该证据记 `pending` 并**保留序列**、逐行标注 `signalp_evidence=pending_reference_layer_has_no_signalp_run` | **排除 0**，276/276 标注 pending |
| lipase box 注释冲突 | 超家族定义为"无 lipase box"；若冻结 mapping manifest 记 `lipase_box_state=supported`（且该超家族不预期 lipase box）**并且**在本 run 内用 Knoll 2009 的规范模式 `Gx1Sx2G` 独立复现 | **排除 3**：`DED_hfam_61_0004`（ABS67380.1，`GLSHG`@42 邻域）、`DED_hfam_61_0009`（EBA09045.1）、`DED_hfam_62_0004`（AAO62349.1）——三者**恰好等于**冻结层标出的 3 条，独立复现零分歧 |
| 长度（完整性） | 项目既有先例 `split_ephaz_seeds.py` 的 `minimum_length_for_core = 200`（策略 `retain_in_broad_and_review`） | **排除 3**：`DED_hfam_61_0043`（94498690，147 aa）、`DED_hfam_67_0044`（167587491，163 aa）、`DED_hfam_64_0011`（AAL03581.1，197 aa）。**保留在 manifest 与输入中**，只是不进 core 训练比对 |
| 重复序列 | 计算集合内重复 | 0 组 |

**结果：276 条 → core 训练集 270 条**（`inputs/training_set.faa`，131,000 B / `de8e03d348f75b37179bcd1ad0b7351867d1e6bc0622a03b575e254889885496`），排除 6 条，逐条理由见 `inputs/training_set_manifest.tsv`（276 行；每条含 accession、家族、保留/排除、排除原因、序列 SHA-256）。

> **同轮证据审计（§7）顺带发现：** 若把审计到的实验证据计入，Cys 家族**仍只有 2 条**不同记录有实验支持（PhaZ1 与 PhaZ2，见 §7），**低于 gate 的 ≥3**。因此本轮的发现层训练**不能**被读作"gate 变得更接近通过"。

---

## 3. 预注册原文（训练与打分之前写入，时间戳 2026-09-17T06:11:11Z / 14:11:11 CST）

文件：`inputs/falsification_preregistration.{tsv,json,md}`（TSV 4 行 G1–G4 + JSON + 逐字 md）。**注册先于训练**：实现脚本 fail-closed —— 若 `results/cys_discovery.hmm`、`results/cys_discovery_scores.tsv` 或 `results/gate_evaluation.json` 已存在，脚本拒绝写入（`RuntimeError`），并有测试钉住。

**钉死的常量：**

- **发现阈值：E-value < 1e-5**；
- **HMMER 比较数 `-Z 109087`**（参考集、混淆集、候选集**同一 E-value 尺度**；否则 724 条集合的 E-value 会因 Z 更小而系统性偏小约 150 倍）；
- **阈值在失败后不得下调**（`threshold_lowering_permitted=false`、`threshold_change_permitted_after_failure=false`）；
- **校准 gate 未变**（`calibration_gate_unchanged=true`，文本见 §1 表末行）；
- 运行前 `pipeline/config/formal_scan_models.tsv` SHA-256 = `8a1c0508…b7a3`。

**四条闸门（逐字要点）：**

| ID | 陈述要点 | 指标 | 阈值 | 总体 |
|---|---|---|---|---|
| **G1** 特异性 | 发现阈值下不得命中胞外参考面板，也不得命中 MCL 解聚酶对照 Q84C08 | `unexplained_hits` | `== 0` | **Q84C08 (1) + 365 条胞外参考**（284 type 1 + 71 type 2 + 6 dPHAMCL + 4 PhaZ7-like）= **366**；29 条 i-nPHASCL-with-lipase-box 若有命中记为 `explained_cross_talk` 并逐条披露，**不计入 G1 失败** |
| **G2** 灵敏度 | 通过完整性过滤的 Cys 参考召回率 ≥ 0.90，**且明确声明是样本内度量** | `recall` | `>= 0.90` | 270 条 core 训练成员 |
| **G3** 分层 | 所有输出标注 `discovery_hmm_uncalibrated`；`formal_scan_models.tsv` 前后逐字节相同 | `sha256_before == sha256_after AND unlabelled_outputs == 0` | 同左 | 本 run 全部输出 |
| **G4** 不产生 family 判定 | 不产生 family 判定、不改任何 `subtype_call`；输出表不得携带 family_call/subtype_call 列 | `subtype_call_rows_changed` | `== 0` | 全部 109,087 条候选 |

---

## 4. 服务器执行（Step 5–6）

**主机：** `grius` / `<SERVER_HOST>`，`<SERVER_USER>`，`${PHB_REMOTE_ROOT}/PHB_gtdb-ds`。SSH 使用 `C:\Windows\System32\OpenSSH\ssh.exe`（Git 自带 ssh 报 `CreateFileMapping` 错误）。

**负载快照（实测）：**

| 时点（UTC） | load average | 备注 |
|---|---|---|
| 训练前 `06:14:39Z` | 0.26 / 0.44 / 1.14 | 80 逻辑核；内存 1007 GB 总 / 992 GB available；2× RTX 4090 均 0% / 1 MiB |
| 训练后 `06:14:43Z` | 0.32 / 0.45 / 1.14 | — |
| 打分前 `06:14:52Z` | 0.30 / 0.44 / 1.14 | — |
| 打分后 `06:19:08Z` | 5.17 / 1.46 / 1.28 | --max 候选扫描期间 |
| 采集 `06:19:36Z` | 3.19 / 1.33 / 1.25 | 磁盘可用 100T |

**工具与命令（全部在服务器 dated `deploy/20260917_phaded_cys_discovery_hmm_01/` 中绑定）：**

| 步骤 | 命令 | 版本 | 线程 | 实测 |
|---|---|---|---|---|
| 比对 | `mafft --auto --thread 40 --anysymbol inputs/training_set.faa > results/cys_discovery_training.aln` | MAFFT **v7.525 (2024/Mar/13)** | 请求 40，**MAFFT 自报 8 thread(s)** | 270 序列 → 932 比对列 |
| 建 HMM | `hmmbuild -n cys_discovery_uncalibrated --cpu 40 --amino results/cys_discovery.hmm results/cys_discovery_training.aln` | HMMER **3.4 (Aug 2023)** | **40 worker threads** | `nseq=270 alen=932 mlen=451 eff_nseq=1.50 re/pos=0.589`；CPU 9.04u / elapsed 0.40s |
| 打分（主） | `hmmsearch --cpu 40 -E 1e-5 -Z 109087 --tblout results/<set>_discovery_threshold.tblout -o logs/<set>_discovery_threshold.txt results/cys_discovery.hmm <set.faa>`（set = candidates / reference / confounders） | HMMER 3.4 | 40 | 候选 7.14s / 参考 1.47s / 混淆 1.85s |
| 打分（补充） | `hmmsearch --cpu 40 --max -E 1e9 --domE 1e9 -Z 109087 --tblout results/<set>_discovery_max.tblout …` | HMMER 3.4 | 40 | 仅用于给**非命中**补数值分数；**命中判定只由主命令决定** |

**产物哈希：** HMM `results/cys_discovery.hmm` = 210,746 B / `70e1b347d8f8d12dd687c51ea907257a71de63d4629cb398f5225d99614f6333`；比对 `results/cys_discovery_training.aln` = 271,045 B / `2084c4ce8de5a84bfaad6befcab2488e92c75e770fd30bbf70f13219c526e459`；全部服务器产物哈希见 `logs/server_artifacts_sha256.txt`，服务器环境快照见 `logs/server_environment_snapshot.json`。

> **⚠️ MAFFT 不可位复现（如实记录）：** 同一条 `--auto --thread 40` 命令在**同一批输入上跑两次得到不同比对** —— 第一次（失败于 hmmbuild 参数错误的那次）`e61612fd876e845b45a2fa027590544a37b77b71469112a62f782ce0c5633314`，第二次（实际使用）`2084c4ce…c526e459`。因此 **HMM 绑定的是 `2084c4ce…` 这张比对**；换机重跑可能得到不同（但等价性未测定）的 HMM。这不是"随机失败"，而是 MAFFT 迭代精修在 8 线程下的非确定性，**必须**在引用本 HMM 时写明。

---

## 5. 打分统计（Step 6–7 实测）

**候选池（109,087 条，全库覆盖，无遗漏行）：**

| 项 | 计数 | 口径 |
|---|---:|---|
| 命中（E < 1e-5，预注册命令） | **31,906**（29.25%） | `discovery_hit=true`，`score_source=pre_registered_threshold_run` |
| 未命中但有数值分数 | 40,909 | 来自 `--max` 补充扫描（`score_source=supplementary_unfiltered_run`） |
| 无论主/补充命令都**没有报告任何比对** | 36,272 | `score_source=not_reported`、`best_E=not_reported`（**不是**"分数很大"，而是"没有可报告的比对"） |
| 合计 | 109,087 | — |

**未分配命中与 Task 6 标记层的重叠：**

| 项 | 计数 |
|---|---:|
| 命中且**未分配**（`candidate_assignment_status=unassigned_PhaDED_like` 或无超家族先验）→ `results/cys_discovery_claims.tsv` | **30,934** |
| 命中但已有分配（不计入 claims） | 972 |
| Task 6 标记层规模（`PF06850 detected ∧ 无 lipase box`，本 run 独立重算） | **29,949**（与 `cys_decision_manifest.json` 的 `cys_type_candidate_evidence_count` 一致） |
| **命中 ∩ 标记层** | **29,924 / 31,906 = 93.79%** |
| **claims ∩ 标记层** | **29,249 / 30,934 = 94.55%** |
| claims 中**不在**标记层 | 1,685 |

**参考面板与混淆集：**

| 集合 | 记录数 | 命中（E < 1e-5） | 备注 |
|---|---:|---:|---|
| 723 条冻结参考 + Q84C08 | 724 | **276** | 命中集**恰好等于** 276 条 Cys 型参考：276/276 全命中；其余 448 条（365 胞外 + 29 with-lipase-box + 52 nPHAMCL + 1 periplasmic + Q84C08）**一条也不命中** |
| x₁ 疏水混淆集（Task 3） | 563 | **443（78.69%）** | ⚠️ **宽度警示**：该发现层不具 x₁ 分辨力；这是"广度"观察，**不是**预注册闸门 |

---

## 6. 预注册闸门结果（`results/gate_evaluation.json`）

| ID | 结果 | 关键实测 |
|---|---|---|
| **G1** 特异性 | ✅ **passed** | 总体 366（365 胞外 + Q84C08）；`unexplained_hits = 0`；lipase-box 超家族交叉 `explained_cross_talk_count = 0`（无命中，因此无需披露清单） |
| **G2** 灵敏度 | ✅ **passed** | `recall = 270/270 = 1.0000`（阈值 0.90）；**`measurement_scope = in_sample_not_held_out`** |
| **G3** 分层 | ✅ **passed** | `formal_scan_models.tsv` 前/后均 `8a1c0508…b7a3`（450 B），逐字节相同；本 run 12 个发现层输出文件**全部**携带标签，`unlabelled_outputs = []` |
| **G4** 不产生 family 判定 | ✅ **passed** | `registered_metric = subtype_call_rows_changed = 0`；`call_columns = []`；边界声明列 `family_call_made` / `new_family_call_made` / `subtype_call_impact` 为"声明无判定"列，**同时**原样列入 `literal_token_columns` 以备审阅 |

**`status = all_gates_passed`、`falsified_gate_ids = []`、`threshold_lowered_after_failure = false`、`calibration_gate_unchanged = true`、`calibration_gate_met = false`。**

### 6.1 G1 的"零命中"是**实测边距**，不是空查询（补充度量）

`results/g1_margin_analysis.json`：在**宽松的 `--max` 扫描**中，366 条 G1 总体里 **209 条有数值分数**，其余 **157 条（含 Q84C08）在任何配置下都没有可报告的比对**。209 条里有数值分数者的**最小 E-value = 230（bit score 2.7）**，即比 1e-5 阈值高 **约 7.4 个数量级**；胞外参考中最接近阈值的是 `CAL93187.1`（E = 230）。**诚实限定：** 对那 157 条，"零命中"的依据是"**没有可报告的比对**"而不是一个数值 E-value；本 run 不把它写成"E-value 极大"。

### 6.2 G2 的诚实限定（最重要的一条）

270 条分母**就是**训练比对的成员，因此 1.00 是**自洽性**度量。冻结台账中**不存在** held-out 的 Cys 阳性（276 条里只有 1 条 `experimental_positive`，且它也在训练集内），所以**本 run 无法测量泛化**。任何引用都必须写成"**在训练成员上 270/270 通过发现阈值**"，**不得**写成敏感度、特异度或家族分辨力。

---

## 7. 证据审计（Step 2）：O87189 / Q7WT48 / Q7WT49

产物：`results/reference_evidence_audit.tsv`（3 行）、`results/cys_evidence_statistics_per_run.json`；原始检索证据（UniProt .txt ×3、PubMed XML ×2、GenBank .gp ×3 + 检索清单）保存于 `inputs/external_evidence/`，**检索日期 2026-09-17，HTTP 全部 200**，逐文件 SHA-256 见该目录的 `retrieval_manifest.txt` 与 `inputs/input_manifest.tsv`。

| accession | 是什么 | 实验证据 | 与冻结台账/既有条目的关系 |
|---|---|---|---|
| **O87189** | C. necator H16 胞内 PHB 解聚酶 **PhaZ1**（`phaZre`，EMBL `BAA33394.1`，AB017612） | **experimental_positive**：PMID **11114905**（Saegusa 2001）——克隆基因在大肠杆菌中表达的粗提物可水解**无定形** PHB 颗粒并释放寡聚/单体 3HB；**不水解结晶 PHB**；产物在 PHB 颗粒上定位；`phaZ` 破坏株在富营养培养基中 PHB 含量高于野生型；且**推导序列缺少经典 lipase box Gly-X-Ser-X-Gly** | **序列与台账 `DED_hfam_65_0001`（`CAJ92291.1`）逐字节相同**（同 419 aa）→ 它是台账**唯一**那条 `experimental_positive` 的**别名**，不是独立记录 |
| **Q7WT49** | H16 胞内 PHB 解聚酶 **PhaZ2**（`phaZ2`，EMBL `AAP74580.1`，AF549808） | **experimental_positive（体内缺失株表型）**：PMID **12813072**（York 2003）——ΔphaZ1/ΔphaZ2/ΔphaZ3 及双、三突变体在 PHB 利用条件下分析，"PhaZ2 is thus suggested to be an intracellular depolymerase"。**注意：不是体外酶活测定** | 序列与台账 **`DED_hfam_61_0029`（`AAP74580.1`）逐字节相同**；与 `DED_hfam_61_0030`（`CAJ93939.1`）为同一蛋白（后者是其 **C 端 404 aa 截短表示**，差 15 个 N 端残基） |
| **Q7WT48** | H16 **PhaZ3**（`phaZ3`，EMBL `AAP74581.1`，AF549809） | **不是** experimental_positive：PMID 12813072 原文为 **"The role of PhaZ3 remains to be established."**（基因真实存在、注释为胞内，但**功能未确立**） | 序列与台账 **`DED_hfam_66_0001`（`AAP74581.1`）逐字节相同**（另有 `DED_hfam_66_0002` = `CAJ95805.1` 同基因基因组注释版） |

**结论与对既有叙述的更正（本 run 只记录，不改冻结台账）：**

1. Task 5 面板把它们记作 `intracellular_non_ephaz_negative`（"胞内定位**阴性**对照"）—— 该标签是**定位轴**的用法，但**三条全部已是冻结台账成员**，其中 **O87189 就是台账唯一实验阳性的别名**、**Q7WT49 有实验支持的解聚酶角色**。因此把它们读作"实验阴性"**与证据不符**；正确的表述是"**胞内定位对照**"，且其 `family_binding_status=unresolved`（既有结论不变：**不能**计入 gate 的 family-resolved 阴性槽位）。
2. **本轮计入 Cys 家族证据统计（仅本 run）：** 有实验支持的不同记录 = **2 条**（`DED_hfam_65_0001` PhaZ1、`DED_hfam_61_0029` PhaZ2），`DED_hfam_66_0001`（PhaZ3）功能未确立。**gate 要求 ≥3 → 仍未达成**（`gate_positive_requirement_met = false`）。冻结台账哈希前后未变。
3. **对本轮训练集的影响 = 0**：训练集严格来自冻结台账的 276 条（O87189/Q7WT49/Q7WT48 的序列已分别以 `CAJ92291.1`/`AAP74580.1`/`AAP74581.1` 形式在台账内），**没有**额外加入任何序列，也没有因此调整排除规则。

---

## 8. 失败与缺陷记录（如实保留，不掩盖）

**(A) 服务器输入被截断（两次），已修复并加 fail-closed 护栏。** 详见 `logs/incident_candidate_fasta_truncation.md`：一条本 agent 的内联远程命令在 Windows PowerShell → `ssh.exe` 传递中**丢失了双引号**，`grep -c "^>" file` 在远端变成 `grep -c ^> file`，其中 `>` 是 **shell 重定向**，把 `inputs/candidate_union.faa` 截断成 2 字节 `0\n`（SHA-256 `9a271f2a…86aa`）。**第一次**在打分启动后被 `03_score.sh` 的"0 records"发现并立即终止任务；**第二次**在复查时复现同一故障。两次都从本地权威副本重新上传并校验，且在 `03_score.sh` 加入**逐文件 SHA-256/字节数/记录数**的 fail-closed 护栏（不符即 exit 2）。**后果：没有任何闸门输入在损坏状态下被评估。**

**(B) 第一轮闸门评估有两个本 agent 自身的实现缺陷（已修，失败证据已留档）。** 修好后由新增测试捕获根因，并可用 `logs/reproduce_pre_fix_gate_evaluation.py` 复现当时结果（写入 `results/gate_evaluation_pre_fix_reproduction.json`，明确标为 `defect_reproduction_not_a_gate_result`）：

1. `parse_tblout` 只以 HMMER 报告的**原始 target 名**为键。参考面板的 FASTA 头是 `>reference_id|accession`，而所有下游表按 **accession** 索引，于是**第一轮的每次 accession 查询都查不到任何东西**：G1 的 "passed" 是**空洞通过**（0 未解释命中只因为 0 次查询命中），G2 报 `recall = 0.0` 的**伪失败**。修复：`parse_tblout` 对含 `|` 的 target 名，**同时**按各字段建别名索引（别名共享同一 entry，不重复计数）。
2. `evaluate_g4` 用子串 `family_call`/`subtype_call` 扫列名，把**声明"没有判定"的边界列**（`family_call_made`、`subtype_call_impact`）也当成判定列，于是第一轮 G4 报**伪失败**。修复：G4 以**预注册的 metric/threshold（`subtype_call_rows_changed == 0`）**判定，同时显式区分 `call_columns` / `boundary_declaration_columns` / `literal_token_columns`，并把字面扫描结果**原样输出**以备审阅。
3. **字段命名更正：** `n_domains` → **`n_reported_rows`**。实测 `hmmsearch --tblout` 每个 target **只输出 1 行**（最佳 domain），把它叫"domain 数"过强（沿用项目对过强字段名 P4/V5 的处置惯例）。

**(C) MAFFT 非位复现**（见 §4 警示框）：两次同命令产生不同比对；本 HMM 绑定 `2084c4ce…`。

---

## 9. 明确不做的事 / 未解决项

**明确不做：** ❌ 不修改 `pipeline/config/`（尤其 `formal_scan_models.tsv`）❌ 不新增 registry 条目 ❌ 不降低校准 gate 或发现阈值 ❌ 不做 `subtype_call` 变更 ❌ 不把命中写成 family 判定或表型 ❌ 不据此删除任何候选 ❌ 不启动 GTDB 全量重扫 ❌ 不删除任何文件 ❌ 不执行 `git commit`/`push`/`clean`/`reset`。

**未解决项（下一轮）：**

1. **该发现层很宽**：x₁ 混淆集命中率 78.69%，claims 里有 1,685 条不在 Task 6 标记层内。要把它当优先级工具，必须先在一个**独立且经验证的阴性集**上测量假阳性率（当前工作区**没有**这样的集合：冻结台账 `experimental_negative = 0`）。
2. **泛化未测**（G2 是样本内）。
3. **Cys 型催化亲核体坐标仍 `pending_reference_annotation`**；本轮 HMM 不含任何残基级锚定，**不得**用于催化三联体论证。
4. **MAFFT 非确定性**使 HMM 不可位复现；若要把该层纳入长期管线，应固定 MAFFT 线程数/种子或记录"等价性容差"。
5. **校准 gate 仍未达成**（Cys 家族实验阳性记录 = 2 < 3）。
6. 上游 gate 计数器只统计绑定到 profile 自身 family 的记录（跨家族轴对照恒为 0）—— **仅报告、不修改**（属其他 agent 的产物）。

---

## 10. 验收与不变量（Step 9 实测）

| 检查 | 结果 |
|---|---|
| 新增聚焦测试 | **38 个**（`test_audit_phaded_cys_reference_evidence` 8 + `test_build_phaded_cys_discovery_training_set` 6 + `test_preregister_phaded_cys_discovery_gates` 6 + `test_parse_phaded_cys_discovery_hmm_scores` 18）；**全部先写失败测试**，红灯日志见 §11。结果：`Ran 38 tests … OK`，exit 0（`logs/step9_focused_tests.log`） |
| 完整测试套件 | `python -m unittest discover -s pipeline/tests` → **`Ran 522 tests in 15.557s` + `OK (skipped=1)`，exit 0**（= 基线 **484** + 本 run 新增 **38**）；日志 `logs/full_test_suite.log` |
| `python -m compileall -q pipeline` | exit 0（`logs/compileall.log`） |
| `git diff --check` | exit 0（无空白错误、无输出；`logs/git_diff_check.log`） |
| run 内产物清单 | `logs/artifact_sha256_manifest.txt`（**96** 条 path/bytes/sha256/mtime；清单文件自身除外）；`logs/verify_artifact_manifest.py` 逐条复验磁盘字节 → **`VERDICT: MATCH`**（mismatch 0、未收录 0） |
| 交付物自洽复验 | `logs/verify_run_outputs.log`：score 表 **109,087 行 / 109,087 唯一 accession**；claims **30,934 行**；审计 **3 行**；训练 manifest **276 行 / 270 条入选**；**边界列违规 0 行**；G1–G4 全 `passed` |
| 本轮作者文件哈希（状态文档、decision record、4 个脚本、4 个测试） | `results/run_integrity_checks.json` 的 `this_run_authored_files` 块 |
| 服务器脚本绑定 | 本地 `logs/server_scripts/` 与服务器 `deploy/20260917_phaded_cys_discovery_hmm_01/` **逐个哈希相同**（`01` `5c2e7f14…`、`02` `b4d3cf9c…`、`03` `57830133…`、`04` `ef249f93…`、`expected_inputs_sha256.tsv` `093b0dbe…`） |
| `pipeline/config/formal_scan_models.tsv` | 450 B / `8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3` → **与运行前一致，未变** |
| 历史 run（`20260901_*`–`20260916_*`） | 扫描 **0** 个文件 mtime ≥ 2026-09-17（`results/run_integrity_checks.json`） |
| 冻结输入 | 台账 / 参考 FASTA / Task 7 表 / Task 2 表 / 候选 FASTA 哈希均记录于 `results/run_integrity_checks.json`，未被本 run 修改 |
| `input_contract.json` | 26 个输入 `status=verified`（path/size/SHA-256），GTDB taxonomy/metadata/tree **三项均 `pending`**（未伪造哈希） |
| `git commit`/`push`/`clean`/`reset` | **未执行**（`git rev-parse HEAD` = `2e44aafe2c91263447f37bac625f97c56ddafc05`，与交接文档基线一致） |

---

## 11. Test-first 证据（红灯日志）

| 日志 | 内容 |
|---|---|
| `logs/step2_test_first_red.log`（33,351 B） | 首轮 31 个测试**全部 error**（实现文件尚不存在：`FileNotFoundError: pipeline/scripts/audit_phaded_cys_reference_evidence.py`）→ `Ran 31 tests … FAILED (errors=31)` |
| `logs/step7_test_first_red2.log` | 别名索引 + G4 边界列缺陷的 5 个新测试 → `failures=4, errors=1` |
| `logs/step7_test_first_red3.log` | `n_reported_rows` 更名 + 混淆集真实记录计的 4 个新测试 → `errors=4` |
| `logs/step2_focused_tests_green_rerun.log`、`logs/step7_focused_tests_green.log` | 对应绿灯 |

**披露：** 首轮红灯日志中出现过一次 `retrieved_http_status` 列名在测试夹具表头里写错（实现正确、夹具漏列），以及一处测试期望排序写错（`sorted` 输出）；两者都是**测试自身**的缺陷，已在同一会话内修正，非实现被放松。

---

*本状态文档为 candidate-only 记录；所有 HMM、profile、domain、motif、定位与序列上下文证据只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。本轮的 HMM 是未校准的发现层模型，不是 registry 条目，也不产生任何 family 判定。*
