# T141 PhaDED 项目全面审查文档（2026-09-21）

> ## ⚠️ 历史文档声明（2026-09-28 追加，**正文一字未改**）
>
> **本文档是 2026-09-21 时点的历史审查记录，其结论已被后续重构部分取代。引用本文任何数字前，请先读 2026-09-28 的两份现行文档：**
>
> | 现行权威 | 路径 |
> |---|---|
> | 设计方案（证据模型） | `docs/superpowers/specs/2026-09-28-phaded-evidence-model-redesign.md` |
> | 实施计划（13 个 Task） | `docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign.md` |
> | 实施状态（v1 权威哈希 + 逐任务记录） | `docs/T141_20260928_phaded_evidence_model_redesign_status.md` |
> | 目标证据运行的授权申请包 | `docs/T141_20260928_phaded_evidence_run_authorization_packets.md` |
>
> **已被取代的口径（本文正文保留原文，不作回写）**：
>
> 1. **二元「高可信」标签**：40,690 / 36,559 / 4,131 是 **v1 冻结口径，仅作 v1 比较输入**。v2 改为每条蛋白唯一 `primary_disposition`（6 类）+ 多 `evidence_flags`，不再用一个二元列混合不同家族的证据标准。
> 2. **`ahsmg` 作为第三种亲核体**：v1 的 `nucleophile_type = cys/ser/ahsmg` 已被取代。**AHSMG 是含 Ser 的 motif class**，v2 拆为 `nucleophile_identity ∈ {ser, cys, unresolved}` × `motif_class`。
> 3. **SignalP 与定位混用**：v1 把 SignalP 预测与「胞内/胞外」混在同一判据里。v2 拆为 `transport_signal_prediction`（预测）与 `localization_evidence`（证据），**SignalP 不再等于定位真相**。
> 4. **SBD 与 type 1/type 2 的绑定**：SBD 不参与 type 1/type 2 的定义，改记 `accessory_domain_architecture`。
> 5. **`model_status = trained` 的单一含义**：v2 引入四层 `model_layer`（`reference_query_only` / `discovery_hmm_uncalibrated` / `sequence_family_hmm_validated` / `calibrated_candidate_model`），并把 `functional_calibration_status` 拆成独立列。**实验阳性不再是构建发现层/序列家族层模型的必要条件**；annotation-only 序列经质量与溯源审计后可训练**发现层** HMM。
> 6. **`hfam_70` 的状态**：本文写「判定 `calibrated_candidate_model`（文档级判定，未 finalize）」。v2 改写为 **`candidate_gate_passed_not_promoted`**，仍**未 finalize、未进 registry、未重扫**。
> 7. **with-lipase「百万假阳性」措辞**：改称「**百万级待判别命中**（million-scale unresolved hits）」，继续 `deferred_structure_review`，**禁删、禁入主候选计数**。
> 8. **生长速率组名**：`degrader / non-degrader` → 「**候选基因携带**」vs「**在规定搜索与质量条件下未检出候选**」；属水平平均差为主估计对象。
>
> **本文仍然权威的部分**：§6 的 11 项错误更正留痕、§8 的 gRodon2 实测数字、§4.2/§4.5 的 v1 冻结计数（作为 v1 比较基线）、§7.2 的六个数字陷阱。
>
> **全局边界（本文原声明继续有效，且被 v2 继承）**：所有 HMM/profile/domain/motif/SignalP/结构/系统发育/分类与生长速率结果，**只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型**。

> **审查对象**：本工作区（`<REPO_ROOT>`，分支 `main`，HEAD `3117287`）自 2026-08 至 2026-09-21 的 PhaDED 主线全部工作。
> **编写依据**：`docs/T141_20260917_project_handoff.md`（交接文档，含 §11 多轮更正）、`docs/T141_20260920_phaded_three_gaps_status.md`（三项缺口 + challenge 验证 + hfam_70 gate 敲定，§④）、`docs/T141_20260920_phaded_grodon_growth_status.md`（生长速率分析）、全部冻结 run 产物与 manifest（本会话逐项实测复核）。
> **全局边界声明（引用本文任何结论前必读）**：本项目所有 HMM、profile、domain、motif、SignalP、结构、系统发育、分类与生长速率预测结果，**只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型**。全文 "高可信度" = 候选证据最强筛选层，非实验确认。

---

## 1. 一句话总览

以 Knoll 2009 PhaDED 架构（8 超家族 / 38 家族）为分类权威，对 GTDB R232（199,923 代表基因组、约 2.92 亿蛋白）完成全库 HMM 扫描与多证据筛选，得到：

- 候选宇宙 **109,087** 条（62,666 基因组）；
- **合并高可信度候选 40,690 条**（池内 36,559 + 池外 strong 通过 4,131，0 重叠）；
- 胞内:胞外 = **76:24**（与文献参考集的 49:51 相反，原因见 §7.3）；
- 校准 gate **尚未解锁**：唯一 `calibrated_candidate_model` 是 `DED_hfam_70`（type 2，六腿全过但未 finalize）；
- 生长速率分析：降解菌相对同属非降解菌预测最大生长速率**略有上移**（Δ=+0.0255 h⁻¹，置换 P=0.004，效应量 r=0.08）；**胞内 vs 胞外无差异**（P=0.48）。

---

## 2. 已完成的工作（按时间线）

### 2.1 分类框架建立（文献层）
- 以 **Knoll 2009**（PhaDED, BMC Bioinformatics 10:89, PMID 19296857）为权威分类：原文收录 **735 序列条目 → 587 个不同蛋白**，分 **8 超家族 + 38 同源家族**；28 条实验验证种子 → 6 个既有超家族 + 2 个新引入（胞内 nPHASCL 有 lipase box、周质）。
- 项目冻结台账 `phaded_reference_ledger.tsv`：**723 条参考**（30 `experimental_positive` + 693 `annotation_only`），8 超家族全部有覆盖、38 家族齐全。
- 分类权威层：**46 profile = 9 trained + 37 reference-only**；8 超家族 = `functional_prior`（`registry_eligible=true`，来自实验种子），38 家族 = `sequence_clustering_2009`（`registry_eligible=false`，禁止进 registry）。

### 2.2 GTDB 全库扫描（scan 13）
- GTDB R232，100 个蛋白分片（约 267 GB / 2.92 亿蛋白），10 模型 × 100 分片 = 1,000 个 hmmsearch 任务（纯 HMMER，无 DIAMOND 预筛——已实测排除历史文档误述）。
- `hits_all.tsv`：6,769,772 命中行 → 6,595,162 唯一蛋白。
- 候选宇宙构建：registry 阈值 + 核心家族仲裁 + 验证门 + tier 阈值 → **109,087 = ePhaZ 70,682 + iPhaZ 38,405**（其余 8 个 registry 模型被家族范围裁剪，属有意设计）。

### 2.3 候选映射与证据层
- 9 个 trained profile 顺序打分 → `phaded_assignment.tsv`（109,087，逐条唯一）。
- 证据层逐条补齐：InterProScan 5.78（108,735 条 supported）、Pfam 架构（`pfam_architecture_01`）、motif 面板（Task 7，109,087 全覆盖）、SignalP 6（定位层）、Foldseek 结构证据（19 条 review/support/conflict）、竞争系统发育审查（42 条）。
- subtype 判定层（`subtype_call`，109,087 行）与 9 个并列证据状态列落版（见 §4.2）。

### 2.4 高可信度筛选与合并（40,690）
- 多证据过滤器：池内 109,087 → 高可信 **36,611** + hold 71,860；核亲体统一（`nucleophile_type`/`nucleophile_family`）。
- 池外 profile 命中层 **35,558**：strong 30,632 / mid 1,503 / weak 3,423；strong 层补全证据后 **4,131** 通过 → 与池内合并 = **40,742**。
- 52 条 SignalP-export Cys 候选 demote（36,611→36,559）→ **合并高可信度 40,690**（本会话 2026-09-20 落地，`runs/20260920_phaded_high_confidence_pool_merge_01`，0 重叠、nucleophile 列保留）。

### 2.5 Cys 型定向工作
- 发现层 HMM（`cys_discovery_uncalibrated`）训练（270 条 Cys 参考，Knoll 2009 方式）与全空间召回：**35,504** 条命中，旧 10 模型 registry 召回率 **99.94%**（35,484/35,504），真新 20 条。
- 分层漏斗审计（本会话，只读）：4,301 条 Cys 命中在 `hits_all.tsv` 却被候选宇宙排除——归因为**有意的家族范围 + E 值 tier 规则**（PhaC 1,970 条 + iPhaZ 4,244 条，非边界 bug）；恢复与否待操作者决策。
- 序列级 Cys 标记：PhaZ1 **Cys183** 锚点（PMID 16233560；C183A/D355A/H388Q 失活、C183S 不失活）→ `[疏水]-C-Q` 模式（−1 疏水来自 Knoll `cysteine-1` 规则）；参考面板 Cys 271/276 vs 非 Cys 22/447（精确 VCQ 层 0.909 vs 0.002）。

### 2.6 校准 gate 进展
- 预注册证伪条件 F1–F6（6 条，4 条数值化全 passed）+ G1–G4（发现层训练闸门全 passed）。
- `DED_hfam_70`（type 2）**六腿全过 → `calibrated_candidate_model`**（唯一）；`hfam_52` 次选（6 篇但含弱证据）。
- `hfam_4`（nPHAMCL）**challenge FAIL**（3/10 假阳性：Moraxella 脂肪酶 5.4e-12、Psychrobacter 脂肪酶 4.1e-11、PcaD γ-内酯酶 6.6e-11）→ 降级 `reference_only_insufficient_panel`，改走结构证据。
- 阳性溯源：30 条台账阳性 → 28 独立 GI → 24 篇一手文献；**严格口径只有 `hfam_52`（6 篇）与 `hfam_70`（7 篇）真正 ≥3 独立阳性**；`hfam_4` 实为"1 真阳性 + 2 条同一基因重复"、`hfam_53` 实为 1。
- **challenge 验证**（5 trained × 10 挑战，`runs/20260920_phaded_heldout_rebuild_01/results/challenge_validation.tsv`）：hfam_52/55/70/8 **10/10 全拒绝**；hfam_4 **3/10 假阳性**（见上）。
- **hfam_4 判别性诊断 + `--pnone` 重建（option b，`runs/20260920_phaded_hfam4_discriminant_01`）**：对齐 5 阳性 + 10 挑战逐列判别——① 催化三联体是**共享的、非判别**（亲核 Ser 与肘部两侧 G 在挑战里 10/10 保守）；② 唯一判别信号是肘部上下文 `G-V-S-W-G` 的 V/W（挑战 0/10 共享）；③ 根因 EFFN=0.378（Dirichlet 先验把 profile 拉向背景）。`--pnone` 重建后 3 假阳性清零但 `eff_nseq=0.00`（阳性 E=0、全挑战 E=999）= **精确匹配过拟合、零召回，假胜利**。**定论：hfam_4 无法做成判别性家族 HMM（序列信号 = 通用 PF00561 折叠），nPHAMCL 改走结构证据；操作者已决定不再修改 HMM**（`hfam4_rebuild_conclusion.md`）。
- **hfam_70 gate 六腿明细**（判定件 `runs/20260920_phaded_heldout_rebuild_01/results/gate_determination_hfam_70.json`）：positive=7 独立 ✅ / heldout=AAB02914.1（E=7.8e-146，max_pairwise_id=0.484）✅ / negative=Q84C08+Q51718 ✅ / challenge=10/10 ✅ / 两个 status=sufficient ✅ → **判定 `calibrated_candidate_model`（文档级判定，未 finalize）**；hfam_52 六腿也全过但含 2 条弱阳性（次选）。
- **澄清（防后续误读）**：本轮**无任何 HMM 实质重建**——challenge 用的 `hfam_70_production.hmm` 与参考 profile 逐字段一致（NSEQ=7、EFFN=0.748535）；hfam_52/55/8 阳性集未变；hfam_4 虽重建但已降级。**故无任何 profile 因本轮工作需要重扫全库**（hfam_70 全库召回早已在 `20260919_phaded_trained_profile_recall_01` 跑完：9 trained profile × 100 shard = 41 分钟，hfam_70 原始命中 39,458 条；单 profile 全库 ≈ 5–10 分钟）。

### 2.7 隔离层
- with-lipase（`DED_hfam_2`）池外命中 **1,206,655** 条 → 归档 deferred 结构验证层（禁删、禁入计数，只经 8YNV Foldseek/结构预测召回）。
- 结构验证试点（8YNV chain A）：抽样 1,001 条 deferred 序列 blastp 实测 **954/987 E<1e-5、最小 E=1.16e-122**；样本可外推性经机器判定（Spearman(index, E 秩)=0.0074）。

**为什么 with-lipase 低可信度、必须走结构确认（2026-09-21 补全）**：

1. **序列信号与普通脂肪酶同源、不可区分**：hfam_2 的判别特征就是「带 lipase box（`G-X-S-X-G`）的胞内 SCL-PHB 解聚酶」，而 lipase box 正是细菌脂肪酶/酯酶共有的通用信号——两者都是 α/β-水解酶折叠 + 同样的肘部基序，**一维序列层面没有可判别的家族特异指纹**。任何序列 HMM 只有两个结局：要么先验主导成「近背景」把脂肪酶全部混进来（海量假阳性），要么收紧成精确匹配（零召回，见 hfam_4 `--pnone` 教训）。这是「α/β-水解酶折叠序列不可分」问题在 SCL 侧（胞内）的对应物，与 nPHAMCL（hfam_4，mcl 侧）同病。
2. **实测噪声量级**：池外命中 **1,206,655** 条（2026-09-19 实测，规则后口径）——全部 46 个 profile 里最大的假阳性源，直接证明「lipase box」这一信号本身无特异性。
3. **实验阳性只有 1 条、无独立复核**：仅 **EAO52570.1** 一条实验阳性，出自 Knoll 2009 单一出处（PMID 19296857）；台账 29 条参考里 28 条是 `annotation_only`。远低于 gate 的「≥3 独立阳性」，也无跨属、跨论文的活性表征。
4. **高可信度 0 通过**：池内 13 条全部 `common_criteria` hold、池外 0 通过 → **0 条进入 40,690**（唯一一个 0 通过的胞内超家族）。
5. **为何「结构」而非「序列」能救它**：解聚酶与脂肪酶的特异性差异在**底物口袋的立体化学（3D 几何）**，不在氨基酸序列。因此只能靠结构证据区分——Foldseek 比对 PDB 8YNV（试点已跑）或结构预测（ESMFold/AlphaFold）后比对底物口袋，从 1.2M 噪声里召回「真带 lipase box 的解聚酶」，通过者再走证据层与校准 gate。**该层禁删、禁入计数；全量漏斗（E<1e-50/1e-30）未授权执行（§9 #6）。**

### 2.8 生长速率分析（本会话 2026-09-20/21）
- 用服务器 gRodon2 工具链，对 40,690 映射的 **30,623 基因组**做预测最大生长速率。
- 同属平衡设计：**8,882 manifest（4,441 对）→ 预测 8,828 ok + 54 failed → 平衡分析 8,776（4,388 对），1,216 属**。
- 复用旧工作区 713 条预测；26,116（85.3%）降解基因组因属内无对照被排除。
- 结果见 §8；完整产物 `runs/20260920_phaded_grodon_growth_01/`。

### 2.9 汇报与图件
- `9.21汇报/gtdb分析/`：管线图集（前序会话，5 图）。
- `9.21汇报/降解速率/`（本会话）：`figure_1_growth_comparison_grodon2.{svg,pdf,tiff,png}` + source_data 6 表 + 数据表 8 张 + qa 3 件 + README + 状态文档。

---

## 3. 数据链与产物清单（审计用）

| 阶段 | 输入 → 输出 | 关键数字 | 权威产物 |
|---|---|---|---|
| 文献冻结 | Knoll 2009 → 台账 | 723 参考 / 30 阳性 | `runs/20260909_phaded_reference_evidence_01/inputs/phaded_reference_ledger.tsv` |
| 权威层 | 台账 → 46 profile | 9 trained + 37 ref-only | `runs/20260917_phaded_classification_authority_01/` |
| 全库扫描 | GTDB R232 → hits_all | 6.77M 行 / 6.60M 蛋白 | 服务器 `runs/20260901_formal_frozen_scan_13/` |
| 候选宇宙 | hits_all → 分层 | 109,087 = ePhaZ 70,682 + iPhaZ 38,405 | `runs/20260905_ephaz_iphaz_layering_01/results/layers_retry_02/protein_layers.tsv` |
| 映射+证据 | 9 profile 打分 + 5 证据源 | 109,087 全对账 | `runs/20260916_phaded_full_library_signalp_01/results/phaded_full_library_localization_subtype_matrix.tsv` |
| 高可信 | 过滤器 | 36,611 + hold 71,860 | `runs/20260919_phaded_nucleophile_unify_01/results/high_confidence_candidates.tsv` |
| 池外 | profile 命中 + 证据补全 | 35,558 → strong 4,131 通过 | `runs/20260919_phaded_pool_external_evidence_01/` |
| demote | SignalP-export 52 条 | 36,611→36,559 | `runs/20260920_phaded_three_gaps_01/results/cys_demotion/` |
| **合并高可信** | 池内 + 池外 | **40,690（0 重叠）** | `runs/20260920_phaded_high_confidence_pool_merge_01/results/merged_high_confidence_candidates.tsv` |
| 发现层 | Cys 发现 HMM 全空间 | 35,504 命中（99.94% 已覆盖） | `runs/20260917_phaded_cys_targeted_recall_01/` |
| 生长速率 | 40,690→30,623 基因组 + gRodon2 | 8,776 平衡 / 1,216 属 | `runs/20260920_phaded_grodon_growth_01/` |

---

## 4. 亚型分层体系（当前完整结构）

### 4.1 分层全景

```
GTDB R232（199,923 基因组）→ 2.92 亿蛋白
 └ scan 13：10 模型 hmmsearch → hits_all.tsv（E≤1e-5，659 万唯一蛋白）
    └【候选宇宙】109,087 = ePhaZ 70,682 + iPhaZ 38,405
       ├【分类权威】8 超家族（可进 registry）+ 38 家族（禁进 registry）
       ├【subtype 判定】7 标签 + 9 证据状态列
       ├【高可信筛选】36,559（池内）+ 4,131（池外 strong）→ 合并 40,690
       │   └【亲核体】cys 29,974 / ser 10,634 / ahsmg 82
       └ hold 71,860（localization 48,488 / architecture 19,562 / common 2,538 / nucleophile 1,272）
 ├【池外 score 层】35,558（strong 30,632 / mid 1,503 / weak 3,423；weak/mid 不升级）
 ├【deferred 结构层】with-lipase 1,206,655（禁删、禁入计数）
 └【发现层 HMM】cys_discovery_uncalibrated（只召回、不筛选、不进 registry）
```

### 4.2 八超家族（文献权威 + 项目台账 + 候选三对照）

| 超家族 | 定位 | 亲核体 | 文献 587 | 台账 723 | 高可信 40,690 |
|---|---|---:|---:|---:|---:|
| intracellular nPHASCL **without lipase box**（Cys 型） | 胞内/天然 SCL | **Cys**（无 GxSxG） | 224（38.2%） | 276（38.2%） | **29,974（73.7%）** |
| extracellular dPHASCL **type 1** | 胞外/变性 SCL | Ser（oxyanion 在 box N 端） | 234（39.9%） | 284（39.3%） | 5,018（12.3%） |
| extracellular dPHASCL type 2 | 胞外/变性 SCL | Ser（oxyanion 在三联体 C 端） | 未单列 | 71（9.8%） | 4,348（10.7%） |
| intracellular nPHAMCL | 胞内/天然 MCL | Ser（带 lid） | 未单列 | 52（7.2%） | 987（2.4%） |
| intracellular nPHASCL **with lipase box** | 胞内/天然 SCL | Ser（x₁=Trp） | 未单列 | 29（4.0%） | **0** |
| extracellular dPHAMCL | 胞外/变性 MCL | Ser（无 SBD） | 未单列 | 6（0.8%） | 280（0.7%） |
| extracellular native-SCL/**PhaZ7-like** | 胞外/天然 SCL | Ser（**AHSMG** 肘部） | 未单列 | 4（0.6%） | 82（0.2%） |
| periplasmic PHA depolymerases | 周质 | Ser（type 2 催化域） | 1（0.2%） | 1（0.1%） | 1（0.0%） |
| **合计** | | | **587** | **723** | **40,690** |

**定位聚合**：文献 49:51（胞内:胞外）→ 台账 49:51 → **高可信 76:24**。

### 4.3 subtype 判定层（109,087，7 标签）

| subtype_call | 条数 |
|---|---:|
| `dPHASCL1_like_candidate` | 67,342 |
| `unassigned_PhaDED_like` | 31,632 |
| `ambiguous_family_within_superfamily` | 9,631 |
| `dPHASCL2_like_candidate` | 356 |
| `hold_architecture_conflict` | 83 |
| `ambiguous_superfamily` | 41 |
| `hold_gene_model_or_structure` | 2 |

### 4.4 九个并列证据状态列（0916 SignalP manifest 实测）

| 维度 | 主要取值 |
|---|---|
| `subtype_confidence` | moderate 67,651 / unresolved 41,304 / hold 85 / low 47 |
| `profile_evidence_status` | trained_hit 67,749 / ambiguous_family 9,664 / unassigned 31,632 / ambiguous_superfamily 42 |
| `domain_evidence_status` | supported_partial 66,632 / generic_interpro 39,999 / explicit_interpro 2,232 / conflict 83 / not_tested 140 / pending 1 |
| `motif_evidence_status` | partial_catalytic_pattern_panel 76,929 / conflict_relative_position 373 / not_tested 31,785 |
| `localization_evidence_status` | export_signal_supported 45,751 / localization_unassigned 62,985 / tool_input_excluded 351 |
| `structure_evidence_status_v2` | **not_tested 109,068** / review 11 / support 5 / conflict 3 |
| `phylogeny_evidence_status_v2` | **not_tested 109,045** / review 22 / support 19 / unresolved 1 |
| `interpro_status` | supported 108,735 / tool_input_excluded 351 / not_reported 1 |

### 4.5 亲核体层（40,690）

| nucleophile_type | 条数 | 依据 |
|---|---:|---|
| `cys` | 29,974 | 无 GxSxG + PF06850 + 强 Cys 家族 profile（`inferred_cys_not_verified`，**未逐条验证**） |
| `ser` | 10,634 | GxSxG（`motif_level`） |
| `ahsmg` | 82 | PhaZ7 型 AHSMG 肘部 |

> 口径提醒：`cys` 数与超家族 "without lipase box" 数完全一致（29,974）；`nucleophile_conflict` 1,272 条被 hold（Cys 带 GxSxG 233 + Ser 缺 GxSxG 1,039）。

### 4.6 校准状态（哪些层"已校准"）

| profile | 状态 | 依据 |
|---|---|---|
| `DED_hfam_70`（type 2） | **calibrated_candidate_model**（唯一） | 6/6 gate 腿；但 **finalize 未执行、未进 registry、未重扫** |
| `DED_hfam_52`（type 1） | 次选 | 6 篇（约 4 篇活性表征 + 2 弱证据） |
| `DED_hfam_55` / `DED_hfam_8` | 不达标 | 各仅 2 篇独立 |
| `DED_hfam_4`（nPHAMCL） | **降级** reference_only_insufficient_panel | challenge 3/10 假阳性；EFFN 0.378 |
| 其余 37 个 | reference_only | 阳性 <3 或 family-resolved 阴性缺失 |

**trained 家族六腿合成 gate 状态（2026-09-20 实测，`three_gaps_status.md` §④）**：

| 家族 | positive≥3 | heldout≥1 | negative≥1 | challenge | 两个 status | gate |
|---|---|---|---|---|---|---|
| hfam_70（type-2） | ✅ 7 | ✅ AAB02914.1（id 0.484） | ✅ Q84C08+Q51718 | ✅ 10/10 | sufficient | ✅ **六腿全过** |
| hfam_52（type-1） | ✅ 6（~4 强+2 弱） | ✅ AAB40611.1（id 0.304） | ✅ Q84C08+Q51718 | ✅ 10/10 | sufficient | ✅ 全过（弱 2 条） |
| hfam_55（type-1） | ❌ 2 | ✅ P52090.1（id 0.701） | ✅ | ✅ 10/10 | sufficient | ❌ 差 1 阳性 |
| hfam_8（dPHAMCL） | ❌ 2 | ✅ AAQ72538.1（id 0.840） | ✅ | ✅ 10/10 | sufficient | ❌ 差 1 阳性 |
| hfam_4（nPHAMCL） | ❌ 1 | ❌ leakage 0.901 | ✅ PcaD | ❌ 3 假阳性 | 负性 sufficient / challenge **not sufficient** | ❌ 多腿失败 |

### 4.7 架构亚型分层规则（判别键，2026-09-21 补充）

| 超家族 | 定位 | 底物 | 催化域/标记 | 亲核体 | lipase box | 配件域(SBD) | 判别要点 |
|---|---|---|---|---|---|---|---|
| type 1 | 胞外(SP/LIPO/TAT) | SCL | PF10503 type-1 | Ser | ✅ | FN3/Ig-like/TSP3/CHB | 唯一带 SBD 配件域 |
| type 2 | 胞外 | SCL | PF10503 type-2 | Ser | ✅ | 无 | 与 type-1 差在配件域 |
| dPHAMCL | 胞外 | mcl | mcl 型（oxyanion 未定） | Ser | ✅ | 无 | 底物 mcl 区分 |
| PhaZ7-like | 胞外 | native-SCL | AHSMG | AHSMG | ❌（AHSMG 替代） | — | 唯一 AHSMG 肘部 |
| nPHAMCL | 胞内(OTHER) | mcl | PF00561（α/β-水解酶） | Ser | ✅ | — | 序列 = 通用 α/β-水解酶 → challenge FAIL |
| with lipase box | 胞内 | SCL | 脂酶型 | Ser | ✅ | — | 与普通脂肪酶不可分 → 0 通过、deferred 结构层 |
| without lipase box（Cys 型） | 胞内 | SCL | PF06850 | **Cys** | ❌ | — | 唯一 Cys 亲核体 |
| periplasmic | 周质 | native PHB | type-2 类 | Ser | ✅ | — | 1 条实验阳性 |

**判别键**：定位（SignalP）→ 底物（SCL/mcl/native）→ 催化域（PF10503/PF00561/PF06850/AHSMG）→ 亲核体（Ser/Cys）→ lipase box 有无 → 配件域（type-1 SBD）。

---

## 5. 证据完备性评估（什么齐了、什么缺了）

### 5.1 已完备

- **motif 层**：109,087 全覆盖（Task 7）。
- **InterPro 层**：108,735 条有命中。
- **定位层（候选）**：45,751 export_signal_supported；29,469 条"胞内"补测 SignalP 实测 OTHER（99.82%）。
- **高可信层**：40,690 全部通过多证据（profile + 架构 + 定位 + 亲核体一致性）。

### 5.2 明确缺口（按严重度排序）

| # | 缺口 | 规模 | 影响 |
|---|---|---|---|
| 1 | **校准 gate 未解锁**：37/37 family-resolved 阴性全 `missing`；37 个 reference-only 无 HMM | 37 家族 | 除 hfam_70 外没有任何 profile 可声称"校准过" |
| 2 | **结构层几乎空白** | 109,068/109,087 = 99.98% `not_tested` | 候选的架构一致性主要靠序列级证据；nPHAMCL 类无法序列层判别，只能走结构证据（8YNV Foldseek/预测，未授权执行） |
| 3 | **系统发育层几乎空白** | 109,045 `not_tested` | 竞争审查只覆盖 42 条；全库候选的进化关系未刻画 |
| 4 | **参考层自身未经定位/架构检验**：284 条 type-1 参考 `domain_architecture` 全 `pending`；整个 723 台账 **0 行提及 signal**；参考层**从未跑过 SignalP** | 723 参考 | "extracellular" 标签是**定义性**的（继承自 DED 超家族名，非实测）→ 与候选层的信号肽门**证据标准不对称**（见 §5.3） |
| 5 | **Cys 催化残基未逐条验证**：`inferred_cys_not_verified` | 28,534（池内） | Cys 判型是推断（无 GxSxG + PF06850 + 家族 profile），未排除"无 GxSxG 的 Ser 型" |
| 6 | **结构断点**：候选 `candidate_ded_alignment_members = 0` | 全部候选 | DED 比对坐标无法传递到候选序列（His/Asp 不可靠） |
| 7 | 池外 mid/weak **4,926** 条未做 SignalP/InterPro（决策：不升级、保留 score_tier） | 4,926 | 若未来要升级需单独授权补算 |
| 8 | 52 条 SignalP-export Cys 已 demote，但**未接入过滤器**（计数变更已由操作者确认执行；`filter_wired=false` 状态未变） | 52 | 未来重跑过滤器须含此规则 |

### 5.3 证据标准不对称（本会话新发现，审查重点）

- 候选层 type-1 门：**必须预测到分泌信号**（否则 hold）；实测 48,038/67,720（71%）因此被 hold。
- 参考层 type-1：**从未被这条标准检验**（无 SignalP、架构全 pending）。
- 文献层：Knoll 2009 自录反例——**gi:74267419**（R. eutropha H16）被归入胞外 type-1 家族、"highly active against artificial amorphous PHB granules, and is **lacking a signal peptide**, a linker domain, and a substrate binding domain"；gi:194292521（C. taiwanensis）GenBank 注胞内、按序列相似性判胞外。
- 被 hold 的 48,038 条中，**11,610 条**同时具备 `type1_verified` 几何 + trained profile 命中 + InterPro 支持 → 仅因缺分泌信号被排除；其中既有"真 type-1 但信号肽未检出/序列缺失"（Knoll 模式）也有"折叠误判"，**当前证据无法裁决**。
- **建议的裁决动作（未执行、需授权）**：① 给 723 参考层补跑 SignalP，标定"已知胞外参考的预测信号肽检出率"；② 给 11,610 条核查 type-1 配件域（FN3/Ig-like/TSP3/CHB，独立于信号肽的胞外证据；`type1_accessory_domain_layer.tsv` 已有 5,018 条）。

### 5.4 低可信度家族归因（本会话实测）

| 归因 | 家族数 | 覆盖参考 | 占比 | 代表 |
|---|---:|---:|---:|---|
| **① 实验阳性不足**（有序列、无锚点） | **21** | **570** | **79%** | hfam_57（83 序列/0 阳性）、hfam_63（64/0）、hfam_65（73/1） |
| ② 种子同质 → profile 统计退化（"序列分不开"） | 1（hfam_4） | — | — | NSEQ=3、EFFN=**0.378**、challenge 3/10 假阳性 |
| ③ 序列量本身太少 | 12 | 16 | 2% | 周质 hfam_3（全库 1 条）等 11 个 1–2 条家族 |
| ④ 比对列不可传递（结构原因） | 1（hfam_70） | — | — | 列保守度 0.0769 < 0.5（type-2 环状置换） |

**结论**：低可信度主因是**实验阳性证据缺口（79%）**，不是"序列互相分不开"、也不是"序列太少"；真正的序列层不可分目前只有 hfam_4 一个实测案例（with-lipase hfam_2 的 0 通过同理——lipase box 与普通脂肪酶信号同源、1,206,655 池外噪声，已 deferred 结构层）。

**hfam_4 判别性诊断补充（2026-09-20，option b 实测）**：阳性肘部 `G-V-S-W-G` 中 V/W 是唯一 0/10 挑战共享的判别位点；催化三联体（Ser + 两侧 G）10/10 挑战共享 = 通用 α/β-水解酶机器、非判别。`--pnone`（去先验）重建虽把 3 假阳性清零，但 `eff_nseq=0.00` = 对 5 条训练序列的精确匹配（阳性 E=0、挑战全 E=999），**零召回、假胜利**——因此 hfam_4 是"序列层不可判别"的实测定论，nPHAMCL（987 条）只能走结构证据。
**trained profile 的 EFFN 体检**：9 个中 8 个 <1（hfam_4 0.378 / hfam_8 0.396 / hfam_55 0.495 / hfam_70 0.749 / hfam_52 0.865 / superfamily type2 0.749 / superfamily dPHAMCL 0.396 / superfamily nPHAMCL 0.378），**唯一 EFFN>1 的是 superfamily type 1（13 种子，1.685）**——即 8/9 的 profile 是先验主导的，即便通过 challenge 者也如此。

---

## 6. 已发现并更正的错误（过程记录，审查时按"是否留痕"检查）

| # | 错误 | 更正状态 |
|---|---|---|
| 1 | PF06850 被当作"胞外 SBD 标记" | 已更正为**胞内 Cys 型 C 端标记**（IPR009656 "degrades PHB granules"；DEFECT_D + three_gaps 任务③）；type-1 催化域 = PF10503（IPR010126，GO extracellular） |
| 2 | hfam_4 训练集"3 阳性" | 溯源为"1 真阳性（P26495）+ 2 条同一基因重复" |
| 3 | hfam_53 "2 阳性" | 溯源为同一 P. stutzeri 基因（实为 1） |
| 4 | 结构验证试点把 **bitscore 当 E 值**（两次复发：初稿 + Fork 3 脚本 `r[4]`） | 已更正（954/987 E<1e-5、最小 1.16e-122），列映射由机器钉死 |
| 5 | "Foldseek 未安装" | 漏查 `tools/`，实为可用（自测 TM 1.00） |
| 6 | patatin "无 GxSxG" | 已更正（patatin 有自己的 GxSxG，Liu 2015 原文 G-X-S47-X-G） |
| 7 | 交接文档 summary 4 处被否决陈述 | amendment A1/A2/A4/A5 已更正 |
| 8 | demote 修订版（36,559）丢失 nucleophile 两列 | 本会话合并时改用 nucleophile_unify 17 列源重建并保留 |
| 9 | 预测表 12 行重启重复 | `dedupe_predictions.py` 去重（8,894→8,882），原始表保留 |
| 10 | gRodon R 多线程失控（data.table 默认并行，load 达 74/80） | 显式 OMP/OPENBLAS/MKL/R_DATATABLE=1 封顶 |
| 11 | 图件审计自身 3 处实现缺陷（缺 text-fill-edge 类、axvline 错坐标变换、直线只测顶点）+ 缺 clipped-artists 检查 | 全部补进审计并经"问题几何回归证明"（详见生长速率状态文档 §9.1 与 汇报 README） |

---

## 7. 关键比例与引用陷阱

### 7.1 比例（引用时写明分母）

| 指标 | 值 |
|---|---|
| 候选宇宙 | 109,087（胞内 iPhaZ 38,405 / 胞外 ePhaZ 70,682）；池基因组 62,666 |
| 合并高可信 | **40,690**（胞内 **76.1%** / 胞外 23.9% / 周质 0.0%） |
| 高可信通过率（池内） | Cys 92.5% / type2 79.5% / dPHAMCL 34.3% / nPHAMCL 14.2% / **type1 7.2%** / with-lipase 0% |
| 文献参考（Knoll 587） | 胞内 48.8% / 胞外 51.0%（type1 39.9% + Cys 38.2%） |
| 发现层 Cys 召回 | 35,504；旧 registry 覆盖 99.94% |

### 7.2 六个必须避开的数字陷阱

1. **8,894 vs 8,882**：生长预测原始表含 12 行重启重复；规范表 `_dedup.tsv`（8,882）。
2. **8,828 ≠ 8,776**：预测 ok 8,828，其中 2 条无有效生长率、50 条因属内平衡剔除。
3. **30,026 vs 29,974**：52 条 demote 前的快照口径（40,742）与现口径（40,690）之别。
4. **SignalP `pending` 是陈旧标签**：29,469 行已实测 OTHER（99.82%）。
5. **with-lipase 用 1,206,655**（trained 覆盖 discovery 规则后），不用 1,238,812（规则前桶）。
6. **文献 vs 台账 vs 候选不是同一分母**：587 蛋白 / 723 记录 / 109,087 候选 / 40,690 高可信，引用比例必须写明层级。

### 7.3 为什么 76:24（本会话结论，须双解释）

胞内:胞外从文献的 49:51 翻转为候选层的 76:24，原因是**判据不对称 + 研究偏向叠加**，不能读作纯生物学比例：
- 判据面：type-1 需要分泌信号 + 几何 + 架构四重证据（71% 卡在分泌信号），Cys 型判据少且松（无 GxSxG + PF06850）→ 通过率 7.2% vs 92.5%。
- 文献面：PhaDED 是策展库，胞外酶易于测活、工业关注高 → 49:51 反映"研究过什么"，不是基因组真实比例。

### 7.4 gate 的「≥3 独立阳性」依据（2026-09-21 澄清）

- **字面出处**：`calibrate_ephaz_subtype_models.py` L109 `len(train) >= 3 and len(genera) >= 3`（承袭 ePhaZ 校准，`runs/20260909_ephaz_core_calibration_01/results/calibration_status.md`）→ PhaDED gate `positive_count >= 3`（`finalize_phaded_reference_panel_calibration.py` L92，fail-closed）。
- **三层依据**：① **HMM 位置特异性谱的最小多样性**——1 条 = 精确匹配零泛化、2 条 = 无法区分"家族保守 vs 个体特异"、≥3 条才形成有意义 profile；② **跨属独立性（防同源泄漏）**——`len(genera) >= 3` 要求 ≥3 属，同一基因簇/同源物（~90% 同一）是冗余进化样本，不重复计数；③ **证据三角化**——≥3 条不同论文/基因/生物的独立实验阳性交叉验证活性。
- **诚实边界**：「3」是**治理启发式**（三角化 + 跨属多样性的最小惯例），项目无「为什么是 3 不是 2/4」的统计推导；其可辩护性靠 ≥3 属硬要求 + held-out/negative/challenge 闭环。
- **实证反例**：hfam_4「3 阳性」实为「1 真阳性 + 2 同基因重复」→ held-out max_id=0.901 leakage_blocked + challenge 3 假阳性 → 降级。这正是该规则存在的直接动机。

---

## 8. gRodon2 生长速率分析（本会话新增，`runs/20260920_phaded_grodon_growth_01`）

### 8.1 设计与账本
- 输入：40,690 → 30,623 降解基因组（纯胞内 22,160 / 纯胞外 5,644 / 兼有 2,819）。
- 同属 1:1 平衡（seed 42；对照 = 同属内不在 62,666 候选池的 GTDB 基因组）；26,116（85.3%）无同属对照被排除。
- 8,882 manifest → 预测 8,828 ok + 54 failed（0.61%，全部 too_few_ribosomal_hits）→ 平衡 8,776（4,388 对）/ 1,216 属。
- 复用旧工作区 713 条预测；R 侧单线程；70 线程为操作者临时授权（已留证）。

### 8.2 结果
| 问题 | 结论 | 关键数 |
|---|---|---|
| ① 降解 vs 同属非降解 | **略有上移，方向一致性弱** | 属水平平均 Δ=**+0.0255 h⁻¹**（CI +0.0065~+0.0446）；Wilcoxon P=0.0033；属内分层置换 P=0.0040；**但** 中位 Δ CI 含 0、符号检验 P=0.115、效应量 r=0.084 |
| ② 胞内 vs 胞外 | **无差异** | Δ=−0.0135 h⁻¹（CI −0.0538~+0.0256）；Mann-Whitney P=0.92；置换 P=0.48 |
| ③ 6 超家族分层 | **全部不显著**（BH 校正后） | q ≥ 0.34；名义最强 type1（+0.0436，P=0.082）与 Cys（+0.0096，P=0.113） |

### 8.3 边界
- gRodon2 = 密码子使用偏好的**预测最大生长速率**，非实测；标签 candidate-only。
- 结论只适用于"降解/非降解混合存在的属"（85.3% 排除的固有覆盖限制）。
- 与旧工作区（旧定义、899 属、无差异 P=0.459）**不可互证**（定义/覆盖不同、集合仅重叠 717 条）。

---

## 9. 待办与需授权项清单

| # | 事项 | 状态 | 需授权 |
|---|---|---|---|
| 1 | 723 参考层补跑 SignalP（标定信号肽门严格度） | 未执行（本会话建议） | ✅ |
| 2 | 11,610 条 SP-less type1 候选配件域核查 | 未执行（本会话建议） | ✅ |
| 3 | 4,301 条 Cys 漏斗丢弃是否恢复 | 待决策（操作者暂缓） | ✅ |
| 4 | `hfam_70` finalize（`finalize_phaded_reference_panel_calibration.py`） | 未执行 | ✅ |
| 5 | 全库重扫（任何新 profile） | 未执行、不需要（99.94% 召回） | ✅ |
| 6 | with-lipase 1.2M deferred 结构验证（Foldseek 漏斗 E<1e-50/1e-30） | 未执行 | ✅ |
| 7 | ESMFold 环境修复（torch 2.12.1+cpu、缺 openfold） | 未执行 | ✅ |
| 8 | 服务器 shard_0001.faa 截断重建 | 未执行（shards_filt 不受影响） | ✅ |
| 9 | 池外 mid/weak 4,926 条升级补算 | 决策为不升级 | ✅（如需） |
| 10 | 52 条 demote 接入过滤器（`filter_wired`） | 计数变更已执行、过滤器规则未接入 | ✅ |
| 11 | push / GitHub Release / Zenodo DOI | 未执行 | ✅ |
| 12 | 52 列 demote 修订版 15 列→17 列修正记录（可选） | 未执行 | 可选 |
| 13 | nPHAMCL（hfam_4）987 条候选走结构证据验证（序列层已定论不可判别；8YNV Foldseek / 结构预测） | 未执行（本会话定论） | ✅ |

---

## 10. 审查要点（如果我是审查者会追问的 10 个问题）

1. **40,690 是"高可信度"还是"证据最强层"？** —— 项目定义是后者；任何论文级表述必须写明候选语义与 gate 未解锁。
2. **为什么胞内 76:24？** —— 必须同时给出判据不对称（§7.3）与文献研究偏向，不得单讲生物学。
3. **type-1 的分泌信号门为什么可以比参考层严？** —— 参考层从未被该标准检验（§5.3），这是当前最大的证据标准缺口。
4. **Cys 型 29,974 条的催化残基是否逐条验证？** —— 否（inferred_cys_not_verified），且 C183S 不失活使"未命中"不能当强阴性。
5. **校准到什么程度？** —— 仅 hfam_70 六腿通过（文档级判定），未 finalize、未进 registry、未重扫；其余 45 profile 均未校准。
6. **发现层 HMM 能用于筛选吗？** —— 不能（x₁ 混淆集命中 78.69%），只召回；AGENTS.md 已固化。
7. **生长速率是实测吗？** —— 否，gRodon2 预测最大生长潜力；结论是"略有上移、效应量小"。
8. **85.3% 的降解基因组为什么不在分析里？** —— 同属设计固有覆盖限制，结论不可外推。
9. **数字是否可复算？** —— 本会话已对 35,504/35,484/4,321/4,301、40,690 合并、8,882/8,776 账本逐项实测复算通过；哈希与 manifest 齐备。
10. **哪些结论写过但已撤回？** —— 见 §6 的 11 项更正记录，引用任何旧文前必须先查 amendment。
11. **为什么 `--pnone` 重建后 challenge 假阳性清零仍是失败？** —— 因 `eff_nseq=0.00` = 对 5 条训练序列的精确匹配（零召回）；challenge 判定必须同时看假阳性清零与真阳性/泛化保留，不能只看前者。
12. **hfam_70 本轮重建过 HMM 吗？** —— 没有；challenge 用的 production HMM 与参考 profile 逐字段一致（NSEQ=7、EFFN=0.748535），gate 敲定是**状态变更**，不产生新序列、不触发重扫。

---

*本文档为 candidate-only 审查记录。所有 HMM/profile/domain/motif/SignalP/结构/系统发育/分类与生长速率证据只表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
