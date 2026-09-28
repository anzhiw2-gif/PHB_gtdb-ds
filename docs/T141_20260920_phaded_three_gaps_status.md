# T141 — 三项证据缺口收口：状态文档

**日期：** 2026-09-20
**run：** `runs/20260920_phaded_three_gaps_01`
**授权：** 操作者授权执行；服务器空闲，SignalP 用 60 线程（动态 `min(60, nproc−load1−10)`）。
**边界：** candidate-only。全程不改 `pipeline/config/`、不改冻结台账、不覆盖历史 `results/`、不 push、不建 DOI。**SignalP 是定位预测、不是实验定位；阳性溯源是对出处（provenance）的更正、不产生新生物学结论。**

三项结论一句话：
1. **定位（①）**：29,521 条"胞内"候选补测 SignalP，**99.82% 被预测确认胞内**；仅 **52 条（0.18%）** 疑似有分泌信号。
2. **阳性（②）**：30 条"实验阳性"的"独立性"**被证伪**——除 `DED_hfam_52` 与 `DED_hfam_70` 外，没有任何家族真正满足"≥3 条独立一手实验阳性"。
3. **SBD（③）**：type-1 高可信度候选"缺 PF06850"是**用错了标记**（coverage_artifact），不是缺底物结合域；**不降级**。

---

## ① SignalP 补测（把"胞内"从推定变实测）

**方法**：从 100 个 `shards_filt` 分片提取 29,521 条序列（66 s，无缺失），SignalP 6.0h（`--organism other --mode fast --torch_num_threads 60 --write_procs 8`，~31 min，returncode 0）。

**对账**：输入 29,521 = SignalP 结果 29,521，`all_matched=true`，无缺失、无多余。

**类别分布**：OTHER **29,469** / SP 49 / LIPO 1 / TAT 2（无 TATLIPO、无 PILIN）。

| 超家族 | 总数 | OTHER（确认胞内） | 会变 export（SP/LIPO/TAT） |
|---|---:|---:|---:|
| Cys（i-nPHASCL 无 lipase box） | 28,534 | 28,482 | **52**（0.182%） |
| nPHAMCL | 987 | 987 | 0 |
| **合计** | 29,521 | 29,469（99.82%） | **52（0.176%）** |

**结论**：此前"≈0.2% 风险"的经验估计被实测证实；**29,469 条（99.82%）的"胞内"从"未测"升级为"预测确认"**；**52 条 Cys 候选疑似携带分泌信号**，已单列 `results/signalp_29521/signalp_export_conflict_52.tsv`（SP 49 / LIPO 1 / TAT 2）。

**待决**（**未擅自执行**）：这 52 条若按实测重跑过滤器，会从 pass 落入 `localization_conflict` hold——**属计数变更，需操作者确认**。其余 29,469 条建议把 `pending` 记为 `OTHER`（实测定位）。

**诚实边界**：SignalP 是**预测**，不是实验定位；本条只把 `pending`（未测）变为"预测标注"，不声称实验验证。

---

## ② 30 条"实验阳性"溯源（gate 的"独立性"被证伪）

**方法**：NCBI esummary/efetch/elink（蛋白 GI → 文献）+ Europe PMC + web_search。

**结果**：30 行 → **28 个独立种子 GI** → **24 篇独立一手文献**（20 篇有活性的 PMID + 4 篇非 PubMed 索引老论文）。按证据类型：克隆/表达+活性 26、纯化+结构 1、纯序列 2、同基因重复 1。

**逐家族"≥3 独立阳性"判定（严格口径）**：

| 家族 | 台账"阳性" | 实际独立一手 | 达标? | 说明 |
|---|---:|---:|---|---|
| `DED_hfam_52`（type1） | 7 | 6 篇 | ✅ | 唯一真正达标的 type1 |
| `DED_hfam_70`（type2） | 7 | 7 篇 | ✅ | 最干净 |
| `DED_hfam_55`（type1） | 3 | **2 篇** | ❌ | 3 基因但 2 篇论文 |
| `DED_hfam_8`（dPHAMCL） | 3 | **2 篇** | ❌ | 3 基因但 2 篇论文 |
| `DED_hfam_4`（nPHAMCL） | 3 | **1** | ❌ | **P26495 真阳性；另 2 条是同一个 P. putida phaZ 基因、纯序列无活性** |
| `DED_hfam_53`（type1） | 2 | **1** | ❌ | **O82950 与 ACG63775.1 是同一 P. stutzeri 基因** |
| `DED_hfam_58` / `65` / `2` / `7` / `3` | 各 1 | 各 1 | ❌ | 均 <3 |

**明确结论**：**gate 的"≥3 条独立实验阳性"只对 `hfam_52` 与 `hfam_70` 成立。** 两个此前被当作"已达标"的家族实为不达标：`hfam_4`（nPHAMCL）训练集 = "1 真阳性 + 2 序列重复"（**实质性溯源缺陷**）；`hfam_53` 离 gate 差 **2** 条（不是 1 条）。`hfam_55`、`hfam_8` 按"不同论文"口径 <3。

**不确定项（如实）**：4 篇非 PubMed 老论文（Zhang&Saito、Yukawa&Uchida、Briese1994、Kobayashi1999）确认存在但未取全文/DOI；BAG32152.1 的 GI 指向 PDB 2D80（结构）；P. putida KT2440 phaZ 未找到专门表征论文。

**产出**：`results/positive_provenance/positive_primary_source_amendment.tsv`（31×17）+ `positive_provenance_summary.json` + `logs/positive_provenance/queries.md`。**未改冻结台账**（amendment 追加制）。

---

## ③ type-1 "缺 SBD" 归因

**verdict = `coverage_artifact`，confidence = high**（已由 InterPro 官方条目独立复核）。

> **溯源更正（诚实声明）：** 这不是新发现——本项目在 **2026-09-17 的 DEFECT_D**（交接文档 §3.2、`reference_residue_mapping_status.md` §DEFECT D）**早已记录**「PF06850 不是胞外 SBD、而是胞内 Cys 型标记（274/276 vs 0/447）」。本任务的增量是：① InterPro 官方条目（`IPR009656` "degrades PHB granules"、`IPR010126` "GO extracellular / 头号引用 PMID 2644188"）的**外部独立确认**；② 补齐 type-1 配件域清单（FN3/Ig-like/TSP3/CHB/Thr-rich/Cad）与"PF10503 = type-1 催化域"的外部出处；③ 修复 `knowledge/family_classification.md` 与 `literature_survey_report.md` 里**与 DEFECT_D 相矛盾的三处残留旧句**（把 PF06850 写成"SBD C 端"）。

- **PF06850 / IPR009656** = "PHB de-polymerase, C-terminal"，描述"degrades **PHB granules**" → 是**胞内（Cys 型）酶的 C 端标记**，不是胞外 type-1 的 SBD。
- **PF10503 / IPR010126** = "Esterase, PHB depolymerase"，GO = **extracellular region**，头号引用 = **A. faecalis 胞外解聚酶（PMID 2644188）** → 是 **type-1/type-2 的催化域**。
- type-1 的配件/连接域 = FN3(PF00041)/Ig-like(PF16403,PF17957)/TSP3(PF02412)/CHB(PF13290)/Thr-rich/Cad。

**证据**：6/13 条 type-1 参考阳性映射到 UniProt 后 **6/6 无 PF06850、全有 PF10503**；候选层 4,884 条池内 type-1 **PF06850=0、PF10503=100%**（FN3 在 InterPro 层 ~14–26%）。

**结论**：5,018 条 type-1 候选**不降级**（100% 有 PF10503 催化域）；SBD 检查应换成 FN3/Ig-like/TSP3/CHB，或把"SBD 必须存在"降级为"记录 caveat"（type-1 的 SBD 是**可选模块**，参考阳性 P52090 等只有催化域）。

**文档更正（已执行）**：`knowledge/family_classification.md` 原第 32 行把 PF06850 写成"胞外 SBD"是张冠李戴，**已更正并保留更正痕迹**（PF06850 归胞内 Cys 型；胞外 SBD = FN3/Ig-like/TSP3/CHB）。

---

## ④ challenge 验证 + hfam_4 判别性重建定论 + hfam_70 gate 敲定（2026-09-20）

**方法**：5 个 trained 家族 HMM × 10 条 challenge 序列（8 条细菌脂肪酶/酯酶 + PcaD γ-内酯酶，见 `challenge_binding_trained_families.json`）跑 `hmmsearch`，阈值 E≤1e-5 判为假阳性。产出：`runs/20260920_phaded_heldout_rebuild_01/results/challenge_validation.tsv`。

**challenge 结果**：

| 家族 | 假阳性 | 判定 |
|---|---|---:|---|
| hfam_52 / 55 / 70 / 8 | 0 | ✅ 10/10 正确拒绝 |
| **hfam_4（nPHAMCL）** | **3** | ❌ **FAIL** |

hfam_4 的 3 个假阳性（全部 ≪1e-5）：P24640（_Moraxella_ 脂肪酶，5.4e-12）、Q02104（_Psychrobacter_ 脂肪酶，4.1e-11）、Q88N36（PcaD，6.6e-11）。

**hfam_4 判别性重建（option b，run `20260920_phaded_hfam4_discriminant_01`）**：对齐 5 阳性 + 10 挑战后逐列判别——① **催化三联体是共享的、非判别**（亲核 Ser 与肘部两侧 G 在挑战里 10/10 保守）；② 判别信号在肘部 `G-V-S-W-G` 的 V/W（挑战 0/10 共享）；③ 根因 `EFFN=0.378`（<1，Dirichlet 先验把 profile 拉向背景）。用 `--pnone` 重建后 3 假阳性消除、但 `eff_nseq=0.00`（阳性 E=0、全挑战 E=999）→ **精确匹配过拟合、零召回，假胜利**。**定论：hfam_4 无法做成判别性家族 HMM**（种子 3~5 条 ~90% 同一 _Pseudomonas_，家族层与超家族层同源，序列信号 = 通用 PF00561）。**操作者已决定暂时不再修改 HMM**；hfam_4 challenge 腿如实记 FAIL，nPHAMCL 改走结构证据层。详见 `runs/20260920_phaded_hfam4_discriminant_01/results/hfam4_rebuild_conclusion.md`。

**trained 家族完整 gate 状态（六腿合成，2026-09-20 实测）**：

| 家族 | positive≥3 | heldout≥1 | negative≥1 | challenge | 负性 status | challenge status | gate |
|---|---|---|---|---|---|---|---|
| **hfam_70（type-2）** | ✅ 7 | ✅ AAB02914.1（id 0.484） | ✅ Q84C08+Q51718 | ✅ 10/10 | sufficient | sufficient | **✅ 六腿全过** |
| hfam_52（type-1） | ✅ 6（~4 强+2 弱） | ✅ AAB40611.1（id 0.304） | ✅ Q84C08+Q51718 | ✅ 10/10 | sufficient | sufficient | ✅ 全过（弱 2 条） |
| hfam_55（type-1） | ❌ 2 | ✅ P52090.1（id 0.701） | ✅ | ✅ 10/10 | sufficient | sufficient | ❌ 差 1 阳性 |
| hfam_8（dPHAMCL） | ❌ 2 | ✅ AAQ72538.1（id 0.840） | ✅ | ✅ 10/10 | sufficient | sufficient | ❌ 差 1 阳性 |
| hfam_4（nPHAMCL） | ❌ 1 | ❌ leakage 0.901 | ✅ PcaD | ❌ 3 假阳性 | sufficient | **not sufficient** | ❌ 多腿失败 |

**hfam_70 gate 敲定（candidate-only determination，非 finalize 脚本执行）**：`family_DED_hfam_70_9818be7f78e3` 六腿全过 → **`calibrated_candidate_model`（判定）**，为 trained 家族中**最干净 gate 候选**（7 条独立强阳性、held-out 干净 id=0.484、challenge 10/10、显式阴性 Q84C08+Q51718）。hfam_52 也六腿全过但含 2 条弱阳性（BAA04986 K1 探针章节、BAA82057 A1「Published Only in Database」），故**次选**。判定件：`runs/20260920_phaded_heldout_rebuild_01/results/gate_determination_hfam_70.json`。**未运行 finalize 脚本、未改 config、未授权全库重扫**（全库重扫仍须单独授权）。

**澄清：本轮 HMM 是否重建（实测，防后续会话误判）**：challenge 用的 `hfam_70_production.hmm` 与参考 profile `family_DED_hfam_70_9818be7f78e3.hmm` **逐字段一致（NSEQ=7、EFFN=0.748535、NAME 相同）**——**hfam_70 未重建**，阳性集 7 条（AAA87070.1/AAB02914.1/AAT09963.1/BAA35137.1/BAA92354.1/BAG32152.1/O05527）自 09-10 起未变。本轮唯一「重建」的临时 6-train HMM 仅用于 held-out 验证 AAB02914.1，非生产 HMM、challenge 未使用。**hfam_52/55/8 阳性集同样未变、未重建；hfam_4 虽重建（5 条修正阳性 + `--pnone`）但已降级 `reference_only_insufficient_panel`、不进正式扫描。** 因此**本轮无任何 profile 需要重扫全库**：hfam_70 的全库召回早已在 `20260919_phaded_trained_profile_recall_01` 跑完（39,458 条原始命中），「gate 敲定」只是给 type-2 候选贴 calibrated 状态，不产生新序列、无需重扫。全库重扫是另一个 deferred 动作，仍须单独授权，且只有 hfam_70 有资格以 calibrated 身份进入。

---

## 影响汇总

1. **type-1 "最可信"的排名不再背着"SBD 全缺"的疑点**——那是工具假象，且项目 DEFECT_D 早已指出 PF06850 非胞外 SBD；5,018 条 100% 有正确催化域标记 PF10503。
2. **gate 的"独立阳性"前提被证伪**（除 52/70）；尤其 **nPHAMCL（hfam_4）训练 profile 的训练集是"1 真阳性 + 2 序列重复"**，这是需要单独处理的实质缺陷。
3. **Cys 30,026 条里 52 条疑似分泌（0.18%）**，其余 99.8% 定位由预测确认。
4. **30 条阳性的二手出处（Knoll 2009）已映射到 24 篇一手文献**（作为 amendment，未覆盖冻结台账）。

## 四项收尾（操作者授权 2026-09-20，均已执行）

1. **✅ 52 条 SignalP-export 的 Cys 候选已 demote**：36,611 → **36,559**（池内），合计高可信度 40,742 → **40,690**。修订版 `results/cys_demotion/high_confidence_candidates_revised.tsv` + `hold_candidates_addendum_52.tsv`（hold_reason = `localization_conflict_signalp_predicted_export`）。**冻结原表未覆盖**；SignalP 是预测，demote 可逆。
2. **✅ `hfam_4` trained profile 已降级**：`family_DED_hfam_4_3c28e8cee1c4` 从 `trained` → `reference_only_insufficient_panel`（`results/profile_status_amendment_hfam_4.json`；理由 = 训练集实为"1 独立阳性 + 2 条同一 P. putida 基因的纯序列重复"）。**未删 HMM、未改 config、未重训**；重训需 ≥2 条新独立阳性 + 授权。gate 账本相应更正：hfam_4 = 1 独立阳性、差 2。
3. **✅ type-1 SBD 判据已闭环**：`results/type1_sbd/type1_sbd_annotation_summary.md` 明确——**过滤器 pass/hold 逻辑无需改**（type-1 判据本就不含 SBD）；PF06850 归属已在知识库更正；可选配件域标注用 FN3/Ig-like/TSP3/CHB。
4. **✅ 4 篇非 PubMed 老论文已补出处**：`results/positive_provenance/non_pubmed_4_papers_retrieval.md`——Kobayashi 1999（Acidovorax TP4，纯化+克隆+Ser20 突变，**强**）、Briese 1994（P. lemoignei 五基因，**基因清单级**）已确认完整出处；R. pickettii A1/K1 两篇定位到 1994 Elsevier 论文集章节，全文/作者卷页**待取**。

## 仍待决（未执行）

1. `hfam_4` 与 nPHAMCL 超家族 profile 的**重训**（需 ≥2 条新独立阳性 + 新 HMM + 授权；当前仅 1 条独立阳性 P26495，**不足够**）。
2. ~~R. pickettii A1/K1 全文~~ **已解决**：A1（Zhang & Saito）= GenBank "Published Only in Database"，**无论文**；K1（Yukawa/Uchida/Kohama/Kurusu 1994）= "Monitoring of polymer biodegradabilities…by a DNA probe method"（书籍章节，DNA 探针方法，非酶活性表征）→ hfam_52 的"6 篇独立论文"应修正为"~4 篇活性表征 + 2 条弱证据"。
3. ~~type-1 配件域标注是否落成正式证据层~~ **已落成**：`type1_accessory_domain_layer.tsv`（5,018 条；4,884 已标注，134 池外 pending；≥1 配件域 1,430=29.3%）。
4. **hfam_70（type-2）已 gate 敲定（六腿全过，见 §④）**，但 **formal finalize（`finalize_phaded_reference_panel_calibration.py`）与全库重扫仍须单独授权**；在此之前 gate 判定仅作文档记录，**不写入 `pipeline/config/formal_scan_models.tsv`、不进入任何全库重扫**。
