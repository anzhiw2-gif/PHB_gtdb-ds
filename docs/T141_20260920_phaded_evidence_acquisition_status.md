# T141 — 证据获取（解锁校准 gate）：状态文档

**Run：** `runs/20260920_phaded_evidence_acquisition_01`
**日期：** 2026-09-20
**状态：** `completed_candidate_only`
**任务来源：** 全面审查报告「接下来该做什么」第 1 项（最高杠杆）。
**性质：** 只读文献/数据库检索 + 判定 + 方案设计。**未修改**冻结台账、`pipeline/config/`、任何历史 run 或 `results/`；未编辑共享交接文档。

---

## 0. 一句话结论

**gate 未解锁，且 2025–2026 新酶路线比原 amendment 预期更弱**：Thermanaeromonas 在**生物体层面确认**（但论文显式 accession 未读到），Nocardiopsis **菌株不符未解决**，C. malaysiensis 的 4 条候选里**有 2 条被证伪为"论文所表征的酶"**（论文的 CmaPHBd 是胞外 type I，而候选 1/2 不是），另 2 条无法区分；P. guguanensis 在 UniProt **无任何记录**。**没有任何家族因此达到 ≥3 阳性**。

## 1. 绑定判定表（逐条，含证据）

| 目标 | UniProt | 判定 | 关键证据 |
|---|---|---|---|
| Thermanaeromonas toyohensis 多域 PHB 水解酶 | A0A1W1W0M1 | **CONFIRMED（生物体层面）** | 614 aa；*T. toyohensis* ToBE；PF10503 + **NF050058 (e_dPHAscl_type1)** + TIGR01840 + PF13290 (CHB_HEX_C_1) + IPR059177 (GH29D-like)；**UniProt 检索证明该物种全基因组只有这 1 条 PF10503 蛋白**；C 端碳水化合物样域与论文原文"C-terminal domains … share similarities with carbohydrate degrading enzymes"吻合 |
| Nocardiopsis dassonvillei PHBD | D7B6U6 | **UNRESOLVED（菌株不符）** | 417 aa；**ATCC 23218（模式株）**；但 PMID 41151231 用 **NCIM 5124** + 体外**密码子优化人工合成基因**；MW ~50 kDa 与带标签构建一致但不构成确认 |
| C. malaysiensis 候选 1 | A0ABN4TJ39 | **REFUTED（非论文之酶）** | RefSeq WP_071070267.1（406 aa）仅注释 "polyhydroxyalkanoate depolymerase"，**不是** "extracellular catalytic domain type 1"；而论文的 CmaPHBd 明确是胞外 type I |
| C. malaysiensis 候选 2 | A0ABM6F3M3 | **REFUTED（非论文之酶）** | RefSeq WP_071037164.1（425 aa，MULTISPECIES Cupriavidus）同样只是泛注释 |
| C. malaysiensis 候选 3 | A0ABM6FE40 | **UNRESOLVED（二选一）** | RefSeq WP_071072601.1（488 aa，BKK80_30885）注释为 "**extracellular catalytic domain type 1** SCL PHA depolymerase"——与论文一致；但无法与候选 4 区分 |
| C. malaysiensis 候选 4 | A0ABM7D8G3 | **UNRESOLVED（二选一）** | RefSeq WP_071072778.1（361 aa，BKK80_32250）同样注释为胞外 type I |
| P. guguanensis PguPHBd | 无 | **UNRESOLVED（多候选）** | UniProt **无该物种任何记录**；NCBI RefSeq 有 6 条胞外 type I 注释（575–576 aa 等），无法判定 |

**论文定位：** PMID **41173112**，DOI `10.1016/j.biortech.2025.133571`，*Bioresour Technol* 441:133571 (2026, e-pub 2025-10-29)，Bodourian et al.——它表征了**两个**酶：**CmaPHBd**（C. malaysiensis）与 **PguPHBd**（P. guguanensis），二者均为 **"extracellular type I depolymerase"**（N 端 CD + C 端 SBD）。

> **对原 amendment 的重要更正**：amendment 曾希望「Cys 家族 hfam_61/hfam_65 各可得一条阳性（来自候选 1/2）」。**该希望不成立**——论文的 CmaPHBd 是胞外 type I，候选 1/2 不是论文之酶；若最终绑定成功，受益的是 **type 1 家族（hfam_45 或 hfam_57）**，与 Cys 家族无关。

## 2. 显式实验阴性清单

**有效阴性（可用）**
1. **Q84C08** — *Aquipseudomonas alcaligenes* LB19 **Mcl-PHA 解聚酶**（278 aa；NCBIfam **NF050057 = e_dPHAmcl**；EC 3.1.1.76；PMID 15995648）。显式"不水解 PHB"。family 可绑定（胞外 dPHAMCL 超家族；具体 hfam 需跑 discovery 层 HMM）。
2. **R. rubrum 周质（原称"胞内"）PHB 解聚酶**（PMID 15489436）— 对**天然 PHB** 特异、**不水解变性 PHB**；属周质家族（PhaDED 家族 3）。**accession 待补**。
3. **P. fluorescens GK13 P(3HO) 解聚酶**（Schirmer 1993 AEM 59:1220 / 1995）— MCL 特异、不水解 SCL-PHA(PHB)。**accession 待补**（UniProt 搜 `P. fluorescens + PF10503` 返回 0 条）。

**被证伪、不得当作阴性（如实列出）**
- `O87189` — 与台账 `DED_hfam_65_0001`（**CAJ92291.1 = PhaZ1**）逐字节相同 → 它是**阳性**。
- `Q7WT49` — = **AAP74580.1 = PhaZ2**，09-19 amendment 已升级为**阳性**。
- `Q7WT48` — = PhaZ3，原文"其作用尚待确定"，**未确定 ≠ 阴性**。
- `Q9AGB6` / `J7K890` — 面板设计性排除（关联文献是合成酶 / 片段低于完整性阈值）。

## 3. Held-out 阳性拆分方案（可执行）

已写入 `results/heldout_split_plan.md`，要点：
1. **触发**：仅当家族独立实验阳性 ≥ 3 时执行。
2. **选点（确定性）**：优先选**跨属**阳性；多条时取 taxonomy_id 最小者；无跨属则取与其余阳性**成对 identity 最低**者；并列取 accession 字典序最小者；优先留出有 PDB/文献锚者。
3. **无泄漏**：held-out **不得进入** hmmbuild 的 MSA；held-out 对训练集**最大成对 identity 必须 < 90%**（否则换下一条，全 ≥90% 则写 `leakage_blocked`）；记录 held-out sha256 + **具体比对哈希**（MAFFT 7.525 不可位复现，须声明）；held-out **不得用于选阈值**。
4. **验收**：held-out 必须被该家族 HMM 检出且 ≥ 阈值；失败结论是"HMM 泛化不足"，**不得**放宽阈值或把 held-out 移回训练。

**可立即执行者**：已 trained 且有 3–7 条阳性的 5 个家族（**hfams_4/52/55/70/8**）——它们已达阳性阈值，held-out 槽位是当前最被低估的可解项。

## 4. Gate 差距表（`results/gate_distance.tsv`）

| 家族 | 现有阳性（含 amendment） | 距 ≥3 还差 | held-out | family-resolved 阴性 | challenge |
|---|---:|---:|---|---|---|
| **DED_hfam_61**（Cys） | 1（PhaZ2） | **2** | 缺 | **缺** | sufficient |
| **DED_hfam_65**（Cys） | 1（PhaZ1） | **2** | 缺 | **缺** | sufficient |
| **DED_hfam_58**（type 1） | 1（BAF35850.1） | **2** | 缺 | **缺** | sufficient |
| **DED_hfam_53**（type 1） | 2 | **1** | — | — | — |
| **DED_hfam_47 / 57 / 45**（type 1） | 0 | **3** | 缺 | **缺** | sufficient |
| hfam_4 / 52 / 55 / 70 / 8（trained） | 3 / 7 / 3 / 7 / 3 | **0（阳性已达标）** | 未记录 | 未记录 | 未记录 |
| Cys 超家族（汇总） | 2（PhaZ1+PhaZ2） | 1（按家族计） | — | — | — |

> **全部 37 个 reference-only profile 的 `heldout_positive = 0`、`family_resolved_negative_status = missing`**（实测 `profile_calibration_readiness.tsv`）；**该表中没有任何家族阳性 ≥3**，故 gate 对 reference-only 层仍然全锁。

## 5. 结论与下一步具体动作

**gate 是否解锁：否。** 没有任何家族达到 `positive_count >= 3`。

**最高价值的三个具体动作（按可执行性排序）**
1. **立即为 5 个 trained 家族（hfams_4/52/55/70/8）执行 held-out 拆分**（阳性已达标，规则已备，不需要新证据）——这能直接点亮 5 个家族的 `heldout >= 1` 槽位。
2. **补齐 family-resolved 阴性**：先固定 `R. rubrum` 周质酶与 GK13 的 accession，再跑 discovery 层 HMM 把 `Q84C08` 与二者绑到具体 hfam；这是 37 个 profile 共同的缺口。
3. **2025–2026 绑定只能靠"读到论文正文/补充材料"才能推进**：本轮因 Wiley/PMC/ScienceDirect 全部拦截（403 / reCAPTCHA / 403）未能取得论文的显式 accession 行。**若拿到 Thermanaeromonas 论文正文的 accession 行 → hfam_58 可 1→2；若拿到 C. malaysiensis 补充材料的序列 → 可在 hfam_45 与 hfam_57 之间定一条（0→1）**。这两者都**不足以单独解锁任何家族**（还差 2 或 1 条）。

**诚实边界**：本轮所有判定均基于可复现的数据库检索；未读到论文正文的显式 accession 处一律写 `unresolved`，**未编造任何 accession、哈希或计数**。

## 6. 产物

- `runs/20260920_phaded_evidence_acquisition_01/input_contract.json`（gtdb 三槽 pending；含输入/输出 SHA-256）
- `results/binding_verification.tsv`（7 行绑定判定 + 证据 + 阻塞原因）
- `results/explicit_negatives.tsv`（3 条有效阴性 + 5 条被证伪/排除项）
- `results/heldout_split_plan.md`（可执行拆分规则 + 无泄漏保证）
- `results/gate_distance.tsv`（逐家族 gate 差距）
- `logs/evidence_queries.md`（全部检索端点、查询与受阻记录）
