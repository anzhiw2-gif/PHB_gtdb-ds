# T141 — 实验阳性证据 amendment（PhaZ2 升级 + 2025–2026 新酶定位）：状态文档

**Run：** `runs/20260919_phazed_experimental_evidence_amendment_01`
**日期：** 2026-09-19
**状态：** `completed_candidate_only`
**关键约束：** **冻结台账未修改**（`phaded_reference_ledger.tsv` 只读核对，SHA-256 = `0281911ca2b963bc0fe5565e8f899d1a7b98d5aa01c1aa59ddcb5b1f1c0c655c`）。本 amendment 是后续使用的**权威解释**。

---

## 1. 结论一：PhaZ2 正式升级为 experimental_positive（已完全确认）

| 项 | 证据 |
|---|---|
| 台账记录 | `DED_hfam_61_0029`，accession **AAP74580.1**（现标 `annotation_only`） |
| UniProt 交叉引用（已取回） | **Q7WT49** = EMBL AF549808 → ProteinId **AAP74580.1**；蛋白名 "**Intracellular PHB depolymerase PhaZ2**"；基因 `phaZ2`；C. necator；419 aa（与台账序列长度一致） |
| 实验文献 | **PMID 12813072**（York et al. 2003, J Bacteriol 185:3788-3794） |
| 实验类型 | 基因缺失株生理实验（ΔphaZ1/ΔphaZ2/ΔphaZ3 单/组合缺失 + PHB 合成与利用测定） |
| 原文结论 | "PhaZ1 and PhaZ2 were **sufficient to account for PHB degradation**… PhaZ2 is thus suggested to be an **intracellular depolymerase**" |
| 判定 | **hfam_61：0 → 1 实验阳性；Cys 超家族：1 → 2**（PhaZ1 CAJ92291.1 + PhaZ2 AAP74580.1） |

**明确不计入**：PhaZ3（AAP74581.1，hfam_66）——原文逐字 "The role of PhaZ3 remains to be established"。**不硬凑。**

**门槛状态**：Cys 超家族 2 条实验阳性，仍 < 3；且 held-out、family-resolved 阴性、challenge 仍缺 → **校准 gate 未解锁**。

## 2. 结论二：2025–2026 新实验酶的家族定位（superfamily 级可靠，accession 绑定未完全确认）

| 酶 | 文献 | 最佳家族（E 值） | 超家族 | 该家族现有实验阳性 | 绑定状态 |
|---|---|---|---|---|---|
| Thermanaeromonas toyohensis 多域水解酶 | PMID 40689726 | **hfam_58**（2.4e-224） | type 1 | 1 | 强候选（菌种匹配，但论文具体 deposit 未核实） |
| Nocardiopsis dassonvillei PHBD | PMID 41151231 | **hfam_47**（1.7e-67） | type 1 | 0 | ⚠️ **菌株不符**：D7B6U6 是 ATCC 23218，论文用 NCIM 5124 |
| C. malaysiensis 候选 1（A0ABN4TJ39） | 2025 Bioresource Tech | **hfam_61**（1.7e-205） | Cys | 0 | ❌ 未绑定（论文具体 accession 未知） |
| C. malaysiensis 候选 2（A0ABM6F3M3） | 同上 | **hfam_65**（1.0e-165） | Cys | 1 | ❌ 未绑定 |
| C. malaysiensis 候选 3（A0ABM6FE40） | 同上 | **hfam_45**（4.0e-30） | type 1 | 0 | ❌ 未绑定 |
| C. malaysiensis 候选 4（A0ABM7D8G3） | 同上 | **hfam_57**（2.6e-125） | type 1 | 0 | ❌ 未绑定 |
| P. guguanensis | 同上 | — | — | — | UniProt 无该菌记录，需 GenBank 查询 |

**判定：这批 2025–2026 酶一律不升级为 experimental_positive**——除 Thermanaeromonas 外，accession 与实验的绑定尚未确认，按"不允许硬凑"原则，只记录定位、不计数。

## 3. 若未来完成绑定确认，各家族可获得的实验阳性

| 家族 | 现有 | 可加 | 绑定前提 |
|---|---:|---:|---|
| hfam_61（Cys） | 0 → 1（PhaZ2） | +1（A0ABN4TJ39） | 确认 2025 论文表征的正是它 |
| hfam_65（Cys） | 1 | +1（A0ABM6F3M3） | 同上 |
| hfam_58（type 1） | 1 | +1（Thermanaeromonas） | 核对论文 deposit |
| hfam_47（type 1） | 0 | +1（Nocardiopsis） | 取得 NCIM 5124 天然序列 |
| hfam_57（type 1） | 0 | +1（A0ABM7D8G3） | 确认 2025 论文绑定 |
| hfam_45（type 1） | 0 | +1（A0ABM6FE40） | 确认 2025 论文绑定 |

**即使全部确认，也无一家族达到 3 条独立实验阳性**——最近的是 hfam_61（可到 2）与 hfam_65（可到 2）、hfam_58（可到 2）。

## 4. 产物

- `results/experimental_evidence_amendment.json`（本 amendment 权威记录）
- `inputs/new_enzymes_2025_2026.faa`（6 条新酶序列，见 `runs/20260919_phazed_new_enzymes_01/inputs/`）
- `input_contract.json`（台账只读哈希锁定）

## 5. 边界

- 冻结台账**未修改**；升级只以 amendment 记录。
- 所有记录仍为 candidate-only；实验阳性只表示"有文献支持的实验证据"，**不等同于**已验证 PHB/PHA 降解表型。
