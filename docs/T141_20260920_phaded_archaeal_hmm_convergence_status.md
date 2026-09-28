# T141 — 古菌 HMM 收敛（审查报告第 4 项）：状态文档

**日期：** 2026-09-20
**性质：** 文档/口径类收敛（**无新计算、无服务器执行、不改任何冻结产物**）
**触发：** 项目全面审查报告「接下来该做什么」第 4 项
**边界：** candidate-only。本文所有数字均来自既有冻结产物的**独立逐行重算**，每条都指向具体路径；未编造任何数值或哈希。

---

## 0. 一句话结论

修正了 `knowledge/family_classification.md` §4 的一处**事实错误**（patatin「非 GxSxG」），并把 `ArchPhaZ_patatin`（及 `ArchPhaZ_hydrolase`）**正式降级为「patatin 折叠召回探针」**——两个古菌 HMM 不得再被用于声称「古菌 PHB 解聚酶」。

## 1. 任务 4.1：事实错误更正（patatin「非 GxSxG」）

**被更正的原文（逐字）：** 「| **ArchPhaZ_patatin**（PhaZh1 型） | patatin/磷脂酶 A₂ 折叠（**非 GxSxG**）；颗粒结合 PGAP；…」

**两条独立证据（均实测，指向具体产物）：**

| # | 证据 | 产物路径 |
|---|---|---|
| 1 | Liu 2015（PMID 25710370）原文：PhaZh1 催化残基「**Ser47 (in a classical lipase box, G-X-S47-X-G)**」——**有** GxSxG | https://pubmed.ncbi.nlm.nih.gov/25710370/ |
| 2 | 2,307 条古菌候选中 `nterm_gxsxg_motif` = `present` **2,118** / `absent` **189** | `runs/20260906_overall_pha_phb_archaea_audit_04/results/archaeal_phazh1_full_audit.tsv` 的 `nterm_gxsxg_motif` 列（本次逐行重算；该 run 的 `results/analysis_manifest.json` 记 `counts.archaeal_records = 2307`） |

**更正后表述：** patatin 折叠**有自己的 GxSxG 序列 motif**，但与 α/β-水解酶折叠的 lipase box **结构语境不同**（独立折叠，催化 Ser 位于 patatin 折叠自身的 GxSxG 中）。推论：**不能**用「无 GxSxG」区分 patatin 型与 α/β-水解酶型；**不能**把 §6.2 的 x₁ 疏水 / lipase box 判据套用到 patatin 候选上打分。

## 2. 任务 4.2：正式降级为「召回探针」

**降级依据（全部实测，各指向具体产物）：**

| # | 依据 | 产物路径 / 数值 |
|---|---|---|
| 1 | 校准**无法建立零假阳性阈值**：支持面板 8/8 命中（E 7.4e-145 … 1.5e-93），但广谱/冲突对照（`negative_control`）**4/4 全部命中**（9.2e-51、6.1e-28、3.7e-23、**5.4e-06**，均 ≤1e-5） | `runs/20260904_archaea_phazh1_candidate_01/results/calibration/calibration.tsv` |
| 2 | 8 种子的 `archaea_PhaZh1_like_candidate.hmm`（20260904 构建）**未注册**进 `pipeline/config/formal_scan_models.tsv`、未启动正式扫描 | `docs/T141_20260904_archaea_phazh1_candidate_01_status.md` |
| 2b | **但需精确区分**：`ArchPhaZ_patatin`（103 种子）与 `ArchPhaZ_hydrolase`（12 种子）**确实在** `pipeline/config/formal_scan_models.tsv` 中注册（该表 10 行模型中的第 6、7 行，阈值 e-5），并已被 run-13 正式扫描使用——**这正是本次降级的要害**：注册模型的**角色**必须限定为「召回探针」，其命中不得被读作「古菌 PHB 解聚酶」，也不得进入细菌 PHB 主结论 | `pipeline/config/formal_scan_models.tsv`（本 run 读取确认，10 行：ePhaZ_curated_core / ePhaZ_broad_discovery / iPhaZ / OH / BdhA / **ArchPhaZ_patatin** / **ArchPhaZ_hydrolase** / PhaJ / phasin / PhaC） |
| 3 | 实验锚点只有 **1 个**：PhaZh1（Liu 2015，体外天然 PHA 颗粒水解；体内角色有限，动员主路为 PhaJ β-氧化，Liu 2016 PMID 27052994） | 同上 + https://pubmed.ncbi.nlm.nih.gov/27052994/ |
| 4 | `ArchPhaZ_hydrolase` HMM 仅 **12 条** Halobacteria 种子、无实验锚点 | `data/hmms/v2/ArchPhaZ_hydrolase.hmm` 头部 `NSEQ  12`（`hmmbuild` 产出） |
| 5 | 唯一非嗜盐古菌旁证（*Sulfolobus acidocaldarius* 脂解酶，43% 相似）**PHA 水解活性从未验证** | Arpigny & Jendrossek 1998, PMID 9785454 |
| 6 | 召回结果绝大多数是 patatin 背景：`review_required` 204 / `exploratory` 2,027 / `excluded` 76 = 2,307 | `runs/20260906_overall_pha_phb_archaea_audit_04/results/analysis_manifest.json` 的 `counts.archaeal_audit` |

**允许的表述**：命名与角色一律写作「patatin 折叠**召回探针**（recall-only，只召回不筛选）」。
**禁止的表述**：「古菌 PHB 解聚酶」「patatin 型解聚酶」等结论性措辞。

## 3. 改动文件清单

| 文件 | 改动 |
|---|---|
| `knowledge/family_classification.md` | §4 表格 patatin 行改为更正后措辞；§4 新增 **§4.1 更正痕迹**（逐字引用原文 + 两条证据 + 更正后表述）与 **§4.2 正式降级声明**（6 条实测依据 + 允许/禁止表述）；§4 `ArchPhaZ_hydrolase` 行加「无实验锚点，见 §4.2」；§6.4 两行用途列改为「召回探针（见 §4.2）」 |
| `docs/T141_20260920_phaded_archaeal_hmm_convergence_status.md` | 本文件 |

**原错误表述未被静默删除**：逐字引用保留在 §4.1 的更正痕迹中。

## 4. 边界与未做事项

- 本次为文档收敛，**未重算任何 HMM、未改动任何 HMM/配置/历史 run 产物**。
- `ArchPhaZ_patatin.hmm` / `ArchPhaZ_hydrolase.hmm` 文件本身**原样保留**（禁删）。
- **未解决**（超出本文档范围，属实验/结构验证）：若要让古菌 HMM 从「召回探针」升级到可用，需补充 PhaZh1 之外的 patatin 候选天然 PHA 颗粒水解活性证据，或 Foldseek 对 PhaZh1 结构比对——HMM 调参无法把「patatin 折叠」与「PHB 解聚活性」区分开（训练集里没有足够真阳性）。
