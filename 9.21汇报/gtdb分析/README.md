# 9.21 阶段性汇报 — PhaDED 候选管线图集（方案 A）

> 绘图规范：`Yuan1z0825/nature-skills` 的 **nature-figure** skill（后端 = Python/matplotlib，排他）
> 生成 run：`runs/20260921_phaded_stage_report_figures_02`（方案 A 图集；run _01 的旧版图 2/3 已被取代，原件完整保留）
> **边界：全部数字为候选同源证据（HMM 同源、结构域、基序、定位、分类学、结构可行性），不等同于已验证的 PHB 降解表型。**
> 按要求**不展示池内/池外来源拆分**。

## 图集逻辑：一张图回答一个问题

| 文件 | 回答的问题 | 面板 |
|---|---|---|
| `figure_1_pipeline_funnel` | 从全库到 40,690 发生了什么 | a 管线示意 ｜ b 数据漏斗 |
| `figure_2_evidence_architecture` | **这 40,690 条靠什么证据成立** | a 支撑强度分层（hero）｜ b 架构证据覆盖矩阵 ｜ c 参考层锚点（行标签带 **Intracellular / Extracellular / Periplasmic** 前缀） |
| `figure_3_calibration_and_controls` | **证据有多强、校准到哪一步** | a 六腿 gate 矩阵（hero）｜ b 对照特异性 ｜ c 主动风险控制 |
| `figure_4_genome_taxonomy` | 这些候选在哪里 | a 门级 ｜ b 属级 top 15（**蓝族单色阶：顶部最深 → 向下渐浅**） |
| `figure_5_gaps_and_next_steps` | **还不知道什么、下一步做什么** | a 缺口量化 ｜ b 结构召回路径 + 待授权清单 |

## 汇报时最该讲的三个数（图 2 / 图 3 的核心）

1. **支撑强度分层（图 2a）**：40,690 条中 **trained profile 5,172（12.7%）**、**家族归属歧义 2,905（7.1%）**、
   **仅发现层 HMM 32,613（80.1%）**。占 73.7% 的 Cys 型 nPHASCL **全部**属最后一档——
   而发现层已被本项目自己判定**无判别力**（对 563 条脂酶/酯酶混淆集命中 443 = 78.69%）。
   → 结论应表述为"**候选证据最强层，但支撑强度分层明显；Cys 型需结构证据补强**"。
2. **校准状态（图 3a）**：五个 trained 家族里**只有 hfam_70（type2）六腿全过**，hfam_52 次选（含 2 条弱阳性）；
   hfam_55/hfam_8 各差 1 条独立阳性；**hfam_4（nPHAMCL）多腿失败**（challenge 3/10 假阳性、
   held-out leakage 0.901，训练集实为"1 真阳性 + 2 条同一基因重复"）→ 已降级并改走结构证据。
   37 个 reference-only 家族**连 HMM 都未建**。
3. **对照特异性（图 3b）**：trained 家族对照假阳性 6.0%（3/50，全部来自 hfam_4）；**发现层在混淆集上 78.7%**；
   Cys `[疏水]-C-Q` 基序在非 Cys 参考上仅 4.9%（阳性对照 271/276 = 98.2%）——该基序是**目前唯一**面向 Cys 型的
   序列级判据，但**尚未接入过滤器**（`filter_wired=false`，无 held-out，需在独立阴性集验证 >5% 则降级）。

## 三个必须避开的数字陷阱

1. **30,026 vs 29,974**：`superfamily_status.md` / `evidence_reliability.json` 的候选层是 **52 条 demote 之前**的
   40,742 快照；本图集一律以 `merged_high_confidence_candidates.tsv`（**40,690**）为计数基准。
2. **SignalP 的 `pending` 是陈旧标签**：29,469 行仍写 `pending`，实际已于 2026-09-20 实测为 OTHER（99.82%）
   → 图 2 按"预测非分泌"计，并在图内声明该升级；否则会严重低估定位层覆盖。
3. **with-lipase deferred 层口径**：采用 `deferred_summary.json` 的 **1,206,655**（trained 覆盖 discovery 规则后）；
   `family_discovery` 的 1,238,812 是规则前桶，未使用。该层**禁止删除、不进任何计数**。

## 质量门（详见 run 内 `results/qa_summary.md`）

| 门 | 结果 |
|---|---|
| 源码静态预检（skill `validate_figure.py`） | 21 pass / 0 warn / 0 fail |
| 渲染期多面板对齐（strict, 1.5 pt） | **5/5 PASS** |
| PDF 字形下限 ≥5 pt | **5/5 PASS** |
| 渲染后碰撞审计（skill `audit_figure_collisions.py`，强制） | **5/5 PASS**（0 fail / 0 warn） |

强制审计本轮又发现并修复了 6 类真实缺陷（面板字母压标题、矩阵旋转表头互相重叠、行内两行标签重叠、
窄条标签越界、待授权清单与路径文字相撞、注记 4.8 pt 低于字形下限）。

## 其他口径提示

- 图 4：属级 top 15 用**蓝族单色阶**（顶部最深 → 向下渐浅，深端锚定声明的 `blue_main`），
  仅表达排名、不引入分类色；面板 a 中 Pseudomonadota 深蓝、其余门中蓝、"Other phyla" 聚合行浅灰。
  （附属事实：top 15 属中有 Streptomyces 1,006 与 Micromonospora 241 来自 Actinomycetota，其余 13 属为 Pseudomonadota。）
- 图 4：n = 30,181 / 30,623 基因组有本地 GTDB R232 分类学（98.6%）；442 个待服务器 metadata，未推断其分类位置。
- 图 5b：结构召回路径**产出为 0**（未执行、待授权）；判据为 Foldseek TM≥0.5 且预测 pLDDT≥70，fail-closed。
- **未入图**：PhaDED 层的基因组邻域与生态 isolation source（覆盖未记录，旧方案 A 口径不可搬用）。
- **本目录只含 GTDB 候选发现相关图**；Grodon 生长速率分析已单独成图，位于同级目录 `../降解速率/`。
- 每图另有 PDF 矢量、SVG 可编辑文本；600-dpi TIFF 与 13 个 source-data TSV 见 run 内 `results/`。
