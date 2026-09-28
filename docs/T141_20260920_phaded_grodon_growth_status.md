# T141 — PhaDED 降解菌 gRodon2 预测最大生长速率分析：状态文档

**Run：** `runs/20260920_phaded_grodon_growth_01`（本地与服务器同名）+ `deploy/20260920_phaded_grodon_growth_01/`
**日期：** 2026-09-20（服务器计算跨至 09-21 00:09 CST）
**状态：** `completed_candidate_only`
**授权：** 操作者确认三项设计（40,690 集 / 同属平衡 + 组间 delta / 复用旧预测）；操作者临时授权 70 线程（超出 AGENTS.md 单任务 40 上限，见 `logs/server_run_notes.md`）

> **边界声明：** gRodon2 输出是基于密码子使用偏好的**预测最大生长速率**，**不是**特定培养条件下的实测生长速率。降解菌标签来自高可信候选分类（candidate-only），**不等同于已验证的 PHB/PHA 降解表型**。

---

## 0. 一句话结论

在 GTDB R232 的**同属平衡**比较中（8,776 基因组 = 4,388 降解 + 4,388 同属非降解，1,216 属），降解菌的预测最大生长速率**略高于**同属非降解菌（属水平平均 Δ = **+0.0255 h⁻¹**，属内分层置换 P = **0.0040**，Wilcoxon P = **0.0033**），但**效应量很小**（r = 0.084）、中位差 CI 含 0、符号检验不显著（P = 0.115）；**胞内与胞外降解菌之间无显著差异**（组间 Δ = −0.0135 h⁻¹，Mann-Whitney P = 0.92，置换 P = 0.48）。

---

## 1. 设计

| 项 | 值 |
|---|---|
| 降解菌集 | 40,690 合并高可信候选 → **30,623 基因组**（`inputs/degrader_genomes.tsv`） |
| 对照 | 同属内、**不在 62,666 候选池基因组**中的 GTDB 基因组，属内 1:1 平衡，seed 42 |
| 分组 | 纯胞内（nPHASCL±lipase box + nPHAMCL）22,160 基因组 / 纯胞外（dPHASCL1/2 + dPHAMCL + PhaZ7-like + periplasmic）5,644 / 兼有 2,819 |
| 胞内 vs 胞外 | 各自做同属平衡 → 属水平 delta → **组间 delta 对比**（Mann-Whitney + 置换 + bootstrap CI），仅用**纯组** |
| 工具 | Pyrodigal CDS → Pfam 核糖体蛋白 HMM（HEG）→ gRodon2 `predictGrowth` |
| 复用 | 旧工作区 `PHB_gtdb` 的 8,788 条预测中，与新 manifest 交集 **713 条**直接复用（同工具同数据同 HMM 路线） |
| 线程 | R 侧强制单线程（`OMP/OPENBLAS/MKL/R_DATATABLE_NUM_THREADS=1`）；worker 数 35→20→40→70（操作者授权） |

**同属设计的固有覆盖限制（必须与结论同读）：** 30,623 个降解基因组中仅 **4,441（14.5%）** 在属内有可配对的非降解基因组，其余 26,116 个因属内"全是候选"被排除；因此本分析覆盖的是**降解菌与非降解菌混合存在的属**，不能外推到全部门类。

---

## 2. 数据账本（`grodon_growth_statistical_tests_newdeg40k.tsv`）

| 项 | 数量 |
|---|---:|
| gRodon2 输入基因组（manifest） | 8,882 |
| 预测成功 | 8,828 |
| 预测失败 | 54（0.61%；全部 `too_few_ribosomal_hits`，涉及 46 属） |
| `status=ok` 但无有效生长率 | 2 |
| 为严格同属平衡额外剔除 | 50 |
| **最终平衡基因组** | **8,776** |
| └ 降解 | **4,388** |
| └ 同属非降解 | **4,388** |
| 最终覆盖属 | **1,216** |

---

## 3. 主结果 A：降解菌 vs 同属非降解菌

**分析单位为属**（属水平 delta = 该属降解菌平均 − 该属对照平均）。

| 指标 | 结果 |
|---|---:|
| 属水平平均 Δ | **+0.02553 h⁻¹**（bootstrap 95% CI +0.00647 ~ +0.04460，**不含 0**） |
| 属水平中位 Δ | +0.00231 h⁻¹（CI −0.00028 ~ +0.00571，**含 0**） |
| 按配对数加权平均 Δ | +0.04369 h⁻¹（CI +0.00251 ~ +0.10872） |
| **Wilcoxon 符号秩检验 P** | **0.00326**（效应量 r = 0.0844） |
| 精确符号检验 P | 0.11470（Δ>0 的属 636 / Δ<0 的属 580） |
| **属内分层置换检验 P** | **0.00399**（非加权）/ **0.00439**（加权） |

**描述性基因组水平**（非独立样本，不作主检验依据）：降解 0.46604 vs 非降解 0.42235 h⁻¹（均值）；0.25885 vs 0.24008（中位）。

**解读**：存在一个**统计上可检出、但幅度很小**的正向位移，且主要由分布上尾驱动（中位差 CI 含 0、符号检验不显著）。**不得**表述为"降解菌生长更快"的强结论；恰当表述为"降解菌相对同属非降解菌的预测最大生长速率略有上移，方向一致性弱（636:580）"。

---

## 4. 主结果 B：胞内 vs 胞外降解菌

| 组 | 属数 | 平均 Δ (h⁻¹) | 中位 Δ | Wilcoxon P | t 检验 P | 符号检验 P | Δ>0 / Δ<0 属 |
|---|---:|---:|---:|---:|---:|---:|---|
| 胞内（纯） | 676 | +0.01990 | +0.00292 | 0.0387 | 0.0214 | 0.233 | 354 / 322 |
| 胞外（纯） | 561 | +0.03344 | +0.00272 | 0.0616 | 0.0665 | 0.353 | 292 / 269 |
| **组间（胞内 − 胞外）** | — | **−0.01354** | — | — | — | — | — |

- **组间差 bootstrap 95% CI：−0.05379 ~ +0.02557（含 0）**
- **Mann-Whitney U P = 0.916；组标签置换检验 P = 0.483**

**结论：胞内与胞外降解菌的预测最大生长速率没有可检出的差异**；两组相对各自同属对照均显示同向的小幅正 delta（胞内 P=0.039、胞外 P=0.062），但组间不区分。

---

## 5. 主结果 C：超家族分层（6 个达到最小属数）

| 超家族 | 属数 | 平均 Δ | 中位 Δ | Wilcoxon P | **BH q** | 方向（+/−） |
|---|---:|---:|---:|---:|---:|---|
| extracellular dPHASCL type 1 | 397 | +0.04361 | +0.00378 | 0.0820 | 0.340 | 210 / 187 |
| intracellular nPHASCL without lipase box（Cys） | 700 | +0.00957 | +0.00154 | 0.1132 | 0.340 | 360 / 340 |
| extracellular dPHASCL type 2 | 203 | +0.06727 | +0.00093 | 0.4777 | 0.638 | 104 / 99 |
| extracellular dPHAMCL | 26 | −0.04392 | −0.00509 | 0.5317 | 0.638 | 11 / 15 |
| extracellular native-SCL/PhaZ7-like | 9 | +0.01526 | +0.01826 | 0.3594 | 0.638 | 7 / 2 |
| intracellular nPHAMCL | 4 | +0.10320 | −0.00734 | 0.8750 | 0.875 | 1 / 3 |

**BH 校正后无任何超家族显著**（q ≥ 0.34）；名义上最接近的是 type 1（P = 0.082）与 Cys 型（P = 0.113）。

---

## 6. 与旧分析（`PHB_gtdb`）的对照

| 项 | 旧分析（旧 phaZ 定义） | 本分析（40,690 新定义） |
|---|---|---|
| 平衡基因组 / 属 | 8,692 / 899 | 8,776 / 1,216 |
| 属水平平均 Δ | −0.00273 h⁻¹ | **+0.02553 h⁻¹** |
| Wilcoxon P | 0.459 | **0.0033** |
| 符号检验 P | 0.386 | 0.115 |
| 分层置换 P | 0.670 / 0.231 | **0.0040 / 0.0044** |
| 亚型分层 | 旧 5 亚型全部不显著（q ≥ 0.25） | 新 6 超家族全部不显著（q ≥ 0.34） |
| 胞内 vs 胞外组间对比 | 旧脚本未做组间检验 | **−0.0135 h⁻¹，P = 0.48（无差异）** |

**如实说明差异原因（不擅自裁决）：** 两次分析的**定义、候选集与属覆盖都不同**（旧 phaZ 定义 vs 新 8 超家族高可信集；899 vs 1,216 属；基因组集合重叠仅 717 条）。旧分析"无差异"与本分析"小幅正 delta"并不互相推翻，但**也不能互相印证**；本分析的效应量很小（r=0.084），任何引用都必须同时给出中位差 CI 含 0 与符号检验不显著这两条限定。

---

## 7. 事故与留证

1. **R 多线程失控（2026-09-20 22:53）**：gRodon 依赖的 data.table 默认按核数并行，35 个并发 R 进程使 load1 达 **74.85/80**，违反"保留约 10 核余量"。处置：显式导出 `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1`（R 侧单线程封顶）后重启；此后即使 70 线程，最高 CPU 进程也只是 hmmsearch 150%，load 稳定 ~70-73。
2. **重启导致重复行**：`read_done` 只把 `status=ok` 记为已完成 → 4 个反复失败的基因组在 4 次重启中各写 1 行（多出 12 行）。处置：新增 `dedupe_predictions.py`（本地 2 单测，先红后绿）去重，原始表**保留**（`grodon_growth_predictions_newdeg40k.tsv`），去重表 `_dedup.tsv` 供下游，去重统计落 `predictions_dedup_stats.json`（8,894 → 8,882，去重 12）。
3. **组间脚本 pandas 索引对齐 bug**：掩码取自过滤子帧时 `df[mask]` 抛 `IndexingError`。处置：`genus_deltas` 内 `pos_mask.reindex(df.index, fill_value=False)` + main 内掩码改由 `df` 构造；新增回归测试 `test_genus_deltas_accepts_mask_from_filtered_subframe`。
4. **pkill 自匹配**：inline ssh 的 pkill 模式串出现在远程 shell 命令行里会自杀该 shell（首次 kill 无输出、exit 1）；改用 `[.]` 正则规避。
5. **scp 通配符不展开**（Windows OpenSSH）：`scp "dir/*"` 静默不传；改为显式文件列表。

---

## 8. 验收

| 检查 | 结果 |
|---|---|
| 新增单测 | `test_build_grodon_manifest`(6) + `test_grodon_group_stats`(6) + `test_dedupe_predictions`(2) = **14 tests OK**（先写失败测试） |
| `python -m compileall -q pipeline deploy/20260920_phaded_grodon_growth_01` | exit 0 |
| `git diff --check` | exit 0 |
| deploy 绑定 | 8 个脚本上传哈希与本地逐一致（`sha256sum` 实测） |
| 预测对账 | manifest 8,882 = 去重后 8,828 ok + 54 failed；平衡 8,776 = 4,388 + 4,388 |
| 旧工作区 | **未修改**（只读取旧预测表与 HMM/脚本；新产物全部在本 run 与 dated deploy） |
| `pipeline/config/`、冻结台账、历史 results | 未改动；未 commit/push |

---

## 9. 产物

| 产物 | 位置（`runs/20260920_phaded_grodon_growth_01/results/`） |
|---|---|
| 同属匹配 manifest（8,882） | `grodon_growth_manifest_newdeg40k.tsv` + `manifest_stats.json` |
| 原始预测（含重启重复行，保留） | `grodon_growth_predictions_newdeg40k.tsv`（8,895 行） |
| **规范预测表** | `grodon_growth_predictions_newdeg40k_dedup.tsv`（8,882）+ `predictions_dedup_stats.json` |
| 同属平衡表 / 汇总 / 审计 / 失败 | `grodon_growth_balanced_by_genus_newdeg40k.tsv`、`..._summary_...tsv`、`grodon_growth_balance_audit_newdeg40k.tsv`、`grodon_failed_*_newdeg40k.tsv` |
| 属内选中计数 | `grodon_growth_balanced_genus_selected_counts_newdeg40k.tsv` |
| **主统计（A）** | `grodon_growth_statistical_tests_newdeg40k.tsv` + `grodon_growth_genus_effects_newdeg40k.tsv` |
| **超家族分层（C）** | `grodon_growth_superfamily_statistical_tests_newdeg40k.tsv` + `..._superfamily_effects_...tsv` |
| **胞内 vs 胞外（B）** | `grodon_growth_group_statistical_tests_newdeg40k.tsv` + `..._group_effects_...tsv` |
| 图（统计脚本自带） | `figures/figure5_grodon_growth_comparison_newdeg40k.{svg,pdf,png}` + `source_data/` |
| **图（nature-figure 风格，6 面板）** | `figures/nature/growth_comparison_nature.{svg,pdf,png}` + `figures/nature/source_data/`（6 张源数据表）；脚本 `pipeline/scripts/plot_phaded_grodon_growth_nature.py` |
| 运行记录 | `logs/server_run_notes.md`（含线程授权、事故、修法） |

### 9.1 nature-figure 风格图（2026-09-21 追加）

沿用项目既有 nature-figure 方法（参考 `D:\16s\...\nature_figure_growth_comparison.py` 与 `D:\PHB_gtdb\scripts\nature_figures_common.py`）：Arial 优先、`svg.fonttype="none"`（**SVG 内 114 个可编辑 `<text>` 元素**）、`pdf.fonttype=42`、7.2 in 双栏宽、固定调色板、面板字母 9 pt 粗体 + 自动碰撞 QC（实测 `overlaps: none`）、SVG/PDF/PNG(600 dpi) 三格式导出并同步写 source data。

六个面板：**a** 基因组水平分布（箱线+抖动，Mann-Whitney P = 0.021、Cliff's δ = +0.028）；**b** 属水平 delta 分布（含 signed-rank / 分层置换 / 符号检验）；**c** 配对属均值（1:1 线、色=delta、大小=配对数）；**d** 胞内 vs 胞外（小提琴 + 组间检验）；**e** 超家族森林图（bootstrap 95% CI + BH q）；**f** 数据账本与结论。

> **注**：本图为 candidate-only 展示；标题与脚注均明写"预测最大生长速率、非实测"与"降解标签为候选"。面板 a/b 的 x 轴按 1st–99th（或 2nd–98th）百分位裁剪，已在图内标注。

---

*本文档为 candidate-only 记录。gRodon2 给出的是密码子使用偏好推断的**预测最大生长速率**，不等于实测生长速率；"降解菌"标签为高可信候选分类，不等同于已验证的 PHB/PHA 降解表型。*
