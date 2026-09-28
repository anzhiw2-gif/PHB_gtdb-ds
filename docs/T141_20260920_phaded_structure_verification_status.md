# T141 — with-lipase 结构验证层启动（PDB 8YNV）：状态文档

**日期：** 2026-09-20
**目标：** 为 `intracellular nPHASCL with lipase box`（deferred 层，1,206,655 条池外命中）建立唯一可行召回路径——结构证据（PDB 8YNV）。
**状态：** `feasibility_measured_and_piloted`（可行性已实测 + 小规模试点完成；**未启动全量计算**）
**边界：** candidate-only。deferred 层**禁止删除**，其任何结果**不得进入候选/高可信度计数**。

---

## 0. 一句话结论

结构验证路径**技术可行**：8YNV 可下载、Foldseek 可用且已验证、ColabFold(AF2/JAX) 可用、deferred 序列可提取。**但必须走「严格 E 值漏斗」**：hfam_2 对 8YNV 的得分是 **E 5.8e-165**，而 1.2M deferred 中 **99.6% 是 E>1e-30 的弱命中**，只有 **E<1e-50 的 1,550 条**（或 E<1e-30 的 4,013 条）值得做结构预测。全量 1.2M 结构预测**不可行**（~324 GPU-日）。

## 1. 可行性实测（服务器，只读）

| 项 | 结论 | 证据 |
|---|---|---|
| **8YNV 可下载** | ✅ | 服务器有外网；`curl https://files.rcsb.org/download/8YNV.pdb` → HTTP 200，1,750,491 B，9,538 ATOM |
| 8YNV 序列 | ✅ | 4 链 A/B/C/D（298/296/297/296 aa）；链 A 序列含 **规范 GxSxG `GWSMG`** + **Asp motif `DRDYVV`** |
| **Foldseek** | ✅ | `${PHB_REMOTE_ROOT}/tools/foldseek_20260912/foldseek/bin/foldseek`，version `463739e0…`；自测 8YNV vs 8YNV → **TM 1.00 / E 9.8e-70** |
| **ColabFold** | ✅ | conda env `colabfold`（colabfold 1.6.1，JAX 0.5.3 见 GPU0）；AF2 单序列 296 aa → **84 s**，产出 PDB |
| ColabFold 模型 | — | 可用：alphafold2 / alphafold2_ptm / multimer v1–v3 / **deepfold_v1**；**无 ESMFold 选项** |
| **ESMFold** | ❌ 环境损坏 | esmfold env 的 torch 是 **`2.12.1+cpu`（CUDA 不可用）**，且**缺 `openfold`** → 无法 folding |
| MMseqs2 / blastp / diamond | ✅ | mmseqs2-18.8cc5c；phb_gtdb env 的 blastp、diamond |
| 本地 colabfold MSA 库 | ❌ 无 | `~/.cache/colabfold/` 只有 `params/`（12 个模型）；MSA 模式需 mmseqs2 互联网服务 |
| **deferred 序列可取** | ✅ | 抽样 1,001 条 → **1,001/1,001 全部取到**，扫 100 分片（268 GB）**耗时 293 s** |
| GPU | ⚠️ | 2× RTX 4090；**GPU0 空闲**（1 MiB, 0%），**GPU1 被他人占用**（23.8 GB, 62–100%） |
| 磁盘 | ✅ | `/home/data` 可用 100 T |

## 2. 试点结果（抽样 1,001 条，已执行）

1. **提取**：1,001/1,001 成功，293 s（成本在扫分片，故全量 1.2M 提取耗时相同）。
2. **序列预筛（blastp vs 8YNV）**：初版记录为「0/1,001 达到 E<1e-5」——**该数值错误，见 §2.2 更正记录**。更正后：987/1,001 在 E≤1e-3 有命中，其中 **E<1e-5 = 954**、最小 **1.16e-122**；pident 中位 24.3%、最大 54.1%。
   → 结论方向保留、证据替换：**对单条参考做 blastp 不具判别力，不能用作过滤器，无命中也不得记为阴性**（详见 §2.2）。
3. **motif 普查**：**78.9%（790/1,001）有规范 GxSxG**；x₁ 疏水 619/790（78%），分布 M 257 / L 167 / F 51 / W 48。
4. **hfam_2 对 8YNV**：**E 5.8e-165、bitscore 533.6** → 真 with-lipase 落在最强桶。
5. **E 值分布（全量 1,206,655）**：`<1e-100` **1,031** / `1e-100–1e-50` **519** / `1e-50–1e-30` 2,463 / `1e-30–1e-10` 676,148 / `1e-10–1e-5` 526,223。
   → **E<1e-50 = 1,550 条；E<1e-30 = 4,013 条**（0.13% / 0.33%）。
6. **⚠️ AF2 单序列预测不可用于结构比对**：用 8YNV 自身序列做 AF2 单序列预测，再与晶体结构做 Foldseek，**TM 仅 0.20**（pLDDT 36.1）——低置信预测使结构比对失效。**这是个硬约束**。

## 2.2 更正记录（correction trace，2026-09-20）

> §2 第 2 点的原始数值「0/1,001 达到 E<1e-5」**已被证伪并更正**。原始值保留在
> `runs/20260919_phaded_pool_external_evidence_01/logs/structure_verification_recon_pilot_20260920.json`
> 的 `pilot.blastp_vs_8ynv`（已标 `superseded: true`），更正值写入同文件 `corrections[0]`。**不覆盖、不删除历史值。**

- **根因**：`plan/pilot_analyze_deferred_sample.sh` 第 33 行用 `-outfmt "6 qseqid pident length evalue bitscore"`（5 列，**index 3 = evalue，index 4 = bitscore**），但计数写成 `float(r[4])` —— **把 E 值阈值套在了 bitscore 列上**。bitscore 最小 25.8，故每个阈值都返回 0。
- **检测方式**：复核时从**同一产物**独立重算，并用不依赖主观判断的**列不变式**钉死映射：该命令带 `-evalue 1e-3`，因此真正的 E 值列必须在**每一行**都 ≤1e-3，而 bitscore 列（25.8–341）不可能满足 —— 列映射由此机器判定，不再靠肉眼读列。
- **更正值**（同产物，E 值取 index 3）：命中 **987/1,001**（E≤1e-3）；**E<1e-5 = 954**、E<1e-10 = 649、E<1e-30 = 2、E<1e-50 = 1、E<1e-100 = 1，最小 **1.16e-122**、中位 1.6e-12；pident 中位 24.254%、最大 54.138%；bitscore 四分位 **44.3 / 52.8 / 62.0**，最大 341。
- **⚠️ 口径警告（关键）**：该 blastp 用 `-subject` 只比对 **1 条序列**（8YNV 链 A，298 aa），E 值按**单序列搜索空间**计算，**比 GTDB 规模搜索宽松约 8.4×10⁵ 倍**（E 随 subject 残基数线性增长，log10 = 5.92）。这些 E 值**不可**与 hfam_2 的 HMM E 值或任何全库 E 值直接比较。换算到 GTDB 规模（2.5e8 残基）：E<1e-5 ↔ 每 1,001 条中 **578** 条，E<1e-10 ↔ **174** 条，E<1e-30 ↔ 1 条。
- **抽样设计已机器判定**（`plan/verify_pilot_blastp_columns.py`，只读，全 PASS 才退出 0）：样本取自 deferred TSV（1,206,655 行，按 genome accession 排序）的**等步长系统抽样**——首行 index 0、末行 1,206,000、**步长恒为 1,206**（1,001 行 × 1,206 ≈ 全表，**无跳区**）。
- **抽样代表性已独立验证**：1,001 条样本在群体中按 `discovery_best_evalue` 的排名 **中位 605,238（≈第 50 百分位）**，进入 top-1000 的仅 1 条（均匀抽样期望 0.83 条）；**文件行序与 E 值轴独立**（Spearman(行序, E 值秩) = **0.0074**，故步长未与 E 值对齐）；accession 前缀 `GCA_000` 占比 **样本 0.0030 / 群体 0.0026**（135/194 个前缀覆盖）→ 样本**代表群体**，**不是**高端尾部抽样，也**不是**早期登入号偏倚。**结论：该样本的 E 值/命中率统计可外推到全层**。同时复核了群体分布 E<1e-50 = 1,550 / E<1e-30 = 4,013（与记录一致）。
- **更正后的科学解读**（操作结论方向保留，**证据替换**）：
  1. deferred 层**并非**与 with-lipase 参考「无序列相似性」——代表性成员的 bitscore 中位 52.8，属**长而高度分歧的 twilight-zone 比对**；「有命中者中 96.7% 落在 E<1e-5」（单序列空间）。
  2. 因此**对单条参考做 blastp 不能作为判别性过滤器**（几乎全部通过最宽松阈值），且**无命中绝不能记为阴性**（远缘同源会被漏掉）。
  3. 换算到 GTDB 规模后仅收窄 ~1.7 倍（1,001→578，E<1e-5）/ ~5.8 倍（1,001→174，E<1e-10）→ **hfam_2 严格 E 值漏斗（E<1e-50 = 1,550 / E<1e-30 = 4,013）仍是唯一有效预筛**。
- **本层任何计数均未改变**；deferred 层仍完全隔离，未删除、未入候选计数。
- **复现与回归**：`python runs/20260919_phaded_pool_external_evidence_01/plan/verify_pilot_blastp_columns.py`（只读，全 PASS 才退出 0）；回归测试 `pipeline/tests/test_pilot_blastp_columns.py`（9 项；**先红后绿**——更正前红色态为 4 项「更正记录」断言失败，5 项列映射/统计断言通过）。日志：`logs/verify_pilot_blastp_columns.log`、`logs/pytest_test_pilot_blastp_columns.log`。
- **复发提示**：同一列错位在本任务中**已发生两次**（父代理 §11.13 point 6；本子任务的记录）。今后凡 blastp 自定格式产物，**先跑上述不变式检查再解读**。

## 3. 方案（三段漏斗）

```
deferred 1,206,655 (E<1e-5)
  └─ Stage 1 严格 E 值预筛（hfam_2 discovery_best_evalue）
       E<1e-50 → 1,550 条   ← 推荐（保守）
       E<1e-30 → 4,013 条   ← 宽口径备选
  └─ Stage 2 结构预测（仅 Stage 1 子集；GPU0）
  └─ Stage 3 Foldseek：8YNV(链A) 为 query vs 预测结构库
       --tmscore-threshold + 预注册判据
  └─ Stage 4 通过者 → 证据层/校准 gate（仍 candidate-only，不入现有计数）
```

**关于序列层**：§2.2 已证明「blastp vs 单条 8YNV」不可用作 Stage 1 之前的过滤器（既无判别力，也无命中≠阴性）。若日后要加序列层，唯一站得住的形式是**数据库规模搜索**（query = deferred 序列 vs 大型参考库，并使用同一搜索空间下的 E 值/bitscore 阈值），且必须先预注册阈值、单独授权。

**脚本位置**（全部脚本文件 scp 执行，禁止内联 ssh）：`runs/20260919_phaded_pool_external_evidence_01/plan/`
- `recon_structure_verification.sh`、`recon_structure_step2.sh`、`recon_final_tools.sh`
- `pilot_extract_deferred_sample.sh`、`pilot_analyze_deferred_sample.sh`、`pilot_step3_motif_models.sh`
- `pilot_esmfold_foldseek.sh`、`pilot_test_colabfold.sh`、`pilot_test_colabfold2.sh`、`pilot_foldseek_validate.sh`
- `verify_pilot_blastp_columns.py`（§2.2 更正的可复现校验脚本，只读）

## 4. 预注册判据（提议，需操作者确认后才执行）

| 项 | 阈值 | 依据 |
|---|---|---|
| 结构通过 | Foldseek **qtmscore ≥ 0.5** 且 **E ≤ 1e-3** | TM≥0.5 = 同一折叠的通用判据；8YNV 自测 1.00 |
| 预测质量门 | 预测 **pLDDT ≥ 70** 才纳入比对 | 试点证明 pLDDT 36 的预测 TM 仅 0.20，会污染结论 |
| fail-closed | 无结构 / pLDDT < 70 / 无 Foldseek 命中 → **not_assessed**（**不得**记为阴性） | 与项目既有 fail-closed 口径一致 |
| 单序列 AF2 结果 | **不得**作为通过依据（仅可作探索） | 试点实测同序列 TM 0.20 |
| 计数隔离 | 通过者**不进入**任何候选/高可信度计数，单列 `structural_validation` 层 | AGENTS.md 禁 deferred 入计数 |

## 5. 算力估算

| 步骤 | E<1e-50（1,550） | E<1e-30（4,013） | 全量 1.2M |
|---|---|---|---|
| 序列提取（100 分片扫描） | ~5 min | ~5 min | ~5 min |
| 结构预测 · AF2 单序列 84 s/条 | ~36 GPU-h（~1.5 天） | ~94 GPU-h（~4 天） | **~324 GPU-日（不可行）** |
| 结构预测 · ESMFold（若修好，~5 s/条） | **~2.2 GPU-h** | ~5.6 GPU-h | ~70 天（仍不可行） |
| 结构预测 · AF2+MSA（需联网 MSA，~2–5 min/条） | ~1–2 GPU-日 | ~3–6 GPU-日 | 不可行 |
| Foldseek 比对 | 分钟级 | 分钟级 | — |

> 线程/GPU 规则：AGENTS.md 单任务 ≤40 线程、`min(40, nproc-load-10)`、保 ~10 核余量。本任务为 GPU 计算，建议**只用 GPU0**（GPU1 被他人占用），CPU 侧 ≤8 线程（IO/后处理）。

## 6. 需要操作者授权的具体动作

1. **修复 ESMFold 环境**（`pip install openfold` + 装 GPU 版 torch）——可将 Stage 2 从 ~36 GPU-h 降到 ~2.2 GPU-h，**但属"安装"操作，AGENTS.md 要求明确授权**。
2. **或**批准用 ColabFold 的 **`deepfold_v1`**（需先测其速度/精度）或 **AF2+MSA**（需联网 MSA）。
3. **批准 Stage 2 计算规模**：建议先跑 **E<1e-50 的 1,550 条**（保守、~36 GPU-h AF2 / ~2.2 GPU-h ESMFold）。
4. **确认预注册判据（§4）**，特别是「pLDDT ≥ 70」的预测质量门与 TM ≥ 0.5 的通过阈值。
5. **GPU1 占用者确认**：GPU1 已被他人用到 23.8 GB，本任务只用 GPU0 不影响对方。

## 7. 未执行说明（诚实声明）

- **未做全量 1.2M 结构预测**（不可行）。
- **未做 MSA 模式的 AF2**（未获授权联网跑 MSA；且需先确认吞吐）。
- **未修复 ESMFold env**（安装操作需授权）。
- 结构通过的候选**尚未产生任何计数**——deferred 层依旧完全隔离。
