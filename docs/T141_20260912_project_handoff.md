# PHB-GTDB / PhaDED 项目交接文档

**交接日期：** 2026-09-12  
**工作区：** `<REPO_ROOT>`  
**分支：** `main`  
**当前 HEAD：** `2e44aafe2c91263447f37bac625f97c56ddafc05`  
**远程：** `https://github.com/anzhiw2-gif/PHB_gtdb-ds`  
**当前项目状态：** PhaDED candidate-only 计算阶段基本完成；ePhaZ 四种子替代模型校准失败于证据不足；正式 registry 修改和 GTDB 全量重扫未启动。

## 0. 新会话先读什么

按以下顺序读取，不要先运行新扫描：

1. `AGENTS.md`
2. 本文件 `docs/T141_20260912_project_handoff.md`
3. `docs/superpowers/plans/2026-09-09-phaded-ephaz-iphaz-next-stage.md`
4. `docs/T141_20260909_phaded_candidate_mapping_01_status.md`
5. `docs/T141_20260911_phaded_priority_review_03_status.md`
6. `docs/T141_20260912_phaded_structure_evidence_review_01_status.md`
7. `docs/T141_20260912_phaded_candidate_architecture_review_01_status.md`
8. `docs/T141_20260912_phaded_competition_phylo_review_01_status.md`
9. `runs/<run_id>/input_contract.json` 和对应结果 manifest

旧的 `docs/T141_20260909_project_handoff.md` 保留为历史交接，不再作为本阶段最新状态来源。

## 1. 不可突破的治理边界

- 所有 HMM、Pfam、InterPro、SignalP、motif、邻域、结构和系统发育结果都只能表示候选同源或功能潜力，不能表示已验证 PHB/PHA 降解表型。
- 不得把 `unassigned_PhaDED_like` 当作生物学阴性。
- 不得把 `PhaC-vs-PhaZ` challenge 候选加入 ePhaZ 阳性训练集。
- 不得把高结构 TM-score、pLDDT、树支持度或 SignalP 结果用于覆盖 Pfam 架构冲突。
- 不得修改 `pipeline/config/formal_scan_models.tsv`、`ePhaZ_curated_core`、历史 HMM、历史运行目录或正式扫描结果。
- 不得启动新的 GTDB 全量扫描，除非先完成独立 registry decision、获得明确授权、创建新的 dated run 和 dated deploy。
- T141 服务器仅可执行对应 `deploy/<run_id>/` 中绑定的源码；服务器线程不得超过 40。
- 保留所有失败重试、overlength exclusion、challenge、negative 和 superseded 结果；不要清理。
- commit、push、安装软件、改变服务器配置、删除文件均需要新的明确授权。

## 2. 主计划任务状态

主计划已更新复选框：`docs/superpowers/plans/2026-09-09-phaded-ephaz-iphaz-next-stage.md`。

| 任务 | 状态 | 说明 |
|---|---|---|
| Task 0 冻结权威输入 | **完成，待 provenance 对账** | 三个 dated run、输入合同、哈希和服务器 deploy 已建立；两个旧 contract 状态仍需单独 amendment |
| Task 1 PhaDED reference ledger | **完成** | 38 个 family 定义、723 条参考序列、ledger/FASTA 校验通过；仅 30 条实验阳性，693 条 annotation-only |
| Task 2 profile | **完成** | 46 条 profile 记录；9 个训练 HMM，37 个 `reference_only` |
| Task 3 candidate mapping | **完成** | 109,087 条候选全部 reconciliation；保留 ambiguous/unassigned |
| Task 4 独立架构/特征证据 | **完成主流程，详细人工复核仍在继续** | Pfam、InterPro、结构、完整性、系统发育已整合；仍有候选 hold 需要逐条解释 |
| Task 5 四种子 ePhaZ 校准 | **完成，结论为不足** | 三个 subtype 均为 `reference_only_insufficient_panel`，没有新模型 |
| Task 6 下游候选解释 | **部分完成** | ranking、summary、代表序列、竞争树、结构审查完成；Figure-only 修订和实验验证层未完成 |
| Task 7 registry/rescan | **未启动/授权阻塞** | 没有 formal-change plan、formal run 或 deploy；这是预期状态 |

计划中的 `[x]` 代表相应执行证据已存在；Task 5 的“完成”不代表模型成功升级。

## 3. 已完成工作与权威结果

### 3.1 PhaDED reference ledger

目录：`runs/20260909_phaded_reference_evidence_01/`

权威文件：

- `inputs/phaded_family_definitions.tsv`
- `inputs/phaded_reference_ledger.tsv`
- `inputs/phaded_reference.faa`
- `results/reference_validation.json`
- `results/profile_build_01/profile_manifest.tsv`
- `results/profile_build_01/profile_audit.json`

验证结果：

- 参考记录：723
- family 定义：38/38
- 实验阳性：30
- annotation-only：693
- 重复 accession：0
- validation failures：0

解释：38 个 family 是第二层参考词汇，不要求每个 family 都有独立 HMM。没有足够实验阳性支持的 family 必须保持 `reference_only`。

### 3.2 Profile 构建

`profile_audit.json` 的权威摘要：

- profile records：46（38 family + 8 superfamily）
- trained profiles：9
- reference-only profiles：37

训练 profile 只包括有足够训练支持的范围。不得为了覆盖 38 个 family 而强行扩展 profile。

### 3.3 109,087 条候选映射

目录：`runs/20260909_phaded_candidate_mapping_01/`

权威文件：

- `input_contract.json`
- `results/phaded_profile_scores.tsv`
- `results/phaded_assignment.tsv`
- `results/phaded_assignment_evidence.tsv`
- `results/phaded_feature_annotation.tsv`
- `results/mapping_summary.json`
- `results/mapping_reconciliation.json`
- `results/hmmsearch_provenance.json`

统计：

| 状态 | 数量 |
|---|---:|
| `assigned` | 67,749 |
| `ambiguous_family` | 9,664 |
| `ambiguous_superfamily` | 42 |
| `unassigned_PhaDED_like` | 31,632 |
| 总候选 | 109,087 |

服务器执行：HMMER 3.4，profile sequential，`--cpu 40`，使用 dated deploy。`unassigned_PhaDED_like` 仅表示当前训练 profile 未报告分数，不是阴性。

### 3.4 Pfam/InterPro/结构/系统发育整合

主要目录：

- `runs/20260911_phaded_pfam_architecture_01/`
- `runs/20260911_phaded_interpro_priority_02/`
- `runs/20260911_phaded_priority_review_03/`
- `runs/20260911_phaded_competition_phylo_review_01/`
- `runs/20260912_phaded_structure_prediction_01/`
- `runs/20260912_phaded_structure_evidence_review_01/`
- `runs/20260912_phaded_candidate_architecture_review_01/`

19 条高信息候选的结构证据整合：

- 5 条 `high_confidence_candidate_only`
- 3 条 `architecture_conflict_manual_review`
- 9 条 `intermediate_structure_manual_review`
- 2 条 `low_confidence_or_truncation_hold`

结构预测：19/19 候选完成，95 个 PDB；T141 使用 GPU 0，CPU 线程不超过 40。

竞争系统发育：4 个 panel，共 42 条候选：

- 4 条 explicit extracellular
- 12 条 PhaC-vs-PhaZ challenge
- 8 条 structure anomaly
- 18 条 ePhaZ competition

Foldseek：7 组 accession-bound 结构比较，服务器使用 4 线程。结果只支持结构同源/架构证据，不改变 Pfam 冲突状态。

### 3.5 关键候选决策

- 3 条 dPHASCL2-like：`P0_architecture_conflict_hold`
- 2 条 source/gene-model 或低结构质量 hold
- 14 条保留为 bounded comparison candidates
- 12 条 PhaC-vs-PhaZ：始终保持 challenge/unresolved
- 18 条 ePhaZ competition：按结构域和树证据保留不确定性

特别注意：3 条 dPHASCL2 候选虽然有高 pLDDT、较强树支持和 Foldseek 相似性，但缺失 dPHASCL2 参考指纹，不能升级为 curated positive。

### 3.6 四种子 ePhaZ 校准

目录：`runs/20260909_ephaz_core_calibration_01/`

权威文件：

- `inputs/subtype_training_panels.tsv`
- `results/calibration/leave_one_out_metrics.tsv`
- `results/calibration/subtype_model_decision.json`
- `results/ephaz_seed_audit.tsv`
- `results/calibration_status.md`
- `results/sha256_manifest.tsv`

结果：

| subtype | train | held-out | genera | negative/challenge hit | decision |
|---|---:|---:|---:|---:|---|
| `e_dPHASCL_type1` | 2 | 5 | 7 | 0 | `reference_only_insufficient_panel` |
| `e_dPHASCL_type2` | 1 | 0 | 1 | 0 | `reference_only_insufficient_panel` |
| `e_dPHAMCL` | 1 | 1 | 2 | 1 | `reference_only_insufficient_panel` |

`Q84C08` 是实验阴性 challenge control，却出现 `negative_hit=true`。这说明当前四种子体系不具备足够特异性，不能升级为广谱分类器。

没有修改正式 registry、历史 HMM 或正训练集。

## 4. 尚未完成的任务

### 4.1 Provenance reconciliation

不要覆盖旧 contract；创建新的 amendment/reconciliation 输出，记录：

1. `runs/20260909_phaded_reference_evidence_01/input_contract.json` 仍显示 `planned_candidate_only`，但 reference validation 和 profile build 已执行。
2. `runs/20260911_phaded_competition_phylo_review_01/input_contract.json` 的 execution 字段仍是 `tree: not_run`，但 `results/phylogeny/` 和 `results/phylogeny/phylogeny_review_manifest.json` 显示树已完成。

需要说明原始合同不可追溯修改，amendment 才是后续使用时的当前解释。

### 4.2 19 条候选的剩余人工结构复核

优先顺序：

1. 3 条 dPHASCL2 PF00756 架构冲突候选；
2. `GCA_030946875.1|JALYXY010000192.1_6` 的 C-terminal partial source record；
3. `GCF_006442545.2|NZ_CP083952.1_5059` 的低结构质量 ORF；
4. 4 条最高 outside-domain low-pLDDT burden 候选；
5. 其余 5 条 intermediate structure 候选。

对应状态：

- `docs/T141_20260912_phaded_candidate_architecture_review_01_status.md`
- `docs/T141_20260912_phaded_dphASCL2_conflict_review_01_status.md`
- `docs/T141_20260912_phaded_intermediate_structure_review_01_status.md`
- `docs/T141_20260912_phaded_hold_source_retrieval_01_status.md`

### 4.3 Figure-only 修订

尚未创建新的 visualization-only run。后续修订必须新建 dated run，不得覆盖 `runs/20260907_progress_report_01/`。

必须修正：

- Figure 4 的单位改为 `genome-family assignments`；
- 明确标记 `Palsa-465` 为 GTDB placeholder，或按记录规则替换为正式属名；
- 图注保留 GTDB 纳入、HMM 检出和注释偏差说明；
- 所有图继续标注 candidate-only 边界。

### 4.4 实验验证层

尚未完成 candidate-specific wet-lab validation plan。后续至少要为每个候选记录：

- accession、克隆蛋白和信号肽/成熟蛋白构建；
- 聚合物底物种类、组成和物理状态；
- 表达、纯化和定位条件；
- 反应条件、时间和温度；
- 空载体、热灭活和阴性对照；
- 产物检测方式；
- 生物学重复和判定标准。

实验结果必须只改变对应候选自己的 experimental-evidence 字段，不能自动验证全部同源序列。

### 4.5 正式 registry/rescan

以下文件不存在，说明 Task 7 尚未启动：

- `docs/superpowers/plans/2026-09-09-phaded-formal-registry-change.md`
- `runs/20260909_phaded_formal_scan_01/`
- `deploy/20260909_phaded_formal_scan_01/`

在新增跨菌属、分型明确的实验阳性序列并重新完成 challenge 校准前，不应创建这些文件。

## 5. 下一会话的执行顺序

### 第一步：只读状态刷新

```powershell
Get-Content -Raw AGENTS.md
Get-Content -Raw docs/T141_20260912_project_handoff.md
git status --short
git branch --show-current
git rev-parse HEAD
```

### 第二步：建立 provenance reconciliation

只读检查以下输入和输出，创建新的 dated reconciliation run 或结果 amendment；不要覆盖旧 contract：

```powershell
Get-Content -Raw runs/20260909_phaded_reference_evidence_01/input_contract.json
Get-Content -Raw runs/20260909_phaded_reference_evidence_01/results/reference_validation.json
Get-Content -Raw runs/20260909_phaded_reference_evidence_01/results/profile_build_01/profile_audit.json
Get-Content -Raw runs/20260911_phaded_competition_phylo_review_01/input_contract.json
Get-Content -Raw runs/20260911_phaded_competition_phylo_review_01/results/phylogeny/phylogeny_review_manifest.json
```

### 第三步：完成 19 条候选的结构/域人工复核

使用已有 PDB、Pfam/InterPro 坐标和系统发育表。不要重新预测全部结构，不要把结构结果写进 ePhaZ 训练集。

### 第四步：准备实验优先级表

优先考虑结构和架构均无硬冲突、序列完整、跨证据一致的 dPHASCL1-like 候选；3 条 dPHASCL2 conflict 和 12 条 PhaC-vs-PhaZ challenge 不得直接进入阳性表达清单。

### 第五步：补充实验阳性并重新校准

只有取得至少 3 个独立训练阳性、至少 3 个属、独立 held-out、预先冻结阈值和零未解释 challenge/negative hit 后，才可考虑 subtype candidate model。

### 第六步：之后才讨论 registry

先写三方案 decision record：

1. 保留当前 registry，仅发布 candidate-only PhaDED mapping；
2. 增加窄范围 subtype model，同时保留历史模型；
3. 用有版本的新模型替代旧模型。

在 decision record 审查和新的明确授权之前，不得修改 `pipeline/config/`、创建正式扫描 run 或执行服务器全库扫描。

## 6. 验证命令

当前代码验证结果：

```text
python -m unittest discover -s pipeline/tests -v
208 passed, 1 skipped
python -m compileall -q pipeline
git diff --check
```

新会话修改代码后至少重新运行：

```powershell
python -m unittest discover -s pipeline/tests -v
python -m compileall -q pipeline
git diff --check
```

## 7. 当前 Git 与发布状态

- 当前分支：`main`
- 未执行 commit、push 或 fetch
- 工作区包含大量未提交和未跟踪的项目脚本、测试、运行状态文档和生成结果
- 不得使用 `git clean`、`git reset --hard` 或删除命令清理工作区
- 发布前必须单独核对本地、GitHub、dated deploy 和 run manifest 的权威关系
- 当前 registry SHA-256：`8a1c05085f882daad9fddb9cf7616def74438d3a7bf55f219e0b2046adc5b7a3`

## 8. 一句话交接结论

PhaDED 的 38-family 参考词汇、9 个可训练 profile、109,087 条候选映射、Pfam/InterPro/结构/系统发育证据和四种子 ePhaZ 校准均已完成候选级审计；当前最重要的后续工作是完成 provenance 对账、收口 19 条高信息候选的人工结构复核、形成实验验证优先级，并补充足够的独立实验阳性序列，而不是修改 registry 或启动 GTDB 全量重扫。
