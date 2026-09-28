# PHB_gtdb-ds 对话—修改—运行过程复盘文档

- 文档日期：2026-09-06
- 工作区：<REPO_ROOT>
- 当前 Git HEAD：`2e44aafe2c91263447f37bac625f97c56ddafc05`
- 当前接受版本：`runs/20260906_overall_pha_phb_archaea_audit_04`
- 文档用途：把本轮对话中提出的问题、形成的判断、实际修改、服务器运行、结果解释和剩余边界集中记录，便于后续复盘与交接。

> 重要口径：本项目输出的是基于 HMM、覆盖度、SignalP、结构域、邻域和分类信息的候选同源/功能潜力。除非另有独立实验、accession 级证据，本项目不把计算命中称为已验证的 PHB 或 MCL-PHA 降解表型。

## 1. 当前最终结论（先看这一节）

1. 当前最可靠的主线结果是“整体 PHA 候选的 ePhaZ/iPhaZ 分层”，不是 PHB 特异降解酶的最终统计。
2. 接受 run `20260906_overall_pha_phb_archaea_audit_04` 共保留 109,087 条去重蛋白候选；全部标为 `PHB_candidate_review`，`phb_supported=0`。这表示没有建立 GTDB 命中与独立实验参考 accession 的直接对应关系，不表示“没有 PHB 酶”。
3. ePhaZ/iPhaZ 的分层已经完成，但各层是证据优先级，不是实验活性等级：ePhaZ core nonsecreted 16,818、ePhaZ SCL secreted 21,856、ePhaZ tier2 review 31,990、ePhaZ competition review 18、iPhaZ intracellular 32,926、iPhaZ tier2 review 5,479。
4. 古菌 PhaZh1-like 采用独立 candidate-only 重建。2,307 条古菌记录中，197 条为 `PhaZh1_like_high`（来自 173 个基因组），但仍需审计和实验验证，不能称为 197 条古菌 PHB 解聚酶。
5. 正式 frozen scan 13 已完成，但本轮没有再次启动正式扫描，也没有修改 `pipeline/config/formal_scan_models.tsv`。当前仍未进行 Git commit 或 push。

## 2. 项目规则和证据层级

本次工作以 `AGENTS.md` 为现场规则来源，核心约束如下：

- 每次运行必须使用新的、可审计的 `runs/<run_id>/`，保留 `logs/`、`inputs/`、`results/` 和 `input_contract.json`；不得覆盖历史结果。
- 服务器只能运行 dated `deploy/<run_id>/` 中已绑定源码，不直接调用服务器项目根目录旧脚本。
- GTDB taxonomy、metadata、tree、HMM、源码和环境必须记录路径、版本、大小和 SHA-256；缺失内容写 `pending`，不能伪造哈希。
- HMM、domain、SignalP、邻域和树结果表示候选同源或功能潜力，不能替代酶活、敲除/互补、纯化蛋白底物实验等表型证据。
- 重算、安装、配置变更、删除、提交、推送和正式扫描都必须在得到明确授权后执行；本次复盘只读核对并新增文档，不执行这些动作。
- 代码改动遵循“先失败测试、后实现”，完成后运行相关测试、`compileall` 和 `git diff --check`。

## 3. 对话主线和阶段性问题

下表按实际推进顺序概括用户问题与项目动作。历史报告中的旧数字仅作为当时状态，当前事实以接受 run、manifest 和最新状态页为准。

| 阶段 | 用户关心的问题 | 实际决策和动作 | 产物/状态 |
|---|---|---|---|
| A. 项目接管 | 阅读审查报告、检查本地全部文件、GitHub 和 T141 服务器 | 核对 `AGENTS.md`、审计报告、源码、测试、run/deploy 结构和 Git 状态；服务器执行边界固定为 dated deploy | 建立“本地仓库—GitHub—T141 run/deploy”权威关系；禁止把旧快照当现状 |
| B. 扫描加速 | 用 60 线程扫描并写监控脚本 | 正式扫描任务按任务级并发运行，HMMER 使用 60 CPU；允许复用已完成输出，但必须写入新的 run 并保留 joblog/manifest | 形成 run 12/13 的可审计运行链，不覆盖历史目录 |
| C. 停止 run 12、启动 run 13 | 停止旧运行，建立 run 13、60 并发、复用完成输出 | 新建 dated run，复用已完成分片，未复用缺少绑定证据的未知输出；要求 1,000/1,000 任务标记和失败清单闭合 | `20260901_formal_frozen_scan_13` 完成；297 个任务复用，703 个实际执行，失败 0 |
| D. run 13 QA 和下游 | 检查是否正常结束、审查 raw 命中、解释 147,690 与 38,741 | 将 raw HMM 命中、registry-threshold、验证后 tier1 和历史 44,814 严格分开；不混用不同输入快照和定义 | raw accepted 6,740,900；严格四家族 union 38,741 genomes；历史 44,814 仅作 Scheme A 比较值 |
| E. ePhaZ 初步治理 | 旧模型是否能证明 PHB；MCL-PHA 漏检是什么 | 识别 ePhaZ HMM 的底物范围并非 PHB 专一；把 PHB、MCL-PHA、胞内动员、非靶标和 challenge 分开 | 形成独立阳性/阴性面板和 candidate-only 规则 |
| F. MCL-PHA 边界和快速优化 | 补充跨属实验阳性、快速调阈值、是否正式扫描 | 增加 Bdellovibrio、Burkholderia、Streptomyces 的完整实验支持序列；固定 profile，做 LOO、coverage/E-value 和 challenge 校准；divergent profile 仍是 singleton，拒绝注册 | `20260902_ephaz_mcl_crossgenus_candidate_07`、`subfamily_candidate_12`、`threshold_candidate_14`；未改正式 registry，未启动正式扫描 |
| G. 家族标准梳理 | 为什么不同家族 HMM 阈值不同；OH、ArchPhaZ、PhaJ/BdhA/PhaC/phasin 如何解释 | 接受“家族特异阈值”是正常做法，要求每家族记录校准依据、用途和证据边界 | OH 保持 E<=1e-5、min_cov=0.6；ArchPhaZ 按 >=0.8/0.6–0.8/<0.6 分层；背景家族不计入解聚酶 union |
| H. OH/ArchPhaZ 审计 | OH 文件损坏、ArchPhaZ 是否真是古菌扩张 | 先保留零字节损坏证据，再从同一 run tier 文件恢复可验证的 `OH_validated.faa`；跨门类抽查；不把低覆盖古菌命中称为古菌扩张 | `20260903_oh_arch_crossphylum_audit_03`；OH 高层全为细菌，ArchPhaZ >=0.8 全为细菌；古菌命中集中低覆盖层 |
| I. 古菌模型重建 | 古菌序列为什么能用于搜索细菌；古菌 PhaZ 如何处理 | 将 Haloferax PhaZh1 作为独立 patatin-like 候选家族，不并入细菌 ePhaZ/iPhaZ PHB 主结论；方案 B 固定 8 条参考序列，建一次 HMM，保留 negative panel | `20260904_archaea_phazh1_candidate_01`；模型 candidate-only |
| J. ePhaZ/iPhaZ 分层 | 是否已经真正分层，能否先做亚型 | 复用 run-13/tier 和 SignalP/competition 证据，按 (family, accession) 去重；保留 secreted、nonsecreted、tier2、competition、intracellular 层 | `20260905_ephaz_iphaz_layering_01`；protein 109,087，去重后无重复键 |
| K. 整体 PHA→保守 PHB→古菌审计 | 先做整体 PHA，后单独挑 PHB；审计 197 条古菌高层 | 新增 overall audit 脚本和测试；保留原始层级，增加 `pha_subtype_evidence`、`phb_focus`、`phb_exclusion_reason` 和古菌 motif/邻域审计字段 | 接受版本为 `audit_04`；PHB 不作强行升级，古菌 197 条全部 review_required/candidate-only |
| L. 当前复盘 | 梳理所有提问、修改和结果，判断是否过度治理 | 本文统一记录事实、解释、pending 和下一步；不删除异常残留、不提交、不推送 | 本文是当前对话过程文档 |

## 4. ePhaZ/iPhaZ：到底改了什么

### 4.1 模型和证据分层

早期 ePhaZ 结果把较宽的同源空间与核心候选混在一起。后续采取 split registry：

- `ePhaZ_curated_core`：实验支持更明确、结构和长度更符合核心定义，用于严格候选层。
- `ePhaZ_broad_discovery`：保留远缘、注释型、短序列和可能的 MCL-PHA 相关空间，用于 discovery，不直接并入严格 PHB/解聚酶统计。
- `iPhaZ`：单独作为胞内候选家族，并保留与 ePhaZ 的竞争审计。
- ePhaZ SignalP 结果进一步区分 secreted、lipoprotein、Tat 和 no-signal；SignalP 只支持定位证据，不证明底物。
- 无信号肽 ePhaZ 与 iPhaZ 的竞争命中保留 `competition_review`，不强行归类。

### 4.2 主要校准经验

- 初始小面板显示 ePhaZ/iPhaZ 在 E<=1e-5 下可与部分阴性分开，但外部面板确认 ePhaZ 不能据此获得 PHB 底物专一性。
- 代表性 MCL-PHA 非 PHB 对照 Q84C08 可被 ePhaZ 命中，证明 ePhaZ fold 空间与胞外 PHA 解聚相关，但不能单独证明 PHB。
- 跨属 MCL-PHA 阳性加入后，经典分支可恢复，lipase-associated 和 Streptomyces 分支的单 profile 泛化失败；因此不能仅靠放宽阈值解决漏检。
- coverage >=0.8 在有限 challenge 面板上有用，但不足以替代独立阳性、阴性和跨属验证；divergent 子 profile 在正式注册前仍需要额外 accession 级阳性。

### 4.3 实际代码修复

本轮接受脚本 `pipeline/scripts/overall_pha_phb_archaea_audit.py` 处理了三类已发现问题：

1. 历史标签 `MCL_like` 与新标签 `MCL_lipase_associated` 的兼容映射。
2. `ePhaZ_tier2_review` 与 `ePhaZ_competition_review` 不再错误压成 `unresolved_pha_evidence`。
3. 字符串 `nonsecreted` 含有子串 `secreted`，原判断会错分；现已调整判断顺序并加回归测试。

## 5. run 11/12/13 与下游关系

### 5.1 失败的 run 11

HMMER 3.4 对目标序列长度超过 100,000 aa 有工具限制。run 11 因此失败；这不是生物学阴性。失败任务、stderr 和部分构建证据必须保留，不能删除或改写为成功。

### 5.2 run 12 到 run 13

run 12 是旧启动快照；随后按用户要求停止旧方案并新建 run 13。run 13 使用 60 线程/任务级并发，复用已有完成输出，但每个复用项仍需记录来源和 manifest。最终：

- 1,000/1,000 任务完成；
- 297 个任务复用父运行输出；
- 703 个任务实际执行；
- 失败任务为 0；
- `hits_all.tsv` 约 6,743,198 行（含表头）；
- 4 条超长序列写入 `overlength_exclusions.tsv`，属于工具限制，不是阴性。

### 5.3 run 13 的三个统计层

- raw HMM 命中：接近 6.74M 行，包含多个蛋白/同一基因组多次命中。
- registry-threshold accepted：6,740,900 条，适用于描述注册表阈值命中。
- tier1：经过序列完整性、验证、覆盖度或家族规则后的高可信计算候选；四家族 union 为 38,741 个基因组。

147,690 是严格四家族 union 之前的 registry-threshold 基因组层，38,741 是后续 tier1/验证层；两者不可混用，也不等同实验阳性。历史 44,814 来自旧 Scheme A 定义，必须保留为历史比较值。

## 6. 各家族当前口径

| 家族 | 当前处理 | 可以说什么 | 不能说什么 |
|---|---|---|---|
| ePhaZ curated/broad | 分开保存；E<=1e-5；暂不设统一 coverage 硬门槛 | PHA 解聚酶样候选、分泌/非分泌候选 | 不能直接称 PHB 特异酶；broad 不可并入严格主统计 |
| iPhaZ | 独立胞内层；保留 competition review | 胞内 PHA 动员/解聚相关候选 | 不能与 ePhaZ 合并成同一证据层 |
| OH | E<=1e-5、min_cov=0.6；>=0.8 可作为高可信优先层 | 经过覆盖度门槛的 OH 候选 | 不能称所有命中实验活性；0.8 不是通用阈值 |
| ArchPhaZ_hydrolase | >=0.8 high、0.6–0.8 review、<0.6 exploratory | 跨域 hydrolase 候选；高层优先复核 | 不能把 12,469 tier1 基因组称为古菌 PHB 解聚酶扩张；高层 499 个基因组全为细菌 |
| ArchPhaZ_patatin | broad patatin/lipid-enzyme fold | 邻域联合的 patatin 候选 | 不能单独证明 PHB 降解 |
| PhaJ | 背景通路指标 | PHB/PHA 动员相关背景 | 不是直接解聚酶阳性 |
| BdhA | 背景代谢家族 | 代谢背景/邻域字段 | 普遍存在，不能单独证明降解 |
| PhaC | PHA 合成/颗粒背景 | 合成相关共现证据 | 不能单独证明降解 |
| phasin | 辅助颗粒蛋白；当前命中很少 | 可作为辅助字段和灵敏度警告 | 不能作为阳性或阴性门槛 |

不同家族阈值不同是正常的，因为 profile 长度、保守区、远缘程度、近邻假阳性和目标用途不同。统一要求不是“统一数字”，而是每家族都要有可追溯阈值、校准集、用途和证据边界。

## 7. 古菌 PhaZh1-like 重建和 197 条审计

### 7.1 为什么古菌序列可以用来搜索细菌（以及为什么本项目没有这样混用）

古菌参考序列可用于建立跨域同源探针或发现保守 fold，但跨域命中不能自动继承古菌实验表型，也不能证明细菌命中同样降解 PHB。Haloferax 的实验依据支持的是 patatin-like PhaZh1 这一古菌候选路线；因此本项目将古菌模型独立保存，不并入细菌 ePhaZ/iPhaZ 的 PHB 主结论。

### 7.2 重建动作

- 固定 8 条古菌 patatin-like 参考序列，构建一次 `archaea_PhaZh1_like_candidate.hmm`。
- 使用支持面板、广谱/冲突 patatin 对照和邻域联合分类，输出 high/review/exploratory/non-PhaZh1 四层。
- 来源服务器为 T141（交互使用 `<SERVER_USER>@<SERVER_HOST>`；不记录任何密钥）；执行源必须是 dated deploy `20260904_archaea_phazh1_candidate_01`。
- 全部结果标记 candidate-only；未修改正式 registry，未启动 GTDB 正式扫描。

### 7.3 197 条 high 的实际审计

总古菌记录 2,307 条：

- `PhaZh1_like_high`：197 条，173 个基因组；
- `PhaZh1_like_review`：7 条；
- `archaeal_patatin_exploratory`：2,027 条；
- `non_PhaZh1_patatin`：76 条。

197 条 high 的审计字段：

- 195 条为 `Halobacteriota; Halobacteria`，2 条为 `Thermoproteota; Nitrososphaeria`；
- FASTA 和 taxonomy 均为 197/197 匹配；
- N 端前 80 aa 的 `GxSxG`：188 条 present、9 条 absent；
- 邻域标记：BdhA 192 条、BdhA+PhaJ 4 条、PhaJ 1 条；
- 独立 domain annotation：全部 `pending`。

因此正确表述是“197 条值得进一步审计/实验验证的古菌 patatin-like 候选”，而不是“197 条已确认古菌 PHB 解聚酶”。古菌候选不进入细菌四家族 PHB 主 union。

## 8. 接受 run 04 的实现、输入和输出

### 8.1 新增源码与测试

- `pipeline/scripts/overall_pha_phb_archaea_audit.py`：整体 PHA 层保留、保守 PHB focus、古菌全量审计、manifest 输出。
- `pipeline/tests/test_overall_pha_phb_archaea_audit.py`：覆盖原始层级保持、PHB 保守分类、历史 MCL 标签、review 层、nonsecreted 判定、古菌排除、FASTA accession、motif 字段和空 results 目录。

### 8.2 接受 run 目录

`runs/20260906_overall_pha_phb_archaea_audit_04/` 已包含 `inputs/`、`logs/`、`results/` 和 `input_contract.json`。

主要输出：

- `results/overall_pha_phb.tsv`：109,087 条去重 ePhaZ/iPhaZ 候选；
- `results/archaeal_phazh1_full_audit.tsv`：2,307 条古菌审计记录；
- `results/analysis_manifest.json`：输入哈希、计数和 phenotype boundary；
- `input_contract.json`：输入来源及 SHA-256 绑定。

当前接受 run 的输入哈希：

| 输入 | SHA-256 |
|---|---|
| protein_layers.tsv | `f4a9d99039ef5dc46fbea5b90dd98fd61adc2c65ae1da700c53976cb5b7df08c` |
| archaeal_candidate_layers.tsv | `276ca77fdfc2a64581bf1a0ee759995a23694527bb6001309d0f90b9f5949dc6` |
| tier1_genome_family.tsv | `35e1f69bca364a0fbc693b0b1abb92a45a765b610bbfe20057be17b0c1964e59` |
| archaeal_cluster_context.tsv | `db4e0f1810ae26dbacf9f8470bab986df9503c26f88e76dbbb11d8a01b704658` |
| ArchPhaZ_patatin_archaea.faa | `6391767db390ea99622ef3188b4447f26f915e2f1ebb316faeaf00e0cc07bd93` |
| dual_profile.tsv | `e08d264efcbb0890bd16fcac7818fbd1641b687c31c052e2fb05eb3109f7144f` |

### 8.3 接受 run 和关键文件哈希

- `overall_pha_phb_archaea_audit.py`：`5F2BE0122C506E3C1756FFBA3918063075CDED54E3C5D84F3B32352A8C019EFA`
- `test_overall_pha_phb_archaea_audit.py`：`B3C5AC2E4AD06CE91148E53BB9489FD602CB45EF72250C14EBB232D52DC82EA4`
- `formal_scan_models.tsv`：`8A1C05085F882DAAD9FDDB9CF7616DEF74438D3A7BF55F219E0B2046ADC5B7A3`
- `overall_pha_phb.tsv`：`2BDBD706A951D3985D422190326F995C2547CD9C87BFEBBD29BF48196F624080`
- `archaeal_phazh1_full_audit.tsv`：`3E546B8C1A10A5701156C1DCBF20C727B0D3D92D6B33604A5FF81400F021EDD6`
- `analysis_manifest.json`：`287B991F3B70D2FA9D6D0D941CDB9BBE9401493675A8BEE1EE93B7F4F919AC9B`
- `input_contract.json`：`09055C7EFC51DC0BA1391BE0B192B371C1433449B82F271437D67F2CA6527754`

### 8.4 验证

最近一次验证：

- 测试：170 passed，1 skipped，6 个 subtests passed；
- `python -m compileall pipeline/scripts`：通过；
- `git diff --check`：通过；
- `formal_scan_models.tsv` 未修改；
- 正式 GTDB scan 未启动。

## 9. 历史 run 关系和必须保留的证据

所有历史 run 都有追溯价值，不应因“结果已被替代”而删除。尤其包括：

- `20260831_formal_frozen_scan_11`：失败证据；超长输入属于工具边界。
- `20260901_formal_frozen_scan_13`：正式 60 线程扫描完成结果。
- `20260902_formal_scan13_tier_processing_02`：run 13 tier1/tier2 下游处理。
- `20260903_run13_qa_02`：运行闭合与计数 QA。
- `20260903_oh_arch_crossphylum_audit_03`：OH 恢复和跨门类抽查。
- `20260904_archaea_phazh1_candidate_01`：古菌候选模型和四层分类。
- `20260905_ephaz_iphaz_layering_01`：ePhaZ/iPhaZ 去重分层。
- `20260906_overall_pha_phb_archaea_audit_01` 至 `_04`：整体 PHA/PHB/古菌审计的迭代证据；`_04` 为接受版本，`_01`–`_03` 保留错误和修复轨迹。
- MCL candidate runs `20260901_ephaz_phb_boundary_candidate_02`、`20260902_ephaz_dual_profile_candidate_06`、`20260902_ephaz_mcl_crossgenus_candidate_07`、`20260902_ephaz_mcl_subfamily_candidate_12`、`20260902_ephaz_mcl_threshold_candidate_14`：记录边界样本、跨属阳性和未注册 profile 的原因。

## 10. 是否存在过度治理：复盘判断

### 10.1 有充分理由的治理

- 早期将 HMM 命中直接接近 PHB 表型，确实存在科学解释风险；分开 PHB/MCL-PHA、胞内/胞外和背景家族是必要的。
- OH 的覆盖度门槛有直接的近邻假阳性依据；ArchPhaZ 的跨门类分层有助于避免把低覆盖古菌命中解释成扩张。
- run 目录、输入契约、复用来源、失败证据和模型哈希绑定是大规模计算可复现性的必要条件。
- 古菌 PhaZh1 与细菌 hydrolase 不是同一证据链，独立建模和 candidate-only 标记是合理的。

### 10.2 可能偏重、且已产生边际收益递减的部分

- ePhaZ/MCL 在 2026-09-02 前后进行了多轮 boundary、dual-profile、cross-genus、subfamily、coverage 和 threshold candidate run。前几轮解决了真实的 MCL 漏检和非 PHB 边界问题；在发现 divergent profile 仍为 singleton、缺少独立 accession 后，继续反复调阈值的收益已经很低。
- coverage 0.8 在有限面板上有效，但不能通过继续加严阈值“制造”跨属泛化；需要新的独立实验阳性，而不是更多同一面板搜索。
- 已接受的 run 13 和 audit_04 已足以进入下游候选汇总、人工抽查和论文准备；再次重建全部 HMM 或正式全库扫描没有当前证据必要性。

### 10.3 结论

治理不是全部多余，但应在当前节点收束：保留现有分层和候选证据，停止没有新失败模式支持的重复阈值循环，把资源转向少量高价值人工/实验验证。

## 11. 当前 pending、out-of-scope 和未完成项

### pending

- 接受 run 的古菌记录独立 domain annotation 尚未完成。
- 109,087 条 GTDB 候选没有 accession 级独立实验确认，因此 PHB-specific support 仍为 0。
- 古菌 197 条 high 的真实 PhaZh1/PHB 酶活尚无对应实验验证。
- 全量 ePhaZ/iPhaZ 重新建树和 HGT 分析仍未完成或被暂停；已有树部分存在 stale/input_not_registered 状态。
- Scheme A 与 run 13 的 accession 级重叠在 QA 中为 pending，不能自行推断。

### out-of-scope（本次没有做）

- 未修改 `formal_scan_models.tsv`。
- 未启动新的正式 GTDB 扫描。
- 未删除损坏证据、历史 run、未跟踪文件或异常残留 `--amino`、`=.8=%d`。
- 未执行 Git commit、push、GitHub Release 或 Zenodo 发布。
- 未把古菌命中并入细菌 PHB 主统计，也未把 PhaJ/BdhA/PhaC/phasin 计作直接解聚酶阳性。

### 未消除的 warning

- `OH.faa` 原始损坏文件无法由 tier1/tier2 完全重建，零字节状态作为证据保留；只恢复了可验证的 `OH_validated.faa`。
- raw、accepted、tier1、genome union 和历史 Scheme A 数字定义不同，后续报告必须同时写明层级。
- 服务器上任何新执行都必须重新绑定 dated deploy、源码和输入哈希，不能依赖项目根目录旧脚本。

## 12. 推荐的低治理推进顺序

1. 以 audit_04 和 run-13 tier 结果作为当前冻结候选基线，不再反复调 ePhaZ 阈值。
2. 完成古菌 197 条的独立 domain annotation，并抽取少量跨门类 high/review 记录做人工邻域和结构域复核。
3. 对 ePhaZ/iPhaZ、OH、ArchPhaZ_hydrolase 只做分层后的去重、门/门类统计和候选优先级，不把它们直接写成 PHB 表型。
4. 为 PHB 特异性建立新的独立 accession 级实验面板；达到跨属、非同源重复和零负对照误报后，再讨论 profile 注册。
5. 需要重新计算时，新建 dated run，先写 input contract 和失败测试；完成后再由用户单独授权正式扫描、registry 变更、commit/push 或清理残留。
6. 发布前同步 README/STATUS、registry、run manifest 和 GitHub，明确“候选层”与“实验验证层”的差异。

## 13. 证据索引

- 当前接受状态：[T141_20260906_overall_pha_phb_archaea_audit_04_status.md](T141_20260906_overall_pha_phb_archaea_audit_04_status.md)
- 古菌候选重建：[T141_20260904_archaea_phazh1_candidate_01_status.md](T141_20260904_archaea_phazh1_candidate_01_status.md)
- ePhaZ/iPhaZ 分层：[T141_20260905_ephaz_iphaz_layering_01_status.md](T141_20260905_ephaz_iphaz_layering_01_status.md)
- OH/ArchPhaZ 抽查：[T141_20260904_oh_arch_crossphylum_audit_status.md](T141_20260904_oh_arch_crossphylum_audit_status.md)
- run 13 QA：[T141_20260903_run13_qa_02_status.md](T141_20260903_run13_qa_02_status.md)
- run 13 下游：[T141_20260902_formal_scan13_tier_processing_02_status.md](T141_20260902_formal_scan13_tier_processing_02_status.md)
- 家族特异标准：[T141_20260903_family_specific_hmm_standards.md](T141_20260903_family_specific_hmm_standards.md)
- 项目历史审计：[PROJECT_AUDIT_REPORT_20260901.md](PROJECT_AUDIT_REPORT_20260901.md)
- 最新主状态：[STATUS.md](STATUS.md)
- 接受 run manifest：[../runs/20260906_overall_pha_phb_archaea_audit_04/results/analysis_manifest.json](../runs/20260906_overall_pha_phb_archaea_audit_04/results/analysis_manifest.json)
- 接受 run 输入契约：[../runs/20260906_overall_pha_phb_archaea_audit_04/input_contract.json](../runs/20260906_overall_pha_phb_archaea_audit_04/input_contract.json)

## 14. 交接检查清单

- [x] 当前权威 run、manifest、状态报告已确定。
- [x] raw/accepted/tier1/历史统计已区分。
- [x] ePhaZ/iPhaZ、OH、ArchPhaZ、PhaJ、BdhA、PhaC、phasin 的证据边界已记录。
- [x] 古菌 197 条 high 的审计结果和 pending domain annotation 已记录。
- [x] 错误修复、测试、compileall 和 diff 检查已记录。
- [x] 未提交、未推送、未启动正式扫描已明确记录。
- [ ] 古菌独立 domain annotation 和新实验面板：后续工作。
- [ ] 用户明确授权后再进行正式扫描、registry 变更、commit/push 或清理残留。

## 15. 种子、HMM 和 PHA/PHB 查询边界核查（2026-09-06）

### 15.1 核查范围与本次动作

本节是对现有种子、HMM 构建脚本、正式模型注册表及已接受结果的只读核查记录。未重建 alignment 或 HMM，未修改 `pipeline/config/formal_scan_models.tsv`，未启动新的 GTDB 扫描，未提交或推送。本节的模型解释遵守项目证据边界：HMM 命中、结构域、信号肽、邻域和系统树均表示候选同源或功能潜力，不能单独证明体内 PHB/PHA 降解表型。

### 15.2 原始种子清单的实际构成

核查对象为 `pipeline/seeds/seeds_manifest.tsv`。其中有 78 条唯一 accession：55 条为 reviewed、23 条为 unreviewed；64 条标为 train、14 条标为 validation。按当前注释归类，ePhaZ 相关 10 条、iPhaZ 相关 17 条、OH 33 条、BdhA 15 条、phasin 3 条；分类域为 Bacteria 63 条、Eukaryota 14 条、Archaea 1 条。

这不是一个只由“已实验证实 PHB 解聚酶”组成的窄种子集。`pipeline/seeds/seeds_annotation.md` 还明确记录了 13 条真核 BDH/BDH2 同源物和 5 条 nylon oligomer hydrolase；它们属于边界、近邻或对照信息，不能叙述为 PHB 解聚酶种子。原始收集脚本 `pipeline/scripts/02_collect_seeds.py` 通过 UniProt REST 的 EC、蛋白名称、`phaZ` 和指定物种等查询词检索，因此该种子池从设计上就是多家族、混合功能空间的参考池。

`pipeline/config/params.yaml` 中的 “curated 78” 应视为文档表述问题：78 是原始 manifest 的总条目数，不能等同于 55 条 reviewed 或已人工确认的核心种子数。本次仅记录该问题，未修改配置；后续如更新该文档，需单独获得授权并同时说明“raw manifest”“reviewed”和“实验支持核心集”的不同口径。

### 15.3 HMM 的构建、版本和权威关系

历史构建脚本 `pipeline/scripts/04_build_hmms.sh` 采用 CD-HIT 95% 去冗余、MAFFT 多序列比对和 HMMER `hmmbuild`；`pipeline/scripts/04b_build_hmms_v2.sh` 改为 CD-HIT 90%，当序列数大于 1,500 时可额外按 80% 处理，再进行 MAFFT 和 `hmmbuild`。因此 `data/hmms/` 与 `data/hmms/v2/` 中共存的资产是历史/开发版本，不能仅凭文件名替代当前正式扫描来源。

正式扫描模型的权威入口是 `pipeline/config/formal_scan_models.tsv`，其中登记 `ePhaZ_curated_core`、`ePhaZ_broad_discovery`、`iPhaZ`、`OH`、`BdhA`、`ArchPhaZ_patatin`、`ArchPhaZ_hydrolase`、`PhaJ`、`phasin` 和 `PhaC`。ePhaZ 的拆分来源于 `runs/20260828_ephaz_split_01/input_contract.json`：curated core 使用 4 条种子，broad discovery 使用 4,454 条种子；HMM 及其 alignment 的 SHA-256 已在该输入契约中绑定。二者必须分开报告，broad discovery 不得并入严格核心统计。

### 15.4 各模型回答的生物学问题

| 模型/标记 | 查询用途 | 可作出的保守表述 | 不可作出的表述 |
|---|---|---|---|
| `ePhaZ_curated_core` | 严格的胞外 PHA 解聚酶样候选层 | 胞外 PHA-depolymerase-like 候选 | PHB 特异酶或已验证 PHB 降解表型 |
| `ePhaZ_broad_discovery` | 远缘、注释、短序列及 MCL-PHA 相关的发现层 | 广义远缘 PHA 相关候选 | 严格主统计或 PHB 特异性证据 |
| `iPhaZ` | 胞内 PHA 动员/解聚相关层 | 胞内 PHA mobilization/depolymerase-like 候选 | 胞外降解或与 ePhaZ 等价的同一机制 |
| `OH` | 3HB oligomer hydrolase 相关候选 | 经 `E<=1e-5`、`min_cov=0.6` 校准的 OH 候选 | 全部命中均有 PHB 解聚酶活性；其近邻 nylon hydrolase 混淆已消除 |
| `BdhA`、`PhaJ`、`PhaC`、phasin | 代谢、动员、合成或颗粒背景 | 通路背景或邻域的联合证据 | 单独作为直接解聚酶阳性 |
| `ArchPhaZ_patatin` | 广义 patatin-fold 候选 | patatin-like 候选，需结合古菌审计 | 已证实古菌 PHB 解聚酶 |
| `ArchPhaZ_hydrolase` | 跨域 hydrolase-like 探索层 | 跨域 hydrolase 候选 | 古菌 PHB 解聚酶扩张的证据 |

### 15.5 本项目究竟在查询 PHA 还是 PHB

PHA 是一类聚羟基脂肪酸酯聚合物的总称；PHB 是其中由 3-hydroxybutyrate 单体构成的短链 PHA。当前模型组合查询的是广义 PHA 相关的同源/功能空间，并按胞外 ePhaZ、胞内 iPhaZ、OH 及背景和古菌探索层分层；它不是一个已验证的 PHB 特异酶模型集合。

因此，`runs/20260906_overall_pha_phb_archaea_audit_04` 的 109,087 条记录保留为 `PHB_candidate_review`，并不等于 PHB 阴性或 PHB 酶不存在；`phb_supported=0` 的含义仅是尚未为这些 GTDB accession 建立独立、accession 级的底物/实验关联。要形成 PHB 特异子集，必须另行建立带 accession、底物组成、实验类型和文献证据的 PHB 面板，并以跨属的独立阳性和近邻阴性进行校准；不能从 HMM 命中直接推断。

### 15.6 支持本次解释的项目记录和文献

- 项目记录：`docs/项目详细审查报告.md`、`docs/T141_20260903_family_specific_hmm_standards.md`、`docs/T141_20260904_archaea_phazh1_candidate_01_status.md`、`docs/T141_20260905_ephaz_iphaz_layering_01_status.md`、`docs/T141_20260906_overall_pha_phb_archaea_audit_04_status.md`。
- Knoll et al. 2009, *PhaDED*, BMC Bioinformatics 10:89, PMCID: PMC2666664：PHA depolymerase 家族的系统性归类资源，支持不把多样 PHA 解聚酶压缩为单一 PHB 特异类别。
- Jendrossek and Handrick 2002, *Annual Review of Microbiology*：PHA 代谢与解聚酶综述，支持区分胞外聚合物降解与胞内颗粒动员。
- Jendrossek et al. 1995, *Paucimonas lemoignei* 胞外解聚酶，DOI: `10.1128/jb.177.3.596-607.1995`：胞外 PHB 解聚酶的实验锚点之一。
- Papaneophytou et al. 2009, *Thermus thermophilus* HB8，DOI: `10.1007/s00253-008-1842-2`：胞内 PHA 动员/解聚相关证据，不能与胞外路线混为一类。
- Liu et al. 2015, archaeal `PhaZh1`，DOI: `10.1128/AEM.04269-14`：支持 Haloferax 的 patatin-like 古菌候选路线；不能外推为当前所有 `ArchPhaZ_hydrolase` 命中均具有古菌 PHB 解聚表型。
