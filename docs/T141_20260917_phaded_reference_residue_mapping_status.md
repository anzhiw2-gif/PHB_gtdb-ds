# T141 PhaDED 参考残基/判别域映射补全 + lid 去循环化（Task 2）

日期：2026-09-17
Run：`runs/20260917_phaded_reference_residue_mapping_01`
计划：`docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` §3 Task 2
脚本：`pipeline/scripts/complete_phaded_reference_residue_mapping.py`；测试：`pipeline/tests/test_complete_phaded_reference_residue_mapping.py`
边界：本轮全部本地完成；没有 ssh、没有连接服务器、没有提交任何服务器任务；没有修改 `pipeline/config/formal_scan_models.tsv`、历史 HMM、`ePhaZ_curated_core`、任何历史 run 或 `results/`；没有删除任何文件；没有执行任何 git 写操作。

> **candidate-only 声明：** 本项目的 HMM、profile、domain、motif、SignalP、结构与系统发育结果只表示候选同源或功能潜力，**不等同于已验证的 PHB/PHA 降解表型**。本 Task 产出的参考残基映射、判别域绑定与校准报告同样只是参考层的序列/结构证据，不构成任何表型结论，也不是催化机制验证。

---

## 0. 一句话结论

参考层 motif 面板从"lipase box + oxyanion hole"升级为**有出处的催化三联体 + 判别域面板**：His 从 **0/723 → 90/723**，Asp 从 **25/723 → 71/723**，SBD（PF06850 绑定）**0 → 274**，lid 从**subtype 先验赋值 → 独立序列检出 39 条**；同时诚实保留 633 条 His、650 条 Asp、467 条 oxyanion、682 条 linker 为 `pending_reference_annotation`。**本 Task 不产生任何新的 family 或 superfamily 判定。**

> **口径补充（2026-09-17，见 §8）：** 上面 His 90 与 Asp 71 是 **resolved（锚定列传递 + 序列模式）** 的总数，其中 **45 / 45** 条只来自保守度门槛下的锚定比对列传递，属推断而非独立检出；独立检出为 **45 / 26**。校准报告原字段名 `independent_detection_supported` 承载的是总数（不准确），本次已新增字段语义澄清文件并修正脚本（原报告文件保持不变）。

---

## 1. 命中率：补全前 vs 补全后（参考集 723 条）

"before" 逐列实测自 `runs/20260915_phaded_motif_reconciliation_01/results/reference_motif_panel.tsv`（只读），"after" 来自本 run 的 `results/reference_mapping_manifest.tsv`。表内数字全部由 `results/implementation_calibration_report.json` 的 `criteria.*.before/after` 直接给出，未手工转录。

| 判据列 | before supported | after supported | 变化 | after 剩余 pending |
|---|---:|---:|---:|---:|
| `lipase_box_state` | 441 | **444** | +3 | 0 |
| `catalytic_ser_cys_state` | 443 | **444** | +1 | 273 |
| `his_state` | **0** | **90** | **+90** | 633 |
| `asp_state` | 25 | **71** | +46 | 650 |
| `oxyanion_hole_state` | 263（另 conflict_relative_position 7） | **256** | −7 | 467 |
| `ahsmg_state` | 2 | **2** | 0（回归保护通过） | 721 |
| `sbd_state` | 0 | **274** | +274 | 0（449 条 `not_detected_in_tested_pfam`） |
| `linker_state` | 0 | 0（41 条为**候选区间**，未做 profile） | 0 | 682 |
| `lid_state` | 0（665 `not_expected_for_subtype` + 58 `pending_reference_annotation`） | **39**（另 10 条 `partial_lid_loop1_only`） | +39 | 0（674 条 `not_detected_pattern`） |

补充口径（`pattern_layer_hit_rates`，分母均为 723）：| 层 | 命中 | 命中率 |
|---|---:|---:|
| lipase box 序列模式 `Gx1Sx2G` | 444 | 0.6141 |
| His 序列模式（8DAJ 校准）`GM[A-Z]H` | 42 | 0.0581 |
| Asp 序列模式（8DAJ 校准）`G[A-Z]{2}DYTV` | 25 | 0.0346 |
| oxyanion 序列模式（8DAJ 校准）`HGC[A-Z]Q` | 268 | 0.3707 |
| AHSMG 字面模式 | 2 | 0.0028 |
| lid 双环序列检出 | 39 | 0.0539 |
| lid 单环（仅 loop1） | 10 | 0.0138 |
| SBD `PF06850` 绑定 | 274 | 0.3790 |
| 任一 residue 走了锚定列传递 | 64 | 0.0885 |
| His 最终 resolved（锚定列 + 模式） | 90 | 0.1245 |
| Asp 最终 resolved（锚定列 + 模式） | 71 | 0.0982 |

> 上表两行"最终 resolved"是**总数**：His 90 中 42 条为序列模式命中、3 条为锚行自身文献坐标、**45 条为锚定列传递**；Asp 71 中 23 / 3 / **45**；oxyanion 256 中 255 / 1 / **0**。独立检出（不含传递）因此是 His **45**、Asp **26**。字段语义与逐判据拆分见 §8 与 `results/implementation_calibration_report_field_notes.json`。

**口径变化说明（不是生物学变化）：** 旧实现的 `lipase_box_state` 有 278 条被先验直接判为 `not_expected_for_subtype`（"无 lipase box"家族），不再扫描。新实现把先验移入独立列 `lipase_box_expectation_from_superfamily`（`not_expected_for_superfamily` 276 条 / `expected_for_superfamily` 447 条），状态本身改为**独立检出**，因此 441 → 444：多出的 3 条正是"无 lipase box"家族中偶然带 `Gx1Sx2G` 的序列。同理 `lid_state` 的先验移入 `lid_expectation_from_superfamily`（`expected_for_superfamily` 52 条 / 其余 671 条）。这两列**只用于事后对照**，不参与任何检出。

---

## 2. 缺陷修正（明确记录）

### DEFECT A：`lid_state` 原本由 subtype 先验派生，构成循环论证（本任务的核心修正）

**原状：** `audit_phaded_motifs.py::_architecture_state` 直接由 `subtype` 字符串赋值——`not_expected_for_subtype`（665 条）与 `pending_reference_annotation`（58 条）。即：用 subtype 派生的字段去检验 subtype，无论数据如何都会"自证"。

**修正：** `detect_lid(sequence, signature=None)` 只接受序列与一个 motif 签名，函数体内**不出现** `subtype` / `superfamily` / `phaded_family` 任何标签；先验被移到独立列 `lid_expectation_from_superfamily`，只用于**事后对照**，不参与检出。该不变量由两条测试硬性强制：

- `test_lid_detector_input_excludes_subtype_and_superfamily_labels`：`inspect.getsource(detect_lid)` 断言禁用词不出现，且签名参数 ⊆ {`sequence`, `signature`}；
- `test_lid_detection_call_sites_never_pass_a_label`：断言所有 `detect_lid(` 调用点行不含标签字段。

**独立检出的结果（未使用任何标签）：**

- 检出 39 条，**全部**落在 `intracellular nPHAMCL`；对先验"预期无 lid"的其余超家族 0 检出；
- 阴性对照 `Q51718.1`（PhaZGK13，胞外 mclPHA 解聚酶）→ `not_detected_pattern`；
- 阳性对照 `AAM63408.1`（PhaZKT，胞内 mclPHA 解聚酶）→ `supported`，`lid_coordinates=34-38;190-195`，`lid_loop_spacing=156`（文献值 190-34=156）；
- 先验对照统计（`lid_prior_agreement`）：`detected_and_expected=39` / `detected_but_not_expected=0` / `not_detected_but_expected=13`。**"detected_but_not_expected=0"是本条最关键的诚实数字**：独立检出没有出现任何与文献先验冲突的阳性。

**结论表述：** 修正后 lid 才是独立判据，**才可用于区分胞内与胞外 mcl 酶**（Task 5 判别轴"定位"依赖此列）。同时必须写明其**不完整性**：先验预期有 lid 的 52 条中只检出 39 条（13 条未通过双环判据，其中 10 条仅有 loop1），因此 lid 是"特异性好、敏感性不足"的判据：`supported` 可作证据，`not_detected_pattern` **不得**读作"无 lid 的结构性阴性"。

### DEFECT B：His 模式从未触发（0/723）——根因是文献consensus 与其自身示例互不一致

**根因（本轮实测）：** 既有 `HIS_PATTERN = GM[A-Z]H[A-Z]{2}P[A-Z]{2}G` 要求第 7 位为 Pro。而 LtPHBase 原文（Protein Sci 2022, DOI 10.1002/pro.4470）写的是"the histidine motif, **GMXHXXPXXG**, is found at amino acids **267–274**（**GMGHAWSSG**）"：所给示例 `GMGHAWSSG` 只有 9 个残基、第 7 位是 Ser 而非 Pro；所述区间 267–274 只有 8 个残基；两者都无法匹配其自述的 10 残基 consensus。**该文献 consensus 与其自身示例互相矛盾**，任何照抄都必然 0 命中。

**实测证据：** 测试 `test_legacy_his_pattern_never_fires_on_any_reference` 断言旧模式在 723 条参考上命中列表为**空**（作为缺陷的可复现回归保护）。

**修正：** 改用**已被坐标验证的观测示例**派生模式 `GM[A-Z]H`（8DAJ His270 位于 `GMGHAWSGG`，P12625.1 的同一模式命中其 `GMAHAWPAG`），并在实现里做 fail-closed 校准：模式必须在 8DAJ 上**精确命中 270 号位**，否则 `calibrate_patterns` 直接抛错、不出结果。校准记录见 `results/implementation_calibration_report.json::pattern_calibration`（`state=calibrated_hits_literature_coordinate`）。

### DEFECT C：三联体判据只对 e-dPHASCL type 1 生效

`GxxDYTV`（Asp）与 `GMxH`（His）都是 type 1 序列特征；type 2 是**环状置换**折叠（Hisano 2006），胞内 mcl 酶的 His 落在 `DDGHL`、胞外 mcl 酶落在 `GGHEW`，都不被任何 `GMxH` 模式覆盖。本 Task 用**锚定列传递**解决可传递的部分，并用**保守性闸门**拒绝不可传递的部分（见 §4）。

### DEFECT D：SBD 从未绑定，且 PF06850 不是胞外 SBD（重要反向发现）

原 `sbd_state = pending_tool_input`（717 条）。按计划要求"用绑定代替新写 motif"，本轮把参考层 hmmscan 输出 `results/raw/reference_pfam.domtblout` 逐条绑定（i-Evalue ≤ 1e-5）：**274 条**得到 `supported`。

**但必须同时记录一个与计划 §3 Task 1 措辞相冲突的实测事实：** `PF06850`（PHB_depo_C）在本参考集中**只出现在 `intracellular nPHASCL without lipase box`（274/276）**，在 `extracellular dPHASCL type 1`（284 条）与 `type 2`（71 条）中**0 命中**。因此：

- `PF06850` 绑定**不足以**代表胞外酶的底物结合域判据；manifest 用独立列 `sbd_role_caveat` 显式标注这一点；
- 胞外超家族的 SBD 判据在本轮**仍未解决**（`not_detected_in_tested_pfam` 不等于"SBD 结构性缺失"）；
- 该结论**只**陈述本参考集与本次 Pfam 子集扫描的观测，未做新的域注释。

### DEFECT E：AHSMG 坐标曾被当作催化残基坐标上报

旧实现把 `catalytic_ser_cys_state=supported` 与 `catalytic_ser_cys_coordinates` 一并取自 AHSMG 匹配（`catalytic_state = "supported" if ahsmg`、`catalytic_coord = _coord(ahsmg)`），等于把"替代催化基序"的位置当成亲核残基位置上报（影响 `2VTVA`、`AAK07742.1` 两条）。现已分开：AHSMG 只出现在 `ahsmg_state/ahsmg_coordinates`；两条 PhaZ7 参考的 lipase box 判为 `not_detected_pattern`（它们本来就没有 lipase box）。这解释 `catalytic_ser_cys_state` 443 → 444 的差值（+3 为"无 lipase box"超家族中偶然的 `Gx1Sx2G`，−2 为上述误标移除）。

---

## 3. 文献残基证据表（版本化）

文件：`runs/20260917_phaded_reference_residue_mapping_01/inputs/literature_residue_reference.tsv`（20 行 × 26 列；size 15306；SHA-256 `387c2641b3d1…`，完整值见 `input_contract.json`）

每行记录：`accession / reference_id / organism / residue_role / reported_coordinate / alignment_coordinate / coordinate_offset / coordinate_source / source_database / source_version / retrieval_date / coordinate_retrieval_method / evidence_type / doi / pmid / pmcid / sequence_length / anchor_sequence_sha256 / coordinate_verified_in_sequence / verification_method / notes`。**未知一律写 `pending`，不推断。**

| 锚点 | 类型 | 文献坐标 | 序列中坐标 | 来源 |
|---|---|---|---|---|
| PDB **8DAJ** chain A（LtPHBase，302 aa，sha256 `1e3bcaf7…`） | 晶体结构 | Ser121 / His270 / Asp197 / oxyanion Cys40 | 同（offset 0） | Protein Sci 2022;31(11):e4470，PMID 36222314，PMC9601781，DOI 10.1002/pro.4470 |
| PDB **2D80** chain A（PfPHBase，318 aa，sha256 `3f23a6f8…`） | 晶体结构 | Ser39 / Asp121 / His155 / oxyanion Cys250 | 19 / 101 / 135 / 230（**offset −20**） | J Mol Biol 2006;356:993，PMID 16405909，DOI 10.1016/j.jmb.2005.12.028（坐标取自该文摘要的 ESTHER 数据库记录） |
| `AAM63408.1`（PhaZKT，285 aa，sha256 `001aa9c3…`） | 定点突变 + 结构 | Ser102 / Asp221 / His248；lid loop1 FNGIG 34–38、loop2 YYWQLF 190–195；lid 域 133–192 | 同（offset 0） | Appl Microbiol Biotechnol 2025;109(1):215，PMID 41055782，PMC12504323，DOI 10.1007/s00253-025-13605-z |
| `Q51718.1`（PhaZGK13，278 aa，sha256 `400721c8…`） | 结构 + 定性陈述 | Ser172 / Asp228 / His260；文献明确胞外酶**无 lid** | 同（offset 0） | 同上（原文称菌株为 *Pseudomonas solani* GK13；DED 账本沿用旧名 *P. fluorescens*，同一 GK13 菌株） |
| Knoll 2009（DED） | 文献陈述 | AHSMG（Phaz7 型替代基序，**无坐标**）；"lipase box 与催化三联体残基已在比对中人工标注" | — | BMC Bioinformatics 2009;10:89，PMID 19296857，PMC2666664，DOI 10.1186/1471-2105-10-89 |

**锚点自身校验（fail-closed，17 行全部通过）：** `results/implementation_calibration_report.json::anchor_verification` 逐行记录"坐标处的残基确实等于该 role 期望残基"；任一不符即抛错、不产出结果。附带的交叉证据：**RCSB 2D80 chain A 序列的 SHA-256 与 DED 账本中参考 `2D80A` 的 `sequence_sha256` 完全相同**（`3f23a6f8…`），即"PDB 2D80 chain A ≡ DED 参考 2D80A"是哈希级证据，不是推断。

### 3.1 未能核实的前提（诚实记录）

Knoll 2009 声称 "Residues of the lipase box and the catalytic triad were **manually annotated** … based on multiple sequence alignments"，但**库内 DED 快照的 38 个 `.aln` 文件是裸 ClustalW 块，没有任何注释行或注释列**（已逐文件检查）。因此本 Task **没有**读取 DED 人工注释，而是**用文献坐标 + 比对列做等价传递**；该未核实前提以 `LIT_KNOL2009_ALIGNMENT_ANNOTATION` 一行（`coordinate_verified_in_sequence=not_verified`）显式留痕，不写成已完成。

---

## 4. 方法：模式层 + 锚定列传递层 + 保守性闸门

1. **模式层（标签无关，全 723 条统一施加）**：4 个 role 的序列模式全部由 8DAJ 的文献坐标校准，并在实现中对 8DAJ 做精确坐标断言；AHSMG 仅作回归对照。
2. **锚定列传递层**：对 38 个 family 比对中 `usable_reference_alignment` 的 29 个，若该比对**含有带文献坐标的参考行**，则把该文献残基所在列（及其无缺口坐标）传给同比对的所有行。共 **3 个比对**满足条件：`DED_hfam_4`（锚 `AAM63408.1`）、`DED_hfam_8`（锚 `Q51718.1`）、`DED_hfam_70`（锚 `2D80A`），覆盖 64 条参考（0.0885）。锚行自身使用**其本文献坐标**，与同比对其他行区分（`literature_anchor:…` vs `anchor_column_transfer:…`）。
3. **保守性闸门（本轮新增的关键安全阀）**：传递列只有在锚定残基于该比对中**保守度 ≥ 0.5** 时才允许外传，否则该 role 对被传递行一律 fail-closed 为 `pending_anchor_column_not_conserved`，**绝不**记为 `not_detected_*`。实测保守度：

| 比对 | 锚 | role | 列 | 观测行数 | 保守度 | 可传递 |
|---|---|---|---:|---:|---:|---|
| DED_hfam_4 | AAM63408.1 | catalytic_histidine | 315 | 42 | **1.0000** | 是 |
| DED_hfam_4 | AAM63408.1 | catalytic_aspartate | 287 | 44 | 0.9545 | 是 |
| DED_hfam_4 | AAM63408.1 | catalytic_serine | 165 | 46 | 0.9565 | 是 |
| DED_hfam_8 | Q51718.1 | 三个 role | 264/232/176 | 6 | **1.0000** | 是 |
| DED_hfam_70 | 2D80A | catalytic_histidine | 150 | 13 | **0.0769** | **否** |
| DED_hfam_70 | 2D80A | catalytic_aspartate | 115 | 13 | **0.0769** | **否** |
| DED_hfam_70 | 2D80A | oxyanion_hole_cysteine | 255 | 13 | **0.0769** | **否** |
| DED_hfam_70 | 2D80A | catalytic_serine | 27 | 2 | 0.5000 | 是 |

`DED_hfam_70` 的低保守度不是脚本缺陷：该 family 的锚 2D80A（真菌）与其余 12 条细菌序列只在中后段共列，传递列上读到的是 `W/P` 等非催化残基——正是 Hisano 2006 所述**环状置换**折叠的直接体现。若不设闸门，这 12 条会被误写成"His/Asp 未检出"的假阴性；设闸门后它们保持 `pending`。该行为由测试 `test_unconserved_transferred_columns_fail_closed_to_pending` 强制。

4. **其他判据**：SBD 用 `PF06850` 逐条绑定（§2 DEFECT D 附 caveat）；linker 用"催化域（PF10503）C 端之后、下一个已注释域之前的区间"**界定候选区间**（41 条），显式标注 `not_profiled` 且**不**给 supported；`oxyanion_hole_state` 的宽松旧模式 `H[A-Z]?G[A-Z]?C[A-Z]?Q` 换成文献 consensus `HGC[A-Z]Q`（8DAJ `HGCTQ`、2D80 `HGCLQ`），type 1/2 的几何归属**按 Task 4 的 fail-closed 规则留给 Task 4**，本 run 只输出 `oxyanion_relative_position_vs_lipase_box`（upstream 191 / downstream 63 / undetermined 469）备用。

### 4.1 逐超家族映射表（`results/superfamily_residue_mapping.tsv`）

| 超家族 | 参考数 | 锚 | His 支持 | His pending | Asp 支持 | lid 支持 | lid 未检出 | PF06850 |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| extracellular dPHAMCL | 6 | Q51718.1 | 6 | 0 | 6 | 0 | 6 | 0 |
| extracellular dPHASCL type 1 | 284 | 无锚 | 42 | 242 | 22 | 0 | 284 | 0 |
| extracellular dPHASCL type 2 | 71 | 2D80A | 1 | 58 | 2 | 0 | 71 | 0 |
| extracellular native-SCL/PhaZ7-like | 4 | 无锚 | 0 | 4 | 0 | 0 | 4 | 0 |
| intracellular nPHAMCL | 52 | AAM63408.1 | 41 | 11 | 41 | **39** | 3 | 0 |
| intracellular nPHASCL with lipase box | 29 | 无锚 | 0 | 29 | 0 | 0 | 29 | 0 |
| intracellular nPHASCL without lipase box | 276 | 无锚 | 0 | 276 | 0 | 0 | 276 | **274** |
| periplasmic PHA depolymerases | 1 | 无锚 | 0 | 1 | 0 | 0 | 1 | 0 |

---

## 5. 本 Task 明确不做的事

1. **不产生任何新的 family / superfamily 判定。** 8DAJ 与 2D80 的外部锚在表中 `phaded_superfamily` 一律为 `pending`；模式层标签无关地全库施加；锚定列传递只在**既有** DED 比对内部进行。没有新增聚类、没有重新划分 family、没有改动 `phaded_family_definitions.tsv`。
2. **不改变候选层。** 本 run 只处理 723 条参考；109,087 条候选的 motif 层重算属 Task 7，需新 dated run。
3. **不把 `pending` 当阴性。** 缺失映射一律 `pending_*`，且 `supported` 必带坐标与 `mapping_basis`（由 `test_supported_states_always_carry_a_coordinate_and_a_basis` 强制）。
4. **不降低校准 gate、不启动全库重扫、不启动服务器任务。**

---

## 6. 仍未完成的部分（诚实清单）

| 缺口 | 数量 | 原因 |
|---|---:|---|
| `his_state` pending | 633 | `intracellular nPHASCL（有/无 lipase box）`、`periplasmic`、`extracellular nPHASCL`、`PhaZ7-like` 四个超家族没有任何带残基坐标的文献锚；type 2 除锚自身外列不保守 |
| `asp_state` pending | 650 | 同上；另有 2 条在锚定列读到非 D 残基（`not_detected_at_mapped_column`，证据级而非表型级） |
| `oxyanion_hole_state` pending | 467 | 3 个有锚比对中只有 `DED_hfam_70` 的锚带 oxyanion 坐标，且该列不保守 |
| `linker_state` pending | 682 | 无 linker 家族 profile 可绑定；41 条只做到"候选区间界定"，未做 profile 判定 |
| 胞外 SBD 判据 | 全部胞外参考 | `PF06850` 不覆盖胞外超家族（§2 DEFECT D），需要真正的胞外 SBD 证据来源 |
| Cys 型亲核残基 | 276 条（`intracellular nPHASCL without lipase box`） | 文献（Knoll 2009 正文）只说该家族用 Cys-His-Asp，未给可核对坐标；库内 `comp_ft_reutropha_phaz.txt` 亦无坐标，故不推断 |
| 参考行未进入任何可用比对 | 147 条（His） | 9 个 family 比对为 header-only，另有部分参考不在其家族比对行中 |

这些缺口**只影响参考层的证据密度，不影响已有结论的方向**：His/Asp 的 supported 集合全部带坐标与出处；pending 既不是阳性也不是阴性。

---

## 7. 产物与可复现命令

| 产物 | 路径 | 说明 |
|---|---|---|
| 文献残基证据表 | `runs/20260917_phaded_reference_residue_mapping_01/inputs/literature_residue_reference.tsv` | 20 行，26 列 |
| 锚点序列（含检索日期） | `runs/20260917_phaded_reference_residue_mapping_01/inputs/literature_anchor_sequences.faa` | 8DAJ / 2D80 chain A，逐序列 SHA-256 |
| 参考映射 manifest | `.../results/reference_mapping_manifest.tsv` | 723 行 × 53 列 |
| 逐超家族映射表 | `.../results/superfamily_residue_mapping.tsv` | 8 行 |
| 实现校准报告 | `.../results/implementation_calibration_report.json` | before/after、5 条缺陷、17 条锚点校验、模式校准、列保守度、hit rates、limitations |
| 输入契约 | `.../input_contract.json` | 10 个输入全部 `verified`（path/size/SHA-256）；GTDB 输入 `pending`（本 Task 不消费） |
| 取序列脚本 | `.../logs/fetch_literature_anchor_sequences.py` | 本 run 唯一的联网步骤（RCSB FASTA），带检索日期 |
| 执行脚本 | `.../logs/run_reference_residue_mapping.py` | 复现整个 run |
| 聚焦测试补跑日志（本次新增，P2） | `.../logs/step2.3_focused_tests_green_rerun.log` | `Ran 26 tests ... OK`，含全部 26 个测试方法；SHA-256 `4af8b5856a0d61e0afa7b7eda14d7e5c38650e0b75400e080c97272e0b604580` |
| 校准报告字段语义澄清（本次新增，P4） | `.../results/implementation_calibration_report_field_notes.json` | 字段语义 + 逐判据 supported/独立检出/锚定列传递拆分（从 manifest 重算） |

```powershell
python runs\20260917_phaded_reference_residue_mapping_01\logs\fetch_literature_anchor_sequences.py   # 仅联网步骤
python runs\20260917_phaded_reference_residue_mapping_01\logs\run_reference_residue_mapping.py
python -m unittest discover -s pipeline/tests -p "test_complete_phaded_reference_residue_mapping.py"
python -m unittest discover -s pipeline/tests
python -m compileall -q pipeline
git diff --check
```

**测试（先红后绿）：** 实现前运行该测试文件：5 个可独立执行的测试中 3 个失败，另 3 个测试类的 `setUpClass` 因实现文件缺失全部报错（`Ran 5 tests ... FAILED (failures=3, errors=3)`，日志 `logs/step2.2_failing_test_red.log`）；实现后 26 个测试全部通过。其中 `test_legacy_his_pattern_never_fires_on_any_reference` 永久保留为缺陷回归保护。

**（2026-09-17 更正，P2）** 上段原先引用绿灯日志 `logs/step2.3_focused_tests_green.log` 并称其记录 `Ran 26 tests ... OK`。**该引用有误，已更正**：见 §7.1。

### 7.1 绿灯日志引用更正与补跑（Phase 1 验证报告 P2）

**问题（P2）：** 本文件 §7 曾称绿灯证据为 `logs/step2.3_focused_tests_green.log` 且记录 `Ran 26 tests ... OK`。实测该日志（size 10918，SHA-256 `ec6d932f29fb2e7988a199da4ad982b472879846ad1dfcc443094797ba92e78c`，UTF-16LE 带 BOM）汇总行为 **`Ran 25 tests in 0.332s`**，且其中只出现 **25** 个测试方法名，**完全不含** `test_unconserved_transferred_columns_fail_closed_to_pending` —— 也就是 §4 第 3 步"保守性闸门"的安全阀测试。被引用的绿灯证据恰好漏掉了它。

**实测复核（本次）：**

| 项目 | 值 |
|---|---|
| 测试文件当前方法数 | **26** |
| 原日志汇总行 | `Ran 25 tests in 0.332s` |
| 原日志中出现的方法数 | **25** |
| 原日志缺失的方法 | `test_unconserved_transferred_columns_fail_closed_to_pending` |
| 日志文件 mtime | 2026-09-17 01:11:00 |
| 测试文件 mtime | 2026-09-17 01:12:27 |

**原因说明（不推断）：** 上述两个时间戳显示测试文件在原日志生成后约 77 秒被再次写入，与该测试方法当时尚未存在于文件中相容；但本次没有任何执行记录能证明该次写入就是新增该测试，**故不判定原因**，只如实记录"原日志少一条测试方法"这一事实。

**本次补跑（覆盖证据）：** 用 `python -m unittest discover -s pipeline/tests -p "test_complete_phaded_reference_residue_mapping.py" -v` 重跑，把完整输出写入**新文件** `logs/step2.3_focused_tests_green_rerun.log`（UTF-8，size 5296，SHA-256 `4af8b5856a0d61e0afa7b7eda14d7e5c38650e0b75400e080c97272e0b604580`），**原日志一字未改，继续作为历史证据保留**。补跑结果：`Ran 26 tests ... OK`，log 中出现的测试方法数 = 26 = 测试文件声明的方法数，缺 0 条、多 0 条（明细见 `runs/20260917_phaded_p2_p4_remediation_01/results/task2_focused_test_rerun.json`）。

> 注意：补跑日志是一次具体执行的原始输出，其中包含本次耗时时长，因此重跑会生成等价但不会逐字节相同的日志；上表 SHA-256 对应本文件写入的这一个实例。

**结论：** 测试实现本身没有缺陷（26/26 通过），问题只在**被引用的证据文件不完整**，现已由覆盖全部 26 个方法的补跑日志接替引用。

---

## 8. 校准报告字段命名澄清（Phase 1 验证报告 P4）

**问题（P4）：** `results/implementation_calibration_report.json` 中 `criteria.<判据>.after.independent_detection_supported` 等于该判据的 **supported 总数**，但其中一部分**不是独立检出**，而是由锚定比对列**传递**得来（保守度 ≥ 0.5 的推断）。本文件（§1、§4）的表述是诚实的（写"resolved（锚定列 + 模式）"），但机器可读字段名会被下游误读为"全部独立检出"。

**本次从 `results/reference_mapping_manifest.tsv` 实际重算的结果**（重算脚本 `runs/20260917_phaded_p2_p4_remediation_01/logs/write_calibration_field_notes.py`，未照抄本文档数字）：

| 判据 | resolved（= supported 总数） | 独立检出（不含传递） | 锚定列传递 | 其中锚自身文献坐标 | 模式/行内直接证据 |
|---|---:|---:|---:|---:|---:|
| `his_state` | 90 | 45 | **45** | 3 | 42 |
| `asp_state` | 71 | 26 | **45** | 3 | 23 |
| `oxyanion_hole_state` | 256 | 256 | **0** | 1 | 255 |
| 其余 6 个判据 | 444 / 444 / 2 / 274 / 0 / 39 | 同左 | 0 | 0 | 同左 |

定义：**独立检出** = 证据直接落在该参考行自身（校准后的序列模式命中、Pfam 域绑定、或该行本身就是锚因而带自己的文献坐标），**不依赖**把别人的比对列传过来；**锚定列传递** = 仅凭"从锚行按保守度门槛把比对列搬到该行"获得坐标，属推断。

**与验证报告 P4 所述数字的差异（如实报告）：** 验证报告写"其中 48、48、1 条来自锚定列传递"。实测：报告里 `supported_anchor_column` 字段确实是 48 / 48 / 1，但这**不是纯传递数**——其中 3 / 3 / 1 条是**锚行自身**（`AAM63408.1`、`Q51718.1`、`2D80A`）携带的本文献坐标，真正由传递得到的只有 **45 / 45 / 0**。其中 oxyanion 的传递数为 **0**：唯一带 oxyanion 坐标的锚 `2D80A` 所在 `DED_hfam_70` 列保守度仅 0.0769，未过闸门，因此该家族其余 12 条保持 `pending`（这与 §4 第 3 步一致）。

**顺带记录的既有报告内部不一致（只记录，本次不修改该文件）：** `alignment_column_transfer.references_with_anchor_column_transfer = 65`，而 manifest 中真正带锚定比对列的参考只有 **64** 行；65 = 三个锚定比对的**比对行数**之和（46 + 6 + 13），并非参考数。

**处理方式（优先不覆盖）：**

1. **原报告文件保持不变**：`results/implementation_calibration_report.json`（size 26467，SHA-256 `8a93f73cb1e66e1eb03611e2686a52e294e221b113f85ac83e602ecdd37c4213`）与 `results/reference_mapping_manifest.tsv`（size 956813，SHA-256 `38441f712e3e52e0cb550984ca442584e09f6cf8e8394dd60ff31163620bfda2`）**本次均未改动**，继续作为 Phase 1 证据。
2. **新增澄清文件**：`results/implementation_calibration_report_field_notes.json`，给出上述字段语义、逐判据拆分、与冻结报告的逐判据对账、与验证报告 P4 数字的差异说明，并推荐字段名 **`resolved_supported`**（= supported 总数）。
3. **脚本已为后续 run 修正**：`pipeline/scripts/complete_phaded_reference_residue_mapping.py` 现在输出 `resolved_supported`、`independent_detection_supported`（**已排除**锚定列传递）、`transferred_from_anchor_column`、`anchor_self_literature_coordinate` 与 `supported_evidence_breakdown`，并新增顶层 `field_semantics` 说明块；`his_state`/`asp_state`/`oxyanion_hole_state` 的 `after_note` 也把三类证据来源写清楚。改代码前先写了失败测试 `pipeline/tests/test_phaded_calibration_field_semantics.py`（首跑 6 个测试 1 失败 4 报错），实现后 8 个测试全部通过；原有 26 个 Task 2 测试仍全部通过。
4. **本 run 的 `input_contract.json` 未改动，但其 `inputs.source_script.sha256` 指向改前版本**（`112fc34603f4b66c4481b35445a8d30e65012213df32f20de2c95a51b8a6aba3`），与当前脚本文件已不一致——契约是历史记录，按"不覆盖历史证据"原则保持原样；改前/改后哈希与出处一并记入 `implementation_calibration_report_field_notes.json` 的 `script_correction.previous_sha256` 与 `runs/20260917_phaded_p2_p4_remediation_01/input_contract.json` 的 `modified_by_this_run.source_script`。

---

*本文件为 candidate-only 研究记录；所有映射、域绑定与检出结果表示候选同源或功能潜力，不等同于已验证的 PHB/PHA 降解表型。*
