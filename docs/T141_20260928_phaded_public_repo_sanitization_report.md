# 公开仓库身份信息净化报告（F16 前置）

- **任务**：`docs/superpowers/plans/2026-09-28-phaded-evidence-model-redesign-followup.md` F16「研究快照入版本控制」的前置净化部分
- **范围**：只做「协调者提交之前必须完成」的两件事 —— 修正公开仓库安全测试的缺陷、净化一切会进入版本控制的文件
- **本文件性质**：普通文档（`docs/` 下），**不是 run 产物**，不写入 `runs/`
- **编写日期**：2026-09-28
- **未执行的动作**：本任务**没有**运行 `git add` / `git commit` / `git push` / `git reset` / `git checkout`，**没有写索引、没有提交**；`git diff --check` 通过，`git status` 中暂存区为空

---

## 1. 安全测试缺陷与修复

### 1.1 缺陷

`pipeline/tests/test_public_repo_safety.py` 用下面的写法构造 Windows 路径模式：

```python
r"(?i)" + re.escape("D:" + "\\\\PHB_gtdb-ds")
```

`"\\\\"` 这个**字面量**在 Python 里求值为**两个**反斜杠字符，`re.escape` 又把这两个字符各自转义，于是最终正则要求**四个**连续反斜杠：

| 阶段 | 取值 |
|---|---|
| 源码文本 | `"\\\\"`（4 个反斜杠字符） |
| 字面量求值后 | 2 个反斜杠字符 |
| `re.escape` 之后 | `\\\\`（正则里的 4 个反斜杠字符） |
| 最终正则 | `(?i)D:\\\\PHB_gtdb\-ds` |

因此该模式**只能**匹配被转义写坏的双反斜杠形式，**匹配不到真实存在的单反斜杠路径**；前斜杠写法同样匹配不到。`C:` 那条规则（`C:" + "\\\\Users\\\\HUAWEI"`）有完全相同的缺陷。

### 1.2 实测证据（旧模式漏报真实命中）

对 `git show HEAD:<file>` 的原始内容求值：

| 文件 | 旧模式 | 新模式 |
|---|---|---|
| `docs/T141_20260917_project_handoff.md` | **NO MATCH（漏报）** | MATCH |
| `docs/T141_20260921_phaded_comprehensive_review.md` | **NO MATCH（漏报）** | MATCH |
| `docs/T141_20260928_phaded_evidence_run_authorization_packets.md` | **NO MATCH（漏报）** | MATCH |

合成样本对照（`BS` = 一个反斜杠字符）：

| 样本 | 旧模式 | 新模式 |
|---|---|---|
| `D:<BS>PHB_gtdb-ds`（单反斜杠） | **False** | **True** |
| `D:<BS><BS>PHB_gtdb-ds`（双反斜杠） | True | True |
| `D:<FS>PHB_gtdb-ds`（前斜杠） | **False** | **True** |
| 4 个反斜杠 | False | False（不是本仓库的实际拼写） |

与此同时，远端用户 `@` 形式与远端根路径两条模式一直正常工作 —— 它们没有走「反斜杠字面量 + `re.escape`」这条路径，**本次修复没有削弱它们**。

### 1.3 修复后的模式集合

模式集合被提取为**模块级常量** `IDENTITY_PATTERNS`，并由 `identity_patterns()` 与 `identity_violations(text)` 暴露，使回归测试直接验证**真实模式**而不是副本。

修复办法是用**正则字符类** `[\\/]{1,2}` 表达「1 个反斜杠、2 个反斜杠或 1 个前斜杠」：

| 项 | 修复前 | 修复后 |
|---|---|---|
| 本地仓库根 | `r"(?i)" + re.escape("D:" + "\\\\PHB_gtdb-ds")` | `r"(?i)" + re.escape("D:") + r"[\\/]{1,2}" + re.escape("PHB_gtdb-ds")` |
| 本地家目录 | `r"(?i)" + re.escape("C:" + "\\\\Users\\\\HUAWEI")` | `r"(?i)" + re.escape("C:") + r"[\\/]{1,2}" + re.escape("Users") + r"[\\/]{1,2}" + re.escape("HUAWEI") + r"(?![0-9A-Za-z_])"` |
| 远端根路径 | `r"(?i)(?<!\$\{PHB_REMOTE_ROOT\})" + re.escape("/home/data/" + "hao"+"yu")` | 追加 `(?<![\\/])` 前瞻保护与 `(?![0-9A-Za-z_])` 边界 |
| 其余 5 条 | — | **逐字未改**（`<SERVER_HOST>`、`<SERVER_USER>@`、私钥头、GitHub token） |

修复过程中新增的回归测试**立刻抓到两个真实的过度匹配**，并已一并收紧：

- `C:<BS>Users<BS>HUAWEI2` 被旧写法误判 → 加 `(?![0-9A-Za-z_])`
- `${PHB_REMOTE_ROOT}2/other` 被误判 → 加 `(?![0-9A-Za-z_])`

收紧只增加前瞻，**不减少任何真实命中**：`${PHB_REMOTE_ROOT}` 后面跟反引号、右括号、斜杠、反斜杠、行尾等所有真实拼写仍然命中（回归测试中逐个断言）。

### 1.4 新增回归测试

`IdentityPatternRegressionTests` 断言（7 个测试全部通过）：

- `D:<BS>PHB_gtdb-ds`（单反斜杠）**必须命中**
- `D:<BS><BS>PHB_gtdb-ds`（双反斜杠）**必须命中**
- `D:<FS>PHB_gtdb-ds`（前斜杠）**必须命中**
- `C:<BS>Users<BS>HUAWEI`（单反斜杠）**必须命中**
- 占位符 `<REPO_ROOT>`、`<LOCAL_HOME>`、`<SERVER_USER>`、`<SERVER_HOST>`、`${PHB_REMOTE_ROOT}` **不得命中**
- 无关词 `D:<BS>projects<BS>other`、`C:<BS>Users<BS>PUBLIC`、`10.16.1.142`、`<SERVER_USER>2@example.org` **不得命中**（防过度匹配）
- 已净化文本（占位符混排）**不得命中**（防止净化反被自己判违规）

同时把 `git ls-files` 调用改为 `git -c core.quotePath=false ls-files -z`：原写法在本机 `core.quotePath` 生效时输出**单个** NUL 结尾的大块内容（`decode().split("\0")` 只得到 1 个元素），会让整条门禁形同虚设；改用 `-z` + `quotePath=false` 后路径被正确 NUL 分隔且非 ASCII 路径不再被八进制转义。

---

## 2. 净化范围与实测数量

### 2.1 范围定义

- **tracked**（`git ls-files`）：585 个文件
- **untracked 且未被忽略**（`git ls-files --others --exclude-standard`）：407 个文件
- 合计待进入版本控制的候选：**992** 个（文本 923 + 二进制 69）
- `runs/`、`deploy/` 已被 `.gitignore` 排除；`results/` 只有少量已跟踪文件（本次扫描未发现命中）
- **未触碰**任何二进制文件、`runs/`、`results/`、`deploy/`

### 2.2 与任务书估计的差异（因旧模式漏报，估计值系统性偏低）

| 模式 | 任务书估计（未跟踪文件数） | 本次实测（未跟踪文件数 / 出现次数） | 差异原因 |
|---|---|---|---|
| 本地仓库根 `D:` + 反斜杠 + 仓库名 | 6 | **12 / 166** | 旧模式漏报单反斜杠拼写 |
| 本地家目录 `C:` + 反斜杠 + Users | 0 | **1 / 25** | 旧模式漏报；命中在 `_session_import/session.jsonl` |
| 远端用户 `<SERVER_USER>@` | 10 | **11 / 79** | 旧模式该项正常，差异来自扫描口径 |
| 远端根 `${PHB_REMOTE_ROOT}` | 19 | **20 / 537** | 同上 |
| 服务器主机 `<SERVER_HOST>` | 24 | **25 / 102** | 同上 |
| 裸用户名独立 token（源码写作 `"hao" + "yu"`） | 未给 | **29 / 636** | 本次新增的独立 token 规则 |
| **有命中的未跟踪文本文件数** | 36 | **39** | 旧模式漏报 3 个 |
| **有命中的已跟踪文件数** | 3 | **3**（+ 测试模块自身 1，按设计自排除） | 一致 |

### 2.3 替换约定

| 禁止串 | 替换为 |
|---|---|
| 本地仓库根（单/双反斜杠、前斜杠三种拼写） | `<REPO_ROOT>` |
| 本地家目录（单/双反斜杠、前斜杠） | `<LOCAL_HOME>` |
| 远端用户 `@` 形式 | `<SERVER_USER>@` |
| 服务器主机 IP | `<SERVER_HOST>` |
| 远端根路径 | `${PHB_REMOTE_ROOT}` |
| 裸用户名 token | `<SERVER_USER>` |

### 2.4 未套用约定的例外（有意保留）

- **`LICENSE`** 含 `Copyright (c) 2026 Haoyu`。这是**作者署名**，不是用户名或路径组件，套用 `<SERVER_USER>` 会篡改许可证的法律含义。**按约定「仅限用户名/路径组件」不予替换**，且因首字母大写、并非全小写 token，规则不会命中它。**请协调者确认是否接受署名保留**。
- `_session_import/session.jsonl` 中的 636 处裸 token，经逐条核对 **626 处**其实是路径/邮箱的一部分（`${PHB_REMOTE_ROOT}`、`<SERVER_USER>@`），只有 **10 处**是真正的独立用户名引用（如 `用户 <SERVER_USER>`、`<SERVER_USER> 无计算进程`、`server_user=<SERVER_USER>`、`account: <SERVER_USER>`），全部按约定替换。

### 2.5 逐文件替换计数与编码不变式

`FFFD` 列为改动后的 U+FFFD 计数；`dLF`/`dCRLF` 为与改动前的差值（必须为 0）；`BOM` 列为「BOM 状态未变」。

#### A. 已跟踪文件（3 个，本任务直接修正门禁失败）

| 文件 | 仓库根 | 家目录 | 用户@ | 远端根 | 主机 | 裸用户 | FFFD | dLF | dCRLF | BOM |
|---|---|---|---|---|---|---|---|---|---|---|
| `docs/T141_20260917_project_handoff.md` | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260921_phaded_comprehensive_review.md` | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260928_phaded_evidence_run_authorization_packets.md` | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 未变 |

#### B. 未跟踪文件（35 个）

| 文件 | 仓库根 | 家目录 | 用户@ | 远端根 | 主机 | 裸用户 | FFFD | dLF | dCRLF | BOM |
|---|---|---|---|---|---|---|---|---|---|---|
| `.superpowers/nature-skills/skills/nature-figure/INSTALL_MANIFEST.md` | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `.superpowers/sdd/task-0-local-server-thread-fix-brief.md` | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 未变 |
| `.superpowers/sdd/task-0-report.md` | 0 | 0 | 0 | 0 | 2 | 2 | 0 | 0 | 0 | 未变 |
| `.superpowers/sdd/task-2-brief.md` | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `_session_import/session.jsonl` | 142 | 25 | 67 | 509 | 73 | 10 | **186（改动前即为 186，未新增）** | 0 | 0 | 未变 |
| `docs/T141_20260904_archaea_phazh1_candidate_01_status.md` | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260905_ephaz_iphaz_layering_01_status.md` | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260906_archaea_phazh1_domain_annotation_03_status.md` | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260906_conversation_process_and_decision_log.md` | 1 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260906_overall_pha_phb_archaea_audit_03_status.md` | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260906_overall_pha_phb_archaea_audit_04_status.md` | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260909_project_handoff.md` | 1 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260911_phaded_pfam_architecture_01_status.md` | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260912_phaded_reference_structure_comparison_01_status.md` | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260912_project_handoff.md` | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260913_gtdb_phb_degrader_catalog_report_amendment_02.md` | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260913_gtdb_phb_degrader_literature_reconciliation_02_status.md` | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260914_phaded_full_evidence_reconciliation_01_status.md` | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260914_phaded_interpro_full_01_status.md` | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260915_phaded_evidence_completion_01_status.md` | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260915_phaded_reference_panel_acquisition_01_status.md` | 6 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260915_phaded_reference_panel_acquisition_02_status.md` | 6 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260916_phaded_full_library_evidence_amendment_01_status.md` | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260916_phaded_full_library_signalp_01_status.md` | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260916_phaded_reference_panel_acquisition_03_status.md` | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260917_phaded_cys_discovery_hmm_status.md` | 0 | 0 | 0 | 1 | 1 | 1 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260917_phaded_housekeeping_amendments_status.md` | 2 | 0 | 0 | 2 | 1 | 1 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260917_phaded_inphascl_cys_decision_status.md` | 0 | 0 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 未变 |
| `docs/T141_20260920_phaded_structure_verification_status.md` | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `docs/superpowers/plans/2026-09-17-phaded-classification-authority-and-architecture-validation.md` | 2 | 0 | 2 | 2 | 2 | 0 | 0 | 0 | 0 | 未变 |
| `docs/superpowers/plans/2026-09-17-phaded-inphascl-cys-decision-record.md` | 1 | 0 | 2 | 1 | 2 | 0 | 0 | 0 | 0 | 未变 |
| `pipeline/scripts/finalize_phaded_mapping_contract.py` | 0 | 0 | 0 | 2 | 1 | 1 | 0 | 0 | 0 | 未变 |
| `pipeline/scripts/finalize_phaded_pfam_run.py` | 0 | 0 | 0 | 3 | 1 | 1 | 0 | 0 | 0 | 未变 |
| `pipeline/scripts/run_phaded_candidate_hmmsearch.py` | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `pipeline/scripts/run_phaded_pfam_scan.py` | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `pipeline/scripts/run_phaded_priority_interproscan.py` | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 未变 |
| `pipeline/tests/test_amend_phaded_provenance.py` | 0 | 0 | 0 | 0 | 3 | 2 | 0 | 0 | 0 | 未变 |
| `pipeline/tests/test_phaded_housekeeping_amendments.py` | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 未变 |

#### C. 汇总

| 指标 | 值 |
|---|---|
| 扫描候选文件 | 41（另 1 个 `tmp_f14_probe.py` 已不存在，跳过） |
| 实际改动文件 | **38**（3 已跟踪 + 35 未跟踪） |
| 替换总数 | **928** |
| 仓库根 / 家目录 / 用户@ / 远端根 / 主机 / 裸用户 | 165 / 25 / 79 / 537 / 102 / 20 |
| 改动行总数 | 221 |
| 总字节变化 | +1306（占位符比原串长的行多于短的，净增） |
| 写回校验 `roundtrip_equal` | 38/38 全部 True |
| 写回后残留命中 `residual_after_write` | 38/38 全部为空 |
| U+FFFD 变化 | **38/38 全部为 0 变化**（见 §3） |
| 行数变化 | **38/38 全部为 0 变化** |
| BOM / CRLF 变化 | **38/38 全部无变化** |

---

## 3. 编码安全证据

净化程序**全程在二进制模式**下对 ASCII 字节序列做替换：从不 `decode` 再 `encode`，因此 UTF-8 与**改动前就已非法的 UTF-8**都不可能被本程序破坏。每个文件写回前逐项校验，写回后重新读回再校验一次，并由一份独立脚本从磁盘复核。

### 3.1 逐文件不变式（改动后从磁盘重读）

- **U+FFFD（U+FFFD 替换字符）变化为 0**：38/38 个文件满足；其中 **37 个文件改动后 U+FFFD 计数为 0**
- **行数变化为 0**：38/38 个文件满足（按字节数 `LF`、`CRLF`、`CR` 三个计数分别比对，不是只看"行数"）
- **BOM 无变化**：38/38 个文件满足（改动前无 BOM 的仍无 BOM，有 BOM 的仍有 BOM）
- **总字节变化 = 每条规则 (替换长度 − 命中长度) 之和**：38/38 精确相等
- **逐行精确复现**：对每个改动文件的每一行，把规则重新作用在**改动前的那一行**上，结果必须与**改动后的那一行**逐字节相同（任何一条规则都不跨行匹配，故这是全文件级别的证明）。38/38 通过，无一行例外 —— 这条检查强制「除命中区间外一个字节都不许动」，排除了重排、重换行、重格式化
- **写回往返一致**：38/38 `roundtrip_equal == True`

### 3.2 唯一需要说明的一处：`_session_import/session.jsonl`

该文件**改动前就含有 186 个 U+FFFD 字节**（导入的会话记录本身在导入阶段已损坏），另有 6 行含孤立 `CR`。本次净化：

- **未新增任何一个 U+FFFD**（改动前 186 → 改动后 186，差值 0）
- 行数、`LF`/`CRLF`/`CR` 计数全部不变（含 6 个孤立 `CR` 原样保留）
- 逐行精确复现通过

即：该文件的历史损坏**不是本次操作造成的**，本次操作也没有加重它。若协调者希望仓库中不出现非合法 UTF-8 文本，应单独决定是否将该 5.3 MB 会话记录排除出提交范围 —— **本任务不擅自删除任何文件**。

### 3.3 全范围编码概览（未经本任务改动的文件，仅报告）

- 范围内文本文件中含 U+FFFD 的：`.superpowers/nature-skills/skills/nature-figure/scripts/audit_figure_collisions.py`（2 处，本任务未改动，属既有损坏）
- 范围内含 UTF-8 BOM 的文本文件共 51 个（多为 `research/` 下的 JSON 与若干 `.diff`），均未被本任务改动
- 范围内含 CRLF 的文本文件若干（如 `pipeline/scripts/audit_discovery_layer_usage.py` 2116 处、`research/` 下 JSON），均未被本任务改动

---

## 4. 净化后仍含身份信息的文件（完整清单）

对范围 995 个文件（文本 926 + 二进制 69）重新全量扫描，仍命中身份模式的文件只有 **1 个**：

| 文件 | 命中的模式 | 说明 |
|---|---|---|
| `pipeline/tests/test_public_repo_safety.py` | 仓库根 10、家目录 10、用户@ 2、远端根 6、主机 3、裸用户 9 | **按设计自排除**，必须保留这些模式本身 |

- 二进制文件中含身份模式的：**0 个**（69 个二进制文件全部扫描，字节未改动）
- 本任务完成后新出现的未跟踪文件中含身份模式的：**0 个**

**结论**：除安全测试模块自身外，版本控制范围内不存在任何身份串。

---

## 5. 凭据扫描结果

对 994 个文件（tracked 585 + untracked-nonignored 409）扫描以下签名：

PEM 私钥起始标记（`BEGIN` / `PRIVATE KEY` 由连字符包裹的五短横线形式）、`OPENSSH PRIVATE KEY`、`RSA PRIVATE KEY`、`EC PRIVATE KEY`、`PGP PRIVATE KEY BLOCK`、`ghp_/gho_/ghu_/ghs_/ghr_`、`github_pat_`、`AKIA…`、`xox[baprs]-…`、`sk-…`、`sk-ant-…`、`AIza…`、`password/passwd/passphrase = …`、OpenSSH 私钥 base64 魔数、base64 长行启发式。

### 5.1 结果：**没有发现任何真实凭据**

| 命中 | 性质 | 判定 |
|---|---|---|
| `pipeline/tests/test_public_repo_safety.py`：三种 PEM 私钥起始标记（RSA / OPENSSH）与其模式定义本身 | 测试自身的**合成样本与模式定义**，不是密钥 | 非凭据；该文件按设计自排除 |
| `tmp_positive_search/{BAA32541.faa,BAD70022.faa,refs.faa}`：base64 长行启发式 | **蛋白质序列**恰好满足 base64 字符集 | 假阳性，非凭据 |

### 5.2 文件名与磁盘检查

- 范围内仅 `.env.example`（示例文件，无真实值）命中文件名关键词；`.env` / `*.pem` / `*.key` / `id_rsa*` / `id_ed25519*` / `.ssh/` 在范围与磁盘上**均不存在**（已被 `.gitignore` 覆盖）

**按任务要求：未发现凭据，因此无需 STOP；也未曾对任何疑似凭据做静默处理。**

---

## 6. 残余风险声明（审计看不到的部分）

本审计的结论**仅限于**它能看到的范围。以下内容**未被验证**，不得被读作「已确认无泄露」：

1. **二进制文件**：69 个二进制文件因含 NUL 字节/二进制后缀被按字节跳过（字节保持原样）。本次对它们做了身份模式与凭据签名的原始字节扫描，**0 命中**，但二进制内部可能以**压缩、编码或嵌入对象**形式承载身份串，字节级正则**看不见**这些形式。具体地，本次未解压、未解码：
   - PDF / PNG / SVG（`9.21汇报/`、`results/figures/` 等）
   - `.docx` / `.xlsx`（如 `docs/PHB_核心文档清单.docx`）
   - `.gz`（含 `data/hmms` 之外的各种归档）
   - `.hmm`、`.faa`、`.fna` 等生物序列文件（按文本处理，但序列内容未被逐一人工核读）
2. **被 gitignore 的树**：`runs/`（实测 1,016 GB）、`deploy/`、`results/` 的绝大部分、`data/*`、`*.log` 等**未被扫描**。它们当前不进版本控制，但：
   - 若日后用 `git add -f` 强制加入，本次结论**不适用**
   - 若日后放宽 `.gitignore`，必须重新执行本审计
   - 服务器侧 `runs/` / `deploy/` 的同名身份串**仍然存在**（本任务只处理本地工作树）
3. **拼写不同的身份串**：本审计只枚举了 6 类拼写（仓库根的三种斜杠形式、家目录的三种形式、`<SERVER_USER>@`、`${PHB_REMOTE_ROOT}`、`IP`、裸 "hao"+"yu" token）。以下**不在**覆盖范围：
   - 其他盘符或路径（`E:`、`F:`、网络共享 `\\host\share`、`/mnt/*`、`/media/*`）
   - 用户名的大小写变体 `Haoyu` / `HAOYU`（**已确认仅出现在 `LICENSE` 的署名处**，按约定保留）
   - 其他缩写或别名（`hy`、`haoy`、`h.y.`）、邮箱域名、机构名、其他主机名/内网 IP 段
   - Base64 / URL 编码 / JSON 转义后的身份串，或跨行折断的串
   - 环境变量名以外的**值**（如 `PHB_REMOTE_ROOT=<实际路径>`）
4. **历史提交**：本审计只针对**当前工作树**。**既有 git 历史**（已提交的版本）可能仍含身份串；本任务**未改写历史**，`git log` 的旧提交内容未被检查。若历史中存在身份串，仅靠本次净化不能消除 —— 需要单独决策（如 `filter-repo` 或接受既成事实），且**必须另行授权**。
5. **测试模块自身的模式串**：`pipeline/tests/test_public_repo_safety.py` 仍以明文保存被禁止的拼写（这正是它的功能）。它在 track list 中被自排除，因此**门禁对它自己的内容不做保证**；这是**有意的设计取舍**，不是遗漏。
6. **并发写入**：本任务运行期间工作树在变化 —— `tmp_f14_probe.py` 在扫描中途消失，另有 `.tmp_debug.py`、`.tmp_debug2.py`、`_audit_pool_in_evidence.py`、`--amino`、`-o`、`=.8=%d` 等无意义命名的未跟踪文件出现（经扫描**均不含身份模式**）。因此：**提交前必须重跑一次全范围扫描**，以覆盖本任务结束后新产生的文件。本报告的数字对应的是本次扫描时点。
7. **审计工具的正确性**：本次净化在编写过程中，**审计脚本自身先后出现过 3 处缺陷**（字节增量算术未考虑规则先后顺序、diff 分段把匹配上下文并入改动块、`str.splitlines` 把孤立 `CR` 当换行），三处均已修正并被最终校验暴露。这说明「净化脚本报告 OK」本身不足以作为结论 —— 因此本报告的主要证据是**逐行精确复现**与**独立脚本从磁盘复核**，而不是净化脚本的自述。
8. **本报告本身也必须净化**：第一版报告为了写清「修复前/修复后」的**精确模式字符串**，直接引用了被禁止的字面量，结果**自己的全量扫描把这份报告判为违规文件**（首轮复扫：`REPO_ROOT` 2、`SERVER_USER_AT` 5、`REMOTE_ROOT` 6、`SERVER_HOST` 2、裸用户名 17）。随后按与安全测试模块相同的做法改写为**拆字面量**记法（如 `"hao" + "yu"`）与符号记法（`<BS>`、`<FS>`），复扫通过。留此记录的原因是：**「写文档说明哪些串被禁止」这一行为本身就会引入被禁止的串** —— 后续任何新增说明文档都必须把安全测试作为门禁跑一遍，不能假定文档天然安全。

---

## 7. 验证命令与输出

工作目录：本仓库根。以下三条命令均已实跑。

### 7.1 `python -m unittest pipeline.tests.test_public_repo_safety -v`

```text
test_all_forbidden_samples_are_flagged (...) ... ok
test_doubled_and_forward_slash_spellings_are_flagged (...) ... ok
test_forward_slash_remote_root_placeholder_is_not_flagged (...) ... ok
test_placeholders_and_benign_text_are_not_flagged (...) ... ok
test_predicate_ignores_already_sanitized_placeholder_forms (...) ... ok
test_single_backslash_windows_paths_are_flagged (...) ... ok
test_tracked_files_do_not_expose_machine_identity_or_credentials (...) ... ok

----------------------------------------------------------------------
Ran 7 tests in 2.512s

OK
```

`EXITCODE=0`

### 7.2 `python -m compileall -q pipeline/scripts pipeline/tests`

```text
（无输出即成功）
```

`EXITCODE=0`

### 7.3 `git diff --check`

```text
（无输出即无空白错误）
```

`EXITCODE=0`

### 7.4 未执行的动作

- 未运行 `git add` / `git commit` / `git push` / `git reset` / `git checkout`
- 未写索引：`git status --short` 中不存在任何已暂存条目
- 提交由协调者负责

---

## 8. 交给协调者的动作项

1. **重新扫描**：本任务结束后如有新文件产生（尤其 `.tmp_debug*.py`、`_audit_pool_in_evidence.py` 之类临时文件），提交前重跑 `python -m unittest pipeline.tests.test_public_repo_safety -v` 与一次全范围身份扫描
2. **决定 `_session_import/session.jsonl` 是否入库**：5.3 MB 会话记录，含既有 186 个 U+FFFD、6 个孤立 `CR`、且是本次替换量最大的单个文件（816 处）。已净化，但其入库价值与体积应由协调者判断
3. **确认 `LICENSE` 署名保留**：`Copyright (c) 2026 Haoyu` 属作者署名，未按约定替换
4. **确认临时/无意义文件是否入库**：`--amino`、`-o`、`=.8=%d`、`.tmp_debug.py`、`.tmp_debug2.py`、`_audit_pool_in_evidence.py`、`.mplconfig/`、`_session_import/`、`9.21汇报/` 等属于工作残留，本任务**未删除**，请决定是否加入 `.gitignore` 或排除出提交
5. **历史提交未处理**：本次只净化当前工作树，未改写 git 历史
