# CHANGELOG — PHB_gtdb-ds 项目变更与事件记录

> 本文档记录本次会话(2026-08-17)对项目的全部改动,以及一次重要的运维事件及其根因。
> 仓库: https://github.com/anzhiw2-gif/PHB_gtdb-ds

---

## 一、本次会话改动总览

按时间顺序,共 11 个 git 提交:

| commit | 类别 | 内容 |
|--------|------|------|
| `8162e38` | 工程 | 初始提交:重写 README、统一 9 家族口径、重写 11_clusters/10_distribution、归档探索脚本、.gitignore/.gitattributes |
| `bab3d0f` | 文档 | 修正 analysis_plan_draft 过时的种子数/家族口径 |
| `55a8ee9` | 结果 | 生态分布表 + 19 个 HMM profiles + tier1 系统发育树入库 |
| `86d8ff3` | 文档 | 更新 README/结果报告状态 |
| `3ccea90` | 文献 | 核心文献清单(Word)+ 171 篇核实表 + 生态/基因簇脚本产出 |
| `a7570df` | 文献 | 核心文献清单更新(16 篇付费全获取 + Thermus Glu 例外) |
| `512000c` | 文献 | 清单更新版(16 篇付费全获取 + 补读要点) |
| `ae47856` | 文献 | 合并清单到正式文件名,删 v2 临时版 |
| `148630a` | 修复 | 催化三联体验证放宽为 Asp/Glu(Ser/Cys 型) |
| `bd24294` | 结果 | Glu 修正重跑 — tier1 ePhaZ 38,692 / iPhaZ 32,926(共 +492 序列) |
| `30eb5b8` | 文档 | 更新结果报告 tier1 数字 |

### 本次会话完成的 4 大类工作

1. **工程修复与 git 化**:重写 `11_clusters.py`(补 ±flank_kb 基因簇共定位)、增强 `10_distribution.py`(补生态关联 + 修 accession 前缀 bug)、归档探索脚本、新增 `.gitignore`/`.gitattributes`/`environment.yml`/`docs/reproducibility.md`/`pipeline/sync_from_server.sh`。
2. **服务器执行**:连 T141 跑通 10_distribution(生态分布)与 11_clusters(patatin 基因簇,620 古菌基因组),同步结果。
3. **文献核实 + 全文补读**:171 篇 PMID 经 Europe PMC + OpenAlex 双重核实(全部真实、0 撤稿);生成 `docs/PHB_核心文献清单.docx`;经南科大机构访问读齐 16 篇付费核心文献。
4. **Glu 三联体修正**:Thermus HB8(S183-E310-H405)揭示催化三联体酸性残基 Asp/Glu 可互换;修正 `08_validate.py` + 两处知识文档,服务器重跑整条验证链。

---

## 二、事件记录:tier1.faa 文件被误写为 "0\n"

### 现象

- 时间:2026-08-17 约 21:16(服务器 T141)
- 文件:`data/screen/tiers/{ePhaZ,iPhaZ,OH,ArchPhaZ_hydrolase,ArchPhaZ_patatin}_tier1.faa`(共 5 个)
- 变化:文件内容变为 `0\n`(2 字节),原始序列全部丢失
- 未受影响:`*_tier1.ids`、`*_tier1.tbl`(底层 hmmsearch 结果)、`*_validated.faa` 均完好

### 根因(已复现验证)

不是外部入侵、不是管线 bug、也不是其它项目脚本——**是我(助手)自己一条有 bug 的诊断命令**:

1. 记录"改动前 baseline"时,我用了一条**内联命令**:
   ```bash
   for f in data/screen/tiers/*_tier1.faa; do ... grep -c "^>" $f; done
   ```
2. 该命令经 **PowerShell → ssh** 传递时,`"^>"` 中的**双引号被剥离**,`>` 变成了 shell 重定向符;
3. `grep -c "^>" $f` 变成 `grep -c ^ > $f` =「统计 `^` 的匹配数,并把结果**重定向写入 $f**」;
4. `grep -c ^`(读空 stdin)输出 `0`,`> $f` 就把每个 tier1.faa **截断并写入 `0\n`**。

复现验证:
```
grep -c "^" t.txt    → 2      (引号在,正常)
grep -c ^ > t.txt    → "0\n"  (引号被吃,文件被覆盖)
```

`*_validated.faa` 未遭殃的原因:命令用了 `&&` 串联,第一个循环里 `grep` 退出码为 1,`&&` 链断裂,第二个循环(validated)未执行。

### 影响与恢复

- **数据无实质丢失**:`tier1.ids`/`tier1.tbl`(hmmsearch 原始结果)完好,tier1.faa 只是从它们派生的序列文件,可随时重建。
- **已恢复**:本次 Glu 修正重跑(08c)已从 `validated.faa` 完整重建 tier1.faa,内容正确。
- **不影响 Glu 对比结论**:"改动前"基准用的是最终报告数字(38,275 等),非被误删的文件。

---

## 三、运维教训(务必遵守)

1. **凡涉及 `>`、`"`、`|`、`$()`、中文括号的命令,一律写成 `.sh` 脚本文件 `scp` 到服务器执行,禁止内联拼字符串**。本次会话因此类引号问题多次出错,最终导致一次误删。
2. **对派生的结果文件(`tier1.faa` 等)做任何批量操作前,先确认底层源文件(`.ids`/`.tbl`/`validated.faa`)完好**——它们是真正的不可再生中间结果,派生文件可重建。
3. **重跑链是幂等的**:`08_validate.py` → `08c_tier_rescore.sh` → `09a_tier1_summary.py` → `09d_patatin_filter.py` → `10_distribution.py` 从输入文件出发可完整重建所有 tier 结果,是数据损坏后的兜底手段。
4. **服务器重跑链记录**:见 `pipeline/scripts/` 各脚本 docstring;数据溯源见 `docs/reproducibility.md`。

---

## 四、patatin 基因簇复筛结论更正(重要)

### 之前的口误

对古菌 patatin 做 ±10kb 基因簇复筛时,曾表述为"1,372 条里仅 **63 条(4.6%)** 邻近 PhaC/phasin = 真 PhaZh1 型"。**该判定标准错误。**

### 更正后的正确结论

- 判定真 PhaZh1 型应看"是否邻近 **PHB 代谢基因簇**",而 **PhaZh1 与 bdhA(降解支路)成簇**(Liu 2015, PMID 25710370),不是与 phaC(合成簇)成簇;
- 正确统计(邻近任一 PhaC/phasin/BdhA/PhaJ):

| 邻近标记 | 位点数 |
|---------|-------|
| BdhA(降解支路) | 274 |
| PhaC(合成簇) | 63 |
| PhaE(合成亚基,PF09712) | 56 |
| PhaJ(动员) | 11 |
| **任一 PHB 代谢基因(并集)** | **340(24.8%)** |
| 无 PHB 上下文 | 1,032(75.2%) |

- 因此正确表述为:1,372 条 patatin 中 **24.8% 处于 PHB 代谢基因上下文**(以 BdhA 降解支路为主),其余 75.2% 为广谱磷脂酶/酯酶背景;
- 注:补充古菌 PHA 合成酶亚基 PhaE(PF09712)与颗粒区蛋白 PHA_gran_rgn(PF09650)两个 marker HMM 后重跑,并集计数由 339 微调至 340,结论不变;
- 更深一层:patatin 是广谱脂质水解酶结构域,且 PhaZh1 体内角色有限、**PhaJ 才是古菌 PHB 动员主通路**——已写入 `final_results_report.md` §2.2 的"生物学 caveat"。
