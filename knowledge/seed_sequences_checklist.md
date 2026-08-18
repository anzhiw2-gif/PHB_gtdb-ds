# PHB 降解基因 — 种子序列收集清单（v0.3，供 HMM 构建）

> v0.3：加入穷尽式 PubMed 调研（728 篇）确认的已表征酶清单与关键催化位点。
> 完整报告：research/pubmed/comprehensive_summary.md；分类体系：
> knowledge/family_classification.md（PhaDED 八大超家族）。

## F. 古菌 PHB 降解相关（v0.2 新增，文献核实）

### F1. patatin 样解聚酶（ArchPhaZ_patatin）— 独立 HMM（细菌 HMM 检不出）
| 来源菌 | 酶 | 文献证据 | UniProt |
|--------|----|---------|---------|
| *Haloferax mediterranei* | **PhaZh1**（patatin 样胞内解聚酶，321aa，Ser47+Asp195，与 bdhA 成簇） | PMID 25710370 | **I3RBH0** |
| *Natronomonas moolapensis* | PHB depolymerase（reviewed，323aa） | UniProt 注释 | **M1XPT2** |
| 筛选 profile | Pfam-A "Patatin"（PF01734 家族）HMM | hmmfetch 提取 | data/hmms/patatin_Patatin.hmm |

### F2. 经典 PHB 解聚酶家族酯酶（α/β 水解酶型，细菌 HMM 可覆盖）
| 来源菌 | 酶 | UniProt |
|--------|----|---------|
| *Haladaptatus paucihalophilus* | Esterase, PHB depolymerase family | A0A1M6SXK9, E7QQJ0 |
| *Halomarina oriensis* | PHB depolymerase family esterase | A0A6B0GML5 |
| *Haloarcula salinisoli* | PHB depolymerase family esterase | A0A8J7YBF6 |
| 等共 13 条 | | data/seeds/archaea_seeds.faa |

### F3. 古菌动员通路（辅助，不建筛选 HMM——ECH 太广谱）
| 基因 | 功能 | 文献 |
|------|------|------|
| PhaJ1-5 | 烯酰-CoA 水合酶（PhaJ1 颗粒相关，PHA 动员） | PMID 27052994（119 条 Haloferax ECH 已收集） |
| BdhA（古菌） | 3HB 脱氢酶（6 条已收集） | — |
| PhaP（古菌 phasin） | 颗粒蛋白（I3R9Z2 已在原种子） | PMID 22247127 |

## G. 穷尽式调研补充的已表征酶（v0.3，PubMed 728 篇核实）

### G1. 关键催化位点（HMM 验证锚点）
| 酶 | 催化位点 | PMID |
|----|---------|------|
| *R. eutropha* PhaZ1（胞内） | **Cys183-Asp355-His388**（Cys 型！） | 16233560 |
| *R. eutropha* PhaZd | Ser190-Asp266-His330 | 16199568 |
| *R. rubrum* PhaZ1（周质） | Ser42-Asp138-His178（type II CD） | 15489436 |
| *B. thuringiensis*（胞内） | G-**W**-S102-M-G；S102A 失活 | 16936025, 39592048 |
| *P. putida* PhaZKT（mcl 胞内） | Ser102-Asp221-His248 + lid | 17170116, 41055782 |
| *A. faecalis* T1（胞外） | Ser139-Asp214-His273 | (见 comp 报告) |
| *T. thermophilus* HB8（胞外） | **Ser183-Glu310-His405**（**Glu 型**！酸性残基用 Glu 而非 Asp） | 19214501 |
| *P. lemoignei* PhaZ7（nPHB） | His47-Ser136-Asp242 + 氧阴离子 His306 | 12855176 |
| *Wautersia/R. eutropha* 3HB 寡聚体水解酶 | S-V-S*-N-G（新基序） | 16233278 |
| *Acidovorax* sp. SA1 i3HBOH | 292 aa，G-X-S-X-G 型 | 12070691 |
| 3HB 脱氢酶（BdhA） | *S. meliloti* 258 aa 短链脱氢酶超家族 | 9922248 |

> 注：催化三联体的"酸性残基"Asp(D)与 Glu(E)可互换（如 *T. thermophilus* HB8 用 Glu310）。
> 验证脚本 `08_validate.py` 已按 Ser/Cys-Asp/Glu-His 放宽，勿再写死为 Asp。

### G2. 补充已表征酶（种子扩充参考）
- **mcl 胞外**：*P. alcaligenes* LB19/M4-7、*Pseudomonas* sp. RY-1、*Streptomyces exfoliatus* K10 PhaZSex2、*Bdellovibrio* PhaZBd（Bd3709）、*T. thermophilus* HB8
- **胞内**：*A. vinelandii* PhbZ1（含硫降解）、*P. denitrificans* i-PhaZ+PhaZc+Hbd、*P. oleovorans*、*B. megaterium* PhaZ1（nPHB 特异）
- **寡聚体水解酶**：*A. faecalis*（68-74 kDa）、*Pseudomonas* sp. A1（72.9 kDa，无 GXSXG）、*Zoogloea ramigera* 二聚体水解酶、*P. lemoignei* PhaZc 同源
- **BdhA**：*R. pickettii* T1 BDH1/2/3、*P. denitrificans* Hbd（四聚体 29kDa×4）
- **真菌**：*P. funiculosum*（结构 2D80）、*P. citrinum/expansum/pinophilum*
- **P. lemoignei 全谱**：phaZ1=C、2=B、3=D、4=PHB、5=A、6=真 PHV 酶、7=nPHB 特异（PMID 7836292, 10742216）

## A. 胞外 PHA 解聚酶（e-PhaZ）— 独立建模

| 来源菌 | 酶 | 文献证据 | UniProt（待填） |
|--------|----|---------|----------------|
| *Comamonas acidovorans* YM1609 | 胞外 PHB 解聚酶 | PMID 9406404 | |
| *Pseudomonas lemoignei* | PhaZ1–7 系列 | Jendrossek 综述 PMID 12213937 | |
| *Alcaligenes faecalis* T1 | 胞外 P(3HB) 解聚酶 | PMID 3942778 | |
| *Ralstonia pickettii* T1 | 胞外 P(3HB) 解聚酶 | PMID 18340545/24146107 | |
| *Thermus thermophilus* HB8 | TTHA0199 胞外 | PMID 19214501 | |
| *Streptomyces ascomycinicus* | 胞外 PHB 解聚酶 | PMID 23951224 | |
| *Streptomyces exfoliatus* K10 | mcl-PHA 解聚酶 | PMID 26156240 | |
| *Undibacterium* sp. KW1/YM2 | PhaZUD | PMID 32369496 | |
| *Alteromonas* sp. | 胞外 P(3HB) 解聚酶 | PMID 40500476 | |
| *Penicillium funiculosum* | 真菌 PHB 解聚酶（晶体结构） | PMID 16405909 | |
| *Lihuaxuella thermophila* | 嗜热 PHB 解聚酶 | PMID 36222314 | |
| *Nocardiopsis dassonvillei* | PHB 解聚酶 | PMID 41151231 | |

## B. 胞内 PHA 解聚酶（i-PhaZ）— 独立建模

| 来源菌 | 酶 | 文献证据 | UniProt（待填） |
|--------|----|---------|----------------|
| *Cupriavidus necator* H16 | PhaZ1/Za1（主要胞内） | PMID 29678915 | |
| *Cupriavidus necator* H16 | PhaZ2/Z3 同工酶 | Jendrossek 综述 | |
| *Rhodospirillum rubrum* | PhaZ2（胞质主酶）/PhaZ3 | PMID 41664065（全文） | |
| *Rhodospirillum rubrum* | PhaZ1（周质型，特殊） | PMID 15489436 | |
| *Bacillus thuringiensis* | 新型胞内 phaZ | PMID 16936025 | |
| *Azospirillum brasilense* | PhaZ | PMID 12898135 | |
| *Sinorhizobium meliloti* | 胞内 PhaZ | PMID 20346169 | |
| *Pseudomonas putida* KT2440/2442 | mcl-PHA 胞内解聚酶 | PMID 17170116/19788655；lid 研究 PMID 41055782 | |

## C. 寡聚体水解酶（OH）

| 来源菌 | 酶 | 文献证据 | UniProt（待填） |
|--------|----|---------|----------------|
| *Cupriavidus/Wautersia eutropha* H16 | 胞内 3HB-寡聚体水解酶 | PMID 16030206/16233278/12070691 | |
| *Alcaligenes faecalis* | 胞外 D(-)-3HB 寡聚体水解酶 | PMID 6626560 | |
| *Paracoccus denitrificans* | D(-)-3HB 寡聚体水解酶 | PMID 11814660 | |
| *Zoogloea ramigera* I-16-M | 3HB-二聚体水解酶 | PMID 7285912/1476778 | |
| *Pseudomonas* sp. | 胞外 3HB 寡聚体水解酶 | PMID 8981982 | |

## D. 3HB 单体代谢（可选，辅助验证）

| 基因 | 功能 | 文献证据 |
|------|------|---------|
| bdhA | 3-羟基丁酸脱氢酶 | Azospirillum 综述 |
| 乙酰乙酸-CoA 合成酶/转移酶 | 乙酰乙酸→乙酰乙酰-CoA | Azospirillum 综述 |
| β-酮硫解酶 | 乙酰乙酰-CoA→2 乙酰-CoA | Azospirillum 综述（与合成共用，需方向区分） |

## E. 颗粒/调控/标记基因（可选，基因簇背景）

| 基因 | 功能 | 结构域 |
|------|------|--------|
| phaP（phasin） | 颗粒涂层蛋白 | phasin_2 = PF09361（已核实） |
| phaR | 转录调控子 | — |
| apdA | 解聚激活蛋白 | — |
| phaC | PHA 合酶（合成方向，簇共定位用） | — |

## 备注

- PhaZ 命名跨菌混乱：以已表征酶的 UniProt 序列为准，不依赖基因名。
- 分类参考：PhaDED（PHA Depolymerase Engineering Database, 2009）。
- 建议每家族 HMM 构建后做催化位点（GXSXG + Ser-Asp-His/Glu）与信号肽
  （SignalP）双重验证。
- ⚠️ **种子库含 18 条"同源物"而非 PHB 降解基因**（13 条真核 BDH1/BDH2 酮体代谢酶 +
  5 条尼龙水解酶 nylB/nylC），保留用于 HMM 序列多样性，但**非 PHB 降解基因**。
  详见 `pipeline/seeds/seeds_annotation.md`。
