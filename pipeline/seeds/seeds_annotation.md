# 种子库同源物标注(非 PHB 降解基因)

> 说明:78 条 curated 种子中,绝大多数是真正的 PHB 降解相关基因;但有 18 条是
> **同源折叠/同源代谢酶**,并非 PHB 降解基因。它们对构建稳健的 HMM 有序列多样性
> 价值,故**保留**,但此处明确标注,避免被误读为"PHB 降解基因种子"。
>
> 依据:UniProt 逐一核实 + 催化反应/生物学背景判断。

## 一、3HB 脱氢酶同源物(13 条,真核酮体代谢,非 PHB 降解)

这 13 条是哺乳动物/脊椎动物的 **3-羟基丁酸脱氢酶(BDH1/BDH2)**,参与**酮体代谢**
(饥饿时把血液 3-羟基丁酸氧化供能),**不参与 PHB 降解**。它们与细菌 bdhA 同源
(共享 EC 1.1.1.30 + 短链脱氢酶折叠 + 线粒体 α-变形菌内共生起源),但生物学
上下文完全不同。

| accession | 物种 | 基因 | 蛋白 | 说明 |
|-----------|------|------|------|------|
| P29147 | Rattus norvegicus | Bdh1 | D-β-hydroxybutyrate dehydrogenase, mitochondrial | 真核酮体代谢 |
| Q02337 | Bos taurus | BDH1 | 同上 | 真核酮体代谢 |
| Q02338 | Homo sapiens | BDH1 | 同上 | 真核酮体代谢 |
| Q80XN0 | Mus musculus | Bdh1 | 同上 | 真核酮体代谢 |
| P86198 | Mesocricetus auratus | BDH1 | 同上 | 真核酮体代谢 |
| Q5ZJZ5 | Gallus gallus | BDH1 | 同上 | 真核酮体代谢 |
| D4A1J4 | Rattus norvegicus | Bdh2 | SDR family member 6(**非 EC 1.1.1.30**) | 广谱短链脱氢酶 |
| Q561X9 | Danio rerio | bdh2 | 同上 | 同上 |
| Q8JZV9 | Mus musculus | Bdh2 | 同上 | 同上 |
| Q9BUT1 | Homo sapiens | BDH2 | 同上 | 同上 |
| C1C4R8 | Aquarana catesbeiana | bdh2 | 同上 | 同上 |
| Q3KPT7 | Xenopus laevis | bdh2 | 同上 | 同上 |
| Q3T046 | Bos taurus | BDH2 | 同上 | 同上 |

> ⚠️ 特别注意:其中 **BDH2** 实为广谱短链脱氢酶(SDR family member 6),
> **连 EC 1.1.1.30 都不是**,是被 UniProt 历史命名"3-hydroxybutyrate dehydrogenase
> type 2"误导进来的。真核 BDH1 才是真正的 EC 1.1.1.30(酮体 3HB 脱氢酶)。

**真正的细菌 PHB 降解 bdhA 只有 2 条:**
- O86034(Sinorhizobium meliloti bdhA)
- Q9X6U2(Cupriavidus necator hbdH1)

## 二、尼龙寡聚体水解酶同源物(5 条,非 PHB 寡聚体水解酶)

这 5 条是 **nylB/nylC(6-氨基己酸寡聚体水解酶)**,降解**尼龙-6 寡聚体**,与
3HB-寡聚体水解酶同属 **Ntn(N-末端亲核)水解酶折叠**,但底物/功能不同。

| accession | 物种 | 基因 | 蛋白 |
|-----------|------|------|------|
| Q79F77 | Paenarthrobacter ureafaciens | nylC | 6-aminohexanoate-oligomer endohydrolase |
| Q1EPR4 | Kocuria sp. KY2 | nylC | 同上 |
| Q1EPR5 | Agromyces sp. KY5R | nylC | 同上 |
| P07061 | Paenarthrobacter ureafaciens | nylB | 6-aminohexanoate-dimer hydrolase |
| P07062 | Paenarthrobacter ureafaciens | nylB' | 同上 |

**真正的 3HB 寡聚体水解酶(Q4W8C9 = phaZc, C. necator H16_A1335)** 已混在同一
`oligomer_hydrolase` 标签下,是正确条目。

## 三、真核 PHB 解聚酶(1 条,正确,非"同源物"混淆)

| accession | 物种 | 说明 |
|-----------|------|------|
| B2NHN2 | Talaromyces funiculosus(真菌) | **真 PHB 解聚酶**,真菌分泌胞外解聚酶降解 PHB,正确 |

> 真菌是唯一含真正 PHB 降解基因的真核谱系;其解聚酶与细菌解聚酶同属 α/β 水解酶折叠。

## 四、结论

- **18 条同源物**(13 真核 BDH + 5 尼龙水解酶)**保留**,用于 HMM 的序列多样性;
- 但在方法/论文中,应表述为"**BdhA/OH 的 HMM 纳入了折叠同源物(真核 3HB 脱氢酶、
  尼龙寡聚体水解酶)以增强稳健性**",而非"全部种子均为 PHB 降解基因";
- 真核 BDH 只影响 HMM 构建(筛选对象是 GTDB 细菌+古菌,无真核),影响有限;
  尼龙水解酶则可能被 OH HMM 检出,使 OH 计数轻微虚高(尼龙降解菌在自然界稀有,
  影响小)。
