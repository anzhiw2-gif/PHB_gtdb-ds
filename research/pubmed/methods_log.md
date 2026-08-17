# PubMed 检索方法学日志（2026-06-01）

工具：NCBI E-utilities（pubmed_api.py 包装脚本，rate-limit 合规）
检索策略：8 组查询，结果取 top-N 后合并去重（union = 155 PMIDs）

| # | 查询词 | PubMed 命中 | PMC 命中 | 说明 |
|---|--------|------------|----------|------|
| 1 | polyhydroxybutyrate depolymerase | 191 | 810 | 核心酶学 |
| 2 | polyhydroxyalkanoate depolymerase | 83 | 145 | PHA 广义 |
| 3 | PHB depolymerase | 194 | 447 | 别名 |
| 4 | phaZ AND depolymerase | 87 | 389 | 基因名 |
| 5 | poly(3-hydroxybutyrate) AND (degradation OR depolymerization) | 2415 | 9783 | 宽泛降解 |
| 6 | polyhydroxyalkanoate AND (genome OR metagenome) AND (degradation OR depolymerase) | 855 | 4873 | 基因组/宏基因组筛选 |
| 7 | PHA degradation AND bacteria | 3304 | 17793 | 细菌降解 |
| 8 | polyhydroxyalkanoate degradation review | 535 | 7396 | 综述 |

要点：
- 核心酶学文献规模 ~200 篇；基因组/宏基因组层面筛选研究 ~855 篇（含共现词），
  说明"全库规模解聚酶挖掘"仍有空间但已有大量工作。
- 合并去重后 155 篇进入摘要筛选阶段（top-N 抽样，非全量）。
- 数据文件：research/pubmed/search_*.json（检索结果）、count_*.json（各库命中数）、
  union_pmids.json（合并 PMID）。
