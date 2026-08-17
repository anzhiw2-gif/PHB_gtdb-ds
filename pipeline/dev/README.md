# pipeline/dev — 探索性/一次性脚本归档

> 这些脚本是调研与种子收集阶段的探索性产物，不属于 01–11 主流程。
> 保留供溯源，不再作为主流程调用。

| 脚本 | 用途 | 说明 |
|------|------|------|
| `probe_dbs.py` / `probe_pazy.py` / `probe_pazy2.py` | 探测 ESTHER / PAZy 数据库 API | 种子收集前的接口验证 |
| `query_archaea2.py` / `query_archaea3.py` / `query_archaea4.py` / `query_archaea_uniprot.py` | 古菌种子查询（UniProt） | 古菌 PhaZ/PhaJ 种子收集的迭代脚本 |
| `check_seeds.py` | 种子 FASTA 质量检查 | 通用 QC 工具，可按需调用 |
| `resume_0611.sh` | 2026-06-11 断点续跑脚本 | 一次性，仅作当时运行记录 |
