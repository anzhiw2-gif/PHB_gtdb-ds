#!/usr/bin/env python3
"""
02_collect_seeds.py — 收集 PHB/PHA 降解相关种子序列（UniProt）

按 EC 号 + 关键词 + 关键物种基因名多路查询 UniProt REST API，
输出：
  seeds/phb_degradation_seeds.faa        （合并 FASTA，含冗余注释）
  seeds/seeds_manifest.tsv               （accession / organism / gene / EC / reviewed / 分类）
  seeds/seeds_curated.faa                （仅 reviewed + 去冗余的精选种子）

用法（HPC/本地均可）：
  python 02_collect_seeds.py [--outdir seeds] [--min-reviewed-only]
依赖：requests, pandas（可 pip install 或 uv run --with requests,pandas）
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import sys
import time

import requests

BASE = "https://rest.uniprot.org/uniprotkb"
STREAM = f"{BASE}/stream"
SEARCH = f"{BASE}/search"

# (查询名, 查询串, 分类标签)
QUERIES = [
    # --- 胞外 PHA/PHB 解聚酶（EC 3.1.1.75/76）---
    ("ePhaZ_ec3175", "ec:3.1.1.75 AND reviewed:true", "e-PhaZ_EC3.1.1.75"),
    ("ePhaZ_ec3176", "ec:3.1.1.76 AND reviewed:true", "e-PhaZ_EC3.1.1.76"),
    ("ePhaZ_pname", 'protein_name:"poly(3-hydroxybutyrate) depolymerase" AND reviewed:true', "e-PhaZ_pname"),
    # --- 胞内/胞外 PHA 解聚酶（关键物种，gene + organism_name 双保险）---
    ("iPhaZ_Cnec", 'gene:phaZ AND organism_name:"Cupriavidus necator"', "i-PhaZ_Cupriavidus_necator"),
    ("iPhaZ_Rrub", 'gene:phaZ AND organism_name:"Rhodospirillum rubrum"', "i-PhaZ_Rhodospirillum_rubrum"),
    ("iPhaZ_Pput", 'gene:phaZ AND organism_name:"Pseudomonas putida"', "i-PhaZ_Pseudomonas_putida"),
    ("ePhaZ_Plemo", 'gene:phaZ AND organism_name:"Paucimonas lemoignei"', "e-PhaZ_Paucimonas_lemoignei"),
    ("ePhaZ_Cacido", 'gene:phaZ AND organism_name:"Comamonas acidovorans"', "e-PhaZ_Comamonas_acidorans"),
    ("ePhaZ_Ttherm", 'gene:phaZ AND organism_name:"Thermus thermophilus"', "e-PhaZ_Thermus_thermophilus"),
    ("ePhaZ_Strep", 'gene:phaZ AND organism_id:1883', "e-PhaZ_Streptomyces"),
    ("ePhaZ_Afacal", 'gene:phaZ AND organism_name:"Alcaligenes faecalis"', "e-PhaZ_Alcaligenes_faecalis"),
    ("ePhaZ_Undib", 'gene:phaZ AND organism_name:"Undibacterium"', "e-PhaZ_Undibacterium"),
    ("ePhaZ_Altero", 'gene:phaZ AND organism_name:"Alteromonas"', "e-PhaZ_Alteromonas"),
    ("iPhaZ_Bth", 'gene:phaZ AND organism_name:"Bacillus thuringiensis"', "i-PhaZ_Bacillus_thuringiensis"),
    ("iPhaZ_Smel", 'gene:phaZ AND organism_name:"Sinorhizobium meliloti"', "i-PhaZ_Sinorhizobium_meliloti"),
    ("iPhaZ_Abras", 'gene:phaZ AND organism_name:"Azospirillum brasilense"', "i-PhaZ_Azospirillum_brasilense"),
    ("ePhaZ_Sascom", 'gene:phaZ AND organism_name:"Streptomyces ascomycinicus"', "e-PhaZ_Streptomyces_ascomycinicus"),
    # --- 3HB 寡聚体水解酶（EC 3.1.1.22 hydroxybutyrate-dimer hydrolase）---
    ("OH_ec31122", "ec:3.1.1.22 AND reviewed:true", "oligomer_hydrolase_EC3.1.1.22"),
    ("OH_Ceutro", 'protein_name:"oligomer hydrolase" AND reviewed:true', "oligomer_hydrolase"),
    # --- 3HB 单体代谢（辅助）---
    ("BdhA", "ec:1.1.1.30 AND reviewed:true", "3HB_dehydrogenase_EC1.1.1.30"),
    # --- 颗粒蛋白 phasin（辅助，Pfam PF09361）---
    ("Phasin", "protein_name:phasin AND reviewed:true", "phasin"),
]

# 注意：entryType 形如 "UniProtKB reviewed (Swiss-Prot)" / "UniProtKB unreviewed (TrEMBL)"，
# 两者都含子串 "reviewed"，不能用 `"reviewed" in entryType` 判断，必须用前缀匹配。
# UniProt REST 字段名：`reviewed` 有效（返回的 JSON 键为 entryType）；`entryType`/`citation`
# 作为 fields 参数无效（会 400）。证据 PMID/DOI 由 rebuild_seeds_manifest.py 用完整 JSON 补。
FIELDS = ("accession,id,protein_name,gene_names,organism_name,lineage,ec,"
          "reviewed,length,sequence")


def fetch_stream(query: str, fields: str, retries: int = 3) -> list[dict]:
    """UniProt stream 接口：一次拉取全部命中（默认 format=json, size=500/页）"""
    params = {"query": query, "fields": fields, "format": "json", "size": "500"}
    url = STREAM
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=120)
            r.raise_for_status()
            data = r.json()
            results = data.get("results", [])
            # stream 接口一次返回全部；若被截断则用 next 指针翻页
            while data.get("next"):
                r2 = requests.get("https://rest.uniprot.org" + data["next"], timeout=120)
                r2.raise_for_status()
                data = r2.json()
                results += data.get("results", [])
            return results
        except Exception as e:
            if attempt == retries - 1:
                print(f"  [WARN] query failed after {retries} tries: {query} -> {e}", file=sys.stderr)
                return []
            time.sleep(2 * (attempt + 1))
    return []


def flatten_ec(ec_field) -> str:
    if isinstance(ec_field, list):
        return ";".join(str(e) for e in ec_field)
    return str(ec_field or "")


def is_reviewed(h) -> bool:
    """判定 UniProt reviewed（Swiss-Prot）状态。优先用 reviewed 字段，否则用 entryType 前缀。
    严禁用 `"reviewed" in entryType`（"unreviewed" 也含 "reviewed" 子串）。"""
    r = h.get("reviewed")
    if isinstance(r, bool):
        return r
    if isinstance(r, str) and r.strip().lower() in ("true", "yes", "1"):
        return True
    et = str(h.get("entryType", "") or "")
    return et.startswith("UniProtKB reviewed")


def extract_evidence(h) -> str:
    """从交叉引用提取 PMID/DOI 证据，返回 "pmid:...;doi:..." 字符串。

    注：stream 接口（fields 模式）不返回 references，此处在 stream 下通常为空；
    权威证据由 rebuild_seeds_manifest.py 用完整 JSON 的 references 补全。"""
    pmids, dois = set(), set()
    for xr in h.get("uniProtKBCrossReferences") or []:
        db = (xr.get("database") or "").lower()
        xid = (xr.get("id") or "")
        if db == "pubmed" and xid:
            pmids.add(xid)
        elif db == "doi" and xid:
            dois.add(xid)
    parts = []
    if pmids:
        parts.append("pmid:" + ";".join(sorted(pmids)))
    if dois:
        parts.append("doi:" + ";".join(sorted(dois)))
    return ";".join(parts)


def assign_split(accession: str) -> str:
    """确定性 train/validation 划分（80/20，按 accession md5 取模，可复现）。"""
    h = int(hashlib.md5(accession.encode("utf-8")).hexdigest(), 16)
    return "train" if h % 10 < 8 else "validation"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="seeds")
    ap.add_argument("--queries", default="all", help="逗号分隔的查询名子集，默认 all")
    ap.add_argument("--min-reviewed-only", action="store_true",
                    help="只保留 reviewed(Swiss-Prot) 条目，丢弃 unreviewed")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    selected = [q for q in QUERIES if args.queries == "all" or q[0] in args.queries.split(",")]
    all_rows = []
    fasta = []
    seen_acc = set()

    for name, query, label in selected:
        print(f"[*] {name}: {query}")
        hits = fetch_stream(query, FIELDS)
        print(f"    -> {len(hits)} hits")
        for h in hits:
            acc = h.get("primaryAccession") or h.get("accession", "")
            if not acc or acc in seen_acc:
                continue
            reviewed_bool = is_reviewed(h)
            if args.min_reviewed_only and not reviewed_bool:
                continue
            seen_acc.add(acc)
            seq = (h.get("sequence") or {}).get("value", "")
            org = (h.get("organism") or {}).get("scientificName", "")
            # gene names: genes: [{geneName:{value}}, {orderedLocusNames:[{value}]}, {orfNames:[...]}]
            genes = []
            for g in h.get("genes") or []:
                if g.get("geneName"):
                    genes.append(g["geneName"].get("value", ""))
                for locus in g.get("orderedLocusNames") or []:
                    genes.append(locus.get("value", ""))
                for orf in g.get("orfNames") or []:
                    genes.append(orf.get("value", ""))
            gene_str = ";".join(dict.fromkeys(genes))  # 保序去重
            # EC: proteinDescription.ecNumbers [{value}]
            pdesc = h.get("proteinDescription") or {}
            ecs = [e.get("value", "") for e in pdesc.get("ecNumbers") or []]
            pname = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
            if not pname and pdesc.get("submissionNames"):
                pname = pdesc["submissionNames"][0].get("fullName", {}).get("value", "")
            rows = {
                "accession": acc,
                "query_group": label,
                "reviewed": "true" if reviewed_bool else "false",
                "organism": org,
                "gene": gene_str,
                "ec": ";".join(ecs),
                "protein_name": pname,
                "length": h.get("sequence", {}).get("length", ""),
                "lineage": ";".join((h.get("organism") or {}).get("lineage", [])),
                "evidence": extract_evidence(h),
                "retrieval_date": dt.date.today().isoformat(),
                "split": assign_split(acc),
            }
            all_rows.append(rows)
            if seq:
                header = f">{acc}|{label}|{org}|{gene_str}|EC:{';'.join(ecs)}|rev:{rows['reviewed']}"
                fasta.append(header)
                fasta.append(seq)
        time.sleep(0.5)  # 礼貌限速

    # 写全量 FASTA + manifest
    faa = os.path.join(args.outdir, "phb_degradation_seeds.faa")
    with open(faa, "w", encoding="utf-8") as f:
        f.write("\n".join(fasta) + "\n")
    man = os.path.join(args.outdir, "seeds_manifest.tsv")
    import csv
    FIELDNAMES = ["accession", "query_group", "reviewed", "organism", "gene", "ec",
                  "protein_name", "length", "lineage", "evidence", "retrieval_date", "split"]
    with open(man, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)

    # 精选：仅 reviewed + 去冗余（按 accession 已去重）
    curated = [r for r in all_rows if r["reviewed"] == "true"]
    curated_fasta = []
    for r in curated:
        acc = r["accession"]
        # 从全量 fasta 中找对应序列
        for i, h in enumerate(fasta):
            if h.startswith(f">{acc}|"):
                curated_fasta.append(h)
                curated_fasta.append(fasta[i + 1])
                break
    cur = os.path.join(args.outdir, "seeds_curated.faa")
    with open(cur, "w", encoding="utf-8") as f:
        f.write("\n".join(curated_fasta) + "\n")

    n_train = sum(1 for r in all_rows if r["split"] == "train")
    n_val = sum(1 for r in all_rows if r["split"] == "validation")
    print(f"\n[DONE] total unique: {len(all_rows)}; curated(reviewed): {len(curated)}")
    print(f"  train/validation: {n_train}/{n_val} (80/20, accession md5 取模)")
    print(f"  FASTA(all):   {faa}")
    print(f"  FASTA(curated): {cur}")
    print(f"  manifest:     {man}")
    json.dump({
        "total": len(all_rows), "curated": len(curated),
        "train": n_train, "validation": n_val,
        "min_reviewed_only": bool(args.min_reviewed_only),
        "queries": [q[0] for q in selected],
    }, open(os.path.join(args.outdir, "seeds_stats.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
