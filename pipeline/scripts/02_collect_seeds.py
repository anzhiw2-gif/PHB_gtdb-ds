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

FIELDS = "accession,id,protein_name,gene_names,organism_name,lineage,ec,reviewed,length,sequence"


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="seeds")
    ap.add_argument("--queries", default="all", help="逗号分隔的查询名子集，默认 all")
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
                "reviewed": "true" if "reviewed" in str(h.get("entryType", "")) else "false",
                "organism": org,
                "gene": gene_str,
                "ec": ";".join(ecs),
                "protein_name": pname,
                "length": h.get("sequence", {}).get("length", ""),
                "lineage": ";".join((h.get("organism") or {}).get("lineage", [])),
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
    with open(man, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()) if all_rows else ["accession"])
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

    print(f"\n[DONE] total unique: {len(all_rows)}; curated(reviewed): {len(curated)}")
    print(f"  FASTA(all):   {faa}")
    print(f"  FASTA(curated): {cur}")
    print(f"  manifest:     {man}")
    json.dump({"total": len(all_rows), "curated": len(curated), "queries": [q[0] for q in selected]},
              open(os.path.join(args.outdir, "seeds_stats.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
