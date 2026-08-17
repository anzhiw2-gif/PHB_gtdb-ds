#!/usr/bin/env python3
"""
02d_collect_seeds_comprehensive.py — 全面种子收集（v2）
比 02/02b/02c 更广：按 EC 全量（含 TrEMBL）、蛋白名模式、关键属、结构域注释，
输出按"家族分类方案"组织的种子库（knowledge/family_classification.md）。

输出: data/seeds/v2/{family}.faa + v2_manifest.tsv + v2_stats.json
家族:
  ePhaZ（胞外，含子型标签）: ePhaZ_T1/T2/T3/T4（按注释尽量分型，无法分型归 ePhaZ_general）
  iPhaZ（胞内）: iPhaZ（+菌属标签）
  OH（寡聚体水解酶）
  BdhA（3HB 脱氢酶）
  PhaJ（烯酰-CoA 水合酶，古菌/动员）
  phasin（颗粒蛋白）
  ArchPhaZ_patatin（古菌 patatin 样）
  PhaC（合成酶，基因簇背景）
"""
import argparse
import csv
import json
import os
import re
import sys
import time
from collections import Counter

import requests

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,id,protein_name,gene_names,organism_name,lineage,ec,reviewed,length,sequence,protein_families,cc_function"

QUERIES = [
    # ---- 胞外 PHA 解聚酶（EC 全量含 TrEMBL）----
    ("ePhaZ_ec75_all", "ec:3.1.1.75", "ePhaZ"),
    ("ePhaZ_ec76_all", "ec:3.1.1.76", "ePhaZ"),
    ("ePhaZ_pname", 'protein_name:"PHA depolymerase" OR protein_name:"poly(3-hydroxybutyrate) depolymerase" OR protein_name:"polyhydroxybutyrate depolymerase"', "ePhaZ"),
    # ---- 胞内 PHA 解聚酶 ----
    ("iPhaZ_pname", 'protein_name:"intracellular polyhydroxyalkanoate depolymerase" OR protein_name:"intracellular PHB depolymerase"', "iPhaZ"),
    ("iPhaZ_genus", 'gene:phaZ AND (organism_name:"Cupriavidus" OR organism_name:"Rhodospirillum" OR organism_name:"Azospirillum" OR organism_name:"Sinorhizobium" OR organism_name:"Pseudomonas putida" OR organism_name:"Bacillus")', "iPhaZ"),
    # ---- 寡聚体水解酶 ----
    ("OH_ec", "ec:3.1.1.22", "OH"),
    ("OH_pname", 'protein_name:"oligomer hydrolase" OR protein_name:"hydroxybutyrate oligomer"', "OH"),
    # ---- 3HB 脱氢酶 ----
    ("BdhA", "ec:1.1.1.30", "BdhA"),
    # ---- 古菌 patatin 样解聚酶 ----
    ("patatin_arch", '(protein_name:"patatin" OR protein_name:"phospholipase") AND taxonomy_id:2157 AND (depolymerase OR polyhydroxy OR PHA OR granule)', "ArchPhaZ_patatin"),
    ("depol_arch", '(protein_name:"depolymerase") AND taxonomy_id:2157', "ArchPhaZ_patatin"),
    # ---- 烯酰-CoA 水合酶（动员，古菌+细菌）----
    ("PhaJ", 'protein_name:"enoyl-CoA hydratase" AND (organism_name:"Haloferax" OR organism_name:"Cupriavidus" OR organism_name:"Rhodospirillum")', "PhaJ"),
    # ---- phasin ----
    ("phasin", "protein_name:phasin", "phasin"),
    # ---- PHA 合酶（基因簇背景）----
    ("PhaC", "ec:2.3.1.- AND protein_name:synthase AND (protein_name:polyhydroxy OR protein_name:PHA)", "PhaC"),
]


def fetch(query: str, size: int = 500) -> list:
    try:
        r = requests.get(STREAM, params={"query": query, "fields": FIELDS,
                                         "format": "json", "size": str(size)}, timeout=90)
        data = r.json()
        res = data.get("results", [])
        while data.get("next") and len(res) < size * 4:
            r2 = requests.get("https://rest.uniprot.org" + data["next"], timeout=90)
            data = r2.json()
            res += data.get("results", [])
        return res
    except Exception as e:
        print(f"  ERR {query[:60]}: {e}", file=sys.stderr)
        return []


def parse(h):
    acc = h.get("primaryAccession", "")
    seq = (h.get("sequence") or {}).get("value", "")
    org = (h.get("organism") or {}).get("scientificName", "")
    pdesc = h.get("proteinDescription") or {}
    pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
    if not pn and pdesc.get("submissionNames"):
        pn = pdesc["submissionNames"][0].get("fullName", {}).get("value", "")
    ecs = [e.get("value", "") for e in pdesc.get("ecNumbers") or []]
    genes = []
    for g in h.get("genes") or []:
        if g.get("geneName"):
            genes.append(g["geneName"].get("value", ""))
    fams = [f.get("value", "") for f in h.get("proteinFamilies") or []]
    rev = "true" if "reviewed" in str(h.get("entryType", "")) else "false"
    return {"accession": acc, "organism": org, "protein_name": pn, "ec": ";".join(ecs),
            "gene": ";".join(dict.fromkeys(genes)), "families": ";".join(fams),
            "reviewed": rev, "sequence": seq, "length": len(seq)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/seeds/v2")
    ap.add_argument("--min-len", type=int, default=100)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fam_seqs = {}   # family -> list of (header, meta)
    seen = set()

    for qname, query, fam in QUERIES:
        hits = fetch(query)
        print(f"[*] {qname}: {len(hits)} hits")
        for h in hits:
            m = parse(h)
            if not m["accession"] or m["accession"] in seen or not m["sequence"]:
                continue
            if len(m["sequence"]) < args.min_len:
                continue
            seen.add(m["accession"])
            header = (f">{m['accession']}|{fam}|{m['organism']}|{m['gene']}|"
                      f"EC:{m['ec']}|rev:{m['reviewed']}|{m['protein_name'][:50]}")
            fam_seqs.setdefault(fam, []).append((header, m))
        time.sleep(0.4)

    # 去冗余（同家族内按序列 100% 去重 + 记录统计）
    manifest = []
    stats = {}
    for fam, items in fam_seqs.items():
        uniq = {}
        for header, m in items:
            uniq.setdefault(m["sequence"], (header, m))
        out = os.path.join(args.outdir, f"{fam}.faa")
        with open(out, "w", encoding="utf-8") as f:
            for seq, (header, m) in uniq.items():
                f.write(header + "\n" + seq + "\n")
                manifest.append({"accession": m["accession"], "family": fam, "organism": m["organism"],
                                 "gene": m["gene"], "ec": m["ec"], "reviewed": m["reviewed"],
                                 "protein_name": m["protein_name"], "families": m["families"]})
        stats[fam] = {"raw": len(items), "unique": len(uniq)}
        print(f"  {fam}: raw={len(items)} unique={len(uniq)} -> {out}")

    with open(os.path.join(args.outdir, "v2_manifest.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "family", "organism", "gene", "ec",
                                          "reviewed", "protein_name", "families"], delimiter="\t")
        w.writeheader()
        w.writerows(manifest)
    json.dump(stats, open(os.path.join(args.outdir, "v2_stats.json"), "w"), indent=2)
    print("\n[DONE]")
    for fam, s in stats.items():
        print(f"  {fam}: {s['unique']} unique")
    print("  total unique:", len(seen))


if __name__ == "__main__":
    main()
