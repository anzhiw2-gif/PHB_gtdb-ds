#!/usr/bin/env python3
"""
02c_collect_archaea_seeds.py — 收集古菌 PHB 降解种子（patatin 样解聚酶 + PhaJ + BdhA）
依据文献（PMID 25710370 PhaZh1; PMID 27052994 PhaJ1; UniProt M1XPT2）：
  - ArchPhaZ: 古菌 patatin 样 PHA 解聚酶（PhaZh1 同源，类磷脂酶折叠）
  - PhaJ: 烯酰-CoA 水合酶（颗粒相关，动员通路）
  - 古菌 BdhA 与 phasin（补充）
输出: data/seeds/archaea_seeds.faa + data/seeds/archaea_manifest.tsv
"""
import argparse
import csv
import os
import sys
import time

import requests

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
ARCHAEA_TAXID = 2157

QUERIES = [
    # 古菌经典 PHB 解聚酶家族酯酶（α/β 水解酶型）
    ("depol_family_archaea", f'(protein_name:"PHB depolymerase" OR protein_name:"polyhydroxyalkanoate depolymerase" OR protein_name:"poly(3-hydroxybutyrate) depolymerase" OR protein_name:"PHB depolymerase family") AND taxonomy_id:{ARCHAEA_TAXID}', "ArchPhaZ"),
    # 古菌 patatin 样解聚酶（类磷脂酶折叠，PhaZh1 型）
    ("patatin_archaea", f'(protein_name:"patatin" AND (protein_name:"depolymerase" OR protein_name:"PHA" OR protein_name:"polyhydroxy")) AND taxonomy_id:{ARCHAEA_TAXID}', "ArchPhaZ_patatin"),
    # 古菌烯酰-CoA 水合酶（动员通路，Haloferax 为代表）
    ("phaj_haloferax", 'protein_name:"enoyl-CoA hydratase" AND organism_name:"Haloferax"', "PhaJ"),
    # 古菌 3HB 脱氢酶
    ("bdh_archaea", f'ec:1.1.1.30 AND taxonomy_id:{ARCHAEA_TAXID}', "BdhA"),
]

# 关键已验证种子（文献/UniProt 直指）
CURATED = {
    "M1XPT2": "ArchPhaZ_patatin",   # Natronomonas moolapensis PHB depolymerase (reviewed)
    "I3RBH0": "ArchPhaZ_patatin",   # Haloferax mediterranei patatin-like phospholipase (=PhaZh1, 321aa)
}


def fetch(query: str, size: int = 200) -> list:
    try:
        r = requests.get(STREAM, params={"query": query, "format": "json", "size": str(size)}, timeout=60)
        return r.json().get("results", [])
    except Exception as e:
        print(f"  ERR {query[:50]}: {e}", file=sys.stderr)
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/seeds")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fasta = []
    manifest = []
    seen = set()

    def add(acc, fam, seq, org, pn, rev):
        if not acc or acc in seen or not seq:
            return
        seen.add(acc)
        fasta.append(f">{acc}|{fam}|{org}|rev:{rev}|{pn[:60]}")
        fasta.append(seq)
        manifest.append({"accession": acc, "family": fam, "organism": org,
                         "reviewed": rev, "protein_name": pn})

    for qname, query, fam in QUERIES:
        hits = fetch(query)
        print(f"[*] {qname}: {len(hits)} hits")
        for h in hits:
            acc = h.get("primaryAccession", "")
            seq = (h.get("sequence") or {}).get("value", "")
            org = (h.get("organism") or {}).get("scientificName", "")
            pdesc = h.get("proteinDescription") or {}
            pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
            if not pn and pdesc.get("submissionNames"):
                pn = pdesc["submissionNames"][0].get("fullName", {}).get("value", "")
            rev = "true" if "reviewed" in str(h.get("entryType", "")) else "false"
            add(acc, fam, seq, org, pn, rev)
        time.sleep(0.4)

    # 关键验证序列（即使重复也不影响）
    for acc, fam in CURATED.items():
        try:
            h = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}", timeout=30).json()
            seq = (h.get("sequence") or {}).get("value", "")
            org = (h.get("organism") or {}).get("scientificName", "")
            pdesc = h.get("proteinDescription") or {}
            pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
            if not pn and pdesc.get("submissionNames"):
                pn = pdesc["submissionNames"][0].get("fullName", {}).get("value", "")
            rev = "true" if "reviewed" in str(h.get("entryType", "")) else "false"
            if acc in seen:
                print(f"  {acc} 已在集合中")
            add(acc, fam, seq, org, pn, rev)
        except Exception as e:
            print(f"  {acc} ERR {e}")

    with open(os.path.join(args.outdir, "archaea_seeds.faa"), "w", encoding="utf-8") as f:
        f.write("\n".join(fasta) + "\n")
    with open(os.path.join(args.outdir, "archaea_manifest.tsv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "family", "organism", "reviewed", "protein_name"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(manifest)

    from collections import Counter
    cnt = Counter(m["family"] for m in manifest)
    print("\n[DONE] archaeal seeds:", len(seen))
    for fam, n in cnt.most_common():
        print(f"  {fam}: {n}")
    print("  ->", os.path.join(args.outdir, "archaea_seeds.faa"))


if __name__ == "__main__":
    main()
