#!/usr/bin/env python3
"""02e_curate_archaea.py — 古菌种子精细整理（v2 修正）
- ArchPhaZ_hydrolase: "PHB depolymerase family esterase" 等经典家族（细菌 HMM 可覆盖）
- ArchPhaZ_patatin: patatin 样解聚酶（PhaZh1 型）——含文献验证序列
  M1XPT2 (Natronomonas) + I3RBH0 (=PhaZh1, Hfx mediterranei) + 古菌 patatin 家族同源
- 剔除磷脂酶 D 型（非 PHA 解聚酶）
输出: data/seeds/v2/ArchPhaZ_hydrolase.faa, ArchPhaZ_patatin.faa（修正）
"""
import os
import re
import requests

SEED = "data/seeds/v2"
PATATIN_DB = "data/hmms/patatin_Patatin.hmm"  # Pfam Patatin（已有）
CURATED = ["M1XPT2", "I3RBH0"]

def load_fasta(path):
    seqs = {}
    with open(path) as f:
        hdr = None; buf = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr: seqs[hdr.split("|")[0]] = (hdr, "".join(buf))
                hdr = line[1:]; buf = []
            else:
                buf.append(line)
        if hdr: seqs[hdr.split("|")[0]] = (hdr, "".join(buf))
    return seqs

def main():
    # 1. 从 v2 ArchPhaZ_patatin 中分流
    src = load_fasta(os.path.join(SEED, "ArchPhaZ_patatin.faa"))
    hydro, patatin, discard = {}, {}, {}
    for acc, (hdr, seq) in src.items():
        name = hdr.split("|")[-1].lower()
        if "phospholipase d" in name:
            discard[acc] = (hdr, seq)
        elif "depolymerase" in name or "phb depolymerase family" in name or "phb depolymerase" in name:
            hydro[acc] = (hdr, seq)
        else:
            patatin[acc] = (hdr, seq)
    print(f"分流: hydrolase型={len(hydro)}, patatin型={len(patatin)}, 剔除PLD={len(discard)}")

    # 2. 加入文献验证序列
    for acc in CURATED:
        try:
            h = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}", timeout=30).json()
            seq = (h.get("sequence") or {}).get("value", "")
            org = (h.get("organism") or {}).get("scientificName", "")
            pdesc = h.get("proteinDescription") or {}
            pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
            hdr = f"{acc}|ArchPhaZ_patatin|{org}||rev:true|{pn}"
            patatin.setdefault(acc, (hdr, seq))
            print(f"  加入验证序列: {acc} ({pn}, {len(seq)}aa)")
        except Exception as e:
            print(f"  {acc} ERR {e}")

    # 3. 古菌 patatin 家族同源（UniProt: patatin-like + archaea，宽一点并过滤长度）
    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                         params={"query": '(protein_name:"patatin") AND taxonomy_id:2157',
                                 "format": "json", "size": "100"}, timeout=60)
        for h in r.json().get("results", []):
            acc = h.get("primaryAccession", "")
            seq = (h.get("sequence") or {}).get("value", "")
            org = (h.get("organism") or {}).get("scientificName", "")
            if 200 <= len(seq) <= 500 and acc not in patatin:
                hdr = f"{acc}|ArchPhaZ_patatin|{org}||rev:true|patatin-like"
                patatin[acc] = (hdr, seq)
    except Exception as e:
        print("patatin 同源查询 ERR:", e)
    print(f"patatin 家族最终: {len(patatin)} 条")

    # 4. 写文件
    for fam, d in [("ArchPhaZ_hydrolase", hydro), ("ArchPhaZ_patatin", patatin)]:
        with open(os.path.join(SEED, f"{fam}.faa"), "w", encoding="utf-8") as f:
            for hdr, seq in d.values():
                f.write(">" + hdr + "\n" + seq + "\n")
        print(f"  {fam}.faa: {len(d)} 条")

if __name__ == "__main__":
    main()
