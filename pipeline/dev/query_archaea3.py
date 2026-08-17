#!/usr/bin/env python3
"""获取古菌 PHA 解聚酶种子：M1XPT2 序列 + NCBI 查 PhaZh1 + 输出古菌种子 FASTA"""
import requests

# 1. M1XPT2 完整信息 + 序列
r = requests.get("https://rest.uniprot.org/uniprotkb/M1XPT2", timeout=30)
h = r.json()
pdesc = h.get("proteinDescription") or {}
pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
if not pn and pdesc.get("submissionNames"):
    pn = pdesc["submissionNames"][0].get("fullName", {}).get("value", "")
print("M1XPT2 protein:", pn)
seq1 = (h.get("sequence") or {}).get("value", "")
print("M1XPT2 len:", len(seq1))
print("M1XPT2 lineage:", ";".join((h.get("organism") or {}).get("lineage", []))[-120:])

# 2. NCBI 查 PhaZh1（Haloferax mediterranei）
print("\n=== NCBI protein: PhaZh1 Haloferax mediterranei ===")
try:
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                     params={"db": "protein", "term": "phaZh1[Gene Name] AND Haloferax mediterranei[Organism]",
                             "retmax": "10", "retmode": "json"}, timeout=30)
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    print("hits:", ids)
    if ids:
        r2 = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                          params={"db": "protein", "id": ",".join(ids), "retmode": "json"}, timeout=30)
        for uid, summ in r2.json().get("result", {}).items():
            if uid == "uids":
                continue
            print(f"  {uid}: {summ.get('title','')} | {summ.get('organism','')} | len={summ.get('slen','')}")
except Exception as e:
    print("NCBI ERR:", e)

# 3. 额外查 Haloferax mediterranei patatin（UniProt 家族注释）
print("\n=== UniProt Hfx mediterranei patatin-like ===")
q = 'organism_id:523841 AND (protein_name:patatin OR protein_name:"phospholipase")'
r = requests.get("https://rest.uniprot.org/uniprotkb/stream", params={"query": q, "format": "json", "size": "20"}, timeout=60)
for h in r.json().get("results", []):
    acc = h.get("primaryAccession", "")
    pn = (h.get("proteinDescription") or {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    print(f"  {acc} | {pn}")
