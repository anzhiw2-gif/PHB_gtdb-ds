#!/usr/bin/env python3
"""获取古菌 PHA 解聚酶关键序列：M1XPT2 + PhaZh1(HFX_6463) + patatin 域古菌分布"""
import requests

def get_entry(acc):
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}", timeout=30)
    h = r.json()
    seq = (h.get("sequence") or {}).get("value", "")
    pdesc = h.get("proteinDescription") or {}
    pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
    org = (h.get("organism") or {}).get("scientificName", "")
    ecs = [e.get("value") for e in (pdesc.get("ecNumbers") or [])]
    print(f"{acc} | {pn} | {org} | len={len(seq)} | EC={ecs}")
    return seq

print("=== M1XPT2 Natronomonas moolapensis PHA depolymerase ===")
seq1 = get_entry("M1XPT2")

print("\n=== Haloferax mediterranei HFX_6463 (PhaZh1, locus tag 查询) ===")
try:
    r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                     params={"query": "locus:HFX_6463 OR locus:HFX6463", "format": "json", "size": "10"}, timeout=30)
    hits = r.json().get("results", [])
    print("hits:", len(hits))
    for h in hits:
        acc = h.get("primaryAccession", "")
        org = (h.get("organism") or {}).get("scientificName", "")
        pn = (h.get("proteinDescription") or {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        print(f"  {acc} | {org} | {pn}")
except Exception as e:
    print("ERR", e)

print("\n=== 古菌 patatin 结构域（Pfam PF01734 patatin）分布 ===")
try:
    r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                     params={"query": "family:\"Patatin-like phospholipase (PF01734)\" AND taxonomy_id:2157", "format": "json", "size": "20"}, timeout=60)
    hits = r.json().get("results", [])
    print("archaeal patatin (PF01734) hits:", len(hits))
    seen = set()
    for h in hits:
        org = (h.get("organism") or {}).get("scientificName", "").split("(")[0].strip()
        if org not in seen:
            seen.add(org)
            acc = h.get("primaryAccession", "")
            print(f"  {acc} | {org}")
except Exception as e:
    print("ERR", e)
