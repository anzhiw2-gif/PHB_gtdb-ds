#!/usr/bin/env python3
"""查询 UniProt 中古菌 PHA 降解相关序列：PhaZh1 (patatin) / PhaJ / BdhA / PhaP"""
import requests

def q(query, size=50, tag=""):
    try:
        r = requests.get("https://rest.uniprot.org/uniprotkb/stream",
                         params={"query": query, "format": "json", "size": str(size)}, timeout=60)
        hits = r.json().get("results", [])
        print(f"\n=== {tag} ({len(hits)} hits) ===")
        print("QUERY:", query)
        for h in hits[:12]:
            acc = h.get("primaryAccession", "")
            org = (h.get("organism") or {}).get("scientificName", "")
            pdesc = h.get("proteinDescription") or {}
            pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
            rev = "rev" if "reviewed" in str(h.get("entryType", "")) else "unrev"
            print(f"  {acc} | {rev} | {org[:45]} | {pn[:65]}")
        return hits
    except Exception as e:
        print(f"ERR {query}: {e}")
        return []

# 1. PhaZh1 / patatin depolymerase（含同源）
q('(protein_name:"patatin" AND protein_name:"depolymerase") OR gene:phaZh1', 50, "PhaZh1 patatin-like PHA depolymerase")

# 2. 古菌 PHA 解聚酶（broad）
q('(protein_name:"polyhydroxyalkanoate depolymerase" OR protein_name:"poly(3-hydroxybutyrate) depolymerase") AND (taxonomy_id:2157)', 50, "Archaeal PHA depolymerase (taxid 2157=Archaea)")

# 3. Haloferax mediterranei 全部 PHA 相关
q('organism_name:"Haloferax mediterranei" AND (depolymerase OR phaZh1 OR enoyl-CoA)', 30, "Hfx mediterranei PHA-related")

# 4. 古菌 enoyl-CoA hydratase PhaJ 同源（颗粒相关 - 难查，用 Haloferax）
q('protein_name:"enoyl-CoA hydratase" AND organism_name:"Haloferax"', 20, "Haloferax enoyl-CoA hydratase")
