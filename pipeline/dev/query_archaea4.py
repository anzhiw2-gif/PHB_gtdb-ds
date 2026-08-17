#!/usr/bin/env python3
"""鉴定 Hfx mediterranei patatin 样蛋白中哪个是 PhaZh1，并生成古菌种子"""
import requests

for acc in ["I3R5B5", "I3RBH0", "I3RBJ3"]:
    try:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}", timeout=30)
        h = r.json()
        pdesc = h.get("proteinDescription") or {}
        pn = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
        if not pn and pdesc.get("submissionNames"):
            pn = pdesc["submissionNames"][0].get("fullName", {}).get("value", "")
        seq = (h.get("sequence") or {}).get("value", "")
        ecs = [e.get("value") for e in (pdesc.get("ecNumbers") or [])]
        genes = [g.get("geneName", {}).get("value", "") for g in (h.get("genes") or [])]
        # 评论/功能
        comments = [c.get("text", [{}])[0].get("value", "") for c in (h.get("comments") or []) if c.get("text")]
        print(f"{acc} | {pn} | len={len(seq)} | genes={genes} | EC={ecs}")
        for c in comments[:2]:
            print(f"    comment: {c[:120]}")
    except Exception as e:
        print(acc, "ERR", e)
