#!/usr/bin/env python3
"""探测 PAZy / ESTHER 数据库结构，为种子扩充做准备"""
import requests

print("=== PAZy API 探测 ===")
base = "https://www.pazy.eu/api"
for path in ["families", "plastics", "proteins", "enzymes", "entries", "proteins?limit=5"]:
    try:
        r = requests.get(f"{base}/{path}", timeout=20)
        txt = r.text
        print(f"GET {base}/{path} -> {r.status_code}, len={len(txt)}")
        if r.status_code == 200 and len(txt) < 1500:
            print("   ", txt[:600].replace("\n", " "))
    except Exception as e:
        print(f"GET {base}/{path} -> ERR {e}")

print()
print("=== ESTHER 家族页探测 ===")
for fam in ["Esterase_phb_PHAZ", "PHAZ7_phb_depolymerase"]:
    try:
        r = requests.get(f"https://bioweb.supagro.inrae.fr/ESTHER/family/{fam}", timeout=20)
        print(f"ESTHER {fam} -> {r.status_code}, len={len(r.text)}")
        if r.status_code == 200:
            # 找家族成员表格/序列链接
            import re
            accs = set(re.findall(r'[A-Z][0-9][A-Z0-9]{3,}', r.text))
            print("   candidate accessions:", list(accs)[:10])
    except Exception as e:
        print(f"ESTHER {fam} -> ERR {e}")
